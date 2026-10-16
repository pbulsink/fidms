library(shiny)
library(bslib)
library(fidms)

ui <- page_sidebar(
  title = "GC-MS/FID Sample Processing",
  theme = bs_theme(version = 5),
  sidebar = sidebar(
    width = 340,
    accordion(
      open = "Files",
      accordion_panel(
        "Files",
        textInput(
          "d_path",
          ".D directory path (optional)",
          placeholder = "e.g. C:/data/diesel.D"
        ),
        helpText(
          "If a .D directory is provided, its raw data.ms is used for ",
          "hydrocarbon-class MS deconvolution, and any FIDpeaks.csv/",
          "MSDpeaks.csv/MSDsearch.csv found inside it are used as a ",
          "fallback for the three uploads below. Uploaded CSVs always take",
          " precedence over files found in the .D directory."
        ),
        fileInput(
          "fid",
          "FID peaks CSV (optional if in .D directory)",
          accept = ".csv"
        ),
        fileInput(
          "msd",
          "MSD peaks CSV (optional if in .D directory)",
          accept = ".csv"
        ),
        fileInput(
          "search",
          "MSD search CSV (optional if in .D directory)",
          accept = ".csv"
        )
      ),
      accordion_panel(
        "Sample",
        textInput("sample_id", "Sample number"),
        textInput("sample_name", "Sample name"),
        textInput("project_code", "Project code"),
        textInput("analyst", "Analyst(s)", placeholder = "nrn/pbulsink"),
        numericInput("density", "Density (g/mL)", 0.85, step = 0.0001),
        numericInput("methanol", "Methanol (wt. %)", 0, step = 0.1),
        numericInput("water", "Water (wt. %)", 0, step = 0.1),
        numericInput("dilution", "Dilution factor", 1, min = 1, step = 0.1)
      ),
      accordion_panel(
        "Processing",
        selectInput(
          "engine",
          "Alignment engine",
          c("Grid search (fast, 2 detectors)" = "grid", "GCalignR" = "gcalignr")
        ),
        selectInput(
          "style",
          "Calibration style",
          c(
            "Slope only" = "slope",
            "Slope + intercept" = "slope_intercept",
            "Quadratic" = "quadratic",
            "Internal standard" = "internal"
          )
        ),
        conditionalPanel(
          "input.style == 'internal'",
          helpText(
            "Internal standardisation normalises the sample to 100% against the reference peak. A peak is suggested once the sample is processed; accept it or select another row in the Peaks tab."
          ),
          uiOutput("suggested")
        ),
        checkboxInput("include_solvents", "Include solvents in summary", FALSE),
        numericInput(
          "expected_shift",
          "Expected MS/FID shift (min)",
          -0.017,
          step = 0.001
        ),
        numericInput(
          "shift_window",
          "Shift search window (min)",
          0.2,
          step = 0.01
        ),
        numericInput("shift_step", "Shift step (min)", 0.001, step = 0.001),
        numericInput(
          "partial_threshold",
          "Match threshold (min)",
          0.03,
          step = 0.005
        ),
        numericInput(
          "merge_max_diff",
          "Peak merge max \u0394RT (min)",
          0.05,
          step = 0.005
        )
      )
    ),
    actionButton("run", "Process sample", class = "btn-primary"),
    downloadButton("dl_xlsx", "Excel report"),
    downloadButton("dl_pdf", "PDF report")
  ),
  navset_card_tab(
    nav_panel("Peaks", DT::DTOutput("peaks")),
    nav_panel(
      "Summary",
      layout_columns(
        col_widths = c(6, 6),
        card(card_header("Compound types"), tableOutput("families")),
        card(card_header("Elements"), tableOutput("elements"))
      ),
      card(card_header("QC"), verbatimTextOutput("qc"))
    ),
    nav_panel(
      "Retention index points",
      p(
        "Edit the RT/RI anchors; the polynomial is refit from whatever is here."
      ),
      DT::DTOutput("rical")
    ),
    nav_panel(
      "MS Class Features",
      p(
        "Auditable detail table of MS class sub-features detected within ",
        "each FID peak's window (shoulders/mixed-class regions). Only ",
        "populated when a .D directory with raw scan data was supplied."
      ),
      DT::DTOutput("class_features")
    )
  )
)

server <- function(input, output, session) {
  ri_points <- reactiveVal(fidms::rical_reference_points)

  observeEvent(input$rical_cell_edit, {
    info <- input$rical_cell_edit
    d <- ri_points()
    d[info$row, info$col] <- DT::coerceValue(info$value, d[[info$col]])
    ri_points(d)
  })

  params <- reactive({
    fidms::fidms_params(
      expected_shift = input$expected_shift,
      shift_window = input$shift_window,
      shift_step = input$shift_step,
      partial_threshold = input$partial_threshold,
      merge_max_diff = input$merge_max_diff
    )
  })

  info <- reactive({
    sample_info(
      sample_id = input$sample_id,
      sample_name = input$sample_name,
      project_code = input$project_code,
      analyst = input$analyst,
      density = input$density,
      methanol_wt_pct = input$methanol,
      water_wt_pct = input$water,
      dilution = input$dilution
    )
  })

  # Resolves one of the three peak CSVs: an upload always wins; otherwise
  # fall back to a same-named file inside the .D directory, if one exists.
  resolved_csv_path <- function(uploaded, d_path, filename) {
    if (!is.null(uploaded)) {
      return(uploaded$datapath)
    }
    if (nzchar(d_path)) {
      candidate <- file.path(d_path, filename)
      if (file.exists(candidate)) {
        return(candidate)
      }
    }
    NULL
  }

  # Stage one: parse, align, identify and interpolate. Kept separate from
  # quantitation so a reference peak can be proposed and reviewed first, the
  # way frmSuggested did before Calibrate ran.
  prepared <- eventReactive(input$run, {
    req(nzchar(input$sample_id))
    d_path <- trimws(input$d_path)
    has_d <- nzchar(d_path)
    if (has_d) {
      validate(need(
        dir.exists(d_path),
        paste0(".D directory not found: ", d_path)
      ))
    }

    fid_path <- resolved_csv_path(input$fid, d_path, "FIDpeaks.csv")
    msd_path <- resolved_csv_path(input$msd, d_path, "MSDpeaks.csv")
    search_path <- resolved_csv_path(input$search, d_path, "MSDsearch.csv")
    validate(
      need(
        !is.null(fid_path),
        "Provide the FID peaks CSV (upload, or a .D directory containing FIDpeaks.csv)."
      ),
      need(
        !is.null(msd_path),
        "Provide the MSD peaks CSV (upload, or a .D directory containing MSDpeaks.csv)."
      ),
      need(
        !is.null(search_path),
        "Provide the MSD search CSV (upload, or a .D directory containing MSDsearch.csv)."
      )
    )

    # Raw .D data is optional and additive: when present, it enables
    # hydrocarbon-class MS deconvolution via a scan matrix + class-ion
    # traces; the CSV-only workflow above is unchanged if omitted.
    class_traces <- NULL
    scan_matrix <- NULL
    if (has_d) {
      scan_matrix <- tryCatch(
        fidms:::.read_ms_scan_matrix(d_path),
        error = function(e) {
          showNotification(
            paste0(
              "Could not read raw MS data from the .D directory: ",
              conditionMessage(e)
            ),
            type = "warning"
          )
          NULL
        }
      )
      if (!is.null(scan_matrix)) {
        class_traces <- fidms:::.build_class_ion_traces(scan_matrix)
      }
    }

    withProgress(message = "Aligning and identifying peaks", {
      fidms::prepare_peaks(
        fid_path,
        msd_path,
        search_path,
        params = params(),
        engine = input$engine,
        ri_reference = ri_points(),
        class_traces = class_traces,
        scan_matrix = scan_matrix
      )
    })
  })

  suggestion <- reactive({
    fidms:::.suggest_reference_peak(prepared()$peaks)
  })

  # Adopt the suggested peak as the initial table selection, which the analyst
  # can override by clicking another row (the frmBrowse path).
  observeEvent(suggestion(), {
    s <- suggestion()
    if (!is.na(s$peak_index)) {
      DT::selectRows(DT::dataTableProxy("peaks"), s$peak_index)
    }
  })

  output$suggested <- renderUI({
    s <- tryCatch(suggestion(), error = function(e) NULL)
    if (is.null(s)) {
      return(helpText("Process the sample to see a suggested peak."))
    }
    if (is.na(s$peak_index)) {
      return(div(
        class = "text-danger",
        "No peak passes the internal-standard quality check; pick one manually."
      ))
    }
    tagList(
      tags$dl(
        tags$dt("Suggested peak"),
        tags$dd(s$peak_index),
        tags$dt("Compound"),
        tags$dd(s$name),
        tags$dt("Area %"),
        tags$dd(sprintf("%.3f", s$area_percent)),
        tags$dt("MS score"),
        tags$dd(s$score),
        tags$dt("Of total peaks"),
        tags$dd(s$n_peaks)
      ),
      if (!s$approved) {
        div(
          class = "text-warning",
          "The largest peak failed the quality check; the largest qualifying peak is suggested instead."
        )
      }
    )
  })

  ref_peak <- reactive({
    sel <- input$peaks_rows_selected
    if (length(sel) == 1) sel else suggestion()$peak_index
  })

  result <- reactive({
    p <- prepared()
    ref <- if (input$style == "internal") {
      r <- ref_peak()
      validate(need(
        !is.na(r),
        "Select a reference peak in the Peaks tab to use internal standardisation."
      ))
      r
    } else {
      NULL
    }
    peaks <- fidms::calibrate(
      p$peaks,
      style = input$style,
      reference_peak = ref,
      dilution = input$dilution,
      density = input$density
    )
    out <- fidms::summarize_report(
      peaks,
      info(),
      input$include_solvents,
      shift = p$shift
    )
    out$shift <- p$shift
    out
  })

  # Peak browser: shows the identified peak list as soon as alignment is done,
  # so a reference peak can be chosen; quantified columns appear once a
  # calibration style has been applied. Class-ion columns are only appended
  # when the run was given a .D directory (raw scan data), so CSV-only runs
  # keep the original nine-column table.
  peak_table <- reactive({
    fidms::total_compound_table(
      result(),
      include_class_ion = "class_primary" %in% names(prepared()$peaks)
    )
  })

  output$peaks <- DT::renderDT(
    peak_table(),
    selection = "single",
    options = list(pageLength = 25, scrollX = TRUE)
  )

  output$rical <- DT::renderDT(
    ri_points(),
    editable = TRUE,
    selection = "none",
    options = list(pageLength = 15)
  )

  output$class_features <- DT::renderDT(
    prepared()$ms_class_features,
    selection = "none",
    options = list(pageLength = 15, scrollX = TRUE)
  )

  output$families <- renderTable({
    r <- result()
    rbind(
      setNames(
        data.frame(r$families$compound_type, round(r$families$wt_pct, 1)),
        c("Compound Type", "Total Amount (wt. %)")
      ),
      setNames(
        data.frame(
          c("Methanol (external)", "Water", "TOTAL"),
          round(c(r$totals$methanol, r$totals$water, r$totals$total), 1)
        ),
        c("Compound Type", "Total Amount (wt. %)")
      )
    )
  })

  output$elements <- renderTable({
    r <- result()
    setNames(
      data.frame(
        c(r$elements$element, "Sum"),
        round(c(r$elements$wt_pct, r$totals$element_sum), 2)
      ),
      c("Element", "Calculated wt. %")
    )
  })

  output$qc <- renderPrint(str(result()$qc))

  output$dl_xlsx <- downloadHandler(
    filename = function() {
      sprintf("%s_%s_GC-MS-FID.xlsx", input$project_code, input$sample_id)
    },
    content = function(file) fidms::write_excel_report(result(), file)
  )

  output$dl_pdf <- downloadHandler(
    filename = function() {
      sprintf("%s_%s_GC-MS-FID.pdf", input$project_code, input$sample_id)
    },
    content = function(file) fidms::write_pdf_report(result(), file)
  )
}

shinyApp(ui, server)
