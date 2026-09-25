# Project memory: fidms

## PDF table extraction (ASTM D2425 / D2789 standards)

The ASTM standard PDFs in the repo reference folder (`D2425-23.pdf`, `D2789-95.pdf`, etc.)
contain complex tables (overlapping/merged cells, multi-row headers) that
`pdftotext -layout` extracts heavily garbled for — column/row alignment is not
reliable for these specific tables, even though it works fine for prose
sections. When a future task needs to re-verify or extend transcribed
table data from these PDFs (e.g. `class_response_matrix`, `sim_masses`),
don't rely solely on `pdftotext -layout`. Steps tried, and their outcomes,
when investigating this on this machine:

- `pdftotext -layout` on the table pages: garbled column alignment for wide
  multi-column tables (confirmed unreliable for Table 3/4/5 and Appendix X1
  matrices specifically).
- Checked for PDF-to-image tools to visually inspect the raw page layout
  instead of relying on text extraction:
  - `pdftoppm`, `pdftocairo`, `mutool`: **not installed** on this machine.
  - `python3 -c "import fitz"` (PyMuPDF): **not available**.
  - ImageMagick (`magick`) **is installed**, but has **no Ghostscript
    delegate**, so `magick -density 200 "file.pdf[N]" out.png` fails --
    rendering PDF pages to images is not currently possible on this machine
    without installing Ghostscript (or PyMuPDF) first.
- Fallback approach that worked reasonably well: cross-reference plain
  `pdftotext` (no `-layout`, natural reading order) against the `-layout`
  output for the same page range (`pdftotext -f N -l N file.pdf -`), using
  the reading-order stream to disambiguate which numbers belong to which
  column when `-layout`'s alignment is ambiguous. This is still
  best-effort, not a substitute for visual page inspection.
- There are `.png` files present in the references folder containing screenshots
  of all tables and table-like structures in the D2425 pdf document. These 
  can be reviewed if required.

**Recommendation for next time:** if a task requires high-confidence
transcription of one of these tables, ask the user whether Ghostscript (or
`poppler-utils`' `pdftoppm`/`pdftocairo`) can be installed, or have them
supply page images/screenshots directly, rather than spending effort on
further text-extraction workarounds.
