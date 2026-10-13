test_that(".classify_family uses names then formula heuristics", {
  peaks <- data.frame(
    name = c(
      "Benzene",
      "Cyclohexane",
      "Hexane",
      "Pentane, 2-methyl-",
      "1-Hexene",
      "Ethanol",
      NA
    ),
    formula = c("C6H6", "C6H12", "C6H14", "C6H14", "C6H12", "C2H6O", "C10H22"),
    stringsAsFactors = FALSE
  )
  peaks <- cbind(peaks, .parse_formula(peaks$formula))
  fam <- .classify_family(peaks)
  expect_equal(fam[1], "Aromatic")
  expect_equal(fam[2], "Naphthene")
  expect_equal(fam[3], "Paraffin")
  expect_equal(fam[4], "Isoparaffin")
  expect_equal(fam[5], "Unsaturate")
  expect_equal(fam[6], "Oxygenate")
  expect_equal(fam[7], "Paraffin")
})

test_that("heteroatomic compounds are flagged", {
  peaks <- cbind(data.frame(name = "Thiophene"), .parse_formula("C4H4S"))
  expect_equal(.classify_family(peaks), "Aromatic")
  peaks2 <- cbind(data.frame(name = "unknown"), .parse_formula("C4H10S"))
  expect_equal(.classify_family(peaks2), "Unknown Heteroatomic")
})

test_that("solvents are detected case- and punctuation-insensitively", {
  expect_equal(
    .is_solvent(data.frame(name = c("Methanol", "methyl-alcohol", "Octane"))),
    c(TRUE, TRUE, FALSE)
  )
})
