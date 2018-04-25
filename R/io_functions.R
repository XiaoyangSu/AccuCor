#' Natural Abundance correction for Carbon labeled samples
#'
#' @param path Path to xlsx file.
#' @param sheet Name of sheet in xlsx file with columns 'compound',
#'      'formula', 'isotopelabel', and one column per sample.
#' @param ColumnsToSkip Specify column heading to skip. All other columns not
#'      named 'compound', 'formula', and 'isotopelabel' will be assumed to be
#'      sample names.
#' @param \dots Pass additional parameters to readxl::read_excel
#' @importFrom rlang .data
#' @return List with two dataframes, "original" which is result of read_excel
#'       and "cleaned", which is a dataframe with columns 'compound', 'formula',
#'       'isotope_label', label_index', followed by columns for each sample.
#' @export
#' @examples
#' \dontrun{
#' read_elmaven_xlsx("ExcelFile", "Sheet1")
#' }
read_elmaven_xlsx <- function(path, sheet = NULL, ColumnsToSkip = NULL, ...) {

  InputDF <- readxl::read_excel(path, sheet = sheet, ...)

  # Remove columns that are not needed
  if(is.null(ColumnsToSkip)) {
    ColumnsToSkip = c(
      "label", "metaGroupId", "groupId", "goodPeakCount", "medMz", "medRt",
      "maxQuality", "isotopeLabel", "compoundId", "expectedRtDiff", "ppmDiff",
      "parent"
    )
  }
  keep_col_nums <- which(!(tolower(names(InputDF)) %in% tolower(ColumnsToSkip)))
  formula_col_num <- which(tolower(names(InputDF)) == 'formula')
  formula_col_name <- names(InputDF)[formula_col_num]
  compound_col_num <- which(tolower(names(InputDF)) == 'compound')
  compound_col_name <- names(InputDF)[compound_col_num]
  sample_col_names <- names(InputDF)[which(!(tolower(names(InputDF)) %in% tolower(c(ColumnsToSkip, "compound", "formula"))))]

  # Generate label column as incremental index
  isotope_label_col_num <- which(tolower(names(InputDF)) == 'isotopelabel')
  isotope_label_col_name <- names(InputDF)[isotope_label_col_num]
  isotope_labels <- dplyr::pull(InputDF, isotope_label_col_num)
  label_counts <- create_label_count(isotope_labels)

  cleaned_dataframe <- dplyr::bind_cols(
    dplyr::select(InputDF, compound = compound_col_num,
                  formula = formula_col_num,
                  isotope_label = isotope_label_col_num),
    label_index = label_counts$label_counts,
    dplyr::select(InputDF, sample_col_names))

  return(list(original = InputDF, cleaned = cleaned_dataframe, isotope = label_counts$isotope))
}



create_label_count <- function(isotope_labels) {
  isotope_info <- determine_isotope(isotope_labels)
  label_counts <- stringr::str_replace_all(isotope_labels, isotope_info$isotope_prefix, "")
  label_counts <- stringr::str_replace_all(label_counts, isotope_info$parent_prefix, "0")
  label_counts <- as.integer(label_counts)
  return(list(label_counts = label_counts, isotope = isotope_info$isotope))
}


determine_isotope <- function(isotope_labels) {
  if(any(grepl("^D", isotope_labels))) {
    parent_prefix = "C12 PARENT"
    isotope_prefix = "D-label-"
    isotope = "D"
  } else if(any(grepl("^C13", isotope_labels))) {
    parent_prefix = "C12 PARENT"
    isotope_prefix = "C13-label-"
    isotope = "C"
  } else {
    stop("Unable to determine isotope from isotopeLabel column")
  }
  return(list(isotope = isotope, parent_prefix = parent_prefix, isotope_prefix = isotope_prefix))
}
