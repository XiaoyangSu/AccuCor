#' Natural Abundance correction for Carbon labeled samples
#'
#' @param InputFile Input xlsx file.
#' @param InputSheetName Name of sheet in xlsx file with columns 'compound',
#'      'formula', 'isotopelabel', and one column per sample.
#' @param ColumnsToSkip Specify column heading to skip. All other columns not
#'      named 'compound', 'formula', and 'isotopelabel' will be assumed to be
#'      sample names.
#' @param \dots Pass additional parameters to readxl::read_excel
#' @importFrom rlang .data
#' @return Dataframe with columns 'compound', 'formula', 'isotope_label',
#'      'label_index', followed by columns for each sample.
#' @examples
#' \dontrun{
#' read_elmaven_xlsx("ExcelFile", "Sheet1")
#' }
read_elmaven_xlsx <- function(InputFile, InputSheetName='Sheet1', ColumnsToSkip=NULL, ...) {

  InputDF <- readxl::read_excel(InputFile, sheet=InputSheetName, ...)

  # Remove columns that are not needed
  if(is.null(ColumnsToSkip)) {
    ColumnsToSkip = c(
      "label", "metaGroupId", "groupId", "goodPeakCount", "medMz", "medRt",
      "maxQuality", "isotopeLabel", "compoundId", "expectedRtDiff", "ppmDiff",
      "parent", "junk"
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

  tmp <- dplyr::mutate(InputDF, label_idx = stringr::str_replace(rlang::UQ(as.name(isotope_label_col_name)), "C13-label-", ""))
  tmp <- dplyr::mutate(tmp, label_idx = stringr::str_replace(.data$label_idx, "C12 PARENT", "0"))
  tmp <- dplyr::mutate(tmp, label_idx = as.numeric(.data$label_idx))
  InputDF <- dplyr::select(tmp, compound = compound_col_num,
                           formula = formula_col_num,
                           isotope_label = isotope_label_col_num,
                           label_index = .data$label_idx,
                           sample_col_names)
  # InputDF <- InputDF %>%
  #   dplyr::mutate(label_idx = stringr::str_replace(rlang::UQ(as.name(isotope_label_col_name)), "C13-label-", "")) %>%
  #   dplyr::mutate(label_idx = stringr::str_replace(label_idx, "C12 PARENT", "0")) %>%
  #   dplyr::mutate(label_idx = as.numeric(label_idx)) %>%
  #   dplyr::select(compound = compound_col_num,
  #                 formula = formula_col_num,
  #                 isotope_label = isotope_label_col_num,
  #                 label_index = .data$label_idx,
  #                 sample_col_names)

  return(InputDF)
}
