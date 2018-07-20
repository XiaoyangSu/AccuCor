#' Natural Abundance correction for Carbon labeled samples
#'
#' @param path Path to input file.
#' @param sheet Name of sheet in xlsx file with columns 'compound',
#'      'formula', 'isotopelabel', and one column per sample.
#' @param compound_database Path to compound database in csv format. Only used
#'   for classic MAVEN style input when formula is not specified.
#' @param filetype Specify file type, default is to determine by file extension.
#' @param columns_to_skip Specify column heading to skip. All other columns not
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
read_elmaven <- function(path, sheet = NULL, compound_database = NULL,
                         columns_to_skip = NULL, filetype = NULL, ...) {

  if (is.null(filetype)) {
    filetype = tools::file_ext(path)
  }

  # TODO Break this out to a function
  # Classic MAVEN style import
  if (!is.null(compound_database)) {
    compounds <- readr::read_csv(compound_database)
    names(compounds)[which(tolower(names(compounds)) == 'compound')] <- "Compound"
    names(compounds)[which(tolower(names(compounds)) == 'formula')] <- "Formula"

    if (filetype %in% c("xls", "xlsx")) {
      InputDF <- readxl::read_excel(path, sheet = sheet, col_names = FALSE, ...)
    } else if (filetype %in% c("csv", "txt")) {
      InputDF <- readr::read_csv(path, col_names = FALSE, ...)
    } else if (filetype %in% c("tsv")) {
      InputDF <- readr::read_tsv(path, col_names = FALSE, ...)
    } else {
      stop(paste("Unsupported input filetype: '", filetype, "'", sep = ""))
    }

    tmpInputDF <- InputDF
    names(tmpInputDF) <- sapply(tmpInputDF[1,], as.character)
    names(tmpInputDF)[1] <- "Compound"
    tmpInputDF$IsotopeLabel <- rep(NA,nrow(tmpInputDF))

    # TODO Fix this so it uses the proper prefixes
    tmpInputDF <- dplyr::mutate(tmpInputDF, IsotopeLabel = .data$Compound)

    for (i in 1:nrow(tmpInputDF)) {

      if(grepl("PARENT$", tmpInputDF[[i,1]])) {
        # tmpInputDF$isotopeLabel[i]=0
        tmpInputDF[[i,1]] <- tmpInputDF[[i-1,1]]
      }
      else if(grepl("^\\w+-label-", tmpInputDF[[i,1]])) {
        # tmpInputDF$isotopeLabel[[i]] <- unlist(strapply(tmpInputDF[[i,1]], "(-)(\\d*)", ~ as.numeric(if (nchar(..2)) ..2)))
        tmpInputDF[[i,1]] <- tmpInputDF[[i-1,1]]  }
    }
    tmpInputDF <- dplyr::left_join(tmpInputDF, dplyr::select(compounds,
                                                   .data$Compound,
                                                   .data$Formula),
                         by = "Compound")
    tmpInputDF <- dplyr::filter(tmpInputDF, grepl("PARENT$|^\\w+-label-",
                                                  .data$IsotopeLabel))
    tmpInputDF_1<- dplyr::select(tmpInputDF, .data$Compound,
                                  .data$Formula, .data$IsotopeLabel)
    tmpInputDF_2 <- dplyr::select(tmpInputDF, -.data$Compound,
                                   -.data$Formula, -.data$IsotopeLabel)
    tmpInputDF_2 <- dplyr::mutate_all(tmpInputDF_2,
                                      dplyr::funs(as.numeric(.data$.)))
    InputDF <- dplyr::bind_cols(tmpInputDF_1, tmpInputDF_2)

    # TODO Convert sample data to numeric
    if (any(is.na(InputDF$Formula))) {
      stop(paste("The formula of",i,"is unknown",sep=" "))
    }

  # El-Maven style import
  } else {
    if (filetype %in% c("xls", "xlsx")) {
      InputDF <- readxl::read_excel(path, sheet = sheet, ...)
    } else if (filetype %in% c("csv", "txt")) {
      InputDF <- readr::read_csv(path, ...)
    } else if (filetype %in% c("tsv")) {
      InputDF <- readr::read_tsv(path, ...)
    } else {
      stop(paste("Unsupported input filetype: '", filetype, "'", sep = ""))
    }
  }

  # Remove empty rows (sometimed introduced by El-MAVEN)
  InputDF <- InputDF[!(apply(InputDF, 1, function(x) all(is.na(x)))), ]

  # Remove columns that are not needed
  if(is.null(columns_to_skip)) {
    columns_to_skip = c(
      "label", "goodPeakCount", "medMz", "medRt",
      "maxQuality", "compoundId", "expectedRtDiff", "ppmDiff", "parent"
    )
  }
  keep_col_nums <- which(!(tolower(names(InputDF)) %in% tolower(columns_to_skip)))

  # Find formula and compound columns
  formula_col_num <- which(tolower(names(InputDF)) == 'formula')
  if (rlang::is_empty(formula_col_num)) {
    stop(paste("Unable to find column 'formula' in file: ", path, sep = ""))
  }
  formula_col_name <- names(InputDF)[formula_col_num]

  # Find compound column
  compound_col_num <- which(tolower(names(InputDF)) == 'compound')
  if (rlang::is_empty(compound_col_num)) {
    stop(paste("Unable to find column 'compound' in file: ", path, sep = ""))
  }
  compound_col_name <- names(InputDF)[compound_col_num]

  # Generate label column as incremental index
  isotope_label_col_num <- which(tolower(names(InputDF)) == 'isotopelabel')
  if (rlang::is_empty(isotope_label_col_num)) {
    stop(paste("Unable to find column 'isotopelabel' in file: ", path, sep = ""))
  }
  isotope_label_col_name <- names(InputDF)[isotope_label_col_num]
  isotope_labels <- dplyr::pull(InputDF, isotope_label_col_num)

  sample_col_names <- names(InputDF)[which(!(tolower(names(InputDF)) %in%
                                               tolower(c(columns_to_skip,
                                                         "compound",
                                                         "formula",
                                                         "isotopelabel",
                                                         "metaGroupId",
                                                         "groupId"))))]
  if (length(sample_col_names) == 0) {
    stop(paste("Unable to find sample columns in file: ", path, sep = ""))
  }

  if ('metaGroupId' %in% names(InputDF) && (! all(InputDF$metaGroupId == 0))) {
    metaGroupId <- InputDF$metaGroupId
  } else {
    # Generate metaGroupId
    tmp <- dplyr::group_by(InputDF, rlang::UQ(as.name(compound_col_name)))
    metaGroupId <- dplyr::group_indices(tmp)
    # Check to be sure there is only one peak group per compound
    tmp <- dplyr::summarise(tmp, peak_group_count =
                              sum(grepl("PARENT", rlang::UQ(as.name(isotope_label_col_name)))))
    tmp <- dplyr::filter(tmp, .data$peak_group_count > 1)
    tmp <- dplyr::pull(tmp, rlang::UQ(as.name(compound_col_name)))
    if (! rlang::is_empty(tmp)) {
      stop(sprintf("Multiple peak groups detected for '%s', use metaGroupId column to distinguish peak groups (El-Maven >= 0.4)", tmp))
    }
  }


  # Create column of counts for isotopes
  # Only one type of isotope at a time is supported
  label_counts <- create_label_count(isotope_labels)

  # Combine cleaned data
  cleaned_dataframe <- dplyr::bind_cols(
    metaGroupId = as.factor(metaGroupId),
    dplyr::select(InputDF,
                  compound = compound_col_num,
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
  if (any(is.na(label_counts))) {
    bad_labels <- is.na(label_counts)
    stop(paste0("Unable to parse isotope labels: ",
               paste(isotope_labels[bad_labels], collapse = ", ")))
  }
  return(list(label_counts = label_counts, isotope = isotope_info$isotope))
}


determine_isotope <- function(isotope_labels) {
  if(any(grepl("^D-", isotope_labels))) {
    parent_prefix = "C12 PARENT"
    isotope_prefix = "D-label-"
    isotope = "D"
  } else if(any(grepl("^D2-", isotope_labels))) {
      parent_prefix = "C12 PARENT"
      isotope_prefix = "D2-label-"
      isotope = "D"
  } else if(any(grepl("^C13", isotope_labels))) {
    parent_prefix = "C12 PARENT"
    isotope_prefix = "C13-label-"
    isotope = "C"
  } else if(any(grepl("^N15", isotope_labels))) {
    parent_prefix = "C12 PARENT"
    isotope_prefix = "N15-label-"
    isotope = "N"
  } else {
    stop("Unable to determine isotope from isotopeLabel column, must be one of (C, D, or N)")
  }
  return(list(isotope = isotope, parent_prefix = parent_prefix, isotope_prefix = isotope_prefix))
}


write_output <- function(dataframe_list, path, filetype = NULL, ...) {
  if(identical(FALSE, path)) {
    return(path)
  }

  if (is.null(path)) {
    return(path)
  }

  if (is.null(filetype)) {
    filetype = tools::file_ext(path)
  }

  if (filetype %in% c("xls", "xlsx")) {
    output_path = paste(tools::file_path_sans_ext(path), "_corrected.", filetype, sep="")
    message(paste("Output written to: '", output_path, "'", sep = ""))
    writexl::write_xlsx(dataframe_list, output_path, ...)
  } else if (filetype %in% c("csv", "txt", "tsv")) {
    if (filetype %in% c("csv", "txt")) {
      write_func <- readr::write_csv
    } else if (filetype %in% c("tsv")) {
      write_func <- readr::write_tsv
    }
    for (sheet in names(dataframe_list)) {
      if (tolower(sheet) == "original") {
        next()
      }
      sheet_path <- paste(tools::file_path_sans_ext(path), "_", tolower(sheet),
                          ".", filetype, sep="")
      message(paste("Output written to: '", sheet_path, "'", sep = ""))
      write_func(dataframe_list[[sheet]], sheet_path, ...)
    }
  } else {
    stop(paste("Unsupported output filetype: '", filetype, "'", sep = ""))
  }
  return(path)
}
