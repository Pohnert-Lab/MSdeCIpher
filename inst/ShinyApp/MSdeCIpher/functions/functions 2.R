finalfunction <- function(elements_limits = list(c("C",0,50),c("H",0,50),
                                               c("N",0,50),c("O",0,50),
                                               c("S",0,50), c("Si",0,50), c("P",0,50)), mass_tolerance, raw_file_check = FALSE, EI_raw_file, CI_raw_file) {
  file_names <- list.files(path = "./annotated spectra/")
  spectra_input <- as.list(NULL)
  for (h in 1:length(file_names)) {
    spectra_input[[h]] <- read.csv(paste("./annotated spectra/", file_names[h], sep =""))
  }
  withProgress(message = "Filtering isotopes", value = 0, {
  temp_env$n <- length(spectra_input)
  results_isotope_filtered <- lapply(spectra_input, filter_isotopes, mass_tolerance = mass_tolerance)
  })
  print("isotope filtering done")
  
  temp_env$n <- length(results_isotope_filtered)
  withProgress(message = "Calculating sum formulas", value = 0, {
  results_annotated <- lapply(results_isotope_filtered, build_sum_formula_tree_old_2, mass_tolerance = mass_tolerance, elements_limits = elements_limits)
  })
  dir.create("./annotated spectra results")
  k <- 1
  for(i in file_names) {
    write.csv(results_annotated[[k]], paste("./annotated spectra results/", i, sep = ""), row.names = FALSE)
    k <- k+1
  }
}

#converts the sum formulas given as character strings in the rcdk S4 object @string into a vector with integers (depicting the number of times each element is present in the sum formula) corresponding to the element_definitions vector.
#input is a list of rcdk sum formula objects, output is a list of integer vectors.
convert_strings_into_atom_counts <- function (formula, element_definitions) {
  formulas_as_atoms <- formula
  for (j in 1:length(formula)) {
    formula_separated <- stringr::str_extract_all(formula[[j]]@string, "([:upper:])([:lower:]?+)([:digit:]*)", simplify = TRUE)
    counts_separated <- stringr::str_extract_all(formula_separated, "[:digit:]*", simplify = TRUE)[,2]
    #try(counts_separated <- stringr::str_extract_all(formula_separated, "[:digit:]+", simplify = TRUE)[,1], silent = TRUE)
    elements_separated <- stringr::str_extract_all(formula_separated, "[:alpha:]+", simplify = TRUE)[,1]
    count_matrix <- rbind(elements_separated, counts_separated)
    atom_counts <- NULL
    for (i in 1:length(element_definitions)) {
      z <- which(element_definitions[i] == count_matrix[1,])
      if (length(z) == 0) {
        atom_counts[i] <- 0
      } else {
        y <- count_matrix[2,z]
        if (y == "") {
          atom_counts[i] <- 1
        } else {
          atom_counts[i] <- y
        }
      }
    }
    formulas_as_atoms[[j]] <- as.integer(atom_counts)
  }
  return(formulas_as_atoms)
}

convert_atom_counts_into_string_formula <- function (atom_vector, element_definitions) {
  output_vector <- NULL
  for (i in 1:length(element_definitions)) {
    output_vector <- c(output_vector, element_definitions[i], atom_vector[i])
  }
  return(stringr::str_c(output_vector, collapse = ""))
}

#calculates possible sum formulas of an mz value given tolerance in ppm and with defined element contraints
#rcdk function generate.formula is wrapped in try expression because it produces an error when no sum formulas match within constraints
#this function returns a list of S4 class objects with the sum formula embedded as a character string in @string
calc_sum_formulas <- function (mz, mass_tolerance, elements_limits) {
  mass_tolerance_window <- mz/1000000*mass_tolerance
  sum_formula <- as.list(NULL)
  try(sum_formula <- rcdk::generate.formula(mz+0.00055, window = mass_tolerance_window, elements = elements_limits, validation=TRUE, charge=1), silent = TRUE)
  return(sum_formula)
}

#removes isotopes from pseudospectra list
#checks each m/z entry in the list if there is a m+1 or m+2 peak present, with tolerances depending on input ppm and (hardcoded) on the nature of the isotopes that can be present
#lower limit for the m+1 peak is defined by the 15N isotope which adds 0.99703, upper limit is defined by the 13C isotope which adds 1.00335‬
#lower limit for the m+2 peak is defined by the 34S isotope which adds 1.99579‬, upper limit is defined by double 13C isotope which adds 2.00670
#all other common isotopes (and combinations) should fall between those limits
#then deletes found isotope peaks if their intensity ("into") is below 50% of that of the original m peak
#input is a single table with columns mz, into, pcgroup
#output is the same table but with less rows
filter_isotopes <- function(table_with_isotopes, mass_tolerance) {
  output <- NULL
  for (i in unique(table_with_isotopes$pcgroup)) {
    reduced_result_table <- table_with_isotopes[which(table_with_isotopes$pcgroup == i),]
    x <- 1
    while (x <= length(reduced_result_table$mz)) {
      lower_limit_isotopes <- (reduced_result_table$mz[x]+0.99703)*(1-mass_tolerance/1000000)
      upper_limit_isotopes <- (reduced_result_table$mz[x]+1.00335)*(1+mass_tolerance/1000000)
      isotope_rows <- which((reduced_result_table$mz > lower_limit_isotopes) & (reduced_result_table$mz < upper_limit_isotopes))
      if (length(isotope_rows) != 0) {
        rows_to_delete <- NULL
        for (c in 1:length(isotope_rows)) {
          if (reduced_result_table[isotope_rows[c],]$into < 0.5*reduced_result_table[x,]$into) {
            rows_to_delete <- c(rows_to_delete, isotope_rows[c])
          }
        }
        if (!is.null(rows_to_delete)) {
          reduced_result_table <- reduced_result_table[-rows_to_delete,]
        }
      }
      x <- x+1
    }
    x <- 1
    while (x <= length(reduced_result_table$mz)) {
      lower_limit_isotopes <- (reduced_result_table$mz[x]+1.99579)*(1-mass_tolerance/1000000)
      upper_limit_isotopes <- (reduced_result_table$mz[x]+2.00670)*(1+mass_tolerance/1000000)
      isotope_rows <- which((reduced_result_table$mz > lower_limit_isotopes) & (reduced_result_table$mz < upper_limit_isotopes))
      if (length(isotope_rows) != 0) {
        rows_to_delete <- NULL
        for (c in 1:length(isotope_rows)) {
          if (reduced_result_table[isotope_rows[c],]$into < 0.5*reduced_result_table[x,]$into) {
            rows_to_delete <- c(rows_to_delete, isotope_rows[c])
          }
        }
        if (!is.null(rows_to_delete)) {
          reduced_result_table <- reduced_result_table[-rows_to_delete,]
        }
      }
      x <- x+1
    }
    output <- rbind(output, reduced_result_table)
  }
  incProgress(1/n)
  return(output)
}



#old version
build_sum_formula_tree_old <- function(result_table, mass_tolerance, elements_limits) {
  spectra_groups <- unique(result_table$pcgroup)
  separated_spectra_groups <- as.list(NULL)
  fragment_tree_list <- as.list(NULL)
  j <- 1
  for (w in spectra_groups) {
    separated_spectra_groups[[j]] <- result_table[result_table$pcgroup == w,]
    j <- j+1
  }
  new_table <- NULL
  for (x in 1:length(separated_spectra_groups)) {
    if (x != 1) {
      fragment_tree_list <- fragment_tree_list_EI_base
    }
    sum_formula_column <- NULL
    probability_column <- NULL
    for (z in 1:length(separated_spectra_groups[[x]]$mz)) {
      mz <- separated_spectra_groups[[x]]$mz[z]
      result <- calc_sum_formulas(mz, mass_tolerance, elements_limits)
      correct_sum_characters <- NULL
      if (length(result) == 0) {
        sum_formula_column[z] <- NA
        probability_column[z] <- NA
      } else {
        result_counts <- convert_strings_into_atom_counts(result, element_definitions)
        if (length(fragment_tree_list) == 0) {
          fragment_tree_list <- c(fragment_tree_list, result_counts)
          for (i in 1:length(result)) {
            correct_sum_characters <- c(correct_sum_characters, result[[i]]@string)
          }
          probability_column[z] <- NA
        } else {
          formula_matches <- NULL
          for (i in 1:length(result_counts)) {
            substracted_list <- lapply(fragment_tree_list, function(x, subt_amnt) subt_amnt - x, subt_amnt = result_counts[[i]])
            logical_list <- lapply(lapply(substracted_list, function(x) x >= 0), function(h) all(h)) == TRUE
            formula_matches <- c(formula_matches, sum(logical_list, na.rm = TRUE))
          }
          correct_sum_formulas <- which(formula_matches == max(formula_matches))
          confidence_values <- formula_matches[correct_sum_formulas]/length(fragment_tree_list)*100
          probability_column[z] <- paste(round(confidence_values[1], 2), "%", sep = "", collapse = " ")
          for (i in correct_sum_formulas) {
			if (x == 1) {
				fragment_tree_list <- c(fragment_tree_list, result_counts[i])
			}
            correct_sum_characters <- c(correct_sum_characters, result[[i]]@string)
          }
        }
        sum_formula_column[z] <- stringr::str_c(correct_sum_characters, sep = " ", collapse = " ")
      }
    }
    if (x == 1) {
      fragment_tree_list_EI_base <- fragment_tree_list
    }
    new_table <- rbind(new_table, cbind(separated_spectra_groups[[x]], sum_formula_column, probability_column))
  }
  return(new_table)
}

#old version with intensity scoring and sum formula display
build_sum_formula_tree_old_2 <- function(result_table, mass_tolerance, elements_limits) {
incProgress(0, detail = paste("spectrum", result_table$pcgroup[1]))
  print(paste("begin computing spectrum", result_table$pcgroup[1], Sys.time(), sep = " "))
  if (any(c("S", "Cl", "Br") %in% element_definitions)) {
    isotope_check_list <- c("S", "Cl", "Br")[which(c("S", "Cl", "Br") %in% element_definitions)]
  } else {
    isotope_check_list <- NA
  }
  spectra_groups <- unique(result_table$pcgroup)
  separated_spectra_groups <- as.list(NULL)
  fragment_tree_list <- as.list(NULL)
  fragment_tree_list_intensities <- NULL
  fragment_tree_list_mz <- NULL
  j <- 1
  for (w in spectra_groups) {
    separated_spectra_groups[[j]] <- result_table[result_table$pcgroup == w,]
    j <- j+1
  }
  new_table <- NULL
  for (x in 1:length(separated_spectra_groups)) {
    sum_formula_column <- NULL
    probability_column <- NULL
    for (z in 1:length(separated_spectra_groups[[x]]$mz)) {
      if(length(fragment_tree_list)==0 & !(x==1)) {
		sum_formula_column[z] <- "no fragments"
        probability_column[z] <- "no fragments"
	  } else {
      mz <- separated_spectra_groups[[x]]$mz[z]
      result <- calc_sum_formulas(mz, mass_tolerance, elements_limits)
      correct_sum_characters <- NULL
      if (length(result) == 0) {
        sum_formula_column[z] <- NA
        probability_column[z] <- NA
      } else if (length(result)>10 & length(fragment_tree_list)==0) {
		sum_formula_column[z] <- "too many possible sum formulas"
        probability_column[z] <- "too many possible sum formulas"
	  } else {
        result_counts <- convert_strings_into_atom_counts(result, element_definitions)
        if (length(fragment_tree_list)==0) {
          fragment_tree_list <- c(fragment_tree_list, result_counts)
          for (i in 1:length(result)) {
            correct_sum_characters <- c(correct_sum_characters, result[[i]]@string)
            sum_formula_column[z] <- stringr::str_c(correct_sum_characters, sep = " ", collapse = " ")
            fragment_tree_list_intensities <- c(fragment_tree_list_intensities, separated_spectra_groups[[x]]$into[z])
            fragment_tree_list_mz <- c(fragment_tree_list_mz, separated_spectra_groups[[x]]$mz[z])
          }
          probability_column[z] <- NA
		} else {
          formula_matches <- NULL
		      mz_scores <- 0
          try(mz_scores <- remove_duplicates_fragment_tree_intensities(fragment_tree_list_intensities, fragment_tree_list_mz), silent = TRUE)
          for (i in 1:length(result_counts)) {
            substracted_list <- lapply(fragment_tree_list, function(x, subt_amnt) subt_amnt - x, subt_amnt = result_counts[[i]])
            logical_list <- lapply(lapply(substracted_list, function(x) x >= 0), function(h) all(h)) == TRUE
            formula_matches <- c(formula_matches, sum(mz_scores[which(logical_list)], na.rm = TRUE))
          }
		  if (raw_file_check == TRUE) {
          check_for_isotopes <- TRUE
          while ((check_for_isotopes == TRUE) & (!all(formula_matches == 0))) {
            check_for_isotopes <- FALSE
            correct_sum_formulas <- which(formula_matches == max(formula_matches, na.rm = TRUE))
            if (!is.na(isotope_check_list)) {
              for (h in 1:length(correct_sum_formulas)) {
                isotope_validity <- NULL
                for (k in isotope_check_list) {
                  isotope_number <- result_counts[[correct_sum_formulas[h]]][which(k == element_definitions)]
                  if (isotope_number > 0) {
                    rt <- separated_spectra_groups[[x]]$rt[z]
					if(x == 1){
					isotope_validity <- c(isotope_validity, check_for_isotope(rt, mz, mass_spec_EI, k, mass_tolerance))
					} else {
					isotope_validity <- c(isotope_validity, check_for_isotope(rt, mz, mass_spec_CI, k, mass_tolerance))
					}
                  }
                }
                if (FALSE %in% isotope_validity) {
                  formula_matches[correct_sum_formulas[h]] <- 0
                  check_for_isotopes <- TRUE
                }
              }
            }
          }
		  }
		  correct_sum_formulas <- which(formula_matches == max(formula_matches, na.rm = TRUE))
          confidence_values <- formula_matches[correct_sum_formulas]/sum(unique(mz_scores))*100
          probability_column[z] <- round(confidence_values[1], 2)
          for (i in correct_sum_formulas) {
            if (x == 1) {
              fragment_tree_list <- c(fragment_tree_list, result_counts[i])
              fragment_tree_list_intensities <- c(fragment_tree_list_intensities, separated_spectra_groups[[x]]$into[z])
              fragment_tree_list_mz <- c(fragment_tree_list_mz, separated_spectra_groups[[x]]$mz[z])
            }
            correct_sum_characters <- c(correct_sum_characters, result[[i]]@string)
          }
          sum_formula_column[z] <- stringr::str_c(correct_sum_characters, sep = " ", collapse = " ")
        }
	  }
      }
    }
    new_table <- rbind(new_table, cbind(separated_spectra_groups[[x]], sum_formula_column, probability_column))
  }
  colnames(new_table) <- c("mz", "rt", "into", "pcgroup", "sum formula", "probability (%)")
  incProgress(1/n, detail = paste("spectrum", result_table$pcgroup[1]))
  print(paste("finished spectrum", result_table$pcgroup[1], Sys.time(), sep = " "))
  return(new_table)
}


remove_duplicates_fragment_tree_intensities <- function(fragment_tree_list_intensities=fragment_tree_list_intensities, fragment_tree_list_mz=fragment_tree_list_mz) {
	dat <- rbind(fragment_tree_list_intensities, fragment_tree_list_mz)
    dat <- dat[,!duplicated(dat[2,])]
	relative_intensities <- fragment_tree_list_intensities/sum(dat[1,])
	return(fragment_tree_list_mz*log(relative_intensities*100+1))
}




#version with top 3 sum formula display and connected probabilities
build_sum_formula_tree_multiple_probabilities <- function(result_table, mass_tolerance, elements_limits) {
  spectra_groups <- unique(result_table$pcgroup)
  separated_spectra_groups <- as.list(NULL)
  fragment_tree_list <- as.list(NULL)
  fragment_tree_list_intensities <- NULL
  j <- 1
  for (w in spectra_groups) {
    separated_spectra_groups[[j]] <- result_table[result_table$pcgroup == w,]
    j <- j+1
  }
  new_table <- NULL
  for (x in 1:length(separated_spectra_groups)) {
    if (x != 1) {
      fragment_tree_list <- fragment_tree_list_EI_base
      fragment_tree_list_intensities <- fragment_tree_list_intensities_base
    }
    sum_formula_column <- NULL
    probability_column <- NULL
    for (z in 1:length(separated_spectra_groups[[x]]$mz)) {
      mz <- separated_spectra_groups[[x]]$mz[z]
      result <- calc_sum_formulas(mz, mass_tolerance, elements_limits)
      correct_sum_characters <- NULL
      if (length(result) == 0) {
        sum_formula_column[z] <- NA
        probability_column[z] <- NA
      } else {
        result_counts <- convert_strings_into_atom_counts(result, element_definitions)
        if (length(fragment_tree_list) == 0) {
          fragment_tree_list <- c(fragment_tree_list, result_counts)
          for (i in 1:length(result)) {
            correct_sum_characters <- c(correct_sum_characters, result[[i]]@string)
            fragment_tree_list_intensities <- c(fragment_tree_list_intensities, separated_spectra_groups[[x]]$into[z])
          }
          probability_column[z] <- NA
        } else {
          formula_matches <- NULL
          for (i in 1:length(result_counts)) {
            substracted_list <- lapply(fragment_tree_list, function(x, subt_amnt) subt_amnt - x, subt_amnt = result_counts[[i]])
            logical_list <- lapply(lapply(substracted_list, function(x) x >= 0), function(h) all(h)) == TRUE
            formula_matches <- c(formula_matches, sum(fragment_tree_list_intensities[which(logical_list)]))
          }
          correct_sum_formulas <- which(formula_matches %in% tail(sort(formula_matches), 3))
          confidence_values <- formula_matches[correct_sum_formulas]/sum(fragment_tree_list_intensities)*100
          probability_column[z] <- paste(round(max(confidence_values), 2), "%", sep = "", collapse = " ")
          for (i in correct_sum_formulas) {
            correct_sum_characters <- c(correct_sum_characters, result[[i]]@string, paste(round(formula_matches[i]/sum(fragment_tree_list_intensities)*100, 2), "%", sep = ""))
			if (x == 1) {
              fragment_tree_list <- c(fragment_tree_list, result_counts[i])
              fragment_tree_list_intensities <- c(fragment_tree_list_intensities, separated_spectra_groups[[x]]$into[z])
            }
          }
        }
        sum_formula_column[z] <- stringr::str_c(correct_sum_characters, sep = " ", collapse = " ")
      }
    }
    if (x == 1) {
      fragment_tree_list_EI_base <- fragment_tree_list
      fragment_tree_list_intensities_base <- fragment_tree_list_intensities
    }
    new_table <- rbind(new_table, cbind(separated_spectra_groups[[x]], sum_formula_column, probability_column))
  }
  return(new_table)
}

#calculates sum formulas for all mz values in a table, starting from the lowest
#when multiple possible sum formulas are returned, chooses the most likely by taking into account all previous sum formulas from the EI pseudospectrum, also returns the % of previous sum formulas supporting the new one
#returns the same table with an additional column (% of explained intensity in the EI spectrum by the possible molecular ion)
build_sum_formula_tree <- function(result_table, mass_tolerance, elements_limits) {
  spectra_groups <- unique(result_table$pcgroup)
  separated_spectra_groups <- as.list(NULL)
  fragment_tree_list <- as.list(NULL)
  fragment_tree_list_intensities <- NULL
  j <- 1
  for (w in spectra_groups) {
    separated_spectra_groups[[j]] <- result_table[result_table$pcgroup == w,]
    j <- j+1
  }
  new_table <- NULL
  for (x in 1:length(separated_spectra_groups)) {
    if (x != 1) {
      fragment_tree_list <- fragment_tree_list_EI_base
      fragment_tree_list_intensities <- fragment_tree_list_intensities_base
    }
    sum_formula_column <- NULL
    probability_column <- NULL
    for (z in 1:length(separated_spectra_groups[[x]]$mz)) {
      mz <- separated_spectra_groups[[x]]$mz[z]
      result <- calc_sum_formulas(mz, mass_tolerance, elements_limits)
      if (length(result) == 0) {
        sum_formula_column[z] <- NA
        probability_column[z] <- NA
      } else {
        result_counts <- convert_strings_into_atom_counts(result, element_definitions)
        if (length(fragment_tree_list) == 0) {
          fragment_tree_list <- c(fragment_tree_list, result_counts)
          for (i in 1:length(result)) {
            fragment_tree_list_intensities <- c(fragment_tree_list_intensities, separated_spectra_groups[[x]]$into[z])
          }
          probability_column[z] <- NA
        } else {
          formula_matches <- NULL
          for (i in 1:length(result_counts)) {
            substracted_list <- lapply(fragment_tree_list, function(x, subt_amnt) subt_amnt - x, subt_amnt = result_counts[[i]])
            logical_list <- lapply(lapply(substracted_list, function(x) x >= 0), function(h) all(h)) == TRUE
            formula_matches <- c(formula_matches, sum(fragment_tree_list_intensities[which(logical_list)]))
          }
          correct_sum_formulas <- which(formula_matches == max(formula_matches, na.rm = TRUE))
		  if (x == 1) {
		  probability_column[z] <- NA
		  } else {
          confidence_values <- formula_matches[correct_sum_formulas]/sum(fragment_tree_list_intensities)*100
          probability_column[z] <- paste(round(confidence_values[1], 2), "%", sep = "", collapse = " ")
		  }
          for (i in correct_sum_formulas) {
            if (x == 1) {
              fragment_tree_list <- c(fragment_tree_list, result_counts[i])
              fragment_tree_list_intensities <- c(fragment_tree_list_intensities, separated_spectra_groups[[x]]$into[z])
            }
          }
        }
      }
    }
    if (x == 1) {
      fragment_tree_list_EI_base <- fragment_tree_list
      fragment_tree_list_intensities_base <- fragment_tree_list_intensities
    }
    new_table <- rbind(new_table, cbind(separated_spectra_groups[[x]], probability_column))
  }
  return(new_table)
}

#checks raw data for abundance of M+2 peak of S, Cl, or Br containing ion
#needs retention time in min, m/z, isotope ("S", "Cl" or "Br"), full exact filename of mzXml file in working directory and mass tolerance in ppm as input
#output is TRUE or FALSE
check_for_isotope <- function(check_rt, check_mz, mass_spec, isotope, mass_tolerance) {
  retention_time_vector <- unlist(lapply(mass_spec, function(x) x$metaData$retentionTime))/60
  scan_to_check <- which(check_rt < retention_time_vector)[1]
  spectrum_to_check <- mass_spec[[scan_to_check]]$spectrum$mass
  if (isotope == "S") {
    isotope_plus_mass <- 1.99579
  }
  else if (isotope == "Cl") {
    isotope_plus_mass <- 1.99705
  }
  else if (isotope == "Br") {
    isotope_plus_mass <- 1.99796
  }
  else {
    print("Error undefined isotope")
  }
  lower_boundary <- check_mz+isotope_plus_mass/1000000*(1000000-mass_tolerance)
  upper_boundary <- check_mz+isotope_plus_mass/1000000*(1000000+mass_tolerance)
  return(any((lower_boundary <= spectrum_to_check)&(upper_boundary >= spectrum_to_check)))
}
