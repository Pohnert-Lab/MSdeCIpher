


#sorts a RT between the RT of two standards and outputs the resulting RT for another run with different retention times
calculateRT <-
function(inputRT, EIRTs, CIRTs) {
  if(inputRT %in% EIRTs) {
    outputRT <- CIRTs[match(inputRT, EIRTs)]
  } else if(inputRT > max(EIRTs) | inputRT < min(EIRTs)) {
    print("Error: Retention time out of standards range")
  }
  else {
    EIRTs <- c(EIRTs, inputRT)
    EIRTs <- sort(EIRTs)
    standardlower <- match(inputRT, EIRTs)-1
    standardhigher <- match(inputRT, EIRTs)
    RTlower <- EIRTs[match(inputRT, EIRTs)-1]
    RThigher <- EIRTs[match(inputRT, EIRTs)+1]
    relativeRT <- (inputRT-RTlower)/(RThigher-RTlower)
    outputRT <- relativeRT*(CIRTs[standardhigher]-CIRTs[standardlower])+CIRTs[standardlower]
    return(outputRT)
  }
}

#converts a RT from EI run to CI run and vice versa given two csv files with retention time standards
rtconversion <-
function(RetentionStandards, inputRT, conversiontype = "EItoCI") {
  
  retentionTable <- read.csv(RetentionStandards)
  if (conversiontype == "EItoCI") {
    EIRTs <- retentionTable[,2]
    CIRTs <- retentionTable[,3]
    rtconversion_output <- calculateRT(inputRT, EIRTs, CIRTs)
	return(rtconversion_output)
  } else {
    if (conversiontype == "CItoEI") {
      EIRTs <- retentionTable[,3]
      CIRTs <- retentionTable[,2]
      rtconversion_output <- calculateRT(inputRT, EIRTs, CIRTs)
	  return(rtconversion_output)
    } else {
      print("Error: Please specify the conversion type")
    }
  }
}

#saves individual EI spectra with RT and intensities of each fragment
separate_EI_spectra <- function(min.clustersize_EI, input_file) {
  EI_spectra <- read.csv(input_file)
  EI_spectra_list_output <- as.list(NULL)
  i <- 1
  for (i in 1:length(unique(sort(EI_spectra$pcgroup)))) {
    z <- unique(sort(EI_spectra$pcgroup))[i]
    blub <- which(EI_spectra$pcgroup %in% z)
    if (length(blub) >= min.clustersize_EI) {
      spectrum <- EI_spectra[blub,]
      EI_spectra_list_output[[z]] <- spectrum
      }
    i+1
  }
  return(EI_spectra_list_output)
}

#saves individual CI spectra with RT and intensities of each fragment
separate_CI_spectra <- function(min.clustersize_CI, input_file) {
  CI_spectra <- read.csv(input_file)
  CI_spectra_list_output <- as.list(NULL)
  for (i in 1:length(unique(sort(CI_spectra$pcgroup)))) {
    z <- unique(sort(CI_spectra$pcgroup))[i]
    blub <- which(CI_spectra$pcgroup %in% z)
    if (length(blub) >= min.clustersize_CI) {
      spectrum <- CI_spectra[blub,]
      CI_spectra_list_output[[z]] <- spectrum
      }
  }
  return(CI_spectra_list_output)
}

#creates and index of retention times for CI clusters
create_CI_index <- function(min.clustersize_CI, input_file) {
CI_spectra <- read.csv(input_file)
CI_index <- data.frame(0,0)
n <- length(unique(sort(CI_spectra$pcgroup)))
for (i in 1:length(unique(sort(CI_spectra$pcgroup)))) {
incProgress((1/n), detail = NULL)
  z <- unique(sort(CI_spectra$pcgroup))[i]
  blub <- which(CI_spectra$pcgroup %in% z)
  if (length(blub) >= min.clustersize_CI) {
    spectrum <- CI_spectra[blub,]
    CI_index[(length(CI_index[,1])+1),1] <- z
    CI_index[length(CI_index[,1]),2] <- mean(spectrum$rt)
  }
}
CI_index <- CI_index[-1,]
colnames(CI_index) <- c("pcgroup", "avg rt")
return(CI_index)
}

#creates and index of retention times for EI clusters
create_EI_index <- function(min.clustersize_EI, input_file) {
EI_spectra <- read.csv(input_file)
EI_index <- data.frame(0,0)
n <- length(unique(sort(EI_spectra$pcgroup)))
for (i in 1:length(unique(sort(EI_spectra$pcgroup)))) {
incProgress((1/n), detail = NULL)
  z <- unique(sort(EI_spectra$pcgroup))[i]
  blub <- which(EI_spectra$pcgroup %in% z)
  if (length(blub) >= min.clustersize_EI) {
    spectrum <- EI_spectra[blub,]
    EI_index[(length(EI_index[,1])+1),1] <- z
    EI_index[length(EI_index[,1]),2] <- mean(spectrum$rt)
  }
}
EI_index <- EI_index[-1,]
colnames(EI_index) <- c("pcgroup", "avg rt")
return(EI_index)
}

#function that outputs the mass difference between m/z values with tolerances (output is table with min mass in first column and max mass in second column), considering an input ppm
mass_diff <- function (low_mass_input, high_mass_input, ppm) {
  mass_difference <- high_mass_input-low_mass_input
  tolerance <- high_mass_input/1000000*ppm
  min_mass_output <- mass_difference-tolerance
  max_mass_output <- mass_difference+tolerance
  return(data.frame(min_mass_output, max_mass_output))
}

#tests if a certain masses with tolerances (table with 2 columns, first the low ends, second the high ends) fit to certain adduct masses, returns the identity of the possible adducts as character vector or NA
isAdduct_old <- function(tolerance_values) {
  output_annotations <- NULL
  adduct_identities <- data.frame(c(-16.03130, 16.03130, 44.06260, 56.06260, 28.03130, 40.03130, 12), c("-CH4", "CH4", "C3H8", "C4H8", "C2H4", "C3H4", "C"))
  colnames(adduct_identities) <- c("m/z", "adduct identity")
  check_vector <- NULL
  for (i in 1:length(adduct_identities$"m/z")) {
    check_matrix  <- data.frame(adduct_identities$"m/z"[i] < tolerance_values[2], adduct_identities$"m/z"[i] > tolerance_values[1])
    check_vector <- apply(check_matrix, 1, all)
    if (any(check_vector)) {
      output_annotations <- c(output_annotations, paste(adduct_identities$"adduct identity"[i]))
    }
  }
  if (sum(c("-CH4", "C3H4", "C2H4") %in% output_annotations) == 3)  {
    output_annotations_final <- "[M+H]+"
  } else if (sum(c("-CH4", "C3H4", "C2H4") %in% output_annotations) == 2) {
    output_annotations_final <- "[M+H]+ ?"
  } else if (sum(c("CH4", "C3H8", "C4H8") %in% output_annotations) == 3) {
    output_annotations_final <- "[M-CH3]+"
  } else if (sum(c("CH4", "C3H8", "C4H8") %in% output_annotations) == 2) {
    output_annotations_final <- "[M-CH3]+ ?"
  } else {
    output_annotations_final <- NA
  }
  if (is.null(output_annotations)) {
    return(c(NA, NA))
  } else {
    return(c(paste(output_annotations, collapse = " "), output_annotations_final))
  }
}

#tests if a certain masses with tolerances (table with 2 columns, first the low ends, second the high ends) fit to certain adduct masses, returns TRUE if at least how_many_must_fit different mass differences have been found, otherwise returns FALSE
isAdduct <- function(deltams_values, search_deltams, how_many_must_fit) {
  output_annotations <- NULL
  for (i in search_deltams) {
    check_matrix  <- data.frame(i < deltams_values[2], i > deltams_values[1])
    check_vector <- apply(check_matrix, 1, all)
    if (any(check_vector)) {
      output_annotations <- c(output_annotations, i)
    }
  }
  if (length(unique(output_annotations)) >= how_many_must_fit) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}


#scans a list of m/z values for the CI methane adduct series outputs a list of spectra with possible molecular ion annotations
#input is the minimum peaks per CI spectrum and a .csv file with all m/z peaks of all spectra with colums mz, rt, into, pcgroup
CI_annotate_old <- function(min.clustersize_CI, input_file, mass_tolerance) {
  CI_list <- separate_CI_spectra(min.clustersize_CI, input_file)
  for (i in 1:length(CI_list)) {
    if (!is.null(CI_list[[i]])) {
      ions_list <- CI_list[[i]]
      mz_values <- unique(sort(ions_list$mz))
      annotation_vectors <- data.frame(File=character(), User=character(), stringsAsFactors=FALSE) 
      colnames(annotation_vectors) <- c("adducts", "identity")
      for (j in 1:length(mz_values)) {
        diff_list <- mass_diff(mz_values[j], mz_values, mass_tolerance)
        annotation_vectors[j,] <- isAdduct(diff_list, search_deltams, how_many_must_fit)
      }
      ions_list$adducts <- annotation_vectors[,1]
      ions_list$identity <- annotation_vectors[,2]
      CI_list[[i]] <- ions_list
    }
  }
  return(CI_list)
}

#scans a list of m/z values for a defined list of mass differences (search_deltams) and outputs a list of spectra with TRUE/FALSE annotations depending on if that ion shows at least (how_many_must_fit) of those mass differences to other ions in the list
#input is the minimum peaks per CI spectrum and a .csv file with all m/z peaks of all spectra with colums mz, rt, into, pcgroup
CI_annotate <- function(min.clustersize_CI, input_file, mass_tolerance, check_raw_data, mass_spec_CI) {
  CI_list <- separate_CI_spectra(min.clustersize_CI, input_file)
  for (i in 1:length(CI_list)) {
    if (!is.null(CI_list[[i]])) {
      ions_list <- CI_list[[i]]
      mz_values <- unique(sort(ions_list$mz))
      identity <- NULL
      if (check_raw_data == TRUE) {
        retention_time_vector <- unlist(lapply(mass_spec_CI, function(x) x$metaData$retentionTime))/60
        min_vector <- which(retention_time_vector >= min(ions_list$rt))
        min_vector <- c(min_vector[1]-1, min_vector)
        max_vector <- which(retention_time_vector <= max(ions_list$rt))
        max_vector <- c(max_vector, max_vector[length(max_vector)]+1)
        scans_to_check <- which(max_vector %in% min_vector)
        spectrum_to_check <- NULL
        for (k in scans_to_check) {
          spectrum_to_check <- c(spectrum_to_check, mass_spec_CI[[k]]$spectrum$mass)
        }
        for (j in 1:length(mz_values)) {
          diff_list <- mass_diff(mz_values[j], unique(spectrum_to_check), mass_tolerance)
          identity[j] <- isAdduct(diff_list, search_deltams, how_many_must_fit)
        }
      } else {
        for (j in 1:length(mz_values)) {
          diff_list <- mass_diff(mz_values[j], mz_values, mass_tolerance)
          identity[j] <- isAdduct(diff_list, search_deltams, how_many_must_fit)
        }
      }
      ions_list <- cbind(ions_list, identity)
      CI_list[[i]] <- ions_list
    }
	incProgress(1/length(CI_list), detail = NULL)
  }
  return(CI_list)
}

#skips the annotation step, just turns out a list similar to CI_annotate, but with all identity values as TRUE
CI_no_annotate <- function(min.clustersize_CI, input_file) {
  CI_list <- separate_CI_spectra(min.clustersize_CI, input_file)
  for (i in 1:length(CI_list)) {
    if (!is.null(CI_list[[i]])) {
      ions_list <- CI_list[[i]]
      mz_values <- unique(sort(ions_list$mz))
      identity <- rep.int(TRUE, length(mz_values))
      ions_list <- cbind(ions_list, identity)
      CI_list[[i]] <- ions_list
    }
    incProgress(1/length(CI_list), detail = NULL)
  }
  return(CI_list)
}


#not in use
sort_CItoEI_spectra  <- function(EI_input_file, CI_input_file, min.clustersize_EI = 20, min.clustersize_CI = 5, RetentionStandards = NULL, rt_tolerance = 0.5) {
  EI_index <- create_EI_index(min.clustersize_EI, EI_input_file)
  CI_index <- create_CI_index(min.clustersize_CI, CI_input_file)
  EI_spectra <- separate_EI_spectra(min.clustersize_EI, EI_input_file)
  CI_spectra <- CI_annotate(min.clustersize_CI, CI_input_file, mass_tolerance, check_raw_data, mass_spec_CI)
  dir.create(paste("./annotated spectra/"))
  for (i in 1:length(EI_index$`avg rt`)) {
    if (!is.null(RetentionStandards)) {
      rt_value <- rtconversion(RetentionStandards, EI_index$`avg rt`[i], "EItoCI")
    } else {
      rt_value <- EI_index$`avg rt`[i]
    }
    if (is.numeric(rt_value)) {
      rt_search_value <- c(rt_value-rt_tolerance, rt_value+rt_tolerance)
      fitting_CI_spectra <- NULL
      for (j in 1:length(CI_index$`avg rt`)) {
        check_vector <- all(c(CI_index$`avg rt`[j] < rt_search_value[2], CI_index$`avg rt`[j] > rt_search_value[1]))
        if (check_vector) {
          fitting_CI_spectra <- c(fitting_CI_spectra, CI_index$pcgroup[j])
        }
        
      }
      
      if (!is.null(fitting_CI_spectra)) {
        dir.create(paste("./annotated spectra/annotated spectrum ", i, sep = ""))
        dir.create(paste("./annotated spectra/annotated spectrum ", i, "/EI/", sep = ""))
        dir.create(paste("./annotated spectra/annotated spectrum ", i, "/CI/", sep = ""))
        write.csv(EI_spectra[[i]], paste("./annotated spectra/annotated spectrum ", i, "/EI/", i, ".csv", sep = ""))
        for (k in 1:length(fitting_CI_spectra)) {
          annotated_spectrum <- fitting_CI_spectra[k]
          write.csv(CI_spectra[[annotated_spectrum]], paste("./annotated spectra/annotated spectrum ", i, "/CI/", annotated_spectrum, ".csv", sep = ""))
        }
        
      }
    }
  }
}

#connects all functions and outputs each EI spectrum with possible molecular ions as a seperate .csv file
annotate_EI_list  <- function(EI_input_file, CI_input_file, min.clustersize_EI = 20, min.clustersize_CI = 5, RetentionStandards = NULL, rt_tolerance = 0.1, mass_tolerance = 3, search_deltams, how_many_must_fit, check_raw_data = FALSE, CI_raw_file) {
  withProgress(message = "Creating index of spectra 1/2", value = 0, {
  EI_index <- create_EI_index(min.clustersize_EI, EI_input_file)
  })
  withProgress(message = "Creating index of spectra 2/2", value = 0, {
  CI_index <- create_CI_index(min.clustersize_CI, CI_input_file)
  })
  EI_spectra <- separate_EI_spectra(min.clustersize_EI, EI_input_file)
  withProgress(message = "Finding ions of interest", value = 0, {
    if (length(search_deltams) == 0) {
      CI_spectra <- CI_no_annotate(min.clustersize_CI, CI_input_file)
    } else {
      CI_spectra <- CI_annotate(min.clustersize_CI, CI_input_file, mass_tolerance, check_raw_data, mass_spec_CI)
    }
  })
  dir.create("./annotated spectra/")
  for (i in 1:length(EI_index$`avg rt`)) {
    x <- EI_index$pcgroup[i]
    if (!is.null(RetentionStandards)) {
      rt_value <- rtconversion(RetentionStandards, EI_index$`avg rt`[i], "EItoCI")
    } else {
      rt_value <- EI_index$`avg rt`[i]
    }
    if (is.numeric(rt_value)) {
      rt_search_value <- c(rt_value-rt_tolerance, rt_value+rt_tolerance)
      fitting_CI_spectra <- NULL
      for (j in 1:length(CI_index$`avg rt`)) {
        check_vector <- all(c(CI_index$`avg rt`[j] < rt_search_value[2], CI_index$`avg rt`[j] > rt_search_value[1]))
        if (check_vector) {
          fitting_CI_spectra <- c(fitting_CI_spectra, CI_index$pcgroup[j])
        }
      }
      if (!is.null(fitting_CI_spectra)) {
        for (k in 1:length(fitting_CI_spectra)) {
          annotated_spectrum <- CI_spectra[[fitting_CI_spectra[k]]]
          annotated_spectrum$pcgroup <- paste(annotated_spectrum$pcgroup, "CI")
          molecular_ions <- which(annotated_spectrum$identity)
          if (length(molecular_ions) > 0) {
            append_ion <- annotated_spectrum[molecular_ions,c("mz", "rt", "into", "pcgroup")]
            EI_spectra[[x]] <- rbind(EI_spectra[[x]], append_ion)
            write.csv(EI_spectra[[x]], paste("./annotated spectra/", x, ".csv", sep = ""), row.names = FALSE)
          }
        }
      }
    }
  }	
}