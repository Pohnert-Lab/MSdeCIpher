
plot_old <- function(z) {
  EI_to_plot <- cbind(EI_spectra[which(EI_spectra$pcgroup == z),]$mz, EI_spectra[which(EI_spectra$pcgroup == z),]$into)
  colnames(EI_to_plot) <- c("mz", "into")
  PlotSpec(EI_to_plot)
  annotated_spectrum <- read.csv(paste("./annotated spectra results/", z, ".csv", sep = ""))
  annotated_spectrum <- annotated_spectrum[which(!annotated_spectrum$pcgroup == z), -(1:2)]
  print("Possible molecular ions from CI:")
  print(annotated_spectrum)
}


search_for_M <- function (mz_value, rt_value, input_file, mz_tolerance = 0.5, rt_tolerance = 0.2) {
  temp_env$EI_spectra <- read.csv(input_file)
  fitting_mz <- which((mz_value + mz_tolerance > EI_spectra$mz) & (mz_value - mz_tolerance < EI_spectra$mz))
  fitting_rt <- which((rt_value + rt_tolerance > EI_spectra$rt) & (rt_value - rt_tolerance < EI_spectra$rt))
  fitting_pcgroups <- unique(EI_spectra$pcgroup[fitting_mz[which(fitting_mz %in% fitting_rt)]])
  annotated_files <- list.files(path = "./annotated spectra results/", pattern = ".csv")
  no_hits <- which(!paste(fitting_pcgroups, ".csv", sep = "") %in% annotated_files)
  hits <- which(paste(fitting_pcgroups, ".csv", sep = "") %in% annotated_files)
  for (i in fitting_pcgroups[no_hits]) {
    print(paste("Contained in EI spectrum #", i, ". No molecular ions found.", sep = ""))
  }
  for (i in fitting_pcgroups[hits]) {
    print(paste("Contained in EI spectrum #", i, ". Molecular ion(s) found. Type plot(", i, ") to plot.", sep = ""))
  }
}

search_for_M_shiny <- function (mz_value, rt_value, input_file, mz_tolerance = 0.5, rt_tolerance = 0.2) {
  fitting_mz <- which((mz_value + mz_tolerance > input_file$mz) & (mz_value - mz_tolerance < input_file$mz))
  fitting_rt <- which((rt_value + rt_tolerance > input_file$rt) & (rt_value - rt_tolerance < input_file$rt))
  fitting_pcgroups <- unique(input_file$pcgroup[fitting_mz[which(fitting_mz %in% fitting_rt)]])
  annotated_files <- list.files(path = "./annotated spectra results/", pattern = ".csv")
  no_hits <- which(!paste(fitting_pcgroups, ".csv", sep = "") %in% annotated_files)
  hits <- which(paste(fitting_pcgroups, ".csv", sep = "") %in% annotated_files)
  string_to_display <- NULL
  choice_numbers <- NULL
  sort_string <- NULL
  if (length(fitting_pcgroups) > 0) {
	for (i in fitting_pcgroups[no_hits]) {
    string_to_display <- c(string_to_display, paste("EI spectrum #", i, " @ RT ", round(mean(input_file[input_file[, "pcgroup"] == i ,]$rt), digits = 2), " \u274C", sep = ""))
	choice_numbers <- c(choice_numbers, i)
	}
	for (i in fitting_pcgroups[hits]) {
		string_to_display <- c(string_to_display, paste("EI spectrum #", i, " @ RT ", round(mean(input_file[input_file[, "pcgroup"] == i ,]$rt), digits = 2), " \u2705", sep = ""))
		choice_numbers <- c(choice_numbers, i)
	}
	names(choice_numbers) <- string_to_display
	choice_numbers <- sort(choice_numbers)
	return(choice_numbers)
  } else {
	return(NULL)
  }
  
}

search_for_M_shiny_old2 <- function (mz_value, rt_value, input_file, mz_tolerance = 0.5, rt_tolerance = 0.2) {
  fitting_mz <- which((mz_value + mz_tolerance > input_file$mz) & (mz_value - mz_tolerance < input_file$mz))
  fitting_rt <- which((rt_value + rt_tolerance > input_file$rt) & (rt_value - rt_tolerance < input_file$rt))
  fitting_pcgroups <- unique(input_file$pcgroup[fitting_mz[which(fitting_mz %in% fitting_rt)]])
  annotated_files <- list.files(path = "./annotated spectra results/", pattern = ".csv")
  no_hits <- which(!paste(fitting_pcgroups, ".csv", sep = "") %in% annotated_files)
  hits <- which(paste(fitting_pcgroups, ".csv", sep = "") %in% annotated_files)
  string_to_display <- NULL

  if (length(fitting_pcgroups) > 0) {
	for (i in fitting_pcgroups[no_hits]) {
    string_to_display <- c(string_to_display, paste("No results in EI spectrum #", i, " @ RT ", round(mean(input_file[input_file[, "pcgroup"] == i ,]$rt), digits = 2), sep = ""))
	}
	for (i in fitting_pcgroups[hits]) {
		string_to_display <- c(string_to_display, paste("EI spectrum #", i, " @ RT ", round(mean(input_file[input_file[, "pcgroup"] == i ,]$rt), digits = 2), sep = ""))
	}
	return(string_to_display)
  } else {
	return(NULL)
  }
  
}

search_for_M_shiny_old <- function (mz_value, rt_value, input_file, mz_tolerance = 0.5, rt_tolerance = 0.2) {
  fitting_mz <- which((mz_value + mz_tolerance > input_file$mz) & (mz_value - mz_tolerance < input_file$mz))
  fitting_rt <- which((rt_value + rt_tolerance > input_file$rt) & (rt_value - rt_tolerance < input_file$rt))
  fitting_pcgroups <- unique(input_file$pcgroup[fitting_mz[which(fitting_mz %in% fitting_rt)]])
  annotated_files <- list.files(path = "./annotated spectra results/", pattern = ".csv")
  no_hits <- which(!paste(fitting_pcgroups, ".csv", sep = "") %in% annotated_files)
  hits <- which(paste(fitting_pcgroups, ".csv", sep = "") %in% annotated_files)
  string_to_display <- NULL
  vector_to_display <- NULL
  for (i in fitting_pcgroups[no_hits]) {
    string_to_display <- c(string_to_display, paste("Contained in EI spectrum #", i, " @ RT ", round(mean(input_file[input_file[, "pcgroup"] == i ,]$rt), digits = 2), ". No molecular ions found.", "<br>", sep = ""))
  }
  for (i in fitting_pcgroups[hits]) {
    string_to_display <- c(string_to_display, paste("Contained in EI spectrum #", i, " @ RT ", round(mean(input_file[input_file[, "pcgroup"] == i ,]$rt), digits = 2), ". Molecular ion(s) found.", "<br>", sep = ""))
	vector_to_display <- c(vector_to_display, i)
  }
  return(list(string_to_display, vector_to_display))
}

plot_shiny <- function(z) {
  EI_to_plot <- cbind(search_EI[which(search_EI$pcgroup == z),]$mz, search_EI[which(search_EI$pcgroup == z),]$into)
  colnames(EI_to_plot) <- c("mz", "into")
  plot(EI_to_plot, type = "h", xlab = "m/z", ylab = "Intensity", yaxs = "i")
}

write_M <- function(z) {
  annotated_spectrum <- read.csv(paste("./annotated spectra results/", z, ".csv", sep = ""))
  annotated_spectrum <- annotated_spectrum[which(!annotated_spectrum$pcgroup == z),]
  annotated_spectrum$rt <- round(annotated_spectrum$rt, 2)
  annotated_spectrum$into <- round(annotated_spectrum$into, 2)
  annotated_spectrum$mz <- round(annotated_spectrum$mz, 5)
  colnames(annotated_spectrum) <- c("m/z", "retention time", "area", "spectrum #", "sum formula", "probability (%)")
  annotated_spectrum <- annotated_spectrum[order(annotated_spectrum[,6], decreasing = TRUE),]
  return(annotated_spectrum)
}