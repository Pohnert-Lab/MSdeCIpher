

temp_env <- new.env()
source("functions/functions.R", local = temp_env)
source("functions/functions 2.R", local = temp_env)
source("functions/functions 3.R", local = temp_env)
ui <- fluidPage(
  titlePanel(img(src="MSdeCIpher.png", width="15%", style="display: block; margin-left: auto; margin-right: auto;"), windowTitle = "MSdeCIpher - R session"),
  tabsetPanel(
    tabPanel("Input",
             sidebarLayout(
               sidebarPanel = sidebarPanel(
                 numericInput("retention_tolerance", "Tolerance when matching retention times between the two runs (min)", 0.1, step = 0.1),
                 checkboxInput("use_restandards", "Use retention time standards to connect the two datasets?"),
                 conditionalPanel(condition = "input.use_restandards == true", fileInput("restandards", "File with retention time standards (.csv)", accept = ".csv")),
                 checkboxInput("raw_for_adduct", "Use raw data for adduct/fragment search?"),
                 checkboxInput("raw_for_isotopes", "Use raw data for sum formula correction? (only relevant for elements S/Cl/Br)"),
                 conditionalPanel(condition = "input.raw_for_adduct == true || input.raw_for_isotopes == true", fileInput("raw_file_for_adducts", "CI/molecular ion dataset raw file for m/z difference search and/or sum formula correction (.mzxml)", accept = ".mzxml")),
                 conditionalPanel(condition = "input.raw_for_isotopes == true", fileInput("raw_file_for_isotopes", "EI/fragment dataset raw file for sum formula correction (.mzxml)", accept = ".mzxml")),
                 conditionalPanel(condition = "input.analysis_go < 1", actionButton("analysis_go", "Start analysis", class = "btn-primary btn-lg", width = "200px")),
                 #conditionalPanel(condition = "input.analysis_go >= 1", actionButton("does_nothing", "Analysis started", class = "btn-primary btn-lg", width = "200px")),
                 uiOutput("download_results"),
                 width = 6),
               mainPanel = mainPanel(
                 column(6,
                        numericInput("mass_tolerance_ppm", "Mass accuracy of your system (ppm)", 3, step = 1),
                        fileInput("EIfile", "EI/fragment dataset (.csv)", accept = ".csv"),
                        numericInput("min_clustersize_EI", "Minimum number of m/z values a spectrum in the EI/fragment dataset must have to be included in the analysis", 20, min = 1, max = NA, step = 1),
                        fileInput("CIfile", "CI/molecular ion dataset (.csv) Adduct/fragment series search will be performed in here", accept = ".csv"),
                        numericInput("min_clustersize_CI", "Minimum number of m/z values a spectrum in the CI/molecular ion dataset must have to be included in the analysis", 20, min = 1, max = NA, step = 1),
                 ),
                 column(6,
                        textAreaInput("mass_diffs", "m/z differences to look for in the CI/molecular ion dataset", value = "-16.03130\n28.03130\n40.03130", height = "70px"),
                        numericInput("how_many_must_fit", "How many of those m/z differences need to be found for an ion to be included in the results?", 2, min = 1),
                        selectInput("additional_filter", "Additional filtering for molecular ion candidates per spectrum based on highest", list("intensity", "m/z")),
                        numericInput("topx_filter", "Only process the top X candidates per molecular ion spectrum based on above criterium", value = 3, min = 1, max = NA, step = 1),
                        textAreaInput("elements", "Element constraints for sum formula calculation (element/min/max)", value = "C 0 50\nH 0 50\nN 0 50\nO 0 50\nS 0 50\nSi 0 50\nP 0 50", height = "130px"),
                 ),
                 width = 6)
               , position = "right")
             ),
    tabPanel("Results Viewer",
             sidebarLayout(
               sidebarPanel = sidebarPanel(
                 actionButton("use_internal_results", "Use results from this session"),
                 textOutput("error_1"),
                 tags$div(
                   tags$h4("OR", align = "center")
                 ),
                 fileInput("use_external_results", "Upload results from previous session (.zip)", accept = ".zip"),
                 fileInput("EIfile2", "Previous EI/fragment dataset (.csv)", accept = ".csv"),
                 tags$h6("-----------------------", align = "center"),
                 numericInput("mz_to_search", "m/z value to search for in EI/fragment dataset", value = 50, min = 0),
                 numericInput("mz_to_search_tolerance", "+/- tolerance", value = 0.5, min = 0.00000001),
                 numericInput("rt_to_search", "retention time of spectrum in EI/fragment dataset (min)", value = 5, min = 0.0000001),
                 numericInput("rt_to_search_tolerance", "+/- tolerance", value = 0.2, min = 0.00000001),
                 actionButton("search_data", "Search", class = "btn-primary btn-lg", width = "200px")
               ),
               mainPanel = mainPanel(
                 conditionalPanel(condition = "input.search_data > 0",
                                  #htmlOutput("search_results"),
                                  #numericInput("plot_to_display", NULL, value = 1, min = 1, step = 1),
                                  #actionButton("plot_button", "Plot"),
                                  uiOutput('choose_result_to_display'),
                                  plotOutput("plot_display"),
                                  #div(
                                  #  tags$b("Possible molecular ions from CI/molecular ion dataset for this spectrum"),
                                  uiOutput("table_title"),
                                  DT::dataTableOutput("molecular_ions")
                                  #)
                                 )
               )
               , position = "left")
             )
             
  )
)



server <- local(function(input, output) {
  
  output$columns = renderUI({
    mydata = get(input$dataset)
    selectInput('columns2', 'Columns', names(mydata))
  })
  
  options(shiny.maxRequestSize=10000*1024^2)
  setwd(tempdir())
  unlink(list.files("./"), recursive = TRUE)
  observeEvent(input$analysis_go, {
  withProgress(message = "Preparing analysis", value = 0, {
    temp_env$EI_input_file <- input$EIfile$datapath
    temp_env$CI_input_file <- input$CIfile$datapath
    temp_env$min.clustersize_EI <- round(input$min_clustersize_EI)
    temp_env$min.clustersize_CI <- round(input$min_clustersize_CI)
    temp_env$RetentionStandards <- input$restandards$datapath
    temp_env$rt_tolerance <- input$retention_tolerance
    if (rt_tolerance == 0) {
      temp_env$rt_tolerance <- 1
    }
    temp_env$mass_tolerance <- input$mass_tolerance_ppm
    if (mass_tolerance == 0) {
      temp_env$mass_tolerance <- 1
    }
    temp_env$additional_filter <- input$additional_filter
    temp_env$topx_filter <- input$topx_filter
    temp_env$search_deltams <- input$mass_diffs
    temp_env$search_deltams <- as.double(stringr::str_extract_all(search_deltams, "[-]?[:digit:]+[:punct:]?[:digit:]+")[[1]])
    temp_env$how_many_must_fit <- round(input$how_many_must_fit)
    if (how_many_must_fit > length(search_deltams)) {
      temp_env$how_many_must_fit <- length(search_deltams)
    }
    elements <- input$elements
    element_symbols <- stringr::str_extract_all(elements, "[:upper:][:lower:]?")[[1]]
    element_restrictions <- stringr::str_extract_all(elements, "[:digit:]+")[[1]]
    temp_env$elements_limits <- as.list(NULL)
    for (i in 1:length(element_symbols)) {
      temp_env$elements_limits[[i]] <- c(element_symbols[i], element_restrictions[i*2-1], element_restrictions[i*2])
    }
	elements_matrix <- matrix( NA, nrow = length(elements_limits), ncol = 3)
	for (i in 1:length(elements_limits)) {
    elements_matrix[i,] <- elements_limits[[i]]
	}
	temp_env$element_definitions <- elements_matrix[,1]
    temp_env$check_raw_data <- input$raw_for_adduct
    temp_env$raw_file_check <- input$raw_for_isotopes
    temp_env$CI_raw_file <- input$raw_file_for_adducts$datapath
    temp_env$EI_raw_file <- input$raw_file_for_isotopes$datapath
	if (raw_file_check == TRUE | check_raw_data == TRUE) {
	temp_env$mass_spec_CI <- readMzXmlData::readMzXmlFile(CI_raw_file)
	} else {
	temp_env$mass_spec_CI <- NULL
	}
	if (raw_file_check == TRUE) {
	temp_env$mass_spec_EI <- readMzXmlData::readMzXmlFile(EI_raw_file)
	} else {
	temp_env$mass_spec_EI <- NULL
	}
	})
    annotate_EI_list(EI_input_file = EI_input_file, CI_input_file = CI_input_file, min.clustersize_EI = min.clustersize_EI, min.clustersize_CI = min.clustersize_CI, RetentionStandards = RetentionStandards, mass_tolerance = mass_tolerance, rt_tolerance = rt_tolerance, search_deltams = search_deltams, how_many_must_fit = how_many_must_fit, check_raw_data = check_raw_data, CI_raw_file = CI_raw_file)
    finalfunction(elements_limits = elements_limits, mass_tolerance = mass_tolerance, raw_file_check = raw_file_check, EI_raw_file = EI_raw_file, CI_raw_file = CI_raw_file)
    output$download_now <- downloadHandler("results.zip", content = function(temppath) {
      zip::zipr(zipfile = temppath, files = list.files("./annotated spectra results", pattern = ".csv", full.names = TRUE))
    }, contentType = "application/zip")
    output$download_results <- renderUI(downloadButton("download_now", "Download results"))
  })
  observeEvent(input$use_internal_results, {
    if(!is.null(input$EIfile$datapath)) {
      output$error_1 <- renderText("Session data found")
      temp_env$search_EI <- read.csv(EI_input_file)
    } else {
      output$error_1 <- renderText("No current session data found")
    }
  })
  x <- NULL
  observeEvent(input$search_data, {
    if (!is.null(input$use_external_results) & !is.null(input$EIfile2) & is.null(x)) {
      temp_env$search_EI <- read.csv(input$EIfile2$datapath)
      zip::unzip(input$use_external_results$datapath)
      dir.create("./annotated spectra results/", showWarnings = FALSE)
      file.copy(list.files("./", pattern = ".csv"), "./annotated spectra results/")
      x <- 1
    }
    mz_to_search <- input$mz_to_search
    rt_to_search <- input$rt_to_search
    mz_to_search_tolerance <- input$mz_to_search_tolerance
    rt_to_search_tolerance <- input$rt_to_search_tolerance
    temp_env$choose_list <- search_for_M_shiny(mz_to_search, rt_to_search, search_EI, mz_to_search_tolerance, rt_to_search_tolerance)
    if (!is.null(choose_list)) {
      output$choose_result_to_display <- renderUI({
        selectInput("choose_list", "Search results", choices = choose_list)
      })
    } else {
      output$choose_result_to_display <- renderUI({
        selectInput("choose_list", "Search results", choices = "No fitting spectra found")
      })
    }
  })
  observeEvent(input$choose_list, {
    if (!is.null(choose_list)) {
      plot_temp <- isolate(input$choose_list)
      output$plot_display <- renderPlot(plot_shiny(plot_temp))
      if (paste(isolate(plot_temp), ".csv", sep = "") %in% list.files("./annotated spectra results", pattern = ".csv")) {
        output$table_title <- renderUI({HTML("<hr> <b>Possible molecular ions for this EI/fragment spectrum:</b> <hr>")})
        output$molecular_ions <- DT::renderDataTable({
          write_M(isolate(plot_temp))
        }, rownames = FALSE)
      }
    } else {
      output$plot_display <- NULL
      output$table_title <- NULL
      output$molecular_ions <- NULL
    }
  })
}, temp_env)

shiny::shinyApp(ui = ui, server = server)
