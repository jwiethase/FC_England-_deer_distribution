required_packages <- c("shiny", "leaflet", "leafem", "raster", "rstudioapi")
for (pkg in required_packages) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
            install.packages(pkg)
      }
      library(pkg, character.only = TRUE)
}

source(file.path('source', 'functions.R'))

shinyApp(ui=fluidPage(
      shinyjs::useShinyjs(), 
      tags$style(HTML("
      .sidebar {
        position: absolute;
        z-index: 100;
        background-color: rgba(255, 255, 255, 0.9);
      }
      .main {
        margin-left: 0;
      }
    ")),
      tags$head(
            tags$style(HTML("
    .footer {
      position: fixed;
      bottom: 0;
      width: 100%;
      background-color: rgba(255, 255, 255, 0.6);
      padding: 1px 0;
      z-index: 1000;
    }
    
    .footer p {
      font-size: 10px;
      margin: 0;
    }
  "))
      ),
      # Add a footer with the funding information
      div(class = "footer",
          fluidRow(
                  column(
                        width = 12,
                        p("Project funded by the Forestry Commission, project FEE/1086", align = "left")
                  ))),
      sidebarLayout(
                  sidebarPanel(width=3,
                               class = "sidebar",
                               h3("England Deer Distribution"),
                               p("Explore the relative abundance of four deer species in England: Roe deer, Red deer, Fallow deer, and Muntjac."),
                               a("More information", href = "#", onclick = "Shiny.onInputChange('moreInfo', Math.random())"), # Add this line
                               selectInput("raster_layer", "Deer species:", 
                                           c("Roe deer", 
                                             "Red deer", 
                                             "Fallow deer", 
                                             "Muntjac")),
                               # checkboxInput("show_median_raster", label = "Median", value = TRUE),
                               # checkboxInput("show_sd_raster", label = "SD", value = FALSE),
                               radioButtons("raster_option", label = "Raster Option:", 
                                            choices = c("Median" = "median", "SD" = "sd"), selected = "median"),
                               actionButton("show_val", "Show validation plots"),
                               actionButton("show_effects", "Show effects plots"),
                               tags$hr(),
                               actionButton("show_more_options", "Display options", 
                                            icon = icon("menu-down", lib = "glyphicon"),
                                            style='padding:4px; font-size:80%'),
                               conditionalPanel(
                                     condition = "input.show_more_options % 2 == 1",
                                     sliderInput("opacity", label = "Opacity:", min = 0, value = 0.8, max = 1),
                                     selectInput("type", "Color palette:", c("turbo", "viridis", "inferno", "spectral", "bam"))
                               ),
                               
                  ),
                  mainPanel(width=12,
                            class = "main",
                            leafletOutput("map", width = validateCssUnit("100%"), height = validateCssUnit('100vh'))
                  )
            )
      ),
      server = function(input, output) {
            map_base <- reactive({
                  base()
            })
            
            output$map <- renderLeaflet({
                  map_base() %>%
                        htmlwidgets::onRender("
          function(el, x) {
            this.on('baselayerchange', function(e) {
              e.layer.bringToBack();
            })
          }
        ")
            })
            
            selected_image <- reactive({
                  switch(input$raster_layer,
                         "Roe deer" = if (input$show_val) "Capreolus_capreolusval_plots.png" else "Capreolus_capreolus_effects_plots.png",
                         "Red deer" = if (input$show_val) "Cervus_elaphusval_plots.png" else "Cervus_elaphus_effects_plots.png",
                         "Fallow deer" = if (input$show_val) "Dama_damaval_plots.png" else "Dama_dama_effects_plots.png",
                         "Muntjac" = if (input$show_val) "Muntiacus_reevesival_plots.png" else "Muntiacus_reevesi_effects_plots.png"
                  )
            })
            
            observe({
                  proxy <- leafletProxy("map")
                  
                  if (input$raster_option == "median") {
                        proxy %>%
                              clearGroup("raster")
                        if(input$raster_layer == "Roe deer"){
                              raster_path <- "Capreolus_capreolus_rast.tif"
                        }
                        if(input$raster_layer == "Red deer"){
                              raster_path <- "Cervus_elaphus_rast.tif"
                        }
                        if(input$raster_layer == "Fallow deer"){
                              raster_path <- "Dama_dama_rast.tif"
                        }
                        if(input$raster_layer == "Muntjac"){
                              raster_path <- "Muntiacus_reevesi_rast.tif"
                        }
                        r <- raster(raster_path)
                        domain_values <- c(min(values(r), na.rm = T), max(values(r), na.rm = T))
                        color_opts <- leafem::colorOptions(
                              palette = get_pal(50, input$type),
                              domain = domain_values,
                              na.color = "transparent"
                        )
                        proxy %>% 
                              addGeotiff(
                                    file = raster_path,
                                    project = TRUE,
                                    opacity = input$opacity,
                                    autozoom = FALSE,
                                    group = "raster",
                                    colorOptions = color_opts) %>% 
                              addLegend_decreasing(
                                    pal = colorNumeric(palette=get_pal(50, input$type), 
                                                       domain = domain_values), 
                                    values = domain_values,
                                    title = "Relative abundance", 
                                    # labFormat = labelFormat(transform = function(x) sort(x, decreasing = TRUE), reverse = TRUE),
                                    layerId = "rasterLegend",
                                    decreasing = TRUE)
                  } else if (input$raster_option == "sd") {
                        proxy %>%
                              clearGroup("raster")
                        if(input$raster_layer == "Roe deer"){
                              raster_path <- "Capreolus_capreolus_sd.tif"
                        }
                        if(input$raster_layer == "Red deer"){
                              raster_path <- "Cervus_elaphus_sd.tif"
                        }
                        if(input$raster_layer == "Fallow deer"){
                              raster_path <- "Dama_dama_sd.tif"
                        }
                        if(input$raster_layer == "Muntjac"){
                              raster_path <- "Muntiacus_reevesi_sd.tif"
                        }
                        r <- raster(raster_path)
                        domain_values <- c(min(values(r), na.rm = T), max(values(r), na.rm = T))
                        color_opts <- leafem::colorOptions(
                              palette = get_pal(50, input$type),
                              domain = domain_values,
                              na.color = "transparent"
                        )
                        proxy %>% 
                              addGeotiff(
                                    file = raster_path,
                                    project = TRUE,
                                    opacity = input$opacity,
                                    autozoom = FALSE,
                                    group = "raster",
                                    colorOptions = color_opts) %>% 
                              addLegend(
                                    pal = colorNumeric(palette=get_pal(50, input$type), 
                                                       domain = domain_values), 
                                    values = domain_values,
                                    title = "Median", 
                                    # labFormat = labelFormat(transform = function(x) sort(x, decreasing = TRUE)),
                                    layerId = "rasterLegend")
                  } else {
                        proxy %>%
                              clearGroup("raster") %>%
                              removeControl("rasterLegend")
                  }
            })
            
            observeEvent(input$show_val, {
                  selected_img <- switch(input$raster_layer,
                                           "Roe deer" = "Capreolus_capreolusval_plots.png",
                                           "Red deer" = "Cervus_elaphusval_plots.png",
                                           "Fallow deer" = "Dama_damaval_plots.png",
                                           "Muntjac" = "Muntiacus_reevesival_plots.png"
                  )

                  showModal(modalDialog(
                        title = "Selected Image",
                        img(src = selected_img, width = "100%"),
                        footer = tagList(
                              modalButton("Close")
                        ),
                        easyClose = TRUE,
                        size = "l",
                        fade = TRUE
                  ))
            })
            
            observeEvent(input$show_effects, {
                  selected_img <- switch(input$raster_layer,
                                           "Roe deer" = "Capreolus_capreolus_effects_plots.png",
                                           "Red deer" = "Cervus_elaphus_effects_plots.png",
                                           "Fallow deer" = "Dama_dama_effects_plots.png",
                                           "Muntjac" = "Muntiacus_reevesi_effects_plots.png"
                  )

                  showModal(modalDialog(
                        title = "Selected Image",
                        img(src = selected_img, width = "100%"),
                        footer = tagList(
                              modalButton("Close")
                        ),
                        easyClose = TRUE,
                        size = "l",
                        fade = TRUE
                  ))
            })
            
            observeEvent(input$moreInfo, {
                  showModal(modalDialog(
                        title = "More information",
                        p("Explore the relative abundance of the four most common deer species 
                        in England: Roe deer, Red deer, Fallow deer, and Muntjac.",
                          tags$br(), "This application was developed by Joris Wiethase, Department of Biology, University of York.
                          Relative abundance estimates are derived from 
                        an integrated INLA Species Distribution Model built using the package 'PointedSDMS' in R (Mostert et al. 2022, see also https://github.com/PhilipMostert/PointedSDMs). 
                         On the relative abundance score, 0 reflects absence of the species, and values closer to 1 represent areas where the species is more abundant. You can choose the deer species and raster option 
                        (Median or SD) to display on the map. Additionally, you can view validation plots and effects plots 
                        for each species by clicking the respective buttons. Customize the map's appearance by adjusting the 
                        opacity and selecting your preferred color palette. 
                        This application is funded by the Forestry Commission, project FEE/1086. Deer observation data were provided by the
                          British Trust for Ornithology, the British Deer Society, Langbein Wildlife Associates and 
                          iNaturalist.",
                          tags$br(),
                          tags$u("References:"),
                          tags$br(),
                          "Mostert, Philip, and Robert O'Hara. PointedSDMs--an R package to help facilitate the construction of integrated species distribution models. bioRxiv (2022): 2022-10."),
                        # p("Relative abundance of four deer species
                        # in England: Roe deer, Red deer, Fallow deer, and Muntjac. Relative abundance estimates are derived from
                        # an integrated INLA Species Distribution Model (closer to 1 means higher deer abundance)."),
                        footer = tagList(
                              modalButton("Close")
                        ),
                        easyClose = TRUE,
                        size = "l",
                        fade = TRUE
                  ))
            })
      }
)
