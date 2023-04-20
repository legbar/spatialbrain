home_UI <- function(id) {
  ns <- NS(id)
  tabPanel("Home",
           tags$head(
             # Note the wrapping of the string in HTML()
             tags$style(HTML("
      img {
      border: 1;
      max-width: 100%;
    }
    element.style {
      width: 75%;
    }"))
           ),
           # sidebarLayout(
           #   position = "left",
           # titlePanel(h1("SpatialBrain", align = 'center')),
           br(),
           fluidRow(
             column(4,
                    offset = 1,
                    # verbatimTextOutput(ns("debug")),
                    h1("Welcome to SpatialBrain"),
                    hr(),
                    h3("A platform to view integrated data from the Wade-Martins Laboratory of Molecular Neurodegeneration."),
                    
                    br(),
                    h4(tags$i("Results and data from Kilfeather, Khoo, et al. 2023 (Manuscript in review):")),
                    # br(),
                    h5("Spatial Transcriptomic Analyses:"), 
                    tags$ul(
                      tags$li(tags$strong("Cell Type Markers")),
                      tags$li(tags$strong("SN/VTA Markers in Dopaminergic Neurons")),
                      tags$li(tags$strong("Cell Number Changes in Age")),
                    ),
                    br(),
                    h5("TRAP Analyses:"), 
                    tags$ul(
                      tags$li(tags$strong("Dopaminergic Markers")),
                      tags$li(tags$strong("Dopaminergic Ageing")),
                      tags$li(tags$strong("Alternative Splicing in Dopaminergic Neurons")),
                    ), 
                    br(),
                    tags$i(h5("Coming soon: Integration with iPS-derived dopaminergic transcriptomic datasets.")),
                    br(),
                    imageOutput(ns("opdc"))
                    ), 
             column(4,
                    offset = 1,
                    # img(src="spatial.gif", align = "center",height='613px',width='800px')
                    # imageOutput(ns("spatial"))
                    # tags$head(tags$style(
                    #   type="text/css",
                    #   "#spatial img {max-width: 100%; width: 100%; height: auto}"
                    # )),
                    # imageOutput(ns("spatial"))
                    HTML('<center><img src="spatial.gif", height="50%"></center>'),
                    br(),
                    br(),
                    h3("Download our Poster:"),
                    br(),
                    tags$a(imageOutput(ns("poster")), href = "https://spatialbrain.s3.eu-west-2.amazonaws.com/POSTER_72ppi.png")
                    
                    
                    # br(),
                    # br(),
                    
             )
           ),
           # sidebarLayout(
           #   position = "left",
           #   sidebarPanel = NULL,
           #   # sidebarPanel(
           #   #   width = 4,
           #   #   h3(helpText("Selected Gene Information")),
           #   #   hr()
           #   # ), 
           #   mainPanel = mainPanel(
           #     use_waiter(),
           #     width = 12,
           #     # wellPanel(plotOutput(ns("ma_plot"))),
           #     # h2("Welcome to SpatialBrain.org!", align = "center"),
           #     wellPanel(
           #       h4(helpText(
           #         "Select Genes to Display More Information"
           #       ), align = "center"),
           #       img(src="spatial.gif", align = "center",height='613px',width='800px')
           #     ))
           # fluidRow(
           #   column(width = 6, 
           #          "Welcome to SpatialBrain.org")
           # ),
           # column(width = 6, 
           #        # helpText("Welcome to SpatialBrain.org"),
           #        # br(),
                  # img(src="spatial.gif", align = "center",height='613px',width='800px')
                  # )
           # mainPanel(use_waiter(),
           #           width = 12,
           #           # wellPanel(plotOutput(ns("ma_plot"))),
           #           wellPanel(h4(
           #             helpText("Welcome to SpatialBrain.org"),
           #             br(),
           #             img(src="spatial.gif", align = "center",height='613px',width='800px')
           #           ))
                     #   plotlyOutput(ns("ma_plot")),
                     #   hr(),
                     #   checkboxInput(ns("show_unchanged"),
                     #                 label = "Show Unchanged",
                     #                 value = T),
                     #   checkboxInput(ns("show_depleted"),
                     #                 label = "Show Depleted",
                     #                 value = T)
                     # ),
                     # # wellPanel(),
                     # wellPanel(verbatimTextOutput(ns("debug"))))
                     # )
                     # )
  # )
  )
}

home_SERVER <- function(id) {
  moduleServer(id, function(input, output, session) {
    # namespace ----
    ns <- session$ns
    
    output$spatial <- renderImage({
      path_to_png <- "input/startup/spatial.gif"
      
      # Get width and height of image output
      width  <- session$clientData$output_image_width
      height <- session$clientData$output_image_height
      
      list(src = path_to_png,
           contentType = "image/gif",
           width = width,
           height = height,
           alt = "Spatial Map Zooming")
      
    }, deleteFile = F)
    
    output$opdc <- renderImage({
      path_to_jpg <- "input/images/opdc-logo_downsized.jpg"
      # Get width and height of image output
      width  <- session$clientData$output_image_width
      height <- session$clientData$output_image_height
      
      list(src = path_to_jpg,
           contentType = "image/jpeg",
           width = width,
           # height = width * (1000/372),
           # height = 100,
           height = height,
           alt = "OPDC Logo")
      
    }, deleteFile = F)
    
    output$poster <- renderImage({
      path_to_jpg <- "input/images/POSTER_small.jpg"
      # Get width and height of image output
      width  <- session$clientData$output_image_width
      height <- session$clientData$output_image_height
      
      list(src = path_to_jpg,
           contentType = "image/jpeg",
           width = width,
           height = height,
           alt = "Kilfeather Poster")
      
    }, deleteFile = F)
    
    # output$debug <- renderPrint(session$clientData$output_image_height)
  })
}
