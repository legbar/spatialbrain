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
                    h1("Welcome to SpatialBrain"),
                    hr(),
                    h4("Results and data from Kilfeather, Hui Khoo, et al. 2022"),
                    br(),
                    h4("Spatial Transcriptomic Analyses:"), 
                    tags$ul(
                      tags$li(tags$strong("Cell Type Markers")),
                      tags$li(tags$strong("SN/VTA Markers in Dopaminergic Neurons")),
                      tags$li(tags$strong("Cell Number Changes in Age")),
                    ),
                    br(),
                    h4("TRAP Analyses:"), 
                    tags$ul(
                      tags$li(tags$strong("Dopaminergic Markers")),
                      tags$li(tags$strong("Dopaminergic Ageing")),
                      tags$li(tags$strong("Alternative Splicing in Dopaminergic Neurons")),

                    )
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
                    HTML('<center><img src="spatial.gif", height="50%"></center>')
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
    
    
  })
}
