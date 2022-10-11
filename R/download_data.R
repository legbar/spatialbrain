download_data_UI <- function(id) {
  ns <- NS(id)
  tabPanel("Download Data",

           titlePanel(h1("Welcome to SpatialBrain", align = 'center')),
           
  )
}

download_data_SERVER <- function(id) {
  moduleServer(id, function(input, output, session) {
    # namespace ----
    ns <- session$ns
    
    
    
  })
}
