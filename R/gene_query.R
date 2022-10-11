gene_query_UI <- function(id) {
  ns <- NS(id)
  tabPanel("Gene Query",

           titlePanel(h1("Welcome to SpatialBrain", align = 'center')),
           
  )
}

gene_query_SERVER <- function(id) {
  moduleServer(id, function(input, output, session) {
    # namespace ----
    ns <- session$ns
    
    
    
  })
}
