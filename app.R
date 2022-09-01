options(repos = BiocManager::repositories())
library(shiny)
library(waiter)
library(shinycssloaders)
library(shinythemes)
library(tidyverse)
library(DT)
library(ggsci)
library(cowplot)
library(shinyWidgets)
library(ggrepel)
library(plotly)


url <- "https://i.postimg.cc/FHshdD0G/spatial-full-1920.jpg"
# url <- "https://www.freecodecamp.org/news/content/images/size/w2000/2020/04/w-qjCHPZbeXCQ-unsplash.jpg"

# Define UI for application that draws a histogram
ui <- tagList(
  # shinythemes::themeSelector(),
  navbarPage(
    title = "SpatialBrain",
    selected = "Home",
    theme = shinytheme("flatly"),
    header = list(
      use_waiter(),
      waiter_preloader(html = tagList(
        spin_fading_circles(),
        h1("Loading SpatialBrain...")
      ), 
      # image = url, 
      fadeout = T)
      
    ),
    home_UI("home"),
    navbarMenu(
      title = "Spatial Transcriptomics",
      "----",
      "All Cell Types",
      markers_UI("markers"),
      sn_vta_UI("sn_vta"),
      tabPanel("Ageing"),
      tabPanel("SNCA-OVX")
    ), 
    navbarMenu(
      title = "TRAP", 
      "----", 
      "Gene-level", 
      trap_enrichment_UI("trap_enrichment"),
      tabPanel("Dopaminergic Enrichment"), 
      tabPanel("Dopaminergic Ageing"), 
      "Transcript-level", 
      tabPanel("Alternative Splicing in Dopaminergic Neurons")
    ), 
    gene_query_UI("gene_query"),
    download_data_UI("download_data")

  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  # waiter_show()

  # Sys.sleep(3)
  
  metadata_all_cells <- readRDS("input/startup/metadata_all_cells.rds")
  cell_type_names <- readRDS("input/startup/cell_type_names.rds")
  
  sn_vta_SERVER("sn_vta")
  
  markers_SERVER("markers", metadata_all_cells, cell_type_names)
  
  trap_enrichment_SERVER("trap_enrichment")
  
  waiter_hide()
  
}

# Run the application
shinyApp(ui = ui, server = server)