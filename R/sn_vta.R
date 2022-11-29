sn_vta_UI <- function(id) {
  ns <- NS(id)
  tabPanel("SN/VTA Markers",
           titlePanel(h1("SN/VTA Markers", align = 'left')),
           br(),
           fluidRow(
             column(width = 5,
                    offset = 1,
                    align = "center",
                    plotOutput(ns("violin_plot"), 
                               height = "450px"
                    )
             ), 
             column(width = 5, 
                    align = "center",
                    plotOutput(ns("spatial_plot"), 
                               height = "450px"
                    )
             )
           ),
           hr(),
           
           fluidRow(
             column(3,
                    h4("Label cells above a threshold:"),
                    br(),
                    sliderTextInput(
                      ns("count_threshold"),
                      label = "Count Display Threshold (Log2)",
                      grid = TRUE, 
                      force_edges = TRUE, 
                      choices = c(0, 1, 2, 3, 4, 5, 6, 7, 8), 
                      selected = 5
                    ),
                    br(),
                    # h4("Data description"),
                    # p("Description of the data and analysis methods"),
                    # tags$ul(
                    #   tags$li("item 1"), 
                    #   tags$li("item 2"),
                    #   tags$li("item 3")
                    # ),
                    style = 'border-right: 1px solid'
             ), 
             column(6, 
                    h4("Click on a gene to view the spatial profile..."),
                    DT::dataTableOutput(ns("sn_vta_markers")),
             ), 
             column(3, 
                    h4("Definitions"),
                    strong("LFC: "), span("The log2 fold-change in abundance between the cell type of interest and all other cells."),
                    br(),
                    br(),
                    strong("FDR-P "), span("The P value, adjusted for multiple comparisons (B&H)."),
                    hr(),
                    h4("Data Download"),
                    p(class = 'text-center', downloadButton(
                      ns('download_table'), 'Download Markers'
                    )),
                    p(class = 'text-center', downloadButton(
                      ns('download_spatial_plot'), 'Download Spatial Plot'
                    )),
                    p(class = 'text-center', downloadButton(
                      ns('download_marker_plot'), 'Download Marker Plot'
                    )),
                    # numericInput(ns('plot_size'), 
                    #              label = "Plot size (pixels)",
                    #              value = 800),
                    # numericInput(ns('plot_height'), 
                    #              label = "Plot height (pixels)",
                    #              value = 400),
                    # verbatimTextOutput(ns("debug")),
                    style = 'border-left: 1px solid'
             )
           )
           )
           
}

# sn_vta server ----
sn_vta_SERVER <- function(id) {
  moduleServer(id, function(input, output, session) {
    # namespace ----
    ns <- session$ns
    
    # metadata <- read_csv("input/sn_vta/da_metadata_spatial_outliers_removed.csv")
    metadata <- readRDS("input/startup/da_metadata.rds")
    # counts <- read_csv(
    #   "input/markers/sn_vs_vta_da_counts.csv",
    #   col_names = c("sample_cell_id", "Gene Symbol", "SCT_count")
    # )
    
    sn_vta_vars <- reactiveValues()
    
    sn_vta_vars$results <- read_csv(
      "input/sn_vta/sn_vta_mast.csv",
      col_names = c("Gene Symbol",
                    "LFC",
                    "FDR-P"),
      skip = 1
    ) %>%
      mutate(across(c(`LFC`, `FDR-P`), ~ signif(.x, 3)))
    
    
    
    # observeEvent(c(input$lfc, input$padj),
    #              {
    #                sn_vta_vars$results_filtered <- results %>%
    #                  # filter(abs(`Log2 Fold-change`) > input$lfc) %>%
    #                  filter(`Adjusted P Value` < input$padj)
    #              })
    
    
    
    # SN_VTA TABLE
    output$sn_vta_markers <- DT::renderDataTable({
      # sn_vta_vars$results_filtered
      sn_vta_vars$results
    },
    selection = "single",
    server = T)
    
    observeEvent(input$sn_vta_markers_rows_selected, {
      req(sn_vta_vars$results)
      if (is.null(input$sn_vta_markers_rows_selected)) {
        sn_vta_vars$gene <- sn_vta_vars$results[1,] %>% pull(`Gene Symbol`)
      } else {
        sn_vta_vars$gene <-
          sn_vta_vars$results[input$sn_vta_markers_rows_selected, ]$`Gene Symbol`
      }
      
      sn_vta_vars$counts <- read_csv(paste0("input/sn_vta/da_counts_per_gene/", sn_vta_vars$gene, ".csv"))
      sn_vta_vars$plot_data <- sn_vta_vars$counts %>%
        inner_join(metadata)
      
    }, ignoreNULL = F)
    
    # observeEvent(input$sn_vta_markers_rows_selected, {
    #   
    #   sn_vta_vars$gene <- sn_vta_vars$results[input$sn_vta_markers_rows_selected,]$`Gene Symbol`
    #   
    #   # mutate(SCT_count = as.double(SCT_count))
    #   
    # })
      
      # SN_VTA VIOLIN PLOT
      output$violin_plot <- renderPlot({
        sn_vta_vars$plot_data %>%
          ggplot(aes(
            x = region,
            y = count,
            fill = region
          )) +
          geom_jitter(width = 0.1, alpha = 0.2) +
          geom_violin() +
          scale_fill_d3() +
          scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
          theme_cowplot() +
          theme(legend.position = "none") +
          labs(x = "Region",
               y = "Count",
               title = sn_vta_vars$gene)
      })
      
      # SN_VTA SPATIAL PLOT
      output$spatial_plot <- renderPlot({
        sn_vta_vars$plot_data %>%
          arrange(count) %>%
          # filter(count > 0) %>%
          ggplot(aes(
            x = x,
            y = y,
            colour = count > input$count_threshold, 
            alpha = count > input$count_threshold
          )) +
          geom_point() +
          facet_wrap(vars(mouse_id_presentation_no_genotype),
                     scales = "free") +
          scale_y_reverse() +
          # scale_color_d3() +
          scale_color_manual(values = c("lightgrey", "red")) +
          scale_alpha_manual(values = c(0.25, 1), guide = "none") +
          theme_cowplot() +
          theme(
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            axis.title = element_blank(),
            legend.position = "top"
          ) +
          panel_border() +
          labs(colour = "Count above threshold", 
               title = sn_vta_vars$gene)
      })
    
    # output$debug <- renderPrint(sn_vta_vars$plot_data)
    
  })
}
