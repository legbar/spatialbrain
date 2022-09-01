side_width <- 6

# sample key sidepanel UI ----
sn_vta_UI <- function(id) {
  ns <- NS(id)
  tabPanel("SN/VTA Markers",
           sidebarLayout(
             position = "left",
             sidebarPanel(
               width = side_width,
               h3(helpText("Click on a gene to view the spatial profile...")),
               hr(),
               DT::dataTableOutput(ns("sn_vta_markers")),
               hr(),
               br(),
               # sliderTextInput(
               #   ns("padj"),
               #   label = "Filter Adj. P Value",
               #   grid = TRUE, 
               #   force_edges = TRUE, 
               #   choices = c(0.0001, 
               #               0.001, 
               #               0.01, 
               #               0.1), 
               #   selected = 0.01
               # ),
               #   
               # )
               # numericInput(
               #   ns("padj"),
               #   label = "Filter Adj. P Value",
               #   value = 0.01,
               #   min = 0,
               #   max = 1, 
               #   width = '35%'
               # ),
               # hr(),
               # br(),
               # sliderTextInput(
               #   ns("lfc"),
               #   label = "Filter Absolute Log2 Fold-change",
               #   grid = TRUE, 
               #   force_edges = TRUE, 
               #   choices = c(0, 1, 2, 3, 4, 5), 
               #   selected = 0
               # ),
               # numericInput(
               #   ns("lfc"),
               #   label = "Filter Log2 Fold-change",
               #   value = 0,
               #   min = NA,
               #   max = NA, 
               #   width = '35%'
               # ),
               # hr(),
               # br(),
               sliderTextInput(
                 ns("count_threshold"),
                 label = "Count Display Threshold (Log2)",
                 grid = TRUE, 
                 force_edges = TRUE, 
                 choices = c(0, 1, 2, 3, 4, 5, 6, 7, 8), 
                 selected = 1
               ),
             #   numericInput(
             #     ns("count_threshold"),
             #     label = "Count Display Threshold (Log2)",
             #     value = 0.5,
             #     min = NA,
             #     max = NA, 
             #     width = '35%'
             #   )
             ),
             mainPanel(
               width = 12 - side_width,
               wellPanel(plotOutput(ns("violin_plot"))),
               wellPanel(plotOutput(ns("spatial_plot"))),
               h4(helpText("Debug")),
               wellPanel(verbatimTextOutput(ns("debug")))
             )
           ))
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
                    "Log2 Fold-change",
                    "Adjusted P Value"),
      skip = 1
    ) %>%
      mutate(across(c(`Adjusted P Value`, `Log2 Fold-change`), ~ signif(.x, 3)))
    
    
    
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
      
      gene <- sn_vta_vars$results[input$sn_vta_markers_rows_selected,]$`Gene Symbol`
      counts <- read_csv(paste0("input/sn_vta/da_counts_per_gene/", gene, ".csv"))
      
      sn_vta_vars$plot_data <- counts %>%
        inner_join(metadata)
        # mutate(SCT_count = as.double(SCT_count))
      
      # SN_VTA VIOLIN PLOT
      output$violin_plot <- renderPlot({
        sn_vta_vars$plot_data %>%
          ggplot(aes(
            x = region,
            y = count,
            fill = region
          )) +
          geom_violin() +
          geom_jitter(width = 0.1) +
          scale_fill_d3() +
          scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
          theme_cowplot() +
          theme(legend.position = "none") +
          labs(x = "Region",
               y = "Count",
               title = paste0(unique(
                 sn_vta_vars$results_filtered[input$sn_vta_markers_rows_selected,]$`Gene Symbol`
               )))
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
          scale_color_d3() +
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
          labs(colour = "Count above threshold")
      })
    })
    
    output$debug <- renderPrint(sn_vta_vars$plot_data)
    
  })
}
