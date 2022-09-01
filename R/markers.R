side_width <- 6

# sample key sidepanel UI ----
markers_UI <- function(id) {
  ns <- NS(id)
  tabPanel("Spatial Markers",
           sidebarLayout(
             position = "left",
             sidebarPanel(
               width = side_width,
               h3(helpText("Select a cell type to view markers...")),
               hr(),
               selectizeInput(ns("cell_type"), 
                              label = "Select cell type", 
                              choices = NULL),
               br(),
               DT::dataTableOutput(ns("markers"))
             ),
             mainPanel(
               width = 12 - side_width,
               wellPanel(plotOutput(ns("spatial_plot"))),
               wellPanel(verbatimTextOutput(ns("debug")))
             )
           ))
}

# markers server ----
markers_SERVER <- function(id, metadata_all_cells, cell_type_names) {
  moduleServer(id, function(input, output, session) {
    # namespace ----
    ns <- session$ns
    
    
    # results <- read_csv(
    #   "input/sn_vta/sn_vta_mast.csv",
    #   col_names = c("Gene Symbol",
    #                 "Log2 Fold-change",
    #                 "Adjusted P Value"),
    #   skip = 1
    # ) %>%
    #   mutate(across(c(`Adjusted P Value`, `Log2 Fold-change`), ~ signif(.x, 3)))
    
    # sn_vta_vars <- reactiveValues()
    # 
    # observeEvent(c(input$cell_type),
    #              {
    #                sn_vta_vars$results_filtered <- results %>%
    #                  # filter(abs(`Log2 Fold-change`) > input$lfc) %>%
    #                  filter(`Adjusted P Value` < input$padj)
    #              })
    
    markers_vars <- reactiveValues()
  
    updateSelectizeInput(getDefaultReactiveDomain(), 
                         "cell_type", 
                         choices = cell_type_names)
    
    observeEvent(input$cell_type, {
      req(input$cell_type)
      markers_vars$marker_table <- readRDS(paste0("input/markers/cell_types/", input$cell_type, ".rds"))
    })
    
    # Cell types TABLE
    output$markers <- DT::renderDataTable({
      # req(input$cell_type)
      # readRDS(paste0("input/markers/cell_types/", input$cell_type, ".rds"))
      markers_vars$marker_table
    },
    selection = "single",
    server = TRUE)
    
    observeEvent(input$markers_rows_selected, {

      
      markers_vars$selected_gene <- markers_vars$marker_table[input$markers_rows_selected,]$gene
      markers_vars$counts <- readRDS(paste0("/active/paper_shiny/input/markers/counts_per_gene/", markers_vars$selected_gene, ".rds")) %>%
        cbind(metadata_all_cells[,-1])
      
    })
    
    output$spatial_plot <- renderPlot({
      req(markers_vars$counts)
      rbind(
        slice_sample(markers_vars$counts[markers_vars$counts$count == 0,], prop = 0.1),
        markers_vars$counts[markers_vars$counts$count > 0,]
      ) %>%
        # filter(count > 0) %>%
        arrange(count) %>%
        # mutate(alpha = ifelse(count == 0, 0.01, 1)) %>%
        
        # slice_sample(prop = 0.1) %>%
        ggplot(aes(x = x,
                   y = y,
                   colour = count,
                   alpha = count,
                   size = count)) +
        geom_point() +
        # geom_point(data = counts[counts$cell_type_publish == selected_cell_type,],
        #            colour = "orange") +
        facet_wrap(vars(mouse_id_presentation_no_genotype),
                   scales = "free") +
        theme_cowplot() +
        panel_border() +
        scale_y_reverse() +
        theme(
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.title = element_blank(),
          # legend.position = "none"
        ) +
        scale_colour_viridis_c(option = "A") +
        # scale_color_gradientn(colours = c("white", "white", "red"),
        #                       breaks = c(0, 3, 12),
        #                       limits = c(0, 12),
        #                       oob = scales::squish) +
        scale_size(range = c(0.01, 1.5)) +
        scale_alpha(range = c(0, 1)) +
        labs(colour = "Count", 
             size = "Count", 
             alpha = "Count")
    })
    
    

    output$debug <- renderPrint(markers_vars$counts)
    
  })
}
