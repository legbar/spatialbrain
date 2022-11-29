spatial_markers_UI <- function(id) {
  ns <- NS(id)
  tabPanel("Spatial Markers",
           titlePanel(h1("Spatial Markers", align = 'center')),
           br(),
           fluidRow(
             column(width = 5,
                    offset = 1,
                    align = "center",
                    # uiOutput(ns("spatial_plot_sized"))
                    plotOutput(ns("spatial_plot"), 
                               height = "600px"
                               # , width = "800px", height = "533px"
                               )
             ), 
             column(width = 5, 
                    # offset = 6, 
                    align = "center",
                    plotOutput(ns("marker_violin_plot"), 
                               height = "600px"
                               )
                    )
           ),
           hr(),
           # fluidRow(
           #   column(12,
           #          plotOutput(ns("marker_violin_plot"))
           # )
           # ),
           fluidRow(
             column(3,
                    h3(helpText("Select a cell type to view markers...")),
                    hr(),
                    selectizeInput(ns("cell_type"), 
                                   label = "Select cell type", 
                                   choices = NULL, 
                                   width = "100%"),
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
                    h4("Select a gene to view its spatial expression"),
                    DT::dataTableOutput(ns("spatial_markers"))
             ), 
             column(3, 
                    h4("Definitions"),
                    strong("LFC: "), span("The log2 fold-change in abundance between the cell type of interest and all other cells"),
                    br(),
                    br(),
                    strong("FDR-P "), span("The P value, adjusted for multiple comparisons (B&H)"),
                    hr(),
                    h4("Download Data"),
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
           ))
  
}

# markers server ----
spatial_markers_SERVER <- function(id, metadata_all_cells, cell_type_names) {
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
                         choices = cell_type_names, 
                         selected = "DA_SN")
    
    observeEvent(input$cell_type, {
      req(input$cell_type)
      markers_vars$marker_table <- readRDS(paste0("input/markers/cell_types/", input$cell_type, ".rds")) %>%
        arrange(desc(LFC))
        # filter(lfc > 0) %>%
        # mutate(across(where(is.numeric), ~ signif(.x, 3)))
      markers_vars$selected_cell_type_public <- names(cell_type_names[cell_type_names == input$cell_type])
    })
    
    # Cell types TABLE
    output$spatial_markers <- DT::renderDataTable({
      # req(input$cell_type)
      # readRDS(paste0("input/markers/cell_types/", input$cell_type, ".rds"))
      markers_vars$marker_table
    },
    selection = "single",
    server = TRUE,
    rownames = FALSE)
    
    observeEvent(c(input$spatial_markers_rows_selected, input$cell_type), {
      req(markers_vars$marker_table)
      if (is.null(input$spatial_markers_rows_selected)) {
        markers_vars$selected_gene <- markers_vars$marker_table[1,] %>% pull(Gene)
      } else {
        markers_vars$selected_gene <-
          markers_vars$marker_table[input$spatial_markers_rows_selected, ]$Gene
      }
      markers_vars$counts <-
        readRDS(paste0(
          "input/markers/counts_per_gene/",
          markers_vars$selected_gene,
          ".rds"
        )) %>%
        cbind(metadata_all_cells[, -1])
    }, ignoreNULL = F)
    
    plot_spatial_func <- function(){
        # filter(count > 0) %>%
        # arrange(count) %>%
        # mutate(alpha = ifelse(count == 0, 0.01, 1)) %>%
        
        # slice_sample(prop = 0.1) %>%
        ggplot(data = slice_sample(markers_vars$counts[markers_vars$counts$count == 0,], prop = 0.1),
               aes(x = x,
                   y = y,
                   # colour = count
                   # alpha = count
                   # size = count
                   )) +
        geom_point(colour = "lightgrey", size = 0.1) +
        geom_point(data = {arrange(markers_vars$counts[markers_vars$counts$count > 0,], count) %>%
            mutate(count = ifelse(count > quantile(count, 0.99), 
                                  quantile(count, 0.99), 
                                  count)) %>%
            mutate(count = ifelse(count < quantile(count, 0.05), 
                                  0, 
                                  count))}, 
                   aes(colour = count, 
                       size = count)) +
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
        scale_color_gradient(low = "lightgrey", high = "red") +
        # scale_colour_viridis_c(option = "A") +
        # scale_color_gradientn(colours = c("white", "white", "red"),
        #                       breaks = c(0, 3, 12),
        #                       limits = c(0, 12),
        #                       oob = scales::squish) +
        scale_size(range = c(0.1, 2)) +
        # scale_alpha(range = c(0, 1)) +
        labs(colour = "Count", 
             size = "Count", 
             alpha = "Count", 
             title = markers_vars$selected_gene)
    }
    
    # observeEvent(input$plot_size, {
    #   markers_vars$plot_width <- paste0(input$plot_size, "px")
    #   markers_vars$plot_height <- paste0((input$plot_size / 1.5), "px")
    # })
    
    # output$spatial_plot_sized <- renderUI({
    #   req(markers_vars$plot_height)
    #   plotOutput(ns("spatial_plot")
    #              # width = markers_vars$plot_width,
    #              # height = markers_vars$plot_height,
    #              # width = paste0(600, "px"), 
    #              # height = paste0(600/1.5, "px")
    #              )
    # })
    
    output$spatial_plot <- renderPlot({
      req(markers_vars$counts)
      plot_spatial_func()
    })
    
    plot_marker_func <- function(){
      readRDS(paste0("input/markers/counts_per_gene/", markers_vars$selected_gene, ".rds")) %>% 
        full_join(metadata_all_cells) %>%
        ggplot(aes(x = cell_type_public, 
                   y = count)) +
        geom_violin(scale = "width") +
        theme_cowplot() +
        coord_flip() +
        labs(y = "Count", title = markers_vars$selected_gene) +
        theme(axis.title.y = element_blank())
    }
    
    output$marker_violin_plot <- renderPlot({
      req(markers_vars$counts)
      plot_marker_func()
    })
    
    # download the filtered data
    output$download_table = downloadHandler(
      filename = function() {
      paste0('Spatial Markers - ', markers_vars$selected_cell_type_public, '.csv')
      },
      content = function(file) {
        write_csv(markers_vars$marker_table, file)
      }
    )
    
    output$download_spatial_plot <- downloadHandler(
      filename = function() {
        paste0("Spatial Plot - ", markers_vars$selected_gene, ".png")
      },
      content = function(file) {
        ggsave(file, plot_spatial_func(),
               width = 800, 
               height = 800/1.5,
               # width = input$plot_size,
               # height = (input$plot_size)/1.5,
               units = "px",
               dpi = 72,
               bg = "white")
      }
    )
    
    output$download_marker_plot <- downloadHandler(
      filename = function() {
        paste0("Marker Plot - ", markers_vars$selected_gene, ".png")
      },
      content = function(file) {
        ggsave(file, plot_marker_func(),
               width = 800, 
               height = 800/1.5,
               # width = input$plot_size,
               # height = (input$plot_size)/1.5,
               units = "px",
               dpi = 72,
               bg = "white")
      }
    )

    # output$debug <- renderPrint(reactiveValuesToList(session$clientData))
    
  })
}
