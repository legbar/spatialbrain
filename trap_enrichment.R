trap_enrichment_UI <- function(id) {
  ns <- NS(id)
  tabPanel("Dopaminergic TRAP Enrichment",
           sidebarLayout(
             position = "left",
             sidebarPanel(
               width = 4,
               wellPanel(
                 h3(helpText("What am I looking at?"))
               ),
               wellPanel(
                 h3(helpText("Selected Gene Information")),
                 hr(),
                 # selectizeInput(ns("gene"),
                 #                label = "Select gene",
                 #                choices = NULL),
                 br(),
                 DT::dataTableOutput(ns("enrichment_table")),
                 br(),
                 p(class = 'text-center', downloadButton(
                   ns('download_table'), 'Download Selected Information'
                 )),
                 # br(),
                 p(class = 'text-center', downloadButton(
                   ns('download_enriched_table'), 'Download All Enriched'
                 )),
                 # br(),
                 p(class = 'text-center', downloadButton(
                   ns('download_all_table'), 'Download All Information'
                 ))
               )
             ),
             mainPanel(
               width = 12 - side_width,
               # wellPanel(plotOutput(ns("ma_plot"))),
               wellPanel(
                 h4(helpText(
                   "Select Genes to Display More Information"
                 )),
                 shinycssloaders::withSpinner(
                   plotlyOutput(ns("ma_plot")), 
                   type = 4
                 ),
                 hr(),
                 checkboxInput(ns("show_unchanged"),
                               label = "Show Unchanged", 
                               value = T),
                 checkboxInput(ns("show_depleted"),
                               label = "Show Depleted", 
                               value = T)
               ),
               # wellPanel(),
               wellPanel(verbatimTextOutput(ns("debug")))
             )
           ))
}

trap_enrichment_SERVER <- function(id) {
  moduleServer(id, function(input, output, session) {
    # namespace ----
    ns <- session$ns
    
    trap_enrichment_vars <- reactiveValues()
    
    MB_FRACTION_META <-
      readRDS("input/startup/MB_FRACTION_META.rds")
    plot_data <- readRDS("input/startup/MA_plot_data.rds") %>%
      mutate(enrichment = factor(enrichment,
                                 levels = c("Unchanged",
                                            "Depleted",
                                            "Enriched")))
    
    observeEvent(input$show_unchanged | input$show_depleted,
                 {
                   req(plot_data)
                   trap_enrichment_vars$plot_data <- plot_data
                   
                   if (input$show_unchanged == F) {
                     trap_enrichment_vars$plot_data <- trap_enrichment_vars$plot_data %>%
                       filter(enrichment != "Unchanged")
                   }
                   if (input$show_depleted == F) {
                     trap_enrichment_vars$plot_data <- trap_enrichment_vars$plot_data %>%
                       filter(enrichment != "Depleted")
                   }
                   
                   trap_enrichment_vars$plot_data <-
                     trap_enrichment_vars$plot_data %>%
                     arrange(enrichment)
                   
                 })
    
    
    # updateSelectizeInput(getDefaultReactiveDomain(),
    #                      "gene",
    #                      choices = MB_FRACTION_META$Gene)
    
    # observeEvent(input$gene, {
    #   trap_enrichment_vars$highlight_markers <- plot_data %>%
    #     # left_join(anno) %>%
    #     filter(external_gene_name == input$gene)
    # })
    
    output$ma_plot <- renderPlotly({
      pal <- c("#C1C1C1",
               "#1F77B4FF",
               "#2CA02CFF")
      pal <- setNames(pal, c("Unchanged", "Depleted", "Enriched"))
      
      fig <- plot_ly(
        data = trap_enrichment_vars$plot_data,
        x = ~ baseMean_C1,
        y = ~ log2FoldChange,
        type = "scatter",
        mode = "markers",
        color = ~ enrichment,
        colors = pal,
        opacity = 0.5,
        text = ~ paste("Gene: ", external_gene_name),
        hoverinfo = "text",
        customdata = ~ external_gene_name
      ) %>%
        config(
          displayModeBar = T,
          toImageButtonOptions = list(
            filename = 'TRAP Enrichment Plot',
            width = 1366,
            height =  768
          )
        ) %>%
        layout(
          dragmode = "select",
          title = "TRAP Enrichment Status",
          xaxis = list(type = "log",
                       title = "Mean Counts", 
                       showline = T, 
                       linewidth = 2, 
                       linecolor = "black"),
          yaxis = list(range = list(-6, 6),
                       title = "Log<sub>2</sub> Fold Change", 
                       # zeroline = T, 
                       showline = T, 
                       linewidth = 2, 
                       linecolor = "black"),
          margin = list(t = 75)
        )
    })
    # Cell types TABLE
    output$enrichment_table <- DT::renderDataTable({
      MB_FRACTION_META[MB_FRACTION_META$Gene %in% event_data("plotly_selected")$customdata, ]
    },
    # selection = "single",
    server = T)
    
    # download the filtered data
    output$download_table = downloadHandler(
      'TRAP Enrichment Information.csv',
      content = function(file) {
        write.csv(MB_FRACTION_META[MB_FRACTION_META$Gene %in% event_data("plotly_selected")$customdata, ], file)
      }
    )
    
    output$download_enriched_table = downloadHandler(
      'TRAP Enrichment Information - All Enriched Genes.csv',
      content = function(file) {
        write.csv(MB_FRACTION_META[MB_FRACTION_META$`FDR-P` < 0.01 &
                                     MB_FRACTION_META$`Log2 Fold Enrichment` > 0, ], file)
      }
    )
    output$download_all_table = downloadHandler(
      'TRAP Enrichment Information - All Genes.csv',
      content = function(file) {
        write.csv(MB_FRACTION_META, file)
      }
    )
    
    output$debug <- renderPrint(input$enrichment_table_rows_all)
    
  })
}
