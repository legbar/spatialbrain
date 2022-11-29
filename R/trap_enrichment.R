trap_enrichment_UI <- function(id) {
  ns <- NS(id)
  tabPanel("Dopaminergic TRAP Enrichment",
           titlePanel(h1("Dopaminergic TRAP Enrichment", align = 'center')),
           br(),
           fluidRow(
             column(8,
                    offset = 2,
                    plotOutput(ns("ma_plot")),
                    # verbatimTextOutput(ns("debug"))
                    )
           ),
           hr(),
           fluidRow(
             column(3,
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
                    h4("Select genes above to display information"),
                    DT::dataTableOutput(ns("enrichment_table")),
                    
                    ), 
             column(3, 
                    h4("Plot settings"), 
                    checkboxInput(ns("show_unchanged"),
                                  label = "Show Unchanged",
                                  value = T),
                    checkboxInput(ns("show_depleted"),
                                  label = "Show Depleted",
                                  value = T),
                    p(tags$em("To reset your selection, double click within the plot area")),
                    br(),
                    p(class = 'text-center', downloadButton(
                      ns('download_table'), 'Download Selected Data'
                    )),
                    # br(),
                    p(class = 'text-center', downloadButton(
                      ns('download_enriched_table'), 'Download All Enriched'
                    )),
                    # br(),
                    p(class = 'text-center', downloadButton(
                      ns('download_all_table'), 'Download All Data'
                    )),
                    style = 'border-left: 1px solid'
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
    
    # Genes TABLE
    output$enrichment_table <- DT::renderDataTable({
      # if (is.null(event_data("plotly_selected")$customdata)) {
      MB_FRACTION_META
      # } else {
      # MB_FRACTION_META[MB_FRACTION_META$Gene %in% event_data("plotly_selected")$customdata, ]
      # }
      
    },
    selection = "single",
    server = TRUE,
    rownames = FALSE)
    
    observeEvent(c(input$enrichment_table_rows_selected), {
      req(MB_FRACTION_META)
      if (is.null(input$enrichment_table_rows_selected)) {
      trap_enrichment_vars$selected_gene <- MB_FRACTION_META[1,] %>% pull(Gene)
      } else {
        trap_enrichment_vars$selected_gene <-
          MB_FRACTION_META[input$enrichment_table_rows_selected, ]$Gene
      }
      trap_enrichment_vars$highlight_markers <- trap_enrichment_vars$plot_data %>%
        filter(external_gene_name == trap_enrichment_vars$selected_gene)
    }, ignoreNULL = F)
    
    output$ma_plot <- renderPlot({
      req(trap_enrichment_vars$highlight_markers)
      trap_enrichment_vars$plot_data %>%
        ggplot(aes(x = baseMean_C1,
                   y = log2FoldChange,
                   colour = enrichment)) +
        theme_cowplot() +
        geom_point(
          size = 0.5,
                   alpha = 0.5
          ) +
        geom_point(data = trap_enrichment_vars$highlight_markers,
                   aes(fill = enrichment),
                   size = 5,
                   shape = 21,
                   colour = "black") +
        geom_label_repel(data = trap_enrichment_vars$highlight_markers,
                         aes(label = external_gene_name),
                         colour = "black",
                         arrow = arrow(length = unit(0.01, "npc")),
                         # box.padding = 1,
                         size = 5
                         # nudge_x = -1
        ) +
        scale_x_log10(breaks = c(1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6),
                      expand = expansion(mult = c(0, 0.05))) +
        scale_y_continuous(limits = c(NA, 6),
                           # breaks = c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5),
                           oob = scales::squish) +
        scale_color_manual(values = c("Depleted" = "#1F77B4FF",
                                      "Enriched" = "#2CA02CFF",
                                      "Unchanged" = "#C1C1C1")) +
        scale_fill_manual(values = c("Depleted" = "#1F77B4FF",
                                     "Enriched" = "#2CA02CFF",
                                     "Unchanged" = "#C1C1C1")) +
        labs(x = "Mean Counts",
             y = "Log2 Fold Change",
             colour = "Enrichment") +
        theme(legend.position = "top",
              # legend.text = element_text(size = 14),
              # legend.title = element_text(size = 14)
        ) +
        guides(color = guide_legend(override.aes = list(size = 3, alpha = 1) ) ,
               fill = FALSE)
      
    })
    
    # output$ma_plot <- renderPlotly({
    #   pal <- c("#C1C1C1",
    #            "#1F77B4FF",
    #            "#2CA02CFF")
    #   pal <- setNames(pal, c("Unchanged", "Depleted", "Enriched"))
    #   
    #   fig <- plot_ly(
    #     data = trap_enrichment_vars$plot_data,
    #     x = ~ baseMean_C1,
    #     y = ~ log2FoldChange,
    #     type = "scatter",
    #     mode = "markers",
    #     color = ~ enrichment,
    #     colors = pal,
    #     opacity = 0.5,
    #     text = ~ paste("Gene: ", external_gene_name),
    #     hoverinfo = "text",
    #     customdata = ~ external_gene_name
    #   ) %>%
    #     config(
    #       displayModeBar = T,
    #       toImageButtonOptions = list(
    #         filename = 'TRAP Enrichment Plot',
    #         width = 1366,
    #         height =  768
    #       )
    #     ) %>%
    #     layout(
    #       dragmode = "select",
    #       # title = list(text = "TRAP Enrichment Status", 
    #       #              xanchor = "left", 
    #       #              x = 0.01),
    #       xaxis = list(type = "log",
    #                    title = "Mean Counts", 
    #                    showline = T, 
    #                    linewidth = 2, 
    #                    linecolor = "black"),
    #       yaxis = list(range = list(-6, 6),
    #                    title = "Log<sub>2</sub> Fold Change", 
    #                    # zeroline = T, 
    #                    showline = T, 
    #                    linewidth = 2, 
    #                    linecolor = "black"),
    #       margin = list(t = 75), 
    #       legend=list(
    #         title=list(text='<b>Enrichment</b>'),
    #                   orientation = 'h', 
    #                   y = 6.5)
    #     )
    # })
    
    
    # download the filtered data
    output$download_table = downloadHandler(
      'TRAP Enrichment Information.csv',
      content = function(file) {
        write_csv(MB_FRACTION_META, file)
        # write_csv(MB_FRACTION_META[MB_FRACTION_META$Gene %in% event_data("plotly_selected")$customdata, ], file)
      }
    )
    
    output$download_enriched_table = downloadHandler(
      'TRAP Enrichment Information - All Enriched Genes.csv',
      content = function(file) {
        write_csv(MB_FRACTION_META[MB_FRACTION_META$`FDR-P` < 0.01 &
                                     MB_FRACTION_META$`Log2 Fold Enrichment` > 0, ], file)
      }
    )
    output$download_all_table = downloadHandler(
      'TRAP Enrichment Information - All Genes.csv',
      content = function(file) {
        write_csv(MB_FRACTION_META, file)
      }
    )
    
    # output$debug <- renderPrint(trap_enrichment_vars$plot_data)
    
  })
}
