ageing_TRAP_UI <- function(id) {
  ns <- NS(id)
  tabPanel(
    "Ageing in Dopaminergic Neurons",
    titlePanel(h1("Dopaminergic Ageing", align = 'center')),
    br(),
    fluidRow(column(
      width = 6,
      offset = 3,
      align = "center",
      plotOutput(ns("plot"),
                 # height = "300px"
                 # , width = "800px", height = "533px")
      ))),
    hr(),
    # fluidRow(
    #   column(12,
    #          plotOutput(ns("marker_violin_plot"))
    # )
    # ),
    fluidRow(
      column(
        3,
        h4(tags$i("Testing for age-related differential expression was performed using DESeq2 on gene-level counts derived from TRAP samples.")),
        br(),
        h4(helpText("Definitions")),
        hr(),
        p(tags$b("TOTAL: "), "Bulk RNA from ventral midbrain"),
        p(tags$b("TRAP: "), "RNA from DAT-TRAP"),
        p(tags$b("YOUNG: "), "Mice aged 3-6 months"),
        p(tags$b("OLD: "), "Mice aged 18-22 months"),
        p(tags$b("LFC: "), "The log2 fold-change in abundance in aged samples, compared to young"),
        p(tags$b("FDR-P: "), "The P value, adjusted for multiple comparisons (B&H)"),
        style = 'border-right: 1px solid'
      ),
      column(
        6,
        h4("Select a gene to view its expression with age"),
        DT::dataTableOutput(ns("table"))
      ),
      column(
        3,
        h4(helpText("Download Data")),
        hr(),
        p(class = 'text-center', downloadButton(
          ns('download_table'), 'Download Ageing Results Table'
        )),
        p(class = 'text-center', downloadButton(
          ns('download_plot'), 'Download Ageing Plot'
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

ageing_TRAP_SERVER <- function(id) {
  moduleServer(id, function(input, output, session) {
    # namespace ----
    ns <- session$ns
    
    ageing_TRAP_vars <- reactiveValues()
    
    ageing_TRAP_vars$table <-
      readRDS("input/startup/MB_AGE_META.rds")
    
    # observeEvent(input$table_rows_selected, {
    #   ageing_TRAP_vars$g <-
    #     ageing_TRAP_vars$table[input$table_rows_selected, ]$Gene
    # })
    
    # Ageing table
    output$table <-
      DT::renderDataTable(
        ageing_TRAP_vars$table,
        selection = "single",
        rownames = FALSE,
        server = T
        )
    
    observeEvent(input$table_rows_selected, {
      req(ageing_TRAP_vars$table)
      if (is.null(input$table_rows_selected)) {
        ageing_TRAP_vars$selected_gene <-
          ageing_TRAP_vars$table[1, ] %>% pull(Gene)
      } else {
        ageing_TRAP_vars$selected_gene <-
          ageing_TRAP_vars$table[input$table_rows_selected,]$Gene
      }
      ageing_TRAP_vars$counts <-
        readRDS(paste0("input/ageing/",
                       ageing_TRAP_vars$selected_gene,
                       ".rds")) %>%
        mutate(Age = factor(Age, levels = c("YOUNG", "OLD")))
    }, ignoreNULL = F)
    
    plot_func <- function() {
      ggplot(data = ageing_TRAP_vars$counts,
             aes(x = Age,
                 y = Count, 
                 fill = Age)) +
        facet_wrap(vars(Fraction), scales = "free_y") +
        geom_boxplot() +
        geom_jitter(shape = 21) +
        # geom_point(shape = 21) +
        theme_cowplot() +
        panel_border() +
        scale_fill_d3() +
        # scale_y_continuous(limits = c(0, 1)) +
        theme(legend.position = "none") +
        labs(
          x = "Age", 
             title = ageing_TRAP_vars$selected_gene)
    }
    #
    # # observeEvent(input$plot_size, {
    # #   markers_vars$plot_width <- paste0(input$plot_size, "px")
    # #   markers_vars$plot_height <- paste0((input$plot_size / 1.5), "px")
    # # })
    #
    # # output$spatial_plot_sized <- renderUI({
    # #   req(markers_vars$plot_height)
    # #   plotOutput(ns("spatial_plot")
    # #              # width = markers_vars$plot_width,
    # #              # height = markers_vars$plot_height,
    # #              # width = paste0(600, "px"),
    # #              # height = paste0(600/1.5, "px")
    # #              )
    # # })
    #
    output$plot <- renderPlot({
      req(ageing_TRAP_vars$counts)
      plot_func()
    })
    
    # plot_marker_func <- function(){
    #   readRDS(paste0("input/markers/counts_per_gene/", markers_vars$selected_gene, ".rds")) %>%
    #     full_join(metadata_all_cells) %>%
    #     ggplot(aes(x = cell_type_public,
    #                y = count)) +
    #     geom_violin(scale = "width") +
    #     theme_cowplot() +
    #     coord_flip() +
    #     labs(y = "Count", title = markers_vars$selected_gene) +
    #     theme(axis.title.y = element_blank())
    # }
    #
    # output$marker_violin_plot <- renderPlot({
    #   req(markers_vars$counts)
    #   plot_marker_func()
    # })
    #
    # download the data
    output$download_table = downloadHandler(
      paste0('Ageing results - TRAP.csv'),
      content = function(file) {
        write_csv(ageing_TRAP_vars$table, file)
      }
    )

    output$download_plot <- downloadHandler(
      filename = function() {
        paste0("Ageing Plot - ", ageing_TRAP_vars$selected_gene, ".png")
      },
      content = function(file) {
        ggsave(file, plot_func(),
               width = 800,
               height = 800/1.5,
               # width = input$plot_size,
               # height = (input$plot_size)/1.5,
               units = "px",
               dpi = 72,
               bg = "white")
      }
    )

    # output$download_marker_plot <- downloadHandler(
    #   filename = function() {
    #     paste0("Marker Plot - ", markers_vars$selected_gene, ".png")
    #   },
    #   content = function(file) {
    #     ggsave(file, plot_marker_func(),
    #            width = 800,
    #            height = 800/1.5,
    #            # width = input$plot_size,
    #            # height = (input$plot_size)/1.5,
    #            units = "px",
    #            dpi = 72,
    #            bg = "white")
    #   }
    # )
    
    # output$debug <- renderPrint(splicing_vars$counts)
    
    
  })
}
