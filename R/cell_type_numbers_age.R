cell_type_numbers_age_UI <- function(id) {
  ns <- NS(id)
  tabPanel("Cell Type Numbers with Age",
           titlePanel(h1("Cell Type Numbers with Age", align = 'center')),
           br(),
           fluidRow(
             column(width = 5,
                    offset = 1,
                    align = "center",
                    # uiOutput(ns("spatial_plot_sized"))
                    plotOutput(ns("volcano_plot"), 
                               # height = "500px"
                               # , width = "800px", height = "533px"
                               )
             ), 
             column(width = 5, 
                    # offset = 6, 
                    align = "center",
                    plotOutput(ns("actual_counts_plot"), 
                               # height = "500px"
                               )
                    )
           ),
           fluidRow(
             column(width = 10,
                    offset = 1,
                    align = "center",
                    plotOutput(ns("xy_plot"),
                               heigh = "250px"
                    )
             )
           ),
           hr(),
           fluidRow(
             column(3,
                    h4(tags$i("Testing for age-related changes in cell type abundance was performed using MASC in Stereo-Seq samples")),
                    br(),
                    h4(helpText("Definitions")),
                    hr(),
                    p(tags$b("YOUNG: "), "Mice aged 3-6 months"),
                    p(tags$b("OLD: "), "Mice aged 18-22 months"),
                    p(tags$b("FDR-P: "), "The P value, adjusted for multiple comparisons (B&H)"),
                    p(tags$b("Estimate: "), "Odds ratio for abundance change in aged brains. Positive values indicate an increase in abundance"),
                    
                    style = 'border-right: 1px solid'
             ), 
             column(6, 
                    h4("Select a cell type to view its abundance across age groups"),
                    DT::dataTableOutput(ns("masc_table")),
                    # selectizeInput(ns("cell_type"),
                    #                label = "Select cell type",
                    #                choices = NULL,
                    #                width = "100%")
             ), 
             column(3, 
                    # h3(helpText("Select a cell type to view markers...")),
                    # hr(),
                    # selectizeInput(ns("cell_type"), 
                    #                label = "Select cell type", 
                    #                choices = NULL, 
                    #                width = "100%"),
                    br(),
                    p(class = 'text-center', downloadButton(
                      ns('download_table'), 'Download Overall Results'
                    )),
                    p(class = 'text-center', downloadButton(
                      ns('download_numbers'), 'Download Cell Type Numbers per Brain'
                    )),
                    # p(class = 'text-center', downloadButton(
                    #   ns('download_marker_plot'), 'Download Marker Plot'
                    # )),
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

# server ----
cell_type_numbers_age_SERVER <- function(id, cell_type_names) {
  moduleServer(id, function(input, output, session) {
    # namespace ----
    ns <- session$ns
    
    cell_type_numbers_age_vars <- reactiveValues()
  
    cell_type_numbers_age_vars$masc <- readRDS("input/ageing_cell_type_numbers/data.rds") %>%
      arrange(fdr) %>%
      mutate(across(where(is.numeric), ~ signif(.x, 3)))
    
    output$masc_table <- DT::renderDataTable({
      cell_type_numbers_age_vars$masc %>%
        select(-cell_type, 
               "Cell Type" = cell_type_full, 
               "Estimate" = estimate, 
               "FDR" = fdr)
    },
    selection = "single",
    server = TRUE,
    rownames = FALSE)
    
    observeEvent(c(input$masc_table_rows_selected, input$cell_type), {
      req(cell_type_numbers_age_vars$masc)
      if (is.null(input$masc_table_rows_selected)) {
        cell_type_numbers_age_vars$selected_cell_type <- cell_type_numbers_age_vars$masc[1,] %>% pull(cell_type)
      } else {
        cell_type_numbers_age_vars$selected_cell_type <-
          cell_type_numbers_age_vars$masc[input$masc_table_rows_selected, ]$cell_type
      }
    }, ignoreNULL = F)
    
    # updateSelectizeInput(getDefaultReactiveDomain(), 
    #                      "cell_type", 
    #                      choices = cell_type_names, 
    #                      selected = "DA_SN")
    
    output$volcano_plot <- renderPlot({
      req(cell_type_numbers_age_vars$selected_cell_type)
      cell_type_numbers_age_vars$masc %>%
        ggplot(aes(x = estimate, 
                   y = -log10(fdr), 
                   label = ifelse(fdr < 0.005, 
                                  cell_type_full, 
                                  ""), 
                   color = ifelse(cell_type == cell_type_numbers_age_vars$selected_cell_type, 
                                  "red", "black"))) +
        geom_point() +
        geom_vline(xintercept = 0, linetype = "dotted") +
        geom_hline(yintercept = -log10(0.01), linetype = "dotted") +
        geom_label_repel() +
        scale_x_continuous(limits = c(-1.5, 1.5), 
                           oob = scales::squish) +
        labs(x = expression(Log[2] ~ OR), 
             y = expression(-Log[10] ~ FDR)) +
        scale_color_manual(values = c("black", "red")) +
        theme_cowplot() +
        theme(legend.position = "none")
    })
    
    actual_numbers <- readRDS("input/ageing_cell_type_numbers/actual_numbers.rds")
    
    output$actual_counts_plot <- renderPlot({
      req(cell_type_numbers_age_vars$selected_cell_type)
    actual_numbers %>%
      filter(cell_type_publish == cell_type_numbers_age_vars$selected_cell_type) %>%
      ggplot(aes(x = age, 
                 y = n)) +
      geom_point() +
      theme_cowplot() +
      labs(x = "Age", y = "Number of Cells", title = names(cell_type_names[cell_type_names == cell_type_numbers_age_vars$selected_cell_type]))
    })
    
    xy <- readRDS("input/ageing_cell_type_numbers/xy.rds")
    
    output$xy_plot <- renderPlot({
      req(xy)
      
      xy %>%
      mutate(
        interest = ifelse(
          cell_type_publish == cell_type_numbers_age_vars$selected_cell_type,
          cell_type_publish,
          "OTHER"
        ),
        size_interest = cell_type_publish == cell_type_numbers_age_vars$selected_cell_type
      ) %>%
        mutate(interest = factor(interest), 
               interest = relevel(interest, ref = cell_type_numbers_age_vars$selected_cell_type)) %>%
        # mutate(interest = ifelse(str_detect(interest, "OTHER"), "OTHER", names(cell_type_names[cell_type_names == cell_type_numbers_age_vars$selected_cell_type]))) %>%
        arrange(desc(interest)) %>%
        ggplot(aes(
          x = x,
          y = -y,
          colour = interest,
          size = size_interest,
          alpha = size_interest
        )) +
        geom_point() +
        scale_color_manual(values = c(pal_d3()(1), "grey")) +
        scale_size_manual(values = c(0.02, 1), guide = 'none') +
        scale_alpha_manual(values = c(0.4, 1), guide = 'none') +
        facet_wrap(vars(mouse_id_label),
                   nrow = 1,
                   scales = "free") +
        theme_minimal() +
        theme(
          rect = element_rect(fill = "transparent"),
          # panel.background = element_rect(fill = 'white', colour = 'white'),
          panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
          strip.text = element_text(face = "bold", size = 14), 
          legend.position = "none"
        ) +
        labs(title = names(cell_type_names[cell_type_names == cell_type_numbers_age_vars$selected_cell_type]),
             color = "")
    })
    
    # download the data
    output$download_table = downloadHandler(
      filename = function() {
        paste0('Cell type numbers with age - MASC Results.csv')
      },
      content = function(file) {
        write_csv({
          cell_type_numbers_age_vars$masc %>%
            select(-cell_type, 
                   "Cell Type" = cell_type_full, 
                   "Estimate" = estimate, 
                   "FDR" = fdr)
        }, file)
      }
    )
    
    output$download_numbers = downloadHandler(
      filename = function() {
        paste0('Cell type numbers per brain.csv')
      },
      content = function(file) {
        write_csv({
          readRDS("input/ageing_cell_type_numbers/actual_numbers_with_ids.rds")
        }, file)
      }
    )

    # output$debug <- renderPrint(cell_type_numbers_age_vars$selected_cell_type)
    
  })
}
