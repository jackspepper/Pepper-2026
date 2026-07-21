library(shiny)
library(ggplot2)
library(dplyr)
library(data.table)
library(scattermore)
library(DT)  # For displaying tables
library(tidyr)  # For pivot_wider()
library(patchwork)

# Increase the maximum upload size to 300MB
options(shiny.maxRequestSize = 300*1024^2)

ui <- fluidPage(
  titlePanel("Scatter Plot with Dynamic Coloring"),
  checkboxInput("show_options", "Show Options", value = TRUE),
  checkboxInput("show_table", "Show Data Table", value = FALSE),  # New checkbox to show/hide table
  sidebarLayout(
    sidebarPanel(
      style = "display: flex; flex-direction: column; height: calc(100vh - 150px);",
      tags$div(
        style = "flex-grow: 1; overflow-y: scroll;",
        uiOutput("group_inputs_ui")
      )
    ),
    mainPanel(
      conditionalPanel(
        condition = "input.show_table == true",  # Condition to show/hide table
        DTOutput("group_table")  # Render the table
      ),
      div(
        style = "position: relative; width: 100%; height: calc(100vh - 150px);",
        plotOutput("scatterPlot", height = "100%", width = "100%")
      )
    )
  ),
  conditionalPanel(
    condition = "input.show_options == true",
    absolutePanel(
      id = "options-panel",
      right = 0, top = 60, width = 300, height = "calc(100vh - 150px)",
      style = "background-color: #f8f9fa; border-left: 1px solid #ddd; padding: 10px; overflow-y: auto;",
      div(
        style = "flex: 1 1 auto; min-width: 100px;",
        fileInput("file", "CSV File", accept = ".csv", width = "100%"),
        fileInput("annotations_file", "Annotations CSV", accept = ".csv", width = "100%"),
        actionButton("update_button", "Update Plot"),
        checkboxInput("update_plot", "Update Plot Automatically", value = FALSE)
      ),
      div(
        style = "flex: 1 1 auto; min-width: 200px;",
        uiOutput("xcol_ui"),
        uiOutput("ycol_ui"),
        uiOutput("clustercol_ui"),
        uiOutput("tissuecol_ui"),
        uiOutput("tissue_checkboxes_ui")
      ),
      div(
        style = "flex: 1 1 auto; min-width: 200px;",
        downloadButton("save_plot", "Save Plot"),
        downloadButton("download_annotations", "Download Annotations"),
        uiOutput("highlight_cluster_ui"),
        numericInput("point_size", "Point Size for Highlighted Cluster", value = 3, min = 0.01)
      )
    )
  )
)

server <- function(input, output, session) {
  # Reactive for loading the main data
  data <- reactive({
    req(input$file)
    fread(input$file$datapath)
  })
  
  annotations <- reactiveVal(NULL)
  
  observeEvent(input$annotations_file, {
    req(input$annotations_file)
    annotations(fread(input$annotations_file$datapath))
  })
  
  output$xcol_ui <- renderUI({
    req(data())
    selectInput("xcol", "X Column", choices = names(data()), selected = if ("umap_x" %in% names(data())) "umap_x" else NULL)
  })
  
  output$ycol_ui <- renderUI({
    req(data())
    selectInput("ycol", "Y Column", choices = names(data()), selected = if ("umap_y" %in% names(data())) "umap_y" else NULL)
  })
  
  output$clustercol_ui <- renderUI({
    req(data())
    selectInput("clustercol", "Cluster Column", choices = names(data()), selected = if ("PhenoGraph" %in% names(data())) "PhenoGraph" else NULL)
  })
  
  output$highlight_cluster_ui <- renderUI({
    req(data())
    cluster_levels <- sort(unique(data()[[input$clustercol]]))
    selectInput("highlight_cluster", "Highlight Cluster", choices = c("None" = "", cluster_levels), selected = "")
  })
  
  output$tissuecol_ui <- renderUI({
    req(data())
    selectInput("tissuecol", "Tissue Column", choices = c("None", names(data())), selected = if ("Tissue" %in% names(data())) "Tissue" else "None")
  })
  
  output$tissue_checkboxes_ui <- renderUI({
    req(input$tissuecol)
    if (input$tissuecol != "None") {
      tissue_levels <- sort(unique(data()[[input$tissuecol]]))
      checkboxGroupInput("selected_tissues", "Select Tissues", choices = tissue_levels, selected = tissue_levels)
    }
  })
  
  output$group_inputs_ui <- renderUI({
    req(input$clustercol)
    cluster_levels <- sort(unique(data()[[input$clustercol]]))
    ann <- annotations()
    lapply(cluster_levels, function(level) {
      group_value <- if (!is.null(ann) && level %in% ann$PhenoGraph) ann$Annotation[ann$PhenoGraph == level] else "TBC"
      notes_value <- if (!is.null(ann) && level %in% ann$PhenoGraph) ann$Notes[ann$PhenoGraph == level] else ""
      tagList(
        textInput(paste0("group_", level), paste("Group for Cluster", level), value = group_value),
        textInput(paste0("notes_", level), paste("Notes for Cluster", level), value = notes_value)
      )
    })
  })
  
  plot_data <- reactive({
    req(input$xcol, input$ycol, input$clustercol)
    plot_data <- data() %>%
      dplyr::select(!!sym(input$xcol), !!sym(input$ycol), !!sym(input$clustercol), !!sym(input$tissuecol), Timepoint) %>%
      mutate(across(all_of(input$clustercol), as.factor),
             Timempoint = factor(Timepoint, levels = c(2, 6, 24, 48, 72, 96, 144)))
    
    if (input$tissuecol != "None" && !is.null(input$selected_tissues)) {
      plot_data <- plot_data %>% filter(.data[[input$tissuecol]] %in% input$selected_tissues)
    }
    
    glimpse(plot_data)
    plot_data
  })
  
  plot_reactive <- reactiveVal(NULL)
  
  observeEvent(input$update_button, {
    plot_reactive({
      cluster_levels <- sort(unique(plot_data()[[input$clustercol]]))
      groups <- sapply(cluster_levels, function(level) input[[paste0("group_", level)]])
      names(groups) <- cluster_levels
      
      plot_data_with_group <- plot_data() %>%
        mutate(Group = factor(groups[as.character(plot_data()[[input$clustercol]])]))
      
      glimpse(plot_data_with_group)
      
      p <- ggplot(plot_data_with_group, aes(x = .data[[input$xcol]], y = .data[[input$ycol]], color = Group)) +
        geom_scattermore(data = plot_data_with_group, pointsize = 1, show.legend = TRUE)
      
      if (input$highlight_cluster != "") {
        p <- p + geom_scattermore(data = plot_data_with_group %>% filter(.data[[input$clustercol]] == input$highlight_cluster), 
                                  aes(x = .data[[input$xcol]], y = .data[[input$ycol]], color = Group), 
                                  pointsize = input$point_size, show.legend = TRUE)
      }
      
      line_plot <- NULL
      print(input$highlight_cluster)
      if (input$highlight_cluster != "") {
        highlighted_data <- plot_data_with_group %>% filter(.data[[input$clustercol]] == input$highlight_cluster)
        count_data <- highlighted_data %>%
          group_by(Timepoint) %>%
          summarise(Count = n())
        
        line_plot <- ggplot(count_data, aes(x = Timepoint, y = Count)) +
          geom_line() +
          geom_point() +
          labs(title = "Count of Highlighted Population Over Time",
               x = "Timepoint",
               y = "Count") +
          theme_minimal()
      }
      
      if (!is.null(line_plot)) {
        print(line_plot)
        combined_plot <- p / line_plot
      } else {
        combined_plot <- p
      }
      
      combined_plot + scale_color_discrete() + theme_minimal() + coord_fixed(ratio = 1)
    })
  })
  
  
  output$scatterPlot <- renderPlot({
    if (input$update_plot) {
      cluster_levels <- sort(unique(plot_data()[[input$clustercol]]))
      groups <- sapply(cluster_levels, function(level) input[[paste0("group_", level)]])
      names(groups) <- cluster_levels
      
      plot_data_with_group <- plot_data() %>%
        mutate(Group = factor(groups[as.character(plot_data()[[input$clustercol]])]))
      
      p <- ggplot(plot_data_with_group, aes(x = .data[[input$xcol]], y = .data[[input$ycol]], color = Group)) +
        geom_scattermore(data = plot_data_with_group, pointsize = 1, show.legend = TRUE)
      
      if (input$highlight_cluster != "") {
        p <- p + geom_scattermore(data = plot_data_with_group %>% filter(.data[[input$clustercol]] == input$highlight_cluster), 
                                  aes(x = .data[[input$xcol]], y = .data[[input$ycol]], color = Group), 
                                  pointsize = input$point_size, show.legend = TRUE)
      }
      
      p + scale_color_discrete() + theme_minimal() + coord_fixed(ratio = 1)
    } else {
      req(plot_reactive())
      plot_reactive()
    }
  }, height = function() {
    session$clientData$output_scatterPlot_width
  })

  
  
  output$save_plot <- downloadHandler(
    filename = function() {
      paste("scatter_plot", Sys.Date(), ".jpeg", sep = "")
    },
    content = function(file) {
      ggsave(file, plot = last_plot(), device = "jpeg", width = 10, height = 10, units = "in", dpi = 300)
    }
  )
  
  output$download_annotations <- downloadHandler(
    filename = function() {
      paste("annotations", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      cluster_levels <- sort(unique(data()[[input$clustercol]]))
      groups <- sapply(cluster_levels, function(level) input[[paste0("group_", level)]])
      notes <- sapply(cluster_levels, function(level) input[[paste0("notes_", level)]])
      names(groups) <- cluster_levels
      names(notes) <- cluster_levels
      
      annotations <- data.frame(
        PhenoGraph = cluster_levels,
        Annotation = groups,
        Notes = notes
      )
      
      write.csv(annotations, file, row.names = FALSE)
    }
  )
  
  output$group_table <- renderDT({
    req(plot_data())
    
    # Create a mapping of cluster levels to group names
    cluster_levels <- sort(unique(plot_data()[[input$clustercol]]))
    groups <- sapply(cluster_levels, function(level) input[[paste0("group_", level)]])
    names(groups) <- cluster_levels
    
    # Use groups in the table data
    table_data <- plot_data() %>%
      mutate(Group = groups[as.character(.data[[input$clustercol]])]) %>%
      group_by(Group, Cluster = .data[[input$clustercol]], Tissue = .data[[input$tissuecol]]) %>%
      summarise(Count = n(), .groups = 'drop') %>%
      group_by(Group, Tissue) %>%
      summarise(Count = sum(Count), Clusters = paste(unique(Cluster), collapse = ", "), .groups = 'drop') %>%
      pivot_wider(names_from = Tissue, values_from = Count, values_fill = list(Count = 0)) %>%
      dplyr::arrange(Group)  # Optional: Sort by Group
    
    datatable(table_data, options = list(pageLength = 5))
  })
  
}

shinyApp(ui, server)