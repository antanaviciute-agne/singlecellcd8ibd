library(shiny)
library(DT)

ExpressionInput <- function(id) {
  # Create a namespace function using the provided id
  ns <- NS(id)
  tagList(
    titlePanel("Single Cell Gene Expression"),
    sidebarPanel(
      
      selectInput(ns("use_data"), label = "Select Compartment", choices = dts, selected = dts[1]),
      selectInput(ns("colour_by"), label = "Colour Cells By", choices = colnames(meta_data), selected = colnames(meta_data)[1]),
      selectInput(ns("facet_by1"), label = "1st Facet Cells By", choices = c(".", colnames(meta_data)), selected = "."),
      selectInput(ns("facet_by2"), label = "2nd Facet Cells By", choices = c(".", colnames(meta_data)), selected = "."),
      selectInput(ns("cell_size"), label = "Cell Size", choices = c(0.5,1:10), selected = "1"),
      selectInput(ns("filter_cells"), label="Filter Cells By", choices = c(colnames(meta_data)), selected = colnames(meta_data)[1]),
      uiOutput(ns("filter_list"))
    ),
    
    mainPanel(
      
      verticalLayout(
        tabsetPanel(type = "tabs",
                    tabPanel("UMAP", plotOutput(ns('umap'))),
                    tabPanel("Boxplot", plotOutput(ns('boxplot')))),
        h6(""),
        tabsetPanel(type ="tabs",
                    tabPanel("Gene Selection Table", dataTableOutput( ns('table'))),
                    tabPanel("CITE-Seq Protein Table", dataTableOutput( ns('protein.table'))))
        
      )
    )
  )
}

#Server
Expression <- function(input, output, session) {
  ns <- session$ns
  
  
  last <- NA
  sr <- NA
  covar <- NA
  observeEvent(input$table_rows_selected, {last <<- "gene"
  })
  
  observeEvent(input$protein.table_rows_selected, {last <<- "protein"
  })
  
  observeEvent(input$colour_by, {last <<- "covar"
  
  
  selectRows( proxy = dataTableProxy('table'), selected = NULL)
  selectRows( proxy = dataTableProxy('protein.table'), selected = NULL)
  
  })
  
  observeEvent(input$use_data, {
    
    last <<- "covar"
    
  })
  
  
  observeEvent(input$use_data, {
    
    sr <<- clusterings[[input$use_data]]
    covar <<- meta_data[rownames(sr), ]
    updateSelectInput(session, "colour_by", selected=colnames(meta_data)[2])
    updateSelectInput(session, "colour_by", selected=colnames(meta_data)[1])
    last <<- "covar"
    
  })
  
  display_cell_ids <- rownames(meta_data)
  
  output$filter_list <- renderUI({
    
    ns <- session$ns
    
    opts <- unique(meta_data[, input$filter_cells])
    
    if ( class(opts) %in% c("factor", "character")){
      
      display_cell_ids <<- rownames(meta_data)
      
      checkboxGroupInput(ns("checkbox1"), "", choices = opts,
                         selected = opts)
      
    }
    else{
      display_cell_ids <<- rownames(meta_data)
    }
    
  })

  observeEvent(input$checkbox1, {
    
    if ( class(meta_data[, input$filter_cells]) %in% c("factor", "character")){
      display_cell_ids <<- rownames(meta_data[which(meta_data[, input$filter_cells] %in% input$checkbox1), ])
      print(length(display_cell_ids))
      
      tmp <- input$colour_by
      
      ##trigger event to refresh, switch back to original selection, bit hacky but works
      updateSelectInput(session, "colour_by", selected=colnames(meta_data)[2])
      updateSelectInput(session, "colour_by", selected=colnames(meta_data)[1])
      updateSelectInput(session, "colour_by", selected=tmp)
      
    }
    else{
      display_cell_ids <<- rownames(meta_data)
    }
  })
  
  output$table <- renderDataTable({
    
    datatable(as.data.frame(gene.table), filter = "top", 
              selection = list(mode = "single", selected = NA))
  })
  
  output$protein.table <- renderDataTable({
    
    datatable(as.data.frame(protein.table), filter = "top", 
              selection = list(mode = "single", selected = NA))
  })
  
  output$umap <- renderPlot({
    
    coords <- data.frame(sr, covar)
    coords <- coords[display_cell_ids[display_cell_ids %in% rownames(coords)], ]
    
    p <- ggplot(coords, aes_string("UMAP_1", "UMAP_2", color=input$colour_by)) + geom_point(size=as.numeric(input$cell_size)) + labs(x="UMAP 1", y="UMAP 2")
    
    gene <- input$table_rows_selected
    
    protein <- input$protein.table_rows_selected
    
     
    
    if(!is.null(gene) & last == "gene"){
       
      
      exp <-expression.data[gene, display_cell_ids[display_cell_ids %in% rownames(coords)] ]
      coords$gene <- exp
      
      p <- ggplot(coords, aes_string("UMAP_1", "UMAP_2", color="gene")) + geom_point(size=as.numeric(input$cell_size)
      ) + labs(x="UMAP 1", y="UMAP 2", color="Expression") + scale_color_viridis_c(option = "D", direction=1)
    }
    
    if(!is.null(protein) & last == "protein"){

      
      
      exp <-protein.data[protein, display_cell_ids[display_cell_ids %in% rownames(coords)] ]
      coords$gene <- exp
      
      p <- ggplot(coords, aes_string("UMAP_1", "UMAP_2", color="gene")) + geom_point(size=as.numeric(input$cell_size)
      ) + labs(x="UMAP 1", y="UMAP 2", color="Expression") + scale_color_viridis_c(option = "D", direction=1)
    }
    
    facets <- paste(input$facet_by1, '~', input$facet_by2)
    
    
    if (facets != '. ~ .')
      p <- p + facet_grid(facets)
    
    print(p + theme_bw())
    
  })
  
  
  output$boxplot <- renderPlot({
    
    ind <- which(rownames(covar) %in% display_cell_ids[display_cell_ids %in% rownames(covar)])
    
    
    gene <- input$table_rows_selected
    protein <- input$protein.table_rows_selected
    
    print(last)
    
    if(!is.null(gene) & last == "gene"){
      
      exp <- as.numeric(expression.data[gene, rownames(covar)])
      
      
      covar$gene <- exp
      
    
      if(class(covar[, input$colour_by]) %in% c("numeric")){ ## numeric covariate
        
        
        p <-ggplot(covar[ind, ], aes_string(input$colour_by, "gene", colour=input$colour_by)) + geom_point() + scale_color_viridis_c(option = "D", direction=1)
      }
      else{##factor covariate
        
        p <- ggplot(covar[ind, ], aes_string(x=input$colour_by, y="gene", fill=input$colour_by)) + geom_boxplot(aes_string(color="gene"), outlier.color = NA) + geom_jitter(height=0, width=0.1, alpha=0.2) + coord_flip()
        
        p <- p + aes(color=exp[ind])+ scale_color_viridis_c(option = "D", direction=1)
      }
      
      if(class(covar[, input$colour_by]) == "numeric"){ ## numeric covariate, use nice viridis scale!
        
        p <- p + scale_color_viridis_c(option="D", direction=1)
      }
      
      facets <- paste(input$facet_by1, '~', input$facet_by2)
      
      
      if (facets != '. ~ .')
        p <- p + facet_grid(facets)
      
      p + theme_bw() + labs(y="Expression", x="", color="Exp.") + guides(fill=FALSE)
      
      
    }
    ##########################################
    
    else if(!is.null(protein) & last == "protein"){
      
      exp <- as.numeric(protein.data[protein, rownames(covar)])
      
      
      covar$gene <- exp
      
      
      if(class(covar[, input$colour_by]) %in% c("numeric")){ ## numeric covariate
        
        
        p <-ggplot(covar[ind, ], aes_string(input$colour_by, "gene", colour=input$colour_by)) + geom_point() + scale_color_viridis_c(option = "D", direction=1)
      }
      else{##factor covariate
        
        p <- ggplot(covar[ind, ], aes_string(x=input$colour_by, y="gene", fill=input$colour_by)) + geom_boxplot(aes_string(color="gene"), outlier.color = NA) + geom_jitter(height=0, width=0.1, alpha=0.2) + coord_flip()
        
        p <- p + aes(color=exp[ind])+ scale_color_viridis_c(option = "D", direction=1)
      }
      
      if(class(covar[, input$colour_by]) == "numeric"){ ## numeric covariate, use nice viridis scale!
        
        p <- p + scale_color_viridis_c(option="D", direction=1)
      }
      
      facets <- paste(input$facet_by1, '~', input$facet_by2)
      
      
      if (facets != '. ~ .')
        p <- p + facet_grid(facets)
      
      p + theme_bw() + labs(y="Expression", x="", color="Exp.") + guides(fill=FALSE)
      
      
    }
    
  })
  
}