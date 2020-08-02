options(repos = BiocManager::repositories())
library(shiny)
library(circlize)
library(ggplot2)
library(ggrepel)
library(shinythemes)
library(Seurat)
library(DT)
library(igraph)
library(dplyr)
library(ggraph)
library(cowplot)

options(stringsAsFactors = FALSE)
Logged = TRUE; ## switch this  and add a username + MD5 password hash below to set a very rudimentary log in
PASSWORD <- data.frame(Username = c( ), Password = c())

shinyServer(function(input, output, session) {
  
  source("www/Login.R",  local = TRUE)
  
  observe({
    if (USER$Logged == TRUE) {
      
      output$obs <- renderUI({
        
        navbarPage("Human Gut CD8+ T-Cell Atlas in UC", theme=shinythemes::shinytheme(theme="simplex"),
                   
                   
                   tabPanel("Single Cell Gene & Protein Expression",
                            ExpressionInput("expression.module"))
                   
                   
        )

      })
      
      
      callModule(Expression, "expression.module")

    }
  
  })
  
})
