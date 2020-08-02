shinyUI(bootstrapPage(
  tagList(
    tags$head(
      tags$link(rel="stylesheet", type="text/css",href="style.css"),
      tags$script(type="text/javascript", src = "md5.js"),
      tags$script(type="text/javascript", src = "passwdInputBinding.js")
    )
  ),
  
  ## Login module;
  div(class = "login",
      uiOutput("uiLogin"),
      textOutput("pass")
  ), 
  
  div(class = "span4", uiOutput("obs")),
  div(class = "span8", plotOutput("distPlot"))
  
))