library(shiny)
library(markdown)
library(ggplot2)
library(grid)

source("functions.R")
source("functions_plot.R")

# Define UI for CNV calculator ----
ui <- shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("CNV Power Calculator"),
  
  # Sidebar with controls to select a dataset and specify the number
  # of observations to view
  sidebarPanel(
    helpText("Enter experimental parameters to calculate power"),
    radioButtons('dis', label = 'Distribution*', c('Poisson' = 'pois','Negative Binomial' = 'nb')),
    numericInput("F", "Tumor Purity (range 0-1)*:", 1,min=0,max=1),
    numericInput("alpha", "Significance Level alpha (range 0-1)*:", 0.05,min=0,max=1),
    numericInput("L", "Length of CNV (in bases)*:", 10000,min=0),
    numericInput("N", "Ploidy of CNV*:", 3,min=0),
    numericInput("W", "Window length (in bases)*:", 100,min=0),
    numericInput("l", "Length of read (bases)*:", 100,min=0),
    numericInput("D", "Haploid Read Depth (integer value)*: ", 30,min=0),
    numericInput("phi","Overdispersion parameter mu*phi (optional):",min=0,max=1,value = '')
  ),
  
  # Show a summary of the dataset and an HTML table with the requested
  # number of observations
  mainPanel(
    tabsetPanel(type = "tabs",
                tabPanel("Input Parameters", verbatimTextOutput("toPrint")),
                tabPanel("Power", verbatimTextOutput("power")),
                tabPanel("Power Curve",verbatimTextOutput("powerCurve_text"),
                         selectInput('tovary',label="First parameter to vary along x-axis", choices=c('None','Read Depth','Read length','Window size','Length of CNV','Sample purity')),
                         numericInput("min1","Minimum value of parameter to be tested",1),
                         numericInput("max1","Maximum value of parameter to be tested",10),

                         selectInput('tovary2',label="Second Parameter to Vary (line color)", choices=c('None','Read Depth','Read length','Window size','Length of CNV','Sample purity')),
                         numericInput("val1","Value 1 to be tested for parameter 2",1),
                         numericInput("val2","Value 2 to be tested for parameter 2",2),
                         numericInput("val3","Value 3 to be tested for parameter 2",3),
                         plotOutput("powerCurve",width="75%"),
                         downloadButton("downloadData", "Download Raw Values")
                ),
                tabPanel("About",includeMarkdown("README.md"))
                
    )
#    print("Power is:"),
 #   textOutput('power'),
  #  tableOutput("view")
  )
))

server = shinyServer(function(input, output) {

    output$summary <- renderPrint({
    dataset <- rnorm(input$W)
    summary(dataset)
  })
  output$description <- renderText({
    sprintf("This tool allows researchers to design experiments for CNV detection from sequencing data. ")
  })
    
  output$toPrint <- renderText({
      sprintf(" L: %s \n D: %s \n W: %s \n N: %s \n l: %s \n F: %s \n" ,input$L, input$D, input$W, input$N, input$l, input$F)
  })
  # Show the first "n" observations
  output$view <- renderTable({

  })
 
  output$power <- renderText({
    if(input$dis =="pois"){
      if(input$F == 1){
        paste("Power is", poisson_pure(as.numeric(input$alpha), as.numeric(input$W), as.numeric(input$l), as.numeric(input$D),as.numeric(input$N), as.numeric(input$L)))
      }else{
        paste("Power is", poisson_impurity(as.numeric(input$alpha), as.numeric(input$W), as.numeric(input$l), as.numeric(input$D),as.numeric(input$N), as.numeric(input$L),input$F))
      }  
    }
    else if(input$dis == "nb"){
      if(is.na(input$phi)) {
        phi = 1/(input$W*input$D/input$l)
      }else{
        phi = input$phi
      }
      negativebinomial_purity(input$alpha,input$W,input$l,input$D,input$N,input$L,phi,input$F)
    }
  })
  ##########################
  ##Generate a power curve##
  output$powerCurve <- renderPlot({plot_power(input)});
  output$powerCurve_text <- renderText({
    if(input$tovary == 'None'){
      sprintf("Select parameter to vary to plot a power curve")
    }
    else{
      sprintf(" L: %s \n D: %s \n W: %s \n N: %s \n l: %s \n" ,input$L, input$D, input$W, input$N, input$l)
    }
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(input$dataset, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(powervalues(input),file, row.names = FALSE)
    }
  )
})

shinyApp(ui, server)
