library(shiny)

# Define UI for miles per gallon app ----
ui <- shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("CNV Power Calculator"),
  
  # Sidebar with controls to select a dataset and specify the number
  # of observations to view
  sidebarPanel(
    numericInput("obs", "Number of observations:", 10),
    numericInput("L", "Length of CNV:", 10000),
    numericInput("N", "Ploidy of CNV:", 3),
    numericInput("W", "Window length:", 3),
    selectInput("r","Overdispersion parameter r:",choices=c(0,0.05,0.1,0.25,0.5,0.75,1))
  ),
  
  # Show a summary of the dataset and an HTML table with the requested
  # number of observations
  mainPanel(
    verbatimTextOutput("summary"),
    
    tableOutput("view")
  )
))

server = shinyServer(function(input, output) {
  
  # Expression that generates a plot of the distribution. The expression
  # is wrapped in a call to renderPlot to indicate that:
  #
  #  1) It is "reactive" and therefore should be automatically 
  #     re-executed when inputs change
  #  2) Its output type is a plot 
  #
#  print(input$obs);

  output$summary <- renderPrint({
    dataset <- rnorm(input$W)
    summary(dataset)
  })
  
  # Show the first "n" observations
  output$view <- renderTable({

  })
})
shinyApp(ui, server)
