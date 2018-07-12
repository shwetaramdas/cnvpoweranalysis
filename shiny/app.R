# Load the required libraries
library(shiny)
library(markdown)
library(ggplot2)
library(grid)

# Source the functions from the coresponding files
source("functions.R")
source("functions_plot.R")

# Define UI for CNV calculator ----
ui <- pageWithSidebar(
  
  # Application title
  headerPanel("CNV Power Calculator"),
  
  # Sidebar for users to select a data distribution and specify the parameters
  sidebarPanel(
    helpText("Enter experimental parameters to calculate power"),
    radioButtons('dis', label = 'Distribution*', c('Poisson' = 'pois','Negative Binomial' = 'nb')),
    numericInput("F", "Tumor Purity (range 0-1)*:", 1,min=0,max=1),
    numericInput("alpha", "Significance Level, alpha (range 0-1)*:", 0.05,min=0,max=1),
    numericInput("L", "Length of CNV (in bases)*:", 10000,min=0),
    numericInput("N", "Ploidy of CNV*:", 3,min=0),
    numericInput("W", "Window length (in bases)*:", 1000,min=0),
    numericInput("l", "Length of read (in bases)*:", 100,min=0),
    numericInput("D", "Haploid Read Depth*: ", 30,min=0),
    numericInput("phi","Overdispersion parameter mu*phi or theta (0 for Poisson Distribution):", min=0, value = 0)
  ),
  
  # Main Panel to display a summary of inputs, the calculated power, and a panel for
  # users to generate their own power curve
  mainPanel(
    tabsetPanel(type = "tabs",
                # To print a summary of the input parameters and the calculated power
                tabPanel("Summary", verbatimTextOutput("toPrint")),
                # To generate a power curve
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
                # To display more information about the website
                tabPanel("About",includeMarkdown("README.md"))
                
    )
  )
)


# Define the Server function for CNV calculator
server = function(input, output) {

  # To print a summary of input parameters and the calcualted power
  output$toPrint <- renderText({
    # To calcualte the power coresponding to the inputs
    calculated_power <- UA_CPA (as.numeric(input$alpha), as.numeric(input$W), as.numeric(input$l),
                                as.numeric(input$D), as.numeric(input$N), as.numeric(input$L),
                                as.numeric(input$phi), as.numeric(input$F))
    # To print out a list of the input parameters and the calculated power
    paste(sprintf("Distribution: %s\nTumor Purity, F: %s\nSignificance Level, alpha: %s\nLength of CNV, L: %s\nPloidy of CNV, N: %s\nWindow length, W: %s\nLength of read, l: %s\nHaploid Read Depth, D: %s\nOverdispersion parameter, theta: %s\n\n",
              input$dis, input$F, input$alpha, input$L, input$N, input$W, input$l, input$D, input$phi), paste("\nThe power for this experiment design is", format(round(calculated_power, 2), nsmall = 2)))
  })

  # Generate a power curve
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
}

shinyApp(ui, server)
