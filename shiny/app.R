# Load the required libraries
library(shiny)
library(markdown)
library(ggplot2)
library(grid)

# Source the functions from the coresponding files
source("functions.R")

# Define UI for CNV calculator ----
ui <- fluidPage(
  
  # Application title
  titlePanel("CNV Detection Power Calculator"),
  
  sidebarLayout(
    # Sidebar for users to select a data distribution and specify the parameters
    sidebarPanel(
      tabsetPanel(
        tabPanel("Power Calculator",
                 helpText(" "),
                 helpText("Enter experimental parameters to calculate power"),
                 radioButtons('dis', label = 'Distribution*', c('Poisson' = 'Poisson','Negative Binomial' = 'Negative Binomial')),
                 numericInput("F", "Tumor Purity (range 0-1)*:", 1,min=0,max=1),
                 numericInput("alpha", "Significance Level, alpha (range 0-1)*:", 0.05,min=0,max=1),
                 numericInput("L", "Length of CNV (in bases)*:", 10000,min=0),
                 numericInput("N", "Ploidy of CNV*:", 3,min=0),
                 numericInput("W", "Window length (in bases)*:", 100,min=0),
                 numericInput("l", "Length of read (in bases)*:", 100,min=0),
                 numericInput("D", "Haploid Read Depth*: ", 0.5,min=0),
                 numericInput("phi","Overdispersion parameter, theta (0 for Poisson Distribution):", min=0, value = 0)
                 ),
        tabPanel("Power Curve",
                 helpText(" "),
                 helpText("Select parameter to vary to generate a power curve"),
                 selectInput('tovary',label="First parameter to vary along x-axis", choices=c('None','Read Depth','Read length','Window size','Length of CNV','Sample purity'),selected = "Length of CNV"),
                 numericInput("min1","Minimum value of parameter to be tested",500),
                 numericInput("max1","Maximum value of parameter to be tested",15000),
                 selectInput('tovary2',label="Second Parameter to Vary (line color)", choices=c('None','Read Depth','Read length','Window size','Length of CNV','Sample purity')),
                 numericInput("val1","Value 1 to be tested for parameter 2",1),
                 numericInput("val2","Value 2 to be tested for parameter 2",2),
                 numericInput("val3","Value 3 to be tested for parameter 2",3),
                 checkboxInput("log_scale", "log Scale for x-Axis",FALSE)
                 )
      )
      ),

    # Main Panel to display a summary of inputs, the calculated power, and a panel for users to generate their own power curve
    mainPanel(
      tabsetPanel(
        # To print a summary of the input parameters, the calculated power, and power curve all in one page
        tabPanel("Summary",
                 helpText(),
                 actionButton("preset_00", "Preset 00"),
                 actionButton("preset_01", "Preset 01"),
                 actionButton("preset_02", "Preset 02"),
                 actionButton("preset_03", "Preset 03"),
                 actionButton("preset_04", "Preset 04"),
                 helpText(),
                 verbatimTextOutput("toPrint"),
                 plotOutput("powerCurve",width="98%"),
                 downloadButton("downloadData", "Download Raw Data"),  # button for download data
                 downloadButton("downloadGraph", "Download Graph")  # button for download graph
        ),
        tabPanel("Tutorial", includeMarkdown("tutorial.md")),
        # To display more information about the website
        tabPanel("About", includeMarkdown("README.md")),
        # To display the patch notes
        tabPanel("Patch Notes", includeMarkdown("patch_notes.md"))
                
    )
  )
)
)


# Define the Server function for CNV calculator
server = function(input, output, session) {

  # To print a summary of input parameters and the calcualted power
  output$toPrint <- renderText({
    # To calcualte the power coresponding to the inputs
    calculated_power <- UA_CPA (as.numeric(input$alpha), as.numeric(input$W), as.numeric(input$l), as.numeric(input$D), as.numeric(input$N), as.numeric(input$L), as.numeric(input$phi), as.numeric(input$F))
    # To print out a list of the input parameters and the calculated power
    paste(sprintf("Distribution: %s\nTumor Purity, F: %s\nSignificance Level, alpha: %s\nLength of CNV, L: %s\nPloidy of CNV, N: %s\nWindow length, W: %s\nLength of read, l: %s\nHaploid Read Depth, D: %s\nOverdispersion parameter, theta: %s\n\n",
              input$dis, input$F, input$alpha, input$L, input$N, input$W, input$l, input$D, input$phi), paste("\nThe power for this experiment design is", format(round(calculated_power, 2), nsmall = 2)))
  })

  # Generate a power curve
  output$powerCurve <- renderPlot({UA_CPAP(input)});

  # Adopt preset 00
  observeEvent(input$preset_00, {
    updateNumericInput(session, "F", value = 1)
    updateNumericInput(session, "alpha", value = 0.05)
    updateNumericInput(session, "L", value = 10000)
    updateNumericInput(session, "N", value = 3)
    updateNumericInput(session, "W", value = 100)
    updateNumericInput(session, "l", value = 100)
    updateNumericInput(session, "D", value = 0.5)
    updateNumericInput(session, "phi", value = 0)
    updateSelectInput(session, "tovary", choices = "Length of CNV")
    updateNumericInput(session, "min1", value = 500)
    updateNumericInput(session, "max1", value = 15000)
    updateSelectInput(session,"tovary2", choices = "None")
    updateNumericInput(session, "val1", value = 1)
    updateNumericInput(session, "val2", value = 2)
    updateNumericInput(session, "val3", value = 3)
    updateCheckboxInput(session, "log_scale", value = FALSE)
  })
  
  # Adopt preset 01
  observeEvent(input$preset_01, {
    updateNumericInput(session, "F", value = 1)
    updateNumericInput(session, "alpha", value = 0.05)
    updateNumericInput(session, "L", value = 10000)
    updateNumericInput(session, "N", value = 3)
    updateNumericInput(session, "W", value = 1000)
    updateNumericInput(session, "l", value = 100)
    updateNumericInput(session, "D", value = 0.1)
    updateNumericInput(session, "phi", value = 0)
    updateSelectInput(session, "tovary", choices = "Length of CNV")
    updateNumericInput(session, "min1", value = 5000)
    updateNumericInput(session, "max1", value = 100000)
    updateSelectInput(session,"tovary2", choices = "Read Depth")
    updateNumericInput(session, "val1", value = 1)
    updateNumericInput(session, "val2", value = 0.5)
    updateNumericInput(session, "val3", value = 0.1)
    updateCheckboxInput(session, "log_scale", value = TRUE)
  })
  
  # Adopt preset 02
  observeEvent(input$preset_02, {
    updateNumericInput(session, "F", value = 1)
    updateNumericInput(session, "alpha", value = 0.05)
    updateNumericInput(session, "L", value = 10000)
    updateNumericInput(session, "N", value = 3)
    updateNumericInput(session, "W", value = 1000)
    updateNumericInput(session, "l", value = 100)
    updateNumericInput(session, "D", value = 0.1)
    updateNumericInput(session, "phi", value = 1)
    updateSelectInput(session, "tovary", choices = "Length of CNV")
    updateNumericInput(session, "min1", value = 5000)
    updateNumericInput(session, "max1", value = 500000)
    updateSelectInput(session,"tovary2", choices = "Read Depth")
    updateNumericInput(session, "val1", value = 1)
    updateNumericInput(session, "val2", value = 0.5)
    updateNumericInput(session, "val3", value = 0.1)
    updateCheckboxInput(session, "log_scale", value = TRUE)
  })
  
  # Adopt preset 03
  observeEvent(input$preset_03, {
    updateNumericInput(session, "F", value = 0.1)
    updateNumericInput(session, "alpha", value = 0.05)
    updateNumericInput(session, "L", value = 10000)
    updateNumericInput(session, "N", value = 3)
    updateNumericInput(session, "W", value = 1000)
    updateNumericInput(session, "l", value = 100)
    updateNumericInput(session, "D", value = 0.1)
    updateNumericInput(session, "phi", value = 0)
    updateSelectInput(session, "tovary", choices = "Length of CNV")
    updateNumericInput(session, "min1", value = 5000)
    updateNumericInput(session, "max1", value = 10000000)
    updateSelectInput(session,"tovary2", choices = "Sample purity")
    updateNumericInput(session, "val1", value = 1)
    updateNumericInput(session, "val2", value = 0.5)
    updateNumericInput(session, "val3", value = 0.1)
    updateCheckboxInput(session, "log_scale", value = TRUE)
  })
  
  # Adopt preset 04
  observeEvent(input$preset_04, {
    updateNumericInput(session, "F", value = 0.1)
    updateNumericInput(session, "alpha", value = 0.05)
    updateNumericInput(session, "L", value = 10000)
    updateNumericInput(session, "N", value = 3)
    updateNumericInput(session, "W", value = 1000)
    updateNumericInput(session, "l", value = 100)
    updateNumericInput(session, "D", value = 0.1)
    updateNumericInput(session, "phi", value = 1)
    updateSelectInput(session, "tovary", choices = "Length of CNV")
    updateNumericInput(session, "min1", value = 5000)
    updateNumericInput(session, "max1", value = 1000000)
    updateSelectInput(session,"tovary2", choices = "Sample purity")
    updateNumericInput(session, "val1", value = 1)
    updateNumericInput(session, "val2", value = 0.5)
    updateNumericInput(session, "val3", value = 0.1)
    updateCheckboxInput(session, "log_scale", value = TRUE)
  })
  
  # Allow the user to download the data
  output$downloadData <- downloadHandler(
    filename = "raw_data.csv",
    content = function(file) {
      write.csv(
        x = df_1, file) # df_1 is a global variable of data frame that is generated in the function UA_CPAP
    }
  )
  
  # Allow the user to download the graph
  output$downloadGraph <- downloadHandler(
    filename = "graph.png",
    content = function(file) {
      ggsave(file)
    }
  )
}

shinyApp(ui, server)
