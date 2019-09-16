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
                 radioButtons('dis', label = 'Distribution*', c('Poisson' = 'Poisson','Negative Binomial with theta' = 'Negative Binomial with theta','Negative Binomial with phi' = 'Negative Binomial with phi')),
                 numericInput("F", "Tumor Purity (range 0-1)*:", 1,min=0,max=1),
                 numericInput("alpha", "Significance Level, alpha (range 0-1)*:", 0.05,min=0,max=1),
                 numericInput("L", "Length of CNV (in bases)*:", 5000,min=0),
                 numericInput("N", "Ploidy of CNV*:", 3,min=0),
                 
                 numericInput("l", "Length of read (in bases)*:", 100,min=0),
                 numericInput("D", "Haploid Read Depth*: ", 1,min=0),
                 numericInput("theta_1","Overdispersion parameter, theta (0 for Poisson Distribution):", min=0, value = 0),
                 numericInput("phi_1","Overdispersion parameter, phi (0 for Poisson Distribution):", min=0, value = 0)
                 ),
        tabPanel("Power Curve",
                 helpText(" "),
                 helpText("Select parameter to vary to generate a power curve"),
                 selectInput('tovary',label="First parameter to vary along x-axis", choices=c('None','Read Depth','Read length','Window size','Length of CNV','Sample purity',"Overdispersion parameter"),selected = "Length of CNV"),
                 numericInput("min1","Minimum value of parameter to be tested",500),
                 numericInput("max1","Maximum value of parameter to be tested",10000),
                 selectInput('tovary2',label="Second Parameter to Vary (line color)", choices=c('None','Read Depth','Read length','Window size','Length of CNV','Sample purity',"Overdispersion parameter")),
                 numericInput("val1","Value 1 to be tested for parameter 2",1),
                 numericInput("val2","Value 2 to be tested for parameter 2",2),
                 numericInput("val3","Value 3 to be tested for parameter 2",4),
                 selectInput('tovary3',label="Third Parameter to Vary (line type)", choices=c('None','Read Depth','Read length','Window size','Length of CNV','Sample purity',"Overdispersion parameter")),
                 numericInput("val31","Value 1 to be tested for parameter 3",0.25),
                 numericInput("val32","Value 2 to be tested for parameter 3",0.5),
                 numericInput("val33","Value 3 to be tested for parameter 3",1)
                 ),
        tabPanel("Advanced",
                 helpText(" "),
                 helpText("Advanced setings for power curve"),
                 checkboxInput("log_scale", "Use log-scale for x-axis",FALSE),
                 checkboxInput("fix_y_range", "Fix range for y-axis to 0-1",TRUE),
                 checkboxInput("show_x_at_0", "Include x=0 in the graph (for continuous-scale)",FALSE),
                 numericInput("W", "Window length (in bases)*:", 1,min=0)
                 )
      )
      ),

    # Main Panel to display a summary of inputs, the calculated power, and a panel for users to generate their own power curve
    mainPanel(
      tabsetPanel(
        # To print a summary of the input parameters, the calculated power, and power curve all in one page
        tabPanel("Summary",
                 helpText(),
                 actionButton("preset_00", "Default"),
                 helpText(),
                 actionButton("preset_01", "Figure 1"),
                 actionButton("preset_02", "Figure 2"),
                 actionButton("preset_03", "Figure 3"),
                 actionButton("preset_04", "Figure 4"),
                 helpText(),
                 actionButton("preset_05", "Use Case 1a"),
                 actionButton("preset_06", "USe Case 1b"),
                 helpText(),
                 actionButton("preset_07", "Use Case 2a"),
                 actionButton("preset_08", "USe Case 2b"),
                 helpText(),
                 textOutput("preset_text"),
                 helpText(),
                 verbatimTextOutput("toPrint"),
                 plotOutput("powerCurve",width="98%",hover = hoverOpts(id = "plot_hover", delayType = "throttle")),
                 verbatimTextOutput("plot_hoverinfo"),
                 downloadButton("downloadData", "Download Raw Data"),  # button for download data
                 downloadButton("downloadGraph", "Download Graph")  # button for download graph
        ),
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
    calculated_power <- UA_CPA (as.numeric(input$alpha), as.numeric(input$W), as.numeric(input$l), as.numeric(input$D), as.numeric(input$N), as.numeric(input$L), as.numeric(input$theta_1), as.numeric(input$F), input)
    # To print out a list of the input parameters and the calculated power
    paste(sprintf("Distribution: %s\nTumor Purity, F: %s\nSignificance Level, alpha: %s\nLength of CNV, L: %s\nPloidy of CNV, N: %s\nWindow length, W: %s\nLength of read, l: %s\nHaploid Read Depth, D: %s\nOverdispersion parameter, theta: %s\nOverdispersion parameter, phi: %s\n\n",
              input$dis, input$F, input$alpha, input$L, input$N, input$W, input$l, input$D, input$theta_1, input$phi_1), paste("\nThe power for this experiment design is", format(round(calculated_power, 2), nsmall = 2)))
  })

  output$plot_hoverinfo <- renderPrint({
    if (is.null(input$plot_hover)) {
      cat("x=NA"," y=NA",sep = "")
    } else {
      cat("x=", round(input$plot_hover$x, 2), " y=", round(input$plot_hover$y, 2),sep = "")
    }
  })
  
  # Generate a power curve
  output$powerCurve <- renderPlot({UA_CPAP(input)});

    # Adopt preset 00
  observeEvent(input$preset_00, {
    output$preset_text <- renderText({  # text to describe the preset
      paste("")
    })
    updateRadioButtons(session, "dis", selected = "Poisson")
    updateNumericInput(session, "F", value = 1)
    updateNumericInput(session, "alpha", value = 0.05)
    updateNumericInput(session, "L", value = 5000)
    updateNumericInput(session, "N", value = 3)
    updateNumericInput(session, "W", value = 1)
    updateNumericInput(session, "l", value = 100)
    updateNumericInput(session, "D", value = 1)
    updateNumericInput(session, "theta_1", value = 0)
    updateSelectInput(session, "tovary", selected = "Length of CNV")
    updateNumericInput(session, "min1", value = 500)
    updateNumericInput(session, "max1", value = 10000)
    updateSelectInput(session,"tovary2", selected = "None")
    updateNumericInput(session, "val1", value = 1)
    updateNumericInput(session, "val2", value = 2)
    updateNumericInput(session, "val3", value = 4)
    updateSelectInput(session,"tovary3", selected = "None")
    updateNumericInput(session, "val31", value = 0.25)
    updateNumericInput(session, "val32", value = 0.5)
    updateNumericInput(session, "val33", value = 1)
    updateCheckboxInput(session, "log_scale", value = FALSE)
    updateCheckboxInput(session, "fix_y_range", value = TRUE)
    updateCheckboxInput(session, "show_x_at_0", value = FALSE)
  })
  
  # Adopt preset 01
  observeEvent(input$preset_01, {
    output$preset_text <- renderText({  # text to describe the preset
    })
    updateRadioButtons(session, "dis", selected = "Poisson")
    updateNumericInput(session, "F", value = 1)
    updateNumericInput(session, "alpha", value = 0.05)
    updateNumericInput(session, "L", value = 5000)
    updateNumericInput(session, "N", value = 3)
    updateNumericInput(session, "W", value = 1)
    updateNumericInput(session, "l", value = 100)
    updateNumericInput(session, "D", value = 1)
    updateNumericInput(session, "theta_1", value = 0)
    updateSelectInput(session, "tovary", selected = "Length of CNV")
    updateNumericInput(session, "min1", value = 5000)
    updateNumericInput(session, "max1", value = 100000)
    updateSelectInput(session,"tovary2", selected = "Read Depth")
    updateNumericInput(session, "val1", value = 0.1)
    updateNumericInput(session, "val2", value = 0.5)
    updateNumericInput(session, "val3", value = 1)
    updateSelectInput(session,"tovary3", selected = "None")
    updateNumericInput(session, "val31", value = 0.25)
    updateNumericInput(session, "val32", value = 0.5)
    updateNumericInput(session, "val33", value = 1)
    updateCheckboxInput(session, "log_scale", value = TRUE)
    updateCheckboxInput(session, "fix_y_range", value = TRUE)
    updateCheckboxInput(session, "show_x_at_0", value = FALSE)
  })
  
  # Adopt preset 02
  observeEvent(input$preset_02, {
    output$preset_text <- renderText({  # text to describe the preset
    })
    updateRadioButtons(session, "dis", selected = "Poisson")
    updateNumericInput(session, "F", value = 1)
    updateNumericInput(session, "alpha", value = 0.05)
    updateNumericInput(session, "L", value = 5000)
    updateNumericInput(session, "N", value = 3)
    updateNumericInput(session, "W", value = 1)
    updateNumericInput(session, "l", value = 100)
    updateNumericInput(session, "D", value = 1)
    updateNumericInput(session, "theta_1", value = 0)
    updateSelectInput(session, "tovary", selected = "Length of CNV")
    updateNumericInput(session, "min1", value = 5000)
    updateNumericInput(session, "max1", value = 500000)
    updateSelectInput(session,"tovary2", selected = "Read Depth")
    updateNumericInput(session, "val1", value = 0.1)
    updateNumericInput(session, "val2", value = 0.5)
    updateNumericInput(session, "val3", value = 1)
    updateSelectInput(session,"tovary3", selected = "Overdispersion parameter")
    updateNumericInput(session, "val31", value = 0)
    updateNumericInput(session, "val32", value = 0.5)
    updateNumericInput(session, "val33", value = 1)
    updateCheckboxInput(session, "log_scale", value = TRUE)
    updateCheckboxInput(session, "fix_y_range", value = TRUE)
    updateCheckboxInput(session, "show_x_at_0", value = FALSE)
  })
  
  # Adopt preset 03
  observeEvent(input$preset_03, {
    output$preset_text <- renderText({  # text to describe the preset
    })
    updateRadioButtons(session, "dis", selected = "Poisson")
    updateNumericInput(session, "F", value = 1)
    updateNumericInput(session, "alpha", value = 0.05)
    updateNumericInput(session, "L", value = 10000)
    updateNumericInput(session, "N", value = 3)
    updateNumericInput(session, "W", value = 1)
    updateNumericInput(session, "l", value = 100)
    updateNumericInput(session, "D", value = 1)
    updateNumericInput(session, "theta_1", value = 0)
    updateSelectInput(session, "tovary", selected = "Length of CNV")
    updateNumericInput(session, "min1", value = 5000)
    updateNumericInput(session, "max1", value = 5000000)
    updateSelectInput(session,"tovary2", selected = "Read Depth")
    updateNumericInput(session, "val1", value = 0.1)
    updateNumericInput(session, "val2", value = 0.5)
    updateNumericInput(session, "val3", value = 1)
    updateSelectInput(session,"tovary3", selected = "Sample purity")
    updateNumericInput(session, "val31", value = 0.1)
    updateNumericInput(session, "val32", value = 0.5)
    updateNumericInput(session, "val33", value = 1)
    updateCheckboxInput(session, "log_scale", value = TRUE)
    updateCheckboxInput(session, "fix_y_range", value = TRUE)
    updateCheckboxInput(session, "show_x_at_0", value = FALSE)
  })
  
  # Adopt preset 04
  observeEvent(input$preset_04, {
    output$preset_text <- renderText({  # text to describe the preset
    })
    updateRadioButtons(session, "dis", selected = "Poisson")
    updateNumericInput(session, "F", value = 1)
    updateNumericInput(session, "alpha", value = 0.05)
    updateNumericInput(session, "L", value = 5000)
    updateNumericInput(session, "N", value = 3)
    updateNumericInput(session, "W", value = 1)
    updateNumericInput(session, "l", value = 100)
    updateNumericInput(session, "D", value = 1)
    updateNumericInput(session, "theta_1", value = 0)
    updateSelectInput(session, "tovary", selected = "Length of CNV")
    updateNumericInput(session, "min1", value = 5000)
    updateNumericInput(session, "max1", value = 1000000)
    updateSelectInput(session,"tovary2", selected = "Sample purity")
    updateNumericInput(session, "val1", value = 0.1)
    updateNumericInput(session, "val2", value = 0.5)
    updateNumericInput(session, "val3", value = 1)
    updateSelectInput(session,"tovary3", selected = "Overdispersion parameter")
    updateNumericInput(session, "val31", value = 0)
    updateNumericInput(session, "val32", value = 1)
    updateNumericInput(session, "val33", value = 0)
    updateCheckboxInput(session, "log_scale", value = TRUE)
    updateCheckboxInput(session, "fix_y_range", value = TRUE)
    updateCheckboxInput(session, "show_x_at_0", value = FALSE)
  })
  
  # Adopt preset 05
  observeEvent(input$preset_05, {
    output$preset_text <- renderText({  # text to describe the preset
    })
    updateRadioButtons(session, "dis", selected = "Poisson")
    updateNumericInput(session, "F", value = 1)
    updateNumericInput(session, "alpha", value = 0.05)
    updateNumericInput(session, "L", value = 10000)
    updateNumericInput(session, "N", value = 3)
    updateNumericInput(session, "W", value = 1)
    updateNumericInput(session, "l", value = 100)
    updateNumericInput(session, "D", value = 1)
    updateNumericInput(session, "theta_1", value = 0)
    updateSelectInput(session, "tovary", selected = "Length of CNV")
    updateNumericInput(session, "min1", value = 5000)
    updateNumericInput(session, "max1", value = 10000000)
    updateSelectInput(session,"tovary2", selected = "Read Depth")
    updateNumericInput(session, "val1", value = 0.1)
    updateNumericInput(session, "val2", value = 0.5)
    updateNumericInput(session, "val3", value = 1)
    updateSelectInput(session,"tovary3", selected = "Sample purity")
    updateNumericInput(session, "val31", value = 0.01)
    updateNumericInput(session, "val32", value = 0.1)
    updateNumericInput(session, "val33", value = 0.5)
    updateCheckboxInput(session, "log_scale", value = TRUE)
    updateCheckboxInput(session, "fix_y_range", value = TRUE)
    updateCheckboxInput(session, "show_x_at_0", value = FALSE)
  })
  
  # Adopt preset 06
  observeEvent(input$preset_06, {
    output$preset_text <- renderText({  # text to describe the preset
    })
    updateRadioButtons(session, "dis", selected = "Negative Binomial with theta")
    updateNumericInput(session, "F", value = 1)
    updateNumericInput(session, "alpha", value = 0.05)
    updateNumericInput(session, "L", value = 10000)
    updateNumericInput(session, "N", value = 3)
    updateNumericInput(session, "W", value = 1)
    updateNumericInput(session, "l", value = 100)
    updateNumericInput(session, "D", value = 1)
    updateNumericInput(session, "theta_1", value = 1)
    updateSelectInput(session, "tovary", selected = "Length of CNV")
    updateNumericInput(session, "min1", value = 5000)
    updateNumericInput(session, "max1", value = 10000000)
    updateSelectInput(session,"tovary2", selected = "Read Depth")
    updateNumericInput(session, "val1", value = 0.1)
    updateNumericInput(session, "val2", value = 0.5)
    updateNumericInput(session, "val3", value = 1)
    updateSelectInput(session,"tovary3", selected = "Sample purity")
    updateNumericInput(session, "val31", value = 0.01)
    updateNumericInput(session, "val32", value = 0.1)
    updateNumericInput(session, "val33", value = 0.5)
    updateCheckboxInput(session, "log_scale", value = TRUE)
    updateCheckboxInput(session, "fix_y_range", value = TRUE)
    updateCheckboxInput(session, "show_x_at_0", value = FALSE)
  })
  
  observeEvent(input$preset_07, {
    output$preset_text <- renderText({  # text to describe the preset
    })
    updateRadioButtons(session, "dis", selected = "Poisson")
    updateNumericInput(session, "F", value = 1)
    updateNumericInput(session, "alpha", value = 0.05)
    updateNumericInput(session, "L", value = 10000)
    updateNumericInput(session, "N", value = 3)
    updateNumericInput(session, "W", value = 1)
    updateNumericInput(session, "l", value = 100)
    updateNumericInput(session, "D", value = 1)
    updateNumericInput(session, "theta_1", value = 0)
    updateSelectInput(session, "tovary", selected = "Sample purity")
    updateNumericInput(session, "min1", value = 0)
    updateNumericInput(session, "max1", value = 0.1)
    updateSelectInput(session,"tovary2", selected = "Read Depth")
    updateNumericInput(session, "val1", value = 10)
    updateNumericInput(session, "val2", value = 100)
    updateNumericInput(session, "val3", value = 1000)
    updateSelectInput(session,"tovary3", selected = "Length of CNV")
    updateNumericInput(session, "val31", value = 10000)
    updateNumericInput(session, "val32", value = 5000)
    updateNumericInput(session, "val33", value = 5000)
    updateCheckboxInput(session, "log_scale", value = FALSE)
    updateCheckboxInput(session, "fix_y_range", value = TRUE)
    updateCheckboxInput(session, "show_x_at_0", value = FALSE)
  })
  
  observeEvent(input$preset_08, {
    output$preset_text <- renderText({  # text to describe the preset
    })
    updateRadioButtons(session, "dis", selected = "Negative Binomial with theta")
    updateNumericInput(session, "F", value = 1)
    updateNumericInput(session, "alpha", value = 0.05)
    updateNumericInput(session, "L", value = 10000)
    updateNumericInput(session, "N", value = 3)
    updateNumericInput(session, "W", value = 1)
    updateNumericInput(session, "l", value = 100)
    updateNumericInput(session, "D", value = 1)
    updateNumericInput(session, "theta_1", value = 1)
    updateSelectInput(session, "tovary", selected = "Sample purity")
    updateNumericInput(session, "min1", value = 0)
    updateNumericInput(session, "max1", value = 0.1)
    updateSelectInput(session,"tovary2", selected = "Read Depth")
    updateNumericInput(session, "val1", value = 10)
    updateNumericInput(session, "val2", value = 100)
    updateNumericInput(session, "val3", value = 1000)
    updateSelectInput(session,"tovary3", selected = "Length of CNV")
    updateNumericInput(session, "val31", value = 10000)
    updateNumericInput(session, "val32", value = 5000)
    updateNumericInput(session, "val33", value = 5000)
    updateCheckboxInput(session, "log_scale", value = FALSE)
    updateCheckboxInput(session, "fix_y_range", value = TRUE)
    updateCheckboxInput(session, "show_x_at_0", value = FALSE)
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
