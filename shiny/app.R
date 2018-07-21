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
                 actionButton("preset_00", "Default"),
                 actionButton("preset_01", "Scenario 01"),
                 actionButton("preset_02", "Scenario 02"),
                 actionButton("preset_03", "Scenario 03"),
                 actionButton("preset_04", "Scenario 04"),
                 actionButton("preset_05", "Scenario 05"),
                 actionButton("preset_06", "Scenario 06"),
                 helpText(),
                 textOutput("preset_text"),
                 helpText(),
                 verbatimTextOutput("toPrint"),
                 plotOutput("powerCurve",width="98%"),
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
    calculated_power <- UA_CPA (as.numeric(input$alpha), as.numeric(input$W), as.numeric(input$l), as.numeric(input$D), as.numeric(input$N), as.numeric(input$L), as.numeric(input$phi), as.numeric(input$F))
    # To print out a list of the input parameters and the calculated power
    paste(sprintf("Distribution: %s\nTumor Purity, F: %s\nSignificance Level, alpha: %s\nLength of CNV, L: %s\nPloidy of CNV, N: %s\nWindow length, W: %s\nLength of read, l: %s\nHaploid Read Depth, D: %s\nOverdispersion parameter, theta: %s\n\n",
              input$dis, input$F, input$alpha, input$L, input$N, input$W, input$l, input$D, input$phi), paste("\nThe power for this experiment design is", format(round(calculated_power, 2), nsmall = 2)))
  })

  # Generate a power curve
  output$powerCurve <- renderPlot({UA_CPAP(input)});

    # Adopt preset 00
  observeEvent(input$preset_00, {
    output$preset_text <- renderText({  # text to describe the preset
      paste("")
    })
    updateNumericInput(session, "F", value = 1)
    updateNumericInput(session, "alpha", value = 0.05)
    updateNumericInput(session, "L", value = 10000)
    updateNumericInput(session, "N", value = 3)
    updateNumericInput(session, "W", value = 100)
    updateNumericInput(session, "l", value = 100)
    updateNumericInput(session, "D", value = 0.5)
    updateNumericInput(session, "phi", value = 0)
    updateSelectInput(session, "tovary", selected = "Length of CNV")
    updateNumericInput(session, "min1", value = 500)
    updateNumericInput(session, "max1", value = 15000)
    updateSelectInput(session,"tovary2", selected = "None")
    updateNumericInput(session, "val1", value = 1)
    updateNumericInput(session, "val2", value = 2)
    updateNumericInput(session, "val3", value = 3)
    updateCheckboxInput(session, "log_scale", value = FALSE)
  })
  
  # Adopt preset 01
  observeEvent(input$preset_01, {
    output$preset_text <- renderText({  # text to describe the preset
      paste("When sequencing tumor-derived DNA, sometimes one faces the situation of low budget or low amounts of tissue material, or one wishes to use shallow sequencing to gain a preliminary glimpse of the tumor. These situations share the common constraint of small D. As we show in the following graph, for N=3, l=100, W=1000, and α=0.05, sequencing with 1X read depth can detect 10kb CNV near 100% of the times, while sequencing with 0.5X read depth can also detect 10kb CNV fairly accurately (power = 0.95). Sequencing with 0.1X read depth will miss most of the 10kb CNV (power = 0.34); however, it can still be used to detect CNV longer than 45kb (power = 0.97).")
    })
    updateNumericInput(session, "F", value = 1)
    updateNumericInput(session, "alpha", value = 0.05)
    updateNumericInput(session, "L", value = 10000)
    updateNumericInput(session, "N", value = 3)
    updateNumericInput(session, "W", value = 1000)
    updateNumericInput(session, "l", value = 100)
    updateNumericInput(session, "D", value = 0.5)
    updateNumericInput(session, "phi", value = 0)
    updateSelectInput(session, "tovary", selected = "Length of CNV")
    updateNumericInput(session, "min1", value = 5000)
    updateNumericInput(session, "max1", value = 100000)
    updateSelectInput(session,"tovary2", selected = "Read Depth")
    updateNumericInput(session, "val1", value = 1)
    updateNumericInput(session, "val2", value = 0.5)
    updateNumericInput(session, "val3", value = 0.1)
    updateCheckboxInput(session, "log_scale", value = TRUE)
  })
  
  # Adopt preset 02
  observeEvent(input$preset_02, {
    output$preset_text <- renderText({  # text to describe the preset
      paste("Sometimes the researcher already knows the regions of the genome to focus on and is interested in detecting very rare cell populations with copy number aberrations in these regions, even if it requires ultra-deep targeted sequencing. Such needs arise in (1) monitoring of potential recurrence of a particular cancer clone during remission and hoping to achieve high sensitivity with a large D; (2) monitoring gene-specific aberrations in circulating tumor cells or circulating tumor DNA when a predefined list of oncogenes is particularly relevant for a particular cancer type (such as BRAF in melanoma). These practical demands share the common feature of small F and large D. As we show in the following graph, for N=3, l=100, W=1000, and α=0.05, sequencing with 1000X read depth can detect 10kb CNV in a sample contains as low as 2% interested cells (power = 0.97), 100X read depth can detect 10kb CNV in a sample contains 6% interested cells (power = 0.96), and 30X read depth can only detect 10kb CNV in a sample contains more than 12% interested cells (power = 0.97).")
    })
    updateNumericInput(session, "F", value = 0.02)
    updateNumericInput(session, "alpha", value = 0.05)
    updateNumericInput(session, "L", value = 10000)
    updateNumericInput(session, "N", value = 3)
    updateNumericInput(session, "W", value = 1000)
    updateNumericInput(session, "l", value = 100)
    updateNumericInput(session, "D", value = 1000)
    updateNumericInput(session, "phi", value = 0)
    updateSelectInput(session, "tovary", selected = "Sample purity")
    updateNumericInput(session, "min1", value = 0)
    updateNumericInput(session, "max1", value = 0.2)
    updateSelectInput(session,"tovary2", selected = "Read Depth")
    updateNumericInput(session, "val1", value = 1000)
    updateNumericInput(session, "val2", value = 100)
    updateNumericInput(session, "val3", value = 30)
    updateCheckboxInput(session, "log_scale", value = FALSE)
  })
  
  # Adopt preset 03
  observeEvent(input$preset_03, {
    output$preset_text <- renderText({  # text to describe the preset
      paste("This example with N=3, l=100, W=1000, and α=0.05 shows power as a function of CNA length L and read depth D under a Poisson distribution. For a CNA with L=10kb, read depth of 0.1 shows very low power (0.34), while a read depth of 1 would increase power significantly (1.00). With today’s sequencing technology, 1X coverage of a mammalian genome costs $30-40. In comparison, a 3-million SNP genotyping platform will measure about 10 SNP loci in this CNA, but at a higher cost.")
    })
    updateNumericInput(session, "F", value = 1)
    updateNumericInput(session, "alpha", value = 0.05)
    updateNumericInput(session, "L", value = 10000)
    updateNumericInput(session, "N", value = 3)
    updateNumericInput(session, "W", value = 1000)
    updateNumericInput(session, "l", value = 100)
    updateNumericInput(session, "D", value = 1)
    updateNumericInput(session, "phi", value = 0)
    updateSelectInput(session, "tovary", selected = "Length of CNV")
    updateNumericInput(session, "min1", value = 5000)
    updateNumericInput(session, "max1", value = 100000)
    updateSelectInput(session,"tovary2", selected = "Read Depth")
    updateNumericInput(session, "val1", value = 1)
    updateNumericInput(session, "val2", value = 0.5)
    updateNumericInput(session, "val3", value = 0.1)
    updateCheckboxInput(session, "log_scale", value = TRUE)
  })
  
  # Adopt preset 04
  observeEvent(input$preset_04, {
    output$preset_text <- renderText({  # text to describe the preset
      paste("This example is exact like our preset 03 except that this example is under the assumption of Negative Binomial model with an overdispersion parameter of 1. For the same CNA with L=10kb, the power decreases from 1.00 to 0.95 at the read depth of 1. The decrease in power here is due to the increase in overall variance and can be compensated by higher sequencing depth if desired.")
    })
    updateNumericInput(session, "F", value = 1)
    updateNumericInput(session, "alpha", value = 0.05)
    updateNumericInput(session, "L", value = 10000)
    updateNumericInput(session, "N", value = 3)
    updateNumericInput(session, "W", value = 1000)
    updateNumericInput(session, "l", value = 100)
    updateNumericInput(session, "D", value = 1)
    updateNumericInput(session, "phi", value = 1)
    updateSelectInput(session, "tovary", selected = "Length of CNV")
    updateNumericInput(session, "min1", value = 5000)
    updateNumericInput(session, "max1", value = 100000)
    updateSelectInput(session,"tovary2", selected = "Read Depth")
    updateNumericInput(session, "val1", value = 1)
    updateNumericInput(session, "val2", value = 0.5)
    updateNumericInput(session, "val3", value = 0.1)
    updateCheckboxInput(session, "log_scale", value = TRUE)
  })
  
  # Adopt preset 05
  observeEvent(input$preset_05, {
    output$preset_text <- renderText({  # text to describe the preset
      paste("This example is corresponding to the figure 3 in our paper. It extends our earlier example from preset 3 by showing the effect on power when the sample purity is reduced from 1 to 0.5 and 0.1. For the same CNA with L=10kb at the read depth of 1, the power decreases from 1.00 to 0.80 when the sample purity is 0.5, to 0.08 when the sample purity is 0.1. Again, such decrease in power can be compensated by higher sequencing depth.")
    })
    updateNumericInput(session, "F", value = 0.5)
    updateNumericInput(session, "alpha", value = 0.05)
    updateNumericInput(session, "L", value = 10000)
    updateNumericInput(session, "N", value = 3)
    updateNumericInput(session, "W", value = 1000)
    updateNumericInput(session, "l", value = 100)
    updateNumericInput(session, "D", value = 1)
    updateNumericInput(session, "phi", value = 0)
    updateSelectInput(session, "tovary", selected = "Length of CNV")
    updateNumericInput(session, "min1", value = 5000)
    updateNumericInput(session, "max1", value = 100000)
    updateSelectInput(session,"tovary2", selected = "Sample purity")
    updateNumericInput(session, "val1", value = 1)
    updateNumericInput(session, "val2", value = 0.5)
    updateNumericInput(session, "val3", value = 0.1)
    updateCheckboxInput(session, "log_scale", value = TRUE)
  })
  
  # Adopt preset 06
  observeEvent(input$preset_06, {
    output$preset_text <- renderText({  # text to describe the preset
      paste("This example is exact like our preset 03 except that this example is under the assumption of Negative Binomial model with an overdispersion parameter of 1. For the same CNA with L=10kb and F=0.1, the power decreases from 0.80 to 0.49 at the read depth of 1. The decrease in power here is due to the increase in overall variance and can be compensated by higher sequencing depth if desired.")
    })
    updateNumericInput(session, "F", value = 0.5)
    updateNumericInput(session, "alpha", value = 0.05)
    updateNumericInput(session, "L", value = 10000)
    updateNumericInput(session, "N", value = 3)
    updateNumericInput(session, "W", value = 1000)
    updateNumericInput(session, "l", value = 100)
    updateNumericInput(session, "D", value = 1)
    updateNumericInput(session, "phi", value = 1)
    updateSelectInput(session, "tovary", selected = "Length of CNV")
    updateNumericInput(session, "min1", value = 5000)
    updateNumericInput(session, "max1", value = 100000)
    updateSelectInput(session,"tovary2", selected = "Sample purity")
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
