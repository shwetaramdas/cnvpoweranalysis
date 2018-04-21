library(shiny)
library(markdown)
source("functions.R")


# Define UI for CNV calculator ----
ui <- shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("CNV Power Calculator"),
  
  # Sidebar with controls to select a dataset and specify the number
  # of observations to view
  sidebarPanel(
    helpText("Enter experimental parameters to calculate power"),
    radioButtons('dis', label = 'Distribution*', c('Poisson' = 'pois','Negative Binomial' = 'nb')),
    numericInput("F", "Tumor Purity*:", 1,min=0,max=1),
    numericInput("alpha", "Significance Level alpha*:", 0.05,min=0,max=1),
    numericInput("L", "Length of CNV*:", 10000,min=0),
    numericInput("N", "Ploidy of CNV*:", 3,min=0),
    numericInput("W", "Window length*:", 100,min=0),
    numericInput("l", "Length of read*:", 100,min=0),
    numericInput("D", "Read Depth*: ", 30,min=0),
    numericInput("phi","Overdispersion parameter phi (optional):",min=0,max=1,value = ''),
    helpText("To generate a power curve, select experimental parameter to vary"),
    selectInput('tovary',label="Parameter to Vary", choices=c('None','Read Depth','Read length','Window size','Length of CNV'))
  ),
  
  # Show a summary of the dataset and an HTML table with the requested
  # number of observations
  mainPanel(
    tabsetPanel(type = "tabs",
                tabPanel("Input Parameters", verbatimTextOutput("toPrint")),
                tabPanel("Power", verbatimTextOutput("power")),
                tabPanel("Power Curve",verbatimTextOutput("powerCurve_text"), plotOutput("powerCurve")),
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
      sprintf(" L: %s \n D: %s \n W: %s \n N: %s \n l: %s \n" ,input$L, input$D, input$W, input$N, input$l)
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
  output$powerCurve <- renderPlot({
    ps = c()
    if(input$tovary == 'None'){
      print('Select a parameter to vary to generate a power curve.')
    }
    else if(input$tovary=='Read Depth'){
      Ds = c(0.1,1,10,20,30,50,60,100)
      for(d in c(0.1,1,10,20,30,50,60,100)){
        if(input$dis == "pois"){
         if(input$F == 1){ 
          ps = c(ps,poisson_pure(as.numeric(input$alpha), as.numeric(input$W), as.numeric(input$l), d,as.numeric(input$N), as.numeric(input$L)))
         }else{
           ps = c(ps,poisson_impurity(as.numeric(input$alpha), as.numeric(input$W), as.numeric(input$l), d,as.numeric(input$N), as.numeric(input$L),input$F)) 
         }
        }else{
          if(is.na(input$phi)){phi = 1/(input$W*input$D/input$l);}
          else{phi = input$phi;}
          if(input$F == 1){
            ps = c(ps,negativebinomial_purity(input$alpha,input$W,input$l,input$D,input$N,input$L,phi,input$F))
          }else{
            ps = c(ps, negativebinomial_purity(input$alpha,input$W,input$l,input$D,input$N,input$L,phi,input$F))
          }
        }  
      }
      plot(Ds, ps,type="l",xlab="Read Depth",ylab="Power")
      points(Ds, ps,pch=19);
    }
    else if(input$tovary == 'Read length'){
      ls = c(50,75,100,150)
      for(l in ls){
        ps = c(ps,poisson_pure(as.numeric(input$alpha), as.numeric(input$W), l, as.numeric(input$D),as.numeric(input$N), as.numeric(input$L)))
      }
      plot(ls, ps,type="l",xlab="Read Length",ylab="Power")
      points(ls, ps,pch=19);
      
    }

  })
  output$powerCurve_text <- renderText({
    if(input$tovary == 'None'){
      sprintf("Select parameter to vary to plot a power curve")
    }
    else{
      donothing = 1;
    }
  })  
})

shinyApp(ui, server)
