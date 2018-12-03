library(shiny)

proteins <- c("HYDROLASE","TRANSFERASE","OXIDOREDUCTASE","LYASE")

shinyUI(fluidPage(
  
  # Application title
  titlePanel("Protein Database Information for Common Enzymes"),

  sidebarLayout(
    sidebarPanel(
      uiOutput('equation'),
      selectizeInput("classification", "Classification", selected = "HYDROLASE", choices = proteins),
      numericInput("MW", "Maximum Molecular Weight", value = 200000),
      uiOutput("means"),
      br(),
      numericInput("MR","Maximum Allowable Resolution",value = 5),
      uiOutput("text9")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Histogram",
                 downloadButton("downloadPlot1", "Download Histogram"),
                 plotOutput("plot1")),
        tabPanel("Information", HTML("<strong>Thanks for stopping by!</strong>"), 
                 uiOutput("infoTab"),uiOutput("tab")),
        tabPanel("Residues and Resolution Scatter Plot",
                 downloadButton("downloadPlot2", "Download Scatter"),
                 plotOutput("plot2", click = "plot_click"),
                 verbatimTextOutput("info"),
                 verbatimTextOutput("click_info")),
        tabPanel("Specified Data Table", downloadButton("downloadData", "Download specified Data"),
                 tableOutput("table")),
        tabPanel("Point Info", downloadButton("downloadClick", "Download Point Data"),
                 tableOutput("table2")),
        tabPanel("Supervised Learning -- Tree", 
                 selectizeInput("treetype", "treetype", choices = c("Classification - Protein Class", "Regression - Matthews","Regression - Residue Count")),
                 sliderInput("minsplit","minsplit", min = 1, max = 30, value = 1, step = 1),
                 plotOutput("plot4")),
        tabPanel("Simple Linear Regression", plotOutput("plot3"))
      )
      
    )
  )
))
#save