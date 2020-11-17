
library("adaptMT")
library("splines")
library(MASS)
source("utils.R")
source("algo.R")
source("plot.R")
library(shiny)

# Define UI for app that draws a histogram ----
ui <- fluidPage(
    
    # App title ----
    titlePanel(em("AdaPT Simulations")),
    
    # Sidebar layout with input and output definitions ----
    sidebarLayout(position = 'left',
                  
                  # Sidebar panel for inputs ----
                  sidebarPanel(
                      helpText("Create a histogram from normal distribution."),
                      sliderInput("N", "number of samples:", min = 100, max = 5000, value = 1000),
                      numericInput("nu1", "nu1", value=2),
                      numericInput("f_param1", "f_param1", value=2), numericInput("f_param2", "f_param2", value=2), 
                      sliderInput("f_param3", "f_param3", min = 0.0, max = 1.0, value = 0.9), 
                      sliderInput("f_param4", "f_param4", min = 0.0, max = 1.0, value = 0.1),
                      sliderInput("rho", "mean:", min = 0, max = 1, value = 0),
                      selectInput("opt", "cov:", choices = list("AR(1)" = 'ar1', "CS" = 'cs'), selected = 1),
                      sliderInput("alpha", "alpha:", min = 0.01, max = 0.30, value = 0.01)
                  ),
                  
                  # Main panel for displaying outputs ----
                  mainPanel(
                      textOutput("selected"),
                      plotOutput("hist")
                  )
    )
)

# Define server logic required to draw a histogram ----
server <- function(input, output) {
    
    output$selected <- renderText({ paste("You are plotting a histogrm of ", input$n," normal samples with mean ", input$mean," and sd ", input$sd) })
    output$hist <- renderPlot({ 
        # generate data
        data <- generate_data(input$N, c(0, input$nu1), 
                              c(input$f_param1, input$f_param2, input$f_param3, input$f_param4))
        
        # create dependence
        Cov <- cor_mat(sum(data$H==1), 0.1, 'ar1')
        
        # generate p value
        data$z <- rnorm(input$N, data$nu)
        data$z[data$H==1] <- mvrnorm(1, data$nu[data$H==1], Cov)
        data$pvals <- 1 - pnorm(data$z)
        
        # run algorithms
        df_BH <- summary_BH(data$pvals, data$H, alphas = seq(0.01, 0.3, 0.01))
        df_storey <- summary_storey(data$pvals, data$H, alphas = seq(0.01, 0.3, 0.01))
        
        formulas <- paste0("ns(x, df = ", 6:7, ")")
        adapt <- adapt_glm(x = data.frame(x = data$x), pvals = data$pvals, pi_formulas = formulas,
                           mu_formulas = formulas,  nfits = 10, alphas = seq(0.01, 0.3, 0.01))
        df_adapt <- summary_adapt(adapt, data$pvals, data$H)
        
        # plot
        plot_s_curve(adapt, data$x, data$pvals, input$alpha)
    })
    
}

shinyApp(ui = ui, server = server)

