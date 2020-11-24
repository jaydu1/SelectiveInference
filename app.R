library("adaptMT")
library("splines")
library(MASS)
source("utils.R")
source("algo.R")
source("plot.R")
library(shiny)

# Define UI for app
ui <- fluidPage(
  titlePanel(em("AdaPT Simulations")),
  
  sidebarLayout(position = 'left',
                
                # Sidebar panel for inputs ----
                sidebarPanel(
                  helpText("Create a histogram from normal distribution."),
                  sliderInput("N", "number of samples:", min = 100, max = 5000, value = 1000, step=100),
                  numericInput("nu1", "nu1", value=2),
                  numericInput("f_param1", "f_param1", value=2), numericInput("f_param2", "f_param2", value=2), 
                  sliderInput("f_param3", "f_param3", min = 0.0, max = 1.0, value = 0.9, step=0.05), 
                  sliderInput("f_param4", "f_param4", min = 0.0, max = 1.0, value = 0.1, step=0.05),
                  sliderInput("rho", "rho:", min = 0, max = 1, value = 0),
                  selectInput("opt", "cov:", choices = list("AR(1)" = 'ar1', "CS" = 'cs'), selected = 1),
                  sliderInput("sparsity", "sparsity:", min = 0, max = 1, value = 0.0, step=0.1),
                  sliderInput("alpha", "alpha:", min = 0.01, max = 0.30, value = 0.05, step=0.01)
                ),
                
                # Main panel for displaying outputs ----
                mainPanel(
                    verticalLayout(
                        plotOutput("p1"), 
                        plotOutput("p2")
                                         )
                )
  )
)

# Define server logic required to draw a histogram ----
server <- function(input, output) {
    dataInput <- reactive({
        # generate data
        data <- generate_data(input$N, c(0, input$nu1), 
                              c(input$f_param1, input$f_param2, input$f_param3, input$f_param4))
        
        # create dependence
        Cov <- cor_mat(sum(data$H==1), input$rho, input$opt, input$sparsity)
        
        # generate p value
        data$z <- rnorm(input$N, data$nu)
        data$z[data$H==1] <- mvrnorm(1, data$nu[data$H==1], Cov)
        data$pvals <- 1 - pnorm(data$z)
        
        # run algorithms
        alphas <- seq(0.01, 0.3, 0.01)
        df_BH <- summary_BH(data$pvals, data$H, alphas = alphas)
        df_storey <- summary_storey(data$pvals, data$H, alphas = alphas)
        
        formulas <- paste0("ns(x, df = ", 6:7, ")")
        adapt <- adapt_glm(x = data.frame(x = data$x), pvals = data$pvals, pi_formulas = formulas,
                           mu_formulas = formulas,  nfits = 10, alphas = alphas,
                           verbose=list(print = FALSE, fit = FALSE, ms = FALSE))
        df_adapt <- summary_adapt(adapt, data$pvals, data$H)
        list(data, df_BH, df_storey, adapt)
    })
  
    # output$selected <- renderText({ paste("You are plotting a histogrm of ", input$n," normal samples with mean ", input$mean," and sd ", input$sd) })
    output$p1 <- renderPlot({
        l = dataInput()
        data = l[[1]]
        df_BH = l[[2]]
        df_storey = l[[3]]
        adapt = l[[4]]
        alphas <- seq(0.01, 0.3, 0.01)
        plot_s_curve(adapt, data$x, data$pvals, input$alpha, data$H,
                     df_BH[abs(alphas-input$alpha)<1e-12,'alpha'], df_storey[abs(alphas-input$alpha)<1e-12,'alpha'])
    })
  
    output$p2 <- renderPlot({
        l = dataInput()
        data = l[[1]]
        df_BH = l[[2]]
        df_storey = l[[3]]
        adapt = l[[4]]
        alphas <- seq(0.01, 0.3, 0.01)
        plot_power(alphas, df_BH, df_storey, df_adapt)
    })
  
}

shinyApp(ui = ui, server = server)

