library("adaptMT")
library("splines")
library(MASS)
source("utils.R")
source("algo.R")
source("plot.R")
library(shiny)

# Define UI for app
ui <- fluidPage(
  navbarPage("AdaPT simulation",
      tabPanel("home",
            fluidRow(
              column(
                
              br(),
              p("We consider the problem of multiple hypothesis testing with generic side information: 
              for each hypothesis Hi we observe both a p-value pi and some predictor xi encoding contextual in- formation about the hypothesis. 
              For large-scale problems, adaptively focusing power on the more promising hypotheses (those more likely to yield discoveries) 
              can lead to much more powerful multiple testing procedures. We propose a general iterative framework for this problem, 
              called the Adaptive p-value Thresholding (AdaPT) procedure, which adaptively estimates a Bayes-optimal p-value rejection threshold and 
              controls the false discovery rate (FDR) in finite samples.",
              style="text-align:justify;color:black;background-color:lavender;padding:15px;border-radius:10px"),
              br(),
                
              p("At each iteration of the procedure, the analyst proposes a rejection threshold and observes partially cen- sored p-values, 
              estimates the false discovery proportion (FDP) below the threshold, and proposes another threshold, 
              until the estimated FDP is below α. Our procedure is adaptive in an unusu- ally strong sense, permitting the analyst 
              to use any statistical or machine learning method she chooses to estimate the optimal threshold, 
              and to switch between different models at each itera- tion as information accrues.",
              style="text-align:justify;color:black;background-color:papayawhip;padding:15px;border-radius:10px"),
              
              width=8),
              
              column(
                br(),
                br(),
                p("For more information please check the",em("Lihua Li and William Fithian's paper"),"page clicking",
                br(),
                a(href="https://arxiv.org/abs/1609.06035", "Here",target="_blank"),
                style="text-align:center;color:black"),
                
                width=2)),
            
            hr(),
            p(em("Developed by"),br("Bowei Kang, Jin-Hong Du, Xinxin Chen"),style="text-align:center; font-family: times")
      ),
  #titlePanel(em("AdaPT Simulations")),
  
  tabPanel('Low Total Variation',
      sidebarLayout(position = 'left',
            # Sidebar panel for inputs ----
            sidebarPanel(
              helpText("Create a histogram from normal distribution."),
              sliderInput("N_f1", "number of samples:", min = 100, max = 5000, value = 200, step=100),
              numericInput("nu1_f1", "nu1", value=2),
              numericInput("f_param1_f1", "f_param1", value=2), 
              numericInput("f_param2_f1", "f_param2", value=2), 
              sliderInput("f_param3_f1", "f_param3", min = 0.0, max = 1.0, value = 0.9, step=0.05), 
              sliderInput("f_param4_f1", "f_param4", min = 0.0, max = 1.0, value = 0.1, step=0.05),
              sliderInput("rho_f1", "rho:", min = 0, max = 1, value = 0),
              selectInput("opt_f1", "cov:", choices = list("AR(1)" = 'ar1', "CS" = 'cs'), selected = 1),
              sliderInput("sparsity_f1", "sparsity:", min = 0, max = 1, value = 0.0, step=0.1),
              sliderInput("alpha_f1", "alpha:", min = 0.01, max = 0.30, value = 0.05, step=0.01)
            ),
            
            # Main panel for displaying outputs ----
            mainPanel(
                verticalLayout(
                    plotOutput("p1_f1"), 
                    plotOutput("p2_f1")
                                     )
            )
        )
      ),
  
  tabPanel('Ordered structure',
      sidebarLayout(position = 'left',
                   
           # Sidebar panel for inputs ----
           sidebarPanel(
             sliderInput("N_f2", "number of samples:", min = 100, max = 5000, value = 300, step=100),
             numericInput("nu1_f2", "nu1", value=2),
             numericInput("f_param1_f2", "f_param1", value=-2), 
             numericInput("f_param2_f2", "f_param2", value=-2), 
             sliderInput("rho_f2", "rho:", min = 0, max = 1, value = 0),
             selectInput("opt_f2", "cov:", choices = list("AR(1)" = 'ar1', "CS" = 'cs'), selected = 1),
             sliderInput("sparsity_f2", "sparsity:", min = 0, max = 1, value = 0.0, step=0.1),
             sliderInput("alpha_f2", "alpha:", min = 0.01, max = 0.30, value = 0.05, step=0.01)
           ),
           
           # Main panel for displaying outputs ----
           mainPanel(
             verticalLayout(
               plotOutput("p1_f2"), 
               plotOutput("p2_f2")
             )
           )
    )
    ),
  
 tabPanel('Group structure',
     sidebarLayout(position = 'left',
          
          # Sidebar panel for inputs ----
          sidebarPanel(
            sliderInput("N_f3", "number of samples:", min = 100, max = 5000, value = 300, step=100),
            numericInput("nu1_f3", "nu1", value=2),
            sliderInput("f_param1_f3", "f_param1", min = 2, max = 10, value = 1, step=1), 
            sliderInput("rho_f3", "rho:", min = 0, max = 1, value = 0),
            selectInput("opt_f3", "cov:", choices = list("AR(1)" = 'ar1', "CS" = 'cs'), selected = 1),
            sliderInput("sparsity_f3", "sparsity:", min = 0, max = 1, value = 0.0, step=0.1),
            sliderInput("alpha_f3", "alpha:", min = 0.01, max = 0.30, value = 0.05, step=0.01)
          ),
          
          # Main panel for displaying outputs ----
          mainPanel(
            verticalLayout(
              plotOutput("p1_f3"), 
              plotOutput("p2_f3")
            )
          )
    )
 )
  ) # Navigate
) # fluidPage


# Define server logic required to draw a histogram ----
server <- function(input, output) {
    dataInput_f1 <- reactive({
        # generate data
        data <- generate_data(input$N_f1, c(0, input$nu1_f1), 
                              c(input$f_param1_f1, input$f_param2_f1, input$f_param3_f1, input$f_param4_f1),
                              1)
        
        # create dependence
        Cov <- cor_mat(sum(data$H==1), input$rho_f1, input$opt_f1, input$sparsity_f1)
        
        # generate p value
        data$z <- rnorm(input$N_f1, data$nu)
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
        list(data, df_BH, df_storey, df_adapt, adapt)
    })
    
    dataInput_f2 <- reactive({
      # generate data
      data <- generate_data(input$N_f2, c(0, input$nu1_f2), 
                            c(input$f_param1_f2, input$f_param2_f2),
                            2)
      
      # create dependence
      Cov <- cor_mat(sum(data$H==1), input$rho_f2, input$opt_f2, input$sparsity_f2)
      
      # generate p value
      data$z <- rnorm(input$N_f2, data$nu)
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
      list(data, df_BH, df_storey, df_adapt, adapt)
    })
  
    dataInput_f3 <- reactive({
      # generate data
      data <- generate_data(input$N_f3, c(0, input$nu1_f3), 
                            c(input$f_param1_f3),
                            3)
      
      # create dependence
      Cov <- cor_mat(sum(data$H==1), input$rho_f3, input$opt_f3, input$sparsity_f3)
      
      # generate p value
      data$z <- rnorm(input$N_f3, data$nu)
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
      list(data, df_BH, df_storey, df_adapt, adapt)
    })
    
    output$p1_f1 <- renderPlot({
      l = dataInput_f1()
      data = l[[1]]
      df_BH = l[[2]]
      df_storey = l[[3]]
      adapt = l[[5]]
      alphas <- seq(0.01, 0.3, 0.01)
      plot_s_curve(adapt, data$x, data$pvals, input$alpha_f1, data$H,
                   df_BH[abs(alphas-input$alpha_f1)<1e-12,'alpha'], df_storey[abs(alphas-input$alpha_f1)<1e-12,'alpha'])
    })
    
    output$p2_f1 <- renderPlot({
      l = dataInput_f1()
      data = l[[1]]
      df_BH = l[[2]]
      df_storey = l[[3]]
      df_adapt = l[[4]]
      adapt = l[[5]]
      alphas <- seq(0.01, 0.3, 0.01)
      plot_power(alphas, df_BH, df_storey, df_adapt)
    })
    
    
    output$p1_f2 <- renderPlot({
        l = dataInput_f2()
        data = l[[1]]
        df_BH = l[[2]]
        df_storey = l[[3]]
        adapt = l[[5]]
        alphas <- seq(0.01, 0.3, 0.01)
        plot_s_curve(adapt, data$x, data$pvals, input$alpha_f2, data$H,
                     df_BH[abs(alphas-input$alpha_f2)<1e-12,'alpha'], df_storey[abs(alphas-input$alpha_f2)<1e-12,'alpha'])
    })
  
    output$p2_f2 <- renderPlot({
        l = dataInput_f2()
        data = l[[1]]
        df_BH = l[[2]]
        df_storey = l[[3]]
        df_adapt = l[[4]]
        adapt = l[[5]]
        alphas <- seq(0.01, 0.3, 0.01)
        plot_power(alphas, df_BH, df_storey, df_adapt)
    })
    
    output$p1_f3 <- renderPlot({
      l = dataInput_f3()
      data = l[[1]]
      df_BH = l[[2]]
      df_storey = l[[3]]
      adapt = l[[5]]
      alphas <- seq(0.01, 0.3, 0.01)
      plot_s_curve(adapt, data$x, data$pvals, input$alpha_f3, data$H,
                   df_BH[abs(alphas-input$alpha_f3)<1e-12,'alpha'], df_storey[abs(alphas-input$alpha_f3)<1e-12,'alpha'])
    })
    
    output$p2_f3 <- renderPlot({
      l = dataInput_f3()
      data = l[[1]]
      df_BH = l[[2]]
      df_storey = l[[3]]
      df_adapt = l[[4]]
      adapt = l[[5]]
      alphas <- seq(0.01, 0.3, 0.01)
      plot_power(alphas, df_BH, df_storey, df_adapt)
    })
  
}

shinyApp(ui = ui, server = server)

