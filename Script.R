# Load required libraries
library(shiny)
library(readxl)
library(dplyr)
library(neuralnet)
library(ggplot2)

ui <- fluidPage(
  titlePanel("Additive Manufacturing Cost Estimator (Pham & Wang Model)"),
  selectInput("input_mode", "Input Mode:", choices = c("Manual Input", "Excel Input")),
  
  conditionalPanel(
    condition = "input.input_mode == 'Excel Input'",
    tabsetPanel(
      tabPanel("Excel Upload",
               fileInput("file", "Upload Excel File (.xlsx)", accept = ".xlsx"),
               actionButton("calc_excel", "Calculate Excel Batch"),
               tableOutput("excel_table")
      )
    )
  ),
  
  conditionalPanel(
    condition = "input.input_mode == 'Manual Input'",
    tabsetPanel(
      tabPanel("Details of part",
               numericInput("n_parts","Number of parts:",1),
               numericInput("V_STL","Effective STL volume (mm³):",888867),
               numericInput("infill","Infill density (0 to 1):",0.68,min=0,max=1,step=0.01),
               numericInput("H","Height (mm):",40),
               numericInput("Wp","Build width Wp (mm):",100)
      ),
      tabPanel("Machine Parameters",
               numericInput("V_s","Scan speed (mm/s):",1257),
               numericInput("V_j","Jump speed (mm/s):",8623),
               numericInput("HS","Scan spacing HS (mm):",0.08),
               numericInput("layer_thk","Layer thickness (mm):",0.1),
               numericInput("V_R","Recoater speed (mm/s):",177),
               numericInput("L_R","Recoater travel distance (mm):",1212),
               numericInput("T_ldelay","T_laser delay (µs):",4596),
               numericInput("T_jdelay","T_jump delay (µs):",500),
               numericInput("alpha","Alpha (α):",1.21),
               numericInput("T_warm","Warm-up time (h):",0),
               numericInput("T_cool","Cool-down time (h):",0),
               numericInput("T_x","Inter-layer delay (s):",5),
               numericInput("T_w","Sinter idle time per layer (s):",0)
      ),
      tabPanel("Material",
               numericInput("mat_cost","Material cost (USD/kg):",60),
               numericInput("mat_density","Material density (g/cm³):",1.2),
               numericInput("support_ratio","Support material factor (0 to 1):",0.2),
               numericInput("recycle","Recycling rate (0 to 1):",0.0)
      ),
      tabPanel("Labor & Overhead",
               numericInput("labor_rate","Labor rate (USD/h):",1.17),
               numericInput("overhead","Overhead (USD/h):",4.60)
      ),
      tabPanel("Post-processing",
               numericInput("post_mat_qty","Post-process material (g):",0),
               numericInput("post_mat_cost","Post-process material cost (USD/kg):",50),
               numericInput("post_labor_time","Post-process labor time (h):",1),
               numericInput("post_labor_rate","Post-process labor rate (USD/h):",1.17)
      ),
      tabPanel("Fuzzy Parameters",
               numericInput("uncertainty","Fuzzy uncertainty (%):",10),
               numericInput("target_cost","Target cost (USD):",20)
      ),
      tabPanel("Results",
               actionButton("calc","Calculate"),
               verbatimTextOutput("results")
      ),
      tabPanel("ML Estimation",
               # Now three separate sliders for each hidden layer
               sliderInput("n_neurons1", "Neurons in hidden layer 1:", min = 1, max = 10, value = 5),
               sliderInput("n_neurons2", "Neurons in hidden layer 2:", min = 1, max = 10, value = 5),
               sliderInput("n_neurons3", "Neurons in hidden layer 3:", min = 1, max = 10, value = 5),
               actionButton("train","Train MLP"),
               verbatimTextOutput("metrics"),
               plotOutput("histogram"),
               plotOutput("scatter")
      )
    )
  )
)

server <- function(input, output, session) {
  observeEvent(input$calc, {
    V_STL   <- input$V_STL; H <- input$H; Wp <- input$Wp
    n_parts <- input$n_parts; rho <- input$infill; thk <- input$layer_thk
    V_s     <- input$V_s; V_j <- input$V_j; HS <- input$HS
    alpha   <- input$alpha
    T_ld    <- input$T_ldelay*1e-6; T_jd <- input$T_jdelay*1e-6
    L_R     <- input$L_R; V_R <- input$V_R
    T_x     <- input$T_x; T_w <- input$T_w
    mat_cost <- input$mat_cost; mat_rho <- input$mat_density
    support  <- input$support_ratio; recycle <- input$recycle
    rate     <- input$labor_rate
    T_warm   <- input$T_warm; T_cool <- input$T_cool
    overhead <- input$overhead
    post_g   <- input$post_mat_qty; post_mat <- input$post_mat_cost
    post_t   <- input$post_labor_time; post_r <- input$post_labor_rate
    delta    <- input$uncertainty/100; C_star <- input$target_cost
    
    # Build Time Calculations
    
    f_rho   <- rho*exp(alpha*(1-rho))
    S       <- V_STL*f_rho/H
    Td      <- (Wp/HS)*(4*T_ld+T_jd)
    V_avg   <- V_s*f_rho+V_j*(1-f_rho)
    N       <- ceiling(H/thk)
    t_pow   <- N*(L_R/V_R+T_x)
    t_sint  <- N*T_w
    t_scan  <- N*(S/(HS*V_avg)+Td)
    warm_s  <- (T_warm+T_cool)*3600
    t_build <- warm_s+t_pow+t_sint+t_scan
    
    # Cost Calculations
    
    V_real   <- V_STL*rho
    usage    <- V_real*(1-recycle)+V_real*support
    C_mat    <- (usage*mat_rho/1e6)*mat_cost
    C_mach_hr<- rate+overhead
    C_setup  <- (T_warm+T_cool)*rate
    C_post   <- (post_g/1000)*post_mat+post_t*post_r
    C_unit   <- C_mat+C_mach_hr*(t_build/3600)+C_setup+rate*post_t+C_post
    C_total  <- C_unit*n_parts
    
    # Fuzzy Arithmetic
    
    C_tri_lo  <- C_unit*(1-delta)
    C_tri_med <- C_unit; 
    C_tri_hi <- C_unit*(1+delta)
    PI <- if(C_star<=C_tri_lo)1 
    else if(C_star>=C_tri_hi)0 
    else 1-(C_star-C_tri_lo)/(C_tri_hi-C_tri_lo)
    mu <- 0.9; 
    entropy <- -mu*log2(mu)-(1-mu)*log2(1-mu); 
    beta <- (1-entropy)*PI
    
    output$results <- renderText({
      paste0(sprintf("Build time: %.2f h",t_build/3600),"\n",
             sprintf("Unit cost: %.2f USD",C_unit),"\n",
             sprintf("Total batch cost: %.2f USD",C_total),"\n",
             sprintf("Cost interval: (%.2f,%.2f,%.2f)",C_tri_lo,C_tri_med,C_tri_hi),"\n",
             sprintf("Possibility (PI): %.3f",PI),"\n",
             sprintf("Entropy: %.3f",entropy),"\n",
             sprintf("Confidence (β): %.3f",beta))
    })
    
  })
  
  observeEvent(input$calc_excel, {
    req(input$file)
    df_in <- read_excel(input$file$datapath)
    df_out <- df_in %>% rowwise() %>% mutate(
      T_ld    = T_ldelay*1e-6,
      T_jd    = T_jdelay*1e-6,
      f_rho   = infill*exp(alpha*(1-infill)),
      S       = V_STL*f_rho/H,
      Td      = (Wp/HS)*(4*T_ld+T_jd),
      V_avg   = V_s*f_rho+V_j*(1-f_rho),
      N       = ceiling(H/layer_thk),
      t_pow   = N*(L_R/V_R+T_x),
      t_sint  = N*T_w,
      t_scan  = N*(S/(HS*V_avg)+Td),
      warm_s  = (T_warm+T_cool)*3600,
      t_build = warm_s+t_pow+t_sint+t_scan,
      usage   = V_STL*infill*(1-recycle)+V_STL*infill*support_ratio,
      C_mat   = (usage*mat_density/1e6)*mat_cost,
      C_mach_hr = labor_rate+overhead,
      C_setup = (T_warm+T_cool)*labor_rate,
      C_post  = (post_mat_qty/1000)*post_mat_cost+post_labor_time*post_labor_rate,
      C_unit  = C_mat+C_mach_hr*(t_build/3600)+C_setup+labor_rate*post_labor_time+C_post,
      C_total = C_unit*n_parts,
      C_tri_lo = C_unit*(1-uncertainty/100), C_tri_med=C_unit, C_tri_hi=C_unit*(1+uncertainty/100),
      PI       = if(target_cost<=C_tri_lo)1 else if(target_cost>=C_tri_hi)0 else 1-(target_cost-C_tri_lo)/(C_tri_hi-C_tri_lo),
      mu       = 0.9, entropy=-mu*log2(mu)-(1-mu)*log2(1-mu), beta=(1-entropy)*PI
    ) %>% ungroup()
    output$excel_table <- renderTable(df_out, digits=3)
  })
  
  
  observeEvent(input$train, {
    set.seed(123)
    n     <- 1000
    rho   <- runif(n, 0.2, 1)
    V_s   <- runif(n, 1000, 1500)
    V_j   <- runif(n, 8000, 9000)
    HS    <- runif(n, 0.05, 0.15)
    H     <- runif(n, 10, 50)
    alpha <- 1.21
    T_ld  <- 4596 * 1e-6
    T_jd  <- 500  * 1e-6
    V_STL <- 150000
    layer_thk <- 0.1
    
    f_rho  <- rho * exp(alpha * (1 - rho))
    V_avg  <- V_s * f_rho + V_j * (1 - f_rho)
    S      <- V_STL * f_rho / H
    Td     <- (100 / HS) * (4 * T_ld + T_jd)
    N      <- ceiling(H / layer_thk)
    t_build<- N * (S / (HS * V_avg) + Td)
    sd_n   <- sd(t_build) * 0.2
    y      <- t_build + rnorm(n, 0, sd_n)
    
    df_raw  <- data.frame(rho, V_s, V_j, HS, H, y)
    df_norm <- as.data.frame(lapply(df_raw, function(x) (x - min(x)) / (max(x) - min(x))))
    
    # Split 80/20
    
    index  <- sample(1:nrow(df_raw), round(0.8 * nrow(df_raw)))
    train_ <- df_norm[index, ]
    test_  <- df_norm[-index, ]
    
    # Train with three independent layer sizes
    
    nn <- neuralnet(
      y ~ rho + V_s + V_j + HS + H,
      data          = train_,
      hidden        = c(input$n_neurons1, input$n_neurons2, input$n_neurons3),
      linear.output = TRUE
    )
    plot(nn)
    
    # Evaluate on test set
    pred_n <- predict(nn, newdata = test_)
    pred   <- pred_n * (max(df_raw$y) - min(df_raw$y)) + min(df_raw$y)
    actual <- df_raw$y[-index]
    err    <- actual - pred
    mse    <- mean(err^2)
    r2     <- 1 - sum(err^2)/sum((actual - mean(actual))^2)
    
    output$metrics <- renderPrint({
      cat(sprintf("Test MSE: %.2f\nTest R²: %.4f", mse, r2))
    })
    output$histogram <- renderPlot({
      hist(err, breaks = 30, main = "Error distribution (test)", xlab = "Error (s)")
    })
    output$scatter <- renderPlot({
      plot(actual, pred, pch = 19,
           xlab = "True build time (s)", ylab = "Predicted build time (s)")
      abline(0, 1, col = "red")
    })
  })
}

shinyApp(ui, server)
