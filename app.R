# Load packages
library("tidyverse")
library("coin")

# Load data
dat <- read_csv("dataset_tktd_exposure.csv",
                na = "NA",
                col_types = cols(uvb_duration_hour = col_factor(),
                                 efv_intensity_µg_L = col_factor(),
                                 exposure_order = col_factor())) %>% 
  mutate(time_hour = time_hour/24) %>% 
  rename(time_day = time_hour)

# Build app
ui <- fluidPage(

  titlePanel(strong(em("Comparison of survival curves using Gehan-Breslow-Wilcoxon method (GBW)"))),
  
  strong(p("Four-factor multistress experiment",
           tags$a(href ="#fn1",
                  id = "#ref1",
                  tags$sup("1")))),
  p("Fully-crossed 2 stress factors including control:",
    br(),
    "1. ESF concentration intensity: 9 + 1 levels, (0.00, 0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64, 1.28, 2.56 μg/L)",
    br(),
    "2.  UVB radiation duration: 5 + 1 levels, (0, 4, 8, 10, 12, 14 hrs)",
    br(),
    "Temporal factors:",
    br(),
    "3.  Order of exposure: 2 levels, (ESF-UVB/UVB-ESF)",
    br(),
    "4.  Periods of exposure: 2 levels, (0, 2 days)"),
  
  sidebarLayout(
    
    sidebarPanel(
      
      selectInput(
        
        "conc",
        label = strong("Exposure concentration of esfenvalerate (ESF) (µg/L)"),
        choices = unique(dat$efv_intensity_µg_L),
        selected = 2.56
        
      ),
      
      selectInput(
        
        "dur",
        label = strong("Exposure duration of UVB radiation (h)"),
        choices = unique(dat$uvb_duration_hour),
        selected = 14
        
      ),
      
      selectInput(
        
        "ord1",
        label = strong("Exposure order 1"),
        choices = unique(dat$exposure_order),
        selected = "E0U"
        
      ),
      
      selectInput(
        
        "ord2",
        label = strong("Exposure order 2"),
        choices = unique(dat$exposure_order),
        selected = "E2U"
        
      )
      
    ),
    
    mainPanel(
    
      plotOutput(outputId = "SurvCurv")

    )
  ),
  
  br(),
  br(),
  br(),
  br(),
  p(tags$sup(id = "fn1", "1"),"[Florian Schunck and Matthias Liess
      Environmental Science & Technology 2022 56 (20), 14660-14667",
      tags$a(href = "https://pubs.acs.org/doi/10.1021/acs.est.2c04345?ref=pdf",
             "DOI: 10.1021/acs.est.2c04345"),
    "]")
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  # Select variables
  pick.var <- reactive({
    
    dat %>% 
      na.omit() %>%
      group_by(exposure_order) %>%
      filter(uvb_duration_hour == input$dur,
             efv_intensity_µg_L == input$conc,
             exposure_order %in% c(input$ord1, input$ord2))
    
  })
    
  output$SurvCurv <- renderPlot({
    
    gbw.test <- logrank_test(Surv(pick.var()$survivors) ~ pick.var()$exposure_order, type = "Gehan-Breslow")
    
    pvalue <- signif(pvalue(gbw.test, method = "Gehan-Breslow"), 2)
    
    ggplot(pick.var()) +
      geom_point(aes(x = time_day,
                     y = survivors,
                     colour = exposure_order),
                 alpha = 0.75,
                 size = 5) +
      geom_line(aes(x = time_day,
                    y = survivors,
                    colour = exposure_order),
                alpha = 0.5,
                linewidth = 2) +
      expand_limits(y = c(0, 4)) +
      scale_x_continuous(breaks = seq(0, 12, 2)) +
      annotate("text",
               x = 8,
               y = 3.5,
               label = paste0("GBW test\n",
                              "p = ",
                              pvalue),
               size = 5) +
      labs(title = paste0(" Combined exposure of ESF = ",
                          input$conc,
                          " µg/L and UVB = ",
                          input$dur,  " h"),
           y = "Number of survivors",
           x = "Time (days)",
           colour = "Exposure order") +
      theme_classic(base_size = 18)
    
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
