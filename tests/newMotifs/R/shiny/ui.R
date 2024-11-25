library(shiny)
library(shinyjs)
library(bslib)

ui<-
  #useShinyjs(),
  #responsive = FALSE,
  page_navbar(title = "Network motifs in transcription networks",
              id = "motif",
    nav_panel("NAR", value = 1,
      card(
        card_header("Negative autoregulation"),
               "Z gets activated when X is on and has expression reduced by Z concentration. Activity of X shown as
               blue shading."),
      fluidRow(
        column(width = 6, 
               sliderInput(inputId = "tmax_NAR", label="Simulation time", min = 1, 
                           max = 50, value = 10, step = 1))),
      fluidRow(  
        column(width = 6,
               sliderInput(inputId = "Xstart_NAR",label = "Signal starts at:", 
                           min = 1, max = 50, value = 1, step = 1)),
        column(width = 6,
               sliderInput(inputId = "Xstop_NAR", label = "Signal stops at:", 
                           min = 1, max = 50, value = 6, step = 1))),
      fluidRow(
        column(width = 4,
               sliderInput(inputId = "base_NAR", label = "Constitutive promoter expression:",
                           min = 0, max = 5, value = 0, step = 0.1)),
        column(width = 4,
               sliderInput(inputId = "aZ_NAR", label = "Z degradation rate:",
                           min = 0, max = 5, value = 1, step = 0.1)),
        column(width = 4,
               sliderInput(inputId = "bZ_NAR", label = "Z production rate:",
                           min = 0, max = 5, value = 1, step = 0.1))),
      fluidRow(
        column(width = 4,
               sliderInput(inputId = "KXZ_NAR", label = "X -> Z regulation strength:",
                           min = 0, max = 5, value = 1, step = 0.1)),
        column(width = 4,
               sliderInput(inputId = "KZ_NAR", label = "Z autoregulation strength:",
                           min = 0, max = 5, value = 1, step = 0.1)),
        column(width = 4,
               sliderInput(inputId = "Hilln_NAR", label = "Hill coefficient:",
                           min = 0, max = 8, value = 1, step = 0.1))),
      fluidRow(
        column(width = 3,
               sliderInput(inputId = "XMult_NAR", label = "Environmental signal strength:",
                           min = 0, max = 5, value = 1, step = 0.1))
      ),
      fluidRow(
        hr(),
        column(width = 12, align = "center",
               plotOutput(outputId = "main_plot_NAR", height = "800px", width = "600px"))
      )   
    ),
    
    nav_panel("PAR", value = 2,
      card(
        card_header("Positive autoregulation"),
        "Z gets activated when X is on and is promoted further by Z concentration. Activity of X shown as
        blue shading."),
      fluidRow(
        column(width = 6, 
               sliderInput(inputId = "tmax_PAR", label="Simulation time", min = 1, 
                           max = 50, value = 10, step = 1))),
      fluidRow(  
        column(width = 6,
               sliderInput(inputId = "Xstart_PAR",label = "Signal starts at:", 
                           min = 1, max = 50, value = 1, step = 1)),
        column(width = 6,
               sliderInput(inputId = "Xstop_PAR", label = "Signal stops at:", 
                           min = 1, max = 50, value = 6, step = 1))),
      fluidRow(
        column(width = 4,
               sliderInput(inputId = "base_PAR", label = "Constitutive promoter expression:",
                           min = 0, max = 5, value = 0.1, step = 0.1)),
        column(width = 4,
               sliderInput(inputId = "aZ_PAR", label = "Z degradation rate:",
                           min = 0, max = 5, value = 1, step = 0.1)),
        column(width = 4,
               sliderInput(inputId = "bZ_PAR", label = "Z production rate:",
                           min = 0, max = 5, value = 1, step = 0.1))),
      fluidRow(
        column(width = 4,
               sliderInput(inputId = "KXZ_PAR", label = "X -> Z regulation strength:",
                           min = 0, max = 5, value = 1, step = 0.1)),
        column(width = 4,
               sliderInput(inputId = "KZ_PAR", label = "Z autoregulation strength:",
                           min = 0, max = 5, value = 1, step = 0.1)),
        column(width = 4,
               sliderInput(inputId = "Hilln_PAR", label = "Hill coefficient:",
                           min = 0, max = 8, value = 1, step = 0.1))),
      fluidRow(
        column(width = 3,
               sliderInput(inputId = "XMult_PAR", label = "Environmental signal strength:",
                           min = 0, max = 5, value = 1, step = 0.1))
      ),
      fluidRow(
        hr(),
        column(width = 12, align = "center",
               plotOutput(outputId = "main_plot_PAR", height = "800px", width = "600px"))
      )
    ),
      nav_panel("FFL-C1", value = 3,
                card(
                  card_header("Feedforward loop (Coherent Type I)"),
                  "Z is activated by X AND Y. Activity of X shown as
               blue shading."),
               
               fluidRow(
                 column(width = 6, 
                        sliderInput(inputId = "tmax_FFLC1", label="Simulation time", min = 1, 
                                    max = 50, value = 10, step = 1))),
               fluidRow(  
                 column(width = 6,
                        sliderInput(inputId = "Xstart_FFLC1",label = "Signal starts at:", 
                                    min = 1, max = 50, value = 1, step = 1)),
                 column(width = 6,
                        sliderInput(inputId = "Xstop_FFLC1", label = "Signal stops at:", 
                                    min = 1, max = 50, value = 6, step = 1))),
               fluidRow(
                 column(width = 4,
                        sliderInput(inputId = "base_FFLC1", label = "Constitutive promoter expression:",
                                    min = 0, max = 5, value = 0, step = 0.1)),
                 column(width = 4,
                        sliderInput(inputId = "aY_FFLC1", label = "Y degradation rate:",
                                    min = 0, max = 5, value = 1, step = 0.1)),
                 column(width = 4,
                        sliderInput(inputId = "bY_FFLC1", label = "Y production rate:",
                                    min = 0, max = 5, value = 1, step = 0.1))),
               fluidRow(
                 column(width = 4,
                        sliderInput(inputId = "KY_FFLC1", label = "X -> Y regulation strength:",
                                    min = 0, max = 5, value = 1, step = 0.1)),
                 column(width = 4,
                        sliderInput(inputId = "KXZ_FFLC1", label = "X -> Z regulation strength:",
                                    min = 0, max = 5, value = 1, step = 0.1)),
                column(width = 4,
                       sliderInput(inputId = "aZ_FFLC1", label = "Z degradation rate:",
                                   min = 0, max = 5, value = 1, step = 0.1))),
              fluidRow(
                column(width = 4,
                       sliderInput(inputId = "bZ_FFLC1", label = "Z production rate:",
                                   min = 0, max = 5, value = 1, step = 0.1)),
                 column(width = 4,
                        sliderInput(inputId = "Hilln_FFLC1", label = "Hill coefficient:",
                                    min = 0, max = 8, value = 1, step = 0.1)),
                 column(width = 4,
                        sliderInput(inputId = "XMult_FFLC1", label = "Environmental signal strength:",
                                    min = 0, max = 5, value = 1, step = 0.1))
               ),
               fluidRow(
                 hr(),
                 column(width = 12, align = "center",
                        plotOutput(outputId = "main_plot_FFLC1", height = "800px", width = "600px"))
               ),   
      ),
    nav_panel("FFL-I1", value = 4,
              card(
                card_header("Feedforward loop (Incoherent Type I)"),
                "Z is activated by X AND deactivated by Y. Activity of X shown as
               blue shading."),
              
              fluidRow(
                column(width = 6, 
                       sliderInput(inputId = "tmax_FFLI1", label="Simulation time", min = 1, 
                                   max = 50, value = 10, step = 1))),
              fluidRow(  
                column(width = 6,
                       sliderInput(inputId = "Xstart_FFLI1",label = "Signal starts at:", 
                                   min = 1, max = 50, value = 1, step = 1)),
                column(width = 6,
                       sliderInput(inputId = "Xstop_FFLI1", label = "Signal stops at:", 
                                   min = 1, max = 50, value = 6, step = 1))),
              fluidRow(
                column(width = 4,
                       sliderInput(inputId = "base_FFLI1", label = "Constitutive promoter expression:",
                                   min = 0, max = 5, value = 0.1, step = 0.1)),
                column(width = 4,
                       sliderInput(inputId = "aY_FFLI1", label = "Y degradation rate:",
                                   min = 0, max = 5, value = 1, step = 0.1)),
                column(width = 4,
                       sliderInput(inputId = "bY_FFLI1", label = "Y production rate:",
                                   min = 0, max = 5, value = 1, step = 0.1))),
              fluidRow(
                column(width = 4,
                       sliderInput(inputId = "KY_FFLI1", label = "X -> Y regulation strength:",
                                   min = 0, max = 5, value = 1, step = 0.1)),
                column(width = 4,
                       sliderInput(inputId = "KXZ_FFLI1", label = "X -> Z regulation strength:",
                                   min = 0, max = 5, value = 1, step = 0.1)),
                column(width = 4,
                       sliderInput(inputId = "aZ_FFLI1", label = "Z degradation rate:",
                                   min = 0, max = 5, value = 1, step = 0.1))),
              fluidRow(
                column(width = 4,
                       sliderInput(inputId = "bZ_FFLI1", label = "Z production rate:",
                                   min = 0, max = 5, value = 1, step = 0.1)),
                column(width = 4,
                       sliderInput(inputId = "Hilln_FFLI1", label = "Hill coefficient:",
                                   min = 0, max = 8, value = 1, step = 0.1)),
                column(width = 4,
                       sliderInput(inputId = "XMult_FFLI1", label = "Environmental signal strength:",
                                   min = 0, max = 5, value = 1, step = 0.1))
              ),
              fluidRow(
                hr(),
                column(width = 12, align = "center",
                       plotOutput(outputId = "main_plot_FFLI1", height = "800px", width = "600px"))
              )
    ),
    nav_panel("FFBH", value = 5,
              card(
                card_header("Feedforward/back hybrid"),
                "Z is activated by X AND Y. X is activated by an environmental signal and Z. Activity of X shown as
               blue shading."),
      
              fluidRow(
                column(width = 6, 
                       sliderInput(inputId = "tmax_FFBH", label="Simulation time", min = 1, 
                                   max = 50, value = 10, step = 1))),
              fluidRow(  
                column(width = 6,
                       sliderInput(inputId = "Xstart_FFBH",label = "Signal starts at:", 
                                   min = 1, max = 50, value = 1, step = 1)),
                column(width = 6,
                       sliderInput(inputId = "Xstop_FFBH", label = "Signal stops at:", 
                                   min = 1, max = 50, value = 6, step = 1))),

              fluidRow(
                column(width = 4,
                       sliderInput(inputId = "base_FFBH", label = "Constitutive promoter expression:",
                                   min = 0, max = 5, value = 0, step = 0.1)),
                column(width = 4,
                       sliderInput(inputId = "aX_FFBH", label = "X degradation rate:",
                                   min = 0, max = 5, value = 1, step = 0.1)),
                column(width = 4,
                       sliderInput(inputId = "KZX_FFBH", label = "Z -> X regulation strength:",
                                   min = 0, max = 5, value = 1, step = 0.1))),
              fluidRow(
                column(width = 4,
                       sliderInput(inputId = "aY_FFBH", label = "Y degradation rate:",
                                   min = 0, max = 5, value = 1, step = 0.1)),
                column(width = 4,
                       sliderInput(inputId = "bY_FFBH", label = "Y production rate:",
                                   min = 0, max = 5, value = 1, step = 0.1)),
                column(width = 4,
                       sliderInput(inputId = "KY_FFBH", label = "X -> Y regulation strength:",
                                   min = 0, max = 5, value = 1, step = 0.1))),
              fluidRow(
                column(width = 4,
                       sliderInput(inputId = "KXZ_FFBH", label = "X -> Z regulation strength:",
                                   min = 0, max = 5, value = 1, step = 0.1)),
                column(width = 4,
                       sliderInput(inputId = "aZ_FFBH", label = "Z degradation rate:",
                                   min = 0, max = 5, value = 1, step = 0.1)),
                column(width = 4,
                       sliderInput(inputId = "bZ_FFBH", label = "Z production rate:",
                                   min = 0, max = 5, value = 1, step = 0.1))),
              fluidRow(
                column(width = 4,
                       sliderInput(inputId = "Hilln_FFBH", label = "Hill coefficient:",
                                   min = 0, max = 8, value = 1, step = 0.1)),
                column(width = 4,
                       sliderInput(inputId = "XMult_FFBH", label = "Environmental signal strength:",
                                   min = 0, max = 5, value = 1, step = 0.1))
              ),
              fluidRow(
                hr(),
                column(width = 12, align = "center",
                       plotOutput(outputId = "main_plot_FFBH", height = "800px", width = "600px"))
              )
    )
  )
  

