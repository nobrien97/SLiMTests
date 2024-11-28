library(deSolve)
library(tidyverse)
library(ggraph)
library(igraph)
library(tidygraph)
library(cowplot)
library(shiny)
library(RColorBrewer)
library(shinyjs)

# Colour scheme:
colX   <- brewer.pal(6, "Paired")[2]
colXbg <- brewer.pal(6, "Paired")[1]
colY   <- brewer.pal(6, "Paired")[6]
colYbg <- brewer.pal(6, "Paired")[5]
colZ   <- brewer.pal(6, "Paired")[4]
colZbg <- brewer.pal(6, "Paired")[3]

colX   <- "#44546A"
colXOutline <- "#222A35"


# plot graphs:

genes <- tibble(name = c(expression("X"), expression("Z")))
activationNAR <- tibble(from = c(1, 2),
                        to =   c(2, 2))
gNAR <- tbl_graph(edges = activationNAR, directed = TRUE, nodes = genes)
plt_NAR <- ggraph(gNAR, layout = "manual", x = c(0,0), y=c(1,0)) +
  geom_node_point(size = 16, color = c(colXOutline)) +
  geom_node_point(size = 15, color = c(colX)) +
  geom_node_text(aes(label = name), size = 8, color = "white", parse = T) +
  geom_edge_arc(aes(start_cap = circle(8, unit = "mm"),
                    end_cap = circle(8, unit = "mm")),
                arrow = arrow(type = "open", angle = 30, length = unit(3, 'mm')),
                strength = 0.0) +
  geom_edge_loop(aes(start_cap = circle(8, unit = "mm"),
                     end_cap = circle(8, unit = "mm"),
                     span = 90, direction = 0, strength = 0.7),
                 arrow = arrow(type = "open", angle = 90, length = unit(3, 'mm'))) +
  scale_x_continuous(limits = c(-0.5,0.5)) +
  scale_y_continuous(limits = c(-0.2,1.1)) +
  coord_fixed() +
  theme_graph()

# PAR
plt_PAR <- ggraph(gNAR, layout = "manual", x = c(0,0), y=c(1,0)) +
  geom_node_point(size = 16, color = c(colXOutline)) +
  geom_node_point(size = 15, color = c(colX)) +
  geom_node_text(aes(label = name), size = 8, color = "white", parse = T) +
  geom_edge_arc(aes(start_cap = circle(8, unit = "mm"),
                    end_cap = circle(8, unit = "mm")),
                arrow = arrow(type = "open", angle = 30, length = unit(3, 'mm')),
                strength = 0.0) +
  geom_edge_loop(aes(start_cap = circle(8, unit = "mm"),
                     end_cap = circle(8, unit = "mm"),
                     span = 90, direction = 0, strength = 0.7),
                 arrow = arrow(type = "open", angle = 30, length = unit(3, 'mm'))) +
  scale_x_continuous(limits = c(-0.5,0.5)) +
  scale_y_continuous(limits = c(-0.2,1.1)) +
  coord_fixed() +
  theme_graph()

genes <- tibble(name = c(expression("X"), expression("Y"), expression("Z")))
activationFFL <- tibble(from = c(1, 1, 2),
                        to =   c(2, 3, 3))
gFFL <- tbl_graph(edges = activationFFL, directed = TRUE, nodes = genes)
plt_FFLC1 <- ggraph(gFFL, layout = "manual", x = c(0,-0.7,0.7), y=c(1,0,0)) +
  geom_node_point(size = 16, color = colXOutline) +
  geom_node_point(size = 15, color = colX) +
  geom_node_text(aes(label = name), size = 8, color = "white", parse = T) +
  geom_edge_arc(aes(start_cap = circle(8, unit = "mm"),
                    end_cap = circle(8, unit = "mm")),
                arrow = arrow(type = "open", angle = 30, length = unit(3, 'mm')),
                strength = c(-0.2, 0.2, -0.2)) +
  scale_x_continuous(limits = c(-1.5,1.5)) +
  scale_y_continuous(limits = c(-0.5,1.5)) +
  coord_fixed() +
  theme_graph()

# I1 FFL
plt_FFLI1 <- ggraph(gFFL, layout = "manual", x = c(0,-0.7,0.7), y=c(1,0,0)) +
  geom_node_point(size = 16, color = colXOutline) +
  geom_node_point(size = 15, color = colX) +
  geom_node_text(aes(label = name), size = 8, color = "white", parse = T) +
  geom_edge_arc(aes(start_cap = circle(8, unit = "mm"),
                    end_cap = circle(8, unit = "mm")),
                arrow = arrow(type = "open", angle = c(30, 30, 90), length = unit(3, 'mm')),
                strength = c(-0.2, 0.2, -0.2)) +
  scale_x_continuous(limits = c(-1.5,1.5)) +
  scale_y_continuous(limits = c(-0.5,1.5)) +
  coord_fixed() +
  theme_graph()

# Arabidopsis FFL with feedback
genes <- tibble(name = c(expression("X"), expression("Y"), expression("Z")))
activationFFLA <- tibble(from = c(1, 1, 2, 3),
                         to =   c(2, 3, 3, 1))
gFFLA <- tbl_graph(edges = activationFFLA, directed = TRUE, nodes = genes)
plt_FFBH <- ggraph(gFFLA, layout = "manual", x = c(0,-0.7,0.7), y=c(1,0,0)) +
  geom_node_point(size = 16, color = colXOutline) +
  geom_node_point(size = 15, color = colX) +
  geom_node_text(aes(label = name), size = 8, color = "white", parse = T) +
  geom_edge_arc(aes(start_cap = circle(8, unit = "mm"),
                    end_cap = circle(8, unit = "mm")),
                arrow = arrow(type = "open", angle = 30, length = unit(3, 'mm')),
                strength = c(-0.2, 0.2, -0.2, 0.2)) +
  scale_x_continuous(limits = c(-1.5,1.5)) +
  scale_y_continuous(limits = c(-0.5,1.5)) +
  coord_fixed() +
  theme_graph()


# ODE system for feedback autoregulation:
ODEs_NAR <- function(t, state, parameters) {
  with (as.list(c(state, parameters)), {
    # step function leads to numerical issues in lsoda:
    #dZ <- bZ * (t > Xstart && t <= Xstop & Z<1) - aZ*Z
    # use Hill function instead:
    X <- XMult * (t > Xstart && t <= Xstop)
    dZ <- base + bZ * (t > Xstart && t <= Xstop) * ((X^Hilln)/(KXZ^Hilln + X^Hilln)) * ((KZ^Hilln)/(KZ^Hilln + Z^Hilln)) - aZ*Z
    return(list(c(dZ)))
  })
}

# ODE system for positive feedback autoregulation:
# When Z = 0, this is a stable point
# Will need to give a basal expression level of Z to activate the circuit
ODEs_PAR <- function(t, state, parameters) {
  with (as.list(c(state, parameters)), {
    # step function leads to numerical issues in lsoda:
    #dZ <- bZ * (t > Xstart && t <= Xstop & Z<1) - aZ*Z
    # use Hill function instead:
    X <- XMult * (t > Xstart && t <= Xstop)
    # Add small amount of Z for basal expression
    # addBase <- ((t > Xstart) && (t < Xstart + 0.01))
    
    dZ <- base + bZ * (X^Hilln/(KXZ^Hilln + X^Hilln)) * ((Z^Hilln)/((KZ^Hilln)+(Z^Hilln))) - aZ*Z
    #dZ <- base + X * bZ *((Z^Hilln)/((KZ^Hilln)+(Z^Hilln))) - aZ*Z
    #dZ <- bZ * (Z^Hilln)/((KZ^Hilln)+(Z^Hilln)) - aZ * Z
    return(list(c(dZ)))
  })
}

# ODE system for C1 FFL:
ODEs_C1_FFL <- function(t, state, parameters) {
  with (as.list(c(state, parameters)), {
    # step function leads to numerical issues in lsoda:
    #dZ <- bZ * (t > Xstart && t <= Xstop & Z<1) - aZ*Z
    # use Hill function instead:
    X <- XMult * (t > Xstart && t <= Xstop)
    dY <- bY * X^Hilln/(KY^Hilln + X^Hilln) - aY*Y
    dZ <- base + bZ * ((X * Y)^Hilln)/((KXZ^Hilln + X^Hilln) * (KY^Hilln + Y^Hilln)) - aZ*Z
    
    return(list(c(dY, dZ)))
  })
}

# I1 FFL
ODEs_I1_FFL <- function(t, state, parameters) {
  with (as.list(c(state, parameters)), {
    # step function leads to numerical issues in lsoda:
    #dZ <- bZ * (t > Xstart && t <= Xstop & Z<1) - aZ*Z
    # use Hill function instead:
    X <- XMult * (t > Xstart && t <= Xstop)
    
    dY <- bY * X^Hilln/(KY^Hilln + X^Hilln) - aY*Y
    dZ <- base + bZ * ((X * KY)^Hilln)/((KXZ^Hilln + X^Hilln) * (KY^Hilln + Y^Hilln)) - aZ*Z
    
    return(list(c(dY, dZ)))
  })
}

# ODE system for Feed forward/back hybrid:
ODEs_FFBH <- function(t, state, parameters) {
  with (as.list(c(state, parameters)), {
    
    # X is step function augmented by Z product
    # baseline X given by environmental cue in Xstart -> Xstop
    # change in X at Xstart is 1, change in X at Xstop is -1
    # Z adds a bit more signal to that
    
    # Manually set X
    X <- XMult * (t >= Xstart && t <= Xstop)
    X <- X + XH
    
    # Hill function component of X, XH
    dXH <- ( Z^Hilln / (KZX^Hilln + Z^Hilln) ) - aX*XH
    
    # Update X
    X <- X + dXH
    
    dY <- bY * X^Hilln/( KY^Hilln + X^Hilln ) - aY*Y
    dZ <- base + bZ *  ((X * Y)^Hilln)/((KXZ^Hilln + X^Hilln) * (KY^Hilln + Y^Hilln)) - aZ*Z
    
    return(list(c(dXH, dY, dZ)))
  })
}


plotDynamics_NAR <- function(Xstart = 1,
                             Xstop = 6,
                             tmax = 10,
                             dt = 0.1,
                             pars) {
  
  #pars <- c(base = 1, aZ = 1, bZ = 1, KXZ = 1, KZ = 1, Hilln = 1, XMult = 1)
  pars <- c(Xstart = Xstart, Xstop = Xstop, tmax = tmax, pars)
  iniState <- c(Z=0)
  times <- seq(0,tmax,by=dt)
  solution <- ode(iniState, times, ODEs_NAR, pars) %>%
    as.data.frame() %>%
    as_tibble() %>%
    mutate(X = ifelse(time > pars["Xstart"] & time <= pars["Xstop"], 1, 0)) %>%
    select(time, X, Z)
  
  plotZ <- ggplot(solution) +
    annotate("rect", xmin = Xstart, xmax = Xstop, ymin = 0, ymax = 1.05,
             alpha = .2, fill = colX) +
    geom_line(aes(time, Z), color = colZ, size = 1.5) +
    scale_y_continuous(limits = c(0,1.05)) +
    theme_bw(base_size = 16)
  
  plotAll <- plot_grid(plt_NAR, plotZ, NULL, nrow = 3, scale = c(1.5,1,1))
  plotAll
}

plotDynamics_PAR <- function(Xstart = 1,
                             Xstop = 6,
                             tmax = 10,
                             dt = 0.1,
                             pars) {
  
  #pars <- c(base = 1, aZ = 1, bZ = 1, KXZ = 1, KZ = 1, Hilln = 1, XMult = 1)
  pars <- c(Xstart = Xstart, Xstop = Xstop, tmax = tmax, pars)
  
  iniState <- c(Z=0)
  times <- seq(0,tmax,by=dt)
  solution <- ode(iniState, times, ODEs_PAR, pars) %>%
    as.data.frame() %>%
    as_tibble() %>%
    mutate(X = ifelse(time > pars["Xstart"] & time <= pars["Xstop"], 1, 0)) %>%
    select(time, X, Z)
  
  plotZ <- ggplot(solution) +
    annotate("rect", xmin = Xstart, xmax = Xstop, ymin = 0, ymax = 1.05,
             alpha = .2, fill = colX) +
    geom_line(aes(time, Z), color = colZ, size = 1.5) +
    scale_y_continuous(limits = c(0,1.05)) +
    theme_bw(base_size = 16)
  
  plotAll <- plot_grid(plt_PAR, plotZ, NULL, nrow = 3, scale = c(1.5,1,1))
  plotAll
}

plotDynamics_FFLC1 <- function(Xstart = 1,
                             Xstop = 6,
                             tmax = 10,
                             dt = 0.1,
                             pars) {
  
  pars <- c(Xstart = Xstart, Xstop = Xstop, tmax = tmax, pars)
  
  #pars <- c(aY = 1, bY = 1, KY = 1, KXZ = 1, aZ = 1, bZ = 1, Hilln = 1, XMult = 1, base = 1)
  iniState <- c(Y=0, Z=0)
  times <- seq(0,tmax,by=dt)
  solution <- ode(iniState, times, ODEs_C1_FFL, pars) %>%
    as.data.frame() %>%
    as_tibble() %>%
    mutate(X = ifelse(time > pars["Xstart"] & time <= pars["Xstop"], 1, 0)) %>%
    select(time, X, Y, Z)
  
  plotZ <- ggplot(solution) +
    annotate("rect", xmin = Xstart, xmax = Xstop, ymin = 0, ymax = 1.05,
             alpha = .2, fill = colX) +
    geom_line(aes(time, Z), color = colZ, size = 1.5) +
    scale_y_continuous(limits = c(0,1.05)) +
    theme_bw(base_size = 16)
  
  plotAll <- plot_grid(plt_FFLC1, plotZ, NULL, nrow = 3, scale = c(1.5,1,1))
  plotAll
}

plotDynamics_FFLI1 <- function(Xstart = 1,
                               Xstop = 6,
                               tmax = 10,
                               dt = 0.1,
                               pars) {
  
  pars <- c(Xstart = Xstart, Xstop = Xstop, tmax = tmax, pars)
  
  #pars <- c(aY = 1, bY = 1, KY = 1, KXZ = 1, aZ = 1, bZ = 1, Hilln = 1, XMult = 1, base = 1)
  iniState <- c(Y=0, Z=0)
  times <- seq(0,tmax,by=dt)
  solution <- ode(iniState, times, ODEs_I1_FFL, pars) %>%
    as.data.frame() %>%
    as_tibble() %>%
    mutate(X = ifelse(time > pars["Xstart"] & time <= pars["Xstop"], 1, 0)) %>%
    select(time, X, Y, Z)
  
  plotZ <- ggplot(solution) +
    annotate("rect", xmin = Xstart, xmax = Xstop, ymin = 0, ymax = 1.05,
             alpha = .2, fill = colX) +
    geom_line(aes(time, Z), color = colZ, size = 1.5) +
    scale_y_continuous(limits = c(0,1.05)) +
    theme_bw(base_size = 16)
  
  plotAll <- plot_grid(plt_FFLI1, plotZ, NULL, nrow = 3, scale = c(1.5,1,1))
  plotAll
}

plotDynamics_FFBH <- function(Xstart = 1,
                               Xstop = 6,
                               tmax = 10,
                               dt = 0.1,
                               pars) {
  
  pars <- c(Xstart = Xstart, Xstop = Xstop, tmax = tmax, pars)
  
  #pars <- c(aX = 1, KZX = 1, aY = 1, bY = 1, KY = 1, KXZ = 1, aZ = 1, bZ = 1, Hilln = 1, XMult = 1, base = 1)
  iniState <- c(XH=0, Y=0, Z=0)
  times <- seq(0,tmax,by=dt)
  solution <- ode(iniState, times, ODEs_FFBH, pars) %>%
    as.data.frame() %>%
    as_tibble() %>%
    mutate(X = ifelse(time > pars["Xstart"] & time <= pars["Xstop"], 1, 0)) %>%
    select(time, X, Y, Z)
  
  plotZ <- ggplot(solution) +
    annotate("rect", xmin = Xstart, xmax = Xstop, ymin = 0, ymax = 1.05,
             alpha = .2, fill = colX) +
    geom_line(aes(time, Z), color = colZ, size = 1.5) +
    scale_y_continuous(limits = c(0,1.05)) +
    theme_bw(base_size = 16)
  
  plotAll <- plot_grid(plt_FFBH, plotZ, NULL, nrow = 3, scale = c(1.5,1,1))
  plotAll
}

motif_strings <- c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH")
server<-function(input, output) {
  
  pars_NAR <- reactive({
    c(base = input$base_NAR, aZ = input$aZ_NAR, bZ = input$bZ_NAR, KXZ = input$KXZ_NAR,
      KZ = input$KZ_NAR, Hilln = input$Hilln_NAR, XMult = input$XMult_NAR)
  })
  
  pars_PAR <- reactive({
    c(base = input$base_PAR, aZ = input$aZ_PAR, bZ = input$bZ_PAR, KXZ = input$KXZ_PAR,
      KZ = input$KZ_PAR, Hilln = input$Hilln_PAR, XMult = input$XMult_PAR)
  })
  
  pars_FFLC1 <- reactive({
    c(base = input$base_FFLC1, aY = input$aY_FFLC1, bY = input$bY_FFLC1, KY = input$KY_FFLC1,
      aZ = input$aZ_FFLC1, bZ = input$bZ_FFLC1,  KXZ = input$KXZ_FFLC1, Hilln = input$Hilln_FFLC1, 
      XMult = input$XMult_FFLC1)
  })
  
  pars_FFLI1 <- reactive({
    c(base = input$base_FFLI1, aY = input$aY_FFLI1, bY = input$bY_FFLI1, KY = input$KY_FFLI1,
      aZ = input$aZ_FFLI1, bZ = input$bZ_FFLI1, KXZ = input$KXZ_FFLI1, Hilln = input$Hilln_FFLI1, 
      XMult = input$XMult_FFLI1)
  })
  
  pars_FFBH <- reactive({
    c(base = input$base_FFBH, aX = input$aX_FFBH, KZX = input$KZX_FFBH, aY = input$aY_FFBH, 
      bY = input$bY_FFBH, KY = input$KY_FFBH, aZ = input$aZ_FFBH, bZ = input$bZ_FFBH, 
      KXZ = input$KXZ_FFBH, Hilln = input$Hilln_FFBH, XMult = input$XMult_FFBH)
  })
  
  
  output$main_plot_NAR <- renderPlot({
    
      plotDynamics_NAR(Xstart = input$Xstart_NAR, Xstop = input$Xstop_NAR,
                       tmax = input$tmax_NAR, pars = pars_NAR())
    
  })
  
  output$main_plot_PAR <- renderPlot({

      plotDynamics_PAR(Xstart = input$Xstart_PAR, Xstop = input$Xstop_PAR,
                       tmax = input$tmax_PAR, pars = pars_PAR())

  })
  
  output$main_plot_FFLC1 <- renderPlot({
      plotDynamics_FFLC1(Xstart = input$Xstart_FFLC1, Xstop = input$Xstop_FFLC1,
                         tmax = input$tmax_FFLC1, pars = pars_FFLC1())

  })
  
  output$main_plot_FFLI1 <- renderPlot({

      plotDynamics_FFLI1(Xstart = input$Xstart_FFLI1, Xstop = input$Xstop_FFLI1,
                         tmax = input$tmax_FFLI1, pars = pars_FFLI1())

  })
  
  output$main_plot_FFBH <- renderPlot({

      plotDynamics_FFBH(Xstart = input$Xstart_FFBH, Xstop = input$Xstop_FFBH, 
                        tmax = input$tmax_FFBH, pars = pars_FFBH())

  })
  
  
}