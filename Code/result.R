## Assumptions ###

#packages
library(RColorBrewer) # for colors
library(copula) # for copulas
library(extraDistr) #for resampling
library(ggplot2) #for plotting
library(dplyr)
library(cowplot)

## Parameters
u_uptail <- seq(from = 0.95, to = 1, length.out = 125)
u_lowtail <- seq(from = 0, to = 0.05, length.out = 125)

B <- 30 # number of replications
n <- 125 # sample size for each replication (not too large, as otherwise the empirical copulas are too close to the true copula)
nu <- 4 #degrees of freedom for student-t

cop.names <- c("Empirical copula", "Empirical beta copula", "Empirical checkerboard copula", "Empirical betabinomial 4", "Empirical beta2 4", "Empirical binomial 4") # empirical copula names

#colour palettes
basecols <- c("#000000", brewer.pal(8, name = "Dark2")[c(8,7,3,1,5,4,2,6)])
mypal <- colorRampPalette(basecols, space = "Lab") # define palette function to interpolate colors
palette(mypal(8)) # get some colors and set them as the new palette

# a palette is a named vector
ColorPalettePlot <- function(myPalette) {
  myPalette = rev(myPalette)
  # create data that references the palette
  colorKey = data.frame(colorName=names(myPalette))
  # plot with ggplot, referencing the palette
  ggplot(data=colorKey, aes(x=1, y = 1:nrow(colorKey), fill=colorName, label=colorName)) +
    geom_tile() +
    scale_fill_manual(values = myPalette) +
    theme_void()+
    theme(legend.position="none") + 
    geom_text(color="white")
}

#color palette plot
myPalette = c("True Copula" = "#000000", "Empirical Copula" = "#71685D", "Empirical Beta Copula" = "#A0734E", "Empirical Checkerboard Copula" = "#5E8599", "EBC-adapted empirical copula with beta-binomial survival margins of rho = 4" = "#53A24A", "EBC-adapted empirical copula with beta survival margins of rho = 4" = "#CA6470", "EBC-adapted empirical copula with binomial survival margins" = "#DC5924")
ColorPalettePlot(myPalette)
## Simulations ## 

#Notes
#Copulas: Gaussian, Student-t, Clayton, Gumbel
#Kendall's Tau: {-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75}
#different dimensions: d = {2, 3, 4, 5}
#Cramer-von-Mises test-statistic for each case

#1. Survival / Exceedance Probabilities

#create matrix for CvM s
t4_CvM_mat_s = matrix(NA, nrow = length(cop.names), ncol = length(2:5)*length(c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)))
N_CvM_mat_s = matrix(NA, nrow = length(cop.names), ncol = length(2:5)*length(c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)))
C_CvM_mat_s = matrix(NA, nrow = length(cop.names), ncol = length(2:5)*length(c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)))
G_CvM_mat_s = matrix(NA, nrow = length(cop.names), ncol = length(2:5)*length(c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)))

colnames(t4_CvM_mat_s) = paste0(expand.grid(2:5, c("neg0.75", "neg0.5", "neg0.25", "0", "0.25", "0.5", "0.75"))[,1], sep = "_", expand.grid(2:5,c("neg0.75", "neg0.5", "neg0.25", "0", "0.25", "0.5", "0.75"))[,2] )
rownames(t4_CvM_mat_s) = cop.names
colnames(N_CvM_mat_s) = paste0(expand.grid(2:5, c("neg0.75", "neg0.5", "neg0.25", "0", "0.25", "0.5", "0.75"))[,1], sep = "_", expand.grid(2:5,c("neg0.75", "neg0.5", "neg0.25", "0", "0.25", "0.5", "0.75"))[,2] )
rownames(N_CvM_mat_s) = cop.names
colnames(C_CvM_mat_s) = paste0(expand.grid(2:5, c("neg0.75", "neg0.5", "neg0.25", "0", "0.25", "0.5", "0.75"))[,1], sep = "_", expand.grid(2:5,c("neg0.75", "neg0.5", "neg0.25", "0", "0.25", "0.5", "0.75"))[,2] )
rownames(C_CvM_mat_s) = cop.names
colnames(G_CvM_mat_s) = paste0(expand.grid(2:5, c("neg0.75", "neg0.5", "neg0.25", "0", "0.25", "0.5", "0.75"))[,1], sep = "_", expand.grid(2:5,c("neg0.75", "neg0.5", "neg0.25", "0", "0.25", "0.5", "0.75"))[,2] )
rownames(G_CvM_mat_s) = cop.names

#s prob
set.seed(2024)
for (d in 2:5) {
  for (tau in c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)) {
    #create the true copula
    if (d <= 2) {tcop <- tCopula(iTau(tCopula(df = nu), tau = tau), dim = d)}
    else {}
    if (d <= 2) {ncop <- normalCopula(iTau(normalCopula(), tau = tau), dim = d)}
    else {}    
    if (d == 2 |(d > 2 & tau >= 0)) {ccop <- claytonCopula(iTau(claytonCopula(), tau = tau), dim = d)}
    else {}
    gcop <- gumbelCopula(iTau(gumbelCopula(), tau = tau), dim = d)
    
    
    #create exceedance probabilities of the empirical version of the true copula
    if (d <= 2) {assign(paste0("res.t4.", d, "d.tau", if (tau < 0) {paste0("neg", abs(tau))} else {tau}, ".s"), emp_survivalprob(tcop))}
    else {}
    if (d <= 2) {assign(paste0("res.N.", d, "d.tau", if (tau < 0) {paste0("neg", abs(tau))} else {tau}, ".s"), emp_survivalprob(ncop))}
    else {}
    if (d == 2 |(d > 2 & tau >= 0)) {assign(paste0("res.C.", d, "d.tau", if (tau < 0) {paste0("neg", abs(tau))} else {tau}, ".s"), emp_survivalprob(ccop))}
    else {}
    assign(paste0("res.G.", d, "d.tau", if (tau < 0) {paste0("neg", abs(tau))} else {tau}, ".s"), emp_survivalprob(gcop))
    
  }
}

#CvM matrix for s
j = -1*length(2:5)
k = 0
for (d in 2:5) {
  j = j + 1
  k = j
  for (tau in c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)) {
    j = j + 4
    t4_CvM_mat_s[,j] = if (d <= 2) {CvM(eval(parse(text =  paste0("res.t4.", d, "d.tau", if (tau < 0) {paste0("neg", abs(tau))} else {tau}, ".s") )), n) } else {NA}
    N_CvM_mat_s[,j] = if (d <= 2) {CvM(eval(parse(text =  paste0("res.N.", d, "d.tau", if (tau < 0) {paste0("neg", abs(tau))} else {tau}, ".s") )), n) } else {NA}
    C_CvM_mat_s[,j] = if (d == 2 |(d > 2 & tau >= 0)) {CvM(eval(parse(text =  paste0("res.C.", d, "d.tau", if (tau < 0) {paste0("neg", abs(tau))} else {tau}, ".s") )), n) } else {NA}
    G_CvM_mat_s[,j] = CvM(eval(parse(text = paste0("res.G.", d, "d.tau", if (tau < 0) {paste0("neg", abs(tau))} else {tau}, ".s") )), n)
    
  }
  j = k
}

########### exceedance probability plots
### print for student's t

#d = 2
par(mfrow = c(3,3))
plot_res_mat(res.t4.2d.tauneg0.75.s, u_uptail, -0.75)
plot_res_mat(res.t4.2d.tauneg0.5.s, u_uptail, -0.5)
plot_res_mat(res.t4.2d.tauneg0.25.s, u_uptail, -0.25)
plot_res_mat(res.t4.2d.tau0.s, u_uptail, 0)
plot_res_mat(res.t4.2d.tau0.25.s, u_uptail, 0.25)
plot_res_mat(res.t4.2d.tau0.5.s, u_uptail, 0.5)
plot_res_mat(res.t4.2d.tau0.75.s, u_uptail, 0.75)

### print for gaussian
#d = 2
par(mfrow = c(3,3))
plot_res_mat(res.N.2d.tauneg0.75.s, u_uptail, -0.75)
plot_res_mat(res.N.2d.tauneg0.5.s, u_uptail, -0.5)
plot_res_mat(res.N.2d.tauneg0.25.s, u_uptail, -0.25)
plot_res_mat(res.N.2d.tau0.s, u_uptail, 0)
plot_res_mat(res.N.2d.tau0.25.s, u_uptail, 0.25)
plot_res_mat(res.N.2d.tau0.5.s, u_uptail, 0.5)
plot_res_mat(res.N.2d.tau0.75.s, u_uptail, 0.75)

### print for clayton

#d = 2
par(mfrow = c(2,2))
plot_res_mat(res.C.2d.tau0.25.s, u_uptail, 0.25)
plot_res_mat(res.C.2d.tau0.5.s, u_uptail, 0.5)
plot_res_mat(res.C.2d.tau0.75.s, u_uptail, 0.75)

#d = 3
par(mfrow = c(2,2))
plot_res_mat(res.C.3d.tau0.25.s, u_uptail, 0.25)
plot_res_mat(res.C.3d.tau0.5.s, u_uptail, 0.5)
plot_res_mat(res.C.3d.tau0.75.s, u_uptail, 0.75)

#d = 4
par(mfrow = c(2,2))
plot_res_mat(res.C.4d.tau0.25.s, u_uptail, 0.25)
plot_res_mat(res.C.4d.tau0.5.s, u_uptail, 0.5)
plot_res_mat(res.C.4d.tau0.75.s, u_uptail, 0.75)

#d = 5
par(mfrow = c(2,2))
plot_res_mat(res.C.5d.tau0.25.s, u_uptail, 0.25)
plot_res_mat(res.C.5d.tau0.5.s, u_uptail, 0.5)
plot_res_mat(res.C.5d.tau0.75.s, u_uptail, 0.75)

### print for gumbel

#d = 2
par(mfrow = c(2,2))
plot_res_mat(res.G.2d.tau0.25.s, u_uptail, 0.25)
plot_res_mat(res.G.2d.tau0.5.s, u_uptail, 0.5)
plot_res_mat(res.G.2d.tau0.75.s, u_uptail, 0.75)

#d = 3
par(mfrow = c(2,2))
plot_res_mat(res.G.3d.tau0.25.s, u_uptail, 0.25)
plot_res_mat(res.G.3d.tau0.5.s, u_uptail, 0.5)
plot_res_mat(res.G.3d.tau0.75.s, u_uptail, 0.75)

#d = 4
par(mfrow = c(2,2))
plot_res_mat(res.G.4d.tau0.25.s, u_uptail, 0.25)
plot_res_mat(res.G.4d.tau0.5.s, u_uptail, 0.5)
plot_res_mat(res.G.4d.tau0.75.s, u_uptail, 0.75)

#d = 5
par(mfrow = c(2,2))
plot_res_mat(res.G.5d.tau0.25.s, u_uptail, 0.25)
plot_res_mat(res.G.5d.tau0.5.s, u_uptail, 0.5)
plot_res_mat(res.G.5d.tau0.75.s, u_uptail, 0.75)

#2. Cumulative Probabilities
#create matrix for CvM c
t4_CvM_mat_c = matrix(NA, nrow = length(cop.names), ncol = length(2:5)*length(c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)))
N_CvM_mat_c = matrix(NA, nrow = length(cop.names), ncol = length(2:5)*length(c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)))
C_CvM_mat_c = matrix(NA, nrow = length(cop.names), ncol = length(2:5)*length(c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)))
G_CvM_mat_c = matrix(NA, nrow = length(cop.names), ncol = length(2:5)*length(c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)))

colnames(t4_CvM_mat_c) = paste0(expand.grid(2:5, c("neg0.75", "neg0.5", "neg0.25", "0", "0.25", "0.5", "0.75"))[,1], sep = "_", expand.grid(2:5,c("neg0.75", "neg0.5", "neg0.25", "0", "0.25", "0.5", "0.75"))[,2] )
rownames(t4_CvM_mat_c) = cop.names
colnames(N_CvM_mat_c) = paste0(expand.grid(2:5, c("neg0.75", "neg0.5", "neg0.25", "0", "0.25", "0.5", "0.75"))[,1], sep = "_", expand.grid(2:5,c("neg0.75", "neg0.5", "neg0.25", "0", "0.25", "0.5", "0.75"))[,2] )
rownames(N_CvM_mat_c) = cop.names
colnames(C_CvM_mat_c) = paste0(expand.grid(2:5, c("neg0.75", "neg0.5", "neg0.25", "0", "0.25", "0.5", "0.75"))[,1], sep = "_", expand.grid(2:5,c("neg0.75", "neg0.5", "neg0.25", "0", "0.25", "0.5", "0.75"))[,2] )
rownames(C_CvM_mat_c) = cop.names
colnames(G_CvM_mat_c) = paste0(expand.grid(2:5, c("neg0.75", "neg0.5", "neg0.25", "0", "0.25", "0.5", "0.75"))[,1], sep = "_", expand.grid(2:5,c("neg0.75", "neg0.5", "neg0.25", "0", "0.25", "0.5", "0.75"))[,2] )
rownames(G_CvM_mat_c) = cop.names



#c prob
set.seed(2024)
for (d in 5) {
  for (tau in c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)) {
    #create the true copula
    if (d <= 2) {tcop <- tCopula(iTau(tCopula(df = nu), tau = tau), dim = d)}
    else {}
    if (d <= 2) {ncop <- normalCopula(iTau(normalCopula(), tau = tau), dim = d)}
    else {}    
    if (d == 2 |(d > 2 & tau >= 0)) {ccop <- claytonCopula(iTau(claytonCopula(), tau = tau), dim = d)}
    else {}
    gcop <- gumbelCopula(iTau(gumbelCopula(), tau = tau), dim = d)
    
    
    #create exceedance probabilities of the empirical version of the true copula
    if (d <= 2) {assign(paste0("res.t4.", d, "d.tau", if (tau < 0) {paste0("neg", abs(tau))} else {tau}, ".c"), emp_cprob(tcop))}
    else {}
    if (d <= 2) {assign(paste0("res.N.", d, "d.tau", if (tau < 0) {paste0("neg", abs(tau))} else {tau}, ".c"), emp_cprob(ncop))}
    else {}
    if (d == 2 |(d > 2 & tau >= 0)) {assign(paste0("res.C.", d, "d.tau", if (tau < 0) {paste0("neg", abs(tau))} else {tau}, ".c"), emp_cprob(ccop))}
    else {}
    assign(paste0("res.G.", d, "d.tau", if (tau < 0) {paste0("neg", abs(tau))} else {tau}, ".c"), emp_cprob(gcop))
    
  }
}



#CvM matrix for c
j = -1*length(2:5)
k = 0
for (d in 2:5) {
  j = j + 1
  k = j
  for (tau in c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)) {
    j = j + 4
    t4_CvM_mat_c[,j] = if (d <= 2) {CvM(eval(parse(text =  paste0("res.t4.", d, "d.tau", if (tau < 0) {paste0("neg", abs(tau))} else {tau}, ".c") )), n) } else {NA}
    N_CvM_mat_c[,j] = if (d <= 2) {CvM(eval(parse(text =  paste0("res.N.", d, "d.tau", if (tau < 0) {paste0("neg", abs(tau))} else {tau}, ".c") )), n) } else {NA}
    C_CvM_mat_c[,j] = if (d == 2 |(d > 2 & tau >= 0)) {CvM(eval(parse(text =  paste0("res.C.", d, "d.tau", if (tau < 0) {paste0("neg", abs(tau))} else {tau}, ".c") )), n) } else {NA}
    G_CvM_mat_c[,j] = CvM(eval(parse(text = paste0("res.G.", d, "d.tau", if (tau < 0) {paste0("neg", abs(tau))} else {tau}, ".c") )), n)
    
  }
  j = k
}

########### cumulative probability plots
### print for student's t

#d = 2
par(mfrow = c(3,3))
plot_res_mat(res.t4.2d.tauneg0.75.c, u_lowtail, -0.75)
plot_res_mat(res.t4.2d.tauneg0.5.c, u_lowtail, -0.5)
plot_res_mat(res.t4.2d.tauneg0.25.c, u_lowtail, -0.25)
plot_res_mat(res.t4.2d.tau0.c, u_lowtail, 0)
plot_res_mat(res.t4.2d.tau0.25.c, u_lowtail, 0.25)
plot_res_mat(res.t4.2d.tau0.5.c, u_lowtail, 0.5)
plot_res_mat(res.t4.2d.tau0.75.c, u_lowtail, 0.75)

### print for gaussian
#d = 2
par(mfrow = c(3,3))
plot_res_mat(res.N.2d.tauneg0.75.c, u_lowtail, -0.75)
plot_res_mat(res.N.2d.tauneg0.5.c, u_lowtail, -0.5)
plot_res_mat(res.N.2d.tauneg0.25.c, u_lowtail, -0.25)
plot_res_mat(res.N.2d.tau0.c, u_lowtail, 0)
plot_res_mat(res.N.2d.tau0.25.c, u_lowtail, 0.25)
plot_res_mat(res.N.2d.tau0.5.c, u_lowtail, 0.5)
plot_res_mat(res.N.2d.tau0.75.c, u_lowtail, 0.75)

### print for clayton

#d = 2
par(mfrow = c(2,2))
plot_res_mat(res.C.2d.tau0.25.c, u_lowtail, 0.25)
plot_res_mat(res.C.2d.tau0.5.c, u_lowtail, 0.5)
plot_res_mat(res.C.2d.tau0.75.c, u_lowtail, 0.75)

#d = 3
par(mfrow = c(2,2))
plot_res_mat(res.C.3d.tau0.25.c, u_lowtail, 0.25)
plot_res_mat(res.C.3d.tau0.5.c, u_lowtail, 0.5)
plot_res_mat(res.C.3d.tau0.75.c, u_lowtail, 0.75)

#d = 4
par(mfrow = c(2,2))
plot_res_mat(res.C.4d.tau0.25.c, u_lowtail, 0.25)
plot_res_mat(res.C.4d.tau0.5.c, u_lowtail, 0.5)
plot_res_mat(res.C.4d.tau0.75.c, u_lowtail, 0.75)

#d = 5
par(mfrow = c(2,2))
plot_res_mat(res.C.5d.tau0.25.c, u_lowtail, 0.25)
plot_res_mat(res.C.5d.tau0.5.c, u_lowtail, 0.5)
plot_res_mat(res.C.5d.tau0.75.c, u_lowtail, 0.75)

### print for gumbel

#d = 2
par(mfrow = c(2,2))
plot_res_mat(res.G.2d.tau0.25.c, u_lowtail, 0.25)
plot_res_mat(res.G.2d.tau0.5.c, u_lowtail, 0.5)
plot_res_mat(res.G.2d.tau0.75.c, u_lowtail, 0.75)

#d = 3
par(mfrow = c(2,2))
plot_res_mat(res.G.3d.tau0.25.c, u_lowtail, 0.25)
plot_res_mat(res.G.3d.tau0.5.c, u_lowtail, 0.5)
plot_res_mat(res.G.3d.tau0.75.c, u_lowtail, 0.75)

#d = 4
par(mfrow = c(2,2))
plot_res_mat(res.G.4d.tau0.25.c, u_lowtail, 0.25)
plot_res_mat(res.G.4d.tau0.5.c, u_lowtail, 0.5)
plot_res_mat(res.G.4d.tau0.75.c, u_lowtail, 0.75)

#d = 5
par(mfrow = c(2,2))
plot_res_mat(res.G.5d.tau0.25.c, u_lowtail, 0.25)
plot_res_mat(res.G.5d.tau0.5.c, u_lowtail, 0.5)
plot_res_mat(res.G.5d.tau0.75.c, u_lowtail, 0.75)


###CVM stuff

#Graphing CVM!!

cop.short.names = c("EC", "EBC", "ECC", "EBB4", "EBeta4", "EBinom")
C_CvM_mat_c_df = data.frame(cbind(cop.short.names,C_CvM_mat_c))
C_CvM_mat_s_df = data.frame(cbind(cop.short.names,C_CvM_mat_s))
t4_CvM_mat_c_df = data.frame(cbind(cop.short.names,t4_CvM_mat_c))
t4_CvM_mat_s_df = data.frame(cbind(cop.short.names,t4_CvM_mat_s))
G_CvM_mat_c_df = data.frame(cbind(cop.short.names,G_CvM_mat_c))
G_CvM_mat_s_df = data.frame(cbind(cop.short.names,G_CvM_mat_s))
N_CvM_mat_c_df = data.frame(cbind(cop.short.names,N_CvM_mat_c))
N_CvM_mat_s_df = data.frame(cbind(cop.short.names,N_CvM_mat_s))


#d = 2, c
theme_update(plot.title = element_text(hjust = 0.5))

t4_2d_neg0.75_c_plot <- ggplot(data = t4_CvM_mat_c_df, aes(x = cop.short.names, y = as.numeric(X2_neg0.75))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 2)") + ylab("CvM") + 
  ggtitle("tau = -0.75")
t4_2d_neg0.5_c_plot <- ggplot(data = t4_CvM_mat_c_df, aes(x = cop.short.names, y = as.numeric(X2_neg0.5))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 2)") + ylab("CvM") + 
  ggtitle("tau = -0.5")
t4_2d_neg0.25_c_plot <- ggplot(data = t4_CvM_mat_c_df, aes(x = cop.short.names, y = as.numeric(X2_neg0.25))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 2)") + ylab("CvM") + 
  ggtitle("tau = -0.25")
t4_2d_0_c_plot <- ggplot(data = t4_CvM_mat_c_df, aes(x = cop.short.names, y = as.numeric(X2_0))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 2)") + ylab("CvM") + 
  ggtitle("tau = 0")
t4_2d_0.25_c_plot <- ggplot(data = t4_CvM_mat_c_df, aes(x = cop.short.names, y = as.numeric(X2_0.25))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 2)") + ylab("CvM") + 
  ggtitle("tau = 0.25")
t4_2d_0.5_c_plot <- ggplot(data = t4_CvM_mat_c_df, aes(x = cop.short.names, y = as.numeric(X2_0.5))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 2)") + ylab("CvM") + 
  ggtitle("tau = 0.5")
t4_2d_0.75_c_plot <- ggplot(data = t4_CvM_mat_c_df, aes(x = cop.short.names, y = as.numeric(X2_0.75))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 2)") + ylab("CvM") + 
  ggtitle("tau = 0.75")

plot_grid(t4_2d_neg0.75_c_plot, t4_2d_neg0.5_c_plot, t4_2d_neg0.25_c_plot, t4_2d_0_c_plot, t4_2d_0.25_c_plot, t4_2d_0.5_c_plot, t4_2d_0.75_c_plot, labels = "AUTO", ncol = 3)

N_2d_neg0.75_c_plot <- ggplot(data = N_CvM_mat_c_df, aes(x = cop.short.names, y = as.numeric(X2_neg0.75))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 2)") + ylab("CvM") + 
  ggtitle("tau = -0.75")
N_2d_neg0.5_c_plot <- ggplot(data = N_CvM_mat_c_df, aes(x = cop.short.names, y = as.numeric(X2_neg0.5))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 2)") + ylab("CvM") + 
  ggtitle("tau = -0.5")
N_2d_neg0.25_c_plot <- ggplot(data = N_CvM_mat_c_df, aes(x = cop.short.names, y = as.numeric(X2_neg0.25))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 2)") + ylab("CvM") + 
  ggtitle("tau = -0.25")
N_2d_0_c_plot <- ggplot(data = N_CvM_mat_c_df, aes(x = cop.short.names, y = as.numeric(X2_0))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 2)") + ylab("CvM") + 
  ggtitle("tau = 0")
N_2d_0.25_c_plot <- ggplot(data = N_CvM_mat_c_df, aes(x = cop.short.names, y = as.numeric(X2_0.25))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 2)") + ylab("CvM") + 
  ggtitle("tau = 0.25")
N_2d_0.5_c_plot <- ggplot(data = N_CvM_mat_c_df, aes(x = cop.short.names, y = as.numeric(X2_0.5))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 2)") + ylab("CvM") + 
  ggtitle("tau = 0.5")
N_2d_0.75_c_plot <- ggplot(data = N_CvM_mat_c_df, aes(x = cop.short.names, y = as.numeric(X2_0.75))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 2)") + ylab("CvM") + 
  ggtitle("tau = 0.75")

plot_grid(N_2d_neg0.75_c_plot, N_2d_neg0.5_c_plot, N_2d_neg0.25_c_plot, N_2d_0_c_plot, N_2d_0.25_c_plot, N_2d_0.5_c_plot, N_2d_0.75_c_plot, labels = "AUTO", ncol = 3)

G_2d_0.25_c_plot <- ggplot(data = G_CvM_mat_c_df, aes(x = cop.short.names, y = as.numeric(X2_0.25))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 2)") + ylab("CvM")  + 
  ggtitle("tau = 0.25")
G_2d_0.5_c_plot <- ggplot(data = G_CvM_mat_c_df, aes(x = cop.short.names, y = as.numeric(X2_0.5))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 2)") + ylab("CvM")  + 
  ggtitle("tau = 0.5")
G_2d_0.75_c_plot <- ggplot(data = G_CvM_mat_c_df, aes(x = cop.short.names, y = as.numeric(X2_0.75))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 2)") + ylab("CvM")  + 
  ggtitle("tau = 0.75")

plot_grid(G_2d_0.25_c_plot, G_2d_0.5_c_plot, G_2d_0.75_c_plot, labels = "AUTO", ncol = 2)

G_3d_0.25_c_plot <- ggplot(data = G_CvM_mat_c_df, aes(x = cop.short.names, y = as.numeric(X3_0.25))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 3)") + ylab("CvM")  + 
  ggtitle("tau = 0.25")
G_3d_0.5_c_plot <- ggplot(data = G_CvM_mat_c_df, aes(x = cop.short.names, y = as.numeric(X3_0.5))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 3)") + ylab("CvM")  + 
  ggtitle("tau = 0.5")
G_3d_0.75_c_plot <- ggplot(data = G_CvM_mat_c_df, aes(x = cop.short.names, y = as.numeric(X3_0.75))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 3)") + ylab("CvM")  + 
  ggtitle("tau = 0.75")

plot_grid(G_3d_0.25_c_plot, G_3d_0.5_c_plot, G_3d_0.75_c_plot, labels = "AUTO", ncol = 2)

G_4d_0.25_c_plot <- ggplot(data = G_CvM_mat_c_df, aes(x = cop.short.names, y = as.numeric(X4_0.25))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 4)") + ylab("CvM")  + 
  ggtitle("tau = 0.25")
G_4d_0.5_c_plot <- ggplot(data = G_CvM_mat_c_df, aes(x = cop.short.names, y = as.numeric(X4_0.5))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 4)") + ylab("CvM")  + 
  ggtitle("tau = 0.5")
G_4d_0.75_c_plot <- ggplot(data = G_CvM_mat_c_df, aes(x = cop.short.names, y = as.numeric(X4_0.75))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 4)") + ylab("CvM")  + 
  ggtitle("tau = 0.75")

plot_grid(G_4d_0.25_c_plot, G_4d_0.5_c_plot, G_4d_0.75_c_plot, labels = "AUTO", ncol = 2)

G_5d_0.25_c_plot <- ggplot(data = G_CvM_mat_c_df, aes(x = cop.short.names, y = as.numeric(X5_0.25))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 5)") + ylab("CvM")  + 
  ggtitle("tau = 0.25")
G_5d_0.5_c_plot <- ggplot(data = G_CvM_mat_c_df, aes(x = cop.short.names, y = as.numeric(X5_0.5))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 5)") + ylab("CvM")  + 
  ggtitle("tau = 0.5")
G_5d_0.75_c_plot <- ggplot(data = G_CvM_mat_c_df, aes(x = cop.short.names, y = as.numeric(X5_0.75))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 5)") + ylab("CvM")  + 
  ggtitle("tau = 0.75")

plot_grid(G_5d_0.25_c_plot, G_5d_0.5_c_plot, G_5d_0.75_c_plot, labels = "AUTO", ncol = 2)

C_2d_0.25_c_plot <- ggplot(data = C_CvM_mat_c_df, aes(x = cop.short.names, y = as.numeric(X2_0.25))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 2)") + ylab("CvM")  + 
  ggtitle("tau = 0.25")
C_2d_0.5_c_plot <- ggplot(data = C_CvM_mat_c_df, aes(x = cop.short.names, y = as.numeric(X2_0.5))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 2)") + ylab("CvM")  + 
  ggtitle("tau = 0.5")
C_2d_0.75_c_plot <- ggplot(data = C_CvM_mat_c_df, aes(x = cop.short.names, y = as.numeric(X2_0.75))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 2)") + ylab("CvM")  + 
  ggtitle("tau = 0.75")

plot_grid(C_2d_0.25_c_plot, C_2d_0.5_c_plot, C_2d_0.75_c_plot, labels = "AUTO", ncol = 2)

C_3d_0.25_c_plot <- ggplot(data = C_CvM_mat_c_df, aes(x = cop.short.names, y = as.numeric(X3_0.25))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 3)") + ylab("CvM")  + 
  ggtitle("tau = 0.25")
C_3d_0.5_c_plot <- ggplot(data = C_CvM_mat_c_df, aes(x = cop.short.names, y = as.numeric(X3_0.5))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 3)") + ylab("CvM")  + 
  ggtitle("tau = 0.5")
C_3d_0.75_c_plot <- ggplot(data = C_CvM_mat_c_df, aes(x = cop.short.names, y = as.numeric(X3_0.75))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 3)") + ylab("CvM")  + 
  ggtitle("tau = 0.75")

plot_grid(C_3d_0.25_c_plot, C_3d_0.5_c_plot, C_3d_0.75_c_plot, labels = "AUTO", ncol = 2)

C_4d_0.25_c_plot <- ggplot(data = C_CvM_mat_c_df, aes(x = cop.short.names, y = as.numeric(X4_0.25))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 4)") + ylab("CvM")  + 
  ggtitle("tau = 0.25")
C_4d_0.5_c_plot <- ggplot(data = C_CvM_mat_c_df, aes(x = cop.short.names, y = as.numeric(X4_0.5))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 4)") + ylab("CvM")  + 
  ggtitle("tau = 0.5")
C_4d_0.75_c_plot <- ggplot(data = C_CvM_mat_c_df, aes(x = cop.short.names, y = as.numeric(X4_0.75))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 4)") + ylab("CvM")  + 
  ggtitle("tau = 0.75")

plot_grid(C_4d_0.25_c_plot, C_4d_0.5_c_plot, C_4d_0.75_c_plot, labels = "AUTO", ncol = 2)

C_5d_0.25_c_plot <- ggplot(data = C_CvM_mat_c_df, aes(x = cop.short.names, y = as.numeric(X5_0.25))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 5)") + ylab("CvM")  + 
  ggtitle("tau = 0.25")
C_5d_0.5_c_plot <- ggplot(data = C_CvM_mat_c_df, aes(x = cop.short.names, y = as.numeric(X5_0.5))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 5)") + ylab("CvM")  + 
  ggtitle("tau = 0.5")
C_5d_0.75_c_plot <- ggplot(data = C_CvM_mat_c_df, aes(x = cop.short.names, y = as.numeric(X5_0.75))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 5)") + ylab("CvM")  + 
  ggtitle("tau = 0.75")

plot_grid(C_5d_0.25_c_plot, C_5d_0.5_c_plot, C_5d_0.75_c_plot, labels = "AUTO", ncol = 2)

#d = 2, s
t4_2d_neg0.75_s_plot <- ggplot(data = t4_CvM_mat_s_df, aes(x = cop.short.names, y = as.numeric(X2_neg0.75))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 2)") + ylab("CvM")  + 
  ggtitle("tau = -0.75")
t4_2d_neg0.5_s_plot <- ggplot(data = t4_CvM_mat_s_df, aes(x = cop.short.names, y = as.numeric(X2_neg0.5))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 2)") + ylab("CvM")  + 
  ggtitle("tau = -0.5")
t4_2d_neg0.25_s_plot <- ggplot(data = t4_CvM_mat_s_df, aes(x = cop.short.names, y = as.numeric(X2_neg0.25))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 2)") + ylab("CvM")  + 
  ggtitle("tau = -0.25")
t4_2d_0_s_plot <- ggplot(data = t4_CvM_mat_s_df, aes(x = cop.short.names, y = as.numeric(X2_0))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 2)") + ylab("CvM")  + 
  ggtitle("tau = 0")
t4_2d_0.25_s_plot <- ggplot(data = t4_CvM_mat_s_df, aes(x = cop.short.names, y = as.numeric(X2_0.25))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 2)") + ylab("CvM")  + 
  ggtitle("tau = 0.25")
t4_2d_0.5_s_plot <- ggplot(data = t4_CvM_mat_s_df, aes(x = cop.short.names, y = as.numeric(X2_0.5))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 2)") + ylab("CvM")  + 
  ggtitle("tau = 0.5")
t4_2d_0.75_s_plot <- ggplot(data = t4_CvM_mat_s_df, aes(x = cop.short.names, y = as.numeric(X2_0.75))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 2)") + ylab("CvM")  + 
  ggtitle("tau = 0.75")

plot_grid(t4_2d_neg0.75_s_plot, t4_2d_neg0.5_s_plot, t4_2d_neg0.25_s_plot, t4_2d_0_s_plot, t4_2d_0.25_s_plot, t4_2d_0.5_s_plot, t4_2d_0.75_s_plot, labels = "AUTO", ncol = 3)

N_2d_neg0.75_s_plot <- ggplot(data = N_CvM_mat_s_df, aes(x = cop.short.names, y = as.numeric(X2_neg0.75))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 2)") + ylab("CvM")  + 
  ggtitle("tau = -0.75")
N_2d_neg0.5_s_plot <- ggplot(data = N_CvM_mat_s_df, aes(x = cop.short.names, y = as.numeric(X2_neg0.5))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 2)") + ylab("CvM")  + 
  ggtitle("tau = -0.5")
N_2d_neg0.25_s_plot <- ggplot(data = N_CvM_mat_s_df, aes(x = cop.short.names, y = as.numeric(X2_neg0.25))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 2)") + ylab("CvM")  + 
  ggtitle("tau = -0.25")
N_2d_0_s_plot <- ggplot(data = N_CvM_mat_s_df, aes(x = cop.short.names, y = as.numeric(X2_0))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 2)") + ylab("CvM")  + 
  ggtitle("tau = 0")
N_2d_0.25_s_plot <- ggplot(data = N_CvM_mat_s_df, aes(x = cop.short.names, y = as.numeric(X2_0.25))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 2)") + ylab("CvM")  + 
  ggtitle("tau = 0.25")
N_2d_0.5_s_plot <- ggplot(data = N_CvM_mat_s_df, aes(x = cop.short.names, y = as.numeric(X2_0.5))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 2)") + ylab("CvM")  + 
  ggtitle("tau = 0.5")
N_2d_0.75_s_plot <- ggplot(data = N_CvM_mat_s_df, aes(x = cop.short.names, y = as.numeric(X2_0.75))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 2)") + ylab("CvM")  + 
  ggtitle("tau = 0.75")

plot_grid(N_2d_neg0.75_s_plot, N_2d_neg0.5_s_plot, N_2d_neg0.25_s_plot, N_2d_0_s_plot, N_2d_0.25_s_plot, N_2d_0.5_s_plot, N_2d_0.75_s_plot, labels = "AUTO", ncol = 3)

G_2d_0.25_s_plot <- ggplot(data = G_CvM_mat_s_df, aes(x = cop.short.names, y = as.numeric(X2_0.25))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 2)") + ylab("CvM")  + 
  ggtitle("tau = 0.25")
G_2d_0.5_s_plot <- ggplot(data = G_CvM_mat_s_df, aes(x = cop.short.names, y = as.numeric(X2_0.5))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 2)") + ylab("CvM")  + 
  ggtitle("tau = 0.5")
G_2d_0.75_s_plot <- ggplot(data = G_CvM_mat_s_df, aes(x = cop.short.names, y = as.numeric(X2_0.75))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 2)") + ylab("CvM")  + 
  ggtitle("tau = 0.75")

plot_grid(G_2d_0.25_s_plot, G_2d_0.5_s_plot, G_2d_0.75_s_plot, labels = "AUTO", ncol = 2)

G_3d_0.25_s_plot <- ggplot(data = G_CvM_mat_s_df, aes(x = cop.short.names, y = as.numeric(X3_0.25))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 3)") + ylab("CvM")  + 
  ggtitle("tau = 0.25")
G_3d_0.5_s_plot <- ggplot(data = G_CvM_mat_s_df, aes(x = cop.short.names, y = as.numeric(X3_0.5))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 3)") + ylab("CvM")  + 
  ggtitle("tau = 0.5")
G_3d_0.75_s_plot <- ggplot(data = G_CvM_mat_s_df, aes(x = cop.short.names, y = as.numeric(X3_0.75))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 3)") + ylab("CvM")  + 
  ggtitle("tau = 0.75")

plot_grid(G_3d_0.25_s_plot, G_3d_0.5_s_plot, G_3d_0.75_s_plot, labels = "AUTO", ncol = 2)

G_4d_0.25_s_plot <- ggplot(data = G_CvM_mat_s_df, aes(x = cop.short.names, y = as.numeric(X4_0.25))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 4)") + ylab("CvM")  + 
  ggtitle("tau = 0.25")
G_4d_0.5_s_plot <- ggplot(data = G_CvM_mat_s_df, aes(x = cop.short.names, y = as.numeric(X4_0.5))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 4)") + ylab("CvM")  + 
  ggtitle("tau = 0.5")
G_4d_0.75_s_plot <- ggplot(data = G_CvM_mat_s_df, aes(x = cop.short.names, y = as.numeric(X4_0.75))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 4)") + ylab("CvM")  + 
  ggtitle("tau = 0.75")

plot_grid(G_4d_0.25_s_plot, G_4d_0.5_s_plot, G_4d_0.75_s_plot, labels = "AUTO", ncol = 2)

G_5d_0.25_s_plot <- ggplot(data = G_CvM_mat_s_df, aes(x = cop.short.names, y = as.numeric(X5_0.25))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 5)") + ylab("CvM")  + 
  ggtitle("tau = 0.25")
G_5d_0.5_s_plot <- ggplot(data = G_CvM_mat_s_df, aes(x = cop.short.names, y = as.numeric(X5_0.5))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 5)") + ylab("CvM")  + 
  ggtitle("tau = 0.5")
G_5d_0.75_s_plot <- ggplot(data = G_CvM_mat_s_df, aes(x = cop.short.names, y = as.numeric(X5_0.75))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 5)") + ylab("CvM")  + 
  ggtitle("tau = 0.75")

plot_grid(G_5d_0.25_s_plot, G_5d_0.5_s_plot, G_5d_0.75_s_plot, labels = "AUTO", ncol = 2)

C_2d_0.25_s_plot <- ggplot(data = C_CvM_mat_s_df, aes(x = cop.short.names, y = as.numeric(X2_0.25))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 2)") + ylab("CvM")  + 
  ggtitle("tau = 0.25")
C_2d_0.5_s_plot <- ggplot(data = C_CvM_mat_s_df, aes(x = cop.short.names, y = as.numeric(X2_0.5))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 2)") + ylab("CvM")  + 
  ggtitle("tau = 0.5")
C_2d_0.75_s_plot <- ggplot(data = C_CvM_mat_s_df, aes(x = cop.short.names, y = as.numeric(X2_0.75))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 2)") + ylab("CvM")  + 
  ggtitle("tau = 0.75")

plot_grid(C_2d_0.25_s_plot, C_2d_0.5_s_plot, C_2d_0.75_s_plot, labels = "AUTO", ncol = 2)

C_3d_0.25_s_plot <- ggplot(data = C_CvM_mat_s_df, aes(x = cop.short.names, y = as.numeric(X3_0.25))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 3)") + ylab("CvM")  + 
  ggtitle("tau = 0.25")
C_3d_0.5_s_plot <- ggplot(data = C_CvM_mat_s_df, aes(x = cop.short.names, y = as.numeric(X3_0.5))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 3)") + ylab("CvM")  + 
  ggtitle("tau = 0.5")
C_3d_0.75_s_plot <- ggplot(data = C_CvM_mat_s_df, aes(x = cop.short.names, y = as.numeric(X3_0.75))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 3)") + ylab("CvM")  + 
  ggtitle("tau = 0.75")

plot_grid(C_3d_0.25_s_plot, C_3d_0.5_s_plot, C_3d_0.75_s_plot, labels = "AUTO", ncol = 2)

C_4d_0.25_s_plot <- ggplot(data = C_CvM_mat_s_df, aes(x = cop.short.names, y = as.numeric(X4_0.25))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 4)") + ylab("CvM")  + 
  ggtitle("tau = 0.25")
C_4d_0.5_s_plot <- ggplot(data = C_CvM_mat_s_df, aes(x = cop.short.names, y = as.numeric(X4_0.5))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 4)") + ylab("CvM")  + 
  ggtitle("tau = 0.5")
C_4d_0.75_s_plot <- ggplot(data = C_CvM_mat_s_df, aes(x = cop.short.names, y = as.numeric(X4_0.75))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 4)") + ylab("CvM")  + 
  ggtitle("tau = 0.75")

plot_grid(C_4d_0.25_s_plot, C_4d_0.5_s_plot, C_4d_0.75_s_plot, labels = "AUTO", ncol = 2)

C_5d_0.25_s_plot <- ggplot(data = C_CvM_mat_s_df, aes(x = cop.short.names, y = as.numeric(X5_0.25))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 5)") + ylab("CvM")  + 
  ggtitle("tau = 0.25")
C_5d_0.5_s_plot <- ggplot(data = C_CvM_mat_s_df, aes(x = cop.short.names, y = as.numeric(X5_0.5))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 5)") + ylab("CvM")  + 
  ggtitle("tau = 0.5")
C_5d_0.75_s_plot <- ggplot(data = C_CvM_mat_s_df, aes(x = cop.short.names, y = as.numeric(X5_0.75))) +
  geom_bar(stat="identity") + xlab("Types of Empirical Copulas (d = 5)") + ylab("CvM")  + 
  ggtitle("tau = 0.75")

plot_grid(C_5d_0.25_s_plot, C_5d_0.5_s_plot, C_5d_0.75_s_plot, labels = "AUTO", ncol = 2)



