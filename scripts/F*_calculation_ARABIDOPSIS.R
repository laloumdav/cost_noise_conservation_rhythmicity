library(ggplot2)
library(reshape2)
library(officer)
library(magrittr)
library(cowplot)
library(Hmisc)
library(plyr)
## Stochastic gene expression (SGE) following Barroso et al. 2018 and Liu et al. 2019
## Raw code from Barroso available : https://figshare.com/articles/R_code_to_analyse_stochastic_gene_expression/4587169
## Code source from Liu et al. 2019 available from its GitHub

## Functions from Liu et al. 2019 :
#####* function used to calculate noise: global adjusted standard deviation *#####
predictSD <- function(expData) {
  ## polynomial model 
  expData$mean<-rowMeans(expData)
  expData$sd<-apply(expData[,-ncol(expData)],1,function(x) sd(x))
  m1 <- lm(sd~mean, expData)
  m2 <- update(m1, .~. + I(mean^2), expData)
  m3 <- update(m2, .~. + I(mean^3), expData)
  m4 <- update(m3, .~. + I(mean^4), expData)
  m5 <- update(m4, .~. + I(mean^5), expData)
  m6 <- update(m5, .~. + I(mean^6), expData)
  m7 <- update(m6, .~. + I(mean^7), expData)
  m8 <- update(m7, .~. + I(mean^8), expData)
  m9 <- update(m8, .~. + I(mean^9), expData)
  m10 <- update(m9, .~. + I(mean^10), expData)
  m11 <- update(m10, .~. + I(mean^11), expData)
  totalM<-list(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11)
  for (i in c(1:10)) {
    modelComp <- anova(totalM[[i]],totalM[[i+1]])
    if (is.na(modelComp$`Pr(>F)`[2]) | modelComp$`Pr(>F)`[2]>0.05) {
      m<-totalM[[i]]
      break
    }
  }
  return(predict(m))
}

noiseAJSD <- function(expData) {
  noise.adj.sd <- apply(expData, 1, FUN = sd) / predictSD(expData)
  return(noise.adj.sd)
}


## ## ## ## ## ## ## ## ## ## 
# From raw dataset
## ## ## ## ## ## ## ## ## ## 
data <- read.csv("~/Documents/cost_theory_workingspace/DATA/arabidopsis/root/single-cell/GSE46226_SC_Expression.csv", header = TRUE, fill=TRUE, row.names = "Locus")
# Deletes all-zero rows:
data <- data[-(which(rowSums(data) == 0)),]
# Evaluates if genes are expressed above a certain treshold (log(x+1) > 1.5)
fBigger <- function(x){
  any(log(x + 1) > 1.5);
}
appExpressed <- apply(data, 1, fBigger)
# Deletes corresponding rows
data <- data[appExpressed,]
subData <- data

# Get mean and variance of expression levels:
data$Mean <- rowMeans(subData)
data$Variance <- apply(subData, 1, var)
data$sd <- apply(subData, 1, sd)
# Get Liu et al. adjusted noise :
data$`noise adjusted sd` <- noiseAJSD(subData)
# Define which polynomial degree was used
## polynomial model 
m1 <- lm(sd~Mean, data)
m2 <- update(m1, .~. + I(Mean^2), data)
m3 <- update(m2, .~. + I(Mean^3), data)
m4 <- update(m3, .~. + I(Mean^4), data)
m5 <- update(m4, .~. + I(Mean^5), data)
m6 <- update(m5, .~. + I(Mean^6), data)
m7 <- update(m6, .~. + I(Mean^7), data)
m8 <- update(m7, .~. + I(Mean^8), data)
m9 <- update(m8, .~. + I(Mean^9), data)
m10 <- update(m9, .~. + I(Mean^10), data)
m11 <- update(m10, .~. + I(Mean^11), data)
totalM<-list(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11)
for (i in c(1:10)){
  modelComp <- anova(totalM[[i]],totalM[[i+1]])
  print(modelComp)
  if (is.na(modelComp$`Pr(>F)`[2]) | modelComp$`Pr(>F)`[2]>0.05) {
    print(paste("=> polynome of ", i, " degree. ANOVA pval = ", modelComp$`Pr(>F)`[2], sep=""))
    break
  }
}
# Model 1: sd ~ Mean + I(Mean^2) + I(Mean^3) + I(Mean^4) + I(Mean^5) + I(Mean^6) + 
#   I(Mean^7) + I(Mean^8) + I(Mean^9)
# Model 2: sd ~ Mean + I(Mean^2) + I(Mean^3) + I(Mean^4) + I(Mean^5) + I(Mean^6) + 
#   I(Mean^7) + I(Mean^8) + I(Mean^9) + I(Mean^10)
# Res.Df    RSS Df Sum of Sq     F Pr(>F)
# 1  13242 3324.0                          
# 2  13241 3323.6  1   0.37802 1.506 0.2198

# Gets Barroso et al. adjusted noise :
m1 <- lm(log(Variance)~log(Mean), data)
m2 <- update(m1, .~. + I(log(Mean)^2), data)
m3 <- update(m2, .~. + I(log(Mean)^3), data)
m4 <- update(m3, .~. + I(log(Mean)^4), data)
m5 <- update(m4, .~. + I(log(Mean)^5), data)

par(mfrow=c(2,2))
plot(log(data$Mean), log(data$Variance))
plot(data$Mean, data$Variance)
plot(data$Mean, data$sd)

# Defines F*
data$Fstar1 <- exp(log(data$Variance) - predict(m1))
data$Fstar2 <- exp(log(data$Variance) - predict(m2))
data$Fstar3 <- exp(log(data$Variance) - predict(m3))
data$Fstar4 <- exp(log(data$Variance) - predict(m4))
data$Fstar5 <- exp(log(data$Variance) - predict(m5))

cor.test(~Fstar1+Mean, data, method = "kendall") # p-value = 0.7805
cor.test(~Fstar2+Mean, data, method = "kendall") # p-value < 2.2e-16
cor.test(~Fstar3+Mean, data, method = "kendall") # p-value = 4.104e-09
cor.test(~Fstar4+Mean, data, method = "kendall") # p-value = 0.116
cor.test(~Fstar5+Mean, data, method = "kendall") # p-value = 0.0693
# Degree 1 or Degree 4 capture most the effect of mean

# Other measures of stochasticity available in the literature:
# Computing fano factor
fFano <- function(x){
  if (all(x==0)) return(NA)
  else return (var(x) / mean(x))
}

# coefficient of variation squared
squaredCV <- function(x) {
  if (all(x==0)) return(NA)
  else return (var(x) / (mean(x))^2)
}

data$'Fano Factor' <- apply(subData, 1, fFano)
data$'Coeff Variation squared' <- apply(subData, 1, squaredCV)
# Distance to median 
library(scran)
data$'distance to median' <- DM(mean = data$Mean, cv2 = data$'Coeff Variation squared', win.size=51)

dat <- data[, c("Mean", "Variance", "sd", "noise adjusted sd", "Fano Factor", "Fstar1", "Fstar2", "Fstar3", "Fstar4", "Fstar5", "Coeff Variation squared", "distance to median")]
#dat <- data.frame("ID"=rownames(dat), dat)
#write.table(dat, "~/Documents/cost_theory_workingspace/DATA/arabidopsis/root/single-cell/noise_data.txt", sep = "\t", quote = F, row.names = F)


#data <- reshape::rename(data, replace = c("noise.adj.sd"="noise adjusted sd", "FanoFactor"="Fano Factor", "Fstar1"="F*_1", "Fstar4"="F*_4"))
dat <- data[, c("Mean", "Variance", "sd", "noise adjusted sd", "Fano Factor", "Fstar1", "Fstar4", "Coeff Variation squared", "distance to median")]

# Measure the slope of the lm
slopeAndRsquared <- function(x) {
  lmtest <- lm(formula = dat$Mean ~ x)
  return(c(signif(lmtest$coef[[2]], 3), signif(summary(lmtest)$r.squared, 3)))
}
slope.and.rsquared <- data.frame(t(apply(dat[, -1], 2, FUN = slopeAndRsquared)))
colnames(slope.and.rsquared) <- c("slope", "Rsquared")
slope.and.rsquared$SGEtype <- rownames(slope.and.rsquared)
# Re-order the plots
slope.and.rsquared$SGEtype <- factor(slope.and.rsquared$SGEtype, 
                                     levels=c("Variance", "sd", "noise adjusted sd", "Fano Factor", "Fstar1", "Fstar4", "Coeff Variation squared", "distance to median"))


dat <- reshape::rename(dat, replace = c("Variance" = paste("Variance", "\n", "slope = ", slope.values[slope.values$SGEtype=="Variance", "slope.values"], sep=""), 
                                        "sd" = paste("sd", "\n", "slope = ", slope.values[slope.values$SGEtype=="sd", "slope.and.rsquared"], sep=""),
                                        "noise adjusted sd" = paste("noise adjusted sd", "\n", "slope = ", slope.values[slope.values$SGEtype=="noise adjusted sd", "slope.values"], sep=""),
                                        "Fano Factor" = paste("Fano Factor", "\n", "slope = ", slope.values[slope.values$SGEtype=="Fano Factor", "slope.values"], sep=""),
                                        "Fstar1" = paste("Fstar1", "\n", "slope = ", slope.values[slope.values$SGEtype=="Fstar1", "slope.values"], sep=""),
                                        "Fstar4" = paste("Fstar4", "\n", "slope = ", slope.values[slope.values$SGEtype=="Fstar4", "slope.values"], sep=""),
                                        "Coeff Variation squared" = paste("Coeff Variation squared", "\n", "slope = ", slope.values[slope.values$SGEtype=="Coeff Variation squared", "slope.values"], sep=""),
                                        "distance to median" = paste("distance to median", "\n", "slope = ", slope.values[slope.values$SGEtype=="distance to median", "slope.values"], sep="")))

# Transform everything in the same scale
#scalingFunction <- function(x) {
#  if (min(x)<0) {
#    x = x + abs(min(x, na.rm = TRUE))
#  }
#  y = x / max(x, na.rm = TRUE)
#  return(y)
#}
#dat[, -1] <- apply(dat[, -1], 2, FUN = scalingFunction)

dat <- melt(data = dat,
            id.vars = "Mean", variable.name = "SGEtype", value.name = "SGEvalue")

p1 <- ggplot(dat, aes(x = Mean, y=SGEvalue)) + 
  geom_point(shape = ".", color = "grey50") + 
  geom_smooth(method = "lm", color = "black", se = FALSE) + 
  facet_wrap(~SGEtype, scales = "free_y", ncol = 2) +
  #scale_x_log10() + scale_y_log10() + 
  theme_bw() +
  theme(axis.title.y = element_blank())
p1

p2bis <- ggplot(dat, aes(SGEvalue)) + 
  geom_histogram(aes(y=..density..), bins = 100, color="grey60", fill="grey60") + 
  facet_wrap(~SGEtype, scales = "free", ncol = 2) + 
  #scale_x_log10() +
  labs(caption="log scaled") +
  theme_bw() 
p2bis

# Alternative representation with discretization:
dat$MeanClass <- cut2(dat$Mean, m = 3000, levels.mean = TRUE)
dat2 <- ddply(.data = dat, .variables = c("MeanClass", "SGEtype"),
              .fun = plyr::summarize,
              median = median(SGEvalue),
              lo95 = quantile(SGEvalue, prob=0.025),
              hi95 = quantile(SGEvalue, prob=0.975))

p2 <- ggplot(dat2, aes(x = as.numeric(as.character(MeanClass)), y=median)) + 
  geom_quantile(data=dat, aes(x=Mean, y=SGEvalue), col = "grey50") + 
  geom_point() + 
  geom_errorbar(aes(ymin = lo95, ymax = hi95)) + 
  facet_wrap(~SGEtype, scales = "free_y", ncol = 2) + 
  scale_x_log10() + #scale_y_log10() + 
  #labs(caption="log scaled") +
  theme_bw() + 
  ylab("SGE measure") + xlab("Mean expression")+
  theme(#aspect.ratio = .9,
    axis.text.y = element_text(size=10), 
    axis.title.y = element_text(size=12), 
    axis.text.x = element_text(size=10), 
    axis.title.x = element_text(size=12), 
    panel.background = element_rect(fill = "grey97", colour = "grey90"), 
    axis.line = element_line(colour = "grey32"),
    panel.spacing = unit(.03, "lines"),
    legend.position = "top",
    plot.caption = element_text(size=12, face = "italic", hjust=0),
    strip.text = element_text(size=11))
p2

# Export as pdf image
#pdf("~/Documents/cost_theory_workingspace/supplementary_information/noise_geomQuantiles_arabidopsis.pdf", height = 6, width = 4.5)
#print(p2)
#dev.off()

##############
# Add the plots to the word document:
read_docx(path = "~/Documents/cost_theory_workingspace/figures_noise.docx") %>%
  cursor_begin() %>% 
  
  cursor_reach(keyword = "##noise_estimations_variance_and_sd_ARABIDOPSIS_Root##") %>%
  body_add_gg(value = p2, width = 5, height = 6) %>%
  cursor_forward() %>%
  body_remove()  %>%
  
  cursor_reach(keyword = "##noise_histograms_variance_and_sd_ARABIDOPSIS_Root##") %>%
  body_add_gg(value = p2bis, width = 5, height = 6) %>%
  cursor_forward() %>%
  body_remove()  %>%
  
  cursor_end() %>%
  print(target = "~/Documents/cost_theory_workingspace/figures_noise.docx")














# Distribution of noise.adj.sd
p2a <- ggplot(data, aes(`noise adjusted sd`)) + 
  geom_histogram(aes(y=..density..), bins = 100, color="grey60", fill="grey60") + 
  #geom_density(alpha=.5, fill="grey30", size = 1) + 
  geom_vline(xintercept = 1, size = 1, linetype = "dotted") + 
  theme_bw() + 
  xlab("noise adjusted sd")
# Distribution of distance to median
p2b <- ggplot(data, aes(`distance to median`)) + 
  geom_histogram(aes(y=..density..), bins = 100, color="grey60", fill="grey60") + 
  #geom_density(alpha=.5, fill="grey30", size = 1) + 
  geom_vline(xintercept = 0, size = 1, linetype = "dotted") + 
  theme_bw() + 
  xlab("distance to median")
# Distribution of F* with first polynomial degree
p2c <- ggplot(data, aes(`F*_1`)) + 
  geom_histogram(aes(y=..density..), bins = 100, color="grey60", fill="grey60") + 
  #geom_density(alpha=.5, fill="grey30", size = 1) + 
  geom_vline(xintercept = 1, size = 1, linetype = "dotted") + 
  theme_bw() + 
  xlab("F* = F*_1")
# Distribution of F* with 4 polynomial degree
p2d <- ggplot(data, aes(`F*_4`)) + 
  geom_histogram(aes(y=..density..), bins = 100, color="grey60", fill="grey60") + 
  #geom_density(alpha=.5, fill="grey30", size = 1) + 
  geom_vline(xintercept = 1, size = 1, linetype = "dotted") + 
  theme_bw() + 
  xlab("F* = F*_4")

# Combine them:
pp <- plot_grid(p2a, p2b, p2c, p2d, labels = c("b", "c", "d", "e"), ncol = 1)
pp <- plot_grid(p1, pp, labels = c("a", ""), rel_widths = c(.7, .3))

##############
# Add the plots to the word document:
read_docx(path = "~/Documents/cost_theory_workingspace/figures_rmd.docx") %>%
  cursor_begin() %>% 
  
  cursor_reach(keyword = "##noise_estimations_variance_and_sd_ARABIDOPSIS##") %>%
  body_add_gg(value = pp, width = 7, height = 7.5) %>%
  cursor_forward() %>%
  body_remove()  %>%
  
  cursor_end() %>%
  print(target = "~/Documents/cost_theory_workingspace/figures_rmd.docx")




# Alternative representation with discretization:
dat$MeanClass <- cut2(dat$Mean, m = 3000, levels.mean = TRUE)
dat2 <- ddply(.data = dat, .variables = c("MeanClass", "SGEtype"),
              .fun = plyr::summarize,
              median = median(SGEvalue),
              lo95 = quantile(SGEvalue, prob=0.025),
              hi95 = quantile(SGEvalue, prob=0.975))

p1a2 <- ggplot(dat2, aes(x = as.numeric(as.character(MeanClass)), y=median)) + 
  geom_quantile(data=dat, aes(x=Mean, y=SGEvalue), col = "grey50") + 
  geom_point() + 
  geom_errorbar(aes(ymin = lo95, ymax = hi95)) + 
  facet_wrap(~SGEtype, scales = "free_y", ncol = 2) + 
  scale_x_log10() + scale_y_log10() + 
  theme_bw() + 
  ylab("SGE measure") + xlab("Mean expression")
p1a2




