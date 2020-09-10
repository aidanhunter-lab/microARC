
# Standardise nutrient data using linear mixed models

library('nlme')

getwd()
dat = read.csv('cleanedData.csv')

head(dat)
dat$scaled_Value = NA
dat$Event = as.factor(dat$Event)

# Data type
levels(dat$Variable)

for(dv in c('PON','POC')){
  # Model POM with an exponential decay over depth
  
  d = dat[dat$Variable==dv,]

  d$log_Value = log(d$Value)
  # create grouped data structure
  gd = groupedData(formula = log_Value ~ Depth | Event, data = d, order.groups = FALSE, FUN = mean)
  
  par(mfrow=c(1,1))
  plot(gd, ylab = paste0('ln(',dv,')'))
  par(mfrow=c(1,2))
  plot(log_Value ~ Depth, data = gd, ylab = paste0('ln(',dv,')'))
  mu = mean(gd$log_Value)
  abline(h=mu)
  plot(log_Value ~ Event, data = gd, ylab = paste0('ln(',dv,')'))
  abline(h=mu)
  
  # simplest model: single intercept
  mod1 = tryCatch(lm(log_Value ~ 1, data = gd), error = function(e) NULL)
  if(!is.null(mod1)){
    summary(mod1)
    par(mfrow = c(1,2))
    plot(resid(mod1) ~ Event, data = gd, horizontal = T, ylab = 'residuals')
    abline(v=0)
    plot(Depth ~ resid(mod1), data = gd, xlab = 'residuals')
    abline(v=0)
  }
  
  # include depth as covariate
  mod2 = tryCatch(lm(log_Value ~ Depth, data = gd), error = function(e) NULL)
  if(!is.null(mod2)){
    summary(mod2)
    par(mfrow = c(1,2))
    plot(resid(mod2) ~ Event, data = gd, horizontal = T, ylab = 'residuals')
    abline(v=0)
    plot(Depth ~ resid(mod2), data = gd, xlab = 'residuals')
    abline(v=0)
  }
  
  # random event-effect on intercept
  mod3 = tryCatch(lme(fixed = log_Value ~ Depth, data = gd, random = ~ 1 | Event, method = 'REML'), error = function(e) NULL)
  if(!is.null(mod3)){
    summary(mod3)
    par(mfrow = c(1,2))
    plot(resid(mod3) ~ Event, data = gd, horizontal = T, ylab = 'residuals')
    abline(v=0)
    plot(Depth ~ resid(mod3), data = gd, xlab = 'residuals')
    abline(v=0)
  }
  
  # random event-effect on gradient
  mod4 = tryCatch(lme(fixed = log_Value ~ Depth, data = gd, random = ~ -1+Depth | Event, method = 'REML'), error = function(e) NULL)
  if(!is.null(mod4)){
    summary(mod4)
    par(mfrow = c(1,2))
    plot(resid(mod4) ~ Event, data = gd, horizontal = T, ylab = 'residuals')
    abline(v=0)
    plot(Depth ~ resid(mod4), data = gd, xlab = 'residuals')
    abline(v=0)
  }
  
  mod5 = tryCatch(lme(fixed = log_Value ~ Depth, data = gd, random = ~ Depth | Event, method = 'ML'), error = function(e) NULL)
  if(!is.null(mod5)){
    summary(mod5)
    par(mfrow = c(1,2))
    plot(resid(mod5) ~ Event, data = gd, horizontal = T, ylab = 'residuals')
    abline(v=0)
    plot(Depth ~ resid(mod5), data = gd, xlab = 'residuals')
    abline(v=0)
  }
  
  # Use AIC to choose between models
  test_aic = rep(NA,5)
  for(i in 1:5){
    im = get(paste0('mod',i))
    if(!is.null(im)){
      test_aic[i]=AIC(im)
    }
  }

  whichMod = paste0('mod', which.min(test_aic))
  mod = get(whichMod)
  rm(list=c('mod1','mod2','mod3','mod4','mod5'))

  # Extract fitted coefficients and residuals
  cf = mod$coefficients$fixed
  cr = mod$coefficients$random$Event
  rsd = resid(mod)
  stdDev = sd(rsd)
  # stdDev = mod$sigma
  
  # Scale the data
  indEv = as.numeric(as.character(gd$Event))
  indEv = indEv - min(indEv) + 1
  
  # if(whichMod == 'mod1'){
  #   mu = {cf[1] + cr[,1][indEv]} + cf[2] * gd$Depth
  # }
  # if(whichMod == 'mod2'){
  #   mu = {cf[1] + cr[,1][indEv]} + cf[2] * gd$Depth
  # }
  
  if(whichMod == 'mod3'){
    mu = {cf[1] + cr[,1][indEv]} + cf[2] * gd$Depth
  }
  if(whichMod == 'mod4'){
    mu = cf[1] + {cf[2] + cr[,2][indEv]} * gd$Depth
  }
  if(whichMod == 'mod5'){
    mu = {cf[1] + cr[,1][indEv]} + {cf[2] + cr[,2][indEv]} * gd$Depth
  }
  
  ys = {gd$log_Value - mu} / stdDev
  
  gd$scaled_Value = ys
  
  # Compare raw data to scaled data
  par(mfcol=c(2,2), mar = c(5,5,1,1), oma = c(0,0,1.5,0))
  
  xlim = range(gd$Depth)
  ylim = range(gd$log_Value)
  xlim_ = xlim + 0.05 * c(-1,1) * diff(xlim)
  plot(NULL,xlim=xlim,ylim=ylim, type='n',xlab = 'depth (m)', ylab = expression(ln(PON) ~ (mmol ~ N ~ m^{-3})))
  mu = mean(gd$log_Value)
  s = sd(gd$log_Value)
  polygon(c(xlim_[1],xlim_[2],xlim_[2],xlim_[1],xlim_[1]), c(mu-s,mu-s,mu+s,mu+s,mu-1), col = gray(0.8))
  abline(h=mu)
  points(log_Value ~ Depth, data = gd)
  box()
  
  xlim = range(gd$Depth)
  ylim = range(gd$scaled_Value)
  xlim_ = xlim + 0.05 * c(-1,1) * diff(xlim)
  plot(NULL,xlim=xlim,ylim=ylim, type='n',xlab = 'depth (m)', ylab = 'scaled data')
  mu = mean(gd$scaled_Value)
  s = sd(gd$scaled_Value)
  polygon(c(xlim_[1],xlim_[2],xlim_[2],xlim_[1],xlim_[1]), c(mu-s,mu-s,mu+s,mu+s,mu-1), col = gray(0.8))
  abline(h=mu)
  points(scaled_Value ~ Depth, data = gd)
  box()
  
  plot(log_Value ~ Event, data = gd, xlab = 'sampling event', ylab = expression(ln(PON) ~ (mmol ~ N ~ m^{-3})))
  mu = mean(gd$log_Value)
  abline(h=mu)
  
  plot(scaled_Value ~ Event, data = gd, xlab = 'sampling event', ylab = 'scaled data')
  mu = mean(gd$scaled_Value)
  abline(h=mu)
  mtext(text=paste0('Compare raw and standardised ',dv,' measurements'),side=3,outer=TRUE)
  
  
  dat$scaled_Value[dat$Variable==dv] = gd$scaled_Value
  
}




# Inorganic nitrogen modelled as linear over ln(depth)
dv = 'N'
d = dat[dat$Variable==dv,]

d$log_Depth = log(d$Depth)
# create grouped data structure
gd = groupedData(formula = Value ~ log_Depth | Event, data = d, order.groups = FALSE, FUN = mean)

plot(gd)
par(mfrow=c(1,2))
plot(Value ~ log_Depth, data = gd, ylab = dv)
mu = mean(gd$Value)
abline(h=mu)
plot(Value ~ Event, data = gd, ylab = dv)
abline(h=mu)

# simplest model: single intercept
mod1 = tryCatch(lm(Value ~ 1, data = gd), error = function(e) NULL)
if(!is.null(mod1)){
  summary(mod1)
  par(mfrow = c(1,2),oma=c(0,0,1,0))
  plot(resid(mod1) ~ Event, data = gd, horizontal = T, ylab = 'residuals')
  abline(v=0)
  plot(Depth ~ resid(mod1), data = gd, xlab = 'residuals')
  abline(v=0)
  mtext(text='model 1',side=3,line=0,outer=T)
}

# include log(depth) as covariate
mod2 = tryCatch(lm(Value ~ log_Depth, data = gd), error = function(e) NULL)
if(!is.null(mod2)){
  summary(mod2)
  par(mfrow = c(1,2),oma=c(0,0,1,0))
  plot(resid(mod2) ~ Event, data = gd, horizontal = T, ylab = 'residuals')
  abline(v=0)
  plot(Depth ~ resid(mod2), data = gd, xlab = 'residuals')
  abline(v=0)
  mtext(text='model 2',side=3,line=0,outer=T)
}

# random event-effect on intercept
mod3 = tryCatch(lme(fixed = Value ~ log_Depth, data = gd, random = ~ 1 | Event, method = 'REML'), error = function(e) NULL)
if(!is.null(mod3)){
  summary(mod3)
  par(mfrow = c(1,2),oma=c(0,0,1,0))
  plot(resid(mod3) ~ Event, data = gd, horizontal = T, ylab = 'residuals')
  abline(v=0)
  plot(Depth ~ resid(mod3), data = gd, xlab = 'residuals')
  abline(v=0)
  mtext(text='model 3',side=3,line=0,outer=T)
}

# random event-effect on gradient
mod4 = tryCatch(lme(fixed = Value ~ log_Depth, data = gd, random = ~ -1 + log_Depth | Event, method = 'REML'), error = function(e) NULL)
if(!is.null(mod4)){
  summary(mod4)
  par(mfrow = c(1,2),oma=c(0,0,1,0))
  plot(resid(mod4) ~ Event, data = gd, horizontal = T, ylab = 'residuals')
  abline(v=0)
  plot(Depth ~ resid(mod4), data = gd, xlab = 'residuals')
  abline(v=0)
  mtext(text='model 4',side=3,line=0,outer=T)
}

# random event-effect on intercept and gradient
mod5 = tryCatch(lme(fixed = Value ~ log_Depth, data = gd, random = ~ log_Depth | Event, method = 'REML'), error = function(e) NULL)
if(!is.null(mod5)){
  summary(mod5)
  par(mfrow = c(1,2),oma=c(0,0,1,0))
  plot(resid(mod5) ~ Event, data = gd, horizontal = T, ylab = 'residuals')
  abline(v=0)
  plot(Depth ~ resid(mod5), data = gd, xlab = 'residuals')
  abline(v=0)
  mtext(text='model 5',side=3,line=0,outer=T)
}

# Use AIC to choose between models
test_aic = rep(NA,5)
for(i in 1:5){
  im = get(paste0('mod',i))
  if(!is.null(im)){
    test_aic[i]=AIC(im)
  }
}

whichMod = paste0('mod', which.min(test_aic))
mod = get(whichMod)
rm(list=c('mod1','mod2','mod3','mod4','mod5'))

# Extract fitted coefficients and residuals
cf = mod$coefficients$fixed
cr = mod$coefficients$random$Event
rsd = resid(mod)
stdDev = sd(rsd)

# Scale the data
indEv = as.numeric(as.character(gd$Event))
indEv = indEv - min(indEv) + 1

if(whichMod=='mod3'){
  mu = {cf[1] + cr[,1][indEv]} + cf[2] * gd$log_Depth
}
if(whichMod=='mod4'){
  mu = cf[1] + {cf[2] + cr[,2][indEv]} * gd$log_Depth
}
if(whichMod=='mod5'){
  mu = {cf[1] + cr[,1][indEv]} + {cf[2] + cr[,2][indEv]} * gd$log_Depth
}

ys = {gd$Value - mu} / stdDev

gd$scaled_Value = ys

# Compare raw data to scaled data
par(mfcol=c(2,2), mar = c(5,5,1,1), oma = c(0,0,1.5,0))

xlim = range(gd$log_Depth)
ylim = range(gd$Value)
xlim_ = xlim + 0.05 * c(-1,1) * diff(xlim)
plot(NULL,xlim=xlim,ylim=ylim, type='n',xlab = 'ln(depth)', ylab = expression(NO[2] + NO[3] ~ (mmol ~ N ~ m^{-3})))
mu = mean(gd$Value)
s = sd(gd$Value)
polygon(c(xlim_[1],xlim_[2],xlim_[2],xlim_[1],xlim_[1]), c(mu-s,mu-s,mu+s,mu+s,mu-1), col = gray(0.8))
abline(h=mu)
points(Value ~ log_Depth, data = gd)
box()

xlim = range(gd$log_Depth)
ylim = range(gd$scaled_Value)
xlim_ = xlim + 0.05 * c(-1,1) * diff(xlim)
plot(NULL,xlim=xlim,ylim=ylim, type='n',xlab = 'ln(depth)', ylab = 'scaled data')
mu = mean(gd$scaled_Value)
s = sd(gd$scaled_Value)
polygon(c(xlim_[1],xlim_[2],xlim_[2],xlim_[1],xlim_[1]), c(mu-s,mu-s,mu+s,mu+s,mu-1), col = gray(0.8))
abline(h=mu)
points(scaled_Value ~ log_Depth, data = gd)
box()

plot(Value ~ Event, data = gd, xlab = 'sampling event', ylab = expression(NO[2] + NO[3] ~ (mmol ~ N ~ m^{-3})))
mu = mean(gd$Value)
abline(h=mu)

plot(scaled_Value ~ Event, data = gd, xlab = 'sampling event', ylab = 'scaled data')
mu = mean(gd$scaled_Value)
abline(h=mu)
mtext(text='Compare raw and standardised N measurements',side=3,outer=TRUE)


dat$scaled_Value[dat$Variable==dv] = gd$scaled_Value




# Plot scaled data of each type side-by-side
par(mfrow=c(1,1))
nbreaks=12
hist_PON = hist(dat$scaled_Value[dat$Variable=='PON'], freq=FALSE, plot=FALSE, breaks=nbreaks)
hist_POC = hist(dat$scaled_Value[dat$Variable=='POC'], freq=FALSE, plot=FALSE, breaks=nbreaks)
hist_N = hist(dat$scaled_Value[dat$Variable=='N'], freq=FALSE, plot=FALSE, breaks=nbreaks)
xl = range(c(hist_PON$breaks,hist_POC$breaks,hist_N$breaks))
yl = c(0,max(c(hist_PON$density,hist_POC$density,hist_N$density)))

hist(dat$scaled_Value[dat$Variable=='PON'], freq=FALSE, breaks=nbreaks, xlim = xl, ylim = yl, xlab = 'scaled data', main = 'Scaled data distributions', col = adjustcolor('green', alpha.f = 0.5))
hist(dat$scaled_Value[dat$Variable=='POC'], freq=FALSE, breaks=nbreaks, add = TRUE, col = adjustcolor('red', alpha.f = 0.5))
hist(dat$scaled_Value[dat$Variable=='N'], freq=FALSE, breaks=nbreaks, add = TRUE, col = adjustcolor('blue', alpha.f = 0.5))
box()
legend('topleft', c('PON','POC','N'), fill = adjustcolor(c('green','red','blue'),alpha.f=0.5), bty = 'n')




par(mfrow=c(1,1))
plot_dat = groupedData(formula = scaled_Value ~ 1 | Event, data = dat, order.groups = FALSE, FUN = mean)
plot(plot_dat, inner = ~ Variable, xlab = 'standardised data', ylab = 'sampling event')

plot_dat = groupedData(formula = scaled_Value ~ Depth | Variable, data = dat, order.groups = FALSE, FUN = mean)
# plot(plot_dat)
plot(plot_dat, inner = ~ Event)

plot_dat = groupedData(formula = scaled_Value ~ Depth | Event, data = dat, order.groups = FALSE, FUN = mean)
# plot(plot_dat)
plot(plot_dat, inner = ~ Variable, ylab = 'standardised data')




