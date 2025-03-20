library("install.load")
install_load("tidyverse", "data.table", "lubridate")

# load selected phq and weather data
phq_weather <- read.csv("/Users/yuezhou/Documents/Work/Data/MDD_Weather/processed_data/analysis_data/phq_weather_selected.csv") %>%
  select(p_id, date, sum_phq8, temp_mean, pressure_mean, humidity_mean, clouds_mean, wind_speed_mean, sun_diff)
colnames(phq_weather)[2] <- "time_str"

# subtypes of weather
summary_coeff <- phq_weather %>%
  dplyr::group_by(p_id) %>%
  dplyr::summarise(
                   temp = cor.test(temp_mean,sum_phq8,method = "spearman")[["estimate"]][["rho"]],
                   pressure= cor.test(pressure_mean,sum_phq8,method = "spearman")[["estimate"]][["rho"]],
                   humidity= cor.test(humidity_mean,sum_phq8,method = "spearman")[["estimate"]][["rho"]],
                   clouds= cor.test(clouds_mean,sum_phq8,method = "spearman")[["estimate"]][["rho"]],
                   wind_speed= cor.test(wind_speed_mean,sum_phq8,method = "spearman")[["estimate"]][["rho"]],
                   sun_diff= cor.test(sun_diff,sum_phq8,method = "spearman")[["estimate"]][["rho"]]
                   
                   
  )
summary_coeff = na.omit(summary_coeff)

# load demographics
demo <- read.csv("/Users/yuezhou/Documents/Work/Data/MDD_Weather/processed_data/analysis_data/demo.csv") %>% select("p_id", "age", "sex", "site","marry_align","edu_age_align", "income_align", "employment")
demo$employment[demo$employment=="Not reported"] = "No"
# load step data
step_feature <- read.csv("/Users/yuezhou/Documents/Work/Data/MDD_Weather/processed_data/analysis_data/phq8_features.csv") %>% select(p_id,phq8,time_str,step_all_mean)
colnames(step_feature)[4] <- "steps"

# merge PHQ and step
phq_weather <- merge(phq_weather, step_feature, by = c("p_id","time_str"))
phq_weather <- na.omit(phq_weather) # remove na
# merge demo and PHQ 
phq_weather <- merge(phq_weather,demo,all.x = TRUE)

#add a confounding variable for National lockdown
phq_weather$time_str = as_datetime(phq_weather$time_str)

phq_weather$lockdown = 0
phq_weather$lockdown[phq_weather$site == "kcl" & phq_weather$time_str >= as.Date("2020-03-23") & phq_weather$time_str < as.Date("2020-06-01")] = 1
phq_weather$lockdown[phq_weather$site == "kcl" & phq_weather$time_str >= as.Date("2020-11-05") & phq_weather$time_str < as.Date("2020-12-02")] = 1
phq_weather$lockdown[phq_weather$site == "kcl" & phq_weather$time_str >= as.Date("2021-01-04") & phq_weather$time_str < as.Date("2021-03-08")] = 1

phq_weather$lockdown[phq_weather$site == "ciber" & phq_weather$time_str >= as.Date("2020-03-14") & phq_weather$time_str < as.Date("2020-06-21")] = 1
phq_weather$lockdown[phq_weather$site == "ciber" & phq_weather$time_str >= as.Date("2020-10-25") & phq_weather$time_str < as.Date("2020-12-09")] = 1
phq_weather$lockdown[phq_weather$site == "ciber" & phq_weather$time_str >= as.Date("2021-01-31") & phq_weather$time_str < as.Date("2021-05-09")] = 1

phq_weather$lockdown[phq_weather$site == "vumc" & phq_weather$time_str >= as.Date("2020-03-15") & phq_weather$time_str < as.Date("2020-05-11")] = 1
phq_weather$lockdown[phq_weather$site == "vumc" & phq_weather$time_str >= as.Date("2020-12-23") & phq_weather$time_str < as.Date("2021-02-09")] = 1
phq_weather$lockdown = as.factor(phq_weather$lockdown)



################################################################################
######################## Mediation Analysis ####################################
################################################################################

library(lme4)
library(mediation)

# divided steps by 100
#phq_weather$steps = phq_weather$steps/100
data <- phq_weather

## PHQ-8 is Mediator
# X = temp, M=phq8, Y=step 
temp_phq8_step <- function(data){
  # model m: M = aX 
  model_m <- lmer(sum_phq8 ~ temp_mean + age+sex+site + marry_align+ edu_age_align + income_align+employment+ (1 | p_id), data = data)
  
  # model y: Y = aX + bM 
  model_y <- lmer(steps~ temp_mean + sum_phq8 + age+sex+site + marry_align+ edu_age_align + income_align+employment+  (1 | p_id), data = data)
  
  # mediation analysis
  med.out <- mediate(model_m, model_y, 
                     treat = "temp_mean", 
                     mediator = "sum_phq8",
                     sims = 500) 
  
  y <- summary(med.out)
  
  result_mat <- matrix(nrow = 5, ncol = 3)
  result_mat = as.data.frame(result_mat)
  colnames(result_mat) <- c("Path", "Est", "P")
  
  result_mat$Path[1] = "Total Effects"
  result_mat$Path[2] = "Direct Effects"
  result_mat$Path[3] = "Indirect Effects"
  result_mat$Path[4] = "a path"
  result_mat$Path[5] = "b path"
  
  
  result_mat$Est[1] =paste(round(y$tau.coef,3),'[', round(y$tau.ci[['2.5%']], 3),',',round(y$tau.ci[['97.5%']],3),']',sep = "")
  result_mat$P[1] = round(y$tau.p, 3)
  
  result_mat$Est[2] =paste(round(y$z0,3),'[', round(y$z0.ci[['2.5%']], 3),',',round(y$z0.ci[['97.5%']],3),']',sep = "")
  result_mat$P[2] = round(y$z0.p, 3)
  
  result_mat$Est[3] =paste(round(y$d0,3),'[', round(y$d0.ci[['2.5%']], 3),',',round(y$d0.ci[['97.5%']],3),']',sep = "")
  result_mat$P[3] = round(y$d0.p, 3)
  
  ci_a = confint(model_m)
  summary_m <- summary(model_m)
  coef_m <- summary_m$coefficients
  p_values_m <- 2 * (1 - pnorm(abs(coef_m[, "t value"])))
  
  result_mat$Est[4] =paste(round(coef_m[2,1],3),'[', round(ci_a[4,][[1]], 3),',',round(ci_a[4,][[2]],3),']',sep = "")
  result_mat$P[4] = round(p_values_m[[2]], 3)
  
  
  ci_b = confint(model_y)
  summary_y <- summary(model_y)
  coef_y <- summary_y$coefficients
  p_values_y <- 2 * (1 - pnorm(abs(coef_y[, "t value"])))
  
  result_mat$Est[5] =paste(round(coef_y[3,1],3),'[', round(ci_b[5,][[1]], 3),',',round(ci_b[5,][[2]],3),']',sep = "")
  result_mat$P[5] = round(p_values_y[[3]], 3)
  return(result_mat)
}
results_1 = temp_phq8_step(phq_weather)
results_2 = temp_phq8_step(phq_weather[phq_weather$p_id %in% summary_coeff$p_id[summary_coeff$temp >= 0.2],])
results_3 = temp_phq8_step(phq_weather[phq_weather$p_id %in% summary_coeff$p_id[summary_coeff$temp <= -0.2],])
results_4 = temp_phq8_step(phq_weather[phq_weather$p_id %in% summary_coeff$p_id[summary_coeff$temp > -0.2 & summary_coeff$temp < 0.2],])
results_combined_all = cbind(results_1, results_2[,2:3],results_3[,2:3],results_4[,2:3])

# X = humidity, M=phq8, Y=step 
humidity_phq8_step <- function(data){
  # model m: M = aX 
  model_m <- lmer(sum_phq8 ~ humidity_mean + age+sex+site + marry_align+ edu_age_align + income_align+employment+ (1 | p_id), data = data)
  
  # model y: Y = aX + bM 
  model_y <- lmer(steps~ humidity_mean + sum_phq8 + age+sex+site + marry_align+ edu_age_align + income_align+employment+  (1 | p_id), data = data)
  
  # mediation analysis
  med.out <- mediate(model_m, model_y, 
                     treat = "humidity_mean", 
                     mediator = "sum_phq8",
                     sims = 500) 
  
  y <- summary(med.out)
  result_mat <- matrix(nrow = 5, ncol = 3)
  result_mat = as.data.frame(result_mat)
  colnames(result_mat) <- c("Path", "Est", "P")
  
  result_mat$Path[1] = "Total Effects"
  result_mat$Path[2] = "Direct Effects"
  result_mat$Path[3] = "Indirect Effects"
  result_mat$Path[4] = "a path"
  result_mat$Path[5] = "b path"
  
  
  result_mat$Est[1] =paste(round(y$tau.coef,3),'[', round(y$tau.ci[['2.5%']], 3),',',round(y$tau.ci[['97.5%']],3),']',sep = "")
  result_mat$P[1] = round(y$tau.p, 3)
  
  result_mat$Est[2] =paste(round(y$z0,3),'[', round(y$z0.ci[['2.5%']], 3),',',round(y$z0.ci[['97.5%']],3),']',sep = "")
  result_mat$P[2] = round(y$z0.p, 3)
  
  result_mat$Est[3] =paste(round(y$d0,3),'[', round(y$d0.ci[['2.5%']], 3),',',round(y$d0.ci[['97.5%']],3),']',sep = "")
  result_mat$P[3] = round(y$d0.p, 3)
  
  ci_a = confint(model_m)
  summary_m <- summary(model_m)
  coef_m <- summary_m$coefficients
  p_values_m <- 2 * (1 - pnorm(abs(coef_m[, "t value"])))
  
  result_mat$Est[4] =paste(round(coef_m[2,1],3),'[', round(ci_a[4,][[1]], 3),',',round(ci_a[4,][[2]],3),']',sep = "")
  result_mat$P[4] = round(p_values_m[[2]], 3)
  
  
  ci_b = confint(model_y)
  summary_y <- summary(model_y)
  coef_y <- summary_y$coefficients
  p_values_y <- 2 * (1 - pnorm(abs(coef_y[, "t value"])))
  
  result_mat$Est[5] =paste(round(coef_y[3,1],3),'[', round(ci_b[5,][[1]], 3),',',round(ci_b[5,][[2]],3),']',sep = "")
  result_mat$P[5] = round(p_values_y[[3]], 3)
  return(result_mat)
}
results_1 = humidity_phq8_step(phq_weather)
results_2 = humidity_phq8_step(phq_weather[phq_weather$p_id %in% summary_coeff$p_id[summary_coeff$humidity >= 0.2],])
results_3 = humidity_phq8_step(phq_weather[phq_weather$p_id %in% summary_coeff$p_id[summary_coeff$humidity <= -0.2],])
results_4 = humidity_phq8_step(phq_weather[phq_weather$p_id %in% summary_coeff$p_id[summary_coeff$humidity > -0.2 & summary_coeff$humidity < 0.2],])
results_combined = cbind(results_1, results_2[,2:3],results_3[,2:3],results_4[,2:3])
results_combined_all = rbind(results_combined_all, results_combined)

# X = Pressure, M=phq8, Y=step 
pressure_phq8_step <- function(data){
  # model m: M = aX 
  model_m <- lmer(sum_phq8 ~ pressure_mean + age+sex+site + marry_align+ edu_age_align + income_align+employment+ (1 | p_id), data = data)
  
  # model y: Y = aX + bM 
  model_y <- lmer(steps~ pressure_mean + sum_phq8 + age+sex+site + marry_align+ edu_age_align + income_align+employment+  (1 | p_id), data = data)
  
  # mediation analysis
  med.out <- mediate(model_m, model_y, 
                     treat = "pressure_mean", 
                     mediator = "sum_phq8",
                     sims = 500) 
  
  y <- summary(med.out)
  result_mat <- matrix(nrow = 5, ncol = 3)
  result_mat = as.data.frame(result_mat)
  colnames(result_mat) <- c("Path", "Est", "P")
  
  result_mat$Path[1] = "Total Effects"
  result_mat$Path[2] = "Direct Effects"
  result_mat$Path[3] = "Indirect Effects"
  result_mat$Path[4] = "a path"
  result_mat$Path[5] = "b path"
  
  
  result_mat$Est[1] =paste(round(y$tau.coef,3),'[', round(y$tau.ci[['2.5%']], 3),',',round(y$tau.ci[['97.5%']],3),']',sep = "")
  result_mat$P[1] = round(y$tau.p, 3)
  
  result_mat$Est[2] =paste(round(y$z0,3),'[', round(y$z0.ci[['2.5%']], 3),',',round(y$z0.ci[['97.5%']],3),']',sep = "")
  result_mat$P[2] = round(y$z0.p, 3)
  
  result_mat$Est[3] =paste(round(y$d0,3),'[', round(y$d0.ci[['2.5%']], 3),',',round(y$d0.ci[['97.5%']],3),']',sep = "")
  result_mat$P[3] = round(y$d0.p, 3)
  
  ci_a = confint(model_m)
  summary_m <- summary(model_m)
  coef_m <- summary_m$coefficients
  p_values_m <- 2 * (1 - pnorm(abs(coef_m[, "t value"])))
  
  result_mat$Est[4] =paste(round(coef_m[2,1],3),'[', round(ci_a[4,][[1]], 3),',',round(ci_a[4,][[2]],3),']',sep = "")
  result_mat$P[4] = round(p_values_m[[2]], 3)
  
  
  ci_b = confint(model_y)
  summary_y <- summary(model_y)
  coef_y <- summary_y$coefficients
  p_values_y <- 2 * (1 - pnorm(abs(coef_y[, "t value"])))
  
  result_mat$Est[5] =paste(round(coef_y[3,1],3),'[', round(ci_b[5,][[1]], 3),',',round(ci_b[5,][[2]],3),']',sep = "")
  result_mat$P[5] = round(p_values_y[[3]], 3)
  return(result_mat)
}
results_1 = pressure_phq8_step(phq_weather)
results_2 = pressure_phq8_step(phq_weather[phq_weather$p_id %in% summary_coeff$p_id[summary_coeff$pressure >= 0.2],])
results_3 = pressure_phq8_step(phq_weather[phq_weather$p_id %in% summary_coeff$p_id[summary_coeff$pressure <= -0.2],])
results_4 = pressure_phq8_step(phq_weather[phq_weather$p_id %in% summary_coeff$p_id[summary_coeff$pressure > -0.2 & summary_coeff$pressure < 0.2],])
results_combined = cbind(results_1, results_2[,2:3],results_3[,2:3],results_4[,2:3])
results_combined_all = rbind(results_combined_all, results_combined)

# X = Clouds, M=phq8, Y=step 
clouds_phq8_step <- function(data){
  # model m: M = aX 
  model_m <- lmer(sum_phq8 ~ clouds_mean + age+sex+site + marry_align+ edu_age_align + income_align+employment+ (1 | p_id), data = data)
  
  # model y: Y = aX + bM 
  model_y <- lmer(steps~ clouds_mean + sum_phq8 + age+sex+site + marry_align+ edu_age_align + income_align+employment+  (1 | p_id), data = data)
  
  # mediation analysis
  med.out <- mediate(model_m, model_y, 
                     treat = "clouds_mean", 
                     mediator = "sum_phq8",
                     sims = 500) 
  
  y <- summary(med.out)
  result_mat <- matrix(nrow = 5, ncol = 3)
  result_mat = as.data.frame(result_mat)
  colnames(result_mat) <- c("Path", "Est", "P")
  
  result_mat$Path[1] = "Total Effects"
  result_mat$Path[2] = "Direct Effects"
  result_mat$Path[3] = "Indirect Effects"
  result_mat$Path[4] = "a path"
  result_mat$Path[5] = "b path"
  
  
  result_mat$Est[1] =paste(round(y$tau.coef,3),'[', round(y$tau.ci[['2.5%']], 3),',',round(y$tau.ci[['97.5%']],3),']',sep = "")
  result_mat$P[1] = round(y$tau.p, 3)
  
  result_mat$Est[2] =paste(round(y$z0,3),'[', round(y$z0.ci[['2.5%']], 3),',',round(y$z0.ci[['97.5%']],3),']',sep = "")
  result_mat$P[2] = round(y$z0.p, 3)
  
  result_mat$Est[3] =paste(round(y$d0,3),'[', round(y$d0.ci[['2.5%']], 3),',',round(y$d0.ci[['97.5%']],3),']',sep = "")
  result_mat$P[3] = round(y$d0.p, 3)
  
  ci_a = confint(model_m)
  summary_m <- summary(model_m)
  coef_m <- summary_m$coefficients
  p_values_m <- 2 * (1 - pnorm(abs(coef_m[, "t value"])))
  
  result_mat$Est[4] =paste(round(coef_m[2,1],3),'[', round(ci_a[4,][[1]], 3),',',round(ci_a[4,][[2]],3),']',sep = "")
  result_mat$P[4] = round(p_values_m[[2]], 3)
  
  
  ci_b = confint(model_y)
  summary_y <- summary(model_y)
  coef_y <- summary_y$coefficients
  p_values_y <- 2 * (1 - pnorm(abs(coef_y[, "t value"])))
  
  result_mat$Est[5] =paste(round(coef_y[3,1],3),'[', round(ci_b[5,][[1]], 3),',',round(ci_b[5,][[2]],3),']',sep = "")
  result_mat$P[5] = round(p_values_y[[3]], 3)
  return(result_mat)
}
results_1 = clouds_phq8_step(phq_weather)
results_2 = clouds_phq8_step(phq_weather[phq_weather$p_id %in% summary_coeff$p_id[summary_coeff$clouds >= 0.2],])
results_3 = clouds_phq8_step(phq_weather[phq_weather$p_id %in% summary_coeff$p_id[summary_coeff$clouds <= -0.2],])
results_4 = clouds_phq8_step(phq_weather[phq_weather$p_id %in% summary_coeff$p_id[summary_coeff$clouds > -0.2 & summary_coeff$clouds < 0.2],])
results_combined = cbind(results_1, results_2[,2:3],results_3[,2:3],results_4[,2:3])
results_combined_all = rbind(results_combined_all, results_combined)

# X = Wind Speed, M=phq8, Y=step 
wind_phq8_step <- function(data){
  # model m: M = aX 
  model_m <- lmer(sum_phq8 ~ wind_speed_mean + age+sex+site + marry_align+ edu_age_align + income_align+employment+ (1 | p_id), data = data)
  
  # model y: Y = aX + bM 
  model_y <- lmer(steps~ wind_speed_mean + sum_phq8 + age+sex+site + marry_align+ edu_age_align + income_align+employment+  (1 | p_id), data = data)
  
  # mediation analysis
  med.out <- mediate(model_m, model_y, 
                     treat = "wind_speed_mean", 
                     mediator = "sum_phq8",
                     sims = 500) 
  
  y <- summary(med.out)
  result_mat <- matrix(nrow = 5, ncol = 3)
  result_mat = as.data.frame(result_mat)
  colnames(result_mat) <- c("Path", "Est", "P")
  
  result_mat$Path[1] = "Total Effects"
  result_mat$Path[2] = "Direct Effects"
  result_mat$Path[3] = "Indirect Effects"
  result_mat$Path[4] = "a path"
  result_mat$Path[5] = "b path"
  
  
  result_mat$Est[1] =paste(round(y$tau.coef,3),'[', round(y$tau.ci[['2.5%']], 3),',',round(y$tau.ci[['97.5%']],3),']',sep = "")
  result_mat$P[1] = round(y$tau.p, 3)
  
  result_mat$Est[2] =paste(round(y$z0,3),'[', round(y$z0.ci[['2.5%']], 3),',',round(y$z0.ci[['97.5%']],3),']',sep = "")
  result_mat$P[2] = round(y$z0.p, 3)
  
  result_mat$Est[3] =paste(round(y$d0,3),'[', round(y$d0.ci[['2.5%']], 3),',',round(y$d0.ci[['97.5%']],3),']',sep = "")
  result_mat$P[3] = round(y$d0.p, 3)
  
  ci_a = confint(model_m)
  summary_m <- summary(model_m)
  coef_m <- summary_m$coefficients
  p_values_m <- 2 * (1 - pnorm(abs(coef_m[, "t value"])))
  
  result_mat$Est[4] =paste(round(coef_m[2,1],3),'[', round(ci_a[4,][[1]], 3),',',round(ci_a[4,][[2]],3),']',sep = "")
  result_mat$P[4] = round(p_values_m[[2]], 3)
  
  
  ci_b = confint(model_y)
  summary_y <- summary(model_y)
  coef_y <- summary_y$coefficients
  p_values_y <- 2 * (1 - pnorm(abs(coef_y[, "t value"])))
  
  result_mat$Est[5] =paste(round(coef_y[3,1],3),'[', round(ci_b[5,][[1]], 3),',',round(ci_b[5,][[2]],3),']',sep = "")
  result_mat$P[5] = round(p_values_y[[3]], 3)
  return(result_mat)
}
results_1 = wind_phq8_step(phq_weather)
results_2 = wind_phq8_step(phq_weather[phq_weather$p_id %in% summary_coeff$p_id[summary_coeff$wind_speed >= 0.2],])
results_3 = wind_phq8_step(phq_weather[phq_weather$p_id %in% summary_coeff$p_id[summary_coeff$wind_speed <= -0.2],])
results_4 = wind_phq8_step(phq_weather[phq_weather$p_id %in% summary_coeff$p_id[summary_coeff$wind_speed > -0.2 & summary_coeff$wind_speed < 0.2],])
results_combined = cbind(results_1, results_2[,2:3],results_3[,2:3],results_4[,2:3])
results_combined_all = rbind(results_combined_all, results_combined)

# X = Sun_diff, M=phq8, Y=step 
sundiff_phq8_step <- function(data){
  # model m: M = aX 
  model_m <- lmer(sum_phq8 ~ sun_diff + age+sex+site + marry_align+ edu_age_align + income_align+employment+ (1 | p_id), data = data)
  
  # model y: Y = aX + bM 
  model_y <- lmer(steps~ sun_diff + sum_phq8 + age+sex+site + marry_align+ edu_age_align + income_align+employment+  (1 | p_id), data = data)
  
  # mediation analysis
  med.out <- mediate(model_m, model_y, 
                     treat = "sun_diff", 
                     mediator = "sum_phq8",
                     sims = 500) 
  
  y <- summary(med.out)
  result_mat <- matrix(nrow = 5, ncol = 3)
  result_mat = as.data.frame(result_mat)
  colnames(result_mat) <- c("Path", "Est", "P")
  
  result_mat$Path[1] = "Total Effects"
  result_mat$Path[2] = "Direct Effects"
  result_mat$Path[3] = "Indirect Effects"
  result_mat$Path[4] = "a path"
  result_mat$Path[5] = "b path"
  
  
  result_mat$Est[1] =paste(round(y$tau.coef,3),'[', round(y$tau.ci[['2.5%']], 3),',',round(y$tau.ci[['97.5%']],3),']',sep = "")
  result_mat$P[1] = round(y$tau.p, 3)
  
  result_mat$Est[2] =paste(round(y$z0,3),'[', round(y$z0.ci[['2.5%']], 3),',',round(y$z0.ci[['97.5%']],3),']',sep = "")
  result_mat$P[2] = round(y$z0.p, 3)
  
  result_mat$Est[3] =paste(round(y$d0,3),'[', round(y$d0.ci[['2.5%']], 3),',',round(y$d0.ci[['97.5%']],3),']',sep = "")
  result_mat$P[3] = round(y$d0.p, 3)
  
  ci_a = confint(model_m)
  summary_m <- summary(model_m)
  coef_m <- summary_m$coefficients
  p_values_m <- 2 * (1 - pnorm(abs(coef_m[, "t value"])))
  
  result_mat$Est[4] =paste(round(coef_m[2,1],3),'[', round(ci_a[4,][[1]], 3),',',round(ci_a[4,][[2]],3),']',sep = "")
  result_mat$P[4] = round(p_values_m[[2]], 3)
  
  
  ci_b = confint(model_y)
  summary_y <- summary(model_y)
  coef_y <- summary_y$coefficients
  p_values_y <- 2 * (1 - pnorm(abs(coef_y[, "t value"])))
  
  result_mat$Est[5] =paste(round(coef_y[3,1],3),'[', round(ci_b[5,][[1]], 3),',',round(ci_b[5,][[2]],3),']',sep = "")
  result_mat$P[5] = round(p_values_y[[3]], 3)
  return(result_mat)
}
results_1 = sundiff_phq8_step(phq_weather)
results_2 = sundiff_phq8_step(phq_weather[phq_weather$p_id %in% summary_coeff$p_id[summary_coeff$sun_diff>= 0.2],])
results_3 = sundiff_phq8_step(phq_weather[phq_weather$p_id %in% summary_coeff$p_id[summary_coeff$sun_diff <= -0.2],])
results_4 = sundiff_phq8_step(phq_weather[phq_weather$p_id %in% summary_coeff$p_id[summary_coeff$sun_diff > -0.2 & summary_coeff$sun_diff < 0.2],])
results_combined = cbind(results_1, results_2[,2:3],results_3[,2:3],results_4[,2:3])
results_combined_all = rbind(results_combined_all, results_combined)

write.csv(results_combined_all, 'mediator_is_phq8_500_new.csv',row.names = F)



## Step is Mediator
# divided steps by 1000
phq_weather$steps = phq_weather$steps/1000

# X = temp
temp_step_phq8 <- function(data){
  # model m: M = aX 
  model_m <- lmer(steps ~ temp_mean + age+sex+site + marry_align+ edu_age_align + income_align+employment+ (1 | p_id), data = data)
  
  # model y: Y = aX + bM 
  model_y <- lmer(sum_phq8~ temp_mean + steps + age+sex+site + marry_align+ edu_age_align + income_align+employment+  (1 | p_id), data = data)
  
  # mediation analysis
  med.out <- mediate(model_m, model_y, 
                     treat = "temp_mean", 
                     mediator = "steps",
                     sims = 500) 
  
  y <- summary(med.out)
  
  result_mat <- matrix(nrow = 5, ncol = 3)
  result_mat = as.data.frame(result_mat)
  colnames(result_mat) <- c("Path", "Est", "P")
  
  result_mat$Path[1] = "Total Effects"
  result_mat$Path[2] = "Direct Effects"
  result_mat$Path[3] = "Indirect Effects"
  result_mat$Path[4] = "a path"
  result_mat$Path[5] = "b path"
  
  
  result_mat$Est[1] =paste(round(y$tau.coef,7),'[', round(y$tau.ci[['2.5%']], 3),',',round(y$tau.ci[['97.5%']],3),']',sep = "")
  
  result_mat$P[1] = round(y$tau.p, 3)
  
  result_mat$Est[2] =paste(round(y$z0,3),'[', round(y$z0.ci[['2.5%']], 3),',',round(y$z0.ci[['97.5%']],3),']',sep = "")
  result_mat$P[2] = round(y$z0.p, 3)
  
  result_mat$Est[3] =paste(round(y$d0,3),'[', round(y$d0.ci[['2.5%']], 3),',',round(y$d0.ci[['97.5%']],3),']',sep = "")
  result_mat$P[3] = round(y$d0.p, 3)
  
  ci_a = confint(model_m)
  summary_m <- summary(model_m)
  coef_m <- summary_m$coefficients
  p_values_m <- 2 * (1 - pnorm(abs(coef_m[, "t value"])))
  
  result_mat$Est[4] =paste(round(coef_m[2,1],3),'[', round(ci_a[4,][[1]], 3),',',round(ci_a[4,][[2]],3),']',sep = "")
  result_mat$P[4] = round(p_values_m[[2]], 3)
  
  
  ci_b = confint(model_y)
  summary_y <- summary(model_y)
  coef_y <- summary_y$coefficients
  p_values_y <- 2 * (1 - pnorm(abs(coef_y[, "t value"])))
  
  result_mat$Est[5] =paste(round(coef_y[3,1],3),'[', round(ci_b[5,][[1]], 3),',',round(ci_b[5,][[2]],3),']',sep = "")
  result_mat$P[5] = round(p_values_y[[3]], 3)
  return(result_mat)
}
results_1 = temp_step_phq8(phq_weather)
results_2 = temp_step_phq8(phq_weather[phq_weather$p_id %in% summary_coeff$p_id[summary_coeff$temp >= 0.2],])
results_3 = temp_step_phq8(phq_weather[phq_weather$p_id %in% summary_coeff$p_id[summary_coeff$temp <= -0.2],])
results_4 = temp_step_phq8(phq_weather[phq_weather$p_id %in% summary_coeff$p_id[summary_coeff$temp > -0.2 & summary_coeff$temp < 0.2],])
results_combined_all = cbind(results_1, results_2[,2:3],results_3[,2:3],results_4[,2:3])


# X = humidity
humidity_step_phq8 <- function(data){
  # model m: M = aX 
  model_m <- lmer(steps ~ humidity_mean + age+sex+site + marry_align+ edu_age_align + income_align+employment+ (1 | p_id), data = data)
  
  # model y: Y = aX + bM 
  model_y <- lmer(sum_phq8~ humidity_mean + steps + age+sex+site + marry_align+ edu_age_align + income_align+employment+  (1 | p_id), data = data)
  
  # mediation analysis
  med.out <- mediate(model_m, model_y, 
                     treat = "humidity_mean", 
                     mediator = "steps",
                     sims = 500) 
  
  y <- summary(med.out)
  
  result_mat <- matrix(nrow = 5, ncol = 3)
  result_mat = as.data.frame(result_mat)
  colnames(result_mat) <- c("Path", "Est", "P")
  
  result_mat$Path[1] = "Total Effects"
  result_mat$Path[2] = "Direct Effects"
  result_mat$Path[3] = "Indirect Effects"
  result_mat$Path[4] = "a path"
  result_mat$Path[5] = "b path"
  
  
  result_mat$Est[1] =paste(round(y$tau.coef,3),'[', round(y$tau.ci[['2.5%']], 3),',',round(y$tau.ci[['97.5%']],3),']',sep = "")
  result_mat$P[1] = round(y$tau.p, 3)
  
  result_mat$Est[2] =paste(round(y$z0,3),'[', round(y$z0.ci[['2.5%']], 3),',',round(y$z0.ci[['97.5%']],3),']',sep = "")
  result_mat$P[2] = round(y$z0.p, 3)
  
  result_mat$Est[3] =paste(round(y$d0,4),'[', round(y$d0.ci[['2.5%']], 3),',',round(y$d0.ci[['97.5%']],3),']',sep = "")
  result_mat$P[3] = round(y$d0.p, 3)
  
  ci_a = confint(model_m)
  summary_m <- summary(model_m)
  coef_m <- summary_m$coefficients
  p_values_m <- 2 * (1 - pnorm(abs(coef_m[, "t value"])))
  
  result_mat$Est[4] =paste(round(coef_m[2,1],3),'[', round(ci_a[4,][[1]], 3),',',round(ci_a[4,][[2]],3),']',sep = "")
  result_mat$P[4] = round(p_values_m[[2]], 3)
  
  
  ci_b = confint(model_y)
  summary_y <- summary(model_y)
  coef_y <- summary_y$coefficients
  p_values_y <- 2 * (1 - pnorm(abs(coef_y[, "t value"])))
  
  result_mat$Est[5] =paste(round(coef_y[3,1],3),'[', round(ci_b[5,][[1]], 3),',',round(ci_b[5,][[2]],3),']',sep = "")
  result_mat$P[5] = round(p_values_y[[3]], 3)
  return(result_mat)
}
results_1 = humidity_step_phq8(phq_weather)
results_2 = humidity_step_phq8(phq_weather[phq_weather$p_id %in% summary_coeff$p_id[summary_coeff$humidity >= 0.2],])
results_3 = humidity_step_phq8(phq_weather[phq_weather$p_id %in% summary_coeff$p_id[summary_coeff$humidity <= -0.2],])
results_4 = humidity_step_phq8(phq_weather[phq_weather$p_id %in% summary_coeff$p_id[summary_coeff$humidity > -0.2 & summary_coeff$humidity < 0.2],])
results_combined = cbind(results_1, results_2[,2:3],results_3[,2:3],results_4[,2:3])
results_combined_all = rbind(results_combined_all, results_combined)

# X = Pressure
pressure_step_phq8 <- function(data){
  # model m: M = aX 
  model_m <- lmer(steps ~ pressure_mean + age+sex+site + marry_align+ edu_age_align + income_align+employment+ (1 | p_id), data = data)
  
  # model y: Y = aX + bM 
  model_y <- lmer(sum_phq8~ pressure_mean + steps + age+sex+site + marry_align+ edu_age_align + income_align+employment+  (1 | p_id), data = data)
  
  # mediation analysis
  med.out <- mediate(model_m, model_y, 
                     treat = "pressure_mean", 
                     mediator = "steps",
                     sims = 500) 
  
  y <- summary(med.out)
  
  result_mat <- matrix(nrow = 5, ncol = 3)
  result_mat = as.data.frame(result_mat)
  colnames(result_mat) <- c("Path", "Est", "P")
  
  result_mat$Path[1] = "Total Effects"
  result_mat$Path[2] = "Direct Effects"
  result_mat$Path[3] = "Indirect Effects"
  result_mat$Path[4] = "a path"
  result_mat$Path[5] = "b path"
  
  
  result_mat$Est[1] =paste(round(y$tau.coef,3),'[', round(y$tau.ci[['2.5%']], 3),',',round(y$tau.ci[['97.5%']],3),']',sep = "")
  result_mat$P[1] = round(y$tau.p, 3)
  
  result_mat$Est[2] =paste(round(y$z0,3),'[', round(y$z0.ci[['2.5%']], 3),',',round(y$z0.ci[['97.5%']],3),']',sep = "")
  result_mat$P[2] = round(y$z0.p, 3)
  
  result_mat$Est[3] =paste(round(y$d0,4),'[', round(y$d0.ci[['2.5%']], 3),',',round(y$d0.ci[['97.5%']],3),']',sep = "")
  result_mat$P[3] = round(y$d0.p, 3)
  
  ci_a = confint(model_m)
  summary_m <- summary(model_m)
  coef_m <- summary_m$coefficients
  p_values_m <- 2 * (1 - pnorm(abs(coef_m[, "t value"])))
  
  result_mat$Est[4] =paste(round(coef_m[2,1],3),'[', round(ci_a[4,][[1]], 3),',',round(ci_a[4,][[2]],3),']',sep = "")
  result_mat$P[4] = round(p_values_m[[2]], 3)
  
  
  ci_b = confint(model_y)
  summary_y <- summary(model_y)
  coef_y <- summary_y$coefficients
  p_values_y <- 2 * (1 - pnorm(abs(coef_y[, "t value"])))
  
  result_mat$Est[5] =paste(round(coef_y[3,1],3),'[', round(ci_b[5,][[1]], 3),',',round(ci_b[5,][[2]],3),']',sep = "")
  result_mat$P[5] = round(p_values_y[[3]], 3)
  return(result_mat)
}
results_1 = pressure_step_phq8(phq_weather)
results_2 = pressure_step_phq8(phq_weather[phq_weather$p_id %in% summary_coeff$p_id[summary_coeff$pressure >= 0.2],])
results_3 = pressure_step_phq8(phq_weather[phq_weather$p_id %in% summary_coeff$p_id[summary_coeff$pressure <= -0.2],])
results_4 = pressure_step_phq8(phq_weather[phq_weather$p_id %in% summary_coeff$p_id[summary_coeff$pressure > -0.2 & summary_coeff$pressure < 0.2],])
results_combined = cbind(results_1, results_2[,2:3],results_3[,2:3],results_4[,2:3])
results_combined_all = rbind(results_combined_all, results_combined)

# X = Clouds
clouds_step_phq8 <- function(data){
  # model m: M = aX 
  model_m <- lmer(steps ~ clouds_mean + age+sex+site + marry_align+ edu_age_align + income_align+employment+ (1 | p_id), data = data)
  
  # model y: Y = aX + bM 
  model_y <- lmer(sum_phq8~ clouds_mean + steps + age+sex+site + marry_align+ edu_age_align + income_align+employment+  (1 | p_id), data = data)
  
  # mediation analysis
  med.out <- mediate(model_m, model_y, 
                     treat = "clouds_mean", 
                     mediator = "steps",
                     sims = 500) 
  
  y <- summary(med.out)
  result_mat <- matrix(nrow = 5, ncol = 3)
  result_mat = as.data.frame(result_mat)
  colnames(result_mat) <- c("Path", "Est", "P")
  
  result_mat$Path[1] = "Total Effects"
  result_mat$Path[2] = "Direct Effects"
  result_mat$Path[3] = "Indirect Effects"
  result_mat$Path[4] = "a path"
  result_mat$Path[5] = "b path"
  
  
  result_mat$Est[1] =paste(round(y$tau.coef,3),'[', round(y$tau.ci[['2.5%']], 3),',',round(y$tau.ci[['97.5%']],3),']',sep = "")
  result_mat$P[1] = round(y$tau.p, 3)
  
  result_mat$Est[2] =paste(round(y$z0,3),'[', round(y$z0.ci[['2.5%']], 3),',',round(y$z0.ci[['97.5%']],3),']',sep = "")
  result_mat$P[2] = round(y$z0.p, 3)
  
  result_mat$Est[3] =paste(round(y$d0,4),'[', round(y$d0.ci[['2.5%']], 3),',',round(y$d0.ci[['97.5%']],3),']',sep = "")
  result_mat$P[3] = round(y$d0.p, 3)
  
  ci_a = confint(model_m)
  summary_m <- summary(model_m)
  coef_m <- summary_m$coefficients
  p_values_m <- 2 * (1 - pnorm(abs(coef_m[, "t value"])))
  
  result_mat$Est[4] =paste(round(coef_m[2,1],3),'[', round(ci_a[4,][[1]], 3),',',round(ci_a[4,][[2]],3),']',sep = "")
  result_mat$P[4] = round(p_values_m[[2]], 3)
  
  
  ci_b = confint(model_y)
  summary_y <- summary(model_y)
  coef_y <- summary_y$coefficients
  p_values_y <- 2 * (1 - pnorm(abs(coef_y[, "t value"])))
  
  result_mat$Est[5] =paste(round(coef_y[3,1],3),'[', round(ci_b[5,][[1]], 3),',',round(ci_b[5,][[2]],3),']',sep = "")
  result_mat$P[5] = round(p_values_y[[3]], 3)
  return(result_mat)
}
results_1 = clouds_step_phq8(phq_weather)
results_2 = clouds_step_phq8(phq_weather[phq_weather$p_id %in% summary_coeff$p_id[summary_coeff$clouds >= 0.2],])
results_3 = clouds_step_phq8(phq_weather[phq_weather$p_id %in% summary_coeff$p_id[summary_coeff$clouds <= -0.2],])
results_4 = clouds_step_phq8(phq_weather[phq_weather$p_id %in% summary_coeff$p_id[summary_coeff$clouds > -0.2 & summary_coeff$clouds < 0.2],])
results_combined = cbind(results_1, results_2[,2:3],results_3[,2:3],results_4[,2:3])
results_combined_all = rbind(results_combined_all, results_combined)

# X = Wind Speed
wind_step_phq8 <- function(data){
  # model m: M = aX 
  model_m <- lmer(steps ~ wind_speed_mean + age+sex+site + marry_align+ edu_age_align + income_align+employment+ (1 | p_id), data = data)
  
  # model y: Y = aX + bM 
  model_y <- lmer(sum_phq8~ wind_speed_mean + steps + age+sex+site + marry_align+ edu_age_align + income_align+employment+  (1 | p_id), data = data)
  
  # mediation analysis
  med.out <- mediate(model_m, model_y, 
                     treat = "wind_speed_mean", 
                     mediator = "steps",
                     sims = 500) 
  
  y <- summary(med.out)
  
  result_mat <- matrix(nrow = 5, ncol = 3)
  result_mat = as.data.frame(result_mat)
  colnames(result_mat) <- c("Path", "Est", "P")
  
  result_mat$Path[1] = "Total Effects"
  result_mat$Path[2] = "Direct Effects"
  result_mat$Path[3] = "Indirect Effects"
  result_mat$Path[4] = "a path"
  result_mat$Path[5] = "b path"
  
  
  result_mat$Est[1] =paste(round(y$tau.coef,3),'[', round(y$tau.ci[['2.5%']], 3),',',round(y$tau.ci[['97.5%']],3),']',sep = "")
  result_mat$P[1] = round(y$tau.p, 3)
  
  result_mat$Est[2] =paste(round(y$z0,3),'[', round(y$z0.ci[['2.5%']], 3),',',round(y$z0.ci[['97.5%']],3),']',sep = "")
  result_mat$P[2] = round(y$z0.p, 3)
  
  result_mat$Est[3] =paste(round(y$d0,4),'[', round(y$d0.ci[['2.5%']], 3),',',round(y$d0.ci[['97.5%']],3),']',sep = "")
  result_mat$P[3] = round(y$d0.p, 3)
  
  ci_a = confint(model_m)
  summary_m <- summary(model_m)
  coef_m <- summary_m$coefficients
  p_values_m <- 2 * (1 - pnorm(abs(coef_m[, "t value"])))
  
  result_mat$Est[4] =paste(round(coef_m[2,1],3),'[', round(ci_a[4,][[1]], 3),',',round(ci_a[4,][[2]],3),']',sep = "")
  result_mat$P[4] = round(p_values_m[[2]], 3)
  
  
  ci_b = confint(model_y)
  summary_y <- summary(model_y)
  coef_y <- summary_y$coefficients
  p_values_y <- 2 * (1 - pnorm(abs(coef_y[, "t value"])))
  
  result_mat$Est[5] =paste(round(coef_y[3,1],3),'[', round(ci_b[5,][[1]], 3),',',round(ci_b[5,][[2]],3),']',sep = "")
  result_mat$P[5] = round(p_values_y[[3]], 3)
  return(result_mat)
}
results_1 = wind_step_phq8(phq_weather)
results_2 = wind_step_phq8(phq_weather[phq_weather$p_id %in% summary_coeff$p_id[summary_coeff$wind_speed >= 0.2],])
results_3 = wind_step_phq8(phq_weather[phq_weather$p_id %in% summary_coeff$p_id[summary_coeff$wind_speed <= -0.2],])
results_4 = wind_step_phq8(phq_weather[phq_weather$p_id %in% summary_coeff$p_id[summary_coeff$wind_speed > -0.2 & summary_coeff$wind_speed < 0.2],])
results_combined = cbind(results_1, results_2[,2:3],results_3[,2:3],results_4[,2:3])
results_combined_all = rbind(results_combined_all, results_combined)

# X = Sun_diff
sundiff_step_phq8 <- function(data){
  # model m: M = aX 
  model_m <- lmer(steps ~ sun_diff + age+sex+site + marry_align+ edu_age_align + income_align+employment+ (1 | p_id), data = data)
  
  # model y: Y = aX + bM 
  model_y <- lmer(sum_phq8~ sun_diff + steps + age+sex+site + marry_align+ edu_age_align + income_align+employment+  (1 | p_id), data = data)
  
  # mediation analysis
  med.out <- mediate(model_m, model_y, 
                     treat = "sun_diff", 
                     mediator = "steps",
                     sims = 500) 
  
  y <- summary(med.out)
  
  result_mat <- matrix(nrow = 5, ncol = 3)
  result_mat = as.data.frame(result_mat)
  colnames(result_mat) <- c("Path", "Est", "P")
  
  result_mat$Path[1] = "Total Effects"
  result_mat$Path[2] = "Direct Effects"
  result_mat$Path[3] = "Indirect Effects"
  result_mat$Path[4] = "a path"
  result_mat$Path[5] = "b path"
  
  
  result_mat$Est[1] =paste(round(y$tau.coef,4),'[', round(y$tau.ci[['2.5%']], 3),',',round(y$tau.ci[['97.5%']],3),']',sep = "")
  result_mat$P[1] = round(y$tau.p, 3)
  
  result_mat$Est[2] =paste(round(y$z0,7),'[', round(y$z0.ci[['2.5%']], 3),',',round(y$z0.ci[['97.5%']],3),']',sep = "")
  result_mat$P[2] = round(y$z0.p, 3)
  
  result_mat$Est[3] =paste(round(y$d0,4),'[', round(y$d0.ci[['2.5%']], 3),',',round(y$d0.ci[['97.5%']],3),']',sep = "")
  result_mat$P[3] = round(y$d0.p, 3)
  
  ci_a = confint(model_m)
  summary_m <- summary(model_m)
  coef_m <- summary_m$coefficients
  p_values_m <- 2 * (1 - pnorm(abs(coef_m[, "t value"])))
  
  result_mat$Est[4] =paste(round(coef_m[2,1],4),'[', round(ci_a[4,][[1]], 3),',',round(ci_a[4,][[2]],3),']',sep = "")
  result_mat$P[4] = round(p_values_m[[2]], 4)
  
  
  ci_b = confint(model_y)
  summary_y <- summary(model_y)
  coef_y <- summary_y$coefficients
  p_values_y <- 2 * (1 - pnorm(abs(coef_y[, "t value"])))
  
  result_mat$Est[5] =paste(round(coef_y[3,1],4),'[', round(ci_b[5,][[1]], 3),',',round(ci_b[5,][[2]],3),']',sep = "")
  result_mat$P[5] = round(p_values_y[[3]], 4)
  return(result_mat)
}
results_1 = sundiff_step_phq8(phq_weather)
results_2 = sundiff_step_phq8(phq_weather[phq_weather$p_id %in% summary_coeff$p_id[summary_coeff$sun_diff>= 0.2],])
results_3 = sundiff_step_phq8(phq_weather[phq_weather$p_id %in% summary_coeff$p_id[summary_coeff$sun_diff <= -0.2],])
results_4 = sundiff_step_phq8(phq_weather[phq_weather$p_id %in% summary_coeff$p_id[summary_coeff$sun_diff > -0.2 & summary_coeff$sun_diff < 0.2],])
results_combined = cbind(results_1, results_2[,2:3],results_3[,2:3],results_4[,2:3])
results_combined_all = rbind(results_combined_all, results_combined)
write.csv(results_combined_all, 'mediator_is_steps_500_new.csv',row.names = F)





