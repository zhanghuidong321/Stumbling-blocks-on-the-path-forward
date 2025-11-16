
rm(list = ls()) #
setwd("D:/R工作路径/CHARLS date")  

###3. library R package ##
library(haven)        # read Stata (.dta) 
library(tidyverse)    #  dplyr, ggplot2  
library(tableone)     # baseline
library(segmented)    # regression analysis
library(splines)      # Spline Regression
library(Hmisc)        # Statistical functions
library(rms)          #  
library(forestplot) # frost
library(stringr)      
library(dplyr)
library(survival)    # COX回归核心包
library(survminer)   # 结果可视化
library(stargazer)   # 表格输出
library(pROC)     
library(ggplot2) 
library(grid)
library(car) 
library(knitr)

# #### 3.import data####

##3.1 data of 2011
demo2011 <- read_dta("2011/demographic_background.dta")
health2011 <- read_dta("2011/health_status_and_functioning.dta")
bio2011 <- read_dta("2011/biomarker.dta")
blood2011 <- read_dta("2011/Blood.dta")

merged2011<-demo2011%>% left_join(health2011,by="ID")%>%
  left_join(bio2011,by="ID")%>%
  left_join(blood2011,by="ID")%>%
  mutate(ID = paste0(householdID, '0', str_sub(ID, -2, -1)))%>% 
  filter(da007_13_ == 2 | is.na(da007_13_))#the first step:filter non_arthritis  


##3.2 data of 2020

health2020<- read_dta('2020/health_status_and_functioning.dta')

arthritis2020 <- health2020 %>% dplyr::select(ID,da003_3_, da002_3_, da003_6_, da002_6_,da002_14_, da003_14_) %>% 
  mutate(
    Arthritis_Status20 = case_when(
      da003_14_ == 1 ~ "arthritis",
      da003_14_ == 2 ~ "non_arthritis",
      is.na(da003_14_) & da002_14_ %in% c(1, 2, 3) ~ "arthritis",
      is.na(da003_14_) & da002_14_ == 99 ~ "non_arthritis",
      TRUE ~ NA_character_
    ),
    liver_dis20 = case_when(
      da003_6_ == 1 ~ 1,
      da003_6_ == 2 ~ 0,
      is.na(da003_6_) & da002_6_ %in% c(1, 2, 3) ~ 1,
      is.na(da003_6_) & da002_6_ == 99 ~ 0,
      TRUE ~ 0 
    ),
    diabetes20 = case_when(
      da003_3_ == 1 ~ 1,
      da003_3_ == 2 ~ 0,
      is.na(da003_3_) & da002_3_ %in% c(1, 2, 3) ~ 1,
      is.na(da003_3_) & da002_3_ == 99 ~ 0,
    ),
    co_ar_dia = case_when(
      Arthritis_Status20 %in% "arthritis" &diabetes20 == 1 ~ "co_ar_dia",
      Arthritis_Status20 %in% "non_arthritis" |diabetes20 == 0 ~ "non_co_ar_dia",
      TRUE ~ NA_character_
    )
  )%>%dplyr::select(-c(da003_3_, da002_3_, da003_6_, da002_6_,da002_14_, da003_14_))%>%
  filter(!is.na(Arthritis_Status20))


status_distribution <- table(arthritis2020$Arthritis_Status20)
print(status_distribution)

#### 4 merged data####

#table(baseline$rgender)

baseline<-merged2011 %>% left_join(arthritis2020 , by= "ID")%>%
  transmute(
    ID,
    age = ifelse(!is.na(ba002_1), 2011 - ba002_1, ba004),
    age_group = ifelse(age < 60, "<60", ">=60"),
    gender = ifelse(rgender == 1, 'Male', 'Female'),
    education = ifelse(bd001 %in% 1:4, 'Primary education',
                       ifelse(bd001 %in% 5:7, "Secondary Education",'Higher Education')),
    marital = ifelse(be001 %in% 1:2, 'married', 'unmarried'),
    location = ifelse(bb006 == 1, 'village', 'city'),
    smoking = ifelse(da059 == 2, 'nonsmoker', 'smoker'),
    drinking = ifelse(da067 == 3, 'nondrinker',"drinker"),
    sleep_time = da049,
    sbp1 = ifelse(qa003 > 300, NA, qa003),
    sbp2 = ifelse(qa007 > 300, NA, qa007),
    sbp3 = ifelse(qa011 > 300, NA, qa011),
    dbp1 = ifelse(qa004 > 200, NA, qa004),
    dbp2 = ifelse(qa008 > 200, NA, qa008),
    dbp3 = ifelse(qa012 > 200, NA, qa012), 
    sbp = (sbp1 + sbp2 + sbp3)/3,
    dbp = (dbp1 + dbp2 + dbp3)/3,
    height = ifelse(qi002 > 300|qi002 < 50, NA, qi002),
    weight = ifelse(ql002 < 30, NA, ql002),
    BMI = weight / ((height / 100)^2),
    WC = ifelse(qm002 > 300, NA, qm002),
    FBG = newglu,
    TG = newtg * 0.011,
    hdl = newhdl * 0.026,
    TC = newcho * 0.011,
    ldl = newldl* 0.011,
    
    VAI = ifelse(gender == 'Male',
                 WC / (39.68 + 1.88 * BMI) * (TG / 1.03) * (1.31 / hdl),
                 WC / (36.58 + 1.89 * BMI) * (TG / 0.81) * (1.52 / hdl)),
    #health status
    Arthritis_Status20,
    diabetes20,
    co_ar_dia ,
    co_disease = ifelse(co_ar_dia%in% 'co_ar_dia',1,0),
    health11 = case_when(
      da001 %in% c(1,2,3)|da002 %in% c(6,7,8)~'Good',
      da001==4|da002==9 ~'Fair',
      da001==5|da002 == 10 ~'Poor'
    ),
    heart_diseases11 = case_when(
      da007_7_ == 1 ~ "heart_diseases",
      da007_7_ == 2 ~ "non_heart_diseases",
      TRUE ~ NA_character_
    ),
    dyslipidemia11 = case_when(
      da007_2_ == 1 ~ "abnormal",
      da007_2_ == 2 ~ "normal",
      TC >= 5.2 ~ "abnormal",  # TC≥5.2 mmol/L
      TG >= 1.7 ~ "abnormal",  # TG≥1.7 mmol/L
      ldl >= 3.4 ~ "abnormal", # LDL-C≥3.4 mmol/L
      gender == "male" & hdl < 1.0 ~ "abnormal",
      gender == "female" & hdl < 1.3 ~ "abnormal",
      TRUE ~ "normal"
    ) 
  ) %>%
  filter(age >= 45,
         #!is.na(education), 
         #!is.na(marital),
         #!is.na(smoking),
         #!is.na(drinking),
         !is.na(VAI), 
         !is.na(Arthritis_Status20))%>%
  mutate(
    VAI_Q = case_when(
      VAI <= quantile(VAI, 0.25, na.rm = TRUE) ~ "Q1",
      VAI <= quantile(VAI, 0.5, na.rm = TRUE) ~ "Q2",
      VAI <= quantile(VAI, 0.75, na.rm = TRUE) ~ "Q3",
      TRUE ~ "Q4"),
    VAI_per_IQR = VAI / IQR(VAI, na.rm = TRUE)
  ) %>%dplyr::select(-c(sbp1:dbp3, height, weight))

write.csv(baseline, file = "2011-2020baseline_dataset.csv", row.names = FALSE, na = "")

#### 5.Descriptive statistics####

table1_arthritis <- CreateTableOne(
  vars = colnames(baseline)[2:29],
  strata = "co_ar_dia",
  data = baseline,
  addOverall = TRUE,
  argsNormal = list(var.equal = TRUE)  
)

# 打印基线资料表
table_final_arthritis <- print(
  table1_arthritis,
  showAllLevels = TRUE,
  exact = NULL,  
  pDigits = 3,   
  catDigits = 1, 
  contDigits = 1,
  smd = FALSE,    # 不计算标准化均数差
  nonnormal = NULL, 
  printToggle = FALSE, 
  quote = FALSE, 
  noSpaces = TRUE 
)

# view table 
print(table_final_arthritis)
View(table_final_arthritis)


####6、Logistic regression ####

# 1.  model 1 Unadjusted

fit1 <- glm(co_disease ~ VAI_per_IQR, data = baseline, family = binomial())
summary(fit1)
# Calculate the OR value and its 95% confidence interval
OR1 <- exp(coef(fit1))
CI1 <- exp(confint(fit1))
cat("Unadjusted OR：", OR1, "\n")
cat("Unadjusted 95%CI：", CI1, "\n")

# Regression Analysis by VAI Percentile

baseline <- baseline %>% 
  group_by(VAI_Q) %>% 
  mutate(VAI_trend = median(VAI, na.rm = TRUE)) %>%
  ungroup()

# Percentile model

fit2 <- glm(co_disease ~ VAI_Q, data = baseline, family = binomial())
summary(fit2)
OR2 <- exp(coef(fit2))
CI2 <- exp(confint(fit2))
cat("Percentile model OR ：", OR2, "\n")
cat("Percentile model 95% CI：", CI2, "\n")

# Trend Testing

fit3 <- glm(co_disease ~ VAI_trend, data = baseline, family = binomial())
summary(fit3)
OR3 <- exp(coef(fit3))
CI3 <- exp(confint(fit3))
cat("Trend  OR ：", OR3, "\n")
cat("Trend  95% CI：", CI3, "\n")

# 2.  model 2 Adjusted Logistic Regression Model

fit4 <- glm(co_disease ~ VAI_per_IQR + age + gender  + location + marital, 
            data = baseline, family = binomial())
summary(fit4)
OR4 <- exp(coef(fit4))
CI4 <- exp(confint(fit4))
cat("Adjusted OR：", OR4, "\n")
cat("Adjusted 95% CI：", CI4, "\n")

# Adjusted Percentile model

fit5 <- glm(co_disease ~ VAI_Q + age + gender + location + marital, 
            data = baseline, family = binomial())
summary(fit5)
OR5 <- exp(coef(fit5))
CI5 <- exp(confint(fit5))
cat("Adjusted Percentile OR ：", OR5, "\n")
cat("Adjusted Percentile 95% CI：", CI5, "\n")

# Adjusted Trend Testing + education

fit6 <- glm(co_disease ~ VAI_trend + age + gender  + location + marital, 
            data = baseline, family = binomial())
summary(fit6)
OR6 <- exp(coef(fit6))
CI6 <- exp(confint(fit6))
cat("Adjusted Trend Testing OR ：", OR6, "\n")
cat("Adjusted Trend Testing 95% CI：", CI6, "\n")

#3.model 3 full adjusted model

fit7 <- glm(co_disease ~ VAI_per_IQR + age + gender  + location + marital+
              smoking + drinking + sbp + dbp + BMI + health11,
            data = baseline, family = binomial())
summary(fit7)
OR7 <- exp(coef(fit7))
CI7 <- exp(confint(fit7))
cat("full adjusted OR ：", OR7, "\n")
cat("full adjusted 95% CI：", CI7, "\n")

# Full Adjusted Percentile model
fit8 <- glm(co_disease ~ VAI_Q + age + gender  + location + marital+
              smoking + drinking + sbp + dbp + BMI + health11,
            data = baseline, family = binomial())
summary(fit8)
OR8 <- exp(coef(fit8))
CI8 <- exp(confint(fit8))
cat("full adjusted percentile OR ：", OR8, "\n")
cat("full adjusted percentile 95% CI：", CI8, "\n")

# full adjusted trend
fit9 <- glm(co_disease ~ VAI_trend + age + gender  + location + marital+
              smoking + drinking + sbp + dbp + BMI + health11,
            data = baseline, family = binomial())
summary(fit9)
OR9 <- exp(coef(fit9))
CI9 <- exp(confint(fit9))
cat("full adjusted trend OR ：", OR9, "\n")
cat("full adjusted trend 95% CI：", CI9, "\n")

####（1）Sensitivity Analysis (Assessing Model Robustness)####

#----6.1.1Method 1: Remove Outliers / High Leverage Points----

baseline_fit9 <- na.omit(baseline[, all.vars(formula(fit9))])
valid_rows <- rownames(baseline_fit9)
n_valid <- nrow(baseline_fit9)

# Calculate the influence statistic (for screening out outliers)

influence_stats <- data.frame(
  hat = hatvalues(fit9),
  cook = cooks.distance(fit9),
  original_row_num = as.integer(valid_rows)
)

# Define the criteria for outliers

k <- length(coef(fit9)) - 1
high_hat <- influence_stats$hat > 2*(k+1)/n_valid
high_cook <- influence_stats$cook > 4/n_valid
outlier_original_rows <- influence_stats$original_row_num[high_hat | high_cook]

cat("Effective sample size for model fit9：", n_valid, "\n")
cat("Detected", length(outlier_original_rows), "outlier/high leverage point\n")

# Data after removing outliers

baseline_no_outlier <- baseline %>%
  filter(row_number() %in% setdiff(1:nrow(baseline), outlier_original_rows))

# Fit the formula using fit7 to fit_sens1
fit_sens1 <- glm(
  formula(fit7),
  data = baseline_no_outlier,
  family = binomial()
)

cat("Coefficient name for fit_sens1：", names(coef(fit_sens1)), "\n")

# Extraction VAI_per_IQR
sens1_OR <- exp(coef(fit_sens1))["VAI_per_IQR"]
sens1_CI <- exp(confint(fit_sens1))["VAI_per_IQR", ]
sens1_p <- summary(fit_sens1)$coefficients["VAI_per_IQR", "Pr(>|z|)"]

# Output comparison results

cat("\n===== Sensitivity Analysis 1: Comparison of Results After Removing Outliers =====\n")
cat("Original Full Model（fit7）OR（VAI per IQR）：", round(exp(coef(fit7)["VAI_per_IQR"]), 3), "\n")
cat("Original Full Model95%CI：", round(exp(confint(fit7))["VAI_per_IQR", ], 3), "，p-value：", round(summary(fit7)$coefficients["VAI_per_IQR", "Pr(>|z|)"], 3), "\n")
cat("After removing outliersOR：", round(sens1_OR, 3), "\n")
cat("After removing outliers95%CI：", round(sens1_CI, 3), "，p-value：", round(sens1_p, 3), "\n")

#----6.1.2 Modify the Covariate Adjustment Set----

# Plan 1：Remove the covariate “health11” and refit the full model
fit_sens2 <- glm(
  co_disease ~ VAI_per_IQR + age + gender + location + marital + smoking + drinking + sbp + dbp + BMI,
  data = baseline, family = binomial()
)

# Extraction Results

sens2_OR <- exp(coef(fit_sens2))["VAI_per_IQR"]
sens2_CI <- exp(confint(fit_sens2))["VAI_per_IQR", ]
sens2_p <- summary(fit_sens2)$coefficients["VAI_per_IQR", "Pr(>|z|)"]

# Output Comparison

cat("Original Full Model OR：", round(exp(coef(fit7)["VAI_per_IQR"]), 3), "\n")
cat("removing “health11” OR：", round(sens2_OR, 3), "，95%CI：", round(sens2_CI, 3), "，p-value：", round(sens2_p, 3), "\n")

# Plan 2：Remove the covariate “drinking” and refit the full model

fit_sens2 <- glm(
  co_disease ~ VAI_per_IQR + age + gender + location + marital + smoking + health11 + sbp + dbp + BMI,
  data = baseline, family = binomial()
)

# Extraction Results
sens2_OR <- exp(coef(fit_sens2))["VAI_per_IQR"]
sens2_CI <- exp(confint(fit_sens2))["VAI_per_IQR", ]
sens2_p <- summary(fit_sens2)$coefficients["VAI_per_IQR", "Pr(>|z|)"]

# Output Comparison
cat("Original Full Model OR：", round(exp(coef(fit7)["VAI_per_IQR"]), 3), "\n")
cat("removing “drinking” OR：", round(sens2_OR, 3), "，95%CI：", round(sens2_CI, 3), "，p-value：", round(sens2_p, 3), "\n")

#----6.1.3 Altering the Definition of Exposure----

# 1. Redefining the VAI Percentile (Replacing the Original Quartile)

baseline <- baseline %>%
  mutate(VAI_Q3 = cut(VAI, breaks = quantile(VAI, c(0, 1/3, 2/3, 1), na.rm = TRUE), 
                      labels = c("Q1", "Q2", "Q3"), include.lowest = TRUE)) %>%
  mutate(VAI_trend_Q3 = median(VAI, na.rm = TRUE), .by = VAI_Q3) 

# 2. Refitting the full model using three-digit numbers 

fit_sens3 <- glm(
  co_disease ~ VAI_per_IQR + age + gender + location + marital + smoking + drinking + sbp + dbp + BMI + health11,
  data = baseline, family = binomial()
)

# Subgroup Analysis of Third-Order Numbers
fit_sens3_q <- glm(
  co_disease ~ VAI_Q3 + age + gender + location + marital + smoking + drinking + sbp + dbp + BMI + health11,
  data = baseline, family = binomial()
)

# Extraction Results
sens3_OR_cont <- round(exp(coef(fit_sens3))["VAI_per_IQR"], 3)  # Continuous OR
sens3_OR_q2 <- round(exp(coef(fit_sens3_q))["VAI_Q3Q2"], 3)  # Q2 vs Q1
sens3_OR_q3 <- round(exp(coef(fit_sens3_q))["VAI_Q3Q3"], 3)  # Q3 vs Q1

# Output Comparison
cat("\n===== Sensitivity Analysis 3: Altering the Definition of Exposure =====\n")
cat("Original quartile Q2/Q1 OR：", round(exp(coef(fit8))["VAI_QQ2"], 3), "，Q4/Q1 OR：", round(exp(coef(fit8))["VAI_QQ4"], 3), "\n")
cat("Third Quartile Q2/Q1 OR：", sens3_OR_q2, "，Q3/Q1 OR：", sens3_OR_q3, "，Continuous OR：", sens3_OR_cont, "\n")

#----6.1.4 Summary of Sensitivity Analysis Results----

# 1. Extract the core results of the original model

original_OR_cont <- round(exp(coef(fit7))["VAI_per_IQR"], 3)
original_CI_cont <- paste0(round(exp(confint(fit7))["VAI_per_IQR", 1], 3), "-", 
                           round(exp(confint(fit7))["VAI_per_IQR", 2], 3))
original_p_cont <- round(summary(fit7)$coefficients["VAI_per_IQR", "Pr(>|z|)"], 3)

# Original Quartile Model

original_OR_q2 <- round(exp(coef(fit8))["VAI_QQ2"], 3)
original_OR_q4 <- round(exp(coef(fit8))["VAI_QQ4"], 3)

# 2. Fit the model “excluding health11”

fit_sens2_health11 <- glm(
  co_disease ~ VAI_per_IQR + age + gender + location + marital + smoking + health11 + sbp + dbp + BMI,
  data = baseline, family = binomial()
)
sens2_health11_OR <- round(exp(coef(fit_sens2_health11))["VAI_per_IQR"], 3)
sens2_health11_CI <- paste0(round(exp(confint(fit_sens2_health11))["VAI_per_IQR", 1], 3), "-", 
                            round(exp(confint(fit_sens2_health11))["VAI_per_IQR", 2], 3))
sens2_health11_p <- round(summary(fit_sens2_health11)$coefficients["VAI_per_IQR", "Pr(>|z|)"], 3)

# 3. Fit the model “excluding drinking”
sens2_drinking_OR <- round(sens2_OR, 3)
sens2_drinking_CI <- paste0(round(sens2_CI[1], 3), "-", round(sens2_CI[2], 3))
sens2_drinking_p <- round(sens2_p, 3)

# 4. Extraction of Results from the Three-Part Numerical Model

sens3_cont_CI <- paste0(round(exp(confint(fit_sens3))["VAI_per_IQR", 1], 3), "-", 
                        round(exp(confint(fit_sens3))["VAI_per_IQR", 2], 3))
sens3_cont_p <- round(summary(fit_sens3)$coefficients["VAI_per_IQR", "Pr(>|z|)"], 3)
sens3_q2_CI <- paste0(round(exp(confint(fit_sens3_q))["VAI_Q3Q2", 1], 3), "-", 
                      round(exp(confint(fit_sens3_q))["VAI_Q3Q2", 2], 3))
sens3_q2_p <- round(summary(fit_sens3_q)$coefficients["VAI_Q3Q2", "Pr(>|z|)"], 3)
sens3_q3_CI <- paste0(round(exp(confint(fit_sens3_q))["VAI_Q3Q3", 1], 3), "-", 
                      round(exp(confint(fit_sens3_q))["VAI_Q3Q3", 2], 3))
sens3_q3_p <- round(summary(fit_sens3_q)$coefficients["VAI_Q3Q3", "Pr(>|z|)"], 3)

# 5. Create a summary data frame

sensitivity_summary <- data.frame(
  Sensitivity_Method = c(
    "Original Model (Reference)",
    "Method1: Remove Outliers/High Leverage Points",
    "Method2: Remove Covariate 'health11'",
    "Method2: Remove Covariate 'drinking'",
    "Method3: Exposure as Tertiles (Continuous)",
    "Method3: Exposure as Tertiles (Q2 vs Q1)",
    "Method3: Exposure as Tertiles (Q3 vs Q1)"
  ),
  Exposure_Type = c(
    "Continuous (VAI per IQR)",
    "Continuous (VAI per IQR)",
    "Continuous (VAI per IQR)",
    "Continuous (VAI per IQR)",
    "Continuous (VAI per IQR)",
    "Categorical (Tertiles)",
    "Categorical (Tertiles)"
  ),
  OR = c(
    original_OR_cont,
    round(sens1_OR, 3),
    sens2_health11_OR,
    sens2_drinking_OR,
    sens3_OR_cont,
    sens3_OR_q2,
    sens3_OR_q3
  ),
  CI95 = c(
    original_CI_cont,
    paste0(round(sens1_CI[1], 3), "-", round(sens1_CI[2], 3)),
    sens2_health11_CI,
    sens2_drinking_CI,
    sens3_cont_CI,
    sens3_q2_CI,
    sens3_q3_CI
  ),
  P_value = c(
    original_p_cont,
    round(sens1_p, 3),
    sens2_health11_p,
    sens2_drinking_p,
    sens3_cont_p,
    sens3_q2_p,
    sens3_q3_p
  ),
  stringsAsFactors = FALSE   
)
# 6. Output Specification Table 3 
kable(
  sensitivity_summary,
  caption = "Sensitivity Analysis Results of the Association Between VAI and Risk of arthritis-diabetes comorbidity",
  align = c("l", "l", "c", "c", "c"),
  booktabs = TRUE   
)


##### （2）ROC and AUC ####

# 1.1 the prediction probabilities for each model

roc_data <- baseline %>%
  dplyr::select(
    ID, co_disease,
    VAI_per_IQR, 
    age, gender, marital, location,  # fit4 
    smoking, drinking, health11, sbp, dbp,BMI  # fit7 
  ) %>% 
  mutate(
    prob_unadj = predict(fit1, newdata = ., type = "response"),  
    prob_adj1 = predict(fit4, newdata = ., type = "response"),  
    prob_adj2 = predict(fit7, newdata = ., type = "response")
  ) %>% 
  filter(
    !is.na(co_disease),   
    !is.na(prob_unadj), 
    !is.na(prob_adj1),   
    !is.na(prob_adj2)  
  )

head(roc_data)
nrow(roc_data)
write.csv(roc_data, file = "2011-2020roc_data.csv", row.names = FALSE, na = "")

# 2.1 model 1  unadjusted （CVAI）-ROC/AUC

roc_unadj <- roc(
  response = roc_data$co_disease,   
  predictor = roc_data$prob_unadj,   
  ci = TRUE,         
  direction = "<",                       
  levels = c(0, 1)    
)

# 2.2 model 2 （VAI + demographic）的ROC和AUC
roc_adj1 <- roc(
  response = roc_data$co_disease,
  predictor = roc_data$prob_adj1,
  ci = TRUE,
  direction = "<",
  levels = c(0, 1)
)

# 2.3 model 3 full adjusted ROC/AUC
roc_adj2 <- roc(
  response = roc_data$co_disease,
  predictor = roc_data$prob_adj2,
  ci = TRUE,
  direction = "<",
  levels = c(0, 1)
)

# 2.4 Output AUC results (including 95% confidence intervals)
cat("=== CVAI co——disease AUC  ===\n")

cat("1. unadjusted model 1 ：AUC =", round(roc_unadj$auc, 3), 
    "，95%CI：", round(roc_unadj$ci[1], 3), "-", round(roc_unadj$ci[3], 3), "\n")
cat("2. adjusted model 2 ：AUC =", round(roc_adj1$auc, 3), 
    "，95%CI：", round(roc_adj1$ci[1], 3), "-", round(roc_adj1$ci[3], 3), "\n")
cat("3. full_adjusted model 3 ：AUC =", round(roc_adj2$auc, 3), 
    "，95%CI：", round(roc_adj2$ci[1], 3), "-", round(roc_adj2$ci[3], 3), "\n")



# 3.1 Extract the ROC coordinates of the model 1

coords_unadj <- coords(
  roc_unadj, 
  x = "all",  
  ret = c("threshold", "specificity", "sensitivity")   
) %>%
  mutate(
    FPR = 1 - specificity,  # false positive rate
    TPR = sensitivity,      # true positive rate
    Model = "Unadjusted
AUC = 0.659 95%CI：0.611 - 0.706"   
  )

# 3.2 Extract the ROC coordinates of the model 2
coords_adj1 <- coords(
  roc_adj1, 
  x = "all",
  ret = c("threshold", "specificity", "sensitivity")
) %>%
  mutate(
    FPR = 1 - specificity,
    TPR = sensitivity,
    Model = "Adjusted
AUC = 0.621 95%CI：0.572 - 0.67"   
  )

# 3.3 Extract the ROC coordinates of the model 3
coords_adj2 <- coords(
  roc_adj2, 
  x = "all",
  ret = c("threshold", "specificity", "sensitivity")
) %>%
  mutate(
    FPR = 1 - specificity,
    TPR = sensitivity,
    Model = "Fully Adjusted 
AUC = 0.708  95%CI：0.66 - 0.755 "   
  )

# 3.4 combine

roc_coords_all <- bind_rows(coords_unadj, coords_adj1, coords_adj2)


# 4.1 Plot the ROC curve

roc_plot <- ggplot() +
  geom_line(
    data = roc_coords_all,
    aes(x = FPR, y = TPR, color = Model, linetype = Model),
    linewidth = 1.2  
  ) +
  
  # 2. Plot a random guess line（AUC=0.5）
  geom_abline(
    slope = 1, 
    intercept = 0, 
    color = "#95A5A6",   
    linetype = "dashed", 
    linewidth = 0.8
  ) +
  
  # 3.Annotate AUC values
  
  #annotate(
  # "text",
  #x = 0.6, y = 0.7,
  #label = paste0("Unadjusted\nAUC = ", round(roc_unadj$auc, 3), "\n95% CI: ", 
  #                   round(roc_unadj$ci[1], 3), "-", round(roc_unadj$ci[3], 3)),
  #size = 4, color = "#E74C3C", fontface = "bold"
  # ) +
  #annotate(
  # "text",
  # x = 0.6, y = 0.6,
  #    label = paste0("Adjusted (Demographics)\nAUC = ", round(roc_adj1$auc, 3), "\n95% CI: ", 
  #                   round(roc_adj1$ci[1], 3), "-", round(roc_adj1$ci[3], 3)),
  #    size = 4, color = "#3498DB", fontface = "bold"
  #  ) +
  #  annotate(
  #    "text",
  #    x = 0.6, y = 0.5,
  #    label = paste0("Fully Adjusted\nAUC = ", round(roc_adj2$auc, 3), "\n95% CI: ", 
  #                   round(roc_adj2$ci[1], 3), "-", round(roc_adj2$ci[3], 3)),
  #    size = 4, color = "#2ECC71", fontface = "bold"
  #  ) 
  
  # 4. Set axis labels and titles 
  labs(
    x = "False Positive Rate (1 - Specificity)",  
    y = "True Positive Rate (Sensitivity)",  
    color = "Model Type",  
    linetype = "Model Type"  
  ) +
  
  # 5. Coordinate axis range
  xlim(0, 1) +
  ylim(0, 1) +
  
  # 6. Theme
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold", color = "#2C3E50"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "#7F8C8D"),
    axis.title.x = element_text(size = 12, face = "bold", color = "#2C3E50"),
    axis.title.y = element_text(size = 12, face = "bold", color = "#2C3E50"),
    axis.text = element_text(size = 10, color = "#2C3E50"),
    # Legend
    legend.position = "bottom",
    legend.title = element_text(size = 11, face = "bold", color = "#2C3E50"),
    legend.text = element_text(size = 10, color = "#2C3E50"),
    panel.grid = element_blank()
  )

# 4.2 Display ROC curve

print(roc_plot)

# 4.3 Save ROC curve 
ggsave(
  filename = "ROC curve for VAI and comorbidity with arthritis and diabetes.tiff",  
  plot = roc_plot,
  width = 10,  
  height = 8, 
  dpi = 300,    
  device = "tiff"  
)


#### （3）plot table2####

# Extract Model 1
model1_vai_or <- exp(coef(fit1))["VAI_per_IQR"]
model1_vai_ci <- exp(confint(fit1))["VAI_per_IQR", ]
model1_vai_p <- coef(summary(fit1))["VAI_per_IQR", "Pr(>|z|)"]

model1_q1_or <- NA  
model1_q1_ci <- NA
model1_q1_p <- NA

model1_q2_or <- exp(coef(fit2))["VAI_QQ2"]
model1_q2_ci <- exp(confint(fit2))["VAI_QQ2", ]
model1_q2_p <- coef(summary(fit2))["VAI_QQ2", "Pr(>|z|)"]

model1_q3_or <- exp(coef(fit2))["VAI_QQ3"]
model1_q3_ci <- exp(confint(fit2))["VAI_QQ3", ]
model1_q3_p <- coef(summary(fit2))["VAI_QQ3", "Pr(>|z|)"]

model1_q4_or <- exp(coef(fit2))["VAI_QQ4"]
model1_q4_ci <- exp(confint(fit2))["VAI_QQ4", ]
model1_q4_p <- coef(summary(fit2))["VAI_QQ4", "Pr(>|z|)"]

model1_trend_p <- coef(summary(fit3))["VAI_trend", "Pr(>|z|)"]

# Extract Model 2 
model2_vai_or <- exp(coef(fit4))["VAI_per_IQR"]
model2_vai_ci <- exp(confint(fit4))["VAI_per_IQR", ]
model2_vai_p <- coef(summary(fit4))["VAI_per_IQR", "Pr(>|z|)"]

model2_q1_or <- NA
model2_q1_ci <- NA
model2_q1_p <- NA

model2_q2_or <- exp(coef(fit5))["VAI_QQ2"]
model2_q2_ci <- exp(confint(fit5))["VAI_QQ2", ]
model2_q2_p <- coef(summary(fit5))["VAI_QQ2", "Pr(>|z|)"]

model2_q3_or <- exp(coef(fit5))["VAI_QQ3"]
model2_q3_ci <- exp(confint(fit5))["VAI_QQ3", ]
model2_q3_p <- coef(summary(fit5))["VAI_QQ3", "Pr(>|z|)"]

model2_q4_or <- exp(coef(fit5))["VAI_QQ4"]
model2_q4_ci <- exp(confint(fit5))["VAI_QQ4", ]
model2_q4_p <- coef(summary(fit5))["VAI_QQ4", "Pr(>|z|)"]

model2_trend_p <- coef(summary(fit6))["VAI_trend", "Pr(>|z|)"]

# Extract Model 3 
model3_vai_or <- exp(coef(fit7))["VAI_per_IQR"]
model3_vai_ci <- exp(confint(fit7))["VAI_per_IQR", ]
model3_vai_p <- coef(summary(fit7))["VAI_per_IQR", "Pr(>|z|)"]

model3_q1_or <- NA
model3_q1_ci <- NA
model3_q1_p <- NA

model3_q2_or <- exp(coef(fit8))["VAI_QQ2"]
model3_q2_ci <- exp(confint(fit8))["VAI_QQ2", ]
model3_q2_p <- coef(summary(fit8))["VAI_QQ2", "Pr(>|z|)"]

model3_q3_or <- exp(coef(fit8))["VAI_QQ3"]
model3_q3_ci <- exp(confint(fit8))["VAI_QQ3", ]
model3_q3_p <- coef(summary(fit8))["VAI_QQ3", "Pr(>|z|)"]

model3_q4_or <- exp(coef(fit8))["VAI_QQ4"]
model3_q4_ci <- exp(confint(fit8))["VAI_QQ4", ]
model3_q4_p <- coef(summary(fit8))["VAI_QQ4", "Pr(>|z|)"]

model3_trend_p <- coef(summary(fit9))["VAI_trend", "Pr(>|z|)"]

# Constructing a Table Data Frame

table_data <- data.frame(
  Variable = c(
    "VAI per IQR",
    "Quartiles of VAI", "Q1", "Q2", "Q3", "Q4", "p for trend"
  ),
  Model_1_OR_CI = c(
    sprintf("%.2f [%.2f, %.2f]", model1_vai_or, model1_vai_ci[1], model1_vai_ci[2]),
    "", "Ref", 
    sprintf("%.2f [%.2f, %.2f]", model1_q2_or, model1_q2_ci[1], model1_q2_ci[2]),
    sprintf("%.2f [%.2f, %.2f]", model1_q3_or, model1_q3_ci[1], model1_q3_ci[2]),
    sprintf("%.2f [%.2f, %.2f]", model1_q4_or, model1_q4_ci[1], model1_q4_ci[2]),
    ""
  ),
  Model_1_p = c(
    ifelse(model1_vai_p < 0.001, "<0.001", sprintf("%.3f", model1_vai_p)),
    "", "", 
    ifelse(model1_q2_p < 0.001, "<0.001", sprintf("%.3f", model1_q2_p)),
    ifelse(model1_q3_p < 0.001, "<0.001", sprintf("%.3f", model1_q3_p)),
    ifelse(model1_q4_p < 0.001, "<0.001", sprintf("%.3f", model1_q4_p)),
    ifelse(model1_trend_p < 0.001, "<0.001", sprintf("%.3f", model1_trend_p))
  ),
  Model_2_OR_CI = c(
    sprintf("%.2f [%.2f, %.2f]", model2_vai_or, model2_vai_ci[1], model2_vai_ci[2]),
    "", "Ref", 
    sprintf("%.2f [%.2f, %.2f]", model2_q2_or, model2_q2_ci[1], model2_q2_ci[2]),
    sprintf("%.2f [%.2f, %.2f]", model2_q3_or, model2_q3_ci[1], model2_q3_ci[2]),
    sprintf("%.2f [%.2f, %.2f]", model2_q4_or, model2_q4_ci[1], model2_q4_ci[2]),
    ""
  ),
  Model_2_p = c(
    ifelse(model2_vai_p < 0.001, "<0.001", sprintf("%.3f", model2_vai_p)),
    "", "", 
    ifelse(model2_q2_p < 0.001, "<0.001", sprintf("%.3f", model2_q2_p)),
    ifelse(model2_q3_p < 0.001, "<0.001", sprintf("%.3f", model2_q3_p)),
    ifelse(model2_q4_p < 0.001, "<0.001", sprintf("%.3f", model2_q4_p)),
    ifelse(model2_trend_p < 0.001, "<0.001", sprintf("%.3f", model2_trend_p))
  ),
  Model_3_OR_CI = c(
    sprintf("%.2f [%.2f, %.2f]", model3_vai_or, model3_vai_ci[1], model3_vai_ci[2]),
    "", "Ref", 
    sprintf("%.2f [%.2f, %.2f]", model3_q2_or, model3_q2_ci[1], model3_q2_ci[2]),
    sprintf("%.2f [%.2f, %.2f]", model3_q3_or, model3_q3_ci[1], model3_q3_ci[2]),
    sprintf("%.2f [%.2f, %.2f]", model3_q4_or, model3_q4_ci[1], model3_q4_ci[2]),
    ""
  ),
  Model_3_p = c(
    ifelse(model3_vai_p < 0.001, "<0.001", sprintf("%.3f", model3_vai_p)),
    "", "", 
    ifelse(model3_q2_p < 0.001, "<0.001", sprintf("%.3f", model3_q2_p)),
    ifelse(model3_q3_p < 0.001, "<0.001", sprintf("%.3f", model3_q3_p)),
    ifelse(model3_q4_p < 0.001, "<0.001", sprintf("%.3f", model3_q4_p)),
    ifelse(model3_trend_p < 0.001, "<0.001", sprintf("%.3f", model3_trend_p))
  )
)

# print table2
knitr::kable(table_data, 
             caption = "Association of VAI with the risk of Co-Arthritis&Diabetes",
             align = c("l", "l", "c", "l", "c", "l", "c"))


# ###7. RCS：Exploring the Nonlinear Relationship Between VAI and Co-Disease ####

ddist <- datadist(baseline)
options(datadist = "ddist")

#  
rcs_fit <- lrm(co_disease ~ rcs(VAI, 4) + age + gender + sbp + dbp + education +
                 location + marital + drinking + smoking +BMI + health11 ,
               data = baseline, x = TRUE, y = TRUE)

print(rcs_fit)

# reference value (median)

ref_value <- median(baseline$VAI, na.rm = TRUE)

# Create drawing data

plot_data <- rms::Predict(rcs_fit, VAI, fun = plogis, ref.zero = TRUE)

# plot RCS

ggplot(plot_data, aes(x = VAI, y = yhat)) +
  geom_line(color = "blue", linewidth = 1.5) +  
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "blue") +
  geom_hline(yintercept = plogis(0), linetype = "dashed", color = "red") +  
  geom_vline(xintercept = ref_value, linetype = "dashed", color = "gray") +
  labs(x = "VAI", 
       y = "OR（95％CI）for Co-disease Risk",
       title = "RCS of the Relationship between VAI and CO-Arthritis & Diabetes") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.caption = element_blank())

ggsave("Fig.3(A) Restrictive Cubic Spline (RCS) for VAI and Arthritis-Diabetes Comorbidity.tif", 
       width = 8, height = 6, dpi = 300, compression = "lzw")



## Testing Nonlinear Relationships（似然比检验Likelihood Ratio Test）
linear_fit <- lrm(co_disease ~ VAI + age + gender + education +
                    location + marital, 
                  data = baseline, x = TRUE, y = TRUE)

# Likelihood Ratio Test

anova_result <- anova(rcs_fit)
print(anova_result)

# Extract p-values for nonlinear tests

nonlinear_p <- anova_result[" Nonlinear", "P"]
cat("nonlinear_p:", nonlinear_p, "\n")

# Model Comparison（AIC（Akaike Information Criterion，赤池信息准则））
AIC_values <- c(
  "Linear model" = AIC(linear_fit),
  "4-Node RCS" = AIC(rcs_fit)
)

cat("AIC values for different models:\n")
print(AIC_values)

# which.min(AIC_values)=best model 

best_model <- names(which.min(AIC_values))
cat("Best Model:", best_model, "\n")

# Calculate OR 95％CI

or_data <- rms::Predict(rcs_fit, VAI, fun = exp, ref.zero = TRUE)

# plot OR RCS

ggplot(or_data, aes(x = VAI, y = yhat)) +
  geom_line(color = "red", linewidth = 1.5) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "red") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  geom_vline(xintercept = ref_value, linetype = "dashed", color = "gray") +
  labs(x = "VAI", 
       y = "OR（95％CI）for Co-disease Risk") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

# Seeking the Inflection Point(拐点)

# 通过计算OR值曲线的导数变化来识别

vai_values <- or_data$VAI
or_values <- or_data$yhat

derivative1 <- diff(or_values) / diff(vai_values)

# 寻找导数变化最大的点
inflection_idx <- which.max(abs(diff(derivative1))) + 1
inflection_point <- vai_values[inflection_idx]

cat("建议拐点位置:", inflection_point, "\n")

# Mark the inflection point

ggplot(or_data, aes(x = VAI, y = yhat)) +
  geom_line(color = "red", size = 1.5) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "red") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  geom_vline(xintercept = inflection_point, linetype = "dashed", color = "blue") +
  annotate("text", x = inflection_point, y = max(or_data$yhat), 
           label = paste("inflection point:", round(inflection_point, 1)), 
           vjust = -1, color = "blue") +
  labs(x = "VAI", 
       y = "OR（95％CI）for Co-disease Risk") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.caption = element_blank())


ggsave("VAI_arthritis_RCS_plot.png", width = 8, height = 6, dpi = 300)
ggsave("Fig.3(B) Restrictive Cubic Spline (RCS) with inflection point.tif", 
       width = 8, height = 6, dpi = 300, compression = "lzw")

# Output Results Summary
cat("=== Summary of Constrained Cubic Spline Analysis Results ===\n")
cat("Reference VAI value:", ref_value, "\n")
cat("Nonlinear Test p-value:", nonlinear_p, "\n")
cat("Best Model:", best_model, "\n")
cat("Recommended inflection point location:", inflection_point, "\n")
cat("The graphic has been saved as: VAI_arthritis_RCS_plot.png 和 Fig.3(B) Restrictive Cubic Spline (RCS) with inflection point.tif\n")


# ###8. Hierarchical analysis ####

#----8.1 Gender ----

# Male subgroup
male_co_disease <- baseline %>% filter(gender == 'Male')
fit_male <- glm(co_disease ~ VAI_per_IQR + age_group +drinking + health11 +
                  smoking + dyslipidemia11 + heart_diseases11,
                data = male_co_disease, family = binomial())
summary(fit_male)
OR_male <- exp(coef(fit_male))
CI_male <- exp(confint(fit_male))
cat("男性亚组OR值：", OR_male, "\n")
cat("男性亚组95%置信区间：", CI_male, "\n")

# Female subgroup
female_aco_disease <- baseline %>% filter(gender == 'Female')
fit_female <- glm(co_disease ~ VAI_per_IQR + age_group +drinking + health11 +
                    smoking + dyslipidemia11 + heart_diseases11,
                  data = female_aco_disease, family = binomial())
summary(fit_female)
OR_female <- exp(coef(fit_female))
CI_female <- exp(confint(fit_female))
cat("女性亚组OR值：", OR_female, "\n")
cat("女性亚组95%置信区间：", CI_female, "\n")

# 检验交互作用（CVAI与性别）
fit_interaction <- glm(co_disease ~ VAI_per_IQR * gender+ age_group +drinking + health11 +
                         smoking + dyslipidemia11 + heart_diseases11, 
                       data = baseline, family = binomial())
summary(fit_interaction)
OR_interaction <- exp(coef(fit_interaction))
CI_interaction <- exp(confint(fit_interaction))
cat("交互作用模型OR值：", OR_interaction, "\n")
cat("交互作用模型95%置信区间：", CI_interaction, "\n")

#---- 8.2 Age_group ----
# <60 
lower_co_disease <- baseline %>% filter(age_group == '<60')
fit_lower <- glm(co_disease ~ VAI_per_IQR + gender + drinking + health11 +
                   smoking + dyslipidemia11 + heart_diseases11,
                 data = lower_co_disease, family = binomial())
summary(fit_lower)
OR_lower <- exp(coef(fit_lower))
CI_lower <- exp(confint(fit_lower))
cat("<60亚组OR值：", OR_lower, "\n")
cat("<60亚组95%置信区间：", CI_lower, "\n")

# >=60
upper_co_disease <- baseline %>% filter(age_group == '>=60')
fit_upper<- glm(co_disease ~ VAI_per_IQR + gender + drinking + health11 +
                  smoking + dyslipidemia11 + heart_diseases11,
                data = upper_co_disease, family = binomial())
summary(fit_upper)
OR_upper <- exp(coef(fit_upper))
CI_upper <- exp(confint(fit_upper))
cat(">=60亚组OR值：", OR_upper, "\n")
cat(">=60亚组95%置信区间：", CI_upper, "\n")

# 检验交互作用（VAI与Age_group）
fit_interaction <- glm(co_disease ~ VAI_per_IQR * age_group + gender + drinking + health11 +
                         smoking + dyslipidemia11 + heart_diseases11,
                       data = baseline, family = binomial())
summary(fit_interaction)
OR_interaction <- exp(coef(fit_interaction))
CI_interaction <- exp(confint(fit_interaction))
cat("交互作用模型OR值：", OR_interaction, "\n")
cat("交互作用模型95%置信区间：", CI_interaction, "\n")

#----8.3 health_status----
# Good
good_co_disease <- baseline %>% filter(health11 == 'Good')
fit_good <- glm(co_disease ~ VAI_per_IQR + gender + age_group +drinking +
                  smoking + dyslipidemia11 + heart_diseases11,
                data = good_co_disease, family = binomial())
summary(fit_good)
OR_good <- exp(coef(fit_good))
CI_good <- exp(confint(fit_good))
cat("good亚组OR值：", OR_good, "\n")
cat("good亚组95%置信区间：", CI_good, "\n")

# Fair
fair_co_disease <- baseline %>% filter(health11 == 'Fair')
fit_fair <- glm(co_disease ~ VAI_per_IQR + gender + age_group +drinking +
                  smoking + dyslipidemia11 + heart_diseases11,  
                data = fair_co_disease, family = binomial())
summary(fit_fair)
OR_fair <- exp(coef(fit_fair))
CI_fair <- exp(confint(fit_fair))
cat("fair亚组OR值：", OR_fair, "\n")
cat("fair亚组95%置信区间：", CI_fair, "\n")

# Poor
poor_co_disease <- baseline %>% filter(health11 == 'Poor')
fit_poor <- glm(co_disease ~ VAI_per_IQR + gender + age_group +drinking +
                  smoking + dyslipidemia11 + heart_diseases11,  
                data = poor_co_disease, family = binomial())
summary(fit_poor)
OR_poor <- exp(coef(fit_poor))
CI_poor <- exp(confint(fit_poor))
cat("poor亚组OR值：", OR_poor, "\n")
cat("poor亚组95%置信区间：", CI_poor, "\n")


# 检验交互作用 
fit_interaction <- glm(co_disease ~ VAI_per_IQR *health11 + age_group + gender + drinking  +
                         smoking + dyslipidemia11 + heart_diseases11,
                       data = baseline, family = binomial())
summary(fit_interaction)
OR_interaction <- exp(coef(fit_interaction))
CI_interaction <- exp(confint(fit_interaction))
cat("交互作用模型OR值：", OR_interaction, "\n")
cat("交互作用模型95%置信区间：", CI_interaction, "\n")

#----8.4 Drinking ----
# Drinker 亚组
drinker_co_disease <- baseline %>% filter(drinking == 'drinker')
fit_drinker <- glm(co_disease ~ VAI_per_IQR + gender + age_group + health11 +
                     smoking + dyslipidemia11 + heart_diseases11,
                   data = drinker_co_disease, family = binomial())
summary(fit_drinker)
OR_drinker <- exp(coef(fit_drinker))
CI_drinker <- exp(confint(fit_drinker))
cat("饮酒者亚组OR值：", OR_drinker, "\n")
cat("饮酒者亚组95%置信区间：", CI_drinker, "\n")

# Nondrinker 亚组
nondrinker_co_disease <- baseline %>% filter(drinking == 'nondrinker')
fit_nondrinker <- glm(co_disease ~ VAI_per_IQR + gender + age_group + health11 +
                        smoking + dyslipidemia11 + heart_diseases11,
                      data = nondrinker_co_disease, family = binomial())
summary(fit_nondrinker)
OR_nondrinker <- exp(coef(fit_nondrinker))
CI_nondrinker <- exp(confint(fit_nondrinker))
cat("非饮酒者亚组OR值：", OR_nondrinker, "\n")
cat("非饮酒者亚组95%置信区间：", CI_nondrinker, "\n")

# 检验交互作用（VAI_per_IQR 与 drinking）
fit_interaction_drinking <- glm(co_disease ~ VAI_per_IQR * drinking + gender + age_group + health11 +
                                  smoking + dyslipidemia11 + heart_diseases11,
                                data = baseline, family = binomial())
summary(fit_interaction_drinking)
OR_interaction_drinking <- exp(coef(fit_interaction_drinking))
CI_interaction_drinking <- exp(confint(fit_interaction_drinking))
cat("饮酒交互作用模型OR值：", OR_interaction_drinking, "\n")
cat("饮酒交互作用模型95%置信区间：", CI_interaction_drinking, "\n")


#---- 8.5 Smoking ----
# Smoker 亚组
smoker_co_disease <- baseline %>% filter(smoking == 'smoker')
fit_smoker <- glm(co_disease ~ VAI_per_IQR + gender + age_group + health11 +
                    drinking + dyslipidemia11 + heart_diseases11,
                  data = smoker_co_disease, family = binomial())
summary(fit_smoker)
OR_smoker <- exp(coef(fit_smoker))
CI_smoker <- exp(confint(fit_smoker))
cat("吸烟者亚组OR值：", OR_smoker, "\n")
cat("吸烟者亚组95%置信区间：", CI_smoker, "\n")

# Nonsmoker 亚组
nonsmoker_co_disease <- baseline %>% filter(smoking == 'nonsmoker')
fit_nonsmoker <- glm(co_disease ~ VAI_per_IQR + gender + age_group + health11 +
                       drinking + dyslipidemia11 + heart_diseases11,
                     data = nonsmoker_co_disease, family = binomial())
summary(fit_nonsmoker)
OR_nonsmoker <- exp(coef(fit_nonsmoker))
CI_nonsmoker <- exp(confint(fit_nonsmoker))
cat("非吸烟者亚组OR值：", OR_nonsmoker, "\n")
cat("非吸烟者亚组95%置信区间：", CI_nonsmoker, "\n")

# 检验交互作用（VAI_per_IQR 与 smoking）
fit_interaction_smoking <- glm(co_disease ~ VAI_per_IQR * smoking + gender + age_group + health11 +
                                 drinking + dyslipidemia11 + heart_diseases11,
                               data = baseline, family = binomial())
summary(fit_interaction_smoking)
OR_interaction_smoking <- exp(coef(fit_interaction_smoking))
CI_interaction_smoking <- exp(confint(fit_interaction_smoking))
cat("吸烟交互作用模型OR值：", OR_interaction_smoking, "\n")
cat("吸烟交互作用模型95%置信区间：", CI_interaction_smoking, "\n")


#----8.6 Dyslipidemia11 ----
# Normal 亚组
normal_co_disease <- baseline %>% filter(dyslipidemia11 == 'normal')
fit_normal <- glm(co_disease ~ VAI_per_IQR + gender + age_group + health11 +
                    drinking + smoking + heart_diseases11,
                  data = normal_co_disease, family = binomial())
summary(fit_normal)
OR_normal <- exp(coef(fit_normal))
CI_normal <- exp(confint(fit_normal))
cat("血脂正常亚组OR值：", OR_normal, "\n")
cat("血脂正常亚组95%置信区间：", CI_normal, "\n")

# Abnormal 亚组
abnormal_co_disease <- baseline %>% filter(dyslipidemia11 == 'abnormal')
fit_abnormal <- glm(co_disease ~ VAI_per_IQR + gender + age_group + health11 +
                      drinking + smoking + heart_diseases11,
                    data = abnormal_co_disease, family = binomial())
summary(fit_abnormal)
OR_abnormal <- exp(coef(fit_abnormal))
CI_abnormal <- exp(confint(fit_abnormal))
cat("血脂异常亚组OR值：", OR_abnormal, "\n")
cat("血脂异常亚组95%置信区间：", CI_abnormal, "\n")

# 检验交互作用（VAI_per_IQR 与 dyslipidemia11）
fit_interaction_dyslipidemia <- glm(co_disease ~ VAI_per_IQR * dyslipidemia11 + gender + age_group + health11 +
                                      drinking + smoking + heart_diseases11,
                                    data = baseline, family = binomial())
summary(fit_interaction_dyslipidemia)
OR_interaction_dyslipidemia <- exp(coef(fit_interaction_dyslipidemia))
CI_interaction_dyslipidemia <- exp(confint(fit_interaction_dyslipidemia))
cat("血脂交互作用模型OR值：", OR_interaction_dyslipidemia, "\n")
cat("血脂交互作用模型95%置信区间：", CI_interaction_dyslipidemia, "\n")


#---- 8.7 Heart_diseases11 ----
# Heart_diseases 亚组
heart_co_disease <- baseline %>% filter(heart_diseases11 == 'heart_diseases')
fit_heart <- glm(co_disease ~ VAI_per_IQR + gender + age_group + health11 +
                   drinking + smoking + dyslipidemia11,
                 data = heart_co_disease, family = binomial())
summary(fit_heart)
OR_heart <- exp(coef(fit_heart))
CI_heart <- exp(confint(fit_heart))
cat("心脏病亚组OR值：", OR_heart, "\n")
cat("心脏病亚组95%置信区间：", CI_heart, "\n")

# Non_heart_diseases 亚组
non_heart_co_disease <- baseline %>% filter(heart_diseases11 == 'non_heart_diseases')
fit_non_heart <- glm(co_disease ~ VAI_per_IQR + gender + age_group + health11 +
                       drinking + smoking + dyslipidemia11,
                     data = non_heart_co_disease, family = binomial())
summary(fit_non_heart)
OR_non_heart <- exp(coef(fit_non_heart))
CI_non_heart <- exp(confint(fit_non_heart))
cat("无心脏病亚组OR值：", OR_non_heart, "\n")
cat("无心脏病亚组95%置信区间：", CI_non_heart, "\n")

# 检验交互作用（VAI_per_IQR 与 heart_diseases11）
fit_interaction_heart <- glm(co_disease ~ VAI_per_IQR * heart_diseases11 + gender + age_group + health11 +
                               drinking + smoking + dyslipidemia11,
                             data = baseline, family = binomial())
summary(fit_interaction_heart)
OR_interaction_heart <- exp(coef(fit_interaction_heart))
CI_interaction_heart <- exp(confint(fit_interaction_heart))
cat("心脏病交互作用模型OR值：", OR_interaction_heart, "\n")
cat("心脏病交互作用模型95%置信区间：", CI_interaction_heart, "\n")



#### 9.Create forest map data ####

#----9.1 Data Preparation ----

# 1.Calculate the overall population effect

fit_total <- glm(co_disease ~ VAI_per_IQR + heart_diseases11 + gender + age_group + health11 +
                   drinking + smoking + dyslipidemia11,
                 data = baseline, family = binomial())
OR_total <- exp(coef(fit_total))["VAI_per_IQR"]
CI_total <- exp(confint(fit_total))["VAI_per_IQR", ]
P_total <- summary(fit_total)$coefficients["VAI_per_IQR", "Pr(>|z|)"]


plot_data <- data.frame(
  # First column：Group Categories
  category = c(
    "General population",
    rep("Stratified by smoking status", 2),
    rep("Stratified by gender", 2),
    rep("Stratified by age_group", 2),
    rep("Stratified by health status", 3),
    rep("Stratified by dyslipidemia ", 2),
    rep("Stratified by drinking", 2),
    rep("Stratified by heart_diseases ", 2)
  ),
  # Second column：subgroup names
  subcategory = c(
    "VAI（每IQR增加）",
    "smoker", "nonsmoker",
    "male", "female",
    "<60", ">=60",
    "Good", "Fair","Poor",
    "abnormal","normal",
    "drinker","nondrinker",
    "heart_diseases","non_heart_diseases"
  ),
  # Third column：OR
  OR = c(
    OR_total,  
    OR_smoker["VAI_per_IQR"], OR_nonsmoker["VAI_per_IQR"],  
    OR_male["VAI_per_IQR"], OR_female["VAI_per_IQR"], 
    OR_lower["VAI_per_IQR"], OR_upper["VAI_per_IQR"], 
    OR_good["VAI_per_IQR"], OR_fair["VAI_per_IQR"], OR_poor["VAI_per_IQR"], 
    OR_abnormal["VAI_per_IQR"], OR_normal["VAI_per_IQR"],  
    OR_drinker["VAI_per_IQR"], OR_nondrinker["VAI_per_IQR"], 
    OR_heart["VAI_per_IQR"], OR_non_heart["VAI_per_IQR"]  
  ),
  # Fourth column：Lower limit of the 95% confidence interval
  lower = c(
    CI_total[1], 
    CI_smoker["VAI_per_IQR", 1], CI_nonsmoker["VAI_per_IQR", 1],  
    CI_male["VAI_per_IQR", 1], CI_female["VAI_per_IQR", 1],  
    CI_lower["VAI_per_IQR", 1], CI_upper["VAI_per_IQR", 1],  
    CI_good["VAI_per_IQR", 1], CI_fair["VAI_per_IQR", 1], CI_poor["VAI_per_IQR", 1],  
    CI_abnormal["VAI_per_IQR", 1], CI_normal["VAI_per_IQR", 1], 
    CI_drinker["VAI_per_IQR", 1], CI_nondrinker["VAI_per_IQR", 1],  
    CI_heart["VAI_per_IQR", 1], CI_non_heart["VAI_per_IQR", 1]  
  ),
  # Fifth column：Upper limit of 95% confidence interval
  upper = c(
    CI_total[2],  
    CI_smoker["VAI_per_IQR", 2], CI_nonsmoker["VAI_per_IQR", 2],  
    CI_male["VAI_per_IQR", 2], CI_female["VAI_per_IQR", 2],  
    CI_lower["VAI_per_IQR", 2], CI_upper["VAI_per_IQR", 2],  
    CI_good["VAI_per_IQR", 2], CI_fair["VAI_per_IQR", 2], CI_poor["VAI_per_IQR", 2],  
    CI_abnormal["VAI_per_IQR", 2], CI_normal["VAI_per_IQR", 2],  
    CI_drinker["VAI_per_IQR", 2], CI_nondrinker["VAI_per_IQR", 2], 
    CI_heart["VAI_per_IQR", 2], CI_non_heart["VAI_per_IQR", 2]  
  ),
  # Sixth column：P value
  P = c(
    P_total,  
    summary(fit_smoker)$coefficients["VAI_per_IQR", "Pr(>|z|)"], 
    summary(fit_nonsmoker)$coefficients["VAI_per_IQR", "Pr(>|z|)"],  
    summary(fit_male)$coefficients["VAI_per_IQR", "Pr(>|z|)"], 
    summary(fit_female)$coefficients["VAI_per_IQR", "Pr(>|z|)"], 
    summary(fit_lower)$coefficients["VAI_per_IQR", "Pr(>|z|)"], 
    summary(fit_upper)$coefficients["VAI_per_IQR", "Pr(>|z|)"],  
    summary(fit_good)$coefficients["VAI_per_IQR", "Pr(>|z|)"], 
    summary(fit_fair)$coefficients["VAI_per_IQR", "Pr(>|z|)"], 
    summary(fit_poor)$coefficients["VAI_per_IQR", "Pr(>|z|)"],  
    summary(fit_abnormal)$coefficients["VAI_per_IQR", "Pr(>|z|)"], 
    summary(fit_normal)$coefficients["VAI_per_IQR", "Pr(>|z|)"], 
    summary(fit_drinker)$coefficients["VAI_per_IQR", "Pr(>|z|)"], 
    summary(fit_nondrinker)$coefficients["VAI_per_IQR", "Pr(>|z|)"],  
    summary(fit_heart)$coefficients["VAI_per_IQR", "Pr(>|z|)"], 
    summary(fit_non_heart)$coefficients["VAI_per_IQR", "Pr(>|z|)"]  
  )
)

print(plot_data)

# ---- 9.2 Drawing a forest map  ----


str(plot_data)

sum(is.na(plot_data$OR))
sum(is.na(plot_data$lower))
sum(is.na(plot_data$upper))

summary(plot_data$OR)
summary(plot_data$lower)
summary(plot_data$upper)

# Construct data
plot_data <- data.frame(
  category = c(
    "Overall",
    # Smoking Status Grouping
    "smoking status",
    "smoker",
    "nonsmoker",
    # Gender Grouping
    "gender",
    "male",
    "female",
    # Age Groups
    "age-group",
    "<60",
    ">=60",
    # Self-Rated Health Status Grouping
    "Self-Rated health status",
    "Good",
    "Fair",
    "Poor",
    # Blood Lipid Level Grouping
    "dyslipidemia level",
    "abnormal",
    "normal",
    # Drinking Status Grouping
    "drinking status",
    "drinker",
    "nondrinker",
    # Cardiac Status Grouping
    "heart disease",
    "heart disease",
    "non_heart disease"
  ),
  number = c(
    NA,  # Overall
    NA, 1876, 2966, 
    NA, 2305, 2542, 
    NA, 2985, 1862,  
    NA, 735, 1138, 500,  
    NA, 503, 4344,  
    NA, 1649, 3193,  
    NA, 416, 4401   
  ),
  OR = c(
    1.142,  # Overall
    NA, 1.099, 1.142,  
    NA, 1.103, 1.157, 
    NA, 1.127, 1.125,  
    NA, 1.137, 1.117, 1.166, 
    NA, 1.139, 1.123, 
    NA, 1.094, 1.160,  
    NA, 1.411, 1.115   
  ),
  lower = c(
    1.062,  # Overall
    NA, 0.925, 1.063,  
    NA, 0.973, 1.065, 
    NA, 1.045, 0.968, 
    NA, 1.012, 1.003, 1.007,  
    NA, 1.020, 1.032, 
    NA, 0.955, 1.068,  
    NA, 1.097, 1.034   
  ),
  upper = c(
    1.227,  # Overall
    NA, 1.243, 1.231, 
    NA, 1.217, 1.265, 
    NA, 1.215, 1.272,  
    NA, 1.296, 1.228, 1.327,  
    NA, 1.303, 1.210, 
    NA, 1.211, 1.270,  
    NA, 1.995, 1.193  
  ),
  P_interaction = c(
    NA,  # Overall
    0.710, NA, NA,  
    0.479, NA, NA, 
    0.926, NA, NA, 
    0.766, NA, NA, NA, 
    0.759, NA, NA, 
    0.393, NA, NA, 
    0.270, NA, NA  
  )
)

plot_data <- plot_data %>%
  mutate(
    n_percent = ifelse(!is.na(number), as.character(number), ""),
    or_ci_text = ifelse(!is.na(OR), sprintf("%.3f (%.3f, %.3f)", OR, lower, upper), ""),
    p_text = ifelse(!is.na(P_interaction), sprintf("%.2f", P_interaction), "")
  )

labeltext <- cbind(
  c("Variable", plot_data$category),       
  c("n", plot_data$n_percent),             
  c("OR (95% CI)", plot_data$or_ci_text),  
  c("P for interaction", plot_data$p_text)  
)

is_summary <- c(TRUE, is.na(plot_data$OR)) 

mean_values <- c(NA, plot_data$OR)          
lower_values <- c(NA, plot_data$lower)      
upper_values <- c(NA, plot_data$upper)      

tiff("forest_plot_with_p_value.tif", width = 7.5, height = 6, units = "in", res = 300)

# Drawing a forest map
forestplot(
  labeltext = labeltext,
  mean = mean_values,
  lower = lower_values,
  upper = upper_values,
  
  graph.pos = 3,              
  align = c("l", "c", "c", "c"),
  
  is.summary = is_summary,     
  col = fpColors(
    box = "royalblue",         
    lines = "darkblue",        
    summary = "royalblue",      
    zero = "black"            
  ),
  
  txt_gp = fpTxtGp(
    label = gpar(cex = 0.7),          
    summary = gpar(cex = 0.7, fontface = "bold"), 
    xlab = gpar(cex = 0.9),             
    ticks = gpar(cex = 0.7),           
    title = gpar(cex = 1.0, fontface = "bold")    
  ),
  xlab = "Odds Ratio (OR)",
  xlog = FALSE,                 
  zero = 1,                     
  clip = c(0.5, 2.0),           
  xticks = c(0.5, 0.7, 1.0, 1.5, 2.0),  
  grid = structure(c(1), gp = gpar(lty = 2, col = "gray")),  
  
  boxsize = 0.2,               
  lineheight = "auto",          
  colgap = unit(6, "mm"),       
  line.margin = unit(1, "mm"),  
)

dev.off()

