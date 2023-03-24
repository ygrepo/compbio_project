#install.packages(c("survival", "lubridate", "ggsurvfit", "gtsummary", "tidycmprsk"))
#remotes::install_github("zabore/condsurv")
# remotes::install_github("zabore/ezfun")
library(survival)
library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)
library(condsurv)
library(dplyr)

lung <- 
  lung %>% 
  mutate(
    status = recode(status, `1` = 0, `2` = 1)
  )

Surv(lung$time, lung$status)[1:10]
s1 <- survfit(Surv(time, status) ~ 1, data = lung)
str(s1)

survfit2(Surv(time, status) ~ 1, data = lung) %>% 
  ggsurvfit() +
  labs(
    x = "Days",
    y = "Overall survival probability"
  ) + 
  add_confidence_interval() +
  add_risktable()

survfit(Surv(time, status) ~ 1, data = lung) %>% 
  tbl_survfit(
    times = 365.25,
    label_header = "**1-year survival (95% CI)**"
  )

survdiff(Surv(time, status) ~ sex, data = lung)

# Cox regression model

coxph(Surv(time, status) ~ sex, data = lung)
coxph(Surv(time, status) ~ sex, data = lung) %>% 
  tbl_regression(exp = TRUE) 

