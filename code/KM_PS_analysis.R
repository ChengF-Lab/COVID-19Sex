setwd("PATH")
library(MatchIt)
library(Matching)
library(gtools)
library(ggplot2)
library(survival)
library(survminer)
library(dplyr)

cov_ps <- read.csv("symptom_hospi.csv")

#=========crude data, cumulative hazard=========
fit_df <- survfit(Surv(hopital_time, hospitalized) ~ gender,conf.int =0.95,data = cov_ps)
#pvalue
res <- pairwise_survdiff(Surv(hopital_time, hospitalized) ~ gender,data = cov_ps)
res
#plot
plt_df<-ggsurvplot(fit_df, data = cov_ps,
                   fun="cumhaz",
                   legend.title = "Gender",
                   legend = "right",
                   ylab = "Cumulative hazard",
                   xlab = "Days since COVID19 symptom/testing onset",
                   title="Hospitalization",
                   legend.labs = c("Male", "Female"),
                   conf.int = TRUE,
                   censor = FALSE,
                   risk.table = TRUE,
                   risk.table.col = "strata",# Risk table color by groups
                   risk.table.y.text = FALSE,
                   break.x.by=15,
                   ggtheme = theme_bw(),
                   palette =c("#DA4453","#3BAFDA"),
                   pval = TRUE,
                   pval.size=8,
                   pval.coord = c(2, 0.05),
                   font.legend=16,
                   font.main = 20,
                   font.x =  18,
                   font.y = 18,
                   font.tickslab = 18)
plt_df

##pdf("./Cumulative_Hazard/symptom_icu_CH.pdf",width=5, height=4)
#print(plt_df$plot) # need print plot first
#dev.off()

#hazard ratio
cov_ps$gender = as.factor(cov_ps$gender)
cov_ps$gender = relevel(cov_ps$gender, ref = "1")
fit_adjust <- coxph(Surv(hopital_time, hospitalized) ~ gender, 
                    data = cov_ps)
sink('./HR/symptom_hopital_raw_summary.txt')
print(res)
print("===============")
summary(fit_adjust)
sink()


#=========PS match, cumulative hazard=========
#PS match
set.seed(1234)
cov_ps_sub<-subset(cov_ps,smoking_now <=1)
cov_ps_sub<-subset(cov_ps_sub,coronary_artery_disease <=1)
cov_ps_sub<-subset(cov_ps_sub,diabetes <=1)
cov_ps_sub<-subset(cov_ps_sub,hypertension <=1)
cov_ps_sub<-subset(cov_ps_sub,copd_emphysema <=1)
#match cov19 positive
match_it<-matchit(hospitalized ~ age+race_other+black+white+smoking_now+
                    coronary_artery_disease+diabetes+hypertension+copd_emphysema,
                  data=cov_ps_sub, method = 'nearest',ratio=1)
summary(match_it)

#Saving the matched samples
df_match <- match.data(match_it)[1:ncol(cov_ps_sub)]

#matchlist
covid_np_mch<-CrossTable(df_match$gender,df_match$hospitalized, 
                         prop.r = F,prop.c = F,prop.t = F,prop.chisq = F)
covid_np_mch_count_list<-c(covid_np_mch$t[4:4],covid_np_mch$t[3:3],covid_np_mch$t[2:2],covid_np_mch$t[1:1])

#cumulative hazard
summary(df_match)

fit_df_match <- survfit(Surv(hopital_time, hospitalized) ~ gender,conf.int =0.95,data = df_match)
#pvalue
res_muh <- pairwise_survdiff(Surv(hopital_time, hospitalized) ~ gender,data = df_match)
res_muh
#plot
plt_df_match<-ggsurvplot(fit_df_match, data = df_match,
                         fun="cumhaz",
                         legend.title = "Gender",
                         legend = "right",
                         ylab = "Cumulative hazard",
                         xlab = "Days since COVID19 symptom onset",
                         title="Hospitalization",
                         legend.labs = c( "Female","Male"),
                         conf.int = TRUE,
                         censor = FALSE,
                         risk.table = TRUE,
                         risk.table.col = "strata",# Risk table color by groups
                         risk.table.y.text = FALSE,
                         break.x.by=15,
                         ggtheme = theme_bw(),
                         palette =c("#3BAFDA","#DA4453"),
                         pval = TRUE,
                         pval.size=8,
                         pval.coord = c(2, 0.05),
                         font.legend=16,
                         font.main = 20,
                         font.x =  18,
                         font.y = 18,
                         font.tickslab = 18)
plt_df_match

#pdf("./Cumulative_Hazard/symptom_icu_CH.pdf",width=5, height=4)
#print(plt_df_match$plot) # need print plot first
#dev.off()

#=======hazard ratio================
df_match$gender = as.factor(df_match$gender)
df_match$gender = relevel(df_match$gender, ref = "1")
fit_adjust <- coxph(Surv(hopital_time, hospitalized) ~ gender, 
                    data = df_match)
sink('./HR/symptom_hospital_PSmch_summary.txt')
print(res_muh)
print("===============")
summary(fit_adjust)
sink()



