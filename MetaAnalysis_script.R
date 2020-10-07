#Copyright (c) 2020 by Alexandra Therond. All Rights Reserved.


#The objective of this code is to:
#Calculate the difference in Hedge's G per observation
#Compute the pooled effect size and variance per study
#Run the model on global cognition
#Conduct the sub-group analysis by cognitive domain
#Conduct the mixed-effects model meta-regression on the four moderators

#required packages
install.packages("metafor")
install.packages("dplyr")
install.packages("forestplot")
library(metafor)
library(dplyr)
library(forestplot)

#the metafor package is employed to compute the Hedge’s G, as well as run the global cognition model, the sub-group analysis by cognitive domain, the mixed-effects model meta-regression and the funnel plot
#the dplyr package is used to more efficiently group the assessments by study and by study/cognitive domain
#the forestplot package is employed to produce the forest plots


#set r correlation coefficient
# we decided to use a conservative threshold of 0.5

r_coef <- 0.5

#import dataset (After publication the .csv file will be freely available upon request at https://www.cranilab.com/open-science)

import_data <- read.csv("~/Desktop/ CRT_depression.csv ", stringsAsFactors=FALSE, fileEncoding="latin1")
R_FINAL <- data.frame(import_data)
head(R_FINAL)

#Loop to add negative if mean decrease is a positive effect (”Dec” in the column Inc_or_Dec)

for (i in 1: nrow(R_FINAL)){
  if(R_FINAL$Inc_or_Dec[[i]] == "Dec"){
    R_FINAL$pre_mean_T[[i]] = - R_FINAL$pre_mean_T[[i]]
    R_FINAL$pre_mean_C[[i]] = - R_FINAL$pre_mean_C[[i]]
    R_FINAL$post_mean_T[[i]] = - R_FINAL$post_mean_T[[i]]
    R_FINAL$post_mean_C[[i]] = - R_FINAL$post_mean_C[[i]]
  }
}

#calculating the control means, SD, ni(sample size), ri (r coef 0.5, rep= replicate) to be entered in the escalc function

raw_data_C <- data.frame(
  m_pre   = c(R_FINAL$pre_mean_C),
  m_post  = c(R_FINAL$post_mean_C),
  sd_pre  = c(R_FINAL$pre_SD_C),
  sd_post = c(R_FINAL$post_SD_C),
  ni      = c(R_FINAL$sample_size_C),
  ri      = c(rep(r_coef, nrow(R_FINAL))))

#calculating the treatment means, SD, ni(sample size), ri (r coef 0.5, rep= replicate) to be entered in the escalc function

raw_data_T <- data.frame(
  m_pre   = c(R_FINAL$pre_mean_T),
  m_post  = c(R_FINAL$post_mean_T),
  sd_pre  = c(R_FINAL$pre_SD_T),
  sd_post = c(R_FINAL$post_SD_T),
  ni      = c(R_FINAL$sample_size_T),
  ri      = c(rep(r_coef, nrow(R_FINAL))))

#calculating the effect size for control and treatment - yi (effect size) and vi (corresponding Variance)

data_C <- escalc(measure="SMCR", m1i=m_post, m2i=m_pre, sd1i=sd_pre, ri=ri, ni=ni, data=raw_data_C)

data_T <- escalc(measure="SMCR", m1i=m_post, m2i=m_pre, sd1i=sd_pre, ri=ri, ni=ni, data=raw_data_T)


#calculating the difference between treatment and control effect size for each assessment of each study (Tyi -Cyi) 
#dat = name of the new dataset with mean and variance  for each assessment


dat <- data.frame(yi = data_T$yi - data_C$yi, vi = data_T$vi + data_C$vi) 

#create dataset with moderators 
#mod_data = name of the new dataset which includes moderators 

mod_data <- data.frame(study_name = R_FINAL$study_name, cog_domain = R_FINAL$cog_domain,
                       Group = R_FINAL$Group, 
                       Duration_hrs = R_FINAL$Duration_hrs, Avg_age = R_FINAL$Avg_age)

mod_data <- cbind(mod_data, dat)

###############################
######## Overall cognition ########
###############################

#Aim 1: Run the model on global cognition

#calculate pooled effect size per study
#create groups per study

groups <- group_by(mod_data, study_name) 

#Using the summarise function: computing the mean of the yi by stydy name, it takes all the effect sizes and computes on mean effect size difference by study

means <- data.frame(summarise(groups,mean_yi = mean(yi)))

#create function to calculate pooled vi (variance), as per Borenstein (2019)
#the equation is stored in pooled_vi_fct in order to reuse it later 

#formula is 1/m^2 * sum(vi + sum r*sqrt(vi)*sqrt(vj) where i not equal j) for i from 1 to m
#m is number of observations within one study

pooled_vi_fct <- function(study_vi){
  if (length(study_vi) > 1){ #compute pooled only if more than 1 observation per study
    out <- 0 #initiate pooled vi to 0
    for(i in 1:length(study_vi)){ #start loop that will go over all vi (1 to length(vi))
      vi <- study_vi[i] #start one vi at a time
      remaining_vi <- study_vi[study_vi != vi] #list of all remaining vi for study
      out <- out + vi #begin by adding value of vi to pooled vi 
      for (j in 1:length(remaining_vi)){ #loop over remaining vi
        out <- out + r_coef*sqrt(vi)*sqrt(remaining_vi[j]) #add r*sqrt(vi)*sqrt(vj) to pooled vi
      }
    }
    pooledvi <- out / (length(study_vi)^2) #take result of sum and divide by m^2
    return(pooledvi) #output of function is final pooled vi for that study
  } else {return (study_vi)} #if only one vi per study, return vi
}

#run pooled variance function over groups
vars <- data.frame(summarise(groups, mean_vi = pooled_vi_fct(vi)))

#combining effect size and vars: so for every study we have an effect size size and variance for the global cognition effect of CRT

avg_study <- merge(means, vars, by="study_name")

#run pooled random effects model on global cognition 

random_avg <- rma(mean_yi, mean_vi, method = "REML", data = avg_study, slab = study_name)
summary(random_avg)

##Vizualize the effect of CRT on global cognition 
#create the meta data to produce the boxes and confidence intervals on the forest plot
meta_data <- 
  structure(list(
    mean  = c(NA, c(random_avg$yi), random_avg$b), 
    lower = c(NA, c(random_avg$yi - 1.96*sqrt(random_avg$vi)), random_avg$ci.lb),
    upper = c(NA, c(random_avg$yi + 1.96*sqrt(random_avg$vi)), random_avg$ci.ub)),
    .Names = c("mean", "lower", "upper"), 
    row.names = c(NA, -11L), 
    class = "data.frame")

#compute the lower and upper bound of the confidence interval to be displayed in the plot
low_CI <- round(c(random_avg$yi - 1.96*sqrt(random_avg$vi), random_avg$ci.lb),2)
high_CI <- round(c(random_avg$yi + 1.96*sqrt(random_avg$vi), random_avg$ci.ub),2)
CI <- paste("[", low_CI, ", ", high_CI, "]", sep="")

#create tabletext, containing the text that will be included in the plot
tabletext<-cbind(
  c("Study", as.character(random_avg$slab), "Random Effects Model"),
  c("Hedge's G [95% CI]", paste(round(random_avg$yi,2), " ", CI[1:9]), 
    paste(round(random_avg$b,2), " ", CI[10])))


#create forest plot
forestplot(tabletext,
           meta_data,
           new_page = TRUE,
           graph.pos = 2,
           is.summary=c(TRUE,rep(FALSE,9),TRUE),
           hrzl_lines = TRUE,
           txt_gp = fpTxtGp(ticks = gpar(cex=0.8)),
           boxsize = 0.25,
           title = "Global Cognition",
           grid = structure(c(random_avg$b), gp = gpar(col = "steelblue", lty=2)),
           col = fpColors(box="grey50",line="black",zero = "grey41", summary = "steelblue"))

#########################
#####Cognitive domains#####
#########################

#Aim 2: Conduct the subgroup analysis by cognitive domains


#calculate average effect size per cognitive domain per study
#group_by: for the following calculations mod_data (defined above, individual effect sizes) will be grouped by study name and cog domain 

groups <- group_by(mod_data, study_name, cog_domain)
means_cog <- data.frame(summarise(groups,mean_yi_cog = mean(yi)))
vars_cog <- data.frame(summarise(groups, mean_vi_cog = pooled_vi_fct(vi)))

#merge into the same dataset: the pooled domain effect size and pooled variance

avg_study_cog <- merge(means_cog, vars_cog, by=c("study_name", "cog_domain"))

#remove observations with IQ because there are not enough studies to proceed with the subgroup analysis (n=2)

avg_study_cog <- avg_study_cog[avg_study_cog$cog_domain != "IQ ",]


##run subgroup analysis for cognitive domains similar as what Cella et al., 2020 did
# run one model for each cognitive domain

#att creating a data set, subseting only attention and speedprocessing
att <- avg_study_cog[avg_study_cog$cog_domain == "Attention_SpeedProcessing",]

#running random effects model for attention
att_avg <- rma(mean_yi_cog, mean_vi_cog, method = "REML", data = att, slab = study_name)
summary(att_avg)

#creating a data set, subseting only executive functioning
exec <- avg_study_cog[avg_study_cog$cog_domain == "Executive_Functioning",]

#running random effects model for executive functioning
exec_avg <- rma(mean_yi_cog, mean_vi_cog, method = "REML", data = exec, slab = study_name)
summary(exec_avg)

#creating a data set, subseting only verbal fluency
verbf <- avg_study_cog[avg_study_cog$cog_domain == "Verbal_Fluency",]

#running random effects model for verbal fluency
verbf_avg <- rma(mean_yi_cog, mean_vi_cog, method = "REML", data = verbf, slab = study_name)
summary(verbf_avg)

#creating a data set, subseting only verbal memory
verbm <- avg_study_cog[avg_study_cog$cog_domain == "Verbal_Memory",]

#running random effects model for verbal memory
verbm_avg <- rma(mean_yi_cog, mean_vi_cog, method = "REML", data = verbm, slab = study_name)
summary(verbm_avg)

#creating a data set, subseting only working memory
work <- avg_study_cog[avg_study_cog$cog_domain == "Working_Memory",]

#running random effects model for working memory
work_avg <- rma(mean_yi_cog, mean_vi_cog, method = "REML", data = work, slab = study_name)
summary(work_avg)

#creating a data set, subseting only visuospatial memory
vis <- avg_study_cog[avg_study_cog$cog_domain == "Visuospatial_Memory",]

#running random effects model for working memory
vis_avg <- rma(mean_yi_cog, mean_vi_cog, method = "REML", data = vis, slab = study_name)
summary(vis_avg)

###Vizualise the effect of CRT on global cognition 

### attention/speed processing ###
#create the meta data to produce the boxes and confidence intervals on the forest plot
meta_data_att <- 
  structure(list(
    mean  = c(NA, c(att_avg$yi), att_avg$b), 
    lower = c(NA, c(att_avg$yi - 1.96*sqrt(att_avg$vi)), att_avg$ci.lb),
    upper = c(NA, c(att_avg$yi + 1.96*sqrt(att_avg$vi)), att_avg$ci.ub)),
    .Names = c("mean", "lower", "upper"), 
    row.names = c(NA, -8L), 
    class = "data.frame")

#compute the lower and upper bound of the confidence interval to be displayed in the plot
low_CI <- round(c(att_avg$yi - 1.96*sqrt(att_avg$vi), att_avg$ci.lb),2)
high_CI <- round(c(att_avg$yi + 1.96*sqrt(att_avg$vi), att_avg$ci.ub),2)
CI <- paste("[", low_CI, ", ", high_CI, "]", sep="")

#create tabletext, containing the text that will be included in the plot
tabletext_att<-cbind(
  c("Study", as.character(att_avg$slab), "Random Effects Model"),
  
  c("Hedge's G [95% CI]", paste(round(att_avg$yi,2), " ", CI[1:6]), 
    paste(round(att_avg$b,2), " ", CI[7])))

#create forest plot
forestplot(tabletext_att,
           meta_data_att,
           new_page = TRUE,
           graph.pos = 2,
           is.summary=c(TRUE,rep(FALSE,6),TRUE),
           hrzl_lines = TRUE,
           txt_gp = own,
           boxsize = 0.25,
           title = "Attention/Processing Speed",
           grid = structure(c(att_avg$b), gp = gpar(col = "steelblue", lty=2)),
           col = fpColors(box="grey50",line="black",zero = "grey41", summary = "steelblue"))

### verbal memory ###
#create the meta data to produce the boxes and confidence intervals on the forest plot
meta_data_verbm <- 
  structure(list(
    mean  = c(NA, c(verbm_avg$yi), verbm_avg$b), 
    lower = c(NA, c(verbm_avg$yi - 1.96*sqrt(verbm_avg$vi)), verbm_avg$ci.lb),
    upper = c(NA, c(verbm_avg$yi + 1.96*sqrt(verbm_avg$vi)), verbm_avg$ci.ub)),
    .Names = c("mean", "lower", "upper"), 
    row.names = c(NA, -9L), 
    class = "data.frame")

#compute the lower and upper bound of the confidence interval to be displayed in the plot
low_CI <- round(c(verbm_avg$yi - 1.96*sqrt(verbm_avg$vi), verbm_avg$ci.lb),2)
high_CI <- round(c(verbm_avg$yi + 1.96*sqrt(verbm_avg$vi), verbm_avg$ci.ub),2)
CI <- paste("[", low_CI, ", ", high_CI, "]", sep="")

#create tabletext, containing the text that will be included in the plot
tabletext_verbm<-cbind(
  c("Study", as.character(verbm_avg$slab), "Random Effects Model"),
  
  c("Hedge's G [95% CI]", paste(round(verbm_avg$yi,2), " ", CI[1:7]), 
    paste(round(verbm_avg$b,2), " ", CI[8])))

#create forest plot
forestplot(tabletext_verbm,
           meta_data_verbm,
           new_page = TRUE,
           graph.pos = 2,
           txt_gp = own,
           is.summary=c(TRUE,rep(FALSE,7),TRUE),
           hrzl_lines = TRUE,
           boxsize = 0.25,
           title = "Verbal Memory",
           grid = structure(c(verbm_avg$b), gp = gpar(col = "steelblue", lty=2)),
           col = fpColors(box="grey50",line="black",zero = "grey41", summary = "steelblue"))

### working memory ###
#create the meta data to produce the boxes and confidence intervals on the forest plot
meta_data_work <- 
  structure(list(
    mean  = c(NA, c(work_avg$yi), work_avg$b), 
    lower = c(NA, c(work_avg$yi - 1.96*sqrt(work_avg$vi)), work_avg$ci.lb),
    upper = c(NA, c(work_avg$yi + 1.96*sqrt(work_avg$vi)), work_avg$ci.ub)),
    .Names = c("mean", "lower", "upper"), 
    row.names = c(NA, -7L), 
    class = "data.frame")

#compute the lower and upper bound of the confidence interval to be displayed in the plot
low_CI <- round(c(work_avg$yi - 1.96*sqrt(work_avg$vi), work_avg$ci.lb),2)
high_CI <- round(c(work_avg$yi + 1.96*sqrt(work_avg$vi), work_avg$ci.ub),2)
CI <- paste("[", low_CI, ", ", high_CI, "]", sep="")

#create tabletext, containing the text that will be included in the plot
tabletext_work<-cbind(
  c("Study", as.character(work_avg$slab), "Random Effects Model"),
  
  c("Hedge's G [95% CI]", paste(round(work_avg$yi,2), " ", CI[1:5]), 
    paste(round(work_avg$b,2), " ", CI[6])))

#create forest plot
forestplot(tabletext_work,
           meta_data_work,
           new_page = TRUE,
           graph.pos = 2,
           txt_gp = own,
           is.summary=c(TRUE,rep(FALSE,5),TRUE),
           hrzl_lines = TRUE,
           boxsize = 0.25,
           title = "Working Memory",
           grid = structure(c(work_avg$b), gp = gpar(col = "steelblue", lty=2)),
           col = fpColors(box="grey50",line="black",zero = "grey41", summary = "steelblue"))

### executive functioning ###
#create the meta data to produce the boxes and confidence intervals on the forest plot
meta_data_exec <- 
  structure(list(
    mean  = c(NA, c(exec_avg$yi), exec_avg$b), 
    lower = c(NA, c(exec_avg$yi - 1.96*sqrt(exec_avg$vi)), exec_avg$ci.lb),
    upper = c(NA, c(exec_avg$yi + 1.96*sqrt(exec_avg$vi)), exec_avg$ci.ub)),
    .Names = c("mean", "lower", "upper"), 
    row.names = c(NA, -7L), 
    class = "data.frame")

#compute the lower and upper bound of the confidence interval to be displayed in the plot
low_CI <- round(c(exec_avg$yi - 1.96*sqrt(exec_avg$vi), exec_avg$ci.lb),2)
high_CI <- round(c(exec_avg$yi + 1.96*sqrt(exec_avg$vi), exec_avg$ci.ub),2)
CI <- paste("[", low_CI, ", ", high_CI, "]", sep="")

#create tabletext, containing the text that will be included in the plot
tabletext_exec<-cbind(
  c("Study", as.character(exec_avg$slab), "Random Effects Model"),
  
  c("Hedge's G [95% CI]", paste(round(exec_avg$yi,2), " ", CI[1:5]), 
    paste(round(exec_avg$b,2), " ", CI[6])))

#create forest plot
forestplot(tabletext_exec,
           meta_data_exec,
           new_page = TRUE,
           graph.pos = 2,
           txt_gp = own,
           is.summary=c(TRUE,rep(FALSE,5),TRUE),
           hrzl_lines = TRUE,
           boxsize = 0.25,
           title = "Executive Functioning",
           grid = structure(c(exec_avg$b), gp = gpar(col = "steelblue", lty=2)),
           col = fpColors(box="grey50",line="black",zero = "grey41", summary = "steelblue"))

### verbal fluency ###
#create the meta data to produce the boxes and confidence intervals on the forest plot
meta_data_verbf <- 
  structure(list(
    mean  = c(NA, c(verbf_avg$yi), verbf_avg$b), 
    lower = c(NA, c(verbf_avg$yi - 1.96*sqrt(verbf_avg$vi)), verbf_avg$ci.lb),
    upper = c(NA, c(verbf_avg$yi + 1.96*sqrt(verbf_avg$vi)), verbf_avg$ci.ub)),
    .Names = c("mean", "lower", "upper"), 
    row.names = c(NA, -5L), 
    class = "data.frame")

#compute the lower and upper bound of the confidence interval to be displayed in the plot
low_CI <- round(c(verbf_avg$yi - 1.96*sqrt(verbf_avg$vi), verbf_avg$ci.lb),2)
high_CI <- round(c(verbf_avg$yi + 1.96*sqrt(verbf_avg$vi), verbf_avg$ci.ub),2)
CI <- paste("[", low_CI, ", ", high_CI, "]", sep="")

#create tabletext, containing the text that will be included in the plot
tabletext_verbf<-cbind(
  c("Study", as.character(verbf_avg$slab), "Random Effects Model"),
  
  c("Hedge's G [95% CI]", paste(round(verbf_avg$yi,2), " ", CI[1:3]), 
    paste(round(verbf_avg$b,2), " ", CI[4])))

#create forest plot
forestplot(tabletext_verbf,
           meta_data_verbf,
           new_page = TRUE,
           graph.pos = 2,
           txt_gp = own,
           is.summary=c(TRUE,rep(FALSE,3),TRUE),
           hrzl_lines = TRUE,
           boxsize = 0.25,
           title = "Verbal Fluency",
           grid = structure(c(verbf_avg$b), gp = gpar(col = "steelblue", lty=2)),
           col = fpColors(box="grey50",line="black",zero = "grey41", summary = "steelblue"))

### visuospatial memory ###
#create the meta data to produce the boxes and confidence intervals on the forest plot
meta_data_vis <- 
  structure(list(
    mean  = c(NA, c(vis_avg$yi), vis_avg$b), 
    lower = c(NA, c(vis_avg$yi - 1.96*sqrt(vis_avg$vi)), vis_avg$ci.lb),
    upper = c(NA, c(vis_avg$yi + 1.96*sqrt(vis_avg$vi)), vis_avg$ci.ub)),
    .Names = c("mean", "lower", "upper"), 
    row.names = c(NA, -5L), 
    class = "data.frame")

#compute the lower and upper bound of the confidence interval to be displayed in the plot
low_CI <- round(c(vis_avg$yi - 1.96*sqrt(vis_avg$vi), vis_avg$ci.lb),2)
high_CI <- round(c(vis_avg$yi + 1.96*sqrt(vis_avg$vi), vis_avg$ci.ub),2)
CI <- paste("[", low_CI, ", ", high_CI, "]", sep="")

#create tabletext, containing the text that will be included in the plot
tabletext_vis<-cbind(
  c("Study", as.character(vis_avg$slab), "Random Effects Model"),
  
  c("Hedge's G [95% CI]", paste(round(vis_avg$yi,2), " ", CI[1:3]), 
    paste(round(vis_avg$b,2), " ", CI[4])))

#create forest plot
forestplot(tabletext_vis,
           meta_data_vis,
           new_page = TRUE,
           graph.pos = 2,
           txt_gp = own,
           is.summary=c(TRUE,rep(FALSE,3),TRUE),
           hrzl_lines = TRUE,
           boxsize = 0.25,
           title = "Visuospatial Memory",
           grid = structure(c(verbf_avg$b), gp = gpar(col = "steelblue", lty=2)),
           col = fpColors(box="grey50",line="black",zero = "grey41", summary = "steelblue"))

#################################
######## MODERATOR ANALYSISModerators analysis #######
################################

#creating a dataset with the moderators by adding moderators to dataset (avg_study)
moderators <- unique(data.frame(study_name = mod_data$study_name, group = mod_data$Group,
                                duration = mod_data$Duration_hrs,
                                age = mod_data$Avg_age))
#Alvarez has two CRT groups with two different mean age that need to be combined
moderators <- unique(data.frame(moderators %>%
                                  group_by(study_name) %>%
                                  mutate(age = mean(age)) %>%
                                  ungroup()))
#adding the mean effect size and variance to the dataset (avg_study)
avg_study <- merge(avg_study, moderators, by = "study_name")

#Aim #3: Conduct the mixed-effects model meta-regression on the four moderators 

#Moderator: GROUP vs INDIVIDUAL

#moderator analysis: hypothesis testing whether each group is sig different from each other
mixed_group2 <- rma(mean_yi, mean_vi, mods = ~ group, method = "REML", data = avg_study)
summary(mixed_group2)

#post-hoc subgroup analysis: hypothesis testing whether each group is sig different than zero
mixed_group1 <- rma(mean_yi, mean_vi, mods = ~ group - 1, method = "REML", data = avg_study)
summary(mixed_group1)

#Moderator: TREATMENT DURATION
mixed_dur <- rma(mean_yi, mean_vi, mods = ~ duration, method = "REML", data = avg_study)
summary(mixed_dur)

#Moderator: PARTICIPANT’S AGE
mixed_age <- rma(mean_yi, mean_vi, mods = ~ age, method = "REML", data = avg_study)
summary(mixed_age)



########################################
########### Creation of funnel plot ##########
########################################

#define colours for funnel plot
colours <-  c("red","dodgerblue3","palevioletred1","springgreen4",
              "steelblue1","gray50","lightgreen","orange","maroon","black")
#produce funnel plot
funnel(random_avg,level=c(95, 99), shade=c("white", "gray88"), col = colours, back="gray96", pch=20, cex.lab=0.85)
#create legend
legend("topright", legend = unique(mod_data$study_name),col=unique(colours), pch = 16, cex = 0.65)


##References
#Borenstein, M. (2019). Common mistakes in meta-analysis and how to avoid them. New Jersey, USA: Biostats Incorporated.


##Aknowledgement
#Annik Gougeon
