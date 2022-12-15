rm(list=ls())
library(dplyr)
library(survival)
library(pubh)
library(ccoptimalmatch)
library(ggplot2)

load(file="drugdata.rds")  #drugdata from the WCE package
source("CXO_funcs.R")
fu = 90

drugdata <- drugdata %>%
  group_by(Id) %>%
  mutate(futime=max(Stop),
         ex=as.numeric(dose > 0)) %>%
  ungroup()

last <- drugdata %>% 
  group_by(Id) %>%
  slice_tail(n=1) %>%
  filter(futime >= fu) %>%
  ungroup()
  
caseids <- last %>% 
  filter(Event == 1)  %>%
  mutate(tc=0) 


##case_time controls

tcids <-  last %>% 
  filter( Event == 0)  %>%
  mutate(tc=1) 
 



##create a CXO study with 90 day time window


cases <- inner_join(drugdata , caseids %>% select(Id), by="Id") %>%
  filter(futime >= fu) %>%
  group_by(Id) %>%
  mutate(minex=min(ex), maxex=max(ex),
         concordant = minex==maxex,
         day=Stop-futime+fu,
         wt=1) %>%
  filter(!concordant, day>0) %>%
  ungroup()



##create case-time-control data

create_subset <- caseids %>%
  arrange( sex) %>%
  distinct( sex,  .keep_all = TRUE) %>%
  mutate(subset = 1:n()) %>%
  select( sex,  subset)


case_with_subset <- full_join(caseids, create_subset)


tc_with_subset <- right_join(tcids, create_subset)

cases_tc <- rbind(case_with_subset,tc_with_subset)

table(cases_tc$tc)



##match time-controls to cases

bdd_controls <- cases_tc %>% 
  filter(tc==1) %>%
  mutate(cluster_case = 0)

bdd_cases <- cases_tc %>% filter(tc==0)
bdd_cases$cluster_case <- paste("case",1:nrow(bdd_cases),sep = "_")


not_processed <- rbind(bdd_cases,bdd_controls)


bdd_temp <- data.frame()
list_p <- unique(bdd_cases$cluster_case)


for(i in 1:length(list_p)){
  temp <- bdd_cases[bdd_cases$cluster_case==list_p[i],]
  subset_identified <- temp$subset
  temp0 <- bdd_controls[bdd_controls$subset==temp$subset,]
  temp_final <- rbind(temp,temp0)
  temp_final$cluster_case <- list_p[i]
  temp_final=temp_final %>%
    group_by(cluster_case) %>%
    mutate(age_diff = abs(age - age[tc==0]),
           fup_diff = futime - futime[tc==0])
  temp_final$age_fup <- ifelse(temp_final$age_diff<=5 & temp_final$fup_diff>= 0,"accept","delete")
  temp_final <- temp_final[temp_final$age_fup=="accept",]
  temp_final$age_fup <- NULL
  bdd_temp <- rbind(bdd_temp,temp_final)
}

table(bdd_temp$tc)


bdd_temp = bdd_temp %>% group_by(cluster_case) %>% mutate(total_control_per_case = n()-1)
bdd_temp$case_ind <- ifelse(bdd_temp$tc==1,1,0)
bdd_temp <- subset(bdd_temp, select=c(cluster_case, Id, tc, case_ind,
                                      sex, age_diff, fup_diff, total_control_per_case))

bdd_temp = bdd_temp %>% group_by(Id) %>% mutate(freq_of_controls = n())

bdd_temp<-bdd_temp[order(bdd_temp$cluster_case,bdd_temp$tc,bdd_temp$fup_diff,
                         bdd_temp$age_diff,bdd_temp$freq_of_controls),]

hist(bdd_temp$freq_of_controls)

final_data <- optimal_matching(bdd_temp, n_con=2, cluster_case, Id, 
                               total_control_per_case, tc, with_replacement = T)

final_data = final_data %>% group_by(cluster_case) %>% 
  mutate(total_control_matched = n()-1) %>%
  ungroup()
table(final_data$tc,final_data$total_control_matched)


##use the IDs to create the final dataset

final_data <- final_data %>% mutate(newId = row_number())

casetimecontrols <- inner_join(drugdata, final_data %>% select(Id, newId, tc), by="Id") %>%
  select(!Id) %>%
  rename(Id=newId) %>%
  filter(futime >= fu) %>%
  group_by(Id) %>%
  mutate(minex=min(ex), maxex=max(ex),
         concordant = minex==maxex,
         day=Stop-futime+fu,
         wt=1) %>%
  filter(!concordant, day>0) %>%
  ungroup()


to_plot <- casetimecontrols %>%
  mutate(tc=factor(tc)) %>%
  group_by(tc, day) %>%
  summarise(ex = mean(ex, na.rm=T), .groups='drop')

to_plot %>%
  ggplot(aes(x=day, y=ex, group=tc, colour=tc)) +
  geom_point()

##save the data
save(drugdata, cases,casetimecontrols, file="CXO.Rdata")

##conditional logistic regression

cfit <- clogit(Event ~ ex + strata(Id) , data=cases,method="efron") #+ offset(wt)
summary(cfit)
exp(cbind(coef(cfit), confint(cfit)))

## M-H estimate


mhor(formula = Event ~ Id/ex, data=cases) 
(bias <- SCL_bias(data=cases, exposure=ex, event=Event, Id=Id))


## simulation data
for (i in 1:8) {
  fn=paste0("Data/sampdata_scenario0", i, ".csv")
  
  temp <- read.csv(fn)
  
  temp <- temp %>% 
    mutate(event = case_period * (1-tc)) %>%
    group_by(Pt_ID) %>%
    mutate(period=row_number()) %>%
    ungroup() 
  
  if (i < 5) {
    temp <- temp %>%
    select(Pt_ID, ex, event, period, tc)
  } else {
    temp <- temp %>%
      select(Pt_ID, ex, event, period, z, tc)
  }
  
  assign(paste0("ctcsim", i), temp)
    
               
}

save(ctcsim1, ctcsim2, ctcsim3, ctcsim4, ctcsim5, ctcsim6, ctcsim7, ctcsim8,
     file="Data/ctcsim.Rdata")



