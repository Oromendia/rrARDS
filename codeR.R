# Clara Oromendia
# March 22, 2017
# Augustine Choi, Ed Schenck, Kevin Ma
# Merge 7 trials to study ARDS: rapid resolving vs persistent


dir = "O:/Main/Research/Research_share/Clara Oromendia/Pulmonology/Rapid Resolving ARDS"
dirF = paste0(dir,"/Reports/Figures")
outTexReport = paste0(dir,"/Reports/Report.tex")
source("O:/Divisional/BIO/mco2004/myFunctions.R")

# library(rpart)
# library(rpart.plot)
# library(rattle)
# library(pROC)
# library(OptimalCutpoints)
library(haven) # read SAS datasets 
col_rrARDS = brewer.pal(3,"Set2")[-3]

#############
# ARMA
#############
cat("\n \\newpage \\section{ARMA Study} \n",file=outTexReport,append=T)
# Baseline demographics and severity of illness 
dat_screen_net1 = data.frame(read_sas(paste0(dir,"/Data_2017_0411/ARMA/data/ARDSNet01/screen.sas7bdat")))
dat_screen_net2 = data.frame(read_sas(paste0(dir,"/Data_2017_0411/ARMA/data/ARDSNet03/screen.sas7bdat")))
dat_screen = merge(dat_screen_net1,dat_screen_net2,all=T)
head(dat_screen)
dim(dat_screen)
colnames(dat_screen) = tolower(colnames(dat_screen))

dat_mine = subset(dat_screen,select=c(ptid,age,gender,ethnic,pneum,sepsis,aspir,trauma,other,multran))


# BMI 
dat_vital = merge(
	data.frame(read_sas(paste0(dir,"/Data_2017_0411/ARMA/data/ARDSNet01/vital.sas7bdat"))),
	data.frame(read_sas(paste0(dir,"/Data_2017_0411/ARMA/data/ARDSNet03/vital.sas7bdat"))),
	all=T)
dat_bmi = subset(dat_vital, VIS==0,select=c(PTID,WEIGHTK,HEIGHTC))
dat_bmi$ptid = dat_bmi$PTID
dat_bmi$bmi = dat_bmi$WEIGHTK/dat_bmi$HEIGHTC

dat_mine = merge(dat_mine,subset(dat_bmi,select=c(PTID,bmi)),all.x=T,by.x="ptid",by.y="PTID")
dim(dat_mine)

# Vasopressor use 
dat_bruss_net1 = data.frame(read_sas(paste0(dir,"/Data_2017_0411/ARMA/data/ARDSNet01/bruss.sas7bdat")))
dat_bruss_net2 = data.frame(read_sas(paste0(dir,"/Data_2017_0411/ARMA/data/ARDSNet03/bruss.sas7bdat")))
dat_bruss = merge(dat_bruss_net1,dat_bruss_net2,all=T)
colnames(dat_bruss) = tolower(colnames(dat_bruss))
dat_bruss_day0 = subset(dat_bruss, vdate==0) # 3 ppl missing day 0 bruss info 

dim(dat_bruss_day0)
head(dat_bruss_day0)
table(dat_bruss$vdate)

# PAFI 
table(dat_bruss_day0$pafi >= 300) # 13 had pafi above 300
length(unique(dat_bruss$ptid[dat_bruss_day0$pafi>= 300]))


##############################################
# Define rapid resolving ARDS
##############################################
# PAFI above 300 day 0 or 1, or extubated day 1 
dat_bruss_day0 = subset(dat_bruss, vdate==0,select=c(ptid,pafi)) 
dat_bruss_day1 = subset(dat_bruss, vdate==1,select=c(ptid,pafi)) 
dat_bruss_maxPafi = merge(dat_bruss_day0,dat_bruss_day1,by="ptid",all=T)
dat_bruss_maxPafi$maxPAFI = pmax(dat_bruss_maxPafi$pafi.x, dat_bruss_maxPafi$pafi.y,na.rm=T)
dat_bruss_maxPafi$maxPAFI_300plus = dat_bruss_maxPafi$maxPAFI >= 300
table(dat_bruss_maxPafi$maxPAFI_300plus)
dat_bruss_maxPafi[is.na(dat_bruss_maxPafi$maxPAFI_300plus),] # 03-106 does not have a PAFI either day 

ptids_pafiMaxAbove300 = dat_bruss_maxPafi$ptid[dat_bruss_maxPafi$maxPAFI_300plus & !is.na(dat_bruss_maxPafi$maxPAFI_300plus)]
ptids_pafi0Above300 = dat_bruss_maxPafi$ptid[dat_bruss_maxPafi$pafi.x >=300 & !is.na(dat_bruss_maxPafi$pafi.x)]
ptids_pafi1Above300 = dat_bruss_maxPafi$ptid[dat_bruss_maxPafi$pafi.y >= 300 & !is.na(dat_bruss_maxPafi$pafi.y)]

# Extubated day 0 or 1?
dat_term_net1 = data.frame(read_sas(paste0(dir,"/Data_2017_0411/ARMA/data/ARDSNet01/term.sas7bdat")))
dat_term_net2 = data.frame(read_sas(paste0(dir,"/Data_2017_0411/ARMA/data/ARDSNet03/term.sas7bdat")))
dat_term = merge(dat_term_net1,dat_term_net2,all=T)
colnames(dat_term) = tolower(colnames(dat_term))
head(dat_term)

# First day of successful extubation is 0 or 1: 
ptids_extubated24hrs = dat_term$ptid[dat_term$unadt %in% c(0,1) & !is.na(dat_term$unadt)]

# rrARDS if either PAFI > 300 or extubation 
ptids_rrARDS = unique(c(ptids_pafiMaxAbove300,ptids_extubated24hrs))
length(ptids_rrARDS) # 75 patients have rrARDS

# Add to dataset 
dat_mine$rrARDS_extubated24hrs = ifelse(dat_mine$ptid %in% ptids_extubated24hrs,"Extubated within 24 hrs","Not")
dat_mine$rrARDS_pafi_max_day01 = ifelse(dat_mine$ptid %in% ptids_pafiMaxAbove300,"PAFI > 300 day 0 or 1","Not")
dat_mine$rrARDS_pafi_day0 = ifelse(dat_mine$ptid %in% ptids_pafi0Above300,"PAFI > 300 day 0","Not")
dat_mine$rrARDS_pafe_day1 = ifelse(dat_mine$ptid %in% ptids_pafi1Above300,"PAFI > 300 day 1","Not ")
dat_mine$rrARDS = ifelse(dat_mine$ptid %in% ptids_rrARDS,"rrARDS","persistent ARDS")
dat_mine$rrARDS = ifelse(dat_mine$ptid %in% ptids_rrARDS,"rrARDS","persistent ARDS")
colnames(dat_mine)




dat_mine = subset(dat_screen,select=c(ptid,age,gender,ethnic,pneum,sepsis,aspir,trauma,other,multran))


# Add variables one at a time 
# Baseline 
dat_mine$age = dat_screen$age
dat_mine$gender = dat_screen$gender
dat_mine$ethnicity = dat_screen$ethnic

# BMI 
dat_mine = merge(dat_mine,subset(dat_bmi,select=c(PTID,bmi)),all.x=T,by.x="ptid",by.y="PTID")

dat_mine$ARDSrisk_pneumonia = NA
dat_mine$ARDSrisk_sepsis = NA
dat_mine$ARDSrisk_aspiration = NA
dat_mine$ARDSrisk_trauma = NA
dat_mine$ARDSrisk_otherCause = NA
dat_mine$ARDSrisk_multiTransfusion = NA
dat_mine$APACHE_III_tot = NA
dat_mine$Vasopressor = NA
dat_mine$organFail_nonPulm = NA
dat_mine$organFail_circulatory = NA
dat_mine$organFail_coag = NA
dat_mine$organFail_renal = NA
dat_mine$organFail_creat = NA
dat_mine$pfRatio_base = NA
dat_mine$ARDS_severity = NA
dat_mine$ARDS_rapidResolve = NA
dat_mine$tidalVol = NA

# Outcomes 
dat_mine$died_90days = NA
dat_mine$nonPulmorganFailFreeDays_in28 = NA
dat_mine$icu_days = NA







vars_baselineTab = c("age","gender","ethnicity","bmi","ARDSrisk_pneumonia","ARDSrisk_sepsis","ARDSrisk_aspiration","ARDSrisk_trauma","ARDSrisk_otherCause","ARDSrisk_multiTransfusion","APACHE_III_tot","Vasopressor","organFail_nonPulm","organFail_circulatory","organFail_coag","organFail_renal","organFail_creat","pfRatio_base","ARDS_cat","tidalVol")

vars_outcome = c("died_90days","nonPulmorganFailFreeDays_in28","icu_days")


# Summary Table
dat_here = dat_mine
grp = dat_here$rrARDS
table_all = rbind(
  c("Overall", paste(table(grp),paste0("[",sprintf("%1.1f",prop.table(table(grp))*100),"\\%]")),NA,length(grp)),
  combine_choices(vars=colnames(dat_here)[-1],data=dat_here,rowPct=T,NA_inPcts=F,nLevelsToBeCont=10,verbose=T)
  # summary_tab(x = dat_here$age , label="Age",type="cont",groups = grp),
  # summary_tab(x = dat_here$gender , label="Gender",type="cat",groups = grp)
  )
head(table_all)
my_printXtable(table_all, name="Baseline Covariates",N=nrow(dat_here),file=outTexReport)

head(dat_mine)

cat("\n  \\end{document} \n",file=outTexReport,append=T)


# For future use 
# Add variables one at a time 
# Baseline 
dat_mine$age = NA
dat_mine$gender = NA
dat_mine$ethnicity = NA
dat_mine$BMI = NA
dat_mine$ARDSrisk_pneumonia = NA
dat_mine$ARDSrisk_sepsis = NA
dat_mine$ARDSrisk_aspiration = NA
dat_mine$ARDSrisk_trauma = NA
dat_mine$ARDSrisk_otherCause = NA
dat_mine$ARDSrisk_multiTransfusion = NA
dat_mine$APACHE_III_tot = NA
dat_mine$Vasopressor = NA
dat_mine$organFail_nonPulm = NA
dat_mine$organFail_circulatory = NA
dat_mine$organFail_coag = NA
dat_mine$organFail_renal = NA
dat_mine$organFail_creat = NA
dat_mine$pfRatio_base = NA
dat_mine$ARDS_severity = NA
dat_mine$ARDS_rapidResolve = NA
dat_mine$tidalVol = NA

# Outcomes 
dat_mine$died_90days = NA
dat_mine$nonPulmorganFailFreeDays_in28 = NA
dat_mine$icu_days = NA

