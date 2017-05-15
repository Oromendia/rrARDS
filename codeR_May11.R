# Clara Oromendia
# March 22, 2017
# Augustine Choi, Ed Schenck, Kevin Ma
# Merge 7 trials to study ARDS: rapid resolving vs persistent


dir = "O:/Main/Research/Research_share/Clara Oromendia/Pulmonology/Rapid Resolving ARDS"
dirF = paste0(dir,"/Reports/Figures")
outTexReport = paste0(dir,"/Reports/Report.tex")
source("O:/Divisional/BIO/mco2004/myFunctions.R")



#########
# How to select row that has the highest value
#########
# oneptid %>% group_by(PTID) %>% top_n(1, PAO2)
#
#
#
#




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

##### Datasets 
dat_screen = merge( # Baseline demographics and severity of illness 
	data.frame(read_sas(paste0(dir,"/Data_2017_0411/ARMA/data/ARDSNet01/screen.sas7bdat"))), # Net 1
	data.frame(read_sas(paste0(dir,"/Data_2017_0411/ARMA/data/ARDSNet03/screen.sas7bdat"))), # Net 2
	all=T)

dat_vital_day0 = subset(merge(
	data.frame(read_sas(paste0(dir,"/Data_2017_0411/ARMA/data/ARDSNet01/vital.sas7bdat"))),
	data.frame(read_sas(paste0(dir,"/Data_2017_0411/ARMA/data/ARDSNet03/vital.sas7bdat"))),
	all=T),VIS==0)

dat_bruss = merge(
	data.frame(read_sas(paste0(dir,"/Data_2017_0411/ARMA/data/ARDSNet01/bruss.sas7bdat"))), # Net1 
	data.frame(read_sas(paste0(dir,"/Data_2017_0411/ARMA/data/ARDSNet03/bruss.sas7bdat"))), # Net 2
	all=T) # 3 ppl missing day 0 
dat_bruss$PAFI[dat_bruss$PAFI == 0 & !is.na(dat_bruss$PAFI == 0)] = NA # PAFI of 0 should be missing
dat_bruss_day0 = subset(dat_bruss, VDATE==0)
dat_bruss_day1 = subset(dat_bruss, VDATE==1)

# Use if bruss is NA 
dat_physio = merge(
	data.frame(read_sas(paste0(dir,"/Data_2017_0411/ARMA/data/ARDSNet01/physio.sas7bdat"))), # Net1 
	data.frame(read_sas(paste0(dir,"/Data_2017_0411/ARMA/data/ARDSNet03/physio.sas7bdat"))), # Net 2
	all=T) # use if bruss missing 
dat_physio$PAFI[dat_physio$PAFI == 0 & !is.na(dat_physio$PAFI == 0)] = NA # PAFI of 0 should be missing
dat_physio_day0 = subset(dat_physio, VIS==0)
dat_physio_day1 = subset(dat_physio, VIS==1)

# Random assignment
dat_random = merge(
	data.frame(read_sas(paste0(dir,"/Data_2017_0411/ARMA/data/ARDSNet01/random.sas7bdat"))), # Net1 
	data.frame(read_sas(paste0(dir,"/Data_2017_0411/ARMA/data/ARDSNet03/random.sas7bdat"))), # Net 2
	all=T)

# Blood abg 
dat_babg_day0 = merge(
	data.frame(read_sas(paste0(dir,"/Data_2017_0411/ARMA/data/ARDSNet01/babg.sas7bdat"))), # Net1 
	data.frame(read_sas(paste0(dir,"/Data_2017_0411/ARMA/data/ARDSNet03/babg.sas7bdat"))), # Net 2
	all=T) # all are day 0, each patient 20 times, 1 40 times
dat_babg_day0 = subset(dat_babg_day0,!( is.na(PAO2) & is.na(FIO2) & is.na(PACO2) & is.na(PH))) # remove if all missing 
dat_babg_day0 = subset(dat_babg_day0,!(PAO2 == 0 & FIO2 == 0 & PACO2 == 0 & PH == 0)) # remove if all 0
dat_babg_day0 = subset(dat_babg_day0,!( is.na(FIO2))) # FIO2 can't be missing at all


# Demo for history
dat_demo = merge( # Baseline demographics and severity of illness 
	data.frame(read_sas(paste0(dir,"/Data_2017_0411/ARMA/data/ARDSNet01/demo.sas7bdat"))), # Net 1
	data.frame(read_sas(paste0(dir,"/Data_2017_0411/ARMA/data/ARDSNet03/demo.sas7bdat"))), # Net 2
	all=T)

# Glascow comma
dat_glasgow = merge( 
	data.frame(read_sas(paste0(dir,"/Data_2017_0411/ARMA/data/ARDSNet01/glasgow.sas7bdat"))), # Net 1
	data.frame(read_sas(paste0(dir,"/Data_2017_0411/ARMA/data/ARDSNet03/glasgow.sas7bdat"))), # Net 2
	all=T)
dat_glasgow_day0 = subset(dat_glasgow, VIS==0)



# Contact points long form
dat_term = merge(
	data.frame(read_sas(paste0(dir,"/Data_2017_0411/ARMA/data/ARDSNet01/term.sas7bdat"))), # Net1 
	data.frame(read_sas(paste0(dir,"/Data_2017_0411/ARMA/data/ARDSNet03/term.sas7bdat"))), # Net 2
	all=T) 
dat_term$PAFI[dat_term$PAFI == 0 & !is.na(dat_term$PAFI == 0)] = NA # PAFI of 0 should be missing
dat_term_day0 = subset(dat_term, VIS==0)
dat_term_day1 = subset(dat_term, VIS==1)


#############################################
## Define rapid resolving ARDS
##############################################
# Rapid resolving if PAFI above 300 day 0 or 1, or extubated day 1 
dat_bruss_maxPafi = merge(dat_bruss_day0,dat_bruss_day1,by="PTID",all=T)
dat_bruss_maxPafi$maxPAFI = pmax(dat_bruss_maxPafi$PAFI.x, dat_bruss_maxPafi$PAFI.y,na.rm=T)
dat_bruss_maxPafi$maxPAFI_300plus = dat_bruss_maxPafi$maxPAFI >= 300

#########
# Look at PAFI (aka pf ratio) from screen, as this should have fit screening criteria
#########
range(dat_screen$PAO2/dat_screen$FIO2,na.rm=T)
dat_screen$PAFI_screen = dat_screen$PAO2/dat_screen$FIO2
mktab_singleCont(var=dat_screen$PAFI_screen, varNm="PAFI from Screen Dataset", d=2,file=outTexReport, notes="Based on screen:PAO2/screen:FIO2")
mktab_singleCont(var=dat_screen$PAO2, varNm="PAO2 from Screen Dataset", d=2,file=outTexReport, notes="")
mktab_singleCont(var=dat_screen$FIO2, varNm="FIO2 from Screen Dataset", d=2,file=outTexReport, notes="")

temp = subset(dat_screen,select=c(PTID,PAFI_screen))
temp = merge(temp,subset(dat_bruss_maxPafi,select=c(PTID,PAFI.x,PAFI.y)),by="PTID",all=T)

pdf(paste0(dirF,"/scatterPlot_PAFIsources.pdf"),width=15,height=6)
par(mfrow=c(1,3))
plot(temp$PAFI_screen,temp$PAFI.x,xlab="PAFI from Screen",ylab="PAFI from Bruss D0",xlim=c(0,550),ylim=c(0,550),main="Screen vs Bruss D0", xaxs="i",yaxs="i"); abline(0,1); abline(h=300,lty=2);abline(v=300,lty=2)
legend("topleft",title="Missing",legend=paste0(c("none:","X only:","Y only:","both:"),table(is.na(temp$PAFI_screen),is.na(temp$PAFI.x))),bty="n")
text(temp$PAFI_screen[which.max(temp$PAFI_screen)],temp$PAFI.x[which.max(temp$PAFI_screen)],temp$PTID[which.max(temp$PAFI_screen)],pos=1)
plot(temp$PAFI_screen,temp$PAFI.y,xlab="PAFI from Screen",ylab="PAFI from Bruss D1",xlim=c(0,550),ylim=c(0,550),main="Screen vs Bruss D1", xaxs="i",yaxs="i"); abline(0,1); abline(h=300,lty=2);abline(v=300,lty=2)
legend("topleft",title="Missing",legend=paste0(c("none:","X only:","Y only:","both:"),table(is.na(temp$PAFI_screen),is.na(temp$PAFI.y))),bty="n")
text(temp$PAFI_screen[which.max(temp$PAFI_screen)],temp$PAFI.x[which.max(temp$PAFI_screen)],temp$PTID[which.max(temp$PAFI_screen)],pos=1)

plot(temp$PAFI.x,temp$PAFI.y,xlab="PAFI from Bruss D0",ylab="PAFI from Bruss D1",xlim=c(0,550),ylim=c(0,550),main="Bruss D0 vs Bruss D1", xaxs="i",yaxs="i"); abline(0,1); abline(h=300,lty=2);abline(v=300,lty=2)
legend("topleft",title="Missing",legend=paste0(c("none:","X only:","Y only:","both:"),table(is.na(temp$PAFI.x),is.na(temp$PAFI.y))),bty="n")
dev.off()
figure_tex("scatterPlot_PAFIsources.pdf",file=outTexReport,caption="Compare inclusion PAFI ratio")

head(dat_bruss_maxPafi)

# Switch up or down 
dat_bruss_maxPafi$PAFI_d0d1 = NA
add_d0 = ifelse(is.na(dat_bruss_maxPafi$PAFI.x),"Missing",ifelse(dat_bruss_maxPafi$PAFI.x >= 300,"High","Low"))
add_d1 = ifelse(is.na(dat_bruss_maxPafi$PAFI.y),"Missing",ifelse(dat_bruss_maxPafi$PAFI.y >= 300,"High","Low"))
dat_bruss_maxPafi$PAFI_d0d1 = paste0(add_d0,"-",add_d1)
table(dat_bruss_maxPafi$PAFI_d0d1)
mktab_singleCate(var=dat_bruss_maxPafi$PAFI_d0d1, varNm="PAFI D0 and D1", d=2,file=outTexReport,notes="Based on Bruss Visits day 0 and 1")

# Patients who had high PAFI 
ptids_pafi0Above300 = dat_bruss_maxPafi$PTID[dat_bruss_maxPafi$PAFI.x >=300 & !is.na(dat_bruss_maxPafi$PAFI.x)]
ptids_pafi1Above300 = dat_bruss_maxPafi$PTID[dat_bruss_maxPafi$PAFI.y >= 300 & !is.na(dat_bruss_maxPafi$PAFI.y)]
ptids_pafiMaxAbove300 = unique(c(ptids_pafi0Above300,ptids_pafi1Above300))

# First day of successful extubation is 0 or 1: 
ptids_extubated24hrs = dat_term$PTID[dat_term$UNADT %in% c(0,1) & !is.na(dat_term$UNADT)]

# rrARDS if either PAFI > 300 or extubation 
ptids_rrARDS = unique(c(ptids_pafiMaxAbove300,ptids_extubated24hrs))
length(ptids_rrARDS) # 75 patients have rrARDS

# Add to dataset 
dat_mine$rrARDS = ifelse(dat_mine$ptid %in% ptids_rrARDS,"rrARDS","persistent ARDS")
dat_mine$rrARDS_type = ""
dat_mine$rrARDS_type[dat_mine$ptid %in% ptids_extubated24hrs] = paste(dat_mine$rrARDS_type[dat_mine$ptid %in% ptids_extubated24hrs], "Extubated wIn 24h",sep=" & ")
dat_mine$rrARDS_type[dat_mine$ptid %in% ptids_pafi0Above300] = paste(dat_mine$rrARDS_type[dat_mine$ptid %in% ptids_pafi0Above300], "PAFI_d0 high",sep=" & ")
dat_mine$rrARDS_type[dat_mine$ptid %in% ptids_pafi1Above300] = paste(dat_mine$rrARDS_type[dat_mine$ptid %in% ptids_pafi1Above300], "PAFI_d1 high",sep=" & ")
dat_mine$rrARDS_type[dat_mine$rrARDS_type == ""] = "Persistent ARDS"
dat_mine$rrARDS_type = trimws(dat_mine$rrARDS_type); dat_mine$rrARDS_type = gsub("^&","",dat_mine$rrARDS_type)

table(dat_mine$rrARDS_type)
mktab_singleCate(var=dat_mine$rrARDS_type, varNm="ARDS designation", d=2,file=outTexReport,notes="Based on Bruss Visits day 0 and 1 and term extubation day")





# APACHE III scoring 
##############
dat_physio_day0$APACHE_tot = NA

# Pulse 
scoreAPACHE_pulse <- function(hrHigh,hrLow,ifMissingMake0 = TRUE) {
	scoreAPACHE_pulse_one <- function(hr){
		score = NA
		score[!is.na(hr) & hr <= 39] = 8
		score[!is.na(hr) & hr >= 40 & hr <= 49] = 5
		score[!is.na(hr) & hr >= 50 & hr <= 99] = 0
		score[!is.na(hr) & hr >= 100 & hr <= 109] = 1
		score[!is.na(hr) & hr >= 110 & hr <= 119] = 5
		score[!is.na(hr) & hr >= 120 & hr <= 139] = 7
		score[!is.na(hr) & hr >= 140 & hr <= 154] = 13
		score[!is.na(hr) & hr >= 155] = 17
		score[is.na(hr)] = 0
		return(score)
	}
	score = pmax(scoreAPACHE_pulse_one(hrHigh),scoreAPACHE_pulse_one(hrLow))
	return(score)
}

# Pulse is worse of high and low heart rate 
dat_physio_day0$APACHE_pulse = scoreAPACHE_pulse(hrHigh = dat_physio_day0$HRATEH,hrLow = dat_physio_day0$HRATEL)
table(dat_physio_day0$APACHE_pulse)

# Blood Pressure
scoreAPACHE_bp <- function(bpHigh,bpLow,ifMissingMake0 = TRUE) {
	scoreAPACHE_bp_one <- function(bp){
		score = NA
		score[!is.na(bp) & bp <= 39] = 23
		score[!is.na(bp) & bp >= 40 & bp <= 59] = 15
		score[!is.na(bp) & bp >= 60 & bp <= 69] = 7
		score[!is.na(bp) & bp >= 70 & bp <= 79] = 6
		score[!is.na(bp) & bp >= 80 & bp <= 99] = 0
		score[!is.na(bp) & bp >= 100 & bp <= 119] = 4
		score[!is.na(bp) & bp >= 120 & bp <= 129] = 7
		score[!is.na(bp) & bp >= 130 & bp <= 139] = 9
		score[!is.na(bp) & bp >= 140] = 10
		score[is.na(bp)] = 0
		return(score)
	}
	score = pmax(scoreAPACHE_bp_one(bpHigh),scoreAPACHE_bp_one(bpLow))
	return(score)
}

# Bp is worse of high and low mean arterial pressure 
dat_physio_day0$APACHE_bp = scoreAPACHE_bp(bpHigh = dat_physio_day0$MEANAPL,bpLow = dat_physio_day0$MEANAPH)
table(dat_physio_day0$APACHE_bp)

# Temperature 
scoreAPACHE_temp <- function(tempHigh,tempLow) {
	scoreAPACHE_temp_one <- function(temp){
		score = NA
		score[!is.na(temp) & temp <= 32.9] = 20
		score[!is.na(temp) & temp >= 33.0 & temp <= 33.4] = 16
		score[!is.na(temp) & temp >= 33.5 & temp <= 33.9] = 13
		score[!is.na(temp) & temp >= 34.0 & temp <= 34.9] = 8
		score[!is.na(temp) & temp >= 35.0 & temp <= 35.9] = 2
		score[!is.na(temp) & temp >= 36.0 & temp <= 39.9] = 0
		score[!is.na(temp) & temp >= 40] = 4
		score[is.na(temp)] = 0
		return(score)
	}
	score = pmax(scoreAPACHE_temp_one(tempHigh),scoreAPACHE_temp_one(tempLow))
	return(score)
}

# Temp is worse of highest and lowest temp 
dat_physio_day0$APACHE_temp = scoreAPACHE_temp(tempHigh = dat_physio_day0$TEMPCL,tempLow = dat_physio_day0$TEMPCH)
table(dat_physio_day0$APACHE_temp)

# Respiratory rate  
scoreAPACHE_resp <- function(respHigh,respLow,onVentilatorwhenLow) {
	scoreAPACHE_resp_one <- function(resp){
		score = NA
		score[!is.na(resp) & resp <= 5] = 17
		score[!is.na(resp) & resp >= 6 & resp <= 11] = 8
		score[!is.na(resp) & resp >= 12 & resp <= 13] = 7
		score[!is.na(resp) & resp >= 14 & resp <= 24] = 0
		score[!is.na(resp) & resp >= 25 & resp <= 34] = 6
		score[!is.na(resp) & resp >= 35 & resp <= 39] = 9
		score[!is.na(resp) & resp >= 40 & resp <= 49] = 11
		score[!is.na(resp) & resp >= 50] = 18
		score[is.na(resp)] = 0
		return(score)
	}
	score_high = scoreAPACHE_resp_one(respHigh)
	score_low = scoreAPACHE_resp_one(respLow)
	score_low[onVentilatorwhenLow == TRUE & !is.na(onVentilatorwhenLow) & respLow >= 6 & respLow <= 13] = 0
	score = pmax(score_high,score_low)
	return(score)
}

dat_physio_day0$APACHE_resp = scoreAPACHE_resp(respHigh = dat_physio_day0$RESPH, respLow = dat_physio_day0$RESPL, onVentilatorwhenLow = dat_physio_day0$LVENT == 1)




# Find babg measurements with worst ABG (arterial blood gas) measurement for that day 

# Find worst abg: emailed Ed, but for now use highest PAO2 
# Use measurement that gives highest PA02_or_AaDO2 score

# PA02_or_AaDO2  
scoreAPACHE_PA02_or_AaDO2 <- function(PAO2,FIO2,PACO2) {
	if (any(is.na(FIO2))) stop("Remove all missing FIO2 entries")
	score_FIO2_50less = NA
	score_FIO2_50less[!is.na(PAO2) & PAO2 <= 49] = 15
	score_FIO2_50less[!is.na(PAO2) & PAO2 >= 50 & PAO2 <= 69] = 5
	score_FIO2_50less[!is.na(PAO2) & PAO2 >= 70 & PAO2 <= 79] = 2
	score_FIO2_50less[!is.na(PAO2) & PAO2 >= 80] = 0
	
	# If FIO2 is more than 50%: 
	AaDo2 = 713*FIO2 - PACO2/0.8 - PAO2
	score_FIO2_50plus = NA
	score_FIO2_50plus[!is.na(AaDo2) & AaDo2 < 100] = 0
	score_FIO2_50plus[!is.na(AaDo2) & AaDo2 >= 100 & AaDo2 < 250] = 7
	score_FIO2_50plus[!is.na(AaDo2) & AaDo2 >= 250 & AaDo2 < 350] = 9
	score_FIO2_50plus[!is.na(AaDo2) & AaDo2 >= 350 & AaDo2 < 500] = 11
	score_FIO2_50plus[!is.na(AaDo2) & AaDo2 >= 500] = 14
	
	score = NA
	score[FIO2 < 0.5] = score_FIO2_50less[FIO2 < 0.5]
	score[FIO2 >= 0.5] = score_FIO2_50plus[FIO2 >= 0.5]
	
	return(score)
}


dat_babg_day0$scoreAPACHE_PAO2orAaDO2 = scoreAPACHE_PA02_or_AaDO2(PAO2 = dat_babg_day0$PAO2,FIO2 = dat_babg_day0$FIO2,PACO2 = dat_babg_day0$PACO2)
# For now calculate every time it was measured, eventually will choose the measurement associated with the highest APACHe score for this part. 


# Hematocrit
scoreAPACHE_hemat <- function(hematHigh,hematLow) {
	scoreAPACHE_hemat_one <- function(hemat){
		score = NA
		score[!is.na(hemat) & hemat <= 40.9] = 3
		score[!is.na(hemat) & hemat >= 41 & hemat <= 49] = 0
		score[!is.na(hemat) & hemat >= 50] = 3
		score[is.na(hemat)] = 0
		return(score)
	}
	score_high = scoreAPACHE_hemat_one(hematHigh)
	score_low = scoreAPACHE_hemat_one(hematLow)
	score = pmax(score_high,score_low)
	return(score)
}

dat_physio_day0$APACHE_hemat = scoreAPACHE_hemat(hematHigh = dat_physio_day0$HCTH, hematLow = dat_physio_day0$HCTL)




# Creatinine
scoreAPACHE_creat <- function(creatHigh,URINE,dialysisChronic) {
	
	hasARF = creatHigh >= 1.5 & URINE < 410 & dialysisChronic == FALSE
	# If they don't have acute renal failure 
	score_creat_ARF = NA
	score_creat_ARF[!is.na(creatHigh) & creatHigh <= 49] = 15
	score_creat_ARF[!is.na(creatHigh) & creatHigh >= 50 & creat <= 69] = 5
	score_creat_ARF[!is.na(creatHigh) & creatHigh >= 70 & creat <= 79] = 2
	score_creat_ARF[!is.na(creatHigh) & creatHigh >= 80] = 0
	
	# If they have acute renal failure 
	score_creat_noARF = NA
	score_creat_noARF[!is.na(creatHigh) & creatHigh < 100] = 0
	score_creat_noARF[!is.na(creatHigh) & creatHigh >= 100 & creat < 250] = 7
	score_creat_noARF[!is.na(creatHigh) & creatHigh >= 250 & creat < 350] = 9
	score_creat_noARF[!is.na(creatHigh) & creatHigh >= 350 & creat < 500] = 11
	score_creat_noARF[!is.na(creatHigh) & creatHigh >= 500] = 14
	
	score = NA
	score[hasARF & !is.na(hasARF)] = score_creat_ARF[hasARF & !is.na(hasARF)]
	score[!hasARF & !is.na(hasARF)] = score_creat_noARF[!hasARF & !is.na(hasARF)]
	
	return(score)
}

# Merge chronic dialysis to ensure ptids match 
dat_physio_day0 = merge(dat_physio_day0, subset(dat_demo,select=c(PTID,DIALY)),all.x=T) 
scoreAPACHE_creat(creatHigh = dat_physio_day0$CREATH,URINE = dat_physio_day0$URINE,dialysisChronic = dat_physio_day0$DIALY == 1)
table(dat_physio_day0$DIALY)
# Note: if missing ARF designation, missing score (not 0)




# Urine output 
scoreAPACHE_urine <- function(urine) {
	score = NA
	score[!is.na(urine) & urine <= 399] = 15
	score[!is.na(urine) & urine >= 400 & urine <= 599] = 8
	score[!is.na(urine) & urine >= 600 & urine <= 899] = 7
	score[!is.na(urine) & urine >= 900 & urine <= 1499] = 5
	score[!is.na(urine) & urine >= 1500 & urine <= 1999] = 4
	score[!is.na(urine) & urine >= 2000 & urine <= 3999] = 0
	score[!is.na(urine) & urine >= 4000] = 1
	score[is.na(urine)] = 0; print(paste(sum(is.na(urine)),"missing made 0"))
	return(score)
}

# Single urine output measurement
dat_physio_day0$APACHE_urine = scoreAPACHE_urine(urine = dat_physio_day0$URINE)
table(dat_physio_day0$APACHE_urine)


# BUN 
scoreAPACHE_bun <- function(bun) {
	score = NA
	score[!is.na(bun) & bun <= 16.9] = 0
	score[!is.na(bun) & bun >= 17 & bun <= 19] = 2
	score[!is.na(bun) & bun >= 20 & bun <= 39] = 7
	score[!is.na(bun) & bun >= 40 & bun <= 79] = 11
	score[!is.na(bun) & bun >= 80] = 12
	score[is.na(bun)] = 0; print(paste(sum(is.na(bun)),"missing made 0"))
	return(score)
}

# Single BUN measurement 
dat_physio_day0$APACHE_bun = scoreAPACHE_bun(bun = dat_physio_day0$BUN)
table(dat_physio_day0$APACHE_bun)



# Sodium
scoreAPACHE_sodium <- function(sodiumHigh,sodiumLow) {
	scoreAPACHE_sodium_one <- function(sodium){
		score = NA
		score[!is.na(sodium) & sodium <= 119] = 3
		score[!is.na(sodium) & sodium >= 120 & sodium <= 134] = 2
		score[!is.na(sodium) & sodium >= 135 & sodium <= 154] = 0
		score[!is.na(sodium) & sodium >= 155] = 4
		score[is.na(sodium)] = 0;  print(paste(sum(is.na(sodium)),"missing made 0"))
		return(score)
	}
	score_high = scoreAPACHE_sodium_one(sodiumHigh)
	score_low = scoreAPACHE_sodium_one(sodiumLow)
	score = pmax(score_high,score_low)
	return(score)
}

dat_physio_day0$APACHE_sodium = scoreAPACHE_sodium(sodiumHigh = dat_physio_day0$SODIUMH, sodiumLow = dat_physio_day0$SODIUML)
table(dat_physio_day0$APACHE_sodium)

# Albumin
scoreAPACHE_albumin <- function(albuminHigh,albuminLow) {
	scoreAPACHE_albumin_one <- function(albumin){
		score = NA
		score[!is.na(albumin) & albumin <= 1.9] = 11
		score[!is.na(albumin) & albumin >= 2.0 & albumin <= 2.4] = 6
		score[!is.na(albumin) & albumin >= 2.5 & albumin <= 4.4] = 0
		score[!is.na(albumin) & albumin >= 4.5] = 4
		score[is.na(albumin)] = 0;  print(paste(sum(is.na(albumin)),"missing made 0"))
		return(score)
	}
	score_high = scoreAPACHE_albumin_one(albuminHigh)
	score_low = scoreAPACHE_albumin_one(albuminLow)
	score = pmax(score_high,score_low)
	return(score)
}

dat_physio_day0$APACHE_albumin = scoreAPACHE_albumin(albuminHigh = dat_physio_day0$ALBUMH, albuminLow = dat_physio_day0$ALBUML)
table(dat_physio_day0$APACHE_albumin)




# Bilirubin
scoreAPACHE_bili <- function(bili) {
	score = NA
	score[!is.na(bili) & bili <= 1.9] = 0
	score[!is.na(bili) & bili >= 2.0 & bili <= 2.9] = 5
	score[!is.na(bili) & bili >= 3.0 & bili <= 4.9] = 6
	score[!is.na(bili) & bili >= 5.0 & bili <= 7.9] = 8
	score[!is.na(bili) & bili >= 8] = 16
	score[is.na(bili)] = 0; print(paste(sum(is.na(bili)),"missing made 0"))
	return(score)
}

# Single BUN measurement 
dat_physio_day0$APACHE_bili = scoreAPACHE_bili(bili = dat_physio_day0$BILI)
table(dat_physio_day0$APACHE_bili)


# Glucose
scoreAPACHE_glucose <- function(glucoseHigh,glucoseLow) {
	scoreAPACHE_glucose_one <- function(glucose){
		score = NA
		score[!is.na(glucose) & glucose <= 39] = 8
		score[!is.na(glucose) & glucose >= 40 & glucose <= 59] = 9
		score[!is.na(glucose) & glucose >= 60 & glucose <= 199] = 0
		score[!is.na(glucose) & glucose >= 200 & glucose <= 349] = 3
		score[!is.na(glucose) & glucose >= 350] = 5
		score[is.na(glucose)] = 0;  print(paste(sum(is.na(glucose)),"missing made 0"))
		return(score)
	}
	score_high = scoreAPACHE_glucose_one(glucoseHigh)
	score_low = scoreAPACHE_glucose_one(glucoseLow)
	score = pmax(score_high,score_low)
	return(score)
}

dat_physio_day0$APACHE_glucose = scoreAPACHE_glucose(glucoseHigh = dat_physio_day0$GLUCH, glucoseLow = dat_physio_day0$GLUCL)
table(dat_physio_day0$APACHE_glucose)


# Age
scoreAPACHE_age <- function(age) {
	score = NA
	score[!is.na(age) & age <= 44] = 0
	score[!is.na(age) & age >= 45 & age <= 59] = 5
	score[!is.na(age) & age >= 60 & age <= 64] = 11
	score[!is.na(age) & age >= 65 & age <= 69] = 13
	score[!is.na(age) & age >= 70 & age <= 74] = 16
	score[!is.na(age) & age >= 75 & age <= 84] = 17
	score[!is.na(age) & age >= 85] = 24
	score[is.na(age)] = 0; print(paste(sum(is.na(age)),"missing made 0"))
	return(score)
}

# Single Age measurement 
dat_physio_day0 = merge(dat_physio_day0,subset(dat_screen,select=c(PTID,AGE)),by="PTID",all.x=T)
dat_physio_day0$APACHE_age = scoreAPACHE_age(age = dat_physio_day0$AGE)
table(dat_physio_day0$APACHE_age)


# Chronic health comorbidities 
scoreAPACHE_chronicHealth <- function(electiveSurg, aids, leukemia, lymphoma, tumor, immune, hepa, cirrhosis) {
	score = NA
	score[is.na(score) & !is.na(aids) & aids == TRUE] = 23
	score[is.na(score) & !is.na(hepa) & hepa == TRUE] = 16
	score[is.na(score) & !is.na(lymphoma) & lymphoma == TRUE] = 13
	score[is.na(score) & !is.na(tumor) & tumor == TRUE] = 11
	score[is.na(score) & !is.na(leukemia) & leukemia == TRUE] = 10
	score[is.na(score) & !is.na(immune) & immune == TRUE] = 10
	score[is.na(score) & !is.na(cirrhosis) & cirrhosis == TRUE] = 4

	# If elective surgery, score is 0 
	score[electiveSurg == T & !is.na(electiveSurg)] = 0
	print(paste(sum(is.na(score)),"missing made 0"))
	score[is.na(score)] = 0
	return(score)
}

# Chronic Health
dat_physio_day0 = merge(dat_physio_day0,subset(dat_demo,select=c(PTID,SURGEL,AIDS,LEUK,LYMPH,TUMOR,IMMUNE,HEPA,CIRR)),by="PTID",all.x=T)
dat_physio_day0$APACHE_chronicH = scoreAPACHE_chronicHealth(
	electiveSurg = dat_physio_day0$SURGEL == 1, 
	aids = dat_physio_day0$AIDS == 1, 
	leukemia = dat_physio_day0$LEUK == 1, 
	lymphoma = dat_physio_day0$LYMPH == 1, 
	tumor = dat_physio_day0$TUMOR == 1, 
	immune = dat_physio_day0$IMMUNE == 1, 
	hepa = dat_physio_day0$HEPA == 1, 
	cirrhosis = dat_physio_day0$CIRR == 1)
table(dat_physio_day0$APACHE_chronicH)



# Neurologic abnormalities 
scoreAPACHE_gcs <- function(eye,verbal,motor) {
	score_yesEye = NA
	score_yesEye[verbal == "Oriented, converses" & motor == "Obeys verbal command"] = 0
	score_yesEye[verbal == "Oriented, converses" & motor == "Localizes pain"] = 3
	score_yesEye[verbal == "Oriented, converses" & motor == "Flexion withdrawl/decorticate rigidity"] = 3
	score_yesEye[verbal == "Oriented, converses" & motor == "Decerebrate rigidity/no response"] = 3
	score_yesEye[verbal == "Confused conversation" & motor == "Obeys verbal command"] = 3
	score_yesEye[verbal == "Confused conversation" & motor == "Localizes pain"] = 8
	score_yesEye[verbal == "Confused conversation" & motor == "Flexion withdrawl/decorticate rigidity"] = 13
	score_yesEye[verbal == "Confused conversation" & motor == "Decerebrate rigidity/no response"] = 13
	score_yesEye[verbal == "Inappropriate Words & Incomprehensible Sounds"& motor == "Obeys verbal command"] = 10
	score_yesEye[verbal == "Inappropriate Words & Incomprehensible Sounds"& motor == "Localizes pain"] = 13
	score_yesEye[verbal == "Inappropriate Words & Incomprehensible Sounds"& motor == "Flexion withdrawl/decorticate rigidity"]= 24
	score_yesEye[verbal == "Inappropriate Words & Incomprehensible Sounds"& motor == "Decerebrate rigidity/no response"] = 29
	score_yesEye[verbal == "No response" & motor == "Obeys verbal command"] = 15
	score_yesEye[verbal == "No response" & motor == "Localizes pain"] = 15
	score_yesEye[verbal == "No response" & motor == "Flexion withdrawl/decorticate rigidity"] = 24
	score_yesEye[verbal == "No response" & motor == "Decerebrate rigidity/no response"] = 29

	score_noEye = NA
	score_noEye[verbal %in%  c("Oriented, converses","Confused conversation")] = 99999 # Should not happen
	score_noEye[verbal == "Inappropriate Words & Incomprehensible Sounds"& motor == "Obeys verbal command"] = 99999
	score_noEye[verbal == "Inappropriate Words & Incomprehensible Sounds"& motor == "Localizes pain"] = 99999
	score_noEye[verbal == "Inappropriate Words & Incomprehensible Sounds"& motor == "Flexion withdrawl/decorticate rigidity"] = 24
	score_noEye[verbal == "Inappropriate Words & Incomprehensible Sounds"& motor == "Decerebrate rigidity/no response"] = 29
	score_noEye[verbal == "No response" & motor == "Obeys verbal command"] = 16
	score_noEye[verbal == "No response" & motor == "Localizes pain"] = 16
	score_noEye[verbal == "No response" & motor == "Flexion withdrawl/decorticate rigidity"] = 33
	score_noEye[verbal == "No response" & motor == "Decerebrate rigidity/no response"] = 48

	score = ifelse(eye,score_yesEye,score_noEye)
table(score)
table(score_yesEye)
table(score_noEye)
	print(paste(sum(is.na(score)),"missing made 0"))
	score[is.na(score)] = 0; 
cat(dat_physio_day0$PTID[which(score == 99999)],sep="\n")

table(dat_physio_day0$MOTOR[which(score == 99999)],dat_physio_day0$VERBAL[which(score == 99999)])
table(dat_physio_day0$MOTOR_cat[which(score == 99999)],dat_physio_day0$VERBAL_cat[which(score == 99999)])
table(dat_physio_day0$MOTOR,dat_physio_day0$VERBAL,dat_physio_day0$EYE)
###
###
###
###
# Emailed Ed 5/12: Issue that 17 have values that don't make sense clinically 
###
###
###
###
table(motor,eye)
	return(score)
}

# Single GCS measurement 
dat_physio_day0 = merge(dat_physio_day0,subset(dat_glasgow_day0,select=c(PTID,EYE,VERBAL,MOTOR)),by="PTID",all.x=T)
dat_physio_day0$VERBAL_cat = NA
dat_physio_day0$VERBAL_cat[dat_physio_day0$VERBAL == 5] = "Oriented, converses"
dat_physio_day0$VERBAL_cat[dat_physio_day0$VERBAL == 4] = "Confused conversation"
dat_physio_day0$VERBAL_cat[dat_physio_day0$VERBAL %in% c(3,2)] = "Inappropriate Words & Incomprehensible Sounds"
dat_physio_day0$VERBAL_cat[dat_physio_day0$VERBAL == 1] = "No response"

dat_physio_day0$MOTOR_cat = NA
dat_physio_day0$MOTOR_cat[dat_physio_day0$MOTOR == 6] = "Obeys verbal command"
dat_physio_day0$MOTOR_cat[dat_physio_day0$MOTOR == 5] = "Localizes pain"
dat_physio_day0$MOTOR_cat[dat_physio_day0$MOTOR %in% c(4,3)] = "Flexion withdrawl/decorticate rigidity"
dat_physio_day0$MOTOR_cat[dat_physio_day0$MOTOR %in% c(2,1)] = "Decerebrate rigidity/no response"

dat_physio_day0$APACHE_gcs = scoreAPACHE_gcs(eye = dat_physio_day0$EYE >= 2, verbal = dat_physio_day0$VERBAL_cat, motor = dat_physio_day0$MOTOR_cat)
table(dat_physio_day0$APACHE_gcs) # 9999 are not clinically credible



# Acid-Base: eventially from same abg are oxygenation  
scoreAPACHE_acidBase <- function(pH, PACO2) {
	score = NA
	score[!is.na(pH) & !is.na(PAO2) & pH < 7.20 & PACO2 < 50] = 12
	score[!is.na(pH) & !is.na(PAO2) & pH < 7.20 & PACO2 >= 50] = 4
	score[!is.na(pH) & !is.na(PAO2) & pH >= 7.20 & pH < 7.35 & PACO2 < 30] = 9
	score[!is.na(pH) & !is.na(PAO2) & pH >= 7.20 & pH < 7.30 & PACO2 >= 30 & PACO2 < 35] = 6
	score[!is.na(pH) & !is.na(PAO2) & pH >= 7.20 & pH < 7.30 & PACO2 >= 35 & PACO2 < 50] = 3
	score[!is.na(pH) & !is.na(PAO2) & pH >= 7.20 & pH < 7.30 & PACO2 >= 50] = 2
	score[!is.na(pH) & !is.na(PAO2) & pH >= 7.35 & pH < 7.50 & PACO2 < 30] = 5
	score[!is.na(pH) & !is.na(PAO2) & pH >= 7.30 & pH < 7.45 & PACO2 >= 30 & PACO2 < 45] = 0
	score[!is.na(pH) & !is.na(PAO2) & pH >= 7.30 & pH < 7.45 & PACO2 >= 45] = 1
	score[!is.na(pH) & !is.na(PAO2) & pH >= 7.45 & pH < 7.50 & PACO2 >= 30 & PACO2 < 35] = 0
	score[!is.na(pH) & !is.na(PAO2) & pH >= 7.45 & pH < 7.50 & PACO2 >= 35 & PACO2 < 45] = 2
	score[!is.na(pH) & !is.na(PAO2) & pH >= 7.50 & pH < 7.65 & PACO2 < 25] = 3
	score[!is.na(pH) & !is.na(PAO2) & pH >= 7.50 & PACO2 >= 25 & PACO2 < 35] = 3
	score[!is.na(pH) & !is.na(PAO2) & pH >= 7.65 & PACO2 < 25] = 0
	score[!is.na(pH) & !is.na(PAO2) & pH >= 7.50 & PACO2 >= 35] = 12
	score[!is.na(pH) & !is.na(PAO2) & pH >= 7.45 & pH < 7.50 & PACO2 >= 45] = 12
	print(paste(sum(is.na(score)),"missing made 0"))
	score[is.na(score)] = 0

	return(score)
}

dat_babg_day0$scoreAPACHE_acidBase = scoreAPACHE_acidBase(pH = dat_babg_day0$PH,PACO2 = dat_babg_day0$PACO2)
table(dat_babg_day0$scoreAPACHE_acidBase)

table(dat_babg_day0$PH == 0)
dat_babg_day0[dat_babg_day0$PH == 0,]

# Declare my dataset for the study
dat_ARMA = data.frame(ptid = dat_screen$PTID)


###################
###################
# Add variables one at a time 
###################
###################

# Baseline
###################
# Age 
dat_screen$age = dat_screen$AGE # format here if necessary 
dat_ARMA = merge(dat_ARMA,subset(dat_screen,select=c(PTID,age)),all.x=T,by.x="ptid",by.y="PTID")

# Gender
dat_screen$gender = dat_screen$GENDER
dat_ARMA = merge(dat_ARMA,subset(dat_screen,select=c(PTID,gender)),all.x=T,by.x="ptid",by.y="PTID")

# Ethnicity 
dat_screen$ethnicity = dat_screen$ETHNIC
dat_ARMA = merge(dat_ARMA,subset(dat_screen,select=c(PTID,ethnicity)),all.x=T,by.x="ptid",by.y="PTID")

# BMI 
	dat_vital_day0$bmi = dat_vital_day0$WEIGHTK/dat_vital_day0$HEIGHTC
	dat_ARMA = merge(dat_ARMA,subset(dat_vital_day0,select=c(PTID,bmi)),all.x=T,by.x="ptid",by.y="PTID")

# ARDS Risk factors (primary only)
	dat_screen$ARDSrisk_pneumonia = (dat_screen$PNEUM == 1)
	dat_ARMA = merge(dat_ARMA,subset(dat_screen,select=c(PTID,ARDSrisk_pneumonia)),all.x=T,by.x="ptid",by.y="PTID")

	dat_screen$ARDSrisk_sepsis = (dat_screen$SEPSIS == 1)
	dat_ARMA = merge(dat_ARMA,subset(dat_screen,select=c(PTID,ARDSrisk_sepsis)),all.x=T,by.x="ptid",by.y="PTID")

	dat_screen$ARDSrisk_aspiration = (dat_screen$ASPIR == 1)
	dat_ARMA = merge(dat_ARMA,subset(dat_screen,select=c(PTID,ARDSrisk_aspiration)),all.x=T,by.x="ptid",by.y="PTID")

	dat_screen$ARDSrisk_trauma = (dat_screen$TRAUMA == 1)
	dat_ARMA = merge(dat_ARMA,subset(dat_screen,select=c(PTID,ARDSrisk_trauma)),all.x=T,by.x="ptid",by.y="PTID")

	dat_screen$ARDSrisk_otherCause = (dat_screen$OTHER == 1)
	dat_ARMA = merge(dat_ARMA,subset(dat_screen,select=c(PTID,ARDSrisk_otherCause)),all.x=T,by.x="ptid",by.y="PTID")

	dat_screen$ARDSrisk_multiTransfusion = (dat_screen$MULTRAN == 1)
	dat_ARMA = merge(dat_ARMA,subset(dat_screen,select=c(PTID,ARDSrisk_multiTransfusion)),all.x=T,by.x="ptid",by.y="PTID")

# APACHE 
dat_ARMA$APACHE_III_tot = NA

# Vasopressors 
dat_bruss_day0$vasopressor = ifelse(dat_bruss_day0$VASO == 1,"Yes",ifelse(dat_bruss_day0$VASO == 2,"No",NA))
dat_ARMA = merge(dat_ARMA,subset(dat_bruss_day0,select=c(PTID,vasopressor)),all.x=T,by.x="ptid",by.y="PTID")


# Organ failure at base 
 # If missing, use physio at visit 0

# Circulatory
dat_bruss_day0$organFail_circulatory = (dat_bruss_day0$SYSBP <= 90 | dat_bruss_day0$VASO == 1)
dat_ARMA = merge(dat_ARMA,subset(dat_bruss_day0,select=c(PTID,organFail_circulatory)),all.x=T,by.x="ptid",by.y="PTID")

# Coagulation 
dat_bruss_day0$organFail_coag = (dat_bruss_day0$PLATE < 80 )
dat_physio_day0$organFail_coag_physio = (dat_physio_day0$PLATE < 80 )
dat_bruss_day0 = merge(dat_bruss_day0,subset(dat_physio_day0,select=c(PTID,organFail_coag_physio)),all=T,by="PTID")
dat_bruss_day0$organFail_coag[is.na(dat_bruss_day0$organFail_coag)] = dat_bruss_day0$organFail_coag_physio[is.na(dat_bruss_day0$organFail_coag)]
dat_ARMA = merge(dat_ARMA,subset(dat_bruss_day0,select=c(PTID,organFail_coag)),all.x=T,by.x="ptid",by.y="PTID")

# Hepatic
dat_bruss_day0$organFail_hepatic = (dat_bruss_day0$BILI >= 2 )
dat_physio_day0$organFail_hepatic_physio = (dat_physio_day0$BILI >= 2 )
dat_bruss_day0 = merge(dat_bruss_day0,subset(dat_physio_day0,select=c(PTID,organFail_hepatic_physio)),all=T,by="PTID")
dat_bruss_day0$organFail_hepatic[is.na(dat_bruss_day0$organFail_hepatic)] = dat_bruss_day0$organFail_hepatic_physio[is.na(dat_bruss_day0$organFail_hepatic)]
dat_ARMA = merge(dat_ARMA,subset(dat_bruss_day0,select=c(PTID,organFail_hepatic)),all.x=T,by.x="ptid",by.y="PTID")

# Renal
dat_bruss_day0$organFail_renal = (dat_bruss_day0$CREAT >= 2 )
dat_ARMA = merge(dat_ARMA,subset(dat_bruss_day0,select=c(PTID,organFail_renal)),all.x=T,by.x="ptid",by.y="PTID")

dat_bruss_day0$organFail_renal = (dat_bruss_day0$CREAT >= 2 )
dat_physio_day0$organFail_renal_physio = (dat_physio_day0$CREATH >= 2 )
dat_bruss_day0 = merge(dat_bruss_day0,subset(dat_physio_day0,select=c(PTID,organFail_renal_physio)),all=T,by="PTID")
dat_bruss_day0$organFail_renal[is.na(dat_bruss_day0$organFail_renal)] = dat_bruss_day0$organFail_renal_physio[is.na(dat_bruss_day0$organFail_renal)]
dat_ARMA = merge(dat_ARMA,subset(dat_bruss_day0,select=c(PTID,organFail_renal)),all.x=T,by.x="ptid",by.y="PTID")


# Count of failures at day 0 
dat_bruss_day0$organFail_numNonPulm = rowSums(cbind(dat_bruss_day0$organFail_circulatory,dat_bruss_day0$organFail_coag,dat_bruss_day0$organFail_hepatic,dat_bruss_day0$organFail_renal),na.rm=T)
dat_ARMA = merge(dat_ARMA,subset(dat_bruss_day0,select=c(PTID,organFail_numNonPulm)),all.x=T,by.x="ptid",by.y="PTID")

# PF ratio
dat_bruss_day0$pfRatio_day0 = dat_bruss_day0$PAFI
dat_ARMA = merge(dat_ARMA,subset(dat_bruss_day0,select=c(PTID,pfRatio_day0)),all.x=T,by.x="ptid",by.y="PTID")

# ARDS severity (based on PF ratio)
dat_ARMA$ARDS_severity_day0 = NA
dat_ARMA$ARDS_severity_day0[dat_ARMA$pfRatio_day0 > 300 & !is.na(dat_ARMA$pfRatio_day0)] = "Day 0 PF Ratio above 300"
dat_ARMA$ARDS_severity_day0[dat_ARMA$pfRatio_day0 >200 & dat_ARMA$pfRatio_day0 <= 300 & !is.na(dat_ARMA$pfRatio_day0)] = "Mild ARDS (200-300]"
dat_ARMA$ARDS_severity_day0[dat_ARMA$pfRatio_day0 >100 & dat_ARMA$pfRatio_day0 <= 200 & !is.na(dat_ARMA$pfRatio_day0)] = "Moderate ARDS (100-200]"
dat_ARMA$ARDS_severity_day0[dat_ARMA$pfRatio_day0 <= 100 & !is.na(dat_ARMA$pfRatio_day0)] = "Severe ARDS (0-100]"


# Treatment assignment and inclusion 
dat_random$tidalVol = ifelse(dat_random$RX == "Randomized: 12 ml/kg","high tidal vol","low tidal vol")
dat_random$include_baseline = "Include" # Include all patients for baseline rrARDS 
dat_random$include_treatmentEffect = "Include"
dat_random$include_treatmentEffect[grep("Assigned",dat_random$RX)] = "Exclude" # some assigned but trial stopped early
dat_ARMA = merge(dat_ARMA,subset(dat_random,select=c(PTID,tidalVol,include_baseline,include_treatmentEffect)),all.x=T,by.x="ptid",by.y="PTID")




dat_ARMA$ARDS_rapidResolve = NA

# Outcomes 
dat_ARMA$died_90days = NA
dat_ARMA$nonPulmorganFailFreeDays_in28 = NA
dat_ARMA$icu_days = NA







vars_baselineTab = c("age","gender","ethnicity","bmi","ARDSrisk_pneumonia","ARDSrisk_sepsis","ARDSrisk_aspiration","ARDSrisk_trauma","ARDSrisk_otherCause","ARDSrisk_multiTransfusion","APACHE_III_tot","Vasopressor","organFail_nonPulm","organFail_circulatory","organFail_coag","organFail_renal","organFail_creat","pfRatio_base","ARDS_cat","tidalVol")

vars_outcome = c("died_90days","nonPulmorganFailFreeDays_in28","icu_days")


# Summary Table
dat_here = dat_ARMA
grp = dat_here$rrARDS
table_all = rbind(
  c("Overall", paste(table(grp),paste0("[",sprintf("%1.1f",prop.table(table(grp))*100),"\\%]")),NA,length(grp)),
  combine_choices(vars=colnames(dat_here)[-1],data=dat_here,rowPct=T,NA_inPcts=F,nLevelsToBeCont=10,verbose=T)
  # summary_tab(x = dat_here$age , label="Age",type="cont",groups = grp),
  # summary_tab(x = dat_here$gender , label="Gender",type="cat",groups = grp)
  )
head(table_all)
my_printXtable(table_all, name="Baseline Covariates",N=nrow(dat_here),file=outTexReport)

head(dat_ARMA)

cat("\n  \\end{document} \n",file=outTexReport,append=T)


# For future use 
# Add variables one at a time 
# Baseline 
dat_ARMA$age = NA
dat_ARMA$gender = NA
dat_ARMA$ethnicity = NA
dat_ARMA$BMI = NA
dat_ARMA$ARDSrisk_pneumonia = NA
dat_ARMA$ARDSrisk_sepsis = NA
dat_ARMA$ARDSrisk_aspiration = NA
dat_ARMA$ARDSrisk_trauma = NA
dat_ARMA$ARDSrisk_otherCause = NA
dat_ARMA$ARDSrisk_multiTransfusion = NA
dat_ARMA$APACHE_III_tot = NA
dat_ARMA$vasopressor = NA
dat_ARMA$organFail_nonPulm = NA
dat_ARMA$organFail_circulatory = NA
dat_ARMA$organFail_coag = NA
dat_ARMA$organFail_renal = NA
dat_ARMA$organFail_creat = NA
dat_ARMA$pfRatio_base = NA
dat_ARMA$ARDS_severity = NA
dat_ARMA$ARDS_rapidResolve = NA
dat_ARMA$tidalVol = NA

# Outcomes 
dat_ARMA$died_90days = NA
dat_ARMA$nonPulmorganFailFreeDays_in28 = NA
dat_ARMA$icu_days = NA

