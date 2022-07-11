# These are all the packages necessary for checking values, creating visualizations, and running FA - you may not need all of them
library(tidyverse)
library(DataExplorer)
library(MVN) # MultiVariate Normality Test
library(readr)
library(GPArotation)
library(caret)
library(mifa)
library(psych)

# Read in data (change to your path) and remove duplicates
new_RWmaster <- read.csv("/Users/emdaly/Downloads/Master file_All_RWcage_data18-21_120522.csv",sep=';')
length(new_RWmaster$Ringnr[duplicated(new_RWmaster$Ringnr)])
# Remove the ~45 duplicates in data
RW_nodupe <- subset(new_RWmaster, new_RWmaster$tested==1)
# length(RW_nodupe$Ringnr[duplicated(RW_nodupe$Ringnr)]) // test and ensure we removed all dupes

# This standardizes the format in which the ring no. ID are in
RW_nodupe$Ringnr=gsub('\\s+', '',RW_nodupe$Ringnr)


##  COLUMN REMOVALS, first pass:  ##
# Missing data (>80%): present_T, latency_T, latency_M, present_M, bill_T, mobs_T, rasp_T
#    mobbing_T, mobs_M, churr_M, mobbing_M, churr_T, present_C, latency_C, bill_C, mobs_C,
#    rasp_C, churr_C, mobbing_C, X.1 (notes column)

# Irrelevant for analysis purposes: nr, year, NestID, Location, tested, species, fieldSex (once PCR sex column added)

# REDUNDANCIES: Total_time? (factored in later), Explored (two other variables say this),
#  Mov_sum (using Movpsec), Breathrate1 and Breathrate2 (using BR which is the average) 

# CHANGE: remove m_sum_start and m_number_movement and create m_Movpsec (m_number_movement/m_sum_start)

# Log10 normal transformation works for some but not all var - still 0 inflation problem as can be seen with:
# a=ggplot(RW_nodupe, aes(P1_time))
# b=ggplot(RW_nodupe, aes(log10(P1_time+1)))

##  COLUMN REMOVALS, second pass:  ##
# ADD: number of perches visited
# REMOVE: Middle_time/m_middle_time (overlapped too much with top and bottom time)
# CHANGE: split mpsec into homeside and exploration side mpsec, make perches into proportions instead

# Read in sequencing data
RAD_samples_and_metadata <- read_delim("Downloads/RAD_samples_and_metadata.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)
# Standardize IDs
RAD_samples_and_metadata$Ring_nro=gsub('\\s+', '',RAD_samples_and_metadata$Ring_nro)
# Merge by ID
rw_master=merge(RW_nodupe,RAD_samples_and_metadata,by.x=c("Ringnr"),by.y=c("Ring_nro"),all.x=TRUE)

# Read in nestling data
nestlingData <- read.csv("/Users/emdaly/Downloads/Master file_Adult_Nestlingdata2018-21 (1).csv",sep=';')
# Standardize IDs
nestlingData$Ringnr=gsub('\\s+', '',nestlingData$Ringnr)
# Merge by date AND ID
rw_master=merge(rw_master, nestlingData, by.x=c("Ringnr","date"), by.y=c("Ringnr","Date"), all.x=TRUE)

# Set up the variables discussed in the column changes section
rw_master$m_Movpsec <- rw_master$m_number_movement/rw_master$m_sum_start

rw_master$P1_pct <- rw_master$P1_time/rw_master$Time_exploration # Percentage of time spent on each perch
rw_master$P2_pct <- rw_master$P2_time/rw_master$Time_exploration
rw_master$P3_pct <- rw_master$P3_time/rw_master$Time_exploration
rw_master$P4_pct <- rw_master$P4_time/rw_master$Time_exploration
rw_master$P5_pct <- rw_master$P5_time/rw_master$Time_exploration

rw_master$n_explored <- apply(rw_master[,24:28], 1, function(x)sum(x!=0)) # Num perches explored (all nonzero entries)

rw_master$m_P1_pct <- rw_master$m_sum_TSP1/rw_master$m_sum_start # Percentage of time spent on each perch, mirror test
rw_master$m_P2_pct <- rw_master$m_sum_TSP2/rw_master$m_sum_start
rw_master$m_P3_pct <- rw_master$m_sum_TSP3/rw_master$m_sum_start
rw_master$m_P4_pct <- rw_master$m_sum_TSP4/rw_master$m_sum_start
rw_master$m_P5_pct <- rw_master$m_sum_TSP5/rw_master$m_sum_start

rw_master$mpsec_homeside <- rw_master$Moves_homeside/rw_master$Time_homeside # Moves per second homeside
rw_master$mpsec_exploration <- rw_master$Moves_exploration/rw_master$Time_exploration # Moves per second exploration side

## CREATING THE MOST RELIABLE SEX COLUMN (sexo) ##
# We start with genetic sex (most reliable) as Sex.y - then Sex.x for those missing - then Sex
# There are a lot of possible entries that do not boil down to M or F so we must coerce them all to 1/0 (see unique(rw_master$Sex.y) for some examples)
for(i in 1:nrow(rw_master)) {
  if (!is.na(rw_master$Sex.y[i])) {
    rw_master$sexo[i]=rw_master$Sex.y[i]
  }
  else if (rw_master$Sex.x[i] != "unknown" && rw_master$Sex.x[i] != "" && !is.na(rw_master$Sex.x[i])) {
    rw_master$sexo[i]=rw_master$Sex.x[i]
  }
  else if (rw_master$Sex[i] != "unknown" && !is.na(rw_master$Sex[i])) {
    rw_master$sexo[i]=rw_master$Sex[i]
  }
  else {
    rw_master$sexo[i]="sex_unknown"
  }
}
for(i in 1:nrow(rw_master)) {
  if (rw_master$sexo[i] == "M" || rw_master$sexo[i] == "male" || rw_master$sexo[i] == "M (in PCA, F field)" || rw_master$sexo[i] == "M (in PCA, F in field)" || rw_master$sexo[i] == "prob M") {
    rw_master$sexo[i]=1
  }
  else if (rw_master$sexo[i] == "F" || rw_master$sexo[i] == "female") {
    rw_master$sexo[i]=0
  }
  else {
    rw_master$sexo[i]=-99
  }
}
rw_master[rw_master=="NaN"]=0
rw_master[rw_master==Inf]=0


## CREATING THE ANALYSIS TABLE WITH ALL NUMERIC AND CONTINUOUS VARIABLES ##
rw_analysis=rw_master[1:300, (colnames(rw_master) %in% c("Ringnr","perches5_time","Bill_open","Bottom_time","Ground_time","Top_time",
                                                    "Latency_time","mpsec_homeside","mpsec_exploration", "Explored", "P1_pct","P2_pct",
                                                    "P3_pct","P4_pct","P5_pct","n_explored","Ruffles",
                                                    "Vocalisations",
                                                    "m_sum_atmirror","m_sum_bill","m_sum_looks","m_sum_low",
                                                    "m_sum_top","m_P1_pct","m_P2_pct","m_P3_pct","m_P4_pct",
                                                    "m_P5_pct","m_number_crown","m_number_display","m_number_flight",
                                                    "m_number_peck","m_number_shake","m_number_snap","m_number_song",
                                                    "m_number_vocal",
                                                    "Net_agg","Hand_agg","m_Movpsec","sexo"))]
# Make "Explored" a numeric factor
rw_analysis$Explored=factor(rw_analysis$Explored, levels=c("n", "y"), labels=c(0, 1))
rw_analysis$Explored=as.numeric(rw_analysis$Explored)-1
# Make any missing sex values NA
rw_analysis$sexo[rw_analysis$sexo==-99]=NA
rw_analysis$sexo=as.numeric(rw_analysis$sexo)

# Remove individuals with more than 40% missing data for this analysis
rowind=(rowMeans(is.na(rw_analysis)))*100
rowind=as.data.frame((rowMeans(is.na(rw_analysis)))*100)
RW_strict<-rw_analysis[-c(which(rowind>40)),]

# THEN, remove variables that still have more than 50% missing data even after individual removal
colind=(colMeans(is.na(RW_strict)))*100
colind=as.data.frame((colMeans(is.na(RW_strict)))*100)
RW_strict<-RW_strict[,-c(which(colind>50))]


# Compute multiple imputation for all the missing values in our analysis.
# See help(mifa) for documentation.

mi <- mifa(
  data  = RW_strict[,2:38],
  ci = "fieller",
  print = TRUE
)

## IGNORE THIS BLOCK - (OLD) MANUAL MULTIPLE IMPUTATION METHOD ##
# sum(nas)
# imp <- mice(RW_strict[,nas], meth='pmm', print=F)
# RW_imp<-complete(imp,1)
# If you have 10% missingness - do 10 MIs, if you have 20% do 20, etc. Cutoff for Caitlin was 25% (same here)
# Always good to do multiple runs (don't just take the first), pay attention to the MI calculation you are using, what does or does not make sense to assume?
# Is there some dependency I am missing?

# now complete dataset (this is poorly written I apologize) of 15 variables that had previously missing values
# maximum likelihood needs to specify distribution more than ordinarily - FA itself does not care as much since it looks at the second moment
# place Ray expects that to matter is in a maximum likelihood direction
# you can specify logit/log etc. link function
# lavaan package - explore
# bias vs variance - how much complexity should you add to the model to make it more robust?

# heywood cases - statistical software spits out best estimate that it knows is technically impossible (you can get some understanding, even if its a little far outside the mathematical bounds)
# structure of items - near-collinearity - too many factors (in EFA)
# RW_final=RW_strict
# RW_final$perches5_time=RW_imp$perches5_time
# RW_final$Latency_time=RW_imp$Latency_time
# RW_final$Explored=RW_imp$Explored
# RW_final$Middle_time=RW_imp$Middle_time
# RW_final$m_sum_low=RW_imp$m_sum_low
# RW_final$m_sum_middle=RW_imp$m_sum_middle
# RW_final$m_sum_top=RW_imp$m_sum_top
# RW_final$m_sum_TSP1=RW_imp$m_sum_TSP1
# RW_final$m_sum_TSP2=RW_imp$m_sum_TSP2
# RW_final$m_sum_TSP3=RW_imp$m_sum_TSP3
# RW_final$m_sum_TSP4=RW_imp$m_sum_TSP4
# RW_final$m_sum_TSP5=RW_imp$m_sum_TSP5
# RW_final$BR=RW_imp$BR
# RW_final$m_Movpsec=RW_imp$m_Movpsec
# RW_final$sexo=RW_imp$sexo
## END OLD BLOCK ##

## EXPLORATORY FACTOR ANALYSIS USING PSYCH PACKAGE ##

# Parallel analysis to determine number of factors
fa.parallel(mi$cov_combined, n.obs=nrow(RW_strict), fm="ml")
# That result can be compared if desired to the Very Simple Structure determination of nfactors
vss(mi$cov_combined, n.obs=nrow(RW_strict), fm="ml")

# This was the one for which I saved the scree plot - has the highest reliability score and no warnings
maxlik=fa(mi$cov_combined, nfactors = 5, n.obs=nrow(RW_strict), rotate = "oblimin", fm="ml")

# Create heatmap of factor loadings
heatmap(maxlik$loadings,Colv=NA, Rowv=NA, col=cm.colors(256))
legend(x="bottomright", cex=0.5,legend=c("neg", "zero", "pos"), fill=(cm.colors(3)))



##  EXPLORATION BEHAVIOR ANALYSIS  ##

# Separate out only exploration analysis variables
rw_analysis=rw_master[1:300, (colnames(rw_master) %in% c("Ringnr","perches5_time","Bill_open","Bottom_time","Ground_time",
                                                    "Latency_time","Explored","P1_pct","P2_pct",
                                                    "P3_pct","P4_pct","P5_pct","Ruffles","Top_time",
                                                    "Vocalisations","sexo","mpsec_homeside","mpsec_exploration","n_explored"))]
# Process variables like before
rw_analysis$Explored=factor(rw_analysis$Explored, levels=c("n", "y"), labels=c(0, 1))
rw_analysis$Explored=as.numeric(rw_analysis$Explored)-1
rw_analysis$sexo[rw_analysis$sexo==-99]=NA
rw_analysis$sexo=as.numeric(rw_analysis$sexo)

# Remove high-missingness entries as before
rowind=(rowMeans(is.na(rw_analysis)))*100
rowind=as.data.frame((rowMeans(is.na(rw_analysis)))*100)
RW_strict<-rw_analysis[-c(which(rowind>40)),]


# Multiple imputation
mi <- mifa(
  data  = RW_strict[,2:19],
  ci = "fieller",
  print = TRUE
)
# RW_strict[is.na(RW_strict)]=0

# OLD MI
# expRW_final=RW_strict
# expRW_final$perches5_time=RW_imp$perches5_time
# expRW_final$Latency_time=RW_imp$Latency_time
# expRW_final$Explored=RW_imp$Explored
# expRW_final$Middle_time=RW_imp$Middle_time
# expRW_final$BR=RW_imp$BR
# expRW_final$sexo=RW_imp$sexo

# expRW_corr=cor(expRW_final[,!(names(expRW_final) %in% c("Ringnr"))])
# highCorr=findCorrelation(expRW_corr,cutoff=0.8,names=T)
# expRW_final=expRW_final[,!(names(expRW_final) %in% highCorr)] # Almost perfectly correlated so needed removal 

# FACTOR ANALYSIS
fa.parallel(mi$cov_combined, n.obs=nrow(RW_strict))
vss(mi$cov_combined, n.obs=nrow(RW_strict))

# This was the one for which I saved the scree plot - has the highest reliability score and no warnings
maxlik=fa(mi$cov_combined, nfactors = 4, n.obs=nrow(RW_strict), rotate = "oblimin", fm="ml")

# Generate heatmap
heatmap(maxlik$loadings,Colv=NA, Rowv=NA, col=cm.colors(256),main="exploration behavior")
legend(x="bottomright", cex=0.5,legend=c("neg", "zero", "pos"), fill=(cm.colors(3)))



## MIRROR BEHAVIOR ANALYSIS ##

# Seperate out only mirror test variables
rw_analysis=rw_master[1:300, (colnames(rw_master) %in% c("Ringnr","m_sum_atmirror","m_sum_bill","m_sum_looks","m_sum_low",
                                                    "m_sum_top","m_P1_pct","m_P2_pct","m_P3_pct","m_P4_pct",
                                                    "m_P5_pct","m_number_crown","m_number_display","m_number_flight",
                                                    "m_number_peck","m_number_shake","m_number_snap","m_number_song",
                                                    "m_number_vocal",
                                                    "m_Movpsec","sexo"))]
# Process variables as before
rw_analysis$sexo[rw_analysis$sexo==-99]=NA
rw_analysis$sexo=as.numeric(rw_analysis$sexo)

# Remove high missingness variables as before
rowind=(rowMeans(is.na(rw_analysis)))*100
rowind=as.data.frame((rowMeans(is.na(rw_analysis)))*100)
RW_strict<-rw_analysis[-c(which(rowind>40)),]

colind=(colMeans(is.na(RW_strict)))*100
colind=as.data.frame((colMeans(is.na(RW_strict)))*100)
RW_strict<-RW_strict[,-c(which(colind>50))]

# Multiple imputation
mi <- mifa(
  data  = rw_analysis[,2:21],
  ci = "fieller",
  print = TRUE
)

# OLD MI
# nas=apply(RW_strict, 2, function(x) any(is.na(x)))
# sum(nas)
# imp <- mice(RW_strict[,nas], meth='pmm', print=F)
# RW_imp<-complete(imp,1)

# mRW_final=RW_strict
# mRW_final$BR=RW_imp$BR
# mRW_final$sexo=RW_imp$sexo

# mRW_corr=cor(mRW_final[,!(names(mRW_final) %in% c("Ringnr"))])
# highCorr=findCorrelation(mRW_corr,cutoff=0.8,names=T)
# mRW_final=mRW_final[,!(names(mRW_final) %in% highCorr)] # Almost perfectly correlated so needed removal 

# cortest.bartlett(mRW_corr,n=232)

fa.parallel(mRW_final[,!(names(mRW_final) %in% c("Ringnr"))])
vss(mRW_final[,!(names(mRW_final) %in% c("Ringnr"))])

# This was the one for which I saved the scree plot - has the highest reliability score and no warnings
maxlik=fa(mi$cov_combined, nfactors = 4, n.obs=nrow(RW_strict), rotate = "oblimin", fm="ml")
heatmap(maxlik$loadings,Colv=NA, Rowv=NA, col=cm.colors(256),main="mirror behavior")
legend(x="bottomright", cex=0.5,legend=c("neg", "zero", "pos"), fill=(cm.colors(3)))



testset_master$m_Movpsec <- testset_master$m_number_movement/testset_master$m_sum_start

testset_master$P1_pct <- testset_master$P1_time/testset_master$Time_exploration # Percentage of time spent on each perch
testset_master$P2_pct <- testset_master$P2_time/testset_master$Time_exploration
testset_master$P3_pct <- testset_master$P3_time/testset_master$Time_exploration
testset_master$P4_pct <- testset_master$P4_time/testset_master$Time_exploration
testset_master$P5_pct <- testset_master$P5_time/testset_master$Time_exploration

testset_master$n_explored <- apply(testset_master[,24:28], 1, function(x)sum(x!=0)) # Num perches explored (all nonzero entries)

testset_master$m_P1_pct <- testset_master$m_sum_TSP1/testset_master$m_sum_start # Percentage of time spent on each perch
testset_master$m_P2_pct <- testset_master$m_sum_TSP2/testset_master$m_sum_start
testset_master$m_P3_pct <- testset_master$m_sum_TSP3/testset_master$m_sum_start
testset_master$m_P4_pct <- testset_master$m_sum_TSP4/testset_master$m_sum_start
testset_master$m_P5_pct <- testset_master$m_sum_TSP5/testset_master$m_sum_start

testset_master$mpsec_homeside <- testset_master$Moves_homeside/testset_master$Time_homeside # Moves per second homeside
testset_master$mpsec_exploration <- testset_master$Moves_exploration/testset_master$Time_exploration # Moves per second exploration side

# Top_time and Latency_time - staying the same but will be added to analysis

# We start with genetic sex (most reliable) - unique(testset_master$Sex.y)
for(i in 1:nrow(testset_master)) {
  if (!is.na(testset_master$Sex.y[i])) {
    testset_master$sexo[i]=testset_master$Sex.y[i]
  }
  else if (testset_master$Sex.x[i] != "unknown" && testset_master$Sex.x[i] != "" && !is.na(testset_master$Sex.x[i])) {
    testset_master$sexo[i]=testset_master$Sex.x[i]
  }
  else if (testset_master$Sex[i] != "unknown" && !is.na(testset_master$Sex[i])) {
    testset_master$sexo[i]=testset_master$Sex[i]
  }
  else {
    testset_master$sexo[i]="sex_unknown"
  }
}
for(i in 1:nrow(testset_master)) {
  if (testset_master$sexo[i] == "M" || testset_master$sexo[i] == "male" || testset_master$sexo[i] == "M (in PCA, F field)" || testset_master$sexo[i] == "M (in PCA, F in field)" || testset_master$sexo[i] == "prob M") {
    testset_master$sexo[i]=1
  }
  else if (testset_master$sexo[i] == "F" || testset_master$sexo[i] == "female") {
    testset_master$sexo[i]=0
  }
  else {
    testset_master$sexo[i]=-99
  }
}
testset_master[testset_master=="NaN"]=0
testset_master[testset_master==Inf]=0

testset_analysis=testset_master[, (colnames(testset_master) %in% c("Ringnr","perches5_time","Bill_open","Bottom_time","Ground_time",
                                                                   "Latency_time","Explored","P1_pct","P2_pct",
                                                                   "P3_pct","P4_pct","P5_pct","Ruffles","Top_time",
                                                                   "Vocalisations","sexo","mpsec_homeside","mpsec_exploration","n_explored"))]
testset_analysis$Explored=factor(testset_analysis$Explored, levels=c("n", "y"), labels=c(0, 1))
testset_analysis$Explored=as.numeric(testset_analysis$Explored)-1
testset_analysis$sexo[testset_analysis$sexo==-99]=NA
testset_analysis$sexo=as.numeric(testset_analysis$sexo)

rowind=(rowMeans(is.na(testset_analysis)))*100
rowind=as.data.frame((rowMeans(is.na(testset_analysis)))*100)
testset_strict<-testset_analysis[-c(which(rowind>40)),]

hs.mod <- 'F1 =~ perches5_time + Latency_time + mpsec_homeside
           F2 =~ Bottom_time + Ground_time + Top_time + P4_pct
           F3 =~ Latency_time + P5_pct + n_explored + mpsec_exploration + P1_pct + P2_pct + P3_pct'

gwasvars=rw_master[, (colnames(rw_master) %in% c("Ringnr","BR","Aggression.hand","Aggression.net","Aggression_total"))]
