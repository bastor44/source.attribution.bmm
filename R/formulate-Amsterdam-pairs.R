require(tidyverse)
require(data.table)
require(ggplot2)
require(rstan)
require(knitr)
require(bayesplot)
library(lubridate)
require(dplyr)
require(ggsci)
#source('R/functions.R')

analysis <- 'analysis_220713'
results <- 'agegps_sensanalysis_210216_MSM-2010_2022'

out.dir <- file.path('/Users/bethanyastor/Library/CloudStorage/OneDrive-ImperialCollegeLondon/Research Project/source.attribution.bmm/out_Amsterdam', results)
clock_model <- '/Users/bethanyastor/Library/CloudStorage/OneDrive-ImperialCollegeLondon/Research Project/source.attribution.bmm/molecular_clock'

args <- list(
  source_dir = '/Users/bethanyastor/Library/CloudStorage/OneDrive-ImperialCollegeLondon/Research Project/source.attribution.bmm',
  indir = 'Users/bethanyastor/Library/CloudStorage/OneDrive-ImperialCollegeLondon/Research Project/source.attrribution.bmm', 
  mig_groups=T, 
  trsm='MSM', 
  sens=T # for sensitivity analysis with patristic distances of 0 
)
source(file.path(args$source_dir, 'R', 'functions_simulation_scenarios.R'))


if(1) dir.create( out.dir )

`%notin%` <- negate(`%in%`)

infile.seq <-	file.path(args$indir, 'Data', 'data_220331/SHM_2201_ROADMAP_220331_tblLAB_seq.rda')
infile.meta <- file.path(args$indir, analysis, 'misc', '220713_sequence_labels.rda')
infile.bas <- file.path(args$indir, 'Data', 'data_220331','SHM_2201_ROADMAP_220331_tblBAS.csv')

load(infile.seq)
load(infile.meta)
dind <- data.table(dind)
dind[, SEQ:= PATIENT %in% ds$PATIENT]
dind <- unique(dind)

## load infection time estimates and metadata ----

dinf <- data.table(read.csv(file.path('data_Ams',analysis,'Infection_date_est_rec.csv')))
setnames(dinf,c("id",'estsctodiagMedian','estsctodiagLL','estsctodiagUL'),c("TO_SEQUENCE_ID",'SER_TO_DIAG','SER_TO_DIAG_LL','SER_TO_DIAG_UL'))
dinf <- unique(dinf)
dinf <- merge(dinf,subset(dind,select=c('PATIENT','CITY','SEQ','TRANSM')),
                                                by.x='TO_SEQUENCE_ID',by.y='PATIENT',all.x=T)
dinf[, SEQ:= TO_SEQUENCE_ID %in% ds$PATIENT]

meta_data <- readRDS(file.path('data_Ams',analysis,'meta_data_mg_country.rds'))
setnames(meta_data,"ID","TO_SEQUENCE_ID")
setnames(meta_data,"bplace","BPLACE")
meta_data <- unique(meta_data)
data_age <- readRDS(file.path('data_Ams',analysis,'data_age.rds'))
data_age <- as.data.table(data_age)
setnames(data_age,c("TO_SEQUENCE_ID","BIRTH_DATE","BIRTH_DATE_DEC"))

# load baseline data
dbas <- data.table(read.csv(infile.bas))

# calculate infection date
dinf[,DIAGNOSIS_DATE:= as.Date(dinf[,hiv_pos_d],format="%Y-%m-%d")]
dinf[,DIAGNOSIS_DATE_N:= hivc.db.Date2numeric(dinf[,DIAGNOSIS_DATE])]
dinf[,EST_INF_DATE:= DIAGNOSIS_DATE_N-SER_TO_DIAG]
dinf[,EST_INF_DATE:= format(date_decimal(EST_INF_DATE), "%Y-%m-%d")]
dinf[,YEAR_OF_INF_EST := year(EST_INF_DATE)]

# add risk group and summarise MSM
cat(paste0('Amsterdam MSM with estimated infection date in 2010-2021: ',nrow(unique(subset(dinf,CITY=='Amsterdam' & TRANSM=='MSM' & YEAR_OF_INF_EST >= 2010 & YEAR_OF_INF_EST<2022)))))
cat(paste0('Amsterdam MSM with estimated infection date in 2010-2021 with a sequence: ',nrow(unique(subset(dinf,CITY=='Amsterdam' & TRANSM=='MSM' &  YEAR_OF_INF_EST >= 2010 & YEAR_OF_INF_EST<2022 & SEQ==T)))))
tmp <- unique(subset(dinf,CITY=='Amsterdam' & TRANSM=='MSM' & YEAR_OF_INF_EST >= 2010 & YEAR_OF_INF_EST<2022))
tmp <- merge(tmp,subset(dbas,select=c('PATIENT','BIRTH_D')),by.x='TO_SEQUENCE_ID',by.y='PATIENT',all.x=T)
tmp$EST_INF_DATE <- as.Date(ISOdate(tmp$EST_INF_DATE, 1, 1))
tmp$BIRTH_D <- as.Date(ISOdate(tmp$BIRTH_D, 1, 1))
tmp[,TO_AGE:= as.numeric(EST_INF_DATE-BIRTH_D)/365]
cat(paste0('Median age of Amsterdam MSM with estimated infection date in 2010-2021: ',median(tmp$TO_AGE)))

# create a table containing all data on recipients
individual_data <- merge(data_age,meta_data,by="TO_SEQUENCE_ID")
individual_data <- merge(individual_data,dinf,by="TO_SEQUENCE_ID")
setnames(individual_data,c("EST_INF_DATE"),c("INFECTION_DATE_EST"))

# format the columns containing date to have date format
individual_data$INFECTION_DATE_EST <- as.Date(individual_data$INFECTION_DATE_EST)
individual_data$DIAGNOSIS_DATE <- as.Date(individual_data$DIAGNOSIS_DATE)

# keep only individuals with complete information of infection date
individual_data_complete <- individual_data[which(!is.na(INFECTION_DATE_EST))]

# add column for year of diagnosis and year of infection #
individual_data_complete[,YEAR_OF_DIAG := year(DIAGNOSIS_DATE)]
individual_data_complete[,YEAR_OF_INF_EST := year(INFECTION_DATE_EST)]

# histogram of risk groups with count #
tmp <- as.data.table(individual_data_complete %>% group_by(RISK_GROUP) %>% summarise(count=n()))
tmp$RISK_GROUP <- factor(tmp$RISK_GROUP,levels=c('HSX','MSM'),labels=c('Heterosexuals', 'Men who have sex \n with men'))

p <- ggplot(tmp,aes(x=RISK_GROUP,y=count,fill=RISK_GROUP))+geom_bar(stat="identity",width=1)+theme_bw()+scale_fill_aaas()+
  labs(x="\n Individual's risk group", y="Sample size \n")+theme(legend.position = "none")+geom_text(aes(label=count),vjust=-0.5,size=3)
p
ggsave(file=file.path(out.dir,'Hist_trsm_gps_complete_individual_data.pdf'), p, w=5.5, h=4)

## subset by risk group ----

# Subset the table to contain only MSM recipients
individual_data_complete <- individual_data_complete[which(individual_data_complete$RISK_GROUP==args$trsm)]

# Plot of individuals per year (cumulative sum)
dat <- as.data.table(individual_data_complete %>% group_by(YEAR_OF_INF_EST)  %>% summarise(count=n()))
dat[order(-YEAR_OF_INF_EST),cumulative:= cumsum(count)]

p <- ggplot(dat,aes(x=as.factor(YEAR_OF_INF_EST),y=cumulative))+geom_bar(stat="identity",fill="dark blue")+theme_bw()+labs(x="\n Starting year for including infected individuals in analysis", y= "Sample size \n")+
  theme(axis.text.x = element_text(angle=45,vjust=0.5))
p
ggsave(file=file.path(out.dir,'Barplot_cum_sum_sample_size.pdf'), p, w=7, h=4)

dat <- as.data.table(individual_data_complete %>% group_by(YEAR_OF_INF_EST, BPLACE_MG)  %>% summarise(count=n()))
dat <- dat[order(BPLACE_MG,YEAR_OF_INF_EST),]
dat <- dat[, list(YEAR_OF_INF_EST=YEAR_OF_INF_EST,cumulative = cumsum(count)),by='BPLACE_MG']
dat[, mwmb:= factor(BPLACE_MG,levels=c('NL','G1','G2','G3','G4','G5','Other'),
                    labels=c('Netherlands','W.Europe, N.America, Oceania','E. & C. Europe', 'Caribbean & S. America','Sub-Saharan Africa','Caribbean & S. America','Other'))]
p <- ggplot(dat,aes(x=YEAR_OF_INF_EST,y=cumulative,fill=mwmb)) +
  geom_bar(stat="identity") +
  scale_fill_npg() +
  facet_wrap(mwmb~.) +
  scale_x_continuous(expand = c(0,0), breaks=seq(1980,2021,2),labels=seq(1980,2021,2)) +
  theme_bw() + labs(x="\n Starting year for including infected individuals in analysis", y= "Sample size \n") +
  theme(axis.text.x = element_text(angle=45,vjust=0.5),
        strip.background=element_blank(), legend.position="none")
p
ggsave(file=file.path(out.dir,'Barplot_cum_sum_sample_size_MG.pdf'), p, w=10, h=6)

## formulate pairs ----

dat <- individual_data_complete
dat[SUBTYPE=='B',CLUSTER_NUMBER:=paste0(CLUSTER_NUMBER,'_',ST_CLADE)]

# summarise cluster size
tmp <- dat[, list(CLU_SIZE=length(TO_SEQUENCE_ID)),by='CLUSTER_NUMBER']
dat <- merge(dat,tmp,by='CLUSTER_NUMBER',all.x=T)
cat(paste0('Number of Amsterdam MSM singletons: ',length(unique(dat$TO_SEQUENCE_ID[dat$CLU_SIZE==1 & dat$YEAR_OF_INF_EST >= 2010 & dat$YEAR_OF_INF_EST<2022]))))

singletons <- unique(dat$TO_SEQUENCE_ID[dat$CLU_SIZE==1 & dat$YEAR_OF_INF_EST >= 2010 & dat$YEAR_OF_INF_EST<2022])

# formulate all pairs first
pairs <- dat %>%
  setNames(paste0(names(.), '_2')) %>%
  tidyr::crossing(dat) %>%
filter(TO_SEQUENCE_ID != TO_SEQUENCE_ID_2 & INFECTION_DATE_EST_2<INFECTION_DATE_EST) # form all pairs ignoring clustering
pairs <- as.data.table(pairs)

# format table for analysis [subset only variables of interest]
pairs <- subset(pairs,select = c("TO_SEQUENCE_ID_2","BIRTH_DATE_2","SUBTYPE_2","CLUSTER_NUMBER_2","LOCATION_2","TREATMENT_DATE_2","SAMPLING_DATE_2","INFECTION_DATE_EST_2","DIAGNOSIS_DATE_2","SER_TO_DIAG_2","SER_TO_DIAG_LL_2","SER_TO_DIAG_UL_2","ORIGIN_2","BIRTH_COUNTRY_2",'BPLACE_MG_2',
                                 "TO_SEQUENCE_ID","BIRTH_DATE","SUBTYPE","CLUSTER_NUMBER","LOCATION","TREATMENT_DATE","SAMPLING_DATE","INFECTION_DATE_EST","DIAGNOSIS_DATE","SER_TO_DIAG","SER_TO_DIAG_LL","SER_TO_DIAG_UL","ORIGIN","BIRTH_COUNTRY",'BPLACE_MG'))
setnames(pairs,names(pairs),c("FROM_SEQUENCE_ID","FROM_BIRTH_DATE","FROM_SUBTYPE","FROM_CLUSTER_NUMBER","FROM_LOCATION","FROM_TREATMENT_DATE","FROM_SAMPLING_DATE","FROM_EST_INFECTION_DATE","FROM_DIAGNOSIS_DATE","FROM_TIME_INF_TO_DIAG","FROM_TIME_INF_TO_DIAG_LL","FROM_INF_TO_DIAG_UL","FROM_ORIGIN","FROM_COUNTRY","FROM_BPLACE",
                              "TO_SEQUENCE_ID","TO_BIRTH_DATE","TO_SUBTYPE","TO_CLUSTER_NUMBER","TO_LOCATION","TO_TREATMENT_DATE","TO_SAMPLING_DATE","TO_EST_INFECTION_DATE","TO_DIAGNOSIS_DATE","TO_INF_TO_DIAG","TO_INF_TO_DIAG_LL","TO_INF_TO_DIAG_UL","TO_ORIGIN","TO_COUNTRY","TO_BPLACE"))

# save all pairs in cohort before selecting based on date
saveRDS(pairs,file=file.path(out.dir, 'all_pairs_all_time.rds'))

### subset to pairs in which recipient was infected since 2010 ----
pairs[, YEAR_OF_INF_EST:= year(TO_EST_INFECTION_DATE)]
pairs <- pairs %>% filter(YEAR_OF_INF_EST>=2010 & YEAR_OF_INF_EST<2022) %>% mutate(YEAR_RANGE="2010-2021")

### remove recipients who are singletons ----
pairs <- subset(pairs, TO_SEQUENCE_ID %notin% singletons)

# number of unique sources
cat(paste0('Number of pairs before exclusions: ',nrow(pairs)))
cat(paste0('Number of sources: ',length(unique(pairs$FROM_SEQUENCE_ID))))
cat(paste0('Number of cases: ',length(unique(pairs$TO_SEQUENCE_ID))))

# save all pairs
saveRDS(pairs,file=file.path(out.dir, 'all_pairs.rds'))

# exclusion criteria ----

### remove any pairs in which transmitter died before infection date of recipient ----
pairs <- merge(pairs,subset(dbas,select=c('PATIENT','MIG_D','MIG_D_aq','DEATH_D')),by.x='FROM_SEQUENCE_ID',by.y='PATIENT',all.x=T)
pairs[, DEATH_D := as.Date(DEATH_D,format="%Y-%m-%d")]

cat(paste0('Number of sources who died before infection date of recipient: ',nrow(subset(pairs,!is.na(DEATH_D) & DEATH_D<TO_EST_INFECTION_DATE))))

pairs <- subset(pairs,DEATH_D>=TO_EST_INFECTION_DATE | is.na(DEATH_D))

### exclude sources who changed residence ----
# remove any pairs in which the transmitter had not yet migrated into NL before infection date of the recipient

pairs[, MIG_D := as.Date(MIG_D,format="%Y-%m-%d")]
pairs[, MIG_D_aq := as.Date(MIG_D_aq,format="%Y-%m-%d")]

# use last possible migration date in window provided, (MIG_D_aq) where date is uncertain
pairs[, FLAG:= 0]
pairs[!is.na(MIG_D) & MIG_D_aq>TO_EST_INFECTION_DATE, FLAG:= 1]
# except for those whose migration date is certain (no uncertainty window), use MIG_D
pairs[!is.na(MIG_D) & is.na(MIG_D_aq) & MIG_D>TO_EST_INFECTION_DATE, FLAG:= 1]

cat(paste0('Number of sources who arrived after infection date of recipient: ',nrow(subset(pairs,FLAG==1))))

pairs <- subset(pairs, FLAG==0)

### exclude suppressed sources ----

infile.rna <-	file.path(indir_data, 'Data', 'data_220331/SHM_2201_ROADMAP_220331_tblLAB_RNA.csv')
dat <- read.csv(infile.rna,header=T)
dat$RNA_D <- as.Date(dat$RNA_D,format=c("%Y-%m-%d"))
dat <- data.table(dat)
# fit model/curve through each patient's VLs to get trajectory and interpolate
# set any values under/over threshold to one below/above the threshold
dat[RNA_V==-1 & !is.na(RNA_L), RNA_V:= RNA_L - 1]
dat[RNA_V==-1 & !is.na(RNA_UL), RNA_V:= RNA_UL + 1]
dat[RNA_V==-1, RNA_V:= 0]
dat[RNA_V==-1000 & RNA_L==-999, RNA_V:= -1]
# copy mean of each person's last 2 measurements to end of follow-up so we don't have missing values
dat2 <- dat %>%
  group_by(PATIENT) %>%
  slice_tail(n = 2)
dat2 <- data.table(dat2)
dat2 <- dat2[, list(RNA_V=mean(RNA_V,na.rm=T)),by='PATIENT']
dat2[, RNA_D:= max(dat$RNA_D)] # add date of last obs to impute until
dat <- merge(dat,dat2,by=c('PATIENT','RNA_D','RNA_V'),all=T)

# just keep patients who are a probable transmitter

dat <- subset(dat,PATIENT %in% pairs$FROM_SEQUENCE_ID)
dat <- merge(dat,subset(pairs,select=c('FROM_SEQUENCE_ID','TO_EST_INFECTION_DATE')),by.x=c('PATIENT','RNA_D'),by.y=c('FROM_SEQUENCE_ID','TO_EST_INFECTION_DATE'),all=T)

# remove patients with less than 4 measurements
tmp <- unique(subset(dat,select=c('PATIENT','RNA_D','RNA_V')))
dn <- tmp[, list(N=length(RNA_V[!is.na(RNA_V)])),by=c('PATIENT')]
dat <- merge(dat,dn,by='PATIENT',all.x=T)
dat <- subset(dat, N>3)
dl <- subset(dat,PATIENT %notin% c('912511','904491','919461')) %>% # exclude  patients with problems fitting loess
  group_by(PATIENT) %>%
  arrange(PATIENT, RNA_D) %>%
  nest() %>%
  mutate(
    pred.response = purrr::map(data, function(x)stats::loess(RNA_V~as.numeric(RNA_D), span= 0.5, data = x,na.action="na.omit") %>%
                                 stats::predict(data.frame(RNA_D = as.numeric(x$RNA_D))))) %>%
  unnest(cols = c(data, pred.response))
dl <- data.table(dl)

# plot fitted loess and predicted VLs at time of infections of all recipients for sample of 10 patients
ggplot(subset(dl, PATIENT %in% unique(dl$PATIENT)[1:10])) +
  geom_line(aes(x=RNA_D,y=RNA_V,col=as.factor(PATIENT))) +
  geom_line(aes(x=RNA_D,y=pred.response,col=as.factor(PATIENT)),linetype=2) +
  labs(x='Date of measurement',y='Viral load',colour='Patient ID') +
  scale_y_continuous(labels = scales::comma) +
  theme_bw(base_size=40) +
  theme(legend.position="bottom")
ggsave(file=file.path(out.dir,'viral_load_trajectories_pred_10pts.png'),w=25,h=12)

# flag patients virally suppressed by putative infection date of recipient
tmp <- data.table(unique(subset(dl,select=c('PATIENT','RNA_D','pred.response'))))
tmp[pred.response<=200, supp:=1] # updated from 100 based on SHM report defn of suppression
tmp[pred.response>200, supp:=0]

pairs <- merge(pairs,tmp,by.x=c('FROM_SEQUENCE_ID','TO_EST_INFECTION_DATE'),by.y=c('PATIENT','RNA_D'),all.x=T)

cat(paste0('Number of sources who were likely durably virally suppressed on infection date of recipient: ',nrow(subset(pairs,supp==1))))

pairs <- subset(pairs,supp==0 | is.na(supp))

cat(paste0('Number of unique pairs: ',length(unique(pairs))))
cat(paste0('Number of unique sources: ',length(unique(pairs$FROM_SEQUENCE_ID))))
cat(paste0('Number of unique cases: ',length(unique(pairs$TO_SEQUENCE_ID))))

### Remove pairs whose time elapsed is larger than 18 year ----

pairs[,TIME_ELAPSED := as.numeric((TO_SAMPLING_DATE-TO_EST_INFECTION_DATE)+abs(FROM_SAMPLING_DATE-TO_EST_INFECTION_DATE))/365]

cat(paste0('Number of pairs with time elapsed > 16 years: ',nrow(subset(pairs,TIME_ELAPSED>=16))))

pairs <- pairs[which(TIME_ELAPSED<16),]

cat(paste0('Number of unique pairs: ',length(unique(pairs))))
cat(paste0('Number of unique sources: ',length(unique(pairs$FROM_SEQUENCE_ID))))
cat(paste0('Number of unique cases: ',length(unique(pairs$TO_SEQUENCE_ID))))

### remove pairs not in same phylo subgraph ----
cat(paste0('Number of sources in different phylo subgraph to recip: ',nrow(pairs[FROM_CLUSTER_NUMBER!=TO_CLUSTER_NUMBER])))
pairs <- subset(pairs,FROM_CLUSTER_NUMBER==TO_CLUSTER_NUMBER)

cat(paste0('Number of unique pairs: ',nrow(pairs)))
cat(paste0('Number of unique sources: ',length(unique(pairs$FROM_SEQUENCE_ID))))
cat(paste0('Number of unique cases: ',length(unique(pairs$TO_SEQUENCE_ID))))

## add genetic distance and calculate time elapsed ----

# Integrate genetic distance of each pair from the distances found from the maximum likelihood phylogenetic tree
gen_dist <- readRDS(file.path('data_Ams',analysis,paste0('pairwise_dist_allSTs_',args$trsm,'.rds')))
pairs <- merge(subset(gen_dist,select = c("FROM_SEQUENCE_ID","TO_SEQUENCE_ID","distance")),pairs,by= c("FROM_SEQUENCE_ID","TO_SEQUENCE_ID"),all.y=T)
setnames(pairs,"distance","GEN_DIST")
# sensitivity analysis: replace distances of 0 with 1 mutation across alignment (1/1302)
if(args$sens==T){
  pairs[is.na(GEN_DIST),GEN_DIST:= 0.00077]
}else{
  pairs <- subset(pairs,!is.na(GEN_DIST))
}

## add extra information ----
# Determine the stage group of the transmitter at the estimated time of infection based on his status
pairs[which(FROM_DIAGNOSIS_DATE<=TO_EST_INFECTION_DATE & !is.na(FROM_TREATMENT_DATE)),TRANS_STAGE:='DIAGNOSED']
pairs[which(FROM_DIAGNOSIS_DATE>TO_EST_INFECTION_DATE & !is.na(FROM_TREATMENT_DATE)),TRANS_STAGE:='UNDIAGNOSED']
# For transmitters with unknown treatment date
pairs[which(is.na(TRANS_STAGE) & FROM_DIAGNOSIS_DATE>TO_EST_INFECTION_DATE),TRANS_STAGE:='UNDIAGNOSED']
pairs[which(FROM_DIAGNOSIS_DATE<=TO_EST_INFECTION_DATE & is.na(TRANS_STAGE)),TRANS_STAGE:='DIAGNOSED']

# Add age of the transmitter at the estimated time of infection
pairs$FROM_BIRTH_DATE <- as.Date(ISOdate(pairs$FROM_BIRTH_DATE, 1, 1))
pairs$TO_BIRTH_DATE <- as.Date(ISOdate(pairs$TO_BIRTH_DATE, 1, 1))
pairs[,FROM_AGE:= as.numeric(TO_EST_INFECTION_DATE-FROM_BIRTH_DATE)/365]
pairs[,TO_AGE:= as.numeric(TO_EST_INFECTION_DATE-TO_BIRTH_DATE)/365]
pairs[, FROM_AGE_GP:= cut(FROM_AGE,breaks=c(15,30,40,50,60,100),include.lowest=T,right=F,
                          labels=c('[15-30)','[30-40)','[40-50)','[50-60)','[60+)'))]
pairs[, TO_AGE_GP:= cut(TO_AGE,breaks=c(15,30,40,50,60,100),include.lowest=T,right=F,
                        labels=c('[15-30)','[30-40)','[40-50)','[50-60)','[60+)'))]
pairs[,PAIR_ID:=seq(1,nrow(pairs))]

## save pairs
cat(paste0('Number of pairs: ',nrow(pairs)))
cat(paste0('Number of unique sources: ',length(unique(pairs$FROM_SEQUENCE_ID))))
cat(paste0('Number of unique cases: ',length(unique(pairs$TO_SEQUENCE_ID))))

saveRDS(pairs,file=file.path(out.dir, paste0(args$trsm,'_pairs.rds')))

cat(paste0('Mean age of sources: ',round(mean(pairs$FROM_AGE),0)))
tmp <- unique(subset(pairs,select=c(FROM_SEQUENCE_ID,FROM_AGE)))
cat(paste0('Median age of sources: ',round(quantile(tmp$FROM_AGE,probs=0.5),0),
           ' [',round(quantile(tmp$FROM_AGE,probs=0.25),0),'-',
           round(quantile(tmp$FROM_AGE,probs=0.75),0),']'))
tmp <- unique(subset(pairs,select=c(TO_SEQUENCE_ID,TO_AGE)))
cat(paste0('Median age of recipients: ',round(quantile(unique(tmp$TO_AGE),probs=0.5),0),
           ' [',round(quantile(unique(tmp$TO_AGE),probs=0.25),0),'-',
           round(quantile(unique(tmp$TO_AGE),probs=0.75),0),']'))
cat(paste0('Median time elapsed (min-max): ',round(quantile(pairs$TIME_ELAPSED,probs=0.5),2),
           ' [',round(min(pairs$TIME_ELAPSED),2),'-',
           round(max(pairs$TIME_ELAPSED),2),']'))
cat(paste0('Median patristic distance (min-max): ',round(quantile(pairs$GEN_DIST,probs=0.5),4),
           ' [',round(min(pairs$GEN_DIST),4),'-',
           round(max(pairs$GEN_DIST),4),']'))

# Boxplot of ages of probable transmitters among all clusters per stage
p <- ggplot(pairs)+geom_boxplot(aes(x=TRANS_STAGE, y=FROM_AGE,fill=as.factor(TRANS_STAGE)))+
  labs(x="\n Stage of the probable transmitter at the \n estimated time of infection", y="Age of the transmitter at the estimated \n time of infection \n")+
  theme_bw()+theme(legend.position = "none")
p
ggsave(file=file.path(out.dir,'Boxplot_age_stage_probable_pairs.pdf'), p, w=6, h=5)

p <- ggplot(pairs,aes(x='All',fill=FROM_BPLACE)) + geom_bar() +
  labs(#x="\n Birth place of recipient at the \n estimated time of infection",
       y="Number of probable transmitters \n",
       fill='Birth place of \nprobable transmitter') +
  theme_bw()+theme(legend.position = "right")
p
ggsave(file=file.path(out.dir,'Boxplot_bplace_transmitter_probable_pairs.pdf'), p, w=6, h=5)

# plot of bplace of probable transmitters among all clusters per stage
p <- ggplot(pairs,aes(x=TRANS_STAGE,fill=FROM_BPLACE)) + geom_bar() +
  labs(x="\n Stage of the probable transmitter at the \n estimated time of infection",
       y="Number of probable transmitters \n",
       fill='Birth place of \nprobable transmitter')+
  theme_bw()+theme(legend.position = "right")
p
ggsave(file=file.path(out.dir,'Boxplot_bplace_stage_probable_pairs.pdf'), p, w=6, h=5)

# plot of bplace of probable transmitters among all clusters per stage
p <- ggplot(pairs,aes(x=TO_BPLACE,fill=FROM_BPLACE)) + geom_bar() +
  labs(x="\n Birth place of recipient at the \n estimated time of infection",
       y="Number of probable transmitters \n",
       fill='Birth place of \nprobable transmitter') +
  theme_bw()+theme(legend.position = "right")
p
ggsave(file=file.path(out.dir,'Boxplot_bplace_to_from_probable_pairs.pdf'), p, w=6, h=5)

tmp <- pairs[, list(N=length(unique(FROM_SEQUENCE_ID))),by='TO_SEQUENCE_ID']

cat(paste0("number of recipients (true pairs) = ", length(unique(pairs$TO_SEQUENCE_ID))))
cat(paste0("mean number of sources per recipient = ", mean(tmp$N)))
cat(paste0("median number of sources per recipient = ", median(tmp$N)))
cat(paste0("number of total pairs = ", nrow(pairs)))
cat(paste0("ratio pairs:non-pairs = ", round(length(unique(pairs$TO_SEQUENCE_ID))/(nrow(pairs)-length(unique(pairs$TO_SEQUENCE_ID))),2)))
cat(paste0("proportion of true pairs = ", round(length(unique(pairs$TO_SEQUENCE_ID))/nrow(pairs)*100),'%'))

write.csv(data.table(table(pairs$FROM_COUNTRY,pairs$FROM_BPLACE)),file=file.path(out.dir, 'from_bplace_includedpairs.csv'))
write.csv(data.table(table(pairs$TO_COUNTRY,pairs$TO_BPLACE)),file=file.path(out.dir, 'to_bplace_includedpairs.csv'))

# save anonymised pairs ----

anon <- subset(pairs,select=c('PAIR_ID','FROM_SEQUENCE_ID','TO_SEQUENCE_ID','GEN_DIST','TIME_ELAPSED','FROM_AGE','TO_AGE','FROM_AGE_GP','TO_AGE_GP'))
write.csv(anon,file=file.path('data_Ams',analysis,paste0(args$trsm,"_anonymised_pairs.csv")))


# summary plots of pairs ----
# plot birthplaces of sources by time period
pairs[YEAR_OF_INF_EST<2016, PERIOD:= '2010-2015']
pairs[YEAR_OF_INF_EST>=2016, PERIOD:= '2016-2021']
tmp <- pairs[, list(N=length(FROM_SEQUENCE_ID)),by=c('PERIOD','FROM_BPLACE')]
tmp[, FROM_BPLACE_LAB:= factor(FROM_BPLACE,levels=c('Overall','NL','G1','G2','G3','Other'),
                          labels=c('Overall','Netherlands','W.Europe,\nN.America,\nOceania','E. & C\nEurope', 'Caribbean &\nS. America','Other'))]
tmp[, TO_BPLACE:= 'Overall']
g1 <- ggplot(tmp) + geom_bar(aes(x=TO_BPLACE,y=N/sum(N),fill=FROM_BPLACE_LAB),stat="identity",position=position_dodge(width=0.9)) +
  theme_bw() + theme(legend.position = "none",
                     axis.text.x = element_text(angle=40, vjust = 0.5)) +
  coord_cartesian(ylim = c(0,1)) +
  scale_y_continuous(labels=scales::label_percent(accuracy = 1L)) +
  labs(x="Birth place of sources",y="",fill="") + scale_fill_npg()
ggsave(file = file.path(out.dir,'birthplaces_included_sources.png'), g1, w = 5, h = 5)


# stratify by bplace of recipient
tmp <- pairs[, list(N=length(FROM_SEQUENCE_ID)),by=c('PERIOD','FROM_BPLACE','TO_BPLACE')]
tmp <- tmp[, list(FROM_BPLACE=FROM_BPLACE,pct=N/sum(N)),by=c('PERIOD','TO_BPLACE')]
tmp[, FROM_BPLACE_LAB:= factor(FROM_BPLACE,levels=c('Overall','NL','G1','G2','G3','Other'),
                               labels=c('Overall','Netherlands','W.Europe,\nN.America,\nOceania','E. & C\nEurope', 'Caribbean &\nS. America','Other'))]
tmp[, TO_BPLACE_LAB:= factor(TO_BPLACE,levels=c('Overall','NL','G1','G2','G3','Other'),
                               labels=c('Overall','Netherlands','W.Europe,\nN.America,\nOceania','E. & C\nEurope', 'Caribbean &\nS. America','Other'))]
g2 <- ggplot(tmp) + geom_bar(aes(x=TO_BPLACE_LAB,y=pct,fill=FROM_BPLACE_LAB),stat="identity",position=position_dodge(width=0.9)) +
  theme_bw() + theme(legend.position = "none",
                     axis.text.x = element_text(angle=40, vjust = 0.5)) +
  coord_cartesian(ylim = c(0,1)) +
  scale_y_continuous(labels=scales::label_percent(accuracy = 1L))  +
  labs(x="Birth place of sources",y="",fill='Birth place of\nrecipient') + scale_fill_npg()
ggsave(file = file.path(out.dir,'birthplaces_included_sources_stratified.png'), g2, w = 5, h = 5)

legend_t <- cowplot::get_legend(g2 + theme(legend.position = "bottom"))

g <- ggarrange(g1 + rremove("xlab"),g2+ rremove("xlab"),ncol=2,widths=c(0.35,0.65),align='hv')
g <- annotate_figure(g, bottom = text_grob("Birth place of recipient",size=18))
g_bplace <- ggarrange(g, legend_t,ncol=1,heights=c(0.8,0.2))
ggsave(file = file.path(out.dir,'birthplaces_included_sources_all_stratified.png'), g_bplace, w = 9, h = 5)
ggsave(file = file.path(out.dir,'birthplaces_included_sources_all_stratified.pdf'), g_bplace, w = 9, h = 5)

# time periods (dutch-born, foreign-born)
pairs[, FROM_BPLACE2:= 'Foreign-born']
pairs[FROM_BPLACE=='NL', FROM_BPLACE2:= 'Dutch-born']
pairs[, TO_BPLACE2:= 'Foreign-born']
pairs[TO_BPLACE=='NL', TO_BPLACE2:= 'Dutch-born']

tmp <- pairs[, list(N=length(FROM_SEQUENCE_ID)),by=c('PERIOD','FROM_BPLACE2')]
tmp <- tmp[, list(FROM_BPLACE2=FROM_BPLACE2,pct=N/sum(N)),by=c('PERIOD')]
pal <- pal_npg("nrc")(4)[c(1,3,4)]

g1 <- ggplot(tmp) + geom_bar(aes(x=FROM_BPLACE2,y=pct,fill=FROM_BPLACE2),stat="identity",width = 1) +
  theme_bw() + theme(legend.position = "none",
                     axis.text.x = element_text(angle=40, vjust = 0.5),
                     strip.background=element_blank()) +
  facet_wrap(PERIOD~.) +
  coord_cartesian(ylim = c(0,1)) +
  scale_y_continuous(labels=scales::label_percent(accuracy = 1L)) +
  labs(x="Birth place of sources",y="",fill="") +
  scale_fill_manual(name="Birth place of\nlikely transmitter",values = c('Dutch-born'=pal[2],'Foreign-born'=pal[3]))
ggsave(file = file.path(out.dir,'birthplaces_time_included_sources.png'), g1, w = 5, h = 5)


# stratify by bplace of recipient
tmp <- pairs[, list(N=length(FROM_SEQUENCE_ID)),by=c('PERIOD','FROM_BPLACE2','TO_BPLACE2')]
tmp <- tmp[, list(FROM_BPLACE2=FROM_BPLACE2,pct=N/sum(N)),by=c('PERIOD','TO_BPLACE2')]
g2 <- ggplot(tmp) + geom_bar(aes(x=TO_BPLACE2,y=pct,fill=FROM_BPLACE2),stat="identity",position=position_dodge(width=0.9)) +
  facet_wrap(PERIOD~.) +
  theme_bw() + theme(legend.position = "bottom",
                     axis.text.x = element_text(angle=40, vjust = 0.5),
                     strip.background=element_blank()) +
  coord_cartesian(ylim = c(0,1)) +
  scale_y_continuous(labels=scales::label_percent(accuracy = 1L))  +
  labs(x="Birth place of recipient",y="",fill='') +
  scale_fill_manual(name="Birth place of\nlikely transmitter",values = c('Dutch-born'=pal[2],'Foreign-born'=pal[3]))
ggsave(file = file.path(out.dir,'birthplaces_time_included_sources_stratified.png'), g2, w = 5, h = 5)

g_time <- ggarrange(g1 + theme(legend.position="none") +  labs(x="",y="",fill="") , g2, ncol=1,heights=c(0.45,0.55),align='hv')

ggsave(file = file.path(out.dir,'birthplaces_time_included_sources_all_stratified.png'), g_time, w = 6, h = 6)
ggsave(file = file.path(out.dir,'birthplaces_time_included_sources_all_stratified.pdf'), g_time, w = 6, h = 6)

g <- ggarrange(g_bplace,g_time,ncol=1,heights=c(0.4,0.6),align='hv',
               labels='AUTO',font.label=list(size=24),vjust=0.8,hjust=0.1)
ggsave(file = file.path(out.dir,'birthplaces_time_included_sources_all_stratified_PANEL.png'), g, w = 8, h = 12)
ggsave(file = file.path(out.dir,'birthplaces_time_included_sources_all_stratified_PANEL.pdf'), g, w = 8, h = 12)


# plot ages
tmp <- pairs[, list(N=length(FROM_SEQUENCE_ID)),by=c('FROM_AGE_GP')]
tmp <- tmp[, list(FROM_AGE_GP=FROM_AGE_GP,pct=N/sum(N))]
pal <- pal_npg("nrc")(4)[c(1,3,4)]

tmp[, TO_AGE_GP:= 'Overall']
g1 <- ggplot(tmp) + geom_bar(aes(x=TO_AGE_GP,y=pct,fill=FROM_AGE_GP),stat="identity",width = 1,position=position_dodge(width=0.9)) +
  theme_bw() + theme(legend.position = "none",
                     axis.text.x = element_text(angle=40, vjust = 0.5),
                     strip.background=element_blank()) +
  coord_cartesian(ylim = c(0,1)) +
  scale_y_continuous(labels=scales::label_percent(accuracy = 1L)) +
  labs(x="Birth place of sources",y="",fill="") +
  scale_fill_npg()
ggsave(file = file.path(out.dir,'ages_included_sources.png'), g1, w = 5, h = 5)

tmp <- pairs[, list(N=length(FROM_SEQUENCE_ID)),by=c('FROM_AGE_GP','TO_AGE_GP')]
tmp <- tmp[, list(FROM_AGE_GP=FROM_AGE_GP,pct=N/sum(N)),by='TO_AGE_GP']

g2 <- ggplot(tmp) + geom_bar(aes(x=TO_AGE_GP,y=pct,fill=FROM_AGE_GP),stat="identity",width = 1,position=position_dodge(width=0.9)) +
  theme_bw() + theme(legend.position = "bottom",
                     axis.text.x = element_text(angle=40, vjust = 0.5),
                     strip.background=element_blank()) +
  coord_cartesian(ylim = c(0,1)) +
  scale_y_continuous(labels=scales::label_percent(accuracy = 1L)) +
  labs(x="Birth place of recipients",y="",fill="Birth place of sources") +
  scale_fill_npg()
ggsave(file = file.path(out.dir,'ages_stratified_included_sources.png'), g2, w = 5, h = 5)

legend_t <- cowplot::get_legend(g2 + theme(legend.position = "bottom"))

g <- ggarrange(g1 + rremove("xlab"),g2+ rremove("xlab") + theme(legend.position="none"),ncol=2,widths=c(0.35,0.65),align='hv')
g <- annotate_figure(g, bottom = text_grob("Birth place of recipients",size=18))
g_bplace <- ggarrange(g, legend_t,ncol=1,heights=c(0.8,0.2))
ggsave(file = file.path(out.dir,'ages_stratified_included_sources_all.png'), g_bplace, w = 9, h = 5)
ggsave(file = file.path(out.dir,'ages_stratified_included_sources_all.pdf'), g_bplace, w = 9, h = 5)


# Check and plot how many points are in the signal cone
pairs$TRANS_STAGE <- factor(pairs$TRANS_STAGE,levels=c('UNDIAGNOSED','DIAGNOSED','ART'),labels=c('Undiagnosed','Diagnosed','On ART (not suppressed)'))


# plot distances

## load estimated molecular clock ----
cat(" \n ------------- \n Load quantiles from fitted molecular clock model \n ------------ \n")

# get quantiles of mean/median distances at prediction points
dpr <- data.table(d_TSeqT = seq(0.1,18,0.1))
dpr[, id := seq_len(nrow(dpr))]

cm <- readRDS(file.path(clock_model,'clock_model_gamma_hier_220315-stan_fit.rds'))
pd <- cm$draws(inc_warmup = FALSE)
po <- list()
tmp <- pd[,,which(grepl('log_alpha1$',dimnames(pd)[[3]]))]
po$log_alpha1 <- as.vector(tmp)
tmp <- pd[,,which(grepl('log_alpha1_pair_sd',dimnames(pd)[[3]]))]
po$log_alpha1_pair_sd <- as.vector(tmp)
tmp <- pd[,,which(grepl('log_phi$',dimnames(pd)[[3]]))]
po$log_phi <- as.vector(tmp)
tmp <- pd[,,which(grepl('log_phi_pair_sd',dimnames(pd)[[3]]))]
po$log_phi_pair_sd <- as.vector(tmp)

## predict distances using posterior medians only
dpo <- data.table(iter = seq_along(po$log_alpha1),
                  log_alpha1 = median(po$log_alpha1),
                  log_alpha1_pair_sd = median(po$log_alpha1_pair_sd),
                  log_phi = median(po$log_phi),
                  log_phi_pair_sd = median(po$log_phi_pair_sd),
                  DUMMY = 1
)
dpo[, log_alpha1_pair := rnorm(nrow(dpo), 0, log_alpha1_pair_sd)]
dpo[, log_phi_pair := rnorm(nrow(dpo), 0, log_phi_pair_sd)]
dpo[, er := exp(log_alpha1 + log_alpha1_pair)]
dpo[, beta := exp( -(log_phi + log_phi_pair))]

dpr[, DUMMY := 1]
dpo <- merge(dpo, dpr, by = 'DUMMY', allow.cartesian = TRUE)
set(dpr, NULL, 'DUMMY', NULL)
dpo[, y_pr := rgamma(nrow(dpo), shape = er * d_TSeqT * beta, rate = beta)]

dps <- dpo[,
           list(
             p = quantile(y_pr, prob = c(0.025, seq(0.1, .9, 0.1), 0.975)),
             qlabel = paste0('q',c(0.025, seq(0.1, .9, 0.1), 0.975)*100)
           ),
           by = c('d_TSeqT')
]
dps_clock <- dcast.data.table(dps, d_TSeqT ~ qlabel, value.var = 'p')
saveRDS(dps_clock,file = file.path(out.dir,'clock_quantiles.rds'))


pal_3 <- pal_npg("nrc")(4)[c(1,3,4)]
p.palette <- RColorBrewer::brewer.pal(5,'Oranges')
p.alpha <- 0.7

p <- ggplot(dps_clock, aes(x = d_TSeqT))+
  geom_ribbon(data = dps_clock, aes(ymin = q2.5, ymax = q10), fill = p.palette[1], alpha = p.alpha) +
  geom_ribbon(data = dps_clock, aes(ymin = q90, ymax = q97.5), fill = p.palette[1], alpha = p.alpha) +
  geom_ribbon(data = dps_clock, aes(ymin = q10, ymax = q20), fill = p.palette[2], alpha = p.alpha) +
  geom_ribbon(data = dps_clock, aes(ymin = q80, ymax = q90), fill = p.palette[2], alpha = p.alpha) +
  geom_ribbon(data = dps_clock, aes(ymin = q20, ymax = q30), fill = p.palette[3], alpha = p.alpha) +
  geom_ribbon(data = dps_clock, aes(ymin = q70, ymax = q80), fill = p.palette[3], alpha = p.alpha) +
  geom_ribbon(data = dps_clock, aes(ymin = q30, ymax = q40), fill = p.palette[4], alpha = p.alpha) +
  geom_ribbon(data = dps_clock, aes(ymin = q60, ymax = q70), fill = p.palette[4], alpha = p.alpha) +
  geom_ribbon(data = dps_clock, aes(ymin = q40, ymax = q60), fill = p.palette[5], alpha = p.alpha) +
  geom_point(data=pairs,aes(shape=TRANS_STAGE,colour=FROM_BPLACE,x=TIME_ELAPSED,y=GEN_DIST))+
  #scale_y_continuous(breaks=seq(0,0.06,0.01),labels=scales::label_percent(accuracy = 1L)) +
  theme_bw(base_size=26)+
  scale_colour_npg()+
  #scale_colour_manual(values=pal_3,name="Stage of probable \n transmitter") +
  labs(x='\n Time elapsed (in years)',y='\n Genetic distance of the pair \n',shape='Stage of probable \n transmitter')+
  theme(legend.position='bottom',
        strip.text=element_text(size=18),
        legend.title=element_text(size=22))+
  coord_cartesian(expand = c(0,0),xlim=c(0,18)) +
  guides(size = "none")
p
ggsave(file=file.path(out.dir,'Ams_all_pairs_genetic_distances_time_elapsed.png'), p, w=12, h=10)

p <- ggplot(dps_clock, aes(x = d_TSeqT))+
  geom_ribbon(data = dps_clock, aes(ymin = q2.5, ymax = q10), fill = p.palette[1], alpha = p.alpha) +
  geom_ribbon(data = dps_clock, aes(ymin = q90, ymax = q97.5), fill = p.palette[1], alpha = p.alpha) +
  geom_ribbon(data = dps_clock, aes(ymin = q10, ymax = q20), fill = p.palette[2], alpha = p.alpha) +
  geom_ribbon(data = dps_clock, aes(ymin = q80, ymax = q90), fill = p.palette[2], alpha = p.alpha) +
  geom_ribbon(data = dps_clock, aes(ymin = q20, ymax = q30), fill = p.palette[3], alpha = p.alpha) +
  geom_ribbon(data = dps_clock, aes(ymin = q70, ymax = q80), fill = p.palette[3], alpha = p.alpha) +
  geom_ribbon(data = dps_clock, aes(ymin = q30, ymax = q40), fill = p.palette[4], alpha = p.alpha) +
  geom_ribbon(data = dps_clock, aes(ymin = q60, ymax = q70), fill = p.palette[4], alpha = p.alpha) +
  geom_ribbon(data = dps_clock, aes(ymin = q40, ymax = q60), fill = p.palette[5], alpha = p.alpha) +
  geom_point(data=pairs,aes(color=TRANS_STAGE,x=TIME_ELAPSED,y=GEN_DIST))+
  scale_y_continuous(breaks=seq(0,0.06,0.01),labels=scales::label_percent(accuracy = 1L)) +
  theme_bw(base_size=26)+
  scale_colour_manual(values=pal_3,name="Stage of probable \n transmitter") +
  labs(x='\n Time elapsed (in years)',y='\n Genetic distance of the pair \n')+
  theme(legend.position='bottom',
        strip.text=element_text(size=18),
        legend.title=element_text(size=22))+
  coord_cartesian(expand = c(0,0),ylim=c(0,0.06),xlim=c(0,10)) +
  guides(size = "none")
p
ggsave(file=file.path(out.dir,'Ams_all_pairs_genetic_distances_time_elapsed_max5pct.png'), p, w=12, h=10)

pairs <- readRDS(file=file.path(out.dir, paste0(args$trsm,'_pairs.rds')))
pairs$TRANS_STAGE <- factor(pairs$TRANS_STAGE,levels=c('UNDIAGNOSED','DIAGNOSED','ART'),labels=c('Undiagnosed','Diagnosed','On ART (not suppressed)'))

p <- ggplot(dps_clock, aes(x = d_TSeqT))+
  geom_ribbon(data = dps_clock, aes(ymin = q2.5, ymax = q10), fill = p.palette[1], alpha = p.alpha) +
  geom_ribbon(data = dps_clock, aes(ymin = q90, ymax = q97.5), fill = p.palette[1], alpha = p.alpha) +
  geom_ribbon(data = dps_clock, aes(ymin = q10, ymax = q20), fill = p.palette[2], alpha = p.alpha) +
  geom_ribbon(data = dps_clock, aes(ymin = q80, ymax = q90), fill = p.palette[2], alpha = p.alpha) +
  geom_ribbon(data = dps_clock, aes(ymin = q20, ymax = q30), fill = p.palette[3], alpha = p.alpha) +
  geom_ribbon(data = dps_clock, aes(ymin = q70, ymax = q80), fill = p.palette[3], alpha = p.alpha) +
  geom_ribbon(data = dps_clock, aes(ymin = q30, ymax = q40), fill = p.palette[4], alpha = p.alpha) +
  geom_ribbon(data = dps_clock, aes(ymin = q60, ymax = q70), fill = p.palette[4], alpha = p.alpha) +
  geom_ribbon(data = dps_clock, aes(ymin = q40, ymax = q60), fill = p.palette[5], alpha = p.alpha) +
  geom_point(data=pairs,aes(color=TRANS_STAGE,x=TIME_ELAPSED,y=GEN_DIST,shape=FROM_BPLACE))+
  scale_y_continuous(breaks=seq(0,0.06,0.01),labels=scales::label_percent(accuracy = 1L)) +
  theme_bw(base_size=26)+
  scale_colour_manual(values=pal_3,name="Stage of probable \n transmitter") +
  labs(x='\n Time elapsed (in years)',y='\n Genetic distance of the pair \n',shape='Birth place of\nprobable transmitter')+
  theme(legend.position='bottom',
        strip.text=element_text(size=18),
        legend.title=element_text(size=22))+
  coord_cartesian(expand = c(0,0),ylim=c(0,0.06),xlim=c(0,10)) +
  guides(size = "none", colour = guide_legend(nrow = 2,title.position='top'),
                                  shape = guide_legend(nrow = 2,title.position='top'), by.col=T)

p
ggsave(file=file.path(out.dir,'Ams_all_pairs_genetic_distances_time_elapsed_max5pct_bplace.png'), p, w=12, h=10)

pairs <- subset(pairs,select = c("PAIR_ID","TIME_ELAPSED","GEN_DIST","TRANS_STAGE","FROM_AGE","TO_AGE"))

# Load data from fitted gamma molecular clock on confirmed pairs from Belgium
set(pairs,NULL,"TIME_ELAPSED",as.factor(round(pairs$TIME_ELAPSED,1)))
#trm.pol.p3 <- readRDS("data_other/signal_cone_Belgium.rds")
trm.pol.p3 <- dps_clock
setnames(trm.pol.p3,"d_TSeqT","TIME_ELAPSED")
set(trm.pol.p3,NULL,"TIME_ELAPSED",as.factor(trm.pol.p3$TIME_ELAPSED))
tmp <- merge(pairs,trm.pol.p3,by = "TIME_ELAPSED")
#tmp[,IN_CONE:= GEN_DIST>=tmp$q1 & GEN_DIST<=tmp$q99]
tmp[,IN_CONE:= GEN_DIST>=tmp$q2.5 & GEN_DIST<=tmp$q97.5]
tmp1 <- as.data.table(tmp  %>% summarise(Pairs=n()))
tmp2 <- as.data.table(tmp  %>% filter(IN_CONE==TRUE) %>% summarise(Pairs_in_signal_cone=n()))
tmp <- cbind(tmp1,tmp2)
tmp <- melt(tmp,measure.vars = c("Pairs","Pairs_in_signal_cone"))
p <- ggplot(tmp,aes(x=as.factor(variable),y=value,fill=as.factor(variable)))+geom_bar(stat="identity",width = 0.2)+theme_bw()+theme(legend.position = "none")+
  labs(x="",y="Total number of probable transmission pairs \n and corresponding pairs in signal cone \n",fill="")
p
ggsave(file=file.path(out.dir,'Pairs_in_cone.pdf'), p, w=5, h=4)


