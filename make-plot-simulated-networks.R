require(data.table)
require(ggplot2)
require(ggsci)
require(ggpubr)
require(igraph)
require(fitdistrplus)
require(GGally)
require(condMVNorm)

args <- list(
  # stanModelFile = 'mm_sigHierG_bgUnif_piReg_221123', # regression on diag stage
  source_dir = '/Users/bethanyastor/Library/CloudStorage/OneDrive-ImperialCollegeLondon/Research Project/source.attribution.bmm',
  clock_model = '/Users/bethanyastor/Library/CloudStorage/OneDrive-ImperialCollegeLondon/Research Project/source.attribution.bmm/molecular_clock',
  hmc_stepsize = 0.02,
  hmc_num_samples = 15,
  hmc_num_warmup = 10,
  seed = 42,
  chain = 1,
  scenario = 1,
  reps = 1,
  npairs = 500,#500
  p_pairs = 0.5,  # 0.5 - 0.01 for non-network sims
  p_trm_undiag = 1, # proportion of true transmitters undiagnosed
  p_nontrm_undiag = 0, # proportion of non-transmitters undiagnosed
  simulate_data = T,
  job_tag = 'simulations_network_500truepairs_prop_subsample100pct',
  local=1,
  networks=T,
  pct_trsm_diag=0,
  pct_trsm_mig=1, # set==1 if all transmitters migrants, 0.8 if 80% are migrants
  sim_bplace='fixed', # or after
  bg = 'unif',
  no_corr=0,
  excl_improb=T,
  all_pairs=F,
  late=F,
  network_data_file='V1.2_patch0_Rand10_Run1_PCseed0_0_transmission_events.csv', # phylogenetic_ZA_seed2_run1.csv
  clu_id=372,
  id_case=32950
)

out.dir <- '/Users/bethanyastor/Library/CloudStorage/OneDrive-ImperialCollegeLondon/Research Project/source.attribution.bmm/figures'
reps <- 1

source('R/functions_plotting.R')

networks <- file.path(args$source_dir,'data_other','transmission_network',args$network_data_file)
infile.inftime <- file.path(args$source_dir,'data_Ams','analysis_220713','roadmap_TSI_estimates_MSM.csv')

seed <- 1234567*reps # set seed for each chain, changing for each replicate

set.seed(seed)

args$reps
args$scenario
cluster_sample=NA
npairs=args$npairs
p_pairs=args$p_pairs
all_pairs=args$all_pairs
bg = args$bg
sim_bplace=args$sim_bplace
pct_trsm_diag=args$pct_trsm_diag
pct_trsm_mig=args$pct_trsm_mig
no_corr=args$no_corr
excl_improb=args$excl_improb

## load estimated molecular clock ----
cat(" \n ------------- \n Load quantiles from fitted molecular clock model \n ------------ \n")

cm <- readRDS(file.path(args$clock_model,'clock_model_gamma_hier_220315-stan_fit.rds'))
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

# specify median posterior of parameters to simulate with (from convergence.csv for gamma clock model)
log_alpha1 <- median(po$log_alpha1) # 0.004653835
log_alpha1_pair_sd <- median(po$log_alpha1_pair_sd)
log_phi <- median(po$log_phi) # 0.00414958
log_phi_pair_sd <- median(po$log_phi_pair_sd)

dpr <- data.table(d_TSeqT = seq(0.1,18,0.1))
dpr[, id := seq_len(nrow(dpr))]
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
dps <- dcast.data.table(dps, d_TSeqT ~ qlabel, value.var = 'p')

## load data ----
dt <- data.table(read.csv(networks,header=T))
setnames(dt,c('TimeOfEvent','InfectedAgeAtEvent','InfectorAgeAtEvent','InfectorARTStatusAtEvent'),c('TimeOfInfection','Infected_age','Infector_age','PartnerARTStatus'))

## define transmission stage of true pairs ----
dt[, TRANS_STAGE:= factor(as.character(PartnerARTStatus),levels=c(-1,0,1,3,5),labels=c('Undiagnosed','Diagnosed','Diagnosed','Diagnosed','Diagnosed'))]
dt[, TRANSMISSION_PAIR:='YES']

dt$IdInfector <- as.character(dt$IdInfector)
dt$IdInfected <- as.character(dt$IdInfected)

bdate <- unique(subset(dt,select=c('IdInfected','TimeOfInfection','Infected_age')))
bdate[, bdate:= TimeOfInfection - Infected_age]

# identify transmission clusters ----
tmp <- subset(dt, select=c(IdInfector, IdInfected))
tmp <- graph.data.frame(tmp, directed=FALSE, vertices=NULL)
rtc <- data.table(ID=V(tmp)$name, CLU=clusters(tmp, mode="weak")$membership)
tmp2 <- rtc[, list(CLU_SIZE=length(ID)), by='CLU'] # max cluster size = 59

pal <- pal_npg("nrc")(9)
p <- ggplot(tmp2) + geom_histogram(aes(x=CLU_SIZE),fill=pal[1],binwidth=2) +
  theme_bw(base_size=26) + labs(x='Cluster size',y='count')

# subsample clusters to achieve a specified number of transmission events
if(!is.na(npairs)){
  
  dt <- merge(dt,rtc,by.x='IdInfected',by.y='ID',all.x=T)
  # re-do to keep most recent npairs events instead
  tab <- dt[order(-TimeOfInfection),]
  ids <- unique(tab[1:500,IdInfected])
  clu <- unique(dt[IdInfected %in% ids, CLU])
}

# get dates
f <- function(x) {
  year <- floor(x);
  month <- floor((x - year) * 12);
  day <- ((x - year) * 12 - month) * 30.42;
  return(sprintf("%i years, %i months, %3.2f days", year, month, day))
}
cat(paste0('from ',f(min(tab[1:500,TimeOfInfection]))))
cat(paste0('to ',f(max(tab[1:500,TimeOfInfection]))))

# keep only clusters which contain one of the 500 most recent incident cases
dt <- subset(dt,CLU %in% clu)
dt[, case:= ifelse(IdInfected %in% ids,1,0)]

## plot network ----
sg_all <- copy(dt)
# sample 50 cases and select their networks
plot_ids <- sample(sg_all$IdInfected,50)
clu_ids <- sg_all[IdInfected %in% plot_ids,'CLU']
sg <- subset(sg_all,select=c('IdInfector','IdInfected','TRANSMISSION_PAIR','case','CLU'),CLU %in% clu_ids$CLU)
network <- graph_from_data_frame(d=sg, directed=T)

# change colour of vertices
col = pal_npg("nrc")(5)
V(network)$cluster<-col[2]
V(network)$cluster[which(V(network)$name %in% sg$IdInfected[sg$case==1])]<-col[1]
V(network)$color=V(network)$cluster

# colour edges
E(network)$color <- 'darkgrey'
E(network)$color[sg$TRANSMISSION_PAIR=='YES'] <- col[1]


colrs <- c(col[1],col[2])
l <- layout_nicely(network)


# only colour for one case
id_src <- dt$IdInfector[dt$IdInfected==args$id_case]
V(network)$cluster<- 'grey50' # mark potential sources
V(network)$cluster[which(V(network)$name %in% sg$IdInfected[sg$CLU==args$clu_id])]<-col[2] # mark phylo sources
V(network)$cluster[which(V(network)$name %in% sg$IdInfector[sg$CLU==args$clu_id])]<-col[2] # mark phylo sources (ensure first case is coloured too)
V(network)$cluster[which(V(network)$name==args$id_case)]<-col[1] # mark case
V(network)$cluster[which(V(network)$name==id_src)]<-col[4] # mark true sources
V(network)$color=V(network)$cluster

V(network)$size <- 1.5
V(network)$size[which(V(network)$name %in% sg$IdInfected[sg$CLU==args$clu_id])] <- 3
V(network)$size[which(V(network)$name %in% sg$IdInfector[sg$CLU==args$clu_id])] <- 3

pal <- c(col[1],'grey50',col[2],col[4])
pdf(file=file.path(out.dir,'simulated_transmission_chains_500cases_ggraph_layout_sample_1case_centred.pdf'),h=30,w=30)
plot(network, layout=layout_nicely, vertex.size=V(network)$size, edge.color = "grey50",
     edge.arrow.size=1, edge.size=3,
     vertex.label=NA)
legend(x=-1, y=-0.93, c("Incident case","Potential sources",
                        "Possible sources in same transmission chain", "True source"), pch=21,
       pt.bg=pal, pt.cex=4, cex=3.5, bty="n", ncol=1)
dev.off()


## plot histogram ----

set.seed(seed)

## load data ----
dt <- data.table(read.csv(networks,header=T))
setnames(dt,c('TimeOfEvent','InfectedAgeAtEvent','InfectorAgeAtEvent','InfectorARTStatusAtEvent'),c('TimeOfInfection','Infected_age','Infector_age','PartnerARTStatus'))

dt$IdInfector <- as.character(dt$IdInfector)
dt$IdInfected <- as.character(dt$IdInfected)

bdate <- unique(subset(dt,select=c('IdInfected','TimeOfInfection','Infected_age')))
bdate[, bdate:= TimeOfInfection - Infected_age]

### identify transmission clusters ----
tmp <- subset(dt, select=c(IdInfector, IdInfected))
tmp <- graph.data.frame(tmp, directed=FALSE, vertices=NULL)
rtc <- data.table(ID=V(tmp)$name, CLU=clusters(tmp, mode="weak")$membership)
tmp2 <- rtc[, list(CLU_SIZE=length(ID)), by='CLU'] # max cluster size = 59

# subsample clusters to achieve a specified number of transmission events
if(!is.na(npairs)){
  
  dt <- merge(dt,rtc,by.x='IdInfected',by.y='ID',all.x=T)
  
  # re-do to keep most recent npairs events instead
  tab <- dt[order(-TimeOfInfection),]
  ids <- unique(tab[1:500,IdInfected])
  real_pairs <- unique(tab[1:500,c('IdInfector','IdInfected')])
  real_pairs[, TRANSMISSION_PAIR:= 'YES']
}

### find all possible-pairs ----
rtc2 <- copy(rtc)
setnames(rtc,'ID','FROM_ID')
setnames(rtc2,'ID','TO_ID')
sg <- expand.grid(FROM_ID=rtc$FROM_ID,TO_ID=rtc2$TO_ID)
sg <- subset(sg, TO_ID %in% ids)
sg <- subset(sg,FROM_ID!=TO_ID)

tmp <- copy(rtc)
setnames(tmp, 'CLU','FROM_CLU')
sg <- merge(sg,tmp,by='FROM_ID',all.x=T)
tmp <- copy(rtc2)
setnames(tmp, 'CLU','TO_CLU')
sg <- merge(sg,tmp,by='TO_ID',all.x=T)
sg <- data.table(sg)

sg <- merge(sg,real_pairs,by.x=c('FROM_ID','TO_ID'),by.y=c('IdInfector','IdInfected'),all.x=T)

# first merge with data on infected (time of infection)
sg <- merge(sg,subset(dt,select=c('IdInfected','TimeOfInfection','TimestepOfInfection','InfectedSPVL','Infected_age','Infector_age')),
            by.x=c('TO_ID'),
            by.y=c('IdInfected'),all.x=T) # just keep the pairs in the subsetted data (otherwise set all=T)

# add date of infection of infector
do <- unique(subset(dt,select=c('IdInfected','TimeOfInfection')))
setnames(do,c('IdInfected','TimeOfInfection'),c('FROM_ID','FROM_TimeOfInfection'))
sg <- merge(sg,do,by=c('FROM_ID'),all=T)
# remove rows with no infected (shouldn't be any)
sg <- subset(sg,!is.na(TO_ID))

# define all non-pairs
sg[is.na(TRANSMISSION_PAIR),TRANSMISSION_PAIR:='NO']

# apply exclusion criteria
# add infection date for those who don't appear in 'TO' column
tmp <- unique(subset(dt,select=c('IdInfector','TimeOfInfection','InfectorTimeElapsedSinceInfection'),!is.na(InfectorTimeElapsedSinceInfection)))
tmp[, FROM_TimeOfInfection_miss:= TimeOfInfection - InfectorTimeElapsedSinceInfection]
tmp <- unique(subset(tmp,select=c('IdInfector','FROM_TimeOfInfection_miss')))
sg <- merge(sg,subset(tmp,select=c('IdInfector','FROM_TimeOfInfection_miss')),by.x='FROM_ID',by.y='IdInfector',all.x=T)
sg[is.na(FROM_TimeOfInfection), FROM_TimeOfInfection:= FROM_TimeOfInfection_miss]
sg[, FROM_TimeOfInfection_miss:= NULL]

if(excl_improb==T){
  # remove pairs where infector was infected after infected
  sg[,EXCLUDE:=ifelse(FROM_TimeOfInfection>TimeOfInfection,1,0)]
}

cat(paste0('Number of pairs where transmitter infected before recipient: ', nrow(sg[EXCLUDE==0]),'\n'))
cat(paste0('Number of sources infected before recipient: ', length(unique(sg$FROM_ID[sg$EXCLUDE==0])),'\n'))

tmp <- subset(sg,select=c('FROM_ID','TO_ID'),EXCLUDE==0)
tmp <- tmp[, list(N_SOURCES=length(unique(FROM_ID))),by='TO_ID']
cat(paste0('Average number of potential sources: ', mean(tmp$N_SOURCES),'\n'))

# exclude those not in same cluster
sg[, PHYLO_EXCLUDE:= ifelse(FROM_CLU!=TO_CLU,1,0)]


### simulate sampling dates ----

# I use the average time from infection to diagnosis to simulate with (assuming most people are sequenced when they are diagnosed)
do <- data.table(read.csv(infile.inftime,header=T))
do <- unique(subset(do,select=c(id,estsctodiagMedian)))

coeff <- fitdist(do$estsctodiagMedian[!is.na(do$estsctodiagMedian)], distr = "weibull", method = "mle")

# get time of infection for all transmitters and recipients
ds <- unique(subset(sg,select=c('TO_ID','TimeOfInfection')))
tmp <- unique(subset(sg,select=c('FROM_ID','FROM_TimeOfInfection')))
setnames(tmp,c('FROM_ID','FROM_TimeOfInfection'),c('TO_ID','TimeOfInfection'))
ds <- unique(rbind(ds,tmp))

# I simulate by adding a rng from weibull dist to the estimated time of infection
ds[, date_s:= TimeOfInfection + rweibull(nrow(ds),shape=coeff$estimate[['shape']],scale=coeff$estimate[['scale']])]

# merge the infector and infected sampling dates
df <- merge(sg,ds,by=c('TO_ID','TimeOfInfection'),all.x=T)
setnames(df,'date_s','TO_DATE_S')
df <- merge(df,subset(ds,select=c('TO_ID','date_s')),by.x=c('FROM_ID'),by.y=c('TO_ID'),all.x=T)
setnames(df,'date_s','FROM_DATE_S')

# calculate time elapsed
df[, TIME_ELAPSED:= abs(TO_DATE_S - TimeOfInfection) + abs(FROM_DATE_S - TimeOfInfection)]

# only keep pairs with time elapsed <16 years
df[TIME_ELAPSED>=16, EXCLUDE:= 1]

cat(paste0('Number of pairs with time elapsed <16yrs: ', nrow(df[df$EXCLUDE==0]),'\n'))
cat(paste0('Number of recipients remaining: ', length(unique(df$TO_ID[df$EXCLUDE==0 & df$TRANSMISSION_PAIR=='YES'])),'\n'))
cat(paste0('% pairs: ', length(unique(df$TO_ID[df$EXCLUDE==0 & df$TRANSMISSION_PAIR=='YES']))/nrow(df[df$EXCLUDE==0])*100,'\n'))

df[FROM_CLU!=TO_CLU, EXCLUDE:= 1]

cat(paste0('Number of phylogenetically possible pairs: ', nrow(df[df$EXCLUDE==0]),'\n'))
cat(paste0('Number of recipients remaining: ', length(unique(df$TO_ID[df$EXCLUDE==0 & df$TRANSMISSION_PAIR=='YES'])),'\n'))
cat(paste0('% pairs: ', length(unique(df$TO_ID[df$EXCLUDE==0 & df$TRANSMISSION_PAIR=='YES']))/nrow(df[df$EXCLUDE==0])*100,'\n'))


### prepare for histograms ----
sg <- copy(df)
tmp <- subset(sg,select=c('FROM_ID','TO_ID'))
tmp2 <- tmp[, list(N_SOURCES=length(unique(FROM_ID))),by='TO_ID']


tmp <- subset(sg,select=c('FROM_ID','TO_ID'), EXCLUDE==0)
tmp3 <- tmp[, list(N_SOURCES=length(unique(FROM_ID))),by='TO_ID']
cat(paste0('Average number of phylogenetically possible sources: ', mean(tmp3$N_SOURCES),'\n'))

# plot cluster sizes
pal <- pal_npg("nrc")(9)
p_hist <- ggplot() +
  geom_histogram(data=tmp3,aes(x=N_SOURCES,fill=pal[1],col=pal[1]),binwidth=1,alpha=0.5) +
  theme_bw(base_size=16) + labs(x='Phylogenetically possible per case',y='Number of incident cases') +
  scale_fill_manual(name="",values = c(pal[1]),
                    labels=c('Phylogenetically possible sources')) +
  scale_linetype_manual(name="",values = c(2),
                        labels=c('Potential sources')) +
  scale_x_continuous(breaks=c(seq(0,60,10)),labels=c(seq(0,60,10))) +
  guides(fill = guide_legend(nrow = 1),
         linetype = guide_legend(nrow = 1),
         colour = 'none') +
  theme(legend.position='bottom', legend.margin=margin(t=-10))
p_hist


# simulate distances ----

# simulate ages
# get age of missing infectors
df <- merge(df,subset(bdate,select=c('IdInfected','bdate')),by.x='FROM_ID',by.y='IdInfected',all.x=T)
df[is.na(Infector_age), Infector_age:= TimeOfInfection - bdate]
# some still missing, so simulate uniformly
df[is.na(Infector_age), Infector_age:= runif(length(which(is.na(Infector_age))), 16, 70)]
setnames(df,c('Infector_age','Infected_age'),c('FROM_AGE','TO_AGE'))

# simulate genetic distances ----
# first simulate for true pairs from gamma dist
df[, log_alpha1_pair := rnorm(nrow(df), 0, log_alpha1_pair_sd)]
df[, log_phi_pair := rnorm(nrow(df), 0, log_phi_pair_sd)]
df[, ER := exp(log_alpha1 + log_alpha1_pair)]
df[, BETA := exp( -(log_phi + log_phi_pair))]

sim_scenario <- as.data.table(df[,c('FROM_ID','TO_ID','TRANSMISSION_PAIR','EXCLUDE','FROM_AGE','TO_AGE','TIME_ELAPSED','ER','BETA')])
sim_scenario[TRANSMISSION_PAIR=="YES", GAMMA_SHAPE:= ER * TIME_ELAPSED * BETA]
sim_scenario[TRANSMISSION_PAIR=="YES", GAMMA_RATE:= BETA]
sim_scenario[TRANSMISSION_PAIR=="YES", GEN_DIST:= rgamma(length(which(sim_scenario$TRANSMISSION_PAIR=="YES")),shape = GAMMA_SHAPE,rate=GAMMA_RATE)]

# Uniform background distribution 
#sim_scenario[TRANSMISSION_PAIR=="NO", GEN_DIST:= runif(length(which(df$TRANSMISSION_PAIR=="NO")),0,0.2)]

# GMM background 
#gmm <- readRDS(file.path(args$source_dir, 'out_Amsterdam/mm_bgGMM_piGP_240711-agegps_sensanalysis_210216_MSM-594613/mm_bgGMM_piGP_240711-agegps_sensanalysis_210216_MSM-GMM_params5clusters.RData'))
gmm <- readRDS(file.path(args$source_dir, 'out_Amsterdam/mm_bgGMM_piGP_240711-agegps_sensanalysis_210216_MSM-1850682-rm-signal/mm_bgGMM_piGP_240711-agegps_sensanalysis_210216_MSM-rm-signal-GMM_params5clusters.RData'))
# assign to cluster 
k <- gmm$nclass
sim_scenario[TRANSMISSION_PAIR=="NO", CLUSTER:=sample(1:k, size=length(which(df$TRANSMISSION_PAIR=="NO")), replace=TRUE, prob=gmm$pi)]
# conditional on TIME_ELAPSED, generate GEN_DIST from mixture model 
sigma <- LTSigma2variance(gmm$LTSigma)
sim_scenario[TRANSMISSION_PAIR=="NO",GEN_DIST:=mapply(CLUSTER, TIME_ELAPSED,FUN=function(cl, t) {rcmvnorm(n=1, mean=gmm$Mu[cl, ], sigma=sigma[,,cl], dependent.ind=2, given.ind=1, X.given=t)})]


# plot distances ----
sim_scenario[, TRSM_PAIR:= 'Competing plausible pair']
sim_scenario[TRANSMISSION_PAIR=='YES', TRSM_PAIR:= 'True pair']
sim_scenario[, TRSM_PAIR:= factor(TRSM_PAIR,levels=c('True pair','Competing plausible pair'))]

# get IQR of quantiles
dpr <- data.table(d_TSeqT = seq(0.00,18,0.01))
dpr[, id := seq_len(nrow(dpr))]
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

iqr <- dpo[,
           list(
             p = quantile(y_pr, prob = c(0.25, 0.75,0.025,0.975)),
             qlabel = paste0('q',c(0.25,0.75,0.025,0.975)*100)
           ),
           by = c('d_TSeqT')
]
iqr <- dcast.data.table(iqr, d_TSeqT ~ qlabel, value.var = 'p')
iqr[, d_TSeqT:= round(d_TSeqT,2)]

sim_scenario[, d_TSeqT:= round(TIME_ELAPSED,2)]
sim_scenario <- merge(sim_scenario,iqr, by='d_TSeqT',all.x=T)

# flag the false positives
sim_scenario[, FALSE_POS:= ifelse(TRANSMISSION_PAIR=='NO' & GEN_DIST>=q2.5 & GEN_DIST<=q97.5,1,0)]
sim_scenario[, cat:= as.factor(ifelse(TRANSMISSION_PAIR=='YES',1,ifelse(FALSE_POS==0,2,3)))]

p.palette <- RColorBrewer::brewer.pal(5,'Oranges')
p.alpha <- 0.7

pal1 <- pal_npg("nrc")(5)

p <- ggplot(data = dps,aes(x=d_TSeqT))+
  geom_ribbon(data = dps, aes(ymin = q2.5, ymax = q10), fill = p.palette[1], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q90, ymax = q97.5), fill = p.palette[1], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q10, ymax = q20), fill = p.palette[2], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q80, ymax = q90), fill = p.palette[2], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q20, ymax = q30), fill = p.palette[3], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q70, ymax = q80), fill = p.palette[3], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q30, ymax = q40), fill = p.palette[4], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q60, ymax = q70), fill = p.palette[4], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q40, ymax = q60), fill = p.palette[5], alpha = p.alpha) +
  geom_point(data=subset(sim_scenario, EXCLUDE==0),aes(color=cat,x=TIME_ELAPSED,y=GEN_DIST,size=cat),alpha=0.7)+
  geom_hline(yintercept=0.015,linetype='dashed',linewidth=1,col='black') +
  theme_bw(base_size=16)+
  scale_color_manual(name="",values = c(pal[8],'grey70',pal[7]),
                     labels=c('Actual transmission pairs','Unlinked pairs','False phylogenetic signal')) +
  scale_size_manual(name="",values = c(2,1,2)) +
  labs(x='Time elapsed (in years)',y='\n Patristic distance of pair\n(nucleotide substitutions per site)\n') +
  guides(colour = guide_legend(nrow = 1),
         size='none',
         by.col=T) +
  theme(legend.position='bottom', legend.margin=margin(t=-10))+
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(xlim=c(0,15),ylim=c(0,0.2))
p
#ggsave(file=file.path(out.dir, 'simulated_data_after_exclusions_falsepos_bgUnif.png'), p, w=7, h=5)
#ggsave(file=file.path(out.dir, 'simulated_data_after_exclusions_falsepos_bgGMM_5.png'), p, w=7, h=5)
ggsave(file=file.path(out.dir, 'simulated_data_after_exclusions_falsepos_bgGMM-rm-signal.png'), p, w=7, h=5)


# plots for poster 
## all pairs coloured
sim_scenario$true_pair <- ifelse(sim_scenario$TRANSMISSION_PAIR=="NO", 0, 1)
p2 <- ggplot(data=dps, aes(x=d_TSeqT)) + 
  # clock signal cone 
  geom_ribbon(data = dps, aes(ymin = q2.5, ymax = q10), fill = p.palette[1], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q90, ymax = q97.5), fill = p.palette[1], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q10, ymax = q20), fill = p.palette[2], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q80, ymax = q90), fill = p.palette[2], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q20, ymax = q30), fill = p.palette[3], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q70, ymax = q80), fill = p.palette[3], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q30, ymax = q40), fill = p.palette[4], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q60, ymax = q70), fill = p.palette[4], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q40, ymax = q60), fill = p.palette[5], alpha = p.alpha) +
  # transmission pair dots
  geom_point(data=subset(sim_scenario, EXCLUDE==0), aes(x=TIME_ELAPSED,y=GEN_DIST, color=TRANSMISSION_PAIR, alpha=TRANSMISSION_PAIR)) +
  scale_color_manual(name="",values = c(pal[8],pal[2]),
                     labels=c('Unlinked pairs', 'Actual transmission pairs'),
                     guide=guide_legend()) +
  scale_alpha_manual(name="", values=c(0.3, 0.8), guide='none') +
  theme_bw(base_size=16)+
  labs(x='Time elapsed (in years)',y='\n Patristic distance of pair\n(nucleotide substitutions per site)\n') +
  theme(legend.position='bottom', legend.margin=margin(t=-10))+
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(xlim=c(0,15),ylim=c(0,0.2))
p2
#ggsave(file=file.path(out.dir, 'simulated_data_bgUnif.png'), p2, w=8, h=6)
#ggsave(file=file.path(out.dir, 'simulated_data_bgGMM_5.png'), p2, w=8, h=6)
ggsave(file=file.path(out.dir, 'simulated_data_bgGMM_5-rm-signal.png'), p2, w=8, h=6)

## only 'signal' pairs coloured
p3 <- ggplot(data = dps,aes(x=d_TSeqT))+
  geom_ribbon(data = dps, aes(ymin = q2.5, ymax = q10), fill = p.palette[1], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q90, ymax = q97.5), fill = p.palette[1], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q10, ymax = q20), fill = p.palette[2], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q80, ymax = q90), fill = p.palette[2], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q20, ymax = q30), fill = p.palette[3], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q70, ymax = q80), fill = p.palette[3], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q30, ymax = q40), fill = p.palette[4], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q60, ymax = q70), fill = p.palette[4], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q40, ymax = q60), fill = p.palette[5], alpha = p.alpha) +
  geom_point(data=subset(sim_scenario, EXCLUDE==0),aes(color=cat,x=TIME_ELAPSED,y=GEN_DIST,alpha=cat))+
  theme_bw(base_size=16)+
  scale_color_manual(name="",values = c(pal[2],'grey70',pal[8]),
                     labels=c('Actual transmission pairs','Unlinked pairs','False phylogenetic signal')) +
  scale_alpha_manual(name="", values=c(0.8, 0.3, 0.3), guide='none')+
  labs(x='Time elapsed (in years)',y='\n Patristic distance of pair\n(nucleotide substitutions per site)\n') +
  guides(colour = guide_legend(nrow = 1),
         size='none',
         by.col=T) +
  theme(legend.position='bottom', legend.margin=margin(t=-10))+
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(xlim=c(0,15),ylim=c(0,0.2))
p3
#ggsave(file=file.path(out.dir,'simulated_data_after_exclusions_falsepos_bgUnif.png'), p3, w=8, h=6)
#ggsave(file=file.path(out.dir, 'simulated_data_after_exclusions_falsepos_bgGMM_5.png'), p3, w=8, h=6)
ggsave(file=file.path(out.dir, 'simulated_data_after_exclusions_falsepos_bgGMM_5-rm-signal.png'), p3, w=8, h=6)

