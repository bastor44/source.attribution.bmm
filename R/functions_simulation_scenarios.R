io_saveRDS <- function(obj, work_dir, out_dir, base_name, check_if_saved_n = 0)
{
  cat('\nSave to file ',file.path(out_dir, base_name),'...')
  tmp <- check_if_saved_n
  repeat
  {
    #   comp_9 = xzfile(tmp, compression = 9)
    # 	saveRDS(fit.gqs, comp_9)
    tryCatch(
      {
        saveRDS(obj, file=file.path(work_dir, base_name))
        if(work_dir!=out_dir)
        {
          file.copy(file.path(work_dir, base_name),
                    file.path(out_dir, base_name),
                    overwrite = TRUE,
                    recursive = FALSE,
                    copy.mode = TRUE,
                    copy.date = TRUE
          )
        }
      }, error = function(err) { warning(err) } )
    if(check_if_saved_n<1)
      break
    check_if_saved <- try(readRDS(file=file.path(out_dir, base_name)))
    if(!'try-error'%in%class(check_if_saved))
      break
    tmp <- tmp-1
    if(tmp<=0)
    {
      stop('Failed to save ',check_if_saved_n,' times')
    }
  }
}

hivc.db.Date2numeric<- function( x )
{
  if(!class(x)%in%c('Date','character'))	return( x )
  x	<- as.POSIXlt(x)
  tmp	<- x$year + 1900
  x	<- tmp + round( x$yday / ifelse((tmp%%4==0 & tmp%%100!=0) | tmp%%400==0,366,365), d=3 )
  x
}


require(mixtools)
require(condMVNorm)


simulate_scenarios_networks <- function(networks,log_alpha1, log_alpha1_pair_sd, 
                                        log_phi, log_phi_pair_sd, ncases, p_pairs, 
                                        all_pairs=F, dps, sim_binary, excl_improb, 
                                        infile.inftime, outfile.base, 
                                        background_distribution=c("Uniform", "GMM"), gmm_params=NULL){

  seed <- 1234567 # set seed for each chain, changing for each replicate

  set.seed(seed)

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
  tmp <- graph_from_data_frame(tmp, directed=FALSE, vertices=NULL)
  rtc <- data.table(ID=V(tmp)$name, CLU=components(tmp, mode="weak")$membership)
  tmp2 <- rtc[, list(CLU_SIZE=length(ID)), by='CLU'] # max cluster size = 59

  pal <- pal_npg("nrc")(9)
  p <- ggplot(tmp2) + geom_histogram(aes(x=CLU_SIZE),fill=pal[1],binwidth=2) +
    theme_bw(base_size=26) + labs(x='Cluster size',y='count')
  ggsave(file=paste0(outfile.base,'-simulated_cluster_sizes.png'), p, w=12, h=8)

  # subsample clusters to achieve a specified number of transmission events
  if(!is.na(ncases)){

    dt <- merge(dt,rtc,by.x='IdInfected',by.y='ID',all.x=T)

    # flag most recent ncases
    tab <- dt[order(-TimeOfInfection),]
    ids <- unique(tab[1:500,IdInfected])

  }

  # find all possible-pairs ----
  rtc2 <- copy(rtc)
  setnames(rtc,'ID','FROM_ID')
  setnames(rtc2,'ID','TO_ID')
  sg <- merge(rtc,rtc2,by=c('CLU'),all=T,allow.cartesian = T)
  if(all_pairs==T){ # if considering all pairs, ignoring clusters
    sg <- data.table(expand.grid(FROM_ID=rtc$FROM_ID,TO_ID=rtc2$TO_ID))
    sg[,FROM_ID:=as.character(FROM_ID)]
    sg[,TO_ID:=as.character(TO_ID)]
  }

  # remove pairs with same FROM_ID and TO_ID (can't infect themselves)
  sg <- subset(sg,FROM_ID!=TO_ID)

  # merge with pairs data ----
  # first merge with data on infected (time of infection)
  sg <- merge(sg,subset(dt,select=c('IdInfected','TimeOfInfection','TimestepOfInfection','InfectedSPVL','Infected_age')),
              by.x=c('TO_ID'),
              by.y=c('IdInfected'),all.x=T) # just keep the pairs in the subsetted data (otherwise set all=T)
  # then add rows to pairs data for all non-pairs
  if(all_pairs==T){
    sg <- merge(sg, dt, by.x=c('FROM_ID','TO_ID','TimeOfInfection','TimestepOfInfection','InfectedSPVL','Infected_age'),
                by.y=c('IdInfector','IdInfected','TimeOfInfection','TimestepOfInfection','InfectedSPVL','Infected_age'),all=T)
  }else{
    sg <- merge(sg, dt, by.x=c('FROM_ID','TO_ID','TimeOfInfection','TimestepOfInfection','InfectedSPVL','Infected_age','CLU'),
                by.y=c('IdInfector','IdInfected','TimeOfInfection','TimestepOfInfection','InfectedSPVL','Infected_age','CLU'),all=T)
  }

  # just keep rows in which infected ID is in the subset of most recent cases
  sg <- subset(sg, TO_ID %in% ids)

  # add date of infection of infector using both infection date and time elapsed since infection columns
  dt[, TimeOfInfection2:= TimeOfInfection - InfectorTimeElapsedSinceInfection]
  tmp <- unique(subset(dt,select=c('IdInfector','TimeOfInfection2')))
  setnames(tmp,'IdInfector','IdInfected')
  do <- unique(subset(dt,select=c('IdInfected','TimeOfInfection')))
  do <- merge(do,tmp,by='IdInfected',all=T)
  do[is.na(TimeOfInfection), TimeOfInfection:=TimeOfInfection2] # no Ids are missing an infection date
  set(do,NULL,'TimeOfInfection2',NULL)

  setnames(do,c('IdInfected','TimeOfInfection'),c('FROM_ID','FROM_TimeOfInfection'))
  sg <- merge(sg,do,by=c('FROM_ID'),all=T)
  # remove rows with no infected
  sg <- subset(sg,!is.na(TO_ID))

  # define all non-pairs
  sg[is.na(TRANSMISSION_PAIR),TRANSMISSION_PAIR:='NO']

  # add infection date for those who don't appear in 'TO' column (shouldn't be any now)
  tmp <- unique(subset(dt,select=c('IdInfector','TimeOfInfection','InfectorTimeElapsedSinceInfection'),!is.na(InfectorTimeElapsedSinceInfection)))
  tmp[, FROM_TimeOfInfection_miss:= TimeOfInfection - InfectorTimeElapsedSinceInfection]
  tmp <- unique(subset(tmp,select=c('IdInfector','FROM_TimeOfInfection_miss')))
  sg <- merge(sg,subset(tmp,select=c('IdInfector','FROM_TimeOfInfection_miss')),by.x='FROM_ID',by.y='IdInfector',all.x=T)
  sg[is.na(FROM_TimeOfInfection), FROM_TimeOfInfection:= FROM_TimeOfInfection_miss]
  sg[, FROM_TimeOfInfection_miss:= NULL]

  # remove pairs where infector was infected after infected (we allow the same date as there are some (8) true pairs where this is the case??)
  if(excl_improb==T){
    sg <- subset(sg,FROM_TimeOfInfection<=TimeOfInfection) # here we lose half of non-pairs as they are highly unlikely
  }

  cat("\n sub-sample non-pairs")
  if(p_pairs!=1){
    N_samp <- sum(sg$TRANSMISSION_PAIR=='NO') - (sum(sg$TRANSMISSION_PAIR=='YES')*(1-p_pairs))/p_pairs
    set(sg,sample(x=which(sg$TRANSMISSION_PAIR=='NO'), size=N_samp),'REMOVE',"1")
    sg <- subset(sg,is.na(REMOVE))
  }

  # simulate sampling dates ----

  # I use the average time from infection to diagnosis to simulate with (assuming most people are sequenced when they are diagnosed)
  do <- data.table(read.csv(infile.inftime,header=T))
  do <- unique(subset(do,select=c(id,estsctodiagMedian)))

  coeff <- fitdist(do$estsctodiagMedian[!is.na(do$estsctodiagMedian)], distr = "weibull", method = "mle")

  # get time of infection for all transmitters and recipients
  ds <- unique(subset(sg,select=c('TO_ID','TimeOfInfection')))
  tmp <- unique(subset(sg,select=c('FROM_ID','FROM_TimeOfInfection')))
  setnames(tmp,c('FROM_ID','FROM_TimeOfInfection'),c('TO_ID','TimeOfInfection'))
  ds <- unique(rbind(ds,tmp))

  # simulate rng from weibull dist and add to the estimated time of infection
  ds[, date_s:= TimeOfInfection + rweibull(nrow(ds),shape=coeff$estimate[['shape']],scale=coeff$estimate[['scale']])]

  # merge the infector and infected sampling dates
  df <- merge(sg,ds,by=c('TO_ID','TimeOfInfection'),all.x=T)
  setnames(df,'date_s','TO_DATE_S')
  df <- merge(df,subset(ds,select=c('TO_ID','date_s')),by.x=c('FROM_ID'),by.y=c('TO_ID'),all.x=T)
  setnames(df,'date_s','FROM_DATE_S')

  # calculate time elapsed
  df[, TIME_ELAPSED:= abs(TO_DATE_S - TimeOfInfection) + abs(FROM_DATE_S - TimeOfInfection)]

  # only keep pairs with time elapsed <16 years
  df <- subset(df,TIME_ELAPSED<16)

  cat(paste0("Earliest infection date of incident cases: ", as.Date(date_decimal(min(df$TimeOfInfection)),format='%d-%m-%Y'),'\n'))
  cat(paste0("Latest infection date of incident cases: ", as.Date(date_decimal(max(df$TimeOfInfection)),format='%d-%m-%Y'),'\n'))

  # simulate ages
  # get age of missing infectors
  df <- merge(df,subset(bdate,select=c('IdInfected','bdate')),by.x='FROM_ID',by.y='IdInfected',all.x=T)
  df[is.na(Infector_age), Infector_age:= TimeOfInfection - bdate]
  # some still missing, so simulate uniformly
  df[is.na(Infector_age), Infector_age:= runif(length(which(is.na(Infector_age))), 16, 70)]
  setnames(df,c('Infector_age','Infected_age'),c('FROM_AGE','TO_AGE'))
  # set any <15 to 15
  df[FROM_AGE<15,FROM_AGE:=15]
  df[TO_AGE<15,TO_AGE:=15]

  # re-simulate ages
  df[which(df$TRANSMISSION_PAIR=='YES'), FROM_AGE := rlnormTrunc(length(which(df$TRANSMISSION_PAIR=='YES')), log(30), log(1.3), min = 16, max = 75)]
  df[which(df$TRANSMISSION_PAIR=='YES'), TO_AGE := rlnormTrunc(length(which(df$TRANSMISSION_PAIR=='YES')), log(FROM_AGE), log(1.25), min = (16), max = (75))]
  #df[which(df$TRANSMISSION_PAIR=='YES'), TO_AGE := rlnormTrunc(length(which(df$TRANSMISSION_PAIR=='YES')), log(FROM_AGE), log(1.8), min = (16), max = (75))] # for less correlated ages

  #
  df[which(df$TRANSMISSION_PAIR=='NO'), FROM_AGE := runif(length(which(df$TRANSMISSION_PAIR == "NO")), 16, 75)]
  df[which(df$TRANSMISSION_PAIR=='NO'), TO_AGE := runif(length(which(df$TRANSMISSION_PAIR == "NO")), 16, 75)]


  # simulate genetic distances ----

  # first simulate for true pairs from gamma dist
  df[, log_alpha1_pair := rnorm(nrow(df), 0, log_alpha1_pair_sd)]
  df[, log_phi_pair := rnorm(nrow(df), 0, log_phi_pair_sd)]
  df[, ER := exp(log_alpha1 + log_alpha1_pair)]
  df[, BETA := exp( -(log_phi + log_phi_pair))]

  sim_scenario <- as.data.table(df[,c('FROM_ID','TO_ID','TRANSMISSION_PAIR','TRANS_STAGE','FROM_AGE','TO_AGE','TIME_ELAPSED','ER','BETA')])
  # simulate binary covariate
  sim_scenario[TRANSMISSION_PAIR=="YES", BIN_COV:= 'cat1']
  sim_scenario[TRANSMISSION_PAIR=="NO", BIN_COV:= 'cat2']
  sim_scenario[TRANSMISSION_PAIR=="YES", GAMMA_SHAPE:= ER * TIME_ELAPSED * BETA]
  sim_scenario[TRANSMISSION_PAIR=="YES", GAMMA_RATE:= BETA]
  sim_scenario[TRANSMISSION_PAIR=="YES", GEN_DIST:= rgamma(length(which(sim_scenario$TRANSMISSION_PAIR=="YES")),shape = GAMMA_SHAPE,rate=GAMMA_RATE)]
  
  # simulate distances for non-transmission pairs (max 20% for uniform)
  if (background_distribution == "Uniform") {
    sim_scenario[TRANSMISSION_PAIR=="NO", GEN_DIST:= runif(length(which(df$TRANSMISSION_PAIR=="NO")),0,0.2)]
  }
  # simulate distances for non-transmission pairs GMM
  if (background_distribution =="GMM") {
    k <- gmm_params$nclass
    p <- gmm_params$p
    #sigma <- LTSigma2variance(gmm_params$LTSigma)
    sigma <- gmm_params$LTSigma
    
    # needs to be conditional on times already generated 
    # assign to cluster 
    sim_scenario[TRANSMISSION_PAIR=="NO", CLUSTER:=sample(1:k, size=length(which(df$TRANSMISSION_PAIR=="NO")), replace=TRUE, prob=gmm_params$pi)]
    # conditional on TIME_ELAPSED, generate GEN_DIST from mixture model 
    sim_scenario[TRANSMISSION_PAIR=="NO", GEN_DIST:=mapply(CLUSTER, TIME_ELAPSED, FUN=function(cl, t) {rcmvnorm(n=1, mean=gmm$Mu[cl, ], sigma=sigma[,,cl], dependent.ind=2, given.ind=1, X.given=t)})]
    #sim_scenario[TRANSMISSION_PAIR=="NO", GEN_DIST:=abs(GEN_DIST)]
  }
  
  
  sim_scenario[, TRANS_STAGE:= 'Diagnosed']
  sim_scenario[TRANSMISSION_PAIR=="YES", TRANS_STAGE:= 'Undiagnosed']

  # prepare for stan data
  sim_scenario$TRANSMISSION_PAIR <- as.factor(sim_scenario$TRANSMISSION_PAIR)
  levels(sim_scenario$TRANSMISSION_PAIR) <- c("No","Yes")
  sim_scenario[,SOURCE_PAIR:= paste(BIN_COV,"-",TRANSMISSION_PAIR)]
  sim_scenario[,LOG_TIME_ELAPSED:=log(TIME_ELAPSED)]
  sim_scenario[,LOG_GEN_DIST:= log(GEN_DIST)]

  tmp <- data.table(TIME_ELAPSED = sort(unique(sim_scenario$TIME_ELAPSED)))
  tmp[, IDX_UNIQUE_TE := seq_len(nrow(tmp))]
  sim_scenario <- merge(sim_scenario, tmp, by = 'TIME_ELAPSED')
  sim_scenario[, PAIR_ID := seq_len(nrow(sim_scenario)) ]

  saveRDS(sim_scenario,file = paste0(outfile.base,'.rds'))

  # plot
  make_plot_simulated_data(dps,sim_scenario,outfile.base)
  plot_logdistance_timeelapsed(sim_scenario,outfile.base)
}


make_plot_simulated_data <- function(dps,sim_scenario,outfile.base){

  sim_scenario[,SOURCE := 'Source category 1']
  sim_scenario[BIN_COV=='cat2', SOURCE := 'Source category 2']

  sim_scenario[SOURCE_PAIR=='Source category 1 - Yes', SOURCE_PAIR := 'Source category 1 - Yes']
  sim_scenario[SOURCE_PAIR=='Source category 2 - No',  SOURCE_PAIR := 'Source category 2 - No']

  p.palette <- RColorBrewer::brewer.pal(5,'Oranges')
  p.alpha <- 0.7

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
    theme_bw(base_size=16)+geom_point(data=sim_scenario,aes(color=SOURCE,x=TIME_ELAPSED,y=GEN_DIST,shape=TRANSMISSION_PAIR))+
    #scale_color_manual(name="Stage of probable \n transmitter & \n Transmission pair",values = c("black","green","yellow2"))+
    scale_colour_npg() +
    labs(x='\n Time elapsed (in years)',y='\n Genetic distance of the pair \n',
         colour = 'Source \ncategory',
         shape = 'Transmission \npair') +
    guides(shape = guide_legend(nrow = 2),
           colour = guide_legend(nrow = 2),by.col=T) +
    theme(legend.position='bottom')+
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    coord_cartesian(xlim=c(0,15),ylim=c(0,0.16))
  p

  ggsave(file=paste0(outfile.base,'-lognorm_gamma.pdf'), p, w=7, h=7)

  if(!is.null(sim_scenario$TRANS_STAGE)){
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
      theme_bw(base_size=16)+geom_point(data=sim_scenario,aes(color=SOURCE_PAIR,x=TIME_ELAPSED,y=GEN_DIST,shape=TRANS_STAGE))+
      #scale_color_manual(name="Stage of probable \n transmitter & \n Transmission pair",values = c("black","green","yellow2"))+
      scale_colour_npg() +
      labs(x='\n Time elapsed (in years)',y='Genetic distance of the pair \n',shape='Stage of transmitter\nat estimated time of\ninfection of recipient',
           colour = 'Source category of\nprobable transmitter\n- Transmission pair')+
      theme(legend.position='bottom',
            legend.title=element_text(size=rel(0.65)),
            legend.text=element_text(size=rel(0.65)))+
      scale_x_continuous(expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0)) +
      coord_cartesian(xlim=c(0,15),ylim=c(0,0.16)) +
    guides(shape = guide_legend(nrow = 2),
             colour = guide_legend(nrow = 2),by.col=T)
    p

    ggsave(file=paste0(outfile.base,'-lognorm_gamma_transm_stage.pdf'), p, w=8, h=7)
  }
}

make_plot_simulated_data_colour_prob_tpair <- function(rep,po,dps_clock,sim_scenario,outfile.base){

  po[SOURCE_PAIR=='Source category 1 - Yes', SOURCE_PAIR := 'Source category 1 - Yes']
  po[SOURCE_PAIR=='Source category 2 - No',  SOURCE_PAIR := 'Source category 2 - No']

  p.palette <- RColorBrewer::brewer.pal(5,'Oranges')
  p.alpha <- 0.7

  p <- ggplot(data = dps_clock,aes(x=d_TSeqT))+
    geom_ribbon(data = dps_clock, aes(ymin = q2.5, ymax = q10), fill = p.palette[1], alpha = p.alpha) +
    geom_ribbon(data = dps_clock, aes(ymin = q90, ymax = q97.5), fill = p.palette[1], alpha = p.alpha) +
    geom_ribbon(data = dps_clock, aes(ymin = q10, ymax = q20), fill = p.palette[2], alpha = p.alpha) +
    geom_ribbon(data = dps_clock, aes(ymin = q80, ymax = q90), fill = p.palette[2], alpha = p.alpha) +
    geom_ribbon(data = dps_clock, aes(ymin = q20, ymax = q30), fill = p.palette[3], alpha = p.alpha) +
    geom_ribbon(data = dps_clock, aes(ymin = q70, ymax = q80), fill = p.palette[3], alpha = p.alpha) +
    geom_ribbon(data = dps_clock, aes(ymin = q30, ymax = q40), fill = p.palette[4], alpha = p.alpha) +
    geom_ribbon(data = dps_clock, aes(ymin = q60, ymax = q70), fill = p.palette[4], alpha = p.alpha) +
    geom_ribbon(data = dps_clock, aes(ymin = q40, ymax = q60), fill = p.palette[5], alpha = p.alpha) +
    theme_bw()+geom_point(data=po,aes(shape=SOURCE_PAIR,x=TIME_ELAPSED,y=GEN_DIST,colour=M))+
    scale_colour_gradient(name="Posterior probability \n of being \n a transmission pair",
                          low = "grey90",
                          high = muted("blue")) +
    #scale_shape_manual(name="Stage of probable \n transmitter & \n Transmission pair",values = c(19,15))+
    labs(x='\n Time elapsed (in years)',y='\n Genetic distance of the pair \n',shape="Source category of \nprobable transmitter -\nTransmission pair")+
    theme(legend.position='bottom')+
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    coord_cartesian(xlim=c(0,15),ylim=c(0,0.2)) +
    guides(shape = guide_legend(nrow = 2),by.col=T)
  p

  ggsave(file=paste0(outfile.base,'-rep_',rep,'-sim_data_prob_tpair.pdf'), p, w=9, h=7)

  p <- ggplot(data = dps_clock,aes(x=d_TSeqT))+
    geom_ribbon(data = dps_clock, aes(ymin = q2.5, ymax = q10), fill = p.palette[1], alpha = p.alpha) +
    geom_ribbon(data = dps_clock, aes(ymin = q90, ymax = q97.5), fill = p.palette[1], alpha = p.alpha) +
    geom_ribbon(data = dps_clock, aes(ymin = q10, ymax = q20), fill = p.palette[2], alpha = p.alpha) +
    geom_ribbon(data = dps_clock, aes(ymin = q80, ymax = q90), fill = p.palette[2], alpha = p.alpha) +
    geom_ribbon(data = dps_clock, aes(ymin = q20, ymax = q30), fill = p.palette[3], alpha = p.alpha) +
    geom_ribbon(data = dps_clock, aes(ymin = q70, ymax = q80), fill = p.palette[3], alpha = p.alpha) +
    geom_ribbon(data = dps_clock, aes(ymin = q30, ymax = q40), fill = p.palette[4], alpha = p.alpha) +
    geom_ribbon(data = dps_clock, aes(ymin = q60, ymax = q70), fill = p.palette[4], alpha = p.alpha) +
    geom_ribbon(data = dps_clock, aes(ymin = q40, ymax = q60), fill = p.palette[5], alpha = p.alpha) +
    theme_bw()+
    geom_point(data=subset(po,TRANSMISSION_PAIR=='No'),aes(shape=SOURCE_PAIR,x=TIME_ELAPSED,y=GEN_DIST,colour=M))+
    scale_colour_gradient(name="Posterior probability \n of being \n a transmission pair",
                          low = "grey90",
                          high = muted("blue")) +
    #scale_shape_manual(name="Stage of probable \n transmitter & \n Transmission pair",values = c(19,15))+
    labs(x='\n Time elapsed (in years)',y='\n Genetic distance of the pair \n',shape="Source category of \nprobable transmitter -\nTransmission pair")+
    theme(legend.position='bottom')+
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    coord_cartesian(xlim=c(0,15)) +
    guides(shape = guide_legend(nrow = 2),by.col=T)
  p
  ggsave(file=paste0(outfile.base,'-rep_',rep,'-sim_data_prob_tpair_nontpairs.pdf'), p, w=9, h=7)

}

make_plot_simulated_data_colour_stage_prob_tpair <- function(rep,po,dps_clock,sim_scenario,outfile.base){

  #sim_scenario <- merge(dps,subset(sim_scenario,select=-TRANSMISSION_PAIR),by=c('PAIR_ID','TRANS_STAGE'),all=T)
  p.palette <- RColorBrewer::brewer.pal(5,'Oranges')
  p.alpha <- 0.7

  p <- ggplot(data = dps_clock,aes(x=d_TSeqT))+
    geom_ribbon(data = dps_clock, aes(ymin = q2.5, ymax = q10), fill = p.palette[1], alpha = p.alpha) +
    geom_ribbon(data = dps_clock, aes(ymin = q90, ymax = q97.5), fill = p.palette[1], alpha = p.alpha) +
    geom_ribbon(data = dps_clock, aes(ymin = q10, ymax = q20), fill = p.palette[2], alpha = p.alpha) +
    geom_ribbon(data = dps_clock, aes(ymin = q80, ymax = q90), fill = p.palette[2], alpha = p.alpha) +
    geom_ribbon(data = dps_clock, aes(ymin = q20, ymax = q30), fill = p.palette[3], alpha = p.alpha) +
    geom_ribbon(data = dps_clock, aes(ymin = q70, ymax = q80), fill = p.palette[3], alpha = p.alpha) +
    geom_ribbon(data = dps_clock, aes(ymin = q30, ymax = q40), fill = p.palette[4], alpha = p.alpha) +
    geom_ribbon(data = dps_clock, aes(ymin = q60, ymax = q70), fill = p.palette[4], alpha = p.alpha) +
    geom_ribbon(data = dps_clock, aes(ymin = q40, ymax = q60), fill = p.palette[5], alpha = p.alpha) +
    theme_bw()+geom_point(data=po,aes(shape=STAGE_PAIR,x=TIME_ELAPSED,y=GEN_DIST,colour=M))+
    scale_colour_gradient(name="Posterior probability \n of being \n a transmission pair",
                          low = "grey90",
                          high = muted("blue")) +
    #scale_shape_manual(name="Stage of probable \n transmitter & \n Transmission pair",values = c(19,15))+
    labs(x='\n Time elapsed (in years)',y='\n Genetic distance of the pair \n',shape="Birth place of\nprobable transmitter\n- Transmission pair")+
    theme(legend.position='bottom')+
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    coord_cartesian(xlim=c(0,15),ylim=c(0,0.2))
  p

  ggsave(file=paste0(outfile.base,'-rep_',rep,'-sim_data_stage_prob_tpair.pdf'), p, w=9, h=7)

  p <- ggplot(data = dps_clock,aes(x=d_TSeqT))+
    geom_ribbon(data = dps_clock, aes(ymin = q2.5, ymax = q10), fill = p.palette[1], alpha = p.alpha) +
    geom_ribbon(data = dps_clock, aes(ymin = q90, ymax = q97.5), fill = p.palette[1], alpha = p.alpha) +
    geom_ribbon(data = dps_clock, aes(ymin = q10, ymax = q20), fill = p.palette[2], alpha = p.alpha) +
    geom_ribbon(data = dps_clock, aes(ymin = q80, ymax = q90), fill = p.palette[2], alpha = p.alpha) +
    geom_ribbon(data = dps_clock, aes(ymin = q20, ymax = q30), fill = p.palette[3], alpha = p.alpha) +
    geom_ribbon(data = dps_clock, aes(ymin = q70, ymax = q80), fill = p.palette[3], alpha = p.alpha) +
    geom_ribbon(data = dps_clock, aes(ymin = q30, ymax = q40), fill = p.palette[4], alpha = p.alpha) +
    geom_ribbon(data = dps_clock, aes(ymin = q60, ymax = q70), fill = p.palette[4], alpha = p.alpha) +
    geom_ribbon(data = dps_clock, aes(ymin = q40, ymax = q60), fill = p.palette[5], alpha = p.alpha) +
    theme_bw()+
    geom_point(data=subset(po,TRANSMISSION_PAIR=='No'),aes(shape=STAGE_PAIR,x=TIME_ELAPSED,y=GEN_DIST,colour=M))+
    scale_colour_gradient(name="Posterior probability \n of being \n a transmission pair",
                          low = "grey90",
                          high = muted("blue")) +
    #scale_shape_manual(name="Stage of probable \n transmitter & \n Transmission pair",values = c(19,15))+
    labs(x='\n Time elapsed (in years)',y='\n Genetic distance of the pair \n',shape="Birth place of\nprobable transmitter\n- Transmission pair")+
    theme(legend.position='bottom')+
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    coord_cartesian(xlim=c(0,15))
  p
  ggsave(file=paste0(outfile.base,'-rep_',rep,'-sim_data_stage_prob_tpair_nontpairs.pdf'), p, w=9, h=7)

}

make_plot_applied_data_colour_stage_prob_tpair <- function(rep,po,dps_clock,sim_scenario,outfile.base){

  #sim_scenario <- merge(dps,subset(sim_scenario,select=-TRANSMISSION_PAIR),by=c('PAIR_ID','TRANS_STAGE'),all=T)
  p.palette <- RColorBrewer::brewer.pal(5,'Oranges')
  p.alpha <- 0.7

  p <- ggplot(data = dps_clock,aes(x=d_TSeqT))+
    geom_ribbon(data = dps_clock, aes(ymin = q2.5, ymax = q10), fill = p.palette[1], alpha = p.alpha) +
    geom_ribbon(data = dps_clock, aes(ymin = q90, ymax = q97.5), fill = p.palette[1], alpha = p.alpha) +
    geom_ribbon(data = dps_clock, aes(ymin = q10, ymax = q20), fill = p.palette[2], alpha = p.alpha) +
    geom_ribbon(data = dps_clock, aes(ymin = q80, ymax = q90), fill = p.palette[2], alpha = p.alpha) +
    geom_ribbon(data = dps_clock, aes(ymin = q20, ymax = q30), fill = p.palette[3], alpha = p.alpha) +
    geom_ribbon(data = dps_clock, aes(ymin = q70, ymax = q80), fill = p.palette[3], alpha = p.alpha) +
    geom_ribbon(data = dps_clock, aes(ymin = q30, ymax = q40), fill = p.palette[4], alpha = p.alpha) +
    geom_ribbon(data = dps_clock, aes(ymin = q60, ymax = q70), fill = p.palette[4], alpha = p.alpha) +
    geom_ribbon(data = dps_clock, aes(ymin = q40, ymax = q60), fill = p.palette[5], alpha = p.alpha) +
    theme_bw()+geom_point(data=po,aes(shape=TRANS_STAGE,x=TIME_ELAPSED,y=GEN_DIST,colour=M))+
    scale_colour_gradient(name="Posterior probability \n of being \n a transmission pair",
                          low = "grey90",
                          high = muted("blue")) +
    #scale_shape_manual(name="Stage of probable \n transmitter & \n Transmission pair",values = c(19,15))+
    labs(x='\n Time elapsed (in years)',y='\n Genetic distance of the pair \n',shape="Birth place of\nprobable transmitter\n- Transmission pair")+
    theme(legend.position='bottom')+
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    coord_cartesian(xlim=c(0,15),ylim=c(0,0.15))
  p

  ggsave(file=paste0(outfile.base,'-rep_',rep,'-app_data_stage_prob_tpair.pdf'), p, w=9, h=7)

}

plot_logdistance_timeelapsed <- function(sim_scenario,outfile.base){

  p <- ggplot(sim_scenario,aes(x=TIME_ELAPSED,y=LOG_GEN_DIST,color=TRANSMISSION_PAIR))+theme_bw()+geom_point()+
    labs(x='\n Time elapsed (in years) ', y='Log of genetic distance of pairs of individuals \n')
  p

  ggsave(file=paste0(outfile.base,'-log_distance_time.pdf'), p, w=8, h=6)

}

plot_posterior_probabilities_trsm_pair <- function(tmp,i){

    pal <- pal_npg("nrc")(4)[3:4]

    p <- ggplot(tmp,aes(x=PAIR_ID,y=M,color=TRANS_STAGE))+geom_point()+
      geom_errorbar(aes(ymin=CL_L, ymax=CL_U)) + facet_grid(~TRANSMISSION_PAIR,scales = "free_x")+
      theme_bw() +
      scale_color_manual(name="Transmission stage",values = pal)+
      labs(x='\n Pair ID ', y='Posterior Probability of being a transmission pair \n with 95% confidence interval \n') +
      coord_cartesian()+scale_y_continuous(labels = scales::percent)
    p
    ggsave(file=file.path(out.dir,paste0('sim_scenario_',i,'-posterior_probability_of_being_transmission_pair.png')), p, w=12, h=7)

    tmp[, flab:= paste0(TRANS_STAGE,' - ',TRANSMISSION_PAIR)]
    p <- ggplot(tmp,aes(x=PAIR_ID,y=M,color=TRANS_STAGE,shape=TRANSMISSION_PAIR))+geom_point()+
      geom_errorbar(aes(ymin=CL_L, ymax=CL_U)) + facet_grid(~flab,scales = "free_x")+
      theme_bw() +
      scale_shape_manual(name="Transmission pair",values=c(20,4)) +
      scale_color_manual(name="Transmission stage",values = pal)+
      labs(x='\n Pair ID ', y='Posterior Probability of being a transmission pair \n with 95% confidence interval \n') +
      coord_cartesian()+scale_y_continuous(labels = scales::percent)
    p
    ggsave(file=file.path(out.dir,paste0('sim_scenario_',i,'-posterior_probability_of_being_transmission_pair_trsmgroup.png')), p, w=12, h=7)

    p <- ggplot(tmp,aes(x=PAIR_ID,y=M,color=TRANS_STAGE,group=PAIR_ID)) +
      geom_boxplot(
        aes(ymin = `20%`, lower = `25%`, middle = M, upper = `75%`, ymax = `80%`),
        stat = "identity"
      ) +
     facet_grid(~flab,scales = "free_x") +
      theme_bw(base_size=20) +
      #theme(strip.background = element_blank()) +
      scale_color_manual(name="Transmission stage",values = pal)+
      labs(x='\n Pair ID ', y='Posterior Probability of being a transmission pair \n with 80% credible interval \n') +
      coord_cartesian()+scale_y_continuous(labels = scales::percent)
    p
    ggsave(file=file.path(out.dir,paste0('-posterior_probability_of_being_transmission_pair_trsmgroup_boxwhisker.png')), p, w=17, h=7)


}

plot_estimated_attributable_fraction <- function(tmp,plot,outfile.base){

  tmp[, SOURCE:= 'Source category 1']
  tmp[BIN_COV=='cat2', SOURCE:= 'Source category 2']

  pal <- pal_npg("nrc")(4)[3:4]

  p <- ggplot(tmp)+geom_point(aes(x=as.factor(SOURCE),y=M,color=as.factor(SOURCE)))+geom_errorbar(aes(x=as.factor(SOURCE),y=M,ymin=CL, ymax=CU,color=as.factor(SOURCE)),width=0.5) +
    scale_color_manual(name="Transmission stage",values = pal) +
    geom_point(aes(x=as.factor(SOURCE),y=true_p)) +
    theme_bw(base_size=22) + theme(legend.position = "none") +
    labs(x='\n Source category of transmitter', y='Attributable proportion of \ntransmission pairs to stage  \n') +
    coord_cartesian(ylim = c(0,1))+scale_y_continuous(labels = scales::label_percent(accuracy = 1L),breaks=seq(0,1,0.1))
  p
  ggsave(file=paste0(outfile.base,'-estimated_attributable_fraction',plot,'.pdf'), p, w=8, h=7)

}

plot_PAF_age <- function(tmp,thrsh,plot,outfile.base){

  tmp[, SOURCE:= factor(FROM_AGE_GP, levels=c("[15-30)", "[30-40)", "[40-50)", "[50-60)", "[60+)"),
                        labels=c("15-29", "30-39", "40-49", "50-59", "60+"  ))]

  tmp <- merge(tmp,thrsh,by='FROM_AGE_GP')
  tmp <- melt(tmp,id.vars=c('FROM_AGE_GP','SOURCE','CL','CU','IL','IU','true_p'),value.vars=c('M','pct_thrsh'))
  tmp[variable=='pct_thrsh', CL:= NA]
  tmp[variable=='pct_thrsh', CU:= NA]

  #  if(is.null(tmp$true_p)){
  #    setnames(tmp,'paf','true_p')
  #  }

  pal1 <- pal_npg("nrc")(4)[1]
  pal2 <- pal_npg("nrc")(4)[2]
  pal3 <- pal_npg("nrc")(4)[3]

  p <- ggplot(tmp) +
    #geom_bar(aes(x=as.factor(SOURCE),y=true_p,color=as.factor(SOURCE)), fill=NA) +
    geom_bar(data=subset(tmp,variable=='M'),aes(x=as.factor(SOURCE),y=true_p,fill='white',color='grey50'),stat='identity',linetype=2) +
    geom_bar(aes(x=as.factor(SOURCE),y=value,fill=as.factor(variable)),stat='identity', position = position_dodge(width=0.9),alpha=0.5) +
    geom_errorbar(aes(x=as.factor(SOURCE),ymin=CL, ymax=CU,color=as.factor(variable)),width=0.3,position = position_dodge(width=0.9)) +
    scale_fill_manual(name="",values = c(pal1,pal2,'white'),labels=c('Estimated flows','1.5% threshold','')) +
    scale_color_manual(name="",values = c('grey50','grey50','grey50'),labels=c('True flows','','')) +
    #geom_point(data=thrsh,aes(x=as.factor(SOURCE),y=value,color=as.factor(SOURCE), shape= as.factor(shp)), position = position_dodge(width=0.5)) +
    #scale_color_manual(name="Estimated proportions",values = c(pal1,pal2),labels=c('Source category 1','Source category 2')) +
    #scale_shape_manual(name="Estimated proportions\nfixed threshold",values = c(17),labels=c('1.5% threshold')) +
    theme_bw(base_size=16) + theme(legend.position = "bottom") +
    labs(x='\n Source category of transmitter', y='Attributable proportion of \ntransmission pairs to source age  \n') +
    coord_cartesian(ylim = c(0,1))+scale_y_continuous(labels = scales::label_percent(accuracy = 1L),breaks=seq(0,1,0.1)) +
    guides(line = guide_legend(nrow = 2),
           fill = guide_legend(nrow = 2,title.position='top',order=1),
           color = guide_legend(nrow = 2,title.position='top'), order=2, by.col=T)

  p

  ggsave(file=paste0(outfile.base,'-estimated_attributable_fraction',plot,'.pdf'), p, w=8, h=7)

}


plot_PAF_age_no_thrsh <- function(tmp, plot, outfile.base){
  pal <- pal_brewer(palette='Paired')(10)[9]
  tmp[, SOURCE:= factor(FROM_AGE_GP, levels=c("[15-30)", "[30-40)", "[40-50)", "[50-60)", "[60+)"),
                        labels=c("15-29", "30-39", "40-49", "50-59", "60+"  ))]
  
  tmp <- melt(tmp,id.vars=c('FROM_AGE_GP','SOURCE','CL','CU','IL','IU','true_p'),value.vars=c('M'))
  
  p <- ggplot(tmp) +
    geom_bar(data=subset(tmp, variable='M'), aes(x=as.factor(SOURCE), y=true_p, colour='grey50', fill='white'),stat='identity', linetype=2) +
    geom_bar(aes(x=as.factor(SOURCE), y=value, color=pal, fill=pal), stat='identity', alpha=0.5) +
    geom_errorbar(aes(x=as.factor(SOURCE), ymin=CL, ymax=CU), color='grey50', width=0.3) +
    scale_fill_manual(name="", values=c(pal, 'white'), labels=c('Estimated flows', '')) +
    scale_color_manual(name="", values=c('grey50', 'grey50'), labels=c('True flows', '')) +
    theme_bw(base_size=14) + theme(legend.position = "bottom") +
    labs(x='\n Source category of transmitter', y='Proportion of transmission pairs \nattributible to source age') +
    coord_cartesian(ylim = c(0,1))+scale_y_continuous(labels = scales::label_percent(accuracy = 1L),breaks=seq(0,1,0.1))
  
  ggsave(file=paste0(outfile.base,'-estimated_attributable_fraction-nothrsh',plot,'.pdf'), p, w=8, h=7)
  return(p)
}

plot_estimated_attributable_fraction_applied <- function(tmp,rep,plot,outfile.base){

  pal <- pal_npg("nrc")(4)[3:4]

  p <- ggplot(tmp)+geom_point(aes(x=as.factor(BIN_COV),y=M,color=as.factor(BIN_COV)))+geom_errorbar(aes(x=as.factor(BIN_COV),y=M,ymin=CL, ymax=CU,color=as.factor(BIN_COV)),width=0.5) +
    scale_color_manual(name="Transmission stage",values = pal) +
    theme_bw(base_size=22) + theme(legend.position = "none") +
    labs(x='\n Birth place of the transmitter', y='Attributable proportion of \ntransmission pairs to stage  \n') +
    coord_cartesian(ylim = c(0,1))+scale_y_continuous(labels = scales::label_percent(accuracy = 1L),breaks=seq(0,1,0.1))
  p
  ggsave(file=paste0(outfile.base,'-rep_',rep,'-estimated_attributable_fraction',plot,'.pdf'), p, w=8, h=7)

}

plot_estimated_attributable_fraction_applied_stratify <- function(tmp,rep,plot,outfile.base){

  pal <- pal_npg("nrc")(4)[3:4]

  p <- ggplot(tmp)+geom_point(aes(x=as.factor(BIN_COV),y=M,color=as.factor(TO_BIN_COV)), position = position_dodge(width=0.9)) +
    geom_errorbar(aes(x=as.factor(BIN_COV),y=M,ymin=CL, ymax=CU,color=as.factor(TO_BIN_COV)),width=0.5,position=position_dodge(width=0.9)) +
    scale_color_manual(name="Birth place of recipient",values = pal) +
    theme_bw(base_size=22) + theme(legend.position = "right") +
    labs(x='\n Birth place of the transmitter', y='Attributable proportion of \ntransmission pairs to stage  \n') +
    coord_cartesian(ylim = c(0,1)) +
    scale_y_continuous(labels = scales::label_percent(accuracy = 1L),breaks=seq(0,1,0.1))
  p
  ggsave(file=paste0(outfile.base,'-rep_',rep,'-estimated_attributable_fraction_bybincov_recipient',plot,'.pdf'), p, w=11, h=7)

}
