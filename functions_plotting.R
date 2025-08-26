
posterior_prob_pair <- function(model_fit,sim_scenario,dps){
  po <- model_fit$draws(inc_warmup = FALSE,
                        format = 'draws_df',
                        variables = 'tpair_prob_w'
  )
  po <- as.data.table(po)
  setnames(po, colnames(po), gsub('^\\.','',colnames(po)))
  po <- melt(po, id.vars = c('chain','iteration','draw'))
  po <- po[,
           list( q = quantile(value, probs = c(0.5, 0.025, 0.25, 0.75, 0.975) ),
                 stat = c('M','CL','IL', 'IU', 'CU')
           ),
           by = c('variable')
  ]
  po <- dcast.data.table(po, variable~stat, value.var = 'q')
  po[, PAIR_ID := as.integer(gsub('tpair_prob_w\\[([0-9]+)\\]','\\1',as.character(variable)))]
  po <- merge(po, sim_scenario, by = 'PAIR_ID')
  tmp <- po[, .(PAIR_ID = PAIR_ID, PAIR_ID2 = seq_along(PAIR_ID)), by = 'TRANSMISSION_PAIR']
  po <- merge(po, tmp, by = c('TRANSMISSION_PAIR','PAIR_ID'))
  
  po[, d_TSeqT:= round(TIME_ELAPSED,1)]
  po <- merge(po,dps,by=c('d_TSeqT'),all.x=T)
  # identify % points in signal cone with trsm prob >50%
  po[, incone:=0]
  po[GEN_DIST>q20 & GEN_DIST<q80, incone:=1]
  po[, prob_pair:=0]
  po[M>0.5, prob_pair:=1]
  return(po)
}

make_plot_simulated_data_colour_prob_tpair <- function(po,dps_clock,sim_scenario,outfile.base){
  
  po[SOURCE_PAIR=='Source category 1 - No', SOURCE_PAIR:='No' ]
  po[SOURCE_PAIR=='Source category 2 - Yes', SOURCE_PAIR:='Yes' ]
  
  p.palette <- RColorBrewer::brewer.pal(5,'Blues')
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
    theme_bw(base_size=16)+geom_point(data=po,aes(shape=TRANSMISSION_PAIR,x=TIME_ELAPSED,y=GEN_DIST,colour=M))+
    scale_colour_gradient(name="Transmission pair\nprobability under model\n(posterior median)",
                          low = "grey90",
                          high = muted("blue"),
                          breaks = seq(0,1,0.25),
                          n.breaks=4,
                          labels = scales::label_number(accuracy = 0.01),
                          guide=guide_colourbar(direction='horizontal',barwidth=15)) +
    labs(x='Time elapsed (in years)',y='\nPatristic distance of\nphylogenetically possible pair',shape="True transmission pair")+
    theme(legend.position='bottom',
          legend.title=element_text(size=rel(0.7)),
          legend.text=element_text(size=rel(0.7)))+
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0),labels = scales::percent, limits = c(0,1)) +
    coord_cartesian(xlim=c(0,15),ylim=c(0,0.2)) +
    guides(shape = guide_legend(nrow = 2,order=1),by.col=T,colour = guide_colorbar(order=2,barwidth=15))
  p
  
}

make_plot_simulated_data_colour_prob_tpair_sizes <- function(po,dps_clock,sim_scenario,outfile.base){
  
  p.palette <- RColorBrewer::brewer.pal(5,'Blues')
  p.alpha <- 0.7
  
  po[TRANSMISSION_PAIR=='No', cat:= ifelse(GEN_DIST>=q2.5 & GEN_DIST<=q97.5,'false_pos',
                                           'true_neg')]
  po[TRANSMISSION_PAIR=='Yes', cat:= ifelse(GEN_DIST<q2.5 | GEN_DIST>q97.5,'false_neg',
                                            'true_pos')]
  po[, cat:=factor(cat,levels=c('true_neg','false_pos','true_pos','false_neg'),
                   labels=c('Unlinked pair with\nno signal',
                            'Unlinked pair with\nfalse signal',
                            'Transmission pair\nwith signal',
                            'Transmission pair\nno signal'))]
  
  #pal <- pal_npg("nrc")(5)
  pal <- pal_brewer(palette='Paired')(10)
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
    #theme_bw(base_size=16)+geom_point(data=po,aes(x=TIME_ELAPSED,y=GEN_DIST,size=M,colour=cat))+
    theme_bw(base_size=14)+geom_point(data=po,aes(x=TIME_ELAPSED,y=GEN_DIST,size=M,colour=TRANSMISSION_PAIR),alpha=0.8)+
    scale_color_manual(name="", values=c(pal[6], pal[3])) +
    #scale_color_manual(name="",values = c(pal[4],pal[5],pal[2],pal[1]),
    #                  labels=c('Unlinked pair with\nno signal',
    #                           'Unlinked pair with\nfalse signal',
    #                           'Transmission pair\nwith signal',
    #                           'Transmission pair\nno signal')) +
    scale_size(range = c(1,3)) +
    labs(x='Time elapsed (in years)',y='Genetic distance of pair \n(nucleotide substitutions per site)',color="True transmission pair",size='Transmission pair probability\nunder model\n(posterior median)')+
    theme(legend.position='bottom',
          legend.title=element_text(size=rel(0.7)),
          legend.text=element_text(size=rel(0.7)))+
    scale_x_continuous(expand = c(0,0)) +
    #scale_y_continuous(expand = c(0,0),labels = scales::percent, limits = c(0,1)) +
    scale_y_continuous(expand=c(0,0)) + 
    coord_cartesian(xlim=c(0,15),ylim=c(0,0.2)) +
    guides(#color = guide_legend(nrow=1,order=1),
      color = guide_legend(nrow=1,order=1),
      size = guide_legend(nrow=1,order=2),
      shape = guide_legend(nrow = 2,order=3),by.col=T)
  #colour = guide_colorbar(order=2,barwidth=15))
  return(p)
}

make_plot_tpair_violin <- function(po,pal){
  
  # plot prob of being a pair
  p <- ggplot(po, aes(x = TRANSMISSION_PAIR, y = M)) +
    geom_jitter(aes(colour=cat), width = 0.3, height = 0, alpha = 0.7) +
    geom_violin(fill = 'transparent') +
    scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
    #ggsci::scale_color_npg() +
    labs(x = 'Transmission pair', y = '\nTransmission pair probability\n(posterior median)') +
    # scale_colour_manual(name = "", values = c(pal[4],pal[5],pal[2],pal[1]),
    #                     labels=c('Unlinked pair, no signal',
    #                              'Unlinked pair, false signal',
    #                              'Transmission pair, with signal',
    #                              'Transmission pair, no signal')) +
    scale_colour_manual(name="", values=c(pal[2], pal[5], pal[3], pal[10]),
                        labels=c('Unlinked pair, no signal',
                                 'Unlinked pair, false signal',
                                 'Transmission pair, with signal',
                                 'Transmission pair, no signal')) + 
    #guide = guide_legend(override.aes = list(color = 'white'))) +
    theme_bw(base_size=14) +
    theme( legend.position = 'bottom',
           #legend.text = element_blank(),
           legend.margin=margin(t=-10,r=-20,b=0,l=-20))+
    guides(col = guide_legend(nrow = 2),by.col=T)
  return(p)
}

plot_estimated_flows_p <- function(tmp,thrsh){
  
  tmp[, SOURCE:= factor(BIN_COV,levels=c('cat1','cat2'),labels=c('Category\n1','Category\n2'))]
  tmp[, type:= 'Mixture model']
  
  thrsh <- melt(thrsh,id.vars=c('p_pairs','SOURCE'))
  setnames(thrsh,'value','M')
  thrsh[, type:= factor(variable,levels=c('p_truth','p_thrsh'),labels=c('Simulated data (truth)','1.5% distance threshold'))]
  
  tmp <- merge(tmp,thrsh,by=c('SOURCE','p_pairs','type','M'),all=T)
  tmp[, type:= factor(type,levels=c('Simulated data (truth)','Mixture model','1.5% distance threshold'))]
  pal3 <- pal_npg("nrc")(1)
  
  pal1 <- pal_npg("nrc")(4)[1]
  pal2 <- pal_npg("nrc")(4)[2]
  pal4 <- pal_npg("nrc")(4)[4]
  
  p <- ggplot(tmp) +
    geom_bar(aes(x=as.factor(SOURCE),y=M,fill=as.factor(type)),stat='identity', position = position_dodge(width=0.9),alpha=0.5, width=0.8) +
    geom_errorbar(aes(x=as.factor(SOURCE),y=M,ymin=CL, ymax=CU,color=as.factor(type)), width=0.3, position = position_dodge(width=0.9), size=1) +
    scale_fill_manual(name="",values = c('grey50',pal2,pal4),labels=c('Simulated data (truth)','Mixture model','1.5% distance threshold')) +
    scale_colour_manual(name="",values = c('grey50','grey50','grey50'),labels=c('Simulated data (truth)','Mixture model','1.5% distance threshold')) +
    scale_shape_manual(name="Estimated proportions\ndistance threshold",values = c(17),labels=c('1.5% threshold')) +
    theme_bw(base_size=16) + theme(legend.position = "bottom", legend.margin=margin(t=-10)) +
    #axis.text.x=element_text(angle=35,vjust=0.7)) +
    labs(x='\n', y='\nPopulation-level\nsources of infection') +
    coord_cartesian(ylim = c(0,1))+scale_y_continuous(labels = scales::label_percent(accuracy = 1L),breaks=seq(0,1,0.1)) +
    guides(colour = "none",
           fill = guide_legend(ncol = 1),
           by.col=T)
  
  p
  
}

plot_MAE_compare_models <- function(tmp,thrsh,ylim=0.4){
  
  setnames(thrsh,'MAE','value')
  
  tmp[, type:= 'Model estimates']
  thrsh[, type:= '1.5% distance threshold']
  tmp[, type:= factor(type,levels=c('Model estimates','1.5% distance threshold'))]
  tmp[, competing_pairs := round(1/p_pairs,2)]
  thrsh[, competing_pairs := round(1/p_pairs,2)]
  
  p2 <- ggplot(tmp,aes(x=competing_pairs),fill=pal3) +
    geom_line(aes(y=M,color=model),size=1) +
    geom_errorbar(aes(y=M,ymin=CL, ymax=CU,color=model), width=0,size=1) +
    geom_point(aes(y=M,color=model), size=2) +
    geom_line(data=thrsh,aes(y=value,color='grey50',linetype = "1.5% distance threshold"), colour = 'grey50', size=1.2) +
    geom_point(data=thrsh,aes(y=value,color='grey50'),color='grey50', size=2) +
    scale_linetype_manual(name = "", values = c(2),
                          guide = guide_legend(override.aes = list(color = c('grey50')),title.position='top',nrow = 1)) +
    scale_colour_npg() +
    theme_bw(base_size=16) + theme(legend.position = "right",
                                   legend.margin=margin(t=-10))+
    labs(x='Average number of phylogenetically\npossible sources per new infection', y='Mean absolute error \n(flows from sources)  \n',
         color='Model') +
    scale_y_continuous(breaks=seq(0,ylim,0.02)) +
    coord_cartesian(ylim = c(0,ylim))+
    guides(line = guide_legend(nrow = 1,order=2),
           colour = guide_legend(nrow = 3,title.position='top',order=1),
           by.col=T)
}

fit_model_2d_data <- function(do,source_dir,outdir){
  
  ggplot(do, aes(x = TIME_ELAPSED, y = GEN_DIST)) +
    geom_point() +
    geom_density2d()
  
  # try to fit a mixture of gaussians using EMCluster package
  data <- matrix(c(do$TIME_ELAPSED,do$GEN_DIST),ncol=2)
  n_clusters <- 3  # Number of clusters
  init_params <- init.EM(data, nclass = n_clusters)
  gmm <- emcluster(data,init_params,assign.class=T)
  print(gmm)
  summary(gmm)
  
  plot(data, col = gmm$class)
  g <- ggplot(data.table(data)) + geom_point(aes(x=V1,y=V2,col=factor(gmm$class))) + 
    scale_color_discrete() + theme_bw() + labs(x='time elapsed',y='patristic distance',col='cluster membership')
  ggsave(file = paste0(outfile.base,'-BNMM_classifications.png'),g, w = 7, h = 5)
  
  png(file=paste0(outfile.base,'-plotEM_BNMM_classifications.png'), width=400, height=350)
  plotem(gmm, data, main = NULL, xlab = NULL, ylab = NULL)
  dev.off()
  
  #plot2d(data, emobj = gmm, k = 3, color.pch = 1,
  #       append.BN = TRUE)
  
  sigma <- LTSigma2variance(gmm$LTSigma)
  for (k in 1:gmm$nclass) {
    plotBN(mu[, k], sigma[, , k], lty = lty, col = col, lwd = lwd, 
           cex = cex)
  }
  
  # store parameters:
  stan_data$K <- n_clusters
  stan_data$mu <- t(gmm$Mu)
  stan_data$gmm_pi <- gmm$pi
  #covariances <- gmm$LTSigma
  stan_data$Sigma <- array(sigma[, , ], dim = c(ncol(data), ncol(data), n_clusters))
  stan_data$Sigma <- aperm(stan_data$Sigma, c(3, 1, 2))
  data_matrix <- array(NA_real_, dim = c(n_clusters,ncol(data),ncol(data)))
  
  options(mc.cores=parallel::detectCores())
  warmup <- 1000
  
  # model using cmdstan
  model = rstan::stan_model(file.path(source_dir,'stan_model_files','sens_analysis_2DMM_bg.stan'))
  model = rstan::stan_model(file.path(source_dir,'stan_model_files','sens_analysis_2DGMM_bg.stan'))
  model = rstan::stan_model(file.path(source_dir,'stan_model_files','sens_analysis_mixlognorm_bg.stan'))
  
  stan_data <- list()
  stan_data$N <- nrow(do)
  stan_data$D <- 2
  stan_data$K <- 2
  stan_data$y <- matrix(NA,nrow=stan_data$N,ncol=2)
  stan_data$y[,1] <- do$TIME_ELAPSED
  stan_data$y[,2] <- do$GEN_DIST
  
  init_fun <- function() {
    list(
      theta = rep(0.5, 2),
      #mu1 = c(1, 0.02),
      #mu2 = c(4, 0.1),
      mu1_base = log(1),
      mu2_base = log(5),
      mu1 = 0,
      mu2 = 0,
      #Sigma1 = diag(2),
      #Sigma2 = diag(2)
      sigma1 = 1,
      sigma2 = 1
    )
  }
  
  fit = rstan::sampling(model,chains=3,data=stan_data,
                        warmup=500,iter=2000,
                        control=list(adapt_delta=.99),
                        #init = init_fun
  )
  #saveRDS(fit,file=paste0(outdir,'stanfit_2D_bg_dist.rds'))
  saveRDS(fit,file=paste0(outdir,'stanfit_MMlnorm_bg_dist.rds'))
  
  # Print the results
  print(fit)
  
  #	examine neff and rhat
  #fit.target.pars <- c('theta','mu1','mu2','Sigma1','Sigma2')
  #fit.target.pars <- c('theta','mu1','mu2','L_Sigma1','L_Sigma2')
  #fit.target.pars <- c('theta','mu1_base','mu2_base','sigma1','sigma2')
  fit.target.pars <- c('theta','mu','L')
  summary <- rstan::monitor(rstan::extract(fit, pars=c(fit.target.pars), permuted = FALSE, inc_warmup = TRUE))
  print(summary,probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
  
  # get samples from the chains
  samples <- rstan::extract(fit, inc_warmup = FALSE)
  #saveRDS(samples,file=paste0(outdir, 'samples_','2D_bg_dist','.rds'))
  saveRDS(samples,file=paste0(outdir, 'samples_','MMlnorm_bg_dist','.rds'))
  
  #	traces
  color_scheme_set("mix-blue-red")
  p <- rstan::traceplot(fit, pars=c(fit.target.pars,'lp__'),inc_warmup=FALSE, ncol = 1)
  #pdf(file=paste0(outdir, 'traces_','2D_bg_dist',".pdf"), w=10, h=20)
  pdf(file=paste0(outdir, 'traces_','MMlnorm_bg_dist',".pdf"), w=10, h=20)
  print(p)
  dev.off()
  
  # Extract and analyze the results
  extract_fit <- extract(fit)
  theta <- extract_fit$theta
  mu1 <- extract_fit$mu1
  mu2 <- extract_fit$mu2
  #Sigma1 <- extract_fit$Sigma1
  #Sigma2 <- extract_fit$Sigma2
  Sigma1 <- extract_fit$L_Sigma1
  Sigma2 <- extract_fit$L_Sigma2
  
  mu1_background <- extract_fit$mu1_background
  mu2_background <- extract_fit$mu2_background
  mu2_background <- melt(mu2_background)
  setnames(mu2_background,c('Var2','value'),c('pair_id','mu2_bg'))
  mu2_background <- mu2_background[, list(mu2_bg=median())]
  sigma1 <- extract_fit$sigma1
  sigma2 <- extract_fit$sigma2
  
  sim <- data.table(
    TIME_ELAPSED=do$TIME_ELAPSED,
    GEN_DIST =rlnorm(nrow(do),median(mu2_background),sigma2)
  )
  
  col.pal <- pal_npg('nrc')(4)
  # Plot the data and the estimated means
  ggplot(do) +
    geom_point(aes(x=TIME_ELAPSED,y=GEN_DIST),col = col.pal[1]) +
    geom_point(aes(x=median(mu1[,1]),y=median(mu1[,2])),col = col.pal[2]) +
    geom_point(aes(x=median(mu2[,1]),y=median(mu2[,2])),col = col.pal[3])
  points(apply(mu1, 2, mean), col = col.pal[1], pch = 19, cex = 2)
  points(apply(mu2, 2, mean), col = col.pal[2], pch = 19, cex = 2)
  
  # simulate data
  require(mnormt)
  x     <- seq(0, 16, 0.1)
  y     <- seq(0, 0.15, 0.01)
  mu1    <- c(median(mu1[,1]), median(mu1[,2]))
  sigma1 <- matrix(c(median(Sigma1[,1,1]), median(Sigma1[,1,2]),
                     median(Sigma1[,2,1]), median(Sigma1[,2,2])), nrow = 2)
  sigma1 <- (sigma1 + t(sigma1)) / 2 # tmp fix to make covariance matrix symmetric
  f     <- function(x, y) dmnorm(cbind(x, y), mu1, sigma1)
  z     <- outer(x, y, f)
  
  contour(x, y, z)
  
  mu2    <- c(median(mu2[,1]), median(mu2[,2]))
  sigma2 <- matrix(c(median(Sigma2[,1,1]), median(Sigma2[,1,2]),
                     median(Sigma2[,2,1]), median(Sigma2[,2,2])), nrow = 2)
  sigma2 <- (sigma2 + t(sigma2)) / 2 # tmp fix to make covariance matrix symmetric
  f     <- function(x, y) dmnorm(cbind(x, y), mu2, sigma2)
  z     <- outer(x, y, f)
  
  contour(x, y, z)
  
  clu1 <- mvrnorm(1000,mu1,sigma1)
  clu2 <- mvrnorm(1000,mu2,sigma2)
  colnames(clu1) <- c('time','gen_dist')
  colnames(clu2) <- c('time','gen_dist')
  ggplot() +  geom_point(data=do,aes(x=TIME_ELAPSED,y=GEN_DIST),col = col.pal[1]) +
    geom_point(data=clu1,aes(x=time,y=gen_dist),col=col.pal[2]) +
    geom_point(data=clu2,aes(x=time,y=gen_dist),col=col.pal[3]) +
    coord_cartesian(xlim=c(0,16),ylim=c(0,0.2)) +
    theme_bw()
  
  # summarise posteriors
  dmu <- data.table(reshape2::melt(samples$mu))
  setnames(dmu,c('iterations','Var2'),c('iter','gp'))
  dmu <- dmu[, list(mu=quantile(value,probs=0.5)),by='gp']
  dmu <- c(dmu[1,2]$mu,dmu[2,2]$mu)
  
  dsig <- data.table(reshape2::melt(samples$Sigma))
  setnames(dsig,c('iterations','Var2','Var3'),c('iter','gp1','gp2'))
  dsig <- dsig[, list(sig=quantile(value,probs=0.5)),by=c('gp1','gp2')]
  dsig <- matrix(c(dsig[1,'sig']$sig,dsig[2,'sig']$sig,dsig[3,'sig']$sig,dsig[4,'sig']$sig),ncol=2)
  
  bivn <- mvrnorm(5000, mu = dmu, Sigma = dsig )  # from Mass package
  bivn.kde <- kde2d(bivn[,1], bivn[,2], n = 50)
  #image(bivn.kde)       # from base graphics package
  contour(bivn.kde, add = TRUE)     # from base graphics package
  
  
  ## other option - fit 2D density
  # Kernel density estimation
  dens <- kde2d(do$TIME_ELAPSED, do$GEN_DIST, n = 50)
  
  # Extract the density estimate for use in another model
  x_grid <- dens$x
  y_grid <- dens$y
  z_density <- dens$z
  
  
  # plot 2D BNMM https://maggielieu.com/2017/03/21/multivariate-gaussian-mixture-model-done-properly/
  params=extract(fit)
  #density plots of the posteriors of the mixture means
  par(mfrow=c(2,2))
  plot(density(params$mu[,1,1]), ylab='', xlab='mu[1]', main='')
  lines(density(params$mu[,1,2]), col=rgb(0,0,0,0.7))
  lines(density(params$mu[,1,3]), col=rgb(0,0,0,0.4))
  lines(density(params$mu[,1,4]), col=rgb(0,0,0,0.1))
  abline(v=c(0), lty='dotted', col='red',lwd=2)
  
  plot(density(params$mu[,2,1]), ylab='', xlab='mu[2]', main='')
  lines(density(params$mu[,2,2]), col=rgb(0,0,0,0.7))
  lines(density(params$mu[,2,3]), col=rgb(0,0,0,0.4))
  lines(density(params$mu[,2,4]), col=rgb(0,0,0,0.1))
  abline(v=c(7), lty='dotted', col='red',lwd=2)
  
  plot(density(params$mu[,3,1]), ylab='', xlab='mu[3]', main='')
  lines(density(params$mu[,3,2]), col=rgb(0,0,0,0.7))
  lines(density(params$mu[,3,3]), col=rgb(0,0,0,0.4))
  lines(density(params$mu[,3,4]), col=rgb(0,0,0,0.1))
  abline(v=c(3), lty='dotted', col='red',lwd=2)
  
}

