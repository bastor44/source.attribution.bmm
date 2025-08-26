require(ggplot2)
require(ggsci)
require(scales)
require(data.table)
require(kableExtra)
require(ggpubr)
require(dplyr)

args <- list(
  indir = '/Users/bethanyastor/Library/CloudStorage/OneDrive-ImperialCollegeLondon/Research Project/source.attribution.bmm', 
  out.dir = '/Users/bethanyastor/Library/CloudStorage/OneDrive-ImperialCollegeLondon/Research Project/source.attribution.bmm/figures', 
  job.name = 'simulations_500truepairs_srcage_exclTE16_subsample', 
  stan.model.unif = 'mm_bgUnif_piGP_221027', 
  stan.model.gmm = 'mm_bgGMM_piGP_240711'
)


source('R/functions_plotting.R')

#tab <- data.table(F=list.dirs('/Users/alexb/Box Sync/Roadmap/source_attribution'))
tab <- data.table(F=list.dirs('/Users/bethanyastor/Library/CloudStorage/OneDrive-ImperialCollegeLondon/Research Project/source.attribution.bmm/out_simulated'))
tab <- tab[ grepl('simulations_500truepairs_srcage_exclTE16_subsample', F) ,]
# tab[, model:= 'vanilla']
# tab[grepl('piReg',F), model:= 'covariate']
# tab[grepl('piGP',F), model:= 'hsgp']
tab[grepl('bgUnif', F), model:='unif']
tab[grepl('bgGMM', F), model:='gmm']
tab[grepl('2stage', F), model:='gmm2']
tab[, prop:= as.numeric(gsub('^.*_subsample([0-9]+)pct-([0-9]+).+','\\1',F))]
tab[prop==0, prop:= 100]
#regex <- '^.*attribution/([A-Za-z0-9_]+)-([A-Za-z0-9_]+)-([0-9]+)'
#regex <- '^.*out_simulated/([A-Za-z0-9_]+)-([A-Za-z0-9_]+)-([0-9]+)'
#regex <- '^(mm_bg[A-Za-z]_piGP_d+[a-z]-simulations_500truepairs_srcage_exclTE16_subsample50pct)'
#tab[, outfile.base:= file.path(F, paste0(gsub(regex, '\\1', F), '-', gsub(regex, '\\2', F)))]
tab

po <- NULL
src <- NULL
thrsh <- NULL
for(i in 1:nrow(tab)){
  tmp <- readRDS(file=paste0(tab[i,outfile.base],'-MAE_5yrage_all_pairs','.RDS'))
  tmp[, p_pairs:= tab[i,prop]/100]
  tmp[p_pairs==1, p_pairs:= 0.08]
  
  tmp[, model:= tab[i,model]]
  
  tmp2 <- readRDS(file = paste0(tab[i,outfile.base],'.rds'))
  tmp3 <- tmp2[, list(N_sources=length(unique(FROM_ID))),by='TO_ID']
  tmp3 <- tmp3[, list(p_pairs=tmp[, p_pairs],
                      N_sources=mean(N_sources))]
  tmp3[, model:= tab[i,model]]
  
  sim_scenario <- readRDS(file=paste0(tab[i,outfile.base],'.RDS'))
  sim_scenario[, FROM_AGE_FIVE:= cut(FROM_AGE,breaks=seq(15,75,5),include.lowest=T,right=F,
                                     labels=c('15-20','20-25','25-30','30-35','35-40','40-45','45-50','50-55','55-60','60-65','65-70','70-75'))]
  
  sim_scenario[, STAGE_thrsh:= 0]
  sim_scenario[GEN_DIST<0.015, STAGE_thrsh:= 1]
  sim_scenario[, p_pairs:= tab[i,prop]/100]
  sim_scenario[p_pairs==1, p_pairs:= 0.08]
  
  tmp4 <- sim_scenario[, list(N=length(FROM_ID[STAGE_thrsh==1]),
                              N_true=length(FROM_ID[TRANSMISSION_PAIR=='Yes'])), by=c('p_pairs','FROM_AGE_FIVE')]
  tmp4 <- tmp4[, list(FROM_AGE_FIVE=FROM_AGE_FIVE,
                      pct_thrsh=N/nrow(sim_scenario[STAGE_thrsh==1]),
                      pct_true=N_true/nrow(sim_scenario[TRANSMISSION_PAIR=='Yes'])), by=c('p_pairs')]
  
  tmp4[, thrsh_error:= abs(pct_true - pct_thrsh)]
  tmp4 <- tmp4[, list(MAE=mean(thrsh_error)),by='p_pairs']
  tmp4[, model:= tab[i,model]]
  
  po <- rbind(po,tmp)
  src <- rbind(src,tmp3)
  thrsh <- rbind(thrsh,tmp4)
}



# po[, model:= factor(model,levels=c('vanilla','covariate','hsgp'),
#                     labels=c('No covariates',
#                              'Discrete covariate',
#                              'Continuous age modelled\nwith random function'))]
po[, model:=factor(model, levels=c('unif', 'gmm', 'gmm_rmsig'), 
                   labels=c('Uniform background', 
                            'GMM background', 
                            'GMM with signal removed'))]
# thrsh[, model:= factor(model,levels=c('vanilla','covariate','hsgp'),
#                        labels=c('No covariates',
#                                 'Discrete covariate',
#                                 'Continuous age modelled\nwith random function'))]
thrsh[, model:=factor(model, levels=c('unif', 'gmm', 'gmm_rmsig'),
                      labels=c('Uniform background', 
                               'GMM background', 
                               'GMM with signal removed'))]

g2 <- plot_MAE_compare_models(po,subset(thrsh,model=='No covariates'),ylim=0.08)

ggsave(g2,filename=file.path(args$out.dir,paste0(args$job.name,'-MAE_models_competingpairs_agesrc_5yr.png')),w=9,h=6)

# save table
po[, L:= paste0(round(M*100,1),'% [',round(CL*100,1),'-',round(CU*100,1),'%]')]
thrsh[, mae_thrsh:= paste0(round(MAE*100,1),'%')]
po <- dcast(po,p_pairs~model,value.var='L')
po <- merge(subset(thrsh,select=c('p_pairs','mae_thrsh'),model=='No covariates'),
            po,
            by=c('p_pairs'))
po <- po[order(-p_pairs)]
saveRDS(po,file=file.path(args$out.dir,paste0(args$job.name,'-compare_MAE_ages_source_5yr.RDS')))

# plot ages ----
#outfile.base <- tab[model=='hsgp' & prop== 50, outfile.base]
outfile.base <- tab[model=='gmm_rmsig', outfile.base][1]

sim_scenario <- readRDS(file = paste0(outfile.base,'.rds'))
sim_scenario[, p_pairs:= 0.5]
sim_scenario[, TRSM:= 'Non-transmission pair']
sim_scenario[TRANSMISSION_PAIR=='Yes', TRSM:= 'True transmission pair']

pal <- pal_npg("nrc")(2)
g1 <- ggplot(sim_scenario) +
  geom_point(aes(x=TO_AGE,y=FROM_AGE,colour=TRANSMISSION_PAIR),alpha=0.7) +
  labs(x='Age of recipient on estimated\ninfection date of recipient',
       y='Age of source on estimated\ninfection date of recipient',
       colour = 'Transmission pair') +
  scale_x_continuous(breaks=seq(15,80,10),labels=seq(15,80,10)) +
  scale_y_continuous(breaks=seq(15,80,10),labels=seq(15,80,10)) +
  scale_colour_manual(values=c(pal[1],pal[2])) +
  theme_bw(base_size=16) +
  theme(strip.background=element_blank(),
        legend.position='bottom',
        legend.margin=margin(t=-10)) +
  guides(colour = guide_legend(nrow = 2,title.position='top'))
g1

g <- ggarrange(g1,g2+theme(legend.position='bottom'),ncol=2,align='hv',labels=c('AUTO'),font.label = list(size=18))
ggsave(g,filename=file.path(args$out.dir,paste0(args$job.name,'-MAE_models_competingpairs_ages_agesrc_5yr_bluepairs_gmm.png')),w=12.5,h=6)
ggsave(g,filename=file.path(args$out.dir,paste0(args$job.name,'-MAE_models_competingpairs_ages_agesrc_5yr_bluepairs_gmm.pdf')),w=12.5,h=6)
