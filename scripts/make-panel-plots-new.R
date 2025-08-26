require(ggplot2)
require(ggsci)
require(scales)
require(data.table)
require(kableExtra)
require(ggpubr)
require(dplyr)
require(gridExtra)

args <- list(
  source_dir = '/Users/bethanyastor/Library/CloudStorage/OneDrive-ImperialCollegeLondon/Research Project/source.attribution.bmm',
  indir = '/Users/bethanyastor/Library/CloudStorage/OneDrive-ImperialCollegeLondon/Research Project/source.attribution.bmm',
  # uniform bg
  #outdir = '/Users/bethanyastor/Library/CloudStorage/OneDrive-ImperialCollegeLondon/Research Project/source.attribution.bmm/out_simulated/mm_bgUnif_piGP_221027b-simulations_500truepairs_srcage_exclTE16_subsample50pct-1052613'
  #outdir = '/Users/bethanyastor/Library/CloudStorage/OneDrive-ImperialCollegeLondon/Research Project/source.attribution.bmm/out_simulated/mm_bgUnif_piGP_221027b-simulations_450truepairs_srcage_exclTE16_subsample14pct-405979'
  # gmm (all) bg
  #outdir = '/Users/bethanyastor/Library/CloudStorage/OneDrive-ImperialCollegeLondon/Research Project/source.attribution.bmm/out_simulated/mm_bgGMM_piGP_240711-simulations_500truepairs_srcage_exclTE16_subsample50pct-53228'
  #outdir = '/Users/bethanyastor/Library/CloudStorage/OneDrive-ImperialCollegeLondon/Research Project/source.attribution.bmm/out_simulated/mm_bgGMM_piGP_240711-simulations_450truepairs_srcage_exclTE16_subsample14pct-323432'
  # gmm (bs) bg
  outdir = '/Users/bethanyastor/Library/CloudStorage/OneDrive-ImperialCollegeLondon/Research Project/source.attribution.bmm/out_simulated/mm_bgGMM_piGP_240711-simulations_500truepairs_srcage_exclTE16_subsample50pct-370047-2stage'
  #outdir = '/Users/bethanyastor/Library/CloudStorage/OneDrive-ImperialCollegeLondon/Research Project/source.attribution.bmm/out_simulated/mm_bgGMM_piGP_240711-simulations_450truepairs_srcage_exclTE16_subsample14pct-847679-2stage'
)
#outfile.base = file.path(args$outdir, 'mm_bgUnif_piGP_221027b-simulations_450truepairs_srcage_exclTE16_subsample14pct')
outfile.base <- file.path(args$outdir, 'mm_bgGMM_piGP_240711-simulations_500truepairs_srcage_exclTE16_subsample50pct')

source(file.path(args$indir, 'R/functions_plotting.R'))
source(file.path(args$indir, 'R/functions_simulation_scenarios.R'))

infile.stanin <- list.files(args$outdir, pattern=paste0('_stanin.RData$'), recursive=TRUE)[1]
stopifnot(length(infile.stanin)>0)
stopifnot(length(infile.stanin)<=1)
tmp <- load(file.path(args$outdir, infile.stanin))
stopifnot(c('args','stan_data')%in%tmp)

tmp <- paste0(outfile.base,'-fitted_stan_model.rds')
cat("\n Read fitted dat ", tmp , "\n")
model_fit <- readRDS(file = tmp)


# load data
sim_scenario <- readRDS(file = paste0(outfile.base,'.rds'))
sim_scenario <- filter(sim_scenario, GEN_DIST >= 0)

# define groupings for model ----
sim_scenario[, FROM_AGE_GP:= cut(FROM_AGE,breaks=c(15,30,40,50,60,100),include.lowest=T,right=F,
                                 labels=c('[15-30)','[30-40)','[40-50)','[50-60)','[60+)'))]
sim_scenario[, TO_AGE_GP:= cut(TO_AGE,breaks=c(15,30,40,50,60,100),include.lowest=T,right=F,
                               labels=c('[15-30)','[30-40)','[40-50)','[50-60)','[60+)'))]

sim_scenario[,FROM_AGE_GP_2:= factor(FROM_AGE_GP,labels=c(1,2,3,4,5))]
sim_scenario[,TO_AGE_GP_2:= factor(TO_AGE_GP,labels=c(1,2,3,4,5))]

if(grepl('GP',args$stanModelFile)){
  # modify ages
  # max age 90
  sim_scenario[TO_AGE>90, TO_AGE:=90]
  sim_scenario[FROM_AGE>90, FROM_AGE:=90]
  
  sim_scenario[,FROM_AGE_INT:= round(FROM_AGE)]
  sim_scenario[,TO_AGE_INT:= round(TO_AGE)]
  
  # get unique pairs ages
  tmp <- data.table(unique(cbind(sim_scenario$FROM_AGE_INT,sim_scenario$TO_AGE_INT)))
  tmp[, IDX_UNIQUE_PAIR := seq_len(nrow(tmp))]
  setnames(tmp, "V1", "FROM_AGE_INT")
  setnames(tmp, "V2", "TO_AGE_INT")
  sim_scenario <- merge(sim_scenario, tmp, by = c("FROM_AGE_INT","TO_AGE_INT"))
  
  # re-order by ages
  sim_scenario <- sim_scenario[order(FROM_AGE_INT,TO_AGE_INT),]
  sim_scenario <- data.table(sim_scenario)
  sim_scenario[, PAIR_ID := seq(1,nrow(sim_scenario)) ]
  
}

if(grepl('GP',args$stanModelFile)){
  sim_scenario[, FROM_AGE_STD:= (FROM_AGE_INT - mean(FROM_AGE_INT))/sd(FROM_AGE_INT)]
  sim_scenario[, TO_AGE_STD:= (TO_AGE_INT - mean(TO_AGE_INT))/sd(TO_AGE_INT)]
  tmp <- data.table(FROM_AGE_STD = sort(unique(sim_scenario$FROM_AGE_STD)))
  tmp[, IDX_UNIQUE_FROM_AGE := seq_len(nrow(tmp))]
  sim_scenario <- merge(sim_scenario,tmp,by='FROM_AGE_STD',all.x=T)
  sim_scenario <- sim_scenario[order(PAIR_ID)]
}

sim_scenario[,FROM_AGE_INT:= round(FROM_AGE)]
sim_scenario[,TO_AGE_INT:= round(TO_AGE)]

# create index for source/recipient categories
sim_scenario[, FROM_BIN_COV:= as.numeric(factor(BIN_COV,levels=c('cat1','cat2')))] # set recipient binary covariate to same as source
sim_scenario[, TO_BIN_COV:= as.numeric(factor(BIN_COV,levels=c('cat1','cat2')))]


# fig A (posterior prob by size) ----
po <- posterior_prob_pair(model_fit,sim_scenario,dps)

g1 <- make_plot_simulated_data_colour_prob_tpair_sizes(po, dps, sim_scenario, outfile.base)
g1
ggsave(filename=paste0(outfile.base, '-posterior-sizes.png'), g1, h=8, w=10)
ggsave(filename=paste0(outfile.base, '-posterior-sizes.pdf'), g1, h=8, w=10)

# fig B (violin plot) ----
# false signal
sim_scenario[, d_TSeqT:= round(TIME_ELAPSED,2)]
dps[, d_TSeqT:= round(d_TSeqT,2)]
sim_scenario <- merge(sim_scenario,subset(dps,select=c('d_TSeqT','q2.5','q97.5')), by='d_TSeqT',all.x=T)

# flag the false positives
sim_scenario[TRANSMISSION_PAIR=='No', cat:= ifelse(GEN_DIST>=q2.5 & GEN_DIST<=q97.5,'false_pos',
                                                   'true_neg')]
sim_scenario[TRANSMISSION_PAIR=='Yes', cat:= ifelse(GEN_DIST<q2.5 | GEN_DIST>q97.5,'false_neg',
                                                    'true_pos')]
sim_scenario[, cat:=factor(cat,levels=c('true_neg','false_pos','true_pos','false_neg'),
                           labels=c('Unlinked pair with\nno signal',
                                    'Unlinked pair with\nfalse signal',
                                    'Transmission pair\nwith signal',
                                    'Transmission pair\nno signal'))]

tmp <- sim_scenario[, list(p=length(PAIR_ID)/nrow(sim_scenario)), by=c('cat')]

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

po[, SOURCE:= 'Source category 1']
po[BIN_COV=='cat2', SOURCE:= 'Source category 2']

#pal <- pal_npg("nrc")(5)
pal <- pal_brewer(palette='Paired')(10)
# plot prob of being a pair
g2 <- make_plot_tpair_violin(po,pal)
g2

g2_leg <- cowplot::get_legend(g2)
ggsave(filename=paste0(outfile.base, '-violin-posterior-prob.png'), g2, h=5, w=6)

# fig C (PAF) ----
po <- model_fit$draws(inc_warmup = FALSE,
                      format = 'draws_df',
                      variables = 'pflows_from'
)
po <- as.data.table(po)
setnames(po, colnames(po), gsub('^\\.','',colnames(po)))

po <- melt(po, id.vars = c('chain','iteration','draw'))

po[, FROM_AGE_GP_2 := as.factor(gsub('pflows_from\\[([0-9]+)\\]','\\1',as.character(variable)))]
po <- merge(po, src_map, by = 'FROM_AGE_GP_2')
po <- po[,
         list( q = quantile(value, probs = c(0.5, 0.025, 0.25, 0.75, 0.975) ),
               stat = c('M','CL','IL', 'IU', 'CU')
         ),
         by = c('FROM_AGE_GP')
]
po <- dcast.data.table(po, FROM_AGE_GP~stat, value.var = 'q')
tmp <- sim_scenario[TRANSMISSION_PAIR == 'Yes',
                    list(value = length(PAIR_ID)),
                    by = 'FROM_AGE_GP'
]

tmp <- model_fit$draws(inc_warmup = FALSE,
                       format = 'draws_df',
                       variables = 'true_flows'
)
tmp <- as.data.table(tmp)
setnames(tmp, colnames(tmp), gsub('^\\.','',colnames(tmp)))
tmp <- melt(tmp, id.vars = c('chain','iteration','draw'))

#age group strat
tmp[, FROM_AGE_GP_2 := as.factor(gsub('true_flows\\[([0-9]+)\\]','\\1',as.character(variable)))]
tmp <- merge(tmp, src_map, by = 'FROM_AGE_GP_2')
tmp <- tmp[,
           list( true_p = quantile(value, probs = c(0.5) ),
                 stat = c('M')
           ),
           by = c('FROM_AGE_GP')
]
po <- merge(po, subset(tmp, select = c('FROM_AGE_GP','true_p')), by = 'FROM_AGE_GP')

plot_PAF_age_no_thrsh(po, plot='-PAF_age_all_pairs', outfile.base)


po[, SOURCE:= factor(FROM_AGE_GP, levels=c("[15-30)", "[30-40)", "[40-50)", "[50-60)", "[60+)"),
                            labels=c("15-29", "30-39", "40-49", "50-59", "60+"  ))]
tmp <- melt(po,id.vars=c('FROM_AGE_GP','SOURCE','CL','CU','IL','IU','true_p'),value.vars=c('M'))

g3 <- ggplot(tmp) +
  geom_bar(data=subset(tmp, variable='M'), aes(x=as.factor(SOURCE), y=true_p, colour='grey50', fill='white'),stat='identity', linetype=2) +
  geom_bar(aes(x=as.factor(SOURCE), y=value, color=pal[9], fill=pal[9]), stat='identity', alpha=0.7) +
  geom_errorbar(aes(x=as.factor(SOURCE), ymin=CL, ymax=CU), color='grey30', width=0.3) +
  scale_fill_manual(name="", values=c(pal[9], 'white'), labels=c('Estimated flows', '')) +
  scale_color_manual(name="", values=c('grey50', 'grey50'), labels=c('True flows', '')) +
  theme_bw(base_size=12) + theme(legend.position = "bottom") +
  labs(x='\n Source category of transmitter', y='Proportion of transmission \npairs attributible \nto source age') +
  coord_cartesian(ylim = c(0,1))+scale_y_continuous(labels = scales::label_percent(accuracy = 1L),breaks=seq(0,1,0.1))
g3


# combine ----
g <- grid.arrange(arrangeGrob(ggarrange(g1+ theme(legend.title=element_text(size=rel(1))),
                                        ggarrange(g2 + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()), g3,
                                                  labels=c('B','C'),font.label = list(size=18), ncol = 2, align='hv'),
                                        labels=c('A','',''),font.label = list(size=18),ncol=1,heights=c(0.55,0.4,0.05))))

ggsave(g,filename=paste0(outfile.base,'-4panelplot.png'),w=11,h=10)
ggsave(g,filename=paste0(outfile.base,'-4panelplot.pdf'),w=11,h=10)
