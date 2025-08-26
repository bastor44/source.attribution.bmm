## preamble ----
require(data.table)  # data mangling
require(bayesplot)
require(hexbin)
require(knitr)
require(ggplot2)
require(ggpubr)
require(rstan)  # run Stan from R
require(cmdstanr)
require(ggsci)
require(scales)
require(grid)
require(viridis)
require(loo)
require(boot)
require(patchwork)
require(RColorBrewer)

if (1)
{
  args <- list(
    source_dir = '/Users/bethanyastor/Library/CloudStorage/OneDrive-ImperialCollegeLondon/Research Project/source.attribution.bmm',
    indir = '/Users/bethanyastor/Library/CloudStorage/OneDrive-ImperialCollegeLondon/Research Project/source.attribution.bmm', 
    #outdir = '/Users/bethanyastor/Library/CloudStorage/OneDrive-ImperialCollegeLondon/Research Project/source.attribution.bmm/out_Amsterdam/mm_bgUnif_piGP_221027b-agegps_sensanalysis_210216_MSM-1115611',
    #outdir = '/Users/bethanyastor/Library/CloudStorage/OneDrive-ImperialCollegeLondon/Research Project/source.attribution.bmm/out_Amsterdam/mm_bgGMM_piGP_240711-agegps_sensanalysis_210216_MSM-632303', 
    #outdir = '/Users/bethanyastor/Library/CloudStorage/OneDrive-ImperialCollegeLondon/Research Project/source.attribution.bmm/out_Amsterdam/mm_bgUnif_piGP_221027b-agegps_sensanalysis_210216_MSM-342783',
    outdir = '/Users/bethanyastor/Library/CloudStorage/OneDrive-ImperialCollegeLondon/Research Project/source.attribution.bmm/out_Amsterdam/mm_bgUnif_piGP_221027b-agegps_sensanalysis_210216_MSM-638302',
    #outdir = '/Users/bethanyastor/Library/CloudStorage/OneDrive-ImperialCollegeLondon/Research Project/source.attribution.bmm/out_Amsterdam/mm_bgGMM_piGP_240711-agegps_sensanalysis_210216_MSM-834222',
    #outdir = '/Users/bethanyastor/Library/CloudStorage/OneDrive-ImperialCollegeLondon/Research Project/source.attribution.bmm/out_Amsterdam/mm_bgGMM_piGP_240711-agegps_sensanalysis_210216_MSM-828325-rm-signal',
    clock_model = '/Users/bethanyastor/Library/CloudStorage/OneDrive-ImperialCollegeLondon/Research Project/source.attribution.bmm/molecular_clock', 
    stanModelFile = 'mm_bgUnif_piGP_221027b',
    #stanModelFile = 'mm_bgGMM_piGP_240711',
    scenario = 15,
    reps = 1,
    rep = 1,
    simulate_data = T,
    job_tag = 'agegps_sensanalysis_210216_MSM',
    weights = T,
    tpair <- 'tpair_prob_w'
  )
}

## command line parsing if any
args_line <-  as.list(commandArgs(trailingOnly=TRUE))
if(length(args_line) > 0)
{
  stopifnot(args_line[[1]]=='-source_dir')
  stopifnot(args_line[[3]]=='-stanModelFile')
  stopifnot(args_line[[5]]=='-indir')
  stopifnot(args_line[[7]]=='-outdir')
  stopifnot(args_line[[9]]=='-job_tag')
  stopifnot(args_line[[11]]=='-scenario')
  stopifnot(args_line[[13]]=='-rep')
  stopifnot(args_line[[15]]=='-weights')
  
  args <- list()
  args[['source_dir']] <- args_line[[2]]
  args[['stanModelFile']] <- args_line[[4]]
  args[['indir']] <- args_line[[6]]
  args[['outdir']] <- args_line[[8]]
  args[['job_tag']] <- args_line[[10]]
  args[['scenario']] <- as.integer(args_line[[12]])
  args[['rep']] <- as.integer(args_line[[14]])
  args[['weights']] <- as.integer(args_line[[16]])
}
args

outfile.base <- file.path(args$outdir, paste0(args$stanModelFile, '-',args$job_tag))
#outfile.base <- file.path(args$outdir, paste0(args$stanModelFile, '-', args$job_tag, '-rm-signal'))

replicate <- args$rep
cat(paste0("rep ", replicate))

## load functions
source(file.path(args$source_dir, 'R', 'functions_simulation_scenarios.R'))
source(file.path(args$source_dir, 'R', 'functions_plotting.R'))
source(file.path(args$source_dir, 'R', 'component_EM.R'))


cat(" \n --------------------------------  with arguments -------------------------------- \n")

## read stanin
cat('\nReading Stan input data...')
infile.stanin <- paste0(outfile.base, '_stanin.RData')
stopifnot(length(infile.stanin)>0)
stopifnot(length(infile.stanin)<=1)
tmp <- load(infile.stanin)
stopifnot(c('args','stan_data')%in%tmp)


# load pairs data ----
#pairs <- fread(file.path('data_Ams',args$analysis,paste0(args$trsm,"_anonymised_pairs.csv")))

pairs <- fread(file.path(args$indir, 'data_Ams',args$analysis,paste0(args$trsm,"_anonymised_pairs.csv")))
tmp <- pairs[, list(N=length(unique(FROM_SEQUENCE_ID))),by='TO_SEQUENCE_ID']

do <- pairs
do[,LOG_TIME_ELAPSED:= log(TIME_ELAPSED)]

# quantiles of time elapsed
quantile(do$TIME_ELAPSED,probs=c(0.5,0.025,0.975))

# for how many pairs is only one direction possible?
dp <- subset(do,select=c('FROM_SEQUENCE_ID','TO_SEQUENCE_ID'))
dp[, FROM_SEQUENCE_ID:= as.numeric(FROM_SEQUENCE_ID)]
dp[, TO_SEQUENCE_ID:= as.numeric(TO_SEQUENCE_ID)]
dp <- dp[, list(pairs= paste0(min(FROM_SEQUENCE_ID,TO_SEQUENCE_ID),'-', max(FROM_SEQUENCE_ID,TO_SEQUENCE_ID))),by=c('FROM_SEQUENCE_ID','TO_SEQUENCE_ID')]
cp <- subset(dp,select=c('pairs'))
cat(paste0('proportion of IDs in which only one direction is possible: ', round(nrow(unique(cp))/nrow(dp)*100,0),'%'))

# define groupings for model ----
do[,FROM_AGE_INT:= round(FROM_AGE)]
do[,TO_AGE_INT:= round(TO_AGE)]
do[,FROM_AGE_GP_2:= factor(FROM_AGE_GP,labels=c(1,2,3,4,5))]
do[,TO_AGE_GP_2:= factor(TO_AGE_GP,labels=c(1,2,3,4,5))]

# get unique time elapsed

tmp <- data.table(TIME_ELAPSED = sort(unique(do$TIME_ELAPSED)))
tmp[, IDX_UNIQUE_TE := seq_len(nrow(tmp))]
do <- merge(do, tmp, by = 'TIME_ELAPSED')
do <- do[order(PAIR_ID),]

# get unique pairs ages
tmp <- data.table(unique(cbind(do$FROM_AGE_INT,do$TO_AGE_INT)))
tmp[, IDX_UNIQUE_PAIR := seq_len(nrow(tmp))]
setnames(tmp, "V1", "FROM_AGE_INT")
setnames(tmp, "V2", "TO_AGE_INT")
do <- merge(do, tmp, by = c("FROM_AGE_INT","TO_AGE_INT"))

# re-order by ages
do <- do[order(FROM_AGE_INT,TO_AGE_INT),]
do <- data.table(do)
do[, PAIR_ID := seq(1,nrow(do)) ]

do <- do[order(PAIR_ID)]

tmp <- paste0(outfile.base, '-fitted_stan_model.rds')
cat("\n Read fitted dat ", tmp , "\n")
model_fit <- readRDS(file = tmp)

cat(" \n -------------------------------- \n Check convergence and divergences \n -------------------------------- \n")
if(grepl('Vanilla',args$stanModelFile)){
  fit.target.pars <- c('logit_y_mix_0','y_mix','log_alpha1_pair[1]','log_phi_pair[1]',"lp__")
}else if(grepl('Reg',args$stanModelFile)){
  fit.target.pars <- c('y_mix[1]','logit_y_mix_0','beta_age_src[1]','beta_age_src[2]','beta_age_src[3]','beta_age_src[4]','beta_age_src[5]',
                       'beta_age_rcp[1]','beta_age_rcp[2]','beta_age_rcp[3]','beta_age_rcp[4]','beta_age_rcp[5]',"lp__")
}else if(grepl('GP',args$stanModelFile)){
  fit.target.pars <- c('logit_y_mix_0','gpscale','lscale[1]','lscale[2]',"lp__")
}
if(grepl('1DGP',args$stanModelFile)){
  if(grepl('_agesrcrec_',args$job_tag)){
    fit.target.pars <- c('logit_y_mix_0',"lp__","gpscale","lscale")
  }else{
    fit.target.pars <- c('logit_y_mix_0',"lp__","gpscale","lscale")
  }
}

po <- model_fit$draws(inc_warmup = TRUE,
                      #format = 'draws_df',
                      variables = fit.target.pars
)
su <- as.data.table(posterior::summarise_draws(po))
#cat(paste("posterior mixing weight: ", round(inv.logit(su$median[1]), 3), "90% CI:", round(inv.logit(su$q5[1]), 3), ',', round(inv.logit(su$q95[1]), 3)))
tmp <- paste0(outfile.base,'-gp_scenario',i,"-convergence.csv")
write.csv(su, file = tmp, row.names = TRUE)
head(su)
su[,min(ess_bulk)]
su[,max(rhat)]


# posterior mixing probability 
po <- model_fit$draws(inc_warmup = FALSE,
                      format = 'draws_df',
                      variables = 'logit_y_mix_0'
)

po <- as.data.table(po)
setnames(po, colnames(po), gsub('^\\.','',colnames(po)))
po <- melt(po, id.vars = c('chain','iteration','draw'))
po[, omega:= inv.logit(value)]
po <- po[,
         list( q = quantile(omega, probs = c(0.5, 0.025, 0.975) ),
               stat = c('M','CL','CU')
         )
]
po <- dcast.data.table(po, .~stat, value.var = 'q')
po[, L:= paste0(round(M*100,0),'% [',round(CL*100,0),'-',round(CU*100,0),'%]')]

po

# args$file_stanModelGQs <- file.path(args$source_dir, 'stan_model_files',paste0(args$stanModelFile,'_GQs','.stan'))
# gq_program <- args$file_stanModelGQs
# mod_gq <- cmdstan_model(gq_program)
# fit_gq <- mod_gq$generate_quantities(model_fit, data = stan_data, seed = 123)
# 
# log_lik <- fit_gq$draws("log_lik", format = "matrix")
# LOO <- loo::loo(log_lik)
# print(LOO)
# saveRDS(LOO, file = paste0(outfile.base, "-LOO.rds"))
# saveRDS(LOO, file = paste0(outfile.base, "-LOO.csv"))

# trace plots ----
color_scheme_set("mix-blue-red")
po <- model_fit$draws(inc_warmup = TRUE,
                      #format = 'draws_df',
                      variables = fit.target.pars
)
p <- bayesplot:::mcmc_trace(po,
                            pars = fit.target.pars,
                            n_warmup = 5e2,
                            facet_args = list(ncol = 1, labeller = label_parsed))+
  bayesplot:::legend_move("bottom")+ggplot2::labs(x="Iteration")+ggplot2::theme(text = element_text(size = 16))+
  xaxis_text(size = 16)+facet_text(size = 16)+legend_text(size = 16)+yaxis_text(size = 16)+xaxis_ticks(size = 14)+
  yaxis_ticks(size = 14)
p
ggsave(file = paste0(outfile.base,'-trace_allpars.pdf'), p, w = 10, h = length(tmp)*10)

tmp <- su$variable[which.min(su$ess_bulk)]
po <- model_fit$draws(inc_warmup = TRUE,
                      #format = 'draws_df',
                      variables = tmp
)
p <- bayesplot:::mcmc_trace(po,
                            pars = tmp,
                            n_warmup = 5e2,
                            facet_args = list(ncol = 1, labeller = label_parsed))+
  bayesplot:::legend_move("bottom")+ggplot2::labs(x="Iteration")+ggplot2::theme(text = element_text(size = 16))+
  xaxis_text(size = 16)+facet_text(size = 16)+legend_text(size = 16)+yaxis_text(size = 16)+xaxis_ticks(size = 14)+
  yaxis_ticks(size = 14)+coord_cartesian(ylim=c(0,3))

p
ggsave(file = paste0(outfile.base,'-trace_lwstneff.pdf'), p, w = 10, h = 5)

# Pairs plots ----
#
cat("\n ----------- make pairs plots: start ----------- \n")
pd <- model_fit$draws(inc_warmup = FALSE,
                      variables = c(fit.target.pars))
bayesplot:::color_scheme_set("mix-blue-pink")
p <- bayesplot:::mcmc_pairs(pd,
                            pars = c(fit.target.pars),
                            diag_fun = "dens",
                            off_diag_fun = "hex"
)
ggsave(p, file = paste0(outfile.base, "-HMC-pairs_transmission_pars.pdf"), w=length(fit.target.pars)*2, h=length(fit.target.pars)*2)
cat("\n ----------- make pairs plots: end ----------- \n")


## GP diagnostics----
fit.target.pars = c("lscale[1]","lscale[2]","gpscale")
if(grepl('1DGP',args$stanModelFile)) fit.target.pars <- c('logit_y_mix_0',"lp__","gpscale","lscale")
po <- model_fit$draws(inc_warmup = FALSE,
                      format = 'draws_df',
                      variables = fit.target.pars
)
po <- as.data.table(po)
setnames(po, colnames(po), gsub('^\\.','',colnames(po)))
po <- melt(po, id.vars = c('chain','iteration','draw'))
po <- po[,
         list( q = quantile(value, probs = c(0.5) ),
               stat = c('M')
         ),
         by = c('variable')
]
po <- dcast.data.table(po, variable~stat, value.var = 'q')

diagnostic <- function(l,l_hat) l_hat + 0.01 >= l
m_QE <- function(c,l,S) ceiling(1.75 * c / (l/S))
l_QE <- function(c,m,S) round(S * 1.75 * c / m, 3)

l_hat1 <- po[1,2]
l_hat2 <- po[2,2]
c(l_hat1, l_hat2)
check1 <- diagnostic(2, l_hat1)
check2 <- diagnostic(2, l_hat2)

print(paste0("Check for lengthscale 1: ", check1))
print(paste0("Check for lengthscale 2: ", check2))



## plot mixing weights----
po <- model_fit$draws(inc_warmup = FALSE,
                      format = 'draws_df',
                      variables = 'y_mix'
)
po <- as.data.table(po)
setnames(po, colnames(po), gsub('^\\.','',colnames(po)))
po <- melt(po, id.vars = c('chain','iteration','draw'))
po[, PAIR_ID := as.integer(gsub('y_mix\\[([0-9]+)\\]','\\1',as.character(variable)))]
po <- merge(po, do, by = 'PAIR_ID')
po <- po[,
         list( q = quantile(value, probs = c(0.5, 0.025, 0.25, 0.75, 0.975) ),
               stat = c('M','CL','IL', 'IU', 'CU')
         ),
         by = c('FROM_AGE_INT')
]
po <- dcast.data.table(po, FROM_AGE_INT~stat, value.var = 'q')

if(grepl('Vanilla|Reg',args$stanModelFile)){
  p <- ggplot(po, aes(x = variable)) +
    geom_point(aes(y=M)) +
    geom_errorbar( aes(ymin = CL, ymax = CU), alpha = 0.5) +
    scale_y_continuous(labels = scales::percent, expand = c(0,0), limit = c(0,1)) +
    labs(x = 'time elapsed (years)', y = 'posterior mixing weight') +
    theme_bw()
  ggsave(file = paste0(outfile.base,'-rep_',replicate,'-postmixweight.png'), p, w = 7, h = 7)
}else{
  
  p1 <- ggplot(po, aes(x = FROM_AGE_INT)) +
    geom_ribbon( aes(ymin = CL, ymax = CU), alpha = 0.5) +
    geom_line( aes(y = M )) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(labels = scales::percent, expand = c(0,0), limit = c(0,1)) +
    labs(x = 'age of source', y = 'posterior mixing weight') +
    theme_bw()
  ggsave(file = paste0(outfile.base,'-rep_',replicate,'-postmixweight_FROMage.png'), p1, w = 7, h = 7)
  
  po <- model_fit$draws(inc_warmup = FALSE,
                        format = 'draws_df',
                        variables = 'y_mix'
  )
  po <- as.data.table(po)
  setnames(po, colnames(po), gsub('^\\.','',colnames(po)))
  po <- melt(po, id.vars = c('chain','iteration','draw'))
  po[, PAIR_ID := as.integer(gsub('y_mix\\[([0-9]+)\\]','\\1',as.character(variable)))]
  po <- merge(po, do, by = 'PAIR_ID')
  po <- po[,
           list( q = quantile(value, probs = c(0.5, 0.025, 0.25, 0.75, 0.975) ),
                 stat = c('M','CL','IL', 'IU', 'CU')
           ),
           by = c('TO_AGE_INT')
  ]
  po <- dcast.data.table(po, TO_AGE_INT~stat, value.var = 'q')
  p2 <- ggplot(po, aes(x = TO_AGE_INT)) +
    geom_ribbon( aes(ymin = CL, ymax = CU), alpha = 0.5) +
    geom_line( aes(y = M )) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(labels = scales::percent, expand = c(0,0), limit = c(0,1)) +
    labs(x = 'age of recipient', y = 'posterior mixing weight') +
    theme_bw()
  ggsave(file = paste0(outfile.base,'-rep_',replicate,'-postmixweight_TOAGE.png'), p2, w = 7, h = 7)
  
  p <- ggarrange(p1, p2, nrow = 2,align='hv')
  ggsave(file=paste0(outfile.base,'-gp-scenario_',i,'mixingweight_posterior_median_ages.png'), p, w=12, h=7)
  ggsave(file=paste0(outfile.base,'-gp-scenario_',i,'mixingweight_posterior_median_ages.pdf'), p, w=12, h=7)
  
}

## GP Random function----
fit.target.pars = c("f")
if(grepl('1DGP',args$stanModelFile)) fit.target.pars = c('f')
po <- model_fit$draws(inc_warmup = FALSE,
                      format = 'draws_df',
                      variables = fit.target.pars
)
po <- as.data.table(po)
setnames(po, colnames(po), gsub('^\\.','',colnames(po)))
po <- melt(po, id.vars = c('chain','iteration','draw'))

if(grepl('1DGP',args$stanModelFile)){
  po[, x_index := as.integer(gsub('f\\[([0-9]+)\\]','\\1',as.character(variable)))]
  if(grepl('_agesrcrec_',args$job_tag)) po[, x_index := as.integer(gsub('hsgp_f_src\\[([0-9]+)\\]','\\1',as.character(variable)))]
  ages1_grid = data.table(x = seq.int(1,stan_data$n,1))
  ages1_grid[, x_index := 1:nrow(ages1_grid)]
  # find all coordinates
  grid = as.data.table( expand.grid(x_index = ages1_grid$x_index) )
  grid = merge(grid, ages1_grid, by = 'x_index')
  po <- merge(po, grid, by = c('x_index'))
  setnames(po,c('x'),c('FROM_AGE_INT'))
}else{
  po[, x_index := as.integer(gsub('f\\[([0-9]+),([0-9]+)\\]','\\1',as.character(variable)))]
  po[, y_index := as.integer(gsub('f\\[([0-9]+),([0-9]+)\\]','\\2',as.character(variable)))]
  po <- merge(po, grid, by = c('x_index','y_index'))
  setnames(po,c('x','y'),c('FROM_AGE_INT','TO_AGE_INT'))
}

po <- po[,
         list( q = quantile((value), probs = c(0.5, 0.025, 0.25, 0.75, 0.975) ),
               stat = c('M','CL','IL', 'IU', 'CU')
         ),
         by = c("FROM_AGE_INT")
]
po <- dcast.data.table(po, FROM_AGE_INT~stat, value.var = 'q')
p1 <- ggplot(subset(po)) + geom_line( aes(x = FROM_AGE_INT, y = M )) +
  geom_ribbon( aes(x = FROM_AGE_INT,ymin = CL, ymax = CU), alpha = 0.3, fill = "light pink") +
  labs(x = 'Age of Source', y = 'Posterior Median of f',
       title = "Posterior Median of HSGP Random Function by Age of Source") +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position='bottom',
        text = element_text(size = 10),plot.title = element_text(face="bold",hjust = 0.5,size=11)
  )
p1
ggsave(file=paste0(outfile.base,'-gp-scenario_',i,'rand_function_posterior_median_FROMage.png'), p1, w=12, h=7)


po <- model_fit$draws(inc_warmup = FALSE,
                      format = 'draws_df',
                      variables = fit.target.pars
)
po <- as.data.table(po)
setnames(po, colnames(po), gsub('^\\.','',colnames(po)))
po <- melt(po, id.vars = c('chain','iteration','draw'))

po[, x_index := as.integer(gsub('f\\[([0-9]+),([0-9]+)\\]','\\1',as.character(variable)))]
po[, y_index := as.integer(gsub('f\\[([0-9]+),([0-9]+)\\]','\\2',as.character(variable)))]
po <- merge(po, grid, by = c('x_index','y_index'))
setnames(po,c('x','y'),c('FROM_AGE_INT','TO_AGE_INT'))

po <- po[,
         list( q = quantile((value), probs = c(0.5, 0.025, 0.25, 0.75, 0.975) ),
               stat = c('M','CL','IL', 'IU', 'CU')
         ),
         by = c("TO_AGE_INT")
]
po <- dcast.data.table(po, TO_AGE_INT~stat, value.var = 'q')

p2 <- ggplot(po) + geom_line( aes(x = TO_AGE_INT, y = M ))+
  geom_ribbon( aes(x = TO_AGE_INT, ymin = CL, ymax = CU), alpha = 0.3, fill = "light pink") +
  labs(x = 'Age of Recipient', y = 'Posterior Median of f',
       title = "Posterior Median of HSGP Random Function by Age of Recipient") +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position='bottom',
        text = element_text(size = 10),plot.title = element_text(face="bold",hjust = 0.5,size=11)
  )
p <- ggarrange(p1, p2, nrow = 2,align='hv')
p
ggsave(file=paste0(outfile.base,'-gp-scenario_',i,'rand_function_posterior_median_ages.png'), p, w=12, h=7)
ggsave(file=paste0(outfile.base,'-gp-scenario_',i,'rand_function_posterior_median_ages.pdf'), p, w=12, h=7)


## plot weighted posterior probabilities of being a pair ----
cat(" \n -------------------------------- \n Plot posterior probabilities of being a pair \n -------------------------------- \n")
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
#po[, PAIR_ID := as.integer(gsub('tpair_prob\\[([0-9]+)\\]','\\1',as.character(variable)))]
po <- merge(po, do, by = 'PAIR_ID')
cat(paste0('Number of pairs with tpair_prob > 95%: ',length(po$PAIR_ID[po$M>=0.95])))
cat(paste0('Number of pairs with tpair_prob > 75%: ',length(po$PAIR_ID[po$M>=0.75])))
cat(paste0('Number of pairs with tpair_prob > 50%: ',length(po$PAIR_ID[po$M>=0.5])))


p <- ggplot(po, aes(x = PAIR_ID, color = FROM_AGE_GP)) +
  geom_point(aes(y = M)) +
  geom_errorbar(aes(ymin = IL, ymax = IU)) +
  scale_y_continuous(labels = scales::percent, expand = c(0,0)) +
  ggsci::scale_color_npg() +
  labs(x = 'pair', y = 'transmission pair probability', colour = 'infector age group') +
  #facet_grid(.~TRANS_STAGE, scales = "free_x") +
  theme_bw() +
  theme( legend.position = 'bottom' )
ggsave(file = paste0(outfile.base,'-rep_',replicate,'-pi_ij_trsmgroup_w.png'), p, w = 12, h = 7)


p <- ggplot(po, aes(x = TIME_ELAPSED, color = FROM_AGE_GP)) +
  geom_point(aes(y = M)) +
  geom_errorbar(aes(ymin = IL, ymax = IU)) +
  scale_y_continuous(labels = scales::percent, expand = c(0,0)) +
  ggsci::scale_color_npg() +
  labs(x = 'time elapsed (years)', y = 'transmission pair probability', colour = 'infector age') +
  #facet_grid(.~TRANS_STAGE, scales = "free_x") +
  theme_bw() +
  theme( legend.position = 'bottom' )
ggsave(file = paste0(outfile.base,'-rep_',replicate,'-pi_ij_timeelapsed_trsmgroup_w.png'), p, w = 12, h = 7)


# plot for each unique recipient ID
p <- ggplot(po, aes(x = TO_SEQUENCE_ID, color = FROM_AGE_GP,y = M)) +
  geom_point() +
  geom_segment(aes(xend=TO_SEQUENCE_ID), yend=0,col="black",alpha=0.2) +
  scale_y_continuous(labels = scales::percent, expand = c(0,0)) +
  ggsci::scale_color_npg() +
  labs(x = 'pair', y = 'transmission pair probability', colour = 'infector age') +
  theme_bw() +
  theme( legend.position = 'bottom' )
ggsave(file = paste0(outfile.base,'-rep_',replicate,'-pi_ij_trsmgroup_recipientID_w.png'), p, w = 20, h = 7)


## read in clock quantiles
dps <- readRDS(file = paste0(outfile.base,'-clock_quantiles.rds'))
po[, d_TSeqT:= round(TIME_ELAPSED,1)]
po <- merge(po,dps,by=c('d_TSeqT'),all.x=T)

p1 <- ggplot(data=dps,aes(x=d_TSeqT)) +
   geom_ribbon(data = dps, aes(ymin = q2.5, ymax = q10), fill = p.palette[1], alpha = p.alpha) +
   geom_ribbon(data = dps, aes(ymin = q90, ymax = q97.5), fill = p.palette[1], alpha = p.alpha) +
   geom_ribbon(data = dps, aes(ymin = q10, ymax = q20), fill = p.palette[2], alpha = p.alpha) +
   geom_ribbon(data = dps, aes(ymin = q80, ymax = q90), fill = p.palette[2], alpha = p.alpha) +
   geom_ribbon(data = dps, aes(ymin = q20, ymax = q30), fill = p.palette[3], alpha = p.alpha) +
   geom_ribbon(data = dps, aes(ymin = q70, ymax = q80), fill = p.palette[3], alpha = p.alpha) +
   geom_ribbon(data = dps, aes(ymin = q30, ymax = q40), fill = p.palette[4], alpha = p.alpha) +
   geom_ribbon(data = dps, aes(ymin = q60, ymax = q70), fill = p.palette[4], alpha = p.alpha) +
   geom_ribbon(data = dps, aes(ymin = q40, ymax = q60), fill = p.palette[5], alpha = p.alpha) +
   geom_point(data=po, aes(x=TIME_ELAPSED,y=GEN_DIST, colour=FROM_AGE_GP, size=M), alpha=p.alpha) +
   theme_bw(base_size=14)+
   scale_colour_brewer(palette = 'Accent')+
   labs(x='Time elapsed \n(years)',y='Genetic distance \n(subst/site)', size='Posterior probability', color='Source age group')+
   #scale_y_continuous(expand = c(0,0), breaks=seq(0,0.18,0.02),labels=scales::label_percent(accuracy = 1L)) +
   scale_y_continuous(expand=c(0,0), breaks=seq(0,0.18,0.02), limits=c(0,0.18))+ 
   scale_x_continuous(expand = c(0,0), breaks=seq(0,16,2), limits=c(0, 16)) + 
   scale_size(range = c(1,3)) + 
   theme(legend.position='bottom') + 
   guides(color=guide_legend(nrow=2))
p1
ggsave(paste0(outfile.base, '-posterior_prob_sizes.png'), p1, w=10, h=8)
ggsave(paste0(outfile.base, '-posterior_prob_sizes.pdf'), p1, w=10, h=8)


points_col <- pal_brewer(palette='Accent')(1)
p2 <- ggplot(data=dps,aes(x=d_TSeqT)) +
  geom_ribbon(data = dps, aes(ymin = q2.5, ymax = q10), fill = p.palette[1], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q90, ymax = q97.5), fill = p.palette[1], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q10, ymax = q20), fill = p.palette[2], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q80, ymax = q90), fill = p.palette[2], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q20, ymax = q30), fill = p.palette[3], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q70, ymax = q80), fill = p.palette[3], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q30, ymax = q40), fill = p.palette[4], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q60, ymax = q70), fill = p.palette[4], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q40, ymax = q60), fill = p.palette[5], alpha = p.alpha) +
  geom_point(data=po, aes(x=TIME_ELAPSED,y=GEN_DIST, size=M), alpha=p.alpha, color=points_col) +
  theme_bw(base_size=12)+
  scale_colour_brewer(palette = 'Accent')+
  labs(x='Time elapsed \n(years)',y='Genetic distance of the pair \n', size='Posterior probability')+
  scale_y_continuous(expand = c(0,0), breaks=seq(0,0.18,0.02),labels=scales::label_percent(accuracy = 1L)) +
  scale_x_continuous(expand = c(0,0), breaks=seq(0,16,2), limits=c(0, 16)) + 
  scale_size(range = c(1,3))
p2
ggsave(paste0(outfile.base, '-posterior_prob_sizes-noage.png'), p2, w=10, h=8)
ggsave(paste0(outfile.base, '-posterior_prob_sizes-noage.pdf'), p2, w=10, h=8)

#GP PAF for each age band----
po <- model_fit$draws(inc_warmup = FALSE,
                       format = 'draws_df',
                      variables = c("pflows_from")
)
po <- as.data.table(po)
setnames(po, colnames(po), gsub('^\\.','',colnames(po)))
po <- melt(po, id.vars = c("chain", "iteration", "draw")) #Wide to Long
po <- po[,
         list( q = quantile(value, probs = c(0.5, 0.025, 0.25, 0.75, 0.975) ),
               stat = c('M','CL','IL', 'IU', 'CU')          ),
         by = c('variable')
]
po <- dcast.data.table(po, variable~stat, value.var = 'q')
po[, AGE_SOURCE_LABEL:=seq(1,5,1)]
po[, FROM_AGE_GP:=unique(do$FROM_AGE_GP)]

po[, L:= paste0(round(M*100,1),'% [',round(CL*100,1),'-',round(CU*100,1),'%]')]

saveRDS(po, file = paste0(outfile.base,"-gp_prob_age_coarse_scenario_",i,".RDS"))


#pal <- pal_npg("nrc")(4)[c(1,3,4)]
po[, TO_AGE_GP:='Overall']
g1 <- ggplot(subset(po)) + geom_bar(aes(x=TO_AGE_GP,y=M,fill=FROM_AGE_GP),stat='identity',position=position_dodge(width=0.9)) +
  geom_errorbar(aes(x=TO_AGE_GP,ymin=CL, ymax=CU,fill=FROM_AGE_GP),width=0.5, colour="black",position=position_dodge(width=0.9))	+
  scale_fill_brewer(palette='Accent') +
  labs(x='Age of recipient', y='Proportion of attributable\ninfections to age group', fill='Source age group') +
  theme_bw(base_size=14) +
  theme(legend.position = 'bottom')+
  #theme(legend.pos='none') +
  coord_cartesian(ylim = c(0,1)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1L),breaks=seq(0,1,0.1))
g1
ggsave(file = paste0(outfile.base,'-rep_',replicate,'-PAF_age_source.png'), g1, w = 10, h = 8)


# stratify by age gp of recipient
po <- model_fit$draws(inc_warmup = FALSE,
                      format = 'draws_df',
                      variables = 'tpair_prob_w'
)
po <- as.data.table(po)
setnames(po, colnames(po), gsub('^\\.','',colnames(po)))
po <- melt(po, id.vars = c('chain','iteration','draw'))
po[, PAIR_ID := as.integer(gsub(paste0('tpair_prob_w','\\[([0-9]+)\\]'),'\\1',as.character(variable)))]
#tmp <- subset(do, select = c('PAIR_ID','FROM_AGE_GP','TO_AGE_GP','TRANS_STAGE'))
tmp <- subset(do, select=c('PAIR_ID', 'FROM_AGE_GP', 'TO_AGE_GP'))
po <- merge(po, tmp, by = 'PAIR_ID')
po <- po[, list(value = sum(value)), by = c('draw','FROM_AGE_GP','TO_AGE_GP')]
tmp <- po[, list(total = sum(value)), by = c('draw','TO_AGE_GP')]
po <- merge(po, tmp, by = c('draw','TO_AGE_GP'))
po[, paf := value/total]

po <- po[,
         list( q = quantile(paf, probs = c(0.5, 0.025, 0.25, 0.75, 0.975) ),
               stat = c('M','CL','IL', 'IU', 'CU')
         ),
         by = c('TO_AGE_GP','FROM_AGE_GP')
]
po <- dcast.data.table(po, TO_AGE_GP+FROM_AGE_GP~stat, value.var = 'q')
po[, L:= paste0(round(M*100,1),'% [',round(CL*100,1),'-',round(CU*100,1),'%]')]

saveRDS(po,file=paste0(outfile.base,'-rep_',replicate,'-PAF_stratify_recipient_agegp_all_pairs','.RDS'))

g2 <- ggplot(subset(po)) + geom_bar(aes(x=TO_AGE_GP,y=M,fill=FROM_AGE_GP),stat='identity',position=position_dodge(width=0.9)) +
  geom_errorbar(aes(x=TO_AGE_GP,ymin=CL, ymax=CU,fill=FROM_AGE_GP),width=0.5, colour="black",position=position_dodge(width=0.9))	+
  scale_fill_brewer(palette = 'Accent') +
  labs(x='Age of recipient', y='',fill='Age of likely source') +
  theme_bw(base_size=14) + 
  #theme(legend.position='bottom') +
  theme(legend.position='none') + 
  coord_cartesian(ylim = c(0,1)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1L),breaks=seq(0,1,0.1))
g2
ggsave(file = paste0(outfile.base,'-rep_',replicate,'-PAF_age_recipient.png'), g2, w = 10, h = 8)

legend_t <- cowplot::get_legend(g2 + theme(legend.position = "bottom"))
#g <- ggarrange(g1 + rremove("xlab"),g2+ theme(legend.position='none'),ncol=2,widths=c(0.35,0.65),align='hv')
#g <- ggarrange(g, legend_t,ncol=1,heights=c(0.8,0.2))
g <- g1 + g2
ggsave(file = paste0(outfile.base,'-rep_',replicate,'-PAF_age_source_rec.png'), g, w = 16, h = 8)
ggsave(file = paste0(outfile.base,'-rep_',replicate,'-PAF_age_source_rec.pdf'), g, w = 16, h = 8)
