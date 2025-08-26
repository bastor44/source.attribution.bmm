require(data.table)
require(ggplot2)
require(ggpubr)
require(cmdstanr)
require(posterior)
require(bayesplot)
require(RColorBrewer)
require(scales)

## Preamble----
if (1)
{
  pkg.dir <- '~/OneDrive - Imperial College London/Research Project'
  indir <- file.path(pkg.dir,'data_other')
  outdir <- '~/OneDrive - Imperial College London/Research Project/source.attribution.bmm/molecular_clock'
  stan.model <- 'clock_model_gamma_hier_220315'
}

outfile.base <- file.path(outdir, "clock_model_gamma_hier_220315")

## load Belgium data----
infile <- '140921_set7_INFO_TRM.R'
file <- paste(indir, '/',infile, sep = '')
load(file)

# exclude B->A, same as A->B
trm.pol.nA <- subset(trm.pol, withA == FALSE, select = c(d_SeqT, d_TSeqT, BRL, FROM, TO))
setkey(trm.pol.nA, d_TSeqT)
trm.pol.nA[, id := seq_len(nrow(trm.pol.nA))]
trm.pol.nA[, pair_id := as.integer(factor(paste0(FROM,'->',TO)))]


# plot raw data
p <- ggplot(trm.pol, aes(x = d_TSeqT, y = BRL)) +
  geom_point(size = 1.2, alpha = 0.5, data = subset(trm.pol, withA == FALSE & BRL > 0.003)) +
  scale_x_continuous(breaks = seq(0,20,2), expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0,0.2,0.02), expand = c(0,0)) +
  coord_trans(ylim = c(0,0.1), xlim = c(0,17) ) +
  theme_bw() +
  labs(x = 'time elapsed\n(years)', y = 'genetic distance\n(subst/site)', title = 'epidemiologically confirmed\ntransmission pairs\n') +
  theme(panel.grid.major = element_line(colour = "grey70", size = 0.4), panel.grid.minor = element_line(colour = "grey70", linewidth = 0.2)) +
  geom_smooth(method = "lm", se = FALSE, color = "black", formula = y ~ x^2)
p
ggsave(file = file.path(outdir,'Genetic_distances_against_time_elapsed_confirmed_pairs.png'), p, w = 6, h = 6)

# plot raw data with individuals
p <- ggplot(trm.pol, aes(x = d_TSeqT, y = BRL)) +
  geom_point(aes(col = paste0(FROM,'->',TO)), size = 1.2, alpha = 0.5, data = subset(trm.pol, withA == FALSE & BRL > 0.003)) +
  scale_x_continuous(breaks = seq(0,20,2), expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0,0.2,0.02), expand = c(0,0)) +
  ggsci::scale_color_aaas() +
  coord_trans(ylim = c(0,0.1), xlim = c(0,17) ) +
  theme_bw() +
  labs( x = 'time elapsed\n(years)',
        y = 'genetic distance\n(subst/site)',
        colour = 'pair',
        title = 'epidemiologically confirmed\ntransmission pairs\n'
  ) +
  theme(panel.grid.major = element_line(colour = "grey70", size = 0.4), panel.grid.minor = element_line(colour = "grey70", size = 0.2)) +
  geom_smooth(method = "lm", se = FALSE, color = "black", formula = y ~ x^2)
p
ggsave(file = file.path(outdir,'Genetic_distances_against_time_elapsed_confirmed_pairs_with_col.png'), p, w = 6, h = 6)



dpr <- data.table(d_TSeqT = seq(0.1,18,0.1))
dpr[, id := seq_len(nrow(dpr))]

cat('---------- fit in stan -----------')

options(mc.cores = parallel::detectCores())
model <- cmdstanr::cmdstan_model(file.path(pkg.dir,'/source.attribution.bmm/stan_model_files',paste0(stan.model,'.stan')),
                                 force_recompile = TRUE
)

# Fit bayesian model in stan----
stan_data <- list(
  N = nrow(trm.pol.nA),
  x = trm.pol.nA$d_TSeqT,
  y = trm.pol.nA$BRL,
  K = length(unique(trm.pol.nA$pair_id)),
  pair_idx = trm.pol.nA$pair_id,
  N_pr = nrow(dpr),
  x_pr = dpr$d_TSeqT
)

stan_init <- list(
  log_alpha1 = log(10^(-2.5)),
  log_alpha1_pair = rep(0, stan_data$K),
  log_alpha1_pair_sd = abs((log(10^(-2.5)) - log(10^(-2.1)))/2),
  log_phi = log(5e-3),
  log_phi_pair = rep(0, stan_data$K)
)

model_fit <- model$sample(
  data = stan_data,
  seed = 42L,
  chains = 4,
  parallel_chains = min(4, parallel::detectCores()),
  refresh = 5e2,
  iter_warmup = 5e2,
  iter_sampling = 2e3,
  save_warmup = TRUE,
  init = list(stan_init, stan_init, stan_init, stan_init)
)

# Save output----
tmp <- paste0(outfile.base,"-stan_fit.rds")
cat("\n Save fitted data to file ", tmp , "\n")
model_fit$save_object(file = tmp)




# Clock model figure for paper 
## load estimated molecular clock ----
cat(" \n ------------- \n Load quantiles from fitted molecular clock model \n ------------ \n")

cm <- readRDS(file.path(outdir,'clock_model_gamma_hier_220315-stan_fit.rds'))
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
# get quantiles of mean/median distances at prediction points
dpr <- data.table(d_TSeqT = seq(0.1,18,0.1))
dpr[, id := seq_len(nrow(dpr))]
dpo <- data.table(iter = seq_along(po$log_alpha1),
                  log_alpha1 = po$log_alpha1,
                  log_alpha1_pair_sd = po$log_alpha1_pair_sd,
                  log_phi = po$log_phi,
                  log_phi_pair_sd = po$log_phi_pair_sd,
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

p.palette <- RColorBrewer::brewer.pal(5,'Blues')
p.alpha <- 0.7


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

dpr <- data.table(d_TSeqT = seq(0.01,40,0.01))
dpr[, id := seq_len(nrow(dpr))]
dpr[, DUMMY := 1]
dpo <- merge(dpo, dpr, by = 'DUMMY', allow.cartesian = TRUE)
set(dpr, NULL, 'DUMMY', NULL)
dpo[, y_pr := rgamma(nrow(dpo), shape = er * d_TSeqT * beta, rate = beta)]

dps <- dpo[,
           list(
             p = quantile(y_pr, prob = c(0.025, seq(0.1, .9, 0.1), 0.975, 0.25, 0.75)),
             qlabel = paste0('q',c(0.025, seq(0.1, .9, 0.1), 0.975, 0.25, 0.75)*100)
           ),
           by = c('d_TSeqT')
]
dps <- dcast.data.table(dps, d_TSeqT ~ qlabel, value.var = 'p')
saveRDS(dps,file = paste0(outfile.base,'-clock_quantiles.rds'))


point_pal <- pal_brewer(palette='Paired')(12)[3:12]
p <- ggplot(dps, aes(x = d_TSeqT)) +
  # clock model quantiles 
  geom_ribbon(data = dps, aes(ymin = q2.5, ymax = q10), fill = p.palette[1], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q90, ymax = q97.5), fill = p.palette[1], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q10, ymax = q20), fill = p.palette[2], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q80, ymax = q90), fill = p.palette[2], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q20, ymax = q30), fill = p.palette[3], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q70, ymax = q80), fill = p.palette[3], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q30, ymax = q40), fill = p.palette[4], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q60, ymax = q70), fill = p.palette[4], alpha = p.alpha) +
  geom_ribbon(data = dps, aes(ymin = q40, ymax = q60), fill = p.palette[5], alpha = p.alpha) +
  geom_line(data = dps, aes(y = q50)) +
  # belgium data pairs 
  geom_point(data = trm.pol.nA, aes(x = d_TSeqT, y = BRL, col = paste0(FROM, '->', TO)), size=1.2, alpha=0.5) + 
  scale_x_continuous(breaks = seq(0,20,2), expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0,0.2,0.02), expand = c(0,0)) +
  #scale_color_brewer(palette='Accent') +
  scale_color_manual(values=point_pal) + 
  coord_trans(ylim = c(0,0.18), xlim = c(0,17) ) +
  theme_bw(base_size = 14) + 
  theme(panel.grid.major = element_line(colour = alpha("grey70", p.alpha), size = 0.4), 
        panel.grid.minor = element_line(colour = alpha("grey70", 0.5), size = 0.2)) + 
  labs(x = 'Time elapsed (years)', y='Genetic distance \n(subst/site)',
       colour='Pair', title='Epidemiologically confirmed transmission pairs') +
  guides(colour = guide_legend(nrow = 1)) + 
  theme(legend.position = 'bottom', legend.title = element_text(size=12), guide_legend(nrow=1))
  
p
ggsave(file=paste0(outfile.base,'-clock-confirmed-pairs-with-colour.png'), p, w=7, h=7)
ggsave(file=paste0(outfile.base,'-clock-confirmed-pairs-with-colour.pdf'), p, w=7, h=6)

