# check component-wise EM function 
require(mvtnorm)
require(EMCluster)
require(ggplot2)
require(ggpubr)
require(patchwork)
require(RColorBrewer)

set.seed(1726)

# source function
args <- list(source_dir = '/Users/bethanyastor/Library/CloudStorage/OneDrive-ImperialCollegeLondon/Research Project/source.attribution.bmm', 
             outdir = '/Users/bethanyastor/Library/CloudStorage/OneDrive-ImperialCollegeLondon/Research Project/source.attribution.bmm/figures'
            )
source(file.path(args$source_dir, 'R', 'component_EM.R'))

# create test data set (well-separated clusters)
test_data <- rbind(rbind(rmvnorm(100, c(-2,2), matrix(c(.1,0.01,0.01,.1), nrow=2)), 
                         rmvnorm(100, c(0,2), matrix(c(.1,0.01,0.01,.1), nrow=2)), 
                         rmvnorm(100, c(2,2), matrix(c(.1,0.01,0.01,.1), nrow=2))))
test_df <- data.frame(x=test_data[,1], y=test_data[,2], c=c(rep(1,100), rep(2,100), rep(3,100)))

# plot test data
test_plot <- ggplot(test_df) +
  geom_point(aes(x=x, y=y)) +
  theme_bw()
test_plot

test_fit <- fit_cEM(test_data, kmin=1, kmax=15, max_iter=2500)
k_best <- test_fit$k
cat(paste('Optimal number of clusters: k =', k_best))

# add gmm density 
x <- seq(-5, 5, 0.1)
y <- seq(-2, 6, 0.1)
grid <- expand.grid(x=x, y=y)

grid$density <- gmm_density(as.matrix(grid), pi=test_fit$pi, mu=test_fit$Mu, sigma=test_fit$C, K=test_fit$k)


line_col <- pal_brewer(palette='Accent')(8)[6]
plot_w_dens_1 <- ggplot(test_df) +
  geom_point(aes(x=x, y=y), size=1.2, alpha=0.7) +
  # gmm density 
  geom_contour(data=grid, aes(x=x, y=y, z=density),bins=5, color=line_col) +
  theme_bw(base_size=16) +
  labs(x='Variable 1', y='Varible 2')
plot_w_dens_1
ggsave(paste0(args$outdir, '/CEM_well_separated.png'), plot_w_dens_1, h=5, w=7)
  


# overlapping clusters
test_data <- rbind(rbind(rmvnorm(100, c(-2,2), matrix(c(1,0.01,0.01,1), nrow=2)), 
                         rmvnorm(100, c(0,2), matrix(c(1,0.01,0.01,1), nrow=2)), 
                         rmvnorm(100, c(2,2), matrix(c(1,0.01,0.01,1), nrow=2))))
test_df <- data.frame(x=test_data[,1], y=test_data[,2], c=c(rep(1,100), rep(2,100), rep(3,100)))
# plot test data
test_plot <- ggplot(test_df) +
  geom_point(aes(x=x, y=y)) +
  theme_bw()
test_plot

test_fit <- fit_cEM(test_data, kmin=1, kmax=15, max_iter=5000)
k_best <- test_fit$k
cat(paste('Optimal number of clusters: k =', k_best))

# add gmm density 
x <- seq(-5, 5, 0.1)
y <- seq(-2, 6, 0.1)
grid <- expand.grid(x=x, y=y)

grid$density <- gmm_density(as.matrix(grid), pi=test_fit$pi, mu=test_fit$Mu, sigma=test_fit$C, K=test_fit$k)

plot_w_dens_2 <- ggplot(test_df) +
  geom_point(aes(x=x, y=y), size=1.2, alpha=0.7) +
  # gmm density 
  geom_contour(data=grid, aes(x=x, y=y, z=density), bins=5, color=line_col) +
  theme_bw(base_size = 16) +
  labs(x='Variable 1', y='Variable 2')
plot_w_dens_2
ggsave(paste0(args$outdir, '/CEM_overlapping.png'), plot_w_dens_2, h=5, w=7)

g <- plot_w_dens_1 + plot_w_dens_2
ggsave(paste0(args$outdir, '/CEM_fit.png'), g, h=5, w=12)
ggsave(paste0(args$outdir, '/CEM-fit.pdf'), g, h=5, w=12)

# another one 
test_data <- rbind(rmvnorm(100, c(-2,-2), matrix(c(1,0.01,0.01,1), nrow=2)), 
                    rmvnorm(200, c(2,2), matrix(c(1,0.01,0.01,1), nrow=2)))

test_df <- data.frame(x=test_data[,1], y=test_data[,2], c=c(rep(1,100), rep(2,200)))

# plot test data
test_plot <- ggplot(test_df) +
  geom_point(aes(x=x, y=y)) +
  theme_bw()
test_plot

test_fit <- fit_cEM(test_data, kmin=1, kmax=15, max_iter=2500)
k_best <- test_fit$k
cat(paste('Optimal number of clusters: k =', k_best))

# add gmm density 
x <- seq(-5, 5, 0.1)
y <- seq(-5, 6, 0.1)
grid <- expand.grid(x=x, y=y)

grid$density <- gmm_density(as.matrix(grid), pi=test_fit$pi, mu=test_fit$Mu, sigma=test_fit$C, K=test_fit$k)

plot_w_dens <- ggplot(test_df) +
  geom_point(aes(x=x, y=y)) +
  # gmm density 
  geom_contour(data=grid, aes(x=x, y=y, z=density)) +
  theme_bw() 
plot_w_dens
