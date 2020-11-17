library("adaptMT")
library("splines")
library(MASS)
source("utils.R")
source("algo.R")
source("plot.R")

N <- 1e3

# generate data
data <- generate_data(N, c(0, 2), c(0, 10, 0.9, 0.1))

# create dependence
Cov <- cor_mat(sum(data$H==1), 0.1, 'ar1')

# generate p value
data$z <- rnorm(N, data$nu)
data$z[data$H==1] <- mvrnorm(1, data$nu[data$H==1], Cov)
data$pvals <- 1 - pnorm(data$z)

# run algorithms
df_BH <- summary_BH(data$pvals, data$H, alphas = seq(0.01, 0.3, 0.01))
df_storey <- summary_storey(data$pvals, data$H, alphas = seq(0.01, 0.3, 0.01))

formulas <- paste0("ns(x, df = ", 6:7, ")")
adapt <- adapt_glm(x = data.frame(x = data$x), pvals = data$pvals, pi_formulas = formulas,
                 mu_formulas = formulas,  nfits = 10, alphas = seq(0.01, 0.3, 0.01))
df_adapt <- summary_adapt(adapt, data$pvals, data$H)

# plot
alpha <- 0.1
plot_s_curve(adapt, data$x, data$pvals, alpha)