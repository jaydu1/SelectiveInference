source("utils.R")
source("algo.R")

# generate data
data <- generate_data(1e4, c(0, 1), c(0,50,0.9,0.1))

# create dependence
Cov <- cor_mat(sum(data$H==1), 0.9, 'ar1')

# run algorithms

# plot