library(Matrix)

f1 <- function(x, parm){ parm[3] * (abs(x-parm[1])<parm[2]) + parm[4]*(abs(x-parm[1])>=parm[2]) }

#parm: beta_0, beta_1 in the model logit(pi) = beta_0 + beta_1*x 
f2 <- function(x, parm){
    x <- scale(x)
    reg <- exp(parm[1]+parm[2]*x)
    return(1/(1+reg))
}

#parm: number of groups, upper bound of pi_1
f3 <- function(x, parm){
    #divide x into different groups based on the values
    #ex: x=1:10 becomes: 1,1,1,1,1,2,2,2,2,2
    group = cut(x, parm[1], labels =FALSE)
    prob = seq(0.01, parm[2],length.out = 11)
    gp_idx = sample(10, parm[1])
    res = sapply(group, 
                 function(gp){return(runif(1, min=prob[gp_idx[gp]], max=prob[gp_idx[gp]+1]))})
    pi_1_range = list()
    for (i in 1:parm[1]) {
        pi_1_range[[i]] = c(prob[gp_idx[i]], prob[gp_idx[i]+1])
    }
    return(list(group = group, pi_1 = res, pi_1_range = pi_1_range))
}


#' Data generating function
#'
#' @params N A number of total observations.
#' @params nu A vector of length 2.
#' @params f A integer specify which type of function to use.
#' @return A list of x, nu and H
#' @export
#' @examples
#' generate_data(1e4, c(0, 1), c(0,50,0.9,0.1))
generate_data <- function(N,
                          nu,
                          f_param,
                          nfun
){
    x_i <- seq(-50, 50, length.out = N)
    
    # To Do: add more types of functions
    # f_pi_1 <- partial(
    #     switch(f, f1, f2), 
    #     parm = f_param)
    
    f_pi_1 <- switch(nfun, f1, f2, f3)
    
    if(nfun<3){
        pi_1 <- f_pi_1(x_i, parm = f_param)    
        H_i <- rbinom(N, 1, pi_1)
        nu_i <- rep(nu[1], N)
        nu_i[H_i==1] <- nu[2]
        
        return(list(x=x_i, nu=nu_i, H=H_i))
    }else{
        res <- f_pi_1(x_i, parm = f_param)
        pi_1 <- res$pi_1
        H_i <- rbinom(N, 1, pi_1)
        nu_i <- rep(nu[1], N)
        nu_i[H_i==1] <- nu[2]
        return(list(x=order(x_i), nu=nu_i, H=H_i, group=res$group, pi_1_range = res$pi_1_range))
    }
}

#' Covariance matrix generating function
#'
#' @params N A number of non nulls.
#' @params rho A number of correlation coefficient.
#' @params opt A string specify which type of dependence to use.
#' @params sparsity A float specify the level of sparsity.
#' @return Sigma A covariance matrix.
#' @export
#' @examples
#' cor_mat(15,0.5,'cs', sparsity=0.8)
cor_mat = function(N, rho=0, opt = c("cs","ar1"), sparsity=1.0){
    if (sparsity==1){
        Sigma = diag(1,N,N) 
    } else {
        n_dense = round((1-sparsity)*N,0)
        order = sample(N)
        
        if (match.arg(opt)=='cs' ) {
            Sigma = matrix(rho, n_dense, n_dense);
            diag(Sigma) = 1
        }
        if (match.arg(opt)=='ar1') {
            exponent = abs(matrix(1:n_dense-1, n_dense, n_dense, byrow=T) - (1:n_dense-1)); 
            Sigma = rho^exponent
        }
        Sigma = bdiag(Diagonal(N-n_dense), Sigma)
        Sigma = Sigma[order,order]
    }
    Sigma = Matrix(Sigma, sparse = TRUE) 
    return(Sigma)
}



