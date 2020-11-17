plot_s_curve <- function(obj, x, pvals,
                           alpha, alpha_BH, alpha_stoery,
                           xlab = "x", xlim = NULL,
                           disp_ymax = 0.2,
                           num_yticks = 3,
                           rand_seed_perturb = NA,
                           ...){
    
    title <- sprintf("Rejection threshold (alpha=%.02f)", alpha)
    n <- length(pvals)
    alphas <- obj$alphas
    if (alpha < min(alphas)){
        s <- rep(0, n)
    } else {
        ind <- max(which(alphas <= alpha))
        s <- obj$s[, ind]
    }
    
    x_ord <- order(x)
    pvals <- pvals[x_ord]
    s <- s[x_ord]
    x <- x[x_ord]
    
    if (!is.null(xlim)){
        inds <- which(x >= xlim[1] & x <= xlim[2])
        x <- x[inds]
        pvals <- pvals[inds]
        s <- s[inds]
        n <- length(inds)
    } else {
        xlim <- c(min(x), max(x))
    }
    
    par(...)
    
    what_type <- ifelse(pvals < s, 1, ifelse(pvals > 1 - s, 2, 3))
    plot(x, (1:n * disp_ymax) / n, type = "n", pch = ".",
         xaxs = "i", yaxs = "i", ylab = "p-values",
         xlab = xlab, xlim = xlim,
         col = c("red", "blue", "black")[what_type], yaxt = "n",
         main = title)
    axis(2, at = seq(0, disp_ymax, length.out = num_yticks),
         labels = seq(0, disp_ymax, length.out = num_yticks))
    hhi <- hlo <- .5
    vlo <- -.01
    vhi <- 1.01
    polygon(x = c(x, rev(x)), y = c(s, rep(vlo, n)),
            col = "#FFDDDD", border = "red")
    polygon(x = c(x, rev(x)), y = c(1 - s, rep(vhi, n)),
            col = "light blue", border = "blue")
    if(!is.na(rand_seed_perturb)) {
        set.seed(rand_seed_perturb)
        points(x, pvals + 0.001 * rnorm(n),
               pch = ".",
               col = c("red", "blue", "black")[what_type])
    }
    points(x, pvals, pch = ".",
           col = c("red", "blue", "black")[what_type])
    abline(a=NULL, b=NULL, h=alpha_BH, v=NULL, col="blue")
    abline(a=NULL, b=NULL, h=alpha_stoery, v=NULL, col="green")
    box()
}