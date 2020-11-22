library(ggplot2)
library(reshape2)
require(gridExtra)

plot_s_curve <- function(obj, x, pvals,
                           alpha, alpha_BH, alpha_stoery,
                           xlab = "x", xlim = NULL,
                           disp_ymax = 0.3,
                           num_yticks = 3,
                           rand_seed_perturb = NA){
    
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

plot_power <- function(alphas, df_BH, df_storey, df_adapt, ...){
    par(...)
    df_res_power <- data.frame(alpha=alphas, AdaPT=df_adapt$power,
                               Storey=df_storey$power, BH=df_BH$power)
    df_res_power <- melt(df_res_power, id.vars=c("alpha"), variable.name='method', value.name="value")
    p1 <- ggplot(data=df_res_power, aes(x=alpha, y=value, color=method)) + 
        geom_line() + 
        scale_y_continuous(

            # Features of the first axis
            name = "Power",
            limits = c(0,1),
            labels = scales::number_format(accuracy = 0.01)
        ) +
        theme_light() + 
        theme(legend.position = "none") + 
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank()) + 
        ggtitle("")
    
    df_res_FDP <- data.frame(alpha=alphas, AdaPT=df_adapt$FDP, 
                             Storey=df_storey$FDP, BH=df_BH$FDP)
    df_res_FDP <- melt(df_res_FDP, id.vars=c("alpha"), variable.name='method', value.name="value")
    p2 <- ggplot(data=df_res_FDP, aes(x=alpha, y=value, color=method)) + 
        geom_line() + 
        scale_y_continuous(
            
            # Features of the first axis
            name = "FDP",
            limits = c(0,0.4),
            labels = scales::number_format(accuracy = 0.01)
        ) +
        theme_light() + 
        theme(legend.position="bottom") 
        ggtitle("")
    grid.arrange(p1, p2,heights=c(0.6, 0.4), ncol=1)
}