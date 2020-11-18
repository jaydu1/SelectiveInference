library(ggplot2)
library(reshape2)

plot_s_curve <- function(obj, x, pvals,
                           alpha, alpha_BH, alpha_stoery,
                           xlab = "x", xlim = NULL,
                           disp_ymax = 0.3,
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

plot_power <- function(alphas, df_BH, df_storey, df_adapt, ...){
    par(...)
    df_res_power <- data.frame(alpha=alphas, AdaPT=df_adapt$power,
                               Storey=df_storey$power, BH=df_BH$power)
    df_res_power <- melt(df_res_power, id.vars=c("alpha"), variable.name='method', value.name="value")
    df_res_power$group <- 'Power'
    
    df_res_FDP <- data.frame(alpha=alphas, AdaPT=df_adapt$FDP, 
                               Storey=df_storey$FDP, BH=df_BH$FDP)
    df_res_FDP <- melt(df_res_FDP, id.vars=c("alpha"), variable.name='method', value.name="value")
    df_res_FDP$group <- 'FDP'
    
    df_res <- rbind(df_res_power, df_res_FDP)
    # df_res$method <- factor(df_res$method, levels=c("AdaPT", "Stoery", "BH"))
    ggplot(data=df_res, aes(x=alpha, y=value, group=interaction(method, group), color=method, linetype=group)) + 
        geom_line() + 
        scale_y_continuous(

            # Features of the first axis
            name = "Power",

            # Add a second axis and specify its features
            sec.axis = sec_axis(~ 1*., name="FDP"),
            limits = c(0,1)
        ) +
        theme(legend.position="bottom") +
        
        ggtitle("")

}