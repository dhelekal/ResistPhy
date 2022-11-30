library(ggplot2)
library(latex2exp)
library(viridis)
library(hexbin)
library(patchwork)

thm1 <- function() {
    thm1 <- theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(), 
        axis.text.x=element_text(size=rel(0.7), angle = 45, hjust=1),
        aspect.ratio=1,
        plot.margin = margin(0, 0, 0, 0, "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(size=rel(0.2), colour = "grey80"),
        plot.title = element_text(hjust = 0.5,size=rel(1.0)))
    return(thm1)
}

#' GGPlot theme 2
thm2 <- function() {
    thm2 <- theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_text(size=rel(0.7), angle = 45, hjust=1),
        axis.text.y=element_text(size=rel(0.7)),
        aspect.ratio=1,
        plot.margin = margin(0, 0, 0, 0, "cm"),
        legend.position="none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(size=rel(0.2), colour = "grey90"),
        plot.title = element_text(hjust = 0.5,size=rel(1.0)))
    return(thm2)
}

blank <- plot_spacer()

gamma_sus_hist <- ggplot(out$draws_df, aes(x=gamma_sus)) + 
                    geom_histogram(aes(y = ..density..), position="identity", bins=30, fill="lightsteelblue3") + 
                    geom_function(fun = function(x) dlnorm(x, meanlog = log(out$gamma_guess), sdlog = out$gamma_log_sd), linetype="longdash", color="steelblue4") + 
                    labs(x="", y="", title=TeX("$\\gamma$")) +
                    geom_vline(xintercept=gamma_sus/365.0,color="red")+
                    theme_minimal() + thm1()
q_u_hist <- ggplot(out$draws_df, aes_(x=as.name("q_u[1]"))) + 
                    geom_histogram(aes(y = ..density..), position="identity", bins=30, fill="lightsteelblue3") + 
                    geom_function(fun = function(x) dlnorm(x, meanlog = 0, sdlog = 0.5), linetype="longdash", color="steelblue4") + 
                    labs(x="", y="", title=TeX("$q_U$")) +
                    geom_vline(xintercept=q_u,color="red")+
                    theme_minimal() + thm1()
q_t_hist <- ggplot(out$draws_df, aes_(x=as.name("q_t[1]"))) + 
                    geom_histogram(aes(y = ..density..), position="identity", bins=30, fill="lightsteelblue3") + 
                    geom_function(fun = function(x) dlnorm(x, meanlog = 0, sdlog = 0.5), linetype="longdash", color="steelblue4") + 
                    labs(x="", y="", title=TeX("$q_T$")) +
                    geom_vline(xintercept=q_t,color="red")+
                    theme_minimal() + thm1()
p1 <- ggplot(out$draws_df, aes_(x=as.name("gamma_sus"), y=as.name("q_u[1]"))) + 
        geom_hex(bins = 20) + 
        scale_fill_viridis() + 
        geom_vline(xintercept=gamma_sus/365.0,color="red")+
        geom_hline(yintercept=q_u,color="red")+
        theme_minimal() +
        thm2()
p2 <- ggplot(out$draws_df, aes_(x=as.name("gamma_sus"), y=as.name("q_t[1]"))) + 
        geom_hex(bins = 20) + 
        scale_fill_viridis() + 
        geom_vline(xintercept=gamma_sus/365.0,color="red")+
        geom_hline(yintercept=q_t, color="red")+
        theme_minimal() +
        thm2()
p3 <- ggplot(out$draws_df, aes_(x=as.name("q_u[1]"), y=as.name("q_t[1]"))) + 
        geom_hex(bins = 20) + 
        scale_fill_viridis() + 
        geom_vline(xintercept=q_u,color="red")+
        geom_hline(yintercept=q_t,color="red")+
        theme_minimal() +
        thm2() 

fig2 <- gamma_sus_hist + blank + blank +
    p1 + q_u_hist + blank +
    p2 + p3 + q_t_hist +
    plot_layout(ncol = 3)

