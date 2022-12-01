#' Interpolate usage data using a step function
#' @param usage usage values
#' @param usage_ts usage time points
#' @export
yearly_usg_stepfunc <- function(usage, usage_ts) {
    k <- length(usage)
    stopifnot(!is.na(usage))
    stopifnot(!is.na(usage_ts))
    stopifnot(k == length(usage_ts))
    out <- function(t) {
        if(t==usage_ts[k]) {
            return(usage[k])
        }
        idx <- min(which(usage_ts > t))-1
        stopifnot(idx >= 1)
        stopifnot(idx <= length(usage))
        return(usage[idx])
    }

    return(out)
}

#' Load synthetic usage data
#' @export
simulated_usage <- function() {
    f <- system.file('extdata',    
                        'synthetic_usage.csv',
                        package='ResistPhy',
                        mustWork = T)
    out <- as.data.frame(read.csv(f, header=T))
    return(out)
}

# GGPlot theme 1
thm1 <- function(...) {
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
        plot.title = element_text(hjust = 0.5,size=rel(1.0)), ...)
    return(thm1)
}

# GGPlot theme 2
thm2 <- function(...) {
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
        plot.title = element_text(hjust = 0.5,size=rel(1.0)), ...)
    return(thm2)
}

# GGPlot theme 3
thm3 <- function(...) {
    thm3 <- theme(
        axis.text.x=element_text(size=rel(0.7), angle = 45, hjust=1),
        axis.text.y=element_text(size=rel(0.7)),
        plot.margin = margin(0, 0, 0, 0, "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(size=rel(0.2), colour = "grey90"),
        plot.title = element_text(hjust = 0.5,size=rel(1.0)), ...)
    return(thm3)
}