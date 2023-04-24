fig5a <- plot_rr_map(out, 1, n_breaks = 100, min_disp_prob=0.5) + ggtitle("Resistant Lineage 1")
fig5a <- fig5a + (fig5a$theme %+replace% theme(legend.position="bottom",
    legend.direction = "horizontal",
    legend.text = element_text(angle=45, vjust=0.5, hjust=0.5, size=rel(1.0)))) +
    guides(fill=guide_coloursteps(label.position = "bottom", title.position="top", label.hjust = 0.8, title.hjust=0.5,nrow=1))

fig5b <- plot_rr_map(out, 2, n_breaks = 100, min_disp_prob=0.5) + ggtitle("Resistant Lineage 2")
fig5b <- fig5b + (fig5b$theme %+replace% theme(axis.ticks.y=element_blank(),
    legend.position="none",
    axis.text.y=element_blank(),
    axis.title.y=element_blank()))

y_min <- max(min(layer_scales(fig5a)$y$range$range), min(layer_scales(fig5b)$y$range$range))
y_max <- min(max(layer_scales(fig5a)$y$range$range), max(layer_scales(fig5b)$y$range$range))

fig5a <- fig5a+ylim(y_min,y_max)
fig5b <- fig5b+ylim(y_min,y_max)

x_min <- max(min(layer_scales(fig5a)$x$range$range), min(layer_scales(fig5b)$x$range$range))
x_max <- min(max(layer_scales(fig5a)$x$range$range), max(layer_scales(fig5b)$x$range$range))

fig5a <- fig5a+xlim(x_min,x_max)
fig5b <- fig5b+xlim(x_min,x_max)

fig5 <- (fig5a | fig5b) + plot_annotation(title = TeX(sprintf("Posterior Probability of $r_s(t) - r_r(t) > c$ given usage level"))) 