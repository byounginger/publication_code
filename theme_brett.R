# Code for theme_brett

theme_brett <- function () {
  theme(axis.text = element_text(color = 'black', size = 6),
        axis.line = element_line(color = 'black'),
        axis.title = element_text(size = 8),
        legend.key.height = unit(0.1, 'in'),
        panel.background = element_blank(),
        strip.text = element_text(size = rel(0.5)),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 5),
        legend.margin = margin(l = 0, r = 0, unit = 'pt'),
        legend.key.width = unit(0.1, unit = 'in'))
}
