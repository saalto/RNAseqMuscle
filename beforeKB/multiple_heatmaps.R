library(ComplexHeatmap)
library(circlize)

col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 1, nc = 2)))

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(Heatmap(MGM, col = col_fun, column_title = "Metabolic Gene Expression", show_heatmap_legend = FALSE), newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
draw(Heatmap(SGM, col = col_fun, column_title = "Other Gene Expression", show_heatmap_legend = FALSE), newpage = FALSE)
upViewport()

lgd = Legend(at = c(-2, -1, 0, 1, 2), col_fun = col_fun, title = "Log Fold Change", title_position = "topleft")

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
grid.draw(lgd)
upViewport()

upViewport()