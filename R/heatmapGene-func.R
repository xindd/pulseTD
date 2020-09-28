heatmapGene <- function(a, c, b){
  colnames(a) = c('h0', 'h1', 'h2', 't0', 't1', 'b')
  a = as.matrix(scale(a))

  colnames(c) = c('h0', 'h1', 'h2', 't0', 't1', 'b')
  c = as.matrix(scale(c))

  colnames(b) = c('h0', 'h1', 'h2', 't0', 't1', 'b')
  b = as.matrix(scale(b))

  ht1 = Heatmap(a, name = "alpha", cluster_columns = FALSE, show_row_names = FALSE)
  ht2 = Heatmap(c, name = "gamma", cluster_columns = FALSE, show_row_names = FALSE)
  ht3 = Heatmap(b, name = "beta", cluster_columns = FALSE)

  ht1 + ht2 + ht3
}



