# Test Plot multiplex adapted from muxviz
#

library(multiweb)
library(muxViz)
library(igraph)
library(dplyr)
library(RColorBrewer)
library(viridis)
# === Step 1: Load data ===
fileName <- system.file("extdata", package = "multiweb")
dn <- list.files(fileName, pattern = "^Kefi2015.*\\.txt$")
g_list <- readNetwork(dn, fileName, skipColumn = 2)
names(g_list) <- c("Negative", "Positive", "Trophic")
lapply(g_list, function(x){igraph::gorder(x)})
# === Step 2: Detect communities using multiweb ===
res <- run_infomap_multi(g_list, layer_names = names(g_list))
communities <- res$communities

# === Paso 3: Crear paletas viridis por módulo ===
modules_multi <- sort(unique(communities$module))
pal.mux <- setNames(sample(viridis(length(modules_multi))), modules_multi)

# Paleta para nodos agregados (módulo dominante por especie)
aggr_df <- communities %>%
  group_by(node, module) %>%
  summarise(flow = sum(flow), .groups = "drop") %>%
  group_by(node) %>%
  slice_max(order_by = flow, n = 1)

modules_aggr <- sort(unique(aggr_df$module))
pal.aggr <- setNames(sample(viridis(length(modules_aggr))), modules_aggr)

# === Paso 4: Matriz de colores por nodo-capa ===
Layers <- length(g_list)
Nodes <- length(unique(communities$node))
node_names <- sort(unique(communities$node))
node_index <- setNames(seq_along(node_names), node_names)

node.colors.matrix <- matrix("#dadada", nrow = Nodes, ncol = Layers)

for (l in seq_len(Layers)) {
  dftmp <- communities %>% filter(layer == names(g_list)[l])
  idxs <- node_index[dftmp$node]
  node.colors.matrix[idxs, l] <- pal.mux[as.character(dftmp$module)]
}

# === Paso 5: Colores agregados por nodo ===
node.colors.aggr <- rep("#dadada", Nodes)
node.colors.aggr[node_index[aggr_df$node]] <- pal.aggr[as.character(aggr_df$module)]

# === Paso 6a: Crear matriz de flujo por nodo-capa ===
node.flow.matrix <- matrix(0, nrow = Nodes, ncol = Layers)

for (l in seq_len(Layers)) {
  dftmp <- communities %>% filter(layer == names(g_list)[l])
  idxs <- node_index[dftmp$node]
  node.flow.matrix[idxs, l] <- dftmp$flow
}

# Optional: Rescale flows for visual clarity (e.g., to radius ~ sqrt(flow))
# Avoid log(0) if flow contains zeros
eps <- 1e-5
scaled.flow.matrix <- sqrt(node.flow.matrix ) *100  # adjust multiplier as needed

# === Paso 6: Layout y visualización ===
lay <- layoutMultiplex(g_list, layout = "drl", ggplot.format = FALSE, box = TRUE)

plot_multi3D(
  g_list,
  layer.layout = lay,
  layer.colors = RColorBrewer::brewer.pal(Layers, "Dark2"),
  layer.shift.x = 0.5,
  layer.space = 2,
  layer.labels = names(g_list),
  layer.labels.cex = 1.5,
  layer.alpha = rep(0.001, length(g_list)),
  node.size.values = scaled.flow.matrix,  # <- flow-based sizing
  node.size.scale = 1,                    # <- no further scaling here  node.size.scale = 0.5,
  node.colors = node.colors.matrix,
  edge.colors = "#FFFFFF",
  edge.alpha = 0.01,
  node.colors.aggr = node.colors.aggr,
  aggr.alpha = 0.01,
  show.aggregate = TRUE,
  show.nodeLabels= TRUE
)

# === Paso 7: Exportar imagen si lo deseas ===
rgl::snapshot3d("multiweb_infomap_3D.png", fmt = "png", width = 2024, height = 2024)


# Adaptar para muxViz
commResult <- convert_infomap_result_to_muxviz(communities)

# Paleta de colores por módulo
library(viridis)
pal.mux <- sample(viridis(commResult$modules.multi))
pal.aggr <- sample(viridis(commResult$modules.aggr))

# Tabla tipo ggplot con modulos por capa y agregados
gplt <- plot_multimodules(commResult, module.colors = pal.mux, show.aggregate = TRUE)


#Calculate PageRank and Degree versatility
sam <- convert_to_supra_adjacency(g_list,use_names = TRUE)$supra_matrix

pr <- GetMultiPageRankCentrality(sam, Layers,Nodes)
deg <- GetMultiDegree(M, Layers,Nodes, isDirected=F)

