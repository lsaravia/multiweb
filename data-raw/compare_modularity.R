# Compara las formas de obtener la modularidad de una particion y que tiene en
# cuenta cada una (pesos / direccion). Todas las columnas de una fila se calculan
# sobre LA MISMA particion, asi que son comparables entre si.
#
# Uso: Rscript data-raw/compare_modularity.R   (desde la raiz del paquete)
# Requiere el binario `infomap` para las filas run_infomap.

suppressMessages({
  library(igraph)
  devtools::load_all(".", quiet = TRUE)
})
has_infomap <- nzchar(Sys.which("infomap"))

# --- Redes: netData[[1]] (dirigida, con 3 lazos) x {sin pesos, con pesos} ------
g_dir <- netData[[1]]
g_und <- as_undirected(g_dir, mode = "collapse")          # no hay aristas mutuas
set.seed(2)
add_w <- function(g) { E(g)$weight <- runif(ecount(g), 0.1, 10); g }
nets <- list(
  "dirigida   sin pesos" = g_dir,
  "dirigida   con pesos" = add_w(g_dir),
  "no dirigida sin pesos" = g_und,
  "no dirigida con pesos" = add_w(g_und)
)

# --- A) modularity(<comunidad>, ...) ignora sus argumentos --------------------
cat("\n=== A) modularity(<objeto comunidad>, ...) ignora directed= y weights= ===\n")
set.seed(1); cm <- cluster_spinglass(g_dir)
print(c(`modularity(cm)` = modularity(cm),
        `directed=TRUE`  = modularity(cm, directed = TRUE),
        `directed=FALSE` = modularity(cm, directed = FALSE),
        `weights=NULL`   = modularity(cm, weights = NULL)), digits = 5)
cat("(es solo max(cm$modularity): el numero que el algoritmo guardo al terminar)\n")

# --- B) Tabla: misma particion, distintas formulas ---------------------------
one_row <- function(red, algo, g, comm) {
  mem <- membership(comm)
  w <- if ("weight" %in% edge_attr_names(g)) E(g)$weight else NULL
  data.frame(
    red = red, algoritmo = algo, grupos = length(unique(mem)),
    guardado        = tryCatch(modularity(comm), error = function(e) NA_real_),
    mod_dir_sinpeso = modularity(g, mem, directed = TRUE),
    mod_und_sinpeso = modularity(g, mem, directed = FALSE),
    mod_dir_conpeso = if (is.null(w)) NA_real_ else modularity(g, mem, weights = w, directed = TRUE),
    mod_und_conpeso = if (is.null(w)) NA_real_ else modularity(g, mem, weights = w, directed = FALSE)
  )
}

rows <- list(); parts <- list()
for (nm in names(nets)) {
  g <- nets[[nm]]
  set.seed(1); sp <- cluster_spinglass(g)
  rows[[length(rows) + 1]] <- one_row(nm, "cluster_spinglass", g, sp)
  parts[[paste(nm, "spinglass")]] <- membership(sp)
  if (has_infomap) {
    d <- tempfile(); dir.create(d)
    ri <- run_infomap(g, output_dir = d, directed = is_directed(g))
    rows[[length(rows) + 1]] <- one_row(nm, "run_infomap", g, ri)
    parts[[paste(nm, "infomap")]] <- membership(ri)
  }
}
tab <- do.call(rbind, rows)

cat("\n=== B) Misma particion, distintas formulas (columnas) ===\n")
cat("guardado = modularity(<comunidad>) | mod_* = modularity(g, membership, ...)\n")
cat("dir/und = directed=TRUE/FALSE | sinpeso/conpeso = weights=NULL / weights=E(g)$weight\n\n")
print(tab, digits = 3, row.names = FALSE)

# --- C) calc_modularity() actual (hace su propia corrida de spinglass) --------
cat("\n=== C) calc_modularity() con cluster_spinglass: directed = FALSE (defecto) / TRUE ===\n")
cm_tab <- sapply(nets, function(g) c(
  `directed=FALSE` = calc_modularity(g, cluster_function = cluster_spinglass)$Modularity,
  `directed=TRUE`  = calc_modularity(g, cluster_function = cluster_spinglass, directed = TRUE)$Modularity))
print(round(cm_tab, 4))
cat("(usa modularity(g, membership, weights = pesos si existen, directed = directed);\n")
cat(" en redes no dirigidas directed no tiene efecto. Cada llamada corre su propio spinglass.)\n")

# --- D) El agrupamiento mismo, usa los pesos y la direccion? ------------------
cat("\n=== D) El AGRUPAMIENTO usa pesos/direccion? (NMI=1: misma particion) ===\n")
nmi <- function(a, b) round(compare(a, b, method = "nmi"), 3)
set.seed(1); sp_u <- cluster_spinglass(g_dir)
set.seed(1); sp_w <- cluster_spinglass(nets[["dirigida   con pesos"]])
set.seed(1); sp_n <- cluster_spinglass(g_und)
cat("spinglass: dirigida vs no dirigida (sin pesos)  NMI =", nmi(membership(sp_u), membership(sp_n)), "\n")
cat("spinglass: sin pesos vs con pesos (dirigida)    NMI =", nmi(membership(sp_u), membership(sp_w)),
    "| grupos:", length(sp_u), "vs", length(sp_w), "\n")
if (has_infomap) {
  run <- function(g, directed) { d <- tempfile(); dir.create(d); membership(run_infomap(g, output_dir = d, directed = directed)) }
  i_dir <- run(g_dir, TRUE); i_und <- run(g_dir, FALSE)
  i_w <- run(nets[["dirigida   con pesos"]], TRUE)
  cat("run_infomap: directed=TRUE vs directed=FALSE    NMI =", nmi(i_dir, i_und),
      "| grupos:", length(unique(i_dir)), "vs", length(unique(i_und)), "\n")
  cat("run_infomap: sin pesos vs con pesos             NMI =", nmi(i_dir, i_w),
      "| grupos:", length(unique(i_dir)), "vs", length(unique(i_w)), "\n")
}
set.seed(1); ci <- cluster_infomap(g_dir)
cat("cluster_infomap (igraph) en netData[[1]]: grupos =", length(ci), "-> modularidad", round(modularity(ci), 4), "\n")
