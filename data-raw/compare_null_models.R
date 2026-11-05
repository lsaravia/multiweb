# Compara los modelos nulos que conservan el grado: curve_ball() y shuffle_network_deg(),
# usando la modularidad (Q) y el codelength de run_infomap() como estadisticos, y verifica
# que se conserva (no se elimina) el canibalismo (autoenlaces / diagonal).
#
# Uso (desde la raiz del paquete; requiere el binario infomap):
#   Rscript data-raw/compare_null_models.R "<nombre en netData>" [nrep] [ks separados por coma]
#   Rscript data-raw/compare_null_models.R "Potter Cove" 30 1,2,5,10,20

args  <- commandArgs(trailingOnly = TRUE)
netnm <- if (length(args) >= 1) args[1] else names(multiweb::netData)[1]
nrep  <- if (length(args) >= 2) as.integer(args[2]) else 100L
ks    <- if (length(args) >= 3) as.numeric(strsplit(args[3], ",")[[1]]) else c(1, 2, 5, 10, 20, 50, 100)
# ks: shuffle_network_deg con k * (numero de aristas) intercambios

suppressMessages({
  library(igraph)
  devtools::load_all(".", quiet = TRUE)
})
stopifnot(nzchar(Sys.which("infomap")))

g  <- netData[[netnm]]
A0 <- as.matrix(as_adjacency_matrix(g, sparse = FALSE))
E0 <- sum(A0)
set.seed(123)

# curve_ball() tiene `for (rep in 1:5*l_hp)`, que R evalua como (1:5)*l_hp (5 iteraciones).
# Copias locales con el bucle corregido a 1:(m*l_hp), sin modificar el paquete.
txt <- paste(deparse(curve_ball), collapse = "\n")
stopifnot(grepl("1:5 * l_hp", txt, fixed = TRUE))
curve_ball_m <- function(m) {
  f <- eval(parse(text = gsub("1:5 * l_hp", sprintf("1:(%d * l_hp)", m), txt, fixed = TRUE)))
  environment(f) <- environment(curve_ball)
  f
}

# --- Estadisticos de cada red aleatoria --------------------------------------
stats_net <- function(x) {
  A <- if (is.matrix(x)) (x > 0) * 1 else as.matrix(as_adjacency_matrix(x, sparse = FALSE))
  gx <- graph_from_adjacency_matrix(A, mode = "directed")       # diag = TRUE: conserva lazos
  d <- tempfile(); dir.create(d)
  r <- run_infomap(gx, output_dir = d, directed = TRUE)
  data.frame(
    Q_dir = modularity(gx, membership(r), directed = TRUE),     # modularidad dirigida
    Q_und = modularity(gx, membership(r), directed = FALSE),    # modularidad no dirigida
    Codelength = r$codelength, Modules = length(r),
    Loops = sum(diag(A)),
    Kept = sum(A & A0) / E0,                                    # fraccion de aristas originales
    DegOK = all(rowSums(A) == rowSums(A0)) && all(colSums(A) == colSums(A0)),
    Comps = components(gx, mode = "weak")$no
  )
}
run_model <- function(nets) do.call(rbind, lapply(nets, stats_net))

res <- list()
cat("Red:", netnm, "| nodos =", nrow(A0), "| aristas =", E0,
    "| autoenlaces (canibalismo) =", sum(diag(A0)), "\n")
res[["original"]] <- stats_net(g)

cat("curve_ball (actual)...\n")
res[["curve_ball actual (5 intentos)"]] <- run_model(curve_ball(g, nsim = nrep))
for (m in c(5, 20)) {
  cat("curve_ball corregido", m, "*n ...\n")
  res[[sprintf("curve_ball corregido (%d*n)", m)]] <- run_model(curve_ball_m(m)(g, nsim = nrep))
}
for (k in ks) {
  cat("shuffle_network_deg k =", k, "...\n")
  nets <- lapply(seq_len(nrep), function(i) shuffle_network_deg(g, delta = k * E0, weighted = FALSE))
  res[[sprintf("shuffle_network_deg k=%g", k)]] <- run_model(nets)
}

# --- Resumen -----------------------------------------------------------------
sm <- function(x) c(
  Q_dir = mean(x$Q_dir), sd_Qd = sd(x$Q_dir), Q_und = mean(x$Q_und),
  Codelen = mean(x$Codelength), sd_Cl = sd(x$Codelength),
  Loops = mean(x$Loops), Kept = mean(x$Kept), DegOK = mean(x$DegOK), Comp1 = mean(x$Comps == 1))
tab <- as.data.frame(t(sapply(res, sm)))
cat("\n=== ", netnm, ": media por modelo (", nrep, " redes) ===\n", sep = "")
cat("Q_dir / Q_und = modularidad dirigida / no dirigida de la particion de run_infomap(directed=TRUE)\n")
cat("Kept = fraccion de aristas originales conservadas | DegOK = grado entrada/salida igual | Comp1 = 1 componente\n\n")
print(round(tab, 3))

# --- Comparacion al converger -------------------------------------------------
kmax <- sprintf("shuffle_network_deg k=%g", max(ks)); cbm <- "curve_ball corregido (20*n)"
cat("\n=== ", kmax, " vs ", cbm, ": KS test ===\n", sep = "")
for (v in c("Q_dir", "Q_und", "Codelength", "Loops", "Kept")) {
  a <- res[[kmax]][[v]]; b <- res[[cbm]][[v]]
  cat(sprintf("%-10s shuffle: %.3f  curve_ball: %.3f  KS p = %.3f\n", v, mean(a), mean(b),
              suppressWarnings(ks.test(a, b)$p.value)))
}
cat("\n=== Autoenlaces (original: ", sum(diag(A0)), ") ===\n", sep = "")
for (nm in names(res)[-1]) cat(sprintf("%-34s media=%.2f  min=%d  max=%d  proporcion con 0 lazos=%.2f\n",
    nm, mean(res[[nm]]$Loops), min(res[[nm]]$Loops), max(res[[nm]]$Loops), mean(res[[nm]]$Loops == 0)))
