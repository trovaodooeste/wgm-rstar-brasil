# ==============================================================================
# RESULTADOS PARA O ORIENTADOR - COM INTERVALOS DE CONFIANCA
# ==============================================================================
# Gera:
# 1. Grafico r* com IC 95% (Hamilton 1986 parametrico)
# 2. Grafico decomposicao g+z com IC
# 3. Tabela completa de parametros com SE
# 4. Tabela r* em datas-chave com IC
# ==============================================================================

library(mFilter)

wgm_path <- "."

# Carregar resultado
resultado <- readRDS(file.path(wgm_path, "outputs/results/wgm_kappa_params_fixos_20260212.rds"))

stage3 <- resultado$stage3
theta <- stage3$theta
param.num <- stage3$param.num
y.data <- stage3$y.data
x.data <- stage3$x.data
lambda.g <- stage3$lambda.g
kappa.vec <- stage3$kappa.vec
c.param <- stage3$c.param
T.obs <- nrow(y.data)
anos.est <- stage3$anos.est
tri.est <- stage3$trimestres.est
datas <- as.Date(paste0(anos.est, "-", (tri.est - 1) * 3 + 1, "-01"))

# ==============================================================================
# RECONSTRUIR xi.00, P.00, bounds
# ==============================================================================

log.output <- as.numeric(resultado$dados$ln_pib)
T.total <- length(log.output) - 8
y.og <- log.output[5:(T.total+8)]
g.pot <- hpfilter(y.og, freq = 36000, type = "lambda", drift = FALSE)$trend
g.pot.diff <- diff(g.pot)
z.initial <- resultado$z.initial

xi.00 <- c(100*g.pot[4], 100*g.pot[3], 0, 0,
            100*mean(g.pot.diff), 100*mean(g.pot.diff),
            z.initial, z.initial)
P.00 <- diag(1, 8)

theta.lb <- rep(-Inf, 12); theta.ub <- rep(Inf, 12)
theta.lb[1] <- theta.ub[1] <- theta[1]
theta.lb[2] <- theta.ub[2] <- theta[2]
theta.lb[3] <- theta.ub[3] <- theta[3]
theta.lb[4] <- theta.ub[4] <- theta[4]
theta.lb[5] <- theta.ub[5] <- theta[5]
theta.lb[6] <- theta.ub[6] <- theta[6]
theta.lb[12] <- theta.ub[12] <- theta[12]
theta.lb[7] <- 0.001; theta.lb[8] <- 0.001
theta.lb[9] <- 0.05; theta.ub[9] <- 0.5
theta.lb[10] <- theta.ub[10] <- lambda.g * theta[9]
theta.lb[11] <- 0.001; theta.ub[11] <- 2.2

params.free <- which(theta.lb != theta.ub)
n.free <- length(params.free)

# ==============================================================================
# KALMAN FILTER COM ll.vec
# ==============================================================================

kalman.filter.with.llvec <- function(theta, y.data, x.data, xi.00, P.00,
                                      param.num, lambda.g, kappa.vec) {
  a.1 <- theta[param.num["a_1"]]; a.2 <- theta[param.num["a_2"]]
  a.r <- theta[param.num["a_r"]]; b.pi1 <- theta[param.num["b_pi1"]]
  b.pi2 <- theta[param.num["b_pi2"]]; b.y <- theta[param.num["b_y"]]
  sigma.is <- theta[param.num["sigma_is"]]; sigma.pi <- theta[param.num["sigma_pi"]]
  sigma.ystar <- theta[param.num["sigma_ystar"]]; c.val <- theta[param.num["c"]]
  sigma.g <- lambda.g * sigma.ystar
  sigma.z <- theta[param.num["sigma_z"]]

  n.states <- 8
  F <- matrix(0, n.states, n.states)
  F[1,1] <- 1; F[1,5] <- 1; F[2,1] <- 1; F[4,3] <- 1
  F[5,5] <- 1; F[6,5] <- 1; F[7,7] <- 1; F[8,7] <- 1
  Q <- matrix(0, n.states, n.states)
  Q[1,1] <- sigma.ystar^2; Q[5,5] <- sigma.g^2; Q[7,7] <- sigma.z^2
  R.base <- diag(c(sigma.is^2, sigma.pi^2))
  T.obs <- nrow(y.data)
  xi.tt <- xi.00; P.tt <- P.00
  ll <- 0; ll.vec <- numeric(T.obs)
  xi.filtered <- matrix(0, T.obs, n.states)

  for (t in 1:T.obs) {
    R <- R.base * (kappa.vec[t]^2)
    xi.ttm1 <- F %*% xi.tt; P.ttm1 <- F %*% P.tt %*% t(F) + Q
    y.t <- y.data[t, ]
    rstar.t <- c.val * xi.ttm1[5] / 100 + xi.ttm1[7]
    rstar.tm1 <- c.val * xi.ttm1[6] / 100 + xi.ttm1[8]
    r.gap <- (x.data[t, 3] + x.data[t, 4]) / 2 - (rstar.t + rstar.tm1) / 2
    y.pred.is <- xi.ttm1[1] + a.1 * (x.data[t, 1] - xi.ttm1[2]) +
                 a.2 * (x.data[t, 2] - xi.ttm1[2]) + a.r * r.gap
    pi.pred <- b.pi1 * x.data[t, 5] + b.pi2 * x.data[t, 6] +
               (1 - b.pi1 - b.pi2) * x.data[t, 7] + b.y * xi.ttm1[3]
    v.t <- y.t - c(y.pred.is, pi.pred)
    H <- matrix(0, n.states, 2)
    H[1, 1] <- 1; H[2, 1] <- -(a.1 + a.2)
    H[5, 1] <- -a.r * c.val / 200; H[6, 1] <- -a.r * c.val / 200
    H[7, 1] <- -a.r / 2; H[8, 1] <- -a.r / 2; H[3, 2] <- b.y
    S.t <- t(H) %*% P.ttm1 %*% H + R + diag(1e-8, 2)
    K.t <- P.ttm1 %*% H %*% solve(S.t)
    xi.tt <- xi.ttm1 + K.t %*% v.t
    P.tt <- P.ttm1 - K.t %*% t(H) %*% P.ttm1
    ll.t <- -0.5 * (log(det(S.t)) + t(v.t) %*% solve(S.t) %*% v.t)
    ll.vec[t] <- as.numeric(ll.t); ll <- ll + ll.t
    xi.filtered[t, ] <- xi.tt
  }
  return(list(ll = as.numeric(ll), ll.vec = ll.vec,
              xi.filtered = xi.filtered, P.tt = P.tt))
}

# ==============================================================================
# MATRIZ DE INFORMACAO
# ==============================================================================

cat("Computando matriz de informacao...\n")
kf.base <- kalman.filter.with.llvec(theta, y.data, x.data, xi.00, P.00,
                                     param.num, lambda.g, kappa.vec)
ll.vec.base <- kf.base$ll.vec
likelihood.gradient <- matrix(NA, T.obs, length(theta))

for (i in params.free) {
  delta <- max(abs(theta[i]) * 1e-6, 1e-6)
  d.theta <- theta; d.theta[i] <- theta[i] + delta
  kf.d <- kalman.filter.with.llvec(d.theta, y.data, x.data, xi.00, P.00,
                                    param.num, lambda.g, kappa.vec)
  likelihood.gradient[, i] <- (kf.d$ll.vec - ll.vec.base) / delta
}

likelihood.gradient.free <- likelihood.gradient[, params.free, drop = FALSE]
info <- tryCatch(
  solve(t(likelihood.gradient.free) %*% likelihood.gradient.free),
  error = function(e) MASS::ginv(t(likelihood.gradient.free) %*% likelihood.gradient.free)
)
bse <- sqrt(diag(info))

# ==============================================================================
# MONTE CARLO
# ==============================================================================

niter <- 5000
cat(sprintf("Monte Carlo: %d iteracoes...\n", niter))

eigenstuff <- eigen(info)
eigenvectors <- eigenstuff$vectors
for (l in 1:n.free) {
  if (eigenvectors[1, l] < 0) eigenvectors[, l] <- -eigenvectors[, l]
}
eigenvalues <- eigenstuff$values
eigenvalues[eigenvalues <= 0] <- 1e-10
hh <- eigenvectors %*% sqrt(diag(eigenvalues, nrow = n.free))

g.base <- kf.base$xi.filtered[, 5] * 4 / 100
z.base <- kf.base$xi.filtered[, 7]
rstar.base <- c.param * g.base * 100 + z.base
ypot.base <- kf.base$xi.filtered[, 1] / 100

cum1 <- matrix(0, T.obs, 4)  # [y*, r*, g, z]
set.seed(50)
good.draws <- 0; excluded.draws <- 0
t.start <- Sys.time()

while (good.draws < niter) {
  perturbation <- rep(0, length(theta))
  perturbation[params.free] <- as.vector(hh %*% rnorm(n.free))
  theta.i <- theta + perturbation

  valid <- TRUE
  for (j in params.free) {
    if (theta.i[j] < theta.lb[j] || theta.i[j] > theta.ub[j]) { valid <- FALSE; break }
  }

  if (valid) {
    kf.i <- tryCatch(
      kalman.filter.with.llvec(theta.i, y.data, x.data, xi.00, P.00,
                                param.num, lambda.g, kappa.vec),
      error = function(e) NULL)
    if (!is.null(kf.i)) {
      g.i <- kf.i$xi.filtered[, 5] * 4 / 100
      z.i <- kf.i$xi.filtered[, 7]
      rstar.i <- theta.i[param.num["c"]] * g.i * 100 + z.i
      ypot.i <- kf.i$xi.filtered[, 1] / 100
      cum1[, 1] <- cum1[, 1] + (ypot.i - ypot.base)^2
      cum1[, 2] <- cum1[, 2] + (rstar.i - rstar.base)^2
      cum1[, 3] <- cum1[, 3] + (g.i * 100 - g.base * 100)^2
      cum1[, 4] <- cum1[, 4] + (z.i - z.base)^2
      good.draws <- good.draws + 1
      if (good.draws %% 1000 == 0) cat(sprintf("  %d/%d (%.0fs)\n", good.draws, niter,
          as.numeric(difftime(Sys.time(), t.start, units="secs"))))
    } else { excluded.draws <- excluded.draws + 1 }
  } else { excluded.draws <- excluded.draws + 1 }
}

elapsed <- as.numeric(difftime(Sys.time(), t.start, units = "secs"))
cat(sprintf("Concluido em %.0fs (%d excluidos)\n", elapsed, excluded.draws))

cum1 <- cum1 / niter
se <- sqrt(cum1)
se.rstar <- se[, 2]; se.g <- se[, 3]; se.z <- se[, 4]

# ==============================================================================
# TABELA 1: PARAMETROS COM IC
# ==============================================================================

cat("\n")
cat("================================================================\n")
cat("TABELA 1: PARAMETROS ESTIMADOS\n")
cat("================================================================\n")
cat("Parametro       | Valor     | SE       | IC 95%              | Status\n")
cat("----------------|-----------|----------|---------------------|--------\n")

# Fixos (do modelo original)
cat(sprintf("a_1             | %9.4f |    -     |         -           | Fixo\n", theta[1]))
cat(sprintf("a_2             | %9.4f |    -     |         -           | Fixo\n", theta[2]))
cat(sprintf("a_r             | %9.4f |    -     |         -           | Fixo\n", theta[3]))
cat(sprintf("b_pi1           | %9.4f |    -     |         -           | Fixo\n", theta[4]))
cat(sprintf("b_pi2           | %9.4f |    -     |         -           | Fixo\n", theta[5]))
cat(sprintf("b_y             | %9.4f |    -     |         -           | Fixo\n", theta[6]))
cat(sprintf("c               | %9.1f |    -     |         -           | Fixo\n", theta[12]))

# Estimados com IC
free.names <- c("sigma_IS", "sigma_pi", "sigma_y*", "sigma_z")
for (k in 1:n.free) {
  i <- params.free[k]
  lo <- theta[i] - 1.96 * bse[k]
  hi <- theta[i] + 1.96 * bse[k]
  cat(sprintf("%-16s| %9.4f | %8.4f | [%6.3f, %6.3f]    | MLE\n",
      free.names[k], theta[i], bse[k], lo, hi))
}

# Kappas
kappas <- resultado$kappas
cat(sprintf("kappa_2020      | %9.2f |    -     |         -           | MLE S1\n", kappas["kappa_2020"]))
cat(sprintf("kappa_2021      | %9.2f |    -     |         -           | MLE S1\n", kappas["kappa_2021"]))
cat(sprintf("kappa_2022      | %9.2f |    -     |         -           | MLE S1\n", kappas["kappa_2022"]))

# ==============================================================================
# TABELA 2: r* EM DATAS-CHAVE COM IC
# ==============================================================================

cat("\n")
cat("================================================================\n")
cat("TABELA 2: r* EM PERIODOS-CHAVE\n")
cat("================================================================\n")
cat("Periodo         | r*       | SE      | IC 95%\n")
cat("----------------|----------|---------|---------------------\n")

periodos <- list(
  "2005-2008" = anos.est >= 2005 & anos.est <= 2008,
  "2009-2013" = anos.est >= 2009 & anos.est <= 2013,
  "2014-2019" = anos.est >= 2014 & anos.est <= 2019,
  "2020-2022" = anos.est >= 2020 & anos.est <= 2022,
  "2023-atual" = anos.est >= 2023,
  "Amostra toda" = rep(TRUE, T.obs)
)

for (p in names(periodos)) {
  idx <- periodos[[p]]
  r.mean <- mean(rstar.base[idx])
  se.mean <- mean(se.rstar[idx])
  cat(sprintf("%-16s| %7.2f%% | %6.2f  | [%.2f%%, %.2f%%]\n",
      p, r.mean, se.mean, r.mean - 1.96*se.mean, r.mean + 1.96*se.mean))
}

# Ultimo trimestre
cat(sprintf("%-16s| %7.2f%% | %6.2f  | [%.2f%%, %.2f%%]\n",
    paste0(tail(anos.est,1), "-Q", tail(tri.est,1)),
    tail(rstar.base, 1), tail(se.rstar, 1),
    tail(rstar.base, 1) - 1.96*tail(se.rstar, 1),
    tail(rstar.base, 1) + 1.96*tail(se.rstar, 1)))

# ==============================================================================
# TABELA 3: DECOMPOSICAO r* = g + z COM IC
# ==============================================================================

cat("\n")
cat("================================================================\n")
cat("TABELA 3: DECOMPOSICAO r* = g + z (valores atuais)\n")
cat("================================================================\n")
cat(sprintf("  g atual:  %.2f%% [IC: %.2f%% - %.2f%%]\n",
    tail(g.base*100, 1),
    tail(g.base*100, 1) - 1.96*tail(se.g, 1),
    tail(g.base*100, 1) + 1.96*tail(se.g, 1)))
cat(sprintf("  z atual:  %.2f%% [IC: %.2f%% - %.2f%%]\n",
    tail(z.base, 1),
    tail(z.base, 1) - 1.96*tail(se.z, 1),
    tail(z.base, 1) + 1.96*tail(se.z, 1)))
cat(sprintf("  r* atual: %.2f%% [IC: %.2f%% - %.2f%%]\n",
    tail(rstar.base, 1),
    tail(rstar.base, 1) - 1.96*tail(se.rstar, 1),
    tail(rstar.base, 1) + 1.96*tail(se.rstar, 1)))

# ==============================================================================
# GRAFICO 1: r* COM IC 95% (principal)
# ==============================================================================

data.arquivo <- format(Sys.Date(), "%Y%m%d")

png(file.path(wgm_path, paste0("outputs/graficos/fig_rstar_ic95_", data.arquivo, ".png")),
    width = 1200, height = 700, res = 130)

par(mar = c(4, 4.5, 3, 1))

rstar.upper <- rstar.base + 1.96 * se.rstar
rstar.lower <- rstar.base - 1.96 * se.rstar

plot(datas, rstar.base, type = "l", col = "red", lwd = 2.5,
     main = expression(paste("Taxa Natural de Juros (", r^"*", ") com IC 95% Parametrico")),
     ylab = "% a.a.", xlab = "",
     ylim = range(c(rstar.lower * 0.95, rstar.upper * 1.05)))

polygon(c(datas, rev(datas)),
        c(rstar.upper, rev(rstar.lower)),
        col = rgb(1, 0, 0, 0.15), border = NA)

lines(datas, rstar.base, col = "red", lwd = 2.5)
lines(datas, rstar.upper, col = "red", lty = 2, lwd = 0.7)
lines(datas, rstar.lower, col = "red", lty = 2, lwd = 0.7)

# Taxa real observada
taxa.real <- as.numeric(resultado$dados$taxa_real_exante_tri[
  match(paste(anos.est, tri.est), paste(resultado$dados$ano, resultado$dados$trimestre))])
lines(datas, taxa.real, col = "gray50", lwd = 1, lty = 1)

# Sombreamento COVID
rect(as.Date("2020-04-01"), par("usr")[3],
     as.Date("2022-12-01"), par("usr")[4],
     col = rgb(0.5, 0.5, 0.5, 0.08), border = NA)

# Anotacoes
text(as.Date("2021-01-01"), par("usr")[4] * 0.95, "COVID", cex = 0.7, col = "gray40")

legend("topright",
       c(expression(r^"*"), "IC 95%", "Taxa real observada"),
       col = c("red", rgb(1, 0, 0, 0.3), "gray50"),
       lwd = c(2.5, 8, 1), lty = c(1, 1, 1), bty = "n", cex = 0.85)

dev.off()
cat(sprintf("\nGrafico 1 salvo: fig_rstar_ic95_%s.png\n", data.arquivo))

# ==============================================================================
# GRAFICO 2: DECOMPOSICAO g + z COM IC
# ==============================================================================

png(file.path(wgm_path, paste0("outputs/graficos/fig_decomposicao_ic95_", data.arquivo, ".png")),
    width = 1400, height = 700, res = 130)

par(mfrow = c(1, 2), mar = c(4, 4.5, 3, 1))

# g com IC
g.pct <- g.base * 100
g.upper <- g.pct + 1.96 * se.g
g.lower <- g.pct - 1.96 * se.g
plot(datas, g.pct, type = "l", col = "blue", lwd = 2.5,
     main = "g (crescimento potencial) com IC 95%",
     ylab = "% a.a.", xlab = "",
     ylim = range(c(g.lower, g.upper)))
polygon(c(datas, rev(datas)), c(g.upper, rev(g.lower)),
        col = rgb(0, 0, 1, 0.12), border = NA)
lines(datas, g.pct, col = "blue", lwd = 2.5)
abline(h = 0, col = "gray", lty = 3)

# z com IC
z.upper <- z.base + 1.96 * se.z
z.lower <- z.base - 1.96 * se.z
plot(datas, z.base, type = "l", col = "darkgreen", lwd = 2.5,
     main = "z (premio de risco e outros) com IC 95%",
     ylab = "% a.a.", xlab = "",
     ylim = range(c(z.lower, z.upper)))
polygon(c(datas, rev(datas)), c(z.upper, rev(z.lower)),
        col = rgb(0, 0.5, 0, 0.12), border = NA)
lines(datas, z.base, col = "darkgreen", lwd = 2.5)

dev.off()
cat(sprintf("Grafico 2 salvo: fig_decomposicao_ic95_%s.png\n", data.arquivo))

# ==============================================================================
# GRAFICO 3: DASHBOARD COMPLETO COM IC
# ==============================================================================

png(file.path(wgm_path, paste0("outputs/graficos/dashboard_ic95_", data.arquivo, ".png")),
    width = 1600, height = 1000, res = 130)

layout(matrix(c(1,1,2,3,4,5), 3, 2, byrow = TRUE), heights = c(1.3, 1, 1))
par(mar = c(4, 4.5, 3, 1))

# Painel 1 (grande): r* com IC
plot(datas, rstar.base, type = "l", col = "red", lwd = 2.5,
     main = expression(paste(r^"*", " com Intervalo de Confianca 95% Parametrico")),
     ylab = "% a.a.", xlab = "",
     ylim = range(c(rstar.lower * 0.95, rstar.upper * 1.05)))
polygon(c(datas, rev(datas)), c(rstar.upper, rev(rstar.lower)),
        col = rgb(1, 0, 0, 0.15), border = NA)
lines(datas, rstar.base, col = "red", lwd = 2.5)
lines(datas, taxa.real, col = "gray50", lwd = 1)
rect(as.Date("2020-04-01"), par("usr")[3],
     as.Date("2022-12-01"), par("usr")[4],
     col = rgb(0.5, 0.5, 0.5, 0.08), border = NA)
legend("topright",
       c(expression(r^"*"), "IC 95%", "Taxa real"),
       col = c("red", rgb(1, 0, 0, 0.3), "gray50"),
       lwd = c(2.5, 8, 1), bty = "n", cex = 0.85)

# Painel 2: g com IC
plot(datas, g.pct, type = "l", col = "blue", lwd = 2,
     main = "g (crescimento potencial)", ylab = "% a.a.", xlab = "",
     ylim = range(c(g.lower, g.upper)))
polygon(c(datas, rev(datas)), c(g.upper, rev(g.lower)),
        col = rgb(0, 0, 1, 0.12), border = NA)
lines(datas, g.pct, col = "blue", lwd = 2)

# Painel 3: z com IC
plot(datas, z.base, type = "l", col = "darkgreen", lwd = 2,
     main = "z (outros fatores)", ylab = "% a.a.", xlab = "",
     ylim = range(c(z.lower, z.upper)))
polygon(c(datas, rev(datas)), c(z.upper, rev(z.lower)),
        col = rgb(0, 0.5, 0, 0.12), border = NA)
lines(datas, z.base, col = "darkgreen", lwd = 2)

# Painel 4: Output gap
gap <- stage3$output.gap.smoothed
plot(datas, gap, type = "l", col = "purple", lwd = 2,
     main = "Hiato do produto", ylab = "%", xlab = "")
abline(h = 0, col = "gray", lty = 2)
rect(as.Date("2020-04-01"), par("usr")[3],
     as.Date("2022-12-01"), par("usr")[4],
     col = rgb(0.5, 0.5, 0.5, 0.08), border = NA)

# Painel 5: SE ao longo do tempo
plot(datas, se.rstar, type = "l", col = "darkred", lwd = 2,
     main = expression(paste("Erro Padrao de ", r^"*")),
     ylab = "SE (pp)", xlab = "")
abline(h = mean(se.rstar), col = "gray", lty = 2)
text(datas[5], mean(se.rstar) * 1.15,
     paste0("Media: ", round(mean(se.rstar), 2), " pp"),
     col = "gray30", cex = 0.8)

dev.off()
cat(sprintf("Grafico 3 salvo: dashboard_ic95_%s.png\n", data.arquivo))

# ==============================================================================
# SALVAR RESULTADOS DOS ICs (para nao ter que re-rodar o MC)
# ==============================================================================

ic.results <- list(
  # SEs por trimestre
  se.rstar = se.rstar,
  se.g = se.g,
  se.z = se.z,
  se.ypot = se[, 1],
  # Estimativas base
  rstar.base = rstar.base,
  g.base = g.base,
  z.base = z.base,
  # ICs 95%
  rstar.upper = rstar.base + 1.96 * se.rstar,
  rstar.lower = rstar.base - 1.96 * se.rstar,
  g.upper = g.base * 100 + 1.96 * se.g,
  g.lower = g.base * 100 - 1.96 * se.g,
  z.upper = z.base + 1.96 * se.z,
  z.lower = z.base - 1.96 * se.z,
  # Datas e metadados
  datas = datas,
  anos.est = anos.est,
  tri.est = tri.est,
  taxa.real = taxa.real,
  gap = stage3$output.gap.smoothed,
  # Info matrix
  bse = bse,
  params.free = params.free,
  info = info,
  # Config
  niter = niter,
  n.excluded = excluded.draws,
  elapsed.secs = elapsed
)

saveRDS(ic.results, file.path(wgm_path, paste0("outputs/results/ic_hamilton_parametrico_", data.arquivo, ".rds")))
cat(sprintf("\nResultados IC salvos: outputs/results/ic_hamilton_parametrico_%s.rds\n", data.arquivo))

cat("\n================================================================\n")
cat("TUDO PRONTO PARA ENVIAR AO ORIENTADOR\n")
cat("================================================================\n")
cat("Graficos:\n")
cat(sprintf("  1. fig_rstar_ic95_%s.png (r* com IC - principal)\n", data.arquivo))
cat(sprintf("  2. fig_decomposicao_ic95_%s.png (g e z com IC)\n", data.arquivo))
cat(sprintf("  3. dashboard_ic95_%s.png (dashboard completo)\n", data.arquivo))
cat("\nTabelas impressas acima.\n")
