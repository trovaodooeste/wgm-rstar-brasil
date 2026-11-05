# ==============================================================================
# HAMILTON (1986) COMPLETO: cum1 + cum2
# ==============================================================================
# Implementa o procedimento Hamilton completo como no HLW (2023):
#   SE = sqrt(cum1/N + cum2/N)
#   cum1 = Var_theta(E[s|theta])  — incerteza parametrica
#   cum2 = E_theta(Var[s|theta])  — incerteza do estado (via P_{t|T})
#
# Diferenca do script anterior (gerar_resultados_orientador_ic.R):
#   - Antigo: só cum1 (incerteza parametrica), usa estados filtrados
#   - Novo: cum1 + cum2, usa estados SUAVIZADOS (RTS smoother)
#
# Para cada draw Monte Carlo:
#   1. Perturba parametros via info matrix
#   2. Roda Kalman filter forward
#   3. Roda RTS smoother backward
#   4. cum1 += (r*_i_suavizado - r*_base_suavizado)^2
#   5. cum2 += entradas relevantes de P_{t|T}
# ==============================================================================

library(mFilter)

wgm_path <- "."
source(file.path(wgm_path, "scripts/kalman_smoother.R"))

# Carregar resultado salvo (para nao re-rodar Stages 1-3)
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

cat("\n================================================================\n")
cat("HAMILTON (1986) COMPLETO: cum1 + cum2\n")
cat("================================================================\n")
cat(sprintf("  Parametros livres: %d (%s)\n", n.free,
    paste(names(param.num)[params.free], collapse=", ")))
cat(sprintf("  Periodos: %d\n", T.obs))

# ==============================================================================
# KALMAN FILTER COMPLETO (salva tudo para o smoother)
# ==============================================================================

kalman.filter.full <- function(theta, y.data, x.data, xi.00, P.00,
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

  # Arrays para filter E smoother
  xi.tt.array <- matrix(0, T.obs, n.states)
  xi.ttm1.array <- matrix(0, T.obs, n.states)
  P.tt.array <- array(0, dim = c(n.states, n.states, T.obs))
  P.ttm1.array <- array(0, dim = c(n.states, n.states, T.obs))

  for (t in 1:T.obs) {
    R <- R.base * (kappa.vec[t]^2)

    # Predict
    xi.ttm1 <- F %*% xi.tt
    P.ttm1 <- F %*% P.tt %*% t(F) + Q

    # Salvar previstos
    xi.ttm1.array[t, ] <- xi.ttm1
    P.ttm1.array[, , t] <- P.ttm1

    # Observacao
    y.t <- y.data[t, ]
    rstar.t <- c.val * xi.ttm1[5] / 100 + xi.ttm1[7]
    rstar.tm1 <- c.val * xi.ttm1[6] / 100 + xi.ttm1[8]
    r.gap <- (x.data[t, 3] + x.data[t, 4]) / 2 - (rstar.t + rstar.tm1) / 2
    y.pred.is <- xi.ttm1[1] + a.1 * (x.data[t, 1] - xi.ttm1[2]) +
                 a.2 * (x.data[t, 2] - xi.ttm1[2]) + a.r * r.gap
    pi.pred <- b.pi1 * x.data[t, 5] + b.pi2 * x.data[t, 6] +
               (1 - b.pi1 - b.pi2) * x.data[t, 7] + b.y * xi.ttm1[3]
    v.t <- y.t - c(y.pred.is, pi.pred)

    # H matrix
    H <- matrix(0, n.states, 2)
    H[1, 1] <- 1; H[2, 1] <- -(a.1 + a.2)
    H[5, 1] <- -a.r * c.val / 200; H[6, 1] <- -a.r * c.val / 200
    H[7, 1] <- -a.r / 2; H[8, 1] <- -a.r / 2; H[3, 2] <- b.y

    # Kalman gain e update
    S.t <- t(H) %*% P.ttm1 %*% H + R + diag(1e-8, 2)
    K.t <- P.ttm1 %*% H %*% solve(S.t)
    xi.tt <- xi.ttm1 + K.t %*% v.t
    P.tt <- P.ttm1 - K.t %*% t(H) %*% P.ttm1

    # Log-likelihood
    ll.t <- -0.5 * (log(det(S.t)) + t(v.t) %*% solve(S.t) %*% v.t)
    ll.vec[t] <- as.numeric(ll.t)
    ll <- ll + ll.t

    # Salvar filtrados
    xi.tt.array[t, ] <- xi.tt
    P.tt.array[, , t] <- P.tt
  }

  return(list(ll = as.numeric(ll), ll.vec = ll.vec,
              xi.tt.array = xi.tt.array, P.tt.array = P.tt.array,
              xi.ttm1.array = xi.ttm1.array, P.ttm1.array = P.ttm1.array,
              F.mat = F, Q = Q))
}

# ==============================================================================
# RODAR FILTER + SMOOTHER BASE
# ==============================================================================

cat("\nRodando filter + smoother base...\n")
kf.base <- kalman.filter.full(theta, y.data, x.data, xi.00, P.00,
                               param.num, lambda.g, kappa.vec)

sm.base <- kalman.smoother.rts(kf.base$xi.tt.array, kf.base$P.tt.array,
                                kf.base$xi.ttm1.array, kf.base$P.ttm1.array,
                                kf.base$F.mat)

# r* base a partir dos estados SUAVIZADOS
g.base.smooth <- sm.base$xi.tT[, 5] * 4 / 100     # fração anual
z.base.smooth <- sm.base$xi.tT[, 7]                 # % anual
rstar.base.smooth <- c.param * g.base.smooth * 100 + z.base.smooth  # % anual
ypot.base.smooth <- sm.base$xi.tT[, 1] / 100

# Comparar filtrado vs suavizado
g.base.filt <- kf.base$xi.tt.array[, 5] * 4 / 100
z.base.filt <- kf.base$xi.tt.array[, 7]
rstar.base.filt <- c.param * g.base.filt * 100 + z.base.filt

cat(sprintf("  r* base (filtrado):  media=%.2f%%, atual=%.2f%%\n",
    mean(rstar.base.filt), tail(rstar.base.filt, 1)))
cat(sprintf("  r* base (suavizado): media=%.2f%%, atual=%.2f%%\n",
    mean(rstar.base.smooth), tail(rstar.base.smooth, 1)))

# Verificar que o smoother reduz P_{t|T} vs P_{t|t}
p77.filt <- sapply(1:T.obs, function(t) kf.base$P.tt.array[7,7,t])
p77.smooth <- sapply(1:T.obs, function(t) sm.base$P.tT[7,7,t])
cat(sprintf("\n  P[7,7] filtrado:  inicio=%.1f, meio=%.1f, fim=%.1f\n",
    p77.filt[1], p77.filt[T.obs %/% 2], p77.filt[T.obs]))
cat(sprintf("  P[7,7] suavizado: inicio=%.1f, meio=%.1f, fim=%.1f\n",
    p77.smooth[1], p77.smooth[T.obs %/% 2], p77.smooth[T.obs]))
cat(sprintf("  Reducao no meio: %.0f%%\n",
    100 * (1 - p77.smooth[T.obs %/% 2] / p77.filt[T.obs %/% 2])))

# ==============================================================================
# MATRIZ DE INFORMACAO
# ==============================================================================

cat("\nComputando matriz de informacao...\n")
ll.vec.base <- kf.base$ll.vec
likelihood.gradient <- matrix(NA, T.obs, length(theta))

for (i in params.free) {
  delta <- max(abs(theta[i]) * 1e-6, 1e-6)
  d.theta <- theta; d.theta[i] <- theta[i] + delta
  kf.d <- kalman.filter.full(d.theta, y.data, x.data, xi.00, P.00,
                              param.num, lambda.g, kappa.vec)
  likelihood.gradient[, i] <- (kf.d$ll.vec - ll.vec.base) / delta
}

likelihood.gradient.free <- likelihood.gradient[, params.free, drop = FALSE]
info <- tryCatch(
  solve(t(likelihood.gradient.free) %*% likelihood.gradient.free),
  error = function(e) MASS::ginv(t(likelihood.gradient.free) %*% likelihood.gradient.free)
)
bse <- sqrt(diag(info))

cat("  Parametros livres com SE:\n")
free.names <- names(param.num)[params.free]
for (k in 1:n.free) {
  cat(sprintf("    %-10s = %.4f (SE=%.4f, t=%.2f)\n",
      free.names[k], theta[params.free[k]], bse[k],
      abs(theta[params.free[k]]) / bse[k]))
}

# ==============================================================================
# MONTE CARLO: HAMILTON COMPLETO (cum1 + cum2)
# ==============================================================================

niter <- 5000
cat(sprintf("\nMonte Carlo Hamilton completo: %d iteracoes...\n", niter))
cat("  Para cada draw: filter forward + smoother backward\n")

eigenstuff <- eigen(info)
eigenvectors <- eigenstuff$vectors
for (l in 1:n.free) {
  if (eigenvectors[1, l] < 0) eigenvectors[, l] <- -eigenvectors[, l]
}
eigenvalues <- eigenstuff$values
eigenvalues[eigenvalues <= 0] <- 1e-10
hh <- eigenvectors %*% sqrt(diag(eigenvalues, nrow = n.free))

# cum1: Var_theta(E[s|theta]) — variancia dos estados suavizados entre draws
# cum2: E_theta(Var[s|theta]) — media da P_{t|T} entre draws
cum1 <- matrix(0, T.obs, 4)  # [y*, r*, g, z]
cum2 <- matrix(0, T.obs, 4)  # [y*, r*, g, z]

set.seed(50)
good.draws <- 0
excluded.draws <- 0
excluded.constraint <- 0
excluded.error <- 0
t.start <- Sys.time()

while (good.draws < niter) {
  # 1. Perturbar parametros
  perturbation <- rep(0, length(theta))
  perturbation[params.free] <- as.vector(hh %*% rnorm(n.free))
  theta.i <- theta + perturbation

  # 2. Checar bounds
  valid <- TRUE
  for (j in params.free) {
    if (theta.i[j] < theta.lb[j] || theta.i[j] > theta.ub[j]) {
      valid <- FALSE; break
    }
  }

  if (!valid) {
    excluded.constraint <- excluded.constraint + 1
    excluded.draws <- excluded.draws + 1
    next
  }

  # 3. Rodar filter forward
  kf.i <- tryCatch(
    kalman.filter.full(theta.i, y.data, x.data, xi.00, P.00,
                        param.num, lambda.g, kappa.vec),
    error = function(e) NULL)

  if (is.null(kf.i)) {
    excluded.error <- excluded.error + 1
    excluded.draws <- excluded.draws + 1
    next
  }

  # 4. Rodar RTS smoother backward
  sm.i <- tryCatch(
    kalman.smoother.rts(kf.i$xi.tt.array, kf.i$P.tt.array,
                         kf.i$xi.ttm1.array, kf.i$P.ttm1.array,
                         kf.i$F.mat),
    error = function(e) NULL)

  if (is.null(sm.i)) {
    excluded.error <- excluded.error + 1
    excluded.draws <- excluded.draws + 1
    next
  }

  # 5. Extrair estados suavizados
  g.i <- sm.i$xi.tT[, 5] * 4 / 100
  z.i <- sm.i$xi.tT[, 7]
  rstar.i <- theta.i[param.num["c"]] * g.i * 100 + z.i
  ypot.i <- sm.i$xi.tT[, 1] / 100

  # 6. Acumular cum1 (incerteza parametrica: desvio dos estados suavizados)
  cum1[, 1] <- cum1[, 1] + (ypot.i - ypot.base.smooth)^2
  cum1[, 2] <- cum1[, 2] + (rstar.i - rstar.base.smooth)^2
  cum1[, 3] <- cum1[, 3] + (g.i * 100 - g.base.smooth * 100)^2
  cum1[, 4] <- cum1[, 4] + (z.i - z.base.smooth)^2

  # 7. Acumular cum2 (incerteza do estado: P_{t|T} suavizada)
  for (t in 1:T.obs) {
    P.tT.t <- sm.i$P.tT[, , t]

    # y*: Var(y*) = P[1,1]
    cum2[t, 1] <- cum2[t, 1] + P.tT.t[1, 1]

    # r* = 4c * xi[5] + xi[7]
    # Var(r*) = (4c)^2 * P[5,5] + P[7,7] + 2*(4c)*P[5,7]
    cum2[t, 2] <- cum2[t, 2] +
      (4 * c.param)^2 * P.tT.t[5, 5] + P.tT.t[7, 7] +
      2 * (4 * c.param) * P.tT.t[5, 7]

    # g (em %): Var(g%) = 16 * P[5,5]
    cum2[t, 3] <- cum2[t, 3] + 16 * P.tT.t[5, 5]

    # z: Var(z) = P[7,7]
    cum2[t, 4] <- cum2[t, 4] + P.tT.t[7, 7]
  }

  good.draws <- good.draws + 1
  if (good.draws %% 500 == 0) {
    elapsed.so.far <- as.numeric(difftime(Sys.time(), t.start, units = "secs"))
    eta <- elapsed.so.far / good.draws * (niter - good.draws)
    cat(sprintf("  %d/%d (%.0fs, ETA %.0fs, excl: %d)\n",
        good.draws, niter, elapsed.so.far, eta, excluded.draws))
  }
}

elapsed <- as.numeric(difftime(Sys.time(), t.start, units = "secs"))
cat(sprintf("\nConcluido em %.0fs\n", elapsed))
cat(sprintf("  Draws aceitos: %d\n", good.draws))
cat(sprintf("  Draws excluidos: %d (constraint: %d, error: %d)\n",
    excluded.draws, excluded.constraint, excluded.error))

# ==============================================================================
# CALCULAR SEs COMPLETOS
# ==============================================================================

cum1 <- cum1 / niter
cum2 <- cum2 / niter

# SE por componente e total
se.cum1 <- sqrt(cum1)   # SE parametrico (sqrt de cada coluna)
se.cum2 <- sqrt(cum2)   # SE estado
se.total <- sqrt(cum1 + cum2)  # SE total

se.rstar <- se.total[, 2]
se.g <- se.total[, 3]
se.z <- se.total[, 4]
se.ypot <- se.total[, 1]

# SE parametrico por variavel
se.rstar.param <- se.cum1[, 2]
se.g.param <- se.cum1[, 3]
se.z.param <- se.cum1[, 4]

# SE estado por variavel
se.rstar.estado <- se.cum2[, 2]
se.g.estado <- se.cum2[, 3]
se.z.estado <- se.cum2[, 4]

# Bandas 95%
rstar.upper <- rstar.base.smooth + 1.96 * se.rstar
rstar.lower <- rstar.base.smooth - 1.96 * se.rstar
g.pct <- g.base.smooth * 100
g.upper <- g.pct + 1.96 * se.g
g.lower <- g.pct - 1.96 * se.g
z.upper <- z.base.smooth + 1.96 * se.z
z.lower <- z.base.smooth - 1.96 * se.z

# Taxa real observada
taxa.real <- as.numeric(resultado$dados$taxa_real_exante_tri[
  match(paste(anos.est, tri.est), paste(resultado$dados$ano, resultado$dados$trimestre))])

# ==============================================================================
# TABELA RESUMO NO CONSOLE
# ==============================================================================

cat("\n================================================================\n")
cat("TABELA RESUMO - HAMILTON (1986) COMPLETO\n")
cat("================================================================\n\n")

cat("                  | SE param  | SE estado | SE total  | % estado\n")
cat("                  | sqrt(cum1)| sqrt(cum2)| sqrt(1+2) |         \n")
cat("------------------|-----------|-----------|-----------|----------\n")
cat(sprintf("r* medio          | %7.2f   | %7.2f   | %7.2f   | %5.0f%%\n",
    mean(se.rstar.param), mean(se.rstar.estado), mean(se.rstar),
    100 * mean(cum2[,2]) / mean(cum1[,2] + cum2[,2])))
cat(sprintf("r* atual          | %7.2f   | %7.2f   | %7.2f   | %5.0f%%\n",
    tail(se.rstar.param, 1), tail(se.rstar.estado, 1), tail(se.rstar, 1),
    100 * tail(cum2[,2], 1) / tail(cum1[,2] + cum2[,2], 1)))
cat(sprintf("g medio           | %7.2f   | %7.2f   | %7.2f   | %5.0f%%\n",
    mean(se.g.param), mean(se.g.estado), mean(se.g),
    100 * mean(cum2[,3]) / mean(cum1[,3] + cum2[,3])))
cat(sprintf("g atual           | %7.2f   | %7.2f   | %7.2f   | %5.0f%%\n",
    tail(se.g.param, 1), tail(se.g.estado, 1), tail(se.g, 1),
    100 * tail(cum2[,3], 1) / tail(cum1[,3] + cum2[,3], 1)))
cat(sprintf("z medio           | %7.2f   | %7.2f   | %7.2f   | %5.0f%%\n",
    mean(se.z.param), mean(se.z.estado), mean(se.z),
    100 * mean(cum2[,4]) / mean(cum1[,4] + cum2[,4])))
cat(sprintf("z atual           | %7.2f   | %7.2f   | %7.2f   | %5.0f%%\n",
    tail(se.z.param, 1), tail(se.z.estado, 1), tail(se.z, 1),
    100 * tail(cum2[,4], 1) / tail(cum1[,4] + cum2[,4], 1)))

cat("\n--- r* COM IC 95% (Hamilton completo) ---\n")
cat(sprintf("  r* atual:  %.2f%% [%.2f%%, %.2f%%]\n",
    tail(rstar.base.smooth, 1),
    tail(rstar.base.smooth, 1) - 1.96 * tail(se.rstar, 1),
    tail(rstar.base.smooth, 1) + 1.96 * tail(se.rstar, 1)))
cat(sprintf("  r* medio:  %.2f%% [%.2f%%, %.2f%%]\n",
    mean(rstar.base.smooth),
    mean(rstar.base.smooth) - 1.96 * mean(se.rstar),
    mean(rstar.base.smooth) + 1.96 * mean(se.rstar)))

cat("\n--- DECOMPOSICAO r* = g + z (valores atuais com IC) ---\n")
cat(sprintf("  g atual:  %.2f%% [%.2f%%, %.2f%%]\n",
    tail(g.pct, 1), tail(g.lower, 1), tail(g.upper, 1)))
cat(sprintf("  z atual:  %.2f%% [%.2f%%, %.2f%%]\n",
    tail(z.base.smooth, 1), tail(z.lower, 1), tail(z.upper, 1)))

# Comparacao com Hamilton parametrico (cum1 only)
ic.param <- tryCatch(
  readRDS(file.path(wgm_path, "outputs/results/ic_hamilton_parametrico_selic_20260213.rds")),
  error = function(e) NULL)

if (!is.null(ic.param)) {
  cat("\n--- COMPARACAO: Hamilton parametrico (cum1) vs completo (cum1+cum2) ---\n")
  cat("                  | Param only | Completo   | Razao\n")
  cat("------------------|------------|------------|------\n")
  cat(sprintf("SE r* medio       | %8.2f   | %8.2f   | %.1fx\n",
      mean(ic.param$se.rstar), mean(se.rstar),
      mean(se.rstar) / mean(ic.param$se.rstar)))
  cat(sprintf("SE r* atual       | %8.2f   | %8.2f   | %.1fx\n",
      tail(ic.param$se.rstar, 1), tail(se.rstar, 1),
      tail(se.rstar, 1) / tail(ic.param$se.rstar, 1)))
  cat(sprintf("SE g medio        | %8.2f   | %8.2f   | %.1fx\n",
      mean(ic.param$se.g), mean(se.g),
      mean(se.g) / mean(ic.param$se.g)))
  cat(sprintf("SE z medio        | %8.2f   | %8.2f   | %.1fx\n",
      mean(ic.param$se.z), mean(se.z),
      mean(se.z) / mean(ic.param$se.z)))
}

# ==============================================================================
# GRAFICO DASHBOARD (4 paineis, como antes)
# ==============================================================================

data.arquivo <- format(Sys.Date(), "%Y%m%d")

png(file.path(wgm_path, paste0("outputs/graficos/fig_rstar_ic_hamilton_completo_", data.arquivo, ".png")),
    width = 1400, height = 900, res = 130)
par(mfrow = c(2, 2), mar = c(4, 4.5, 3, 1))

# 1. r* com IC completo
plot(datas, rstar.base.smooth, type = "l", col = "red", lwd = 2.5,
     main = expression(paste(r^"*", " - Hamilton completo (cum1 + cum2)")),
     ylab = "% a.a.", xlab = "",
     ylim = range(c(rstar.lower, rstar.upper), na.rm = TRUE))
polygon(c(datas, rev(datas)), c(rstar.upper, rev(rstar.lower)),
        col = rgb(1, 0, 0, 0.15), border = NA)
lines(datas, rstar.base.smooth, col = "red", lwd = 2.5)
lines(datas, taxa.real, col = "gray50", lwd = 0.8)
rect(as.Date("2020-04-01"), par("usr")[3],
     as.Date("2022-12-01"), par("usr")[4],
     col = rgb(0.5, 0.5, 0.5, 0.08), border = NA)
legend("topright", c(expression(r^"*"), "IC 95%", "Taxa real"),
       col = c("red", rgb(1,0,0,0.3), "gray50"), lwd = c(2.5, 8, 0.8), bty = "n")

# 2. Decomposicao SE: cum1 vs cum2
plot(datas, se.rstar.param, type = "l", col = "blue", lwd = 2,
     main = "Decomposicao SE(r*): param vs estado",
     ylab = "SE (pp)", xlab = "",
     ylim = c(0, max(se.rstar) * 1.1))
lines(datas, se.rstar.estado, col = "orange", lwd = 2)
lines(datas, se.rstar, col = "red", lwd = 2, lty = 2)
legend("topleft",
       c("Parametro (cum1)", "Estado (cum2)", "Total"),
       col = c("blue", "orange", "red"), lwd = 2, lty = c(1,1,2), bty = "n")

# 3. g com IC
plot(datas, g.pct, type = "l", col = "blue", lwd = 2,
     main = "g (crescimento potencial) com IC 95%",
     ylab = "% a.a.", xlab = "",
     ylim = range(c(g.lower, g.upper)))
polygon(c(datas, rev(datas)), c(g.upper, rev(g.lower)),
        col = rgb(0, 0, 1, 0.12), border = NA)
lines(datas, g.pct, col = "blue", lwd = 2)

# 4. z com IC
plot(datas, z.base.smooth, type = "l", col = "darkgreen", lwd = 2,
     main = "z (outros fatores) com IC 95%",
     ylab = "% a.a.", xlab = "",
     ylim = range(c(z.lower, z.upper)))
polygon(c(datas, rev(datas)), c(z.upper, rev(z.lower)),
        col = rgb(0, 0.5, 0, 0.12), border = NA)
lines(datas, z.base.smooth, col = "darkgreen", lwd = 2)

dev.off()
cat(sprintf("\nDashboard: fig_rstar_ic_hamilton_completo_%s.png\n", data.arquivo))

# ==============================================================================
# GRAFICO (A): r* suavizado sem bandas (fig_rstar_base.png)
# ==============================================================================

png(file.path(wgm_path, "outputs/graficos/fig_rstar_base.png"),
    width = 1200, height = 650, res = 130)
par(mar = c(4, 4.5, 3, 1))

plot(datas, rstar.base.smooth, type = "l", col = "red", lwd = 2.5,
     main = expression(paste("Taxa Natural de Juros (", r^"*", ") - Estimativa Suavizada")),
     ylab = "% a.a.", xlab = "",
     ylim = range(c(rstar.base.smooth, taxa.real), na.rm = TRUE) * c(0.9, 1.1))
lines(datas, taxa.real, col = "gray50", lwd = 1)
rect(as.Date("2020-04-01"), par("usr")[3],
     as.Date("2022-12-01"), par("usr")[4],
     col = rgb(0.5, 0.5, 0.5, 0.08), border = NA)
abline(h = mean(rstar.base.smooth), col = "red", lty = 3, lwd = 0.8)
legend("topright",
       c(expression(r^"*" ~ suavizado), "Taxa real observada",
         paste0("Media r* = ", round(mean(rstar.base.smooth), 1), "%")),
       col = c("red", "gray50", "red"),
       lwd = c(2.5, 1, 0.8), lty = c(1, 1, 3), bty = "n", cex = 0.85)

dev.off()
cat("Grafico (A): fig_rstar_base.png\n")

# ==============================================================================
# GRAFICO (B): g com banda total (fig_g_ic_total.png)
# ==============================================================================

png(file.path(wgm_path, "outputs/graficos/fig_g_ic_total.png"),
    width = 1200, height = 650, res = 130)
par(mar = c(4, 4.5, 3, 1))

plot(datas, g.pct, type = "l", col = "blue", lwd = 2.5,
     main = "Crescimento Potencial (g) com IC 95% - Hamilton Completo",
     ylab = "% a.a.", xlab = "",
     ylim = range(c(g.lower, g.upper)))
polygon(c(datas, rev(datas)), c(g.upper, rev(g.lower)),
        col = rgb(0, 0, 1, 0.12), border = NA)
lines(datas, g.pct, col = "blue", lwd = 2.5)
lines(datas, g.upper, col = "blue", lty = 2, lwd = 0.7)
lines(datas, g.lower, col = "blue", lty = 2, lwd = 0.7)
abline(h = 0, col = "gray", lty = 3)
rect(as.Date("2020-04-01"), par("usr")[3],
     as.Date("2022-12-01"), par("usr")[4],
     col = rgb(0.5, 0.5, 0.5, 0.08), border = NA)
legend("topright",
       c("g (crescimento potencial)", "IC 95%"),
       col = c("blue", rgb(0, 0, 1, 0.3)),
       lwd = c(2.5, 8), bty = "n", cex = 0.85)

dev.off()
cat("Grafico (B): fig_g_ic_total.png\n")

# ==============================================================================
# GRAFICO (C): decomposicao da incerteza do r* (fig_decomp_se_rstar.png)
# ==============================================================================

png(file.path(wgm_path, "outputs/graficos/fig_decomp_se_rstar.png"),
    width = 1200, height = 650, res = 130)
par(mar = c(4, 4.5, 3, 1))

plot(datas, se.rstar, type = "l", col = "red", lwd = 2.5,
     main = expression(paste("Decomposicao do Erro Padrao de ", r^"*")),
     ylab = "SE (pp)", xlab = "",
     ylim = c(0, max(se.rstar) * 1.1))
lines(datas, se.rstar.param, col = "blue", lwd = 2)
lines(datas, se.rstar.estado, col = "orange", lwd = 2)
abline(h = mean(se.rstar), col = "red", lty = 3, lwd = 0.8)
abline(h = mean(se.rstar.param), col = "blue", lty = 3, lwd = 0.8)
rect(as.Date("2020-04-01"), par("usr")[3],
     as.Date("2022-12-01"), par("usr")[4],
     col = rgb(0.5, 0.5, 0.5, 0.08), border = NA)
legend("topleft",
       c(paste0("Total (media: ", round(mean(se.rstar), 1), " pp)"),
         paste0("Parametro - cum1 (media: ", round(mean(se.rstar.param), 1), " pp)"),
         paste0("Estado - cum2 (media: ", round(mean(se.rstar.estado), 1), " pp)")),
       col = c("red", "blue", "orange"),
       lwd = c(2.5, 2, 2), bty = "n", cex = 0.85)

dev.off()
cat("Grafico (C): fig_decomp_se_rstar.png\n")

# ==============================================================================
# SALVAR RESULTADOS (.rds)
# ==============================================================================

ic.completo <- list(
  # Series base (suavizadas)
  rstar.base = rstar.base.smooth,
  g.base = g.base.smooth,        # fracao (ex: 0.021)
  g.base.pct = g.pct,            # percentual (ex: 2.1)
  z.base = z.base.smooth,
  ypot.base = ypot.base.smooth,

  # cum1 por periodo [T x 4: y*, r*, g, z]
  cum1 = cum1,
  # cum2 por periodo [T x 4: y*, r*, g, z]
  cum2 = cum2,

  # SE parametrico = sqrt(cum1) por variavel
  se.rstar.param = se.rstar.param,
  se.g.param = se.g.param,
  se.z.param = se.z.param,

  # SE estado = sqrt(cum2) por variavel
  se.rstar.estado = se.rstar.estado,
  se.g.estado = se.g.estado,
  se.z.estado = se.z.estado,

  # SE total = sqrt(cum1 + cum2) por variavel
  se.rstar = se.rstar,
  se.g = se.g,
  se.z = se.z,
  se.ypot = se.ypot,

  # Matrizes completas de SE [T x 4]
  se.cum1.mat = se.cum1,
  se.cum2.mat = se.cum2,
  se.total.mat = se.total,

  # Bandas 95% (total)
  rstar.upper = rstar.upper,
  rstar.lower = rstar.lower,
  g.upper = g.upper,
  g.lower = g.lower,
  z.upper = z.upper,
  z.lower = z.lower,

  # Metadados
  datas = datas,
  anos.est = anos.est,
  tri.est = tri.est,
  taxa.real = taxa.real,
  gap = stage3$output.gap.smoothed,

  # Info matrix e t-stats
  bse = bse,
  info = info,
  params.free = params.free,

  # MC config
  niter = niter,
  n.accepted = good.draws,
  n.excluded = excluded.draws,
  n.excluded.constraint = excluded.constraint,
  n.excluded.error = excluded.error,
  elapsed.secs = elapsed,

  metodo = "Hamilton (1986) completo: cum1 (param) + cum2 (estado via RTS smoother P.tT)"
)

rds.path <- file.path(wgm_path,
    paste0("outputs/results/ic_hamilton_completo_selic_", data.arquivo, ".rds"))
saveRDS(ic.completo, rds.path)
cat(sprintf("\nResultados salvos: %s\n", basename(rds.path)))

# ==============================================================================
# LISTA DE ARQUIVOS SALVOS
# ==============================================================================

cat("\n================================================================\n")
cat("ARQUIVOS SALVOS\n")
cat("================================================================\n")
cat(sprintf("  .rds:  outputs/results/ic_hamilton_completo_selic_%s.rds\n", data.arquivo))
cat(sprintf("  (A)    outputs/graficos/fig_rstar_base.png\n"))
cat(sprintf("  (B)    outputs/graficos/fig_g_ic_total.png\n"))
cat(sprintf("  (C)    outputs/graficos/fig_decomp_se_rstar.png\n"))
cat(sprintf("  (D)    outputs/graficos/fig_rstar_ic_hamilton_completo_%s.png\n", data.arquivo))

cat("\n================================================================\n")
cat("HAMILTON COMPLETO FINALIZADO\n")
cat("================================================================\n")
cat("\nExemplo de uso do .rds:\n")
cat("  ic <- readRDS('outputs/results/ic_hamilton_completo_selic_YYYYMMDD.rds')\n")
cat("  ic$rstar.base       # serie suavizada r*\n")
cat("  ic$se.rstar         # SE total por periodo\n")
cat("  ic$se.rstar.param   # SE parametrico (sqrt cum1)\n")
cat("  ic$se.rstar.estado  # SE estado (sqrt cum2)\n")
cat("  ic$cum1             # cum1 bruto [T x 4]\n")
cat("  ic$cum2             # cum2 bruto [T x 4]\n")
cat("  ic$rstar.upper      # banda superior 95%\n")
cat("  ic$rstar.lower      # banda inferior 95%\n")
