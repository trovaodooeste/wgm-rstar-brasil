# ==============================================================================
# ICs PARA MODELO DI-PRE (mesmo procedimento do Selic)
# ==============================================================================

library(mFilter)

wgm_path <- "."

cat("=== IC HAMILTON PARAMETRICO - MODELO DI-PRE ===\n\n")

# Carregar resultado DI-Pre
resultado <- readRDS(file.path(wgm_path, "outputs/results/wgm_dipre_params_fixos_20260212.rds"))

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

cat(sprintf("T = %d trimestres\n", T.obs))

# Reconstruir xi.00, P.00, bounds
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

cat(sprintf("Parametros livres: %d\n", n.free))

# Kalman filter (identico ao Selic)
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

# Matriz de informacao
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

# Monte Carlo
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

cum1 <- matrix(0, T.obs, 4)
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

# Resultados
cat(sprintf("\n  r* atual: %.2f%% [IC: %.2f%% - %.2f%%]\n",
    tail(rstar.base, 1),
    tail(rstar.base, 1) - 1.96*tail(se.rstar, 1),
    tail(rstar.base, 1) + 1.96*tail(se.rstar, 1)))
cat(sprintf("  r* medio: %.2f%% [IC: %.2f%% - %.2f%%]\n",
    mean(rstar.base),
    mean(rstar.base) - 1.96*mean(se.rstar),
    mean(rstar.base) + 1.96*mean(se.rstar)))
cat(sprintf("  SE medio r*: %.2f pp\n", mean(se.rstar)))

# Salvar
taxa.real <- as.numeric(resultado$dados$taxa_real_exante_tri[
  match(paste(anos.est, tri.est), paste(resultado$dados$ano, resultado$dados$trimestre))])

data.arquivo <- format(Sys.Date(), "%Y%m%d")

ic.results <- list(
  se.rstar = se.rstar, se.g = se.g, se.z = se.z, se.ypot = se[, 1],
  rstar.base = rstar.base, g.base = g.base, z.base = z.base,
  rstar.upper = rstar.base + 1.96 * se.rstar,
  rstar.lower = rstar.base - 1.96 * se.rstar,
  g.upper = g.base * 100 + 1.96 * se.g,
  g.lower = g.base * 100 - 1.96 * se.g,
  z.upper = z.base + 1.96 * se.z,
  z.lower = z.base - 1.96 * se.z,
  datas = datas, anos.est = anos.est, tri.est = tri.est,
  taxa.real = taxa.real,
  gap = stage3$output.gap.smoothed,
  bse = bse, params.free = params.free, info = info,
  niter = niter, n.excluded = excluded.draws, elapsed.secs = elapsed
)

saveRDS(ic.results, file.path(wgm_path, paste0("outputs/results/ic_hamilton_parametrico_dipre_", data.arquivo, ".rds")))
cat(sprintf("\nSalvo: outputs/results/ic_hamilton_parametrico_dipre_%s.rds\n", data.arquivo))
