##------------------------------------------------------------------------------##
## File:        rstar_stage3_com_kappa.R
##
## Description: Stage 3 WGM (Wellington Gubiani Maiberg)
##              Baseado em HLW 2023, adaptado para Brasil
##              COM KAPPA (volatilidade variavel no tempo)
##              Decomposicao: r* = c*g + z
##
## Curva IS:
##   y_t = y*_t + a_1*(y_{t-1} - y*_{t-1}) + a_2*(y_{t-2} - y*_{t-2})
##       + a_r/2 * [(r_{t-1} - g_{t-1} - z_{t-1}) + (r_{t-2} - g_{t-2} - z_{t-2})]
##
## Curva de Phillips:
##   pi_t = b_pi1*pi_{t-1} + b_pi2*pi_bar_{t-2:4} + (1-b_pi1-b_pi2)*pi_bar_{t-5:8}
##        + b_y * gap_{t-1}
##
## Estados:
##   y*_t = y*_{t-1} + g_{t-1} + e_y*
##   g_t = g_{t-1} + e_g
##   z_t = z_{t-1} + e_z
##
## r* = c*g + z (com c=1)
##
## KAPPA: volatilidade variavel durante COVID (2020 Q2-Q4, 2021, 2022)
##        kappa^2 aplicado a matriz R (erro de medicao)
##
## Data: 2026-02-08 (Adicionado KAPPA)
##------------------------------------------------------------------------------##

# Pacotes
if (!require("nloptr")) install.packages("nloptr")
if (!require("numDeriv")) install.packages("numDeriv")
if (!require("mFilter")) install.packages("mFilter")

library(nloptr)
library(numDeriv)
library(mFilter)

##------------------------------------------------------------------------------##
## Funcao principal: rstar.stage3.com.kappa
##------------------------------------------------------------------------------##
rstar.stage3.com.kappa <- function(log.output,
                                    inflation,
                                    real.interest.rate,
                                    stage2.results,
                                    anos,                     # Vetor de anos
                                    trimestres,               # Vetor de trimestres (1-4)
                                    c.param = 1,              # c na eq r* = c*g + z
                                    lambda.g = NA,            # de Stage 2
                                    lambda.z = NA,            # se NA, usa 0.05 (apenas p/ valor inicial)
                                    z.initial = 2.2,          # Moreira & Portugal: z0 = 2.2
                                    sigma.ystar.constraint = 0.5,
                                    sigma.ystar.lower = 0.001,  # Floor para sigma_y*
                                    sigma.z.upper = 2.2,      # Teto sigma_z (MLE livre com bound)
                                    xi.00 = NA,
                                    P.00 = NA) {

  cat("\n=================================================================\n")
  cat("STAGE 3 WGM - Decomposicao r* = c*g + z\n")
  cat("=================================================================\n\n")

  ##--------------------------------------------------------------------------##
  ## 1. Extrair parametros FIXOS dos Stages 1 e 2
  ##--------------------------------------------------------------------------##

  cat("1. Fixando parametros dos Stages 1 e 2...\n")

  # Do Stage 2:
  a.1.fixed    <- stage2.results$theta[1]
  a.2.fixed    <- stage2.results$theta[2]
  a.r.fixed    <- stage2.results$theta[3]
  b.pi1.fixed  <- stage2.results$theta[4]
  b.pi2.fixed  <- stage2.results$theta[5]
  b.y.fixed    <- stage2.results$theta[6]

  cat(sprintf("   a_1 = %.4f (fixo)\n", a.1.fixed))
  cat(sprintf("   a_2 = %.4f (fixo)\n", a.2.fixed))
  cat(sprintf("   a_r = %.6f (fixo)\n", a.r.fixed))
  cat(sprintf("   b_pi1 = %.4f (fixo)\n", b.pi1.fixed))
  cat(sprintf("   b_pi2 = %.4f (fixo)\n", b.pi2.fixed))
  cat(sprintf("   b_y = %.4f (fixo)\n", b.y.fixed))
  cat(sprintf("   c = %.1f (fixo)\n", c.param))

  # Extrair kappas do Stage 2 (se disponíveis)
  estimar.kappa <- stage2.results$estimar.kappa
  if (is.null(estimar.kappa)) estimar.kappa <- FALSE

  # Kappas APENAS para COVID (2008 usa dummy, não kappa)
  if (estimar.kappa && !is.null(stage2.results$kappas)) {
    kappas <- stage2.results$kappas
    # Verificar se os kappas são válidos (não NA)
    if (!any(is.na(kappas))) {
      kappa_2020 <- as.numeric(kappas["kappa_2020"])
      kappa_2021 <- as.numeric(kappas["kappa_2021"])
      kappa_2022 <- as.numeric(kappas["kappa_2022"])
      cat("\n   Kappas do Stage 2 (fixos, apenas COVID):\n")
      cat(sprintf("   kappa_2020 = %.2f\n", kappa_2020))
      cat(sprintf("   kappa_2021 = %.2f\n", kappa_2021))
      cat(sprintf("   kappa_2022 = %.2f\n", kappa_2022))
      cat("   (2008 usa dummy, não kappa)\n")
    } else {
      kappa_2020 <- 1.0
      kappa_2021 <- 1.0
      kappa_2022 <- 1.0
      estimar.kappa <- FALSE
      cat("\n   Kappas = 1.0 (kappas do Stage 2 são NA)\n")
    }
  } else {
    kappa_2020 <- 1.0
    kappa_2021 <- 1.0
    kappa_2022 <- 1.0
    cat("\n   Kappas = 1.0 (sem volatilidade variável)\n")
  }

  # Lambda_g do Stage 2
  if (is.na(lambda.g)) {
    lambda.g <- stage2.results$lambda.g
  }
  cat(sprintf("   lambda_g = %.4f (do Stage 2)\n", lambda.g))

  # Lambda_z
  if (is.na(lambda.z)) {
    lambda.z <- 0.05  # Valor tipico HLW
  }
  cat(sprintf("   lambda_z = %.4f\n", lambda.z))

  ##--------------------------------------------------------------------------##
  ## 2. Preparar dados
  ##--------------------------------------------------------------------------##

  T <- length(log.output) - 8

  cat(sprintf("\n2. Preparando dados (T = %d trimestres)...\n", T))

  # y.data: [100*log(y), pi]
  y.data <- cbind(
    100 * log.output[9:(T+8)],
    inflation[9:(T+8)]
  )

  # x.data
  x.data <- cbind(
    100 * log.output[8:(T+7)],   # 1: y_{t-1}
    100 * log.output[7:(T+6)],   # 2: y_{t-2}
    real.interest.rate[8:(T+7)], # 3: r_{t-1}
    real.interest.rate[7:(T+6)], # 4: r_{t-2}
    inflation[8:(T+7)],          # 5: pi_{t-1}
    (inflation[7:(T+6)] + inflation[6:(T+5)] + inflation[5:(T+4)]) / 3,   # 6: pi_bar_{t-2:t-4}
    (inflation[4:(T+3)] + inflation[3:(T+2)] + inflation[2:(T+1)] + inflation[1:T]) / 4   # 7: pi_bar_{t-5:t-8}
  )

  # Criar vetor de kappa para cada observacao (T trimestres)
  anos.est <- anos[9:(T+8)]
  trimestres.est <- trimestres[9:(T+8)]

  # Kappa apenas para COVID (2020 Q2-Q4, 2021, 2022)
  # 2008 usa dummy (phi_2008), não kappa
  kappa.vec <- rep(1.0, T)
  for (t.i in 1:T) {
    ano <- anos.est[t.i]
    tri <- trimestres.est[t.i]
    if (ano == 2020 && tri >= 2) {
      kappa.vec[t.i] <- kappa_2020
    } else if (ano == 2021) {
      kappa.vec[t.i] <- kappa_2021
    } else if (ano == 2022) {
      kappa.vec[t.i] <- kappa_2022
    }
  }

  n.crises <- sum(kappa.vec > 1.0)
  cat(sprintf("   Trimestres com kappa > 1: %d\n", n.crises))

  ##--------------------------------------------------------------------------##
  ## 3. Parametros Stage 3
  ##--------------------------------------------------------------------------##

  cat("\n3. Configurando parametros para Stage 3...\n")

  # Stage 3 param vector (12 params):
  # [a_1, a_2, a_r, b_pi1, b_pi2, b_y, sigma_is, sigma_pi, sigma_ystar, sigma_g, sigma_z, c]

  param.num <- c("a_1"=1, "a_2"=2, "a_r"=3, "b_pi1"=4, "b_pi2"=5,
                 "b_y"=6, "sigma_is"=7, "sigma_pi"=8,
                 "sigma_ystar"=9, "sigma_g"=10, "sigma_z"=11, "c"=12)

  # Valores iniciais do Stage 2
  sigma.is.s2 <- stage2.results$theta[8]
  sigma.pi.s2 <- stage2.results$theta[9]
  sigma.ystar.s2 <- stage2.results$theta[10]

  # sigma_g e sigma_z via lambdas
  sigma.g.initial <- lambda.g * sigma.ystar.s2
  sigma.z.initial <- lambda.z * sigma.is.s2 / abs(a.r.fixed)  # lambda_z = a_r * sigma_z / sigma_is

  initial.parameters <- c(
    a.1.fixed,        # 1: a_1 (FIXO)
    a.2.fixed,        # 2: a_2 (FIXO)
    a.r.fixed,        # 3: a_r (FIXO)
    b.pi1.fixed,      # 4: b_pi1 (FIXO)
    b.pi2.fixed,      # 5: b_pi2 (FIXO)
    b.y.fixed,        # 6: b_y (FIXO)
    sigma.is.s2,      # 7: sigma_IS (estimar)
    sigma.pi.s2,      # 8: sigma_pi (estimar)
    sigma.ystar.s2,   # 9: sigma_y* (estimar)
    sigma.g.initial,  # 10: sigma_g (via lambda_g)
    sigma.z.initial,  # 11: sigma_z (via lambda_z)
    c.param           # 12: c (FIXO = 1)
  )

  cat("\n   Parametros iniciais:\n")
  for (i in 1:length(initial.parameters)) {
    status <- if (i %in% c(1,2,3,4,5,6,12)) "FIXO" else "estimar"
    cat(sprintf("   [%2d] %-12s = %10.6f (%s)\n",
                i, names(param.num)[i], initial.parameters[i], status))
  }

  ##--------------------------------------------------------------------------##
  ## 4. Bounds
  ##--------------------------------------------------------------------------##

  theta.lb <- rep(-Inf, 12)
  theta.ub <- rep(Inf, 12)

  # FIXAR parametros estruturais
  theta.lb[1] <- a.1.fixed; theta.ub[1] <- a.1.fixed
  theta.lb[2] <- a.2.fixed; theta.ub[2] <- a.2.fixed
  theta.lb[3] <- a.r.fixed; theta.ub[3] <- a.r.fixed
  theta.lb[4] <- b.pi1.fixed; theta.ub[4] <- b.pi1.fixed
  theta.lb[5] <- b.pi2.fixed; theta.ub[5] <- b.pi2.fixed
  theta.lb[6] <- b.y.fixed; theta.ub[6] <- b.y.fixed
  theta.lb[12] <- c.param; theta.ub[12] <- c.param

  # Variancias positivas
  theta.lb[7] <- 0.001               # sigma_IS
  theta.lb[8] <- 0.001               # sigma_pi
  theta.lb[9] <- sigma.ystar.lower   # sigma_y* (configurable floor)
  theta.lb[10] <- 0.001              # sigma_g
  theta.lb[11] <- 0.001              # sigma_z

  # Ajustar valor inicial de sigma_y* se abaixo do floor
  if (initial.parameters[9] < sigma.ystar.lower) {
    initial.parameters[9] <- sigma.ystar.lower * 1.5
    cat(sprintf("   Ajustando sigma_y* inicial para %.4f (acima do floor %.4f)\n",
                initial.parameters[9], sigma.ystar.lower))
  }

  # Restricao sigma_y*
  if (!is.na(sigma.ystar.constraint)) {
    theta.ub[9] <- sigma.ystar.constraint
    initial.parameters[9] <- min(initial.parameters[9], sigma.ystar.constraint * 0.9)
    cat(sprintf("   Restricao: %.4f <= sigma_y* <= %.4f\n", sigma.ystar.lower, sigma.ystar.constraint))
  }

  # ===========================================================================
  # SIGMA_Z: LIVRE por MLE com teto (upper bound)
  # MLE decide o valor otimo; teto impede explosao
  # Justificativa: lambda_z fixo matava variacao de z (SD=0.03)
  # ===========================================================================
  theta.ub[11] <- sigma.z.upper
  if (initial.parameters[11] > sigma.z.upper * 0.9) {
    initial.parameters[11] <- sigma.z.upper * 0.9
  }
  cat(sprintf("   sigma_z: LIVRE por MLE, teto = %.2f (valor inicial = %.4f)\n",
              sigma.z.upper, initial.parameters[11]))

  ##--------------------------------------------------------------------------##
  ## 5. Estados iniciais
  ##--------------------------------------------------------------------------##

  cat(sprintf("\n4. Inicializando estados (z_0 = %.1f, Moreira & Portugal)...\n", z.initial))

  if (any(is.na(xi.00))) {
    # HP filter para tendencia inicial
    y.og <- log.output[5:(T+8)]
    g.pot <- hpfilter(y.og, freq = 36000, type = "lambda", drift = FALSE)$trend
    g.pot.diff <- diff(g.pot)

    # Stage 3 states: [y*_t, y*_{t-1}, gap_t, gap_{t-1}, g_t, g_{t-1}, z_t, z_{t-1}]
    xi.00 <- c(
      100 * g.pot[4],           # y*_t
      100 * g.pot[3],           # y*_{t-1}
      0,                        # gap_t
      0,                        # gap_{t-1}
      100 * mean(g.pot.diff),   # g_t
      100 * mean(g.pot.diff),   # g_{t-1}
      z.initial,                # z_t (Moreira & Portugal: 2.2)
      z.initial                 # z_{t-1}
    )
    cat(sprintf("   xi.00: y*=%.1f, g=%.2f, z=%.1f\n", xi.00[1], xi.00[5], xi.00[7]))
  }

  if (any(is.na(P.00))) {
    P.00 <- diag(1, 8)
  }

  ##--------------------------------------------------------------------------##
  ## 6. Kalman Filter Stage 3
  ##--------------------------------------------------------------------------##

  kalman.filter.stage3 <- function(theta, y.data, x.data, xi.00, P.00, param.num, lambda.g, lambda.z, kappa.vec) {

    a.1 <- theta[param.num["a_1"]]
    a.2 <- theta[param.num["a_2"]]
    a.r <- theta[param.num["a_r"]]
    b.pi1 <- theta[param.num["b_pi1"]]
    b.pi2 <- theta[param.num["b_pi2"]]
    b.y <- theta[param.num["b_y"]]
    sigma.is <- theta[param.num["sigma_is"]]
    sigma.pi <- theta[param.num["sigma_pi"]]
    sigma.ystar <- theta[param.num["sigma_ystar"]]
    c.val <- theta[param.num["c"]]

    # sigma_g constrainado via lambda_g; sigma_z LIVRE (do theta, MLE decide)
    sigma.g <- lambda.g * sigma.ystar
    sigma.z <- theta[param.num["sigma_z"]]  # LIVRE por MLE (com teto)

    # Matrizes do sistema
    # Estado: [y*_t, y*_{t-1}, gap_t, gap_{t-1}, g_t, g_{t-1}, z_t, z_{t-1}]

    n.states <- 8

    # F: transicao
    F <- matrix(0, n.states, n.states)
    F[1,1] <- 1; F[1,5] <- 1   # y*_t = y*_{t-1} + g_{t-1}
    F[2,1] <- 1                 # y*_{t-1} lag
    F[3,3] <- 0                 # gap_t (calculado)
    F[4,3] <- 1                 # gap_{t-1} lag
    F[5,5] <- 1                 # g_t = g_{t-1}
    F[6,5] <- 1                 # g_{t-1} lag
    F[7,7] <- 1                 # z_t = z_{t-1}
    F[8,7] <- 1                 # z_{t-1} lag

    # Q: covariancia dos choques (CONSTANTE - não varia com kappa)
    Q <- matrix(0, n.states, n.states)
    Q[1,1] <- sigma.ystar^2
    Q[5,5] <- sigma.g^2
    Q[7,7] <- sigma.z^2

    # R base: erros de observacao (será multiplicado por kappa^2)
    R.base <- diag(c(sigma.is^2, sigma.pi^2))

    T.obs <- nrow(y.data)
    xi.tt <- xi.00
    P.tt <- P.00

    ll <- 0
    xi.filtered <- matrix(0, T.obs, n.states)
    P.tt.array <- array(0, dim = c(n.states, n.states, T.obs))  # Covariancia filtrada por periodo

    for (t in 1:T.obs) {
      # Aplicar kappa^2 ao R (conforme HLW 2023)
      kappa.t <- kappa.vec[t]
      R <- R.base * (kappa.t^2)

      # Predict
      xi.ttm1 <- F %*% xi.tt
      P.ttm1 <- F %*% P.tt %*% t(F) + Q

      # Observacao
      y.t <- y.data[t, ]

      # r* = c*g + z
      rstar.t <- c.val * xi.ttm1[5] / 100 + xi.ttm1[7]  # g em % tri, z em % aa
      rstar.tm1 <- c.val * xi.ttm1[6] / 100 + xi.ttm1[8]

      # Hiato dos juros: (r - r*)
      r.gap <- (x.data[t, 3] + x.data[t, 4]) / 2 - (rstar.t + rstar.tm1) / 2

      # IS: y_t = y*_t + a_1*(y_{t-1} - y*_{t-1}) + a_2*(y_{t-2} - y*_{t-2}) + a_r*r_gap
      y.pred.is <- xi.ttm1[1] + a.1 * (x.data[t, 1] - xi.ttm1[2]) +
                   a.2 * (x.data[t, 2] - xi.ttm1[2]) + a.r * r.gap

      # Phillips
      pi.pred <- b.pi1 * x.data[t, 5] + b.pi2 * x.data[t, 6] +
                 (1 - b.pi1 - b.pi2) * x.data[t, 7] + b.y * xi.ttm1[3]

      y.pred <- c(y.pred.is, pi.pred)
      v.t <- y.t - y.pred

      # H matrix: ∂(observacao)/∂(estado)
      # Estados: [y*_t, y*_{t-1}, gap_t, gap_{t-1}, g_t, g_{t-1}, z_t, z_{t-1}]
      # Obs: [y_t (IS), pi_t (Phillips)]
      H <- matrix(0, n.states, 2)

      # IS curve: y_t = y*_t + a_1*(y_{t-1} - y*_{t-1}) + a_2*(y_{t-2} - y*_{t-2}) + a_r*r_gap
      # r_gap = r_avg - (r*_t + r*_{t-1})/2
      # r*_t = c*g_t/100 + z_t
      H[1, 1] <- 1              # ∂y/∂y*_t = 1
      H[2, 1] <- -(a.1 + a.2)   # ∂y/∂y*_{t-1} = -(a_1 + a_2) [y_{t-2} usa y*_{t-1} como proxy]
      H[5, 1] <- -a.r * c.val / 200  # ∂y/∂g_t = a_r * ∂r_gap/∂g_t = a_r * (-c/200)
      H[6, 1] <- -a.r * c.val / 200  # ∂y/∂g_{t-1}
      H[7, 1] <- -a.r / 2            # ∂y/∂z_t = a_r * ∂r_gap/∂z_t = a_r * (-1/2)
      H[8, 1] <- -a.r / 2            # ∂y/∂z_{t-1}

      # Phillips curve: pi_t = ... + b_y * gap_{t-1}
      H[3, 2] <- b.y            # ∂pi/∂gap_t (usando gap filtrado)

      # Kalman gain
      S.t <- t(H) %*% P.ttm1 %*% H + R
      S.t <- S.t + diag(1e-8, 2)

      K.t <- P.ttm1 %*% H %*% solve(S.t)

      # Update
      xi.tt <- xi.ttm1 + K.t %*% v.t
      P.tt <- P.ttm1 - K.t %*% t(H) %*% P.ttm1

      # Log-likelihood
      ll <- ll - 0.5 * (log(det(S.t)) + t(v.t) %*% solve(S.t) %*% v.t)

      xi.filtered[t, ] <- xi.tt
      P.tt.array[,,t] <- P.tt
    }

    return(list(ll = as.numeric(ll), xi.filtered = xi.filtered,
                P.tt = P.tt, xi.tt = xi.tt, P.tt.array = P.tt.array))
  }

  ##--------------------------------------------------------------------------##
  ## 7. Otimizacao MLE
  ##--------------------------------------------------------------------------##

  cat("\n5. Iniciando otimizacao MLE Stage 3...\n")

  f <- function(theta) {
    result <- tryCatch({
      kf <- kalman.filter.stage3(theta, y.data, x.data, xi.00, P.00, param.num, lambda.g, lambda.z, kappa.vec)
      -kf$ll
    }, error = function(e) {
      return(1e10)
    })
    return(result)
  }

  nloptr.out <- nloptr(
    x0 = initial.parameters,
    eval_f = f,
    eval_grad_f = function(x) { grad(f, x) },
    lb = theta.lb,
    ub = theta.ub,
    opts = list(
      algorithm = "NLOPT_LD_LBFGS",
      xtol_rel = 1e-8,
      maxeval = 5000
    )
  )

  theta <- nloptr.out$solution

  cat(sprintf("\n   Status: %s\n", nloptr.out$message))

  ##--------------------------------------------------------------------------##
  ## 8. Extrair resultados
  ##--------------------------------------------------------------------------##

  kf.final <- kalman.filter.stage3(theta, y.data, x.data, xi.00, P.00, param.num, lambda.g, lambda.z, kappa.vec)

  # Estados filtrados
  # [y*_t, y*_{t-1}, gap_t, gap_{t-1}, g_t, g_{t-1}, z_t, z_{t-1}]

  # PIB potencial
  potential.smoothed <- kf.final$xi.filtered[, 1] / 100

  # Output gap
  output.gap.smoothed <- y.data[, 1] - kf.final$xi.filtered[, 1]

  # g (tendencia) - converter para % anual
  g.smoothed <- kf.final$xi.filtered[, 5] * 4 / 100  # tri -> anual

  # z (outros fatores) - ja em % anual
  z.smoothed <- kf.final$xi.filtered[, 7]

  # r* = c*g + z
  rstar.smoothed <- c.param * g.smoothed * 100 + z.smoothed  # Em % anual

  ##--------------------------------------------------------------------------##
  ## 9. Resultados
  ##--------------------------------------------------------------------------##

  cat("\n=================================================================\n")
  cat("STAGE 3 WGM CONCLUIDO!\n")
  cat("=================================================================\n\n")

  cat("Parametros estimados:\n")
  cat(sprintf("  a_1 = %.4f (fixo)\n", theta[1]))
  cat(sprintf("  a_2 = %.4f (fixo)\n", theta[2]))
  cat(sprintf("  a_r = %.6f (fixo)\n", theta[3]))
  cat(sprintf("  b_pi1 = %.4f (fixo)\n", theta[4]))
  cat(sprintf("  b_pi2 = %.4f (fixo)\n", theta[5]))
  cat(sprintf("  b_y = %.4f (fixo)\n", theta[6]))
  cat(sprintf("  sigma_IS = %.4f\n", theta[7]))
  cat(sprintf("  sigma_pi = %.4f\n", theta[8]))
  cat(sprintf("  sigma_y* = %.4f\n", theta[9]))
  cat(sprintf("  sigma_g = %.4f (via lambda_g)\n", lambda.g * theta[9]))
  cat(sprintf("  sigma_z = %.4f (livre MLE, teto=%.2f)\n", theta[11], sigma.z.upper))
  cat(sprintf("  c = %.1f (fixo)\n", theta[12]))
  cat(sprintf("  lambda_g = %.4f\n", lambda.g))
  cat(sprintf("  lambda_z = %.4f\n", lambda.z))

  cat(sprintf("\nLog-likelihood = %.4f\n", kf.final$ll))

  cat(sprintf("\n=== DECOMPOSICAO r* = c*g + z ===\n"))
  cat(sprintf("  g (tendencia):    Media = %.2f%% aa, Ultimo = %.2f%% aa\n",
              mean(g.smoothed) * 100, tail(g.smoothed, 1) * 100))
  cat(sprintf("  z (outros):       Media = %.2f%% aa, Ultimo = %.2f%% aa\n",
              mean(z.smoothed), tail(z.smoothed, 1)))
  cat(sprintf("  r* = c*g + z:     Media = %.2f%% aa, Ultimo = %.2f%% aa\n",
              mean(rstar.smoothed), tail(rstar.smoothed, 1)))

  cat(sprintf("\nOutput Gap:\n"))
  cat(sprintf("  Media: %.2f%%\n", mean(output.gap.smoothed)))
  cat(sprintf("  Range: [%.2f%%, %.2f%%]\n",
              min(output.gap.smoothed), max(output.gap.smoothed)))

  return(list(
    theta = theta,
    param.num = param.num,
    log.likelihood = kf.final$ll,
    potential.smoothed = potential.smoothed,
    output.gap.smoothed = output.gap.smoothed,
    g.smoothed = g.smoothed,
    z.smoothed = z.smoothed,
    rstar.smoothed = rstar.smoothed,
    xi.filtered = kf.final$xi.filtered,
    P.tt.array = kf.final$P.tt.array,
    lambda.g = lambda.g,
    lambda.z = lambda.z,
    c.param = c.param,
    out.mle = nloptr.out,
    y.data = y.data,
    x.data = x.data,
    T = T,
    estimar.kappa = estimar.kappa,
    kappas = c(kappa_2020 = kappa_2020, kappa_2021 = kappa_2021, kappa_2022 = kappa_2022),
    kappa.vec = kappa.vec,
    anos.est = anos.est,
    trimestres.est = trimestres.est
  ))
}

cat("\n")
cat("=================================================================\n")
cat("rstar_stage3_com_kappa.R carregado\n")
cat("=================================================================\n")
cat("\n")
cat("METODOLOGIA HLW 2023 COM KAPPA:\n")
cat("  - Kappas fixos do Stage 2\n")
cat("  - kappa^2 aplicado a matriz R (erro de medicao)\n")
cat("  - Decomposicao: r* = c*g + z\n")
cat("\n")
cat("Uso:\n")
cat("  result <- rstar.stage3.com.kappa(log.output, inflation, real.interest.rate,\n")
cat("                                    stage2.results, anos, trimestres)\n")
cat("\n")
