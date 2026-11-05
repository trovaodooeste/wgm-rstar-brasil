##------------------------------------------------------------------------------##
## File:        rstar_stage2_com_kappa.R
##
## Description: Stage 2 WGM (Wellington Gubiani Maiberg)
##              Baseado em HLW 2023, adaptado para Brasil
##              COM KAPPA (volatilidade variavel no tempo)
##
## Curva IS (com taxa de juros):
##   ỹ_t = a_1·ỹ_{t-1} + a_2·ỹ_{t-2} + a_r/2·(r_{t-1} + r_{t-2} - 2·r*) + ε_IS
##
## Curva de Phillips (mantida do Stage 1):
##   π_t = b_π1·π_{t-1} + b_π2·π̄_{t-2:t-4} + (1-b_π1-b_π2)·π̄_{t-5:t-8} + b_y·ỹ_{t-1} + ε_π
##
## Estados:
##   y*_t = y*_{t-1} + g_{t-1}
##   g_t = g_{t-1} + ε_g
##
## r* = 4·g (anualizado, c=1 assumido)
##
## KAPPA: volatilidade variavel durante COVID (2020 Q2-Q4, 2021, 2022)
##        kappa^2 aplicado a matriz R (erro de medicao)
##
## Parâmetros Stage 2 (11 total):
##   a_1, a_2, a_r, b_π1, b_π2, b_y, g, σ_IS, σ_π, σ_y*, σ_g
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
## Função principal: rstar.stage2.com.kappa
##------------------------------------------------------------------------------##
rstar.stage2.com.kappa <- function(log.output,
                                    inflation,
                                    real.interest.rate,  # Taxa real (% aa)
                                    stage1.results,       # Resultados do Stage 1 (inclui kappas)
                                    anos,                 # Vetor de anos para cada observacao
                                    trimestres,           # Vetor de trimestres (1-4)
                                    b.y.constraint = 0.25,
                                    a.r.constraint = NA,  # Upper bound para a_r (se quiser forçar negativo)
                                    a.r.lower = NA,       # Lower bound para a_r (default: sem limite)
                                    sigma.ystar.constraint = 0.5,
                                    sigma.ystar.lower = 0.001,  # Floor para sigma_y*
                                    lambda.g = NA,        # sigma_g/sigma_y* ratio (se NA, estima)
                                    xi.00 = NA,
                                    P.00 = NA) {

  cat("\n=================================================================\n")
  cat("STAGE 2 WGM - Taxa Natural de Juros para Brasil\n")
  cat("=================================================================\n\n")

  ##--------------------------------------------------------------------------##
  ## 1. Extrair parâmetros FIXOS do Stage 1
  ##--------------------------------------------------------------------------##

  cat("1. Fixando parâmetros do Stage 1...\n")

  # Parâmetros do Stage 1 simples (9 params):
  # [a_1, a_2, b_π1, b_π2, b_y, g, σ_IS, σ_π, σ_y*]
  a.1.fixed    <- stage1.results$theta[1]
  a.2.fixed    <- stage1.results$theta[2]
  b.pi1.fixed  <- stage1.results$theta[3]
  b.pi2.fixed  <- stage1.results$theta[4]
  b.y.fixed    <- stage1.results$theta[5]
  g.initial    <- stage1.results$theta[6]
  sigma.is.s1  <- stage1.results$theta[7]
  sigma.pi.s1  <- stage1.results$theta[8]
  sigma.ystar.s1 <- stage1.results$theta[9]

  cat(sprintf("   a_1 = %.4f (fixo)\n", a.1.fixed))
  cat(sprintf("   a_2 = %.4f (fixo)\n", a.2.fixed))
  cat(sprintf("   b_π1 = %.4f (fixo)\n", b.pi1.fixed))
  cat(sprintf("   b_π2 = %.4f (fixo)\n", b.pi2.fixed))
  cat(sprintf("   b_y = %.4f (fixo)\n", b.y.fixed))
  cat(sprintf("   g_inicial = %.4f\n", g.initial))

  # Extrair kappas do Stage 1 (se disponíveis)
  estimar.kappa <- stage1.results$estimar.kappa
  if (is.null(estimar.kappa)) estimar.kappa <- FALSE

  # Kappas APENAS para COVID (2008 usa dummy, não kappa)
  if (estimar.kappa && !is.null(stage1.results$kappas)) {
    kappas <- stage1.results$kappas
    kappa_2020 <- as.numeric(kappas["kappa_2020"])
    kappa_2021 <- as.numeric(kappas["kappa_2021"])
    kappa_2022 <- as.numeric(kappas["kappa_2022"])
    cat("\n   Kappas do Stage 1 (fixos, apenas COVID):\n")
    cat(sprintf("   kappa_2020 = %.2f\n", kappa_2020))
    cat(sprintf("   kappa_2021 = %.2f\n", kappa_2021))
    cat(sprintf("   kappa_2022 = %.2f\n", kappa_2022))
    cat("   (2008 usa dummy, não kappa)\n")
  } else {
    kappa_2020 <- 1.0
    kappa_2021 <- 1.0
    kappa_2022 <- 1.0
    cat("\n   Kappas = 1.0 (sem volatilidade variável)\n")
  }

  ##--------------------------------------------------------------------------##
  ## 2. Preparar dados
  ##--------------------------------------------------------------------------##

  T <- length(log.output) - 8

  cat(sprintf("\n2. Preparando dados (T = %d trimestres)...\n", T))

  # y.data: [100×log(y), π]
  y.data <- cbind(
    100 * log.output[9:(T+8)],
    inflation[9:(T+8)]
  )

  # x.data para Stage 2 (inclui taxa de juros)
  x.data <- cbind(
    100 * log.output[8:(T+7)],   # 1: y_{t-1}
    100 * log.output[7:(T+6)],   # 2: y_{t-2}
    real.interest.rate[8:(T+7)], # 3: r_{t-1}
    real.interest.rate[7:(T+6)], # 4: r_{t-2}
    inflation[8:(T+7)],          # 5: π_{t-1}
    (inflation[7:(T+6)] + inflation[6:(T+5)] + inflation[5:(T+4)]) / 3,   # 6: π̄_{t-2:t-4}
    (inflation[4:(T+3)] + inflation[3:(T+2)] + inflation[2:(T+1)] + inflation[1:T]) / 4   # 7: π̄_{t-5:t-8}
  )

  cat(sprintf("   y.data: %d x %d\n", nrow(y.data), ncol(y.data)))
  cat(sprintf("   x.data: %d x %d\n", nrow(x.data), ncol(x.data)))

  # Criar vetor de kappa para cada observacao (T trimestres)
  anos.est <- anos[9:(T+8)]
  trimestres.est <- trimestres[9:(T+8)]

  # Kappa apenas para COVID (2008 usa dummy)
  kappa.vec <- rep(1.0, T)
  for (t.i in 1:T) {
    ano <- anos.est[t.i]
    tri <- trimestres.est[t.i]
    # 2008 NAO usa kappa, usa dummy via phi_2008
    if (ano == 2020 && tri >= 2) {
      kappa.vec[t.i] <- kappa_2020
    } else if (ano == 2021) {
      kappa.vec[t.i] <- kappa_2021
    } else if (ano == 2022) {
      kappa.vec[t.i] <- kappa_2022
    }
  }

  n.covid <- sum(kappa.vec > 1.0)
  cat(sprintf("   Trimestres com kappa > 1: %d\n", n.covid))

  ##--------------------------------------------------------------------------##
  ## 3. Estimação inicial de a_r por OLS
  ##--------------------------------------------------------------------------##

  cat("\n3. Estimando a_r (sensibilidade à taxa real) via OLS...\n")

  # Output gap do Stage 1
  gap.s1 <- stage1.results$output.gap.smoothed

  # Ajustar tamanho
  T.gap <- length(gap.s1)
  if (T.gap < T) {
    # Gap do Stage 1 começa depois
    offset <- T - T.gap
    gap.for.is <- gap.s1
    r.for.is <- (real.interest.rate[(8+offset):(T+7)] + real.interest.rate[(7+offset):(T+6)]) / 2
    r.for.is <- r.for.is[1:T.gap]
  } else {
    gap.for.is <- gap.s1[1:T]
    r.for.is <- (real.interest.rate[8:(T+7)] + real.interest.rate[7:(T+6)]) / 2
  }

  # IS curve: ỹ_t = a_1·ỹ_{t-1} + a_2·ỹ_{t-2} + a_r·r_{avg} + ε
  # Resíduo após remover a_1 e a_2:
  T.use <- min(length(gap.for.is) - 2, length(r.for.is) - 2)
  y.is <- gap.for.is[3:(T.use+2)]
  y.is.l1 <- gap.for.is[2:(T.use+1)]
  y.is.l2 <- gap.for.is[1:T.use]
  r.avg <- r.for.is[3:(T.use+2)]

  y.is.adj <- y.is - a.1.fixed * y.is.l1 - a.2.fixed * y.is.l2
  a.r.ols <- sum(y.is.adj * r.avg) / sum(r.avg^2)

  cat(sprintf("   a_r (OLS) = %.6f\n", a.r.ols))

  # Para Brasil, a_r frequentemente é negativo ou muito pequeno
  # Isso é documentado na literatura (Torres 2023, etc)
  if (a.r.ols > 0) {
    cat("   NOTA: a_r positivo (incomum, mas possível)\n")
  } else {
    cat("   NOTA: a_r negativo (comum em emergentes - 'puzzle dos juros')\n")
  }

  ##--------------------------------------------------------------------------##
  ## 4. Parâmetros Stage 2
  ##--------------------------------------------------------------------------##

  cat("\n4. Configurando parâmetros para Stage 2...\n")

  # Stage 2 param vector (11 params):
  # [a_1, a_2, a_r, b_π1, b_π2, b_y, g, σ_IS, σ_π, σ_y*, σ_g]

  param.num <- c("a_1"=1, "a_2"=2, "a_r"=3, "b_pi1"=4, "b_pi2"=5,
                 "b_y"=6, "g"=7, "sigma_is"=8, "sigma_pi"=9,
                 "sigma_ystar"=10, "sigma_g"=11)

  # Lambda_g: ratio sigma_g / sigma_y*
  if (is.na(lambda.g)) {
    lambda.g <- 0.05  # Valor inicial típico HLW
    cat(sprintf("   lambda_g = %.4f (valor inicial)\n", lambda.g))
  } else {
    cat(sprintf("   lambda_g = %.4f (fornecido)\n", lambda.g))
  }

  # sigma_g = lambda_g * sigma_y*
  sigma.g.initial <- lambda.g * sigma.ystar.s1

  initial.parameters <- c(
    a.1.fixed,        # 1: a_1 (FIXO)
    a.2.fixed,        # 2: a_2 (FIXO)
    a.r.ols,          # 3: a_r (ESTIMAR)
    b.pi1.fixed,      # 4: b_π1 (FIXO)
    b.pi2.fixed,      # 5: b_π2 (FIXO)
    b.y.fixed,        # 6: b_y (FIXO)
    g.initial,        # 7: g (ESTIMAR)
    sigma.is.s1,      # 8: σ_IS (ESTIMAR)
    sigma.pi.s1,      # 9: σ_π (ESTIMAR)
    sigma.ystar.s1,   # 10: σ_y* (ESTIMAR com restrição)
    sigma.g.initial   # 11: σ_g (ESTIMAR via lambda_g)
  )

  cat("\n   Parâmetros iniciais:\n")
  for (i in 1:length(initial.parameters)) {
    status <- if (i %in% c(1,2,4,5,6)) "FIXO" else "estimar"
    cat(sprintf("   [%2d] %-12s = %10.6f (%s)\n",
                i, names(param.num)[i], initial.parameters[i], status))
  }

  ##--------------------------------------------------------------------------##
  ## 5. Bounds
  ##--------------------------------------------------------------------------##

  cat(sprintf("\n   DEBUG: a.r.constraint = %s (is.na = %s)\n",
              as.character(a.r.constraint), is.na(a.r.constraint)))

  theta.lb <- rep(-Inf, 11)
  theta.ub <- rep(Inf, 11)

  # FIXAR parâmetros do Stage 1
  theta.lb[1] <- a.1.fixed; theta.ub[1] <- a.1.fixed  # a_1
  theta.lb[2] <- a.2.fixed; theta.ub[2] <- a.2.fixed  # a_2
  theta.lb[4] <- b.pi1.fixed; theta.ub[4] <- b.pi1.fixed  # b_π1
  theta.lb[5] <- b.pi2.fixed; theta.ub[5] <- b.pi2.fixed  # b_π2
  theta.lb[6] <- b.y.fixed; theta.ub[6] <- b.y.fixed  # b_y

  # a_r: permitir negativo (comum em emergentes)
  # Moreira & Portugal: a_r <= -0.0025
  if (!is.na(a.r.constraint)) {
    theta.ub[3] <- a.r.constraint
    # Ajustar valor inicial se estiver acima do bound
    # Nota: a.r.constraint é negativo, então precisamos de valor MAIS negativo
    if (initial.parameters[3] > a.r.constraint) {
      initial.parameters[3] <- a.r.constraint * 1.5  # Mais negativo que o limite
      cat(sprintf("\n   Ajustando a_r inicial para %.6f (dentro do bound)\n", initial.parameters[3]))
    }
    cat(sprintf("   Restrição: a_r <= %.4f (upper bound)\n", a.r.constraint))
  }

  # a_r lower bound (limitar quão negativo pode ser)
  if (!is.na(a.r.lower)) {
    theta.lb[3] <- a.r.lower
    # Ajustar valor inicial se estiver abaixo do bound
    if (initial.parameters[3] < a.r.lower) {
      initial.parameters[3] <- (a.r.lower + ifelse(!is.na(a.r.constraint), a.r.constraint, -0.0025)) / 2
      cat(sprintf("\n   Ajustando a_r inicial para %.6f (dentro do lower bound)\n", initial.parameters[3]))
    }
    cat(sprintf("   Restrição: a_r >= %.4f (lower bound)\n", a.r.lower))
  }

  # Variâncias positivas
  theta.lb[8] <- 0.001               # σ_IS
  theta.lb[9] <- 0.001               # σ_π
  theta.lb[10] <- sigma.ystar.lower  # σ_y* (configurable floor)
  theta.lb[11] <- 0.001              # σ_g

  # Ajustar valor inicial de sigma_y* se abaixo do floor
  if (initial.parameters[10] < sigma.ystar.lower) {
    initial.parameters[10] <- sigma.ystar.lower * 1.5
    cat(sprintf("   Ajustando σ_y* inicial para %.4f (acima do floor %.4f)\n",
                initial.parameters[10], sigma.ystar.lower))
  }

  # Restrição sigma_y*
  if (!is.na(sigma.ystar.constraint)) {
    theta.ub[10] <- sigma.ystar.constraint
    initial.parameters[10] <- min(initial.parameters[10], sigma.ystar.constraint * 0.9)
    cat(sprintf("   Restrição: %.4f <= σ_y* <= %.4f\n", sigma.ystar.lower, sigma.ystar.constraint))
  }

  ##--------------------------------------------------------------------------##
  ## 6. Estados iniciais
  ##--------------------------------------------------------------------------##

  if (any(is.na(xi.00))) {
    # HP filter para tendência inicial
    y.og <- log.output[5:(T+8)]
    g.pot <- hpfilter(y.og, freq = 36000, type = "lambda", drift = FALSE)$trend
    g.pot.diff <- diff(g.pot)

    # Stage 2 states: [y*_t, y*_{t-1}, ỹ_t, ỹ_{t-1}, g_t]
    xi.00 <- c(100 * g.pot[4], 100 * g.pot[3], 0, 0, 100 * mean(g.pot.diff))
    cat(sprintf("\n   xi.00 (HP filter): %.2f %.2f %.2f %.2f %.2f\n",
                xi.00[1], xi.00[2], xi.00[3], xi.00[4], xi.00[5]))
  }

  if (any(is.na(P.00))) {
    P.00 <- diag(1, 5)
  }

  ##--------------------------------------------------------------------------##
  ## 7. Kalman Filter Stage 2
  ##--------------------------------------------------------------------------##

  kalman.filter.stage2 <- function(theta, y.data, x.data, xi.00, P.00, param.num, lambda.g, kappa.vec) {

    a.1 <- theta[param.num["a_1"]]
    a.2 <- theta[param.num["a_2"]]
    a.r <- theta[param.num["a_r"]]
    b.pi1 <- theta[param.num["b_pi1"]]
    b.pi2 <- theta[param.num["b_pi2"]]
    b.y <- theta[param.num["b_y"]]
    g <- theta[param.num["g"]]
    sigma.is <- theta[param.num["sigma_is"]]
    sigma.pi <- theta[param.num["sigma_pi"]]
    sigma.ystar <- theta[param.num["sigma_ystar"]]
    sigma.g <- lambda.g * sigma.ystar  # Forçar relação lambda_g

    # Matrizes do sistema de estados
    # Estado: [y*_t, y*_{t-1}, ỹ_t, ỹ_{t-1}, g_t]

    # F: transição
    F <- matrix(0, 5, 5)
    F[1,1] <- 1; F[1,5] <- 1  # y*_t = y*_{t-1} + g_{t-1}
    F[2,1] <- 1               # y*_{t-1} = y*_{t-1}
    F[3,3] <- 0               # ỹ_t (observado)
    F[4,3] <- 1               # ỹ_{t-1} = ỹ_t
    F[5,5] <- 1               # g_t = g_{t-1} + ε_g

    # Q: covariância dos choques (CONSTANTE - não varia com kappa)
    Q <- matrix(0, 5, 5)
    Q[1,1] <- sigma.ystar^2
    Q[5,5] <- sigma.g^2

    # R base: covariância dos erros de observação (será multiplicado por kappa^2)
    R.base <- diag(c(sigma.is^2, sigma.pi^2))

    T.obs <- nrow(y.data)
    n <- length(xi.00)

    xi.tt <- xi.00
    P.tt <- P.00

    ll <- 0
    xi.filtered <- matrix(0, T.obs, n)

    for (t in 1:T.obs) {
      # Aplicar kappa^2 ao R (conforme HLW 2023)
      kappa.t <- kappa.vec[t]
      R <- R.base * (kappa.t^2)

      # Predict
      xi.ttm1 <- F %*% xi.tt
      P.ttm1 <- F %*% P.tt %*% t(F) + Q

      # Observação
      y.t <- y.data[t, ]

      # Predição de y
      # IS: y_t = y*_t + a_1*(y_{t-1} - y*_{t-1}) + a_2*(y_{t-2} - y*_{t-2})
      #         + a_r/2*(r_{t-1} + r_{t-2} - 2*4*g)
      r.gap <- (x.data[t, 3] + x.data[t, 4]) / 2 - 4 * xi.ttm1[5] / 100  # r* = 4*g anualizado

      y.pred.is <- xi.ttm1[1] + a.1 * (x.data[t, 1] - xi.ttm1[2]) +
                   a.2 * (x.data[t, 2] - xi.ttm1[2]) + a.r * r.gap

      # Phillips: π_t = b_π1*π_{t-1} + b_π2*π̄_{t-2:4} + (1-b_π1-b_π2)*π̄_{t-5:8} + b_y*ỹ_{t-1}
      pi.pred <- b.pi1 * x.data[t, 5] + b.pi2 * x.data[t, 6] +
                 (1 - b.pi1 - b.pi2) * x.data[t, 7] + b.y * xi.ttm1[3]

      y.pred <- c(y.pred.is, pi.pred)
      v.t <- y.t - y.pred

      # Matriz H (como estados afetam observações)
      H <- matrix(0, 5, 2)
      H[1, 1] <- 1                  # y*_t entra em y_t
      H[2, 1] <- -a.1               # -a_1 * y*_{t-1}
      H[3, 2] <- b.y                # b_y * ỹ_{t-1} para Phillips

      # Kalman gain
      S.t <- t(H) %*% P.ttm1 %*% H + R
      S.t <- S.t + diag(1e-8, 2)  # Regularização

      K.t <- P.ttm1 %*% H %*% solve(S.t)

      # Update
      xi.tt <- xi.ttm1 + K.t %*% v.t
      P.tt <- P.ttm1 - K.t %*% t(H) %*% P.ttm1

      # Log-likelihood
      ll <- ll - 0.5 * (log(det(S.t)) + t(v.t) %*% solve(S.t) %*% v.t)

      xi.filtered[t, ] <- xi.tt
    }

    return(list(ll = as.numeric(ll), xi.filtered = xi.filtered,
                P.tt = P.tt, xi.tt = xi.tt))
  }

  ##--------------------------------------------------------------------------##
  ## 8. Otimização MLE
  ##--------------------------------------------------------------------------##

  cat("\n5. Iniciando otimização MLE Stage 2...\n")

  f <- function(theta) {
    result <- tryCatch({
      kf <- kalman.filter.stage2(theta, y.data, x.data, xi.00, P.00, param.num, lambda.g, kappa.vec)
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
  ## 9. Extrair resultados
  ##--------------------------------------------------------------------------##

  kf.final <- kalman.filter.stage2(theta, y.data, x.data, xi.00, P.00, param.num, lambda.g, kappa.vec)

  # PIB potencial
  potential.smoothed <- kf.final$xi.filtered[, 1] / 100

  # Output gap
  output.gap.smoothed <- y.data[, 1] - kf.final$xi.filtered[, 1]

  # g (tendência trimestral)
  g.smoothed <- kf.final$xi.filtered[, 5]

  # r* preliminar = 4 * g (anualizado)
  rstar.preliminary <- g.smoothed * 4 / 100  # Em % anual

  ##--------------------------------------------------------------------------##
  ## 10. Resultados
  ##--------------------------------------------------------------------------##

  cat("\n=================================================================\n")
  cat("STAGE 2 WGM CONCLUÍDO!\n")
  cat("=================================================================\n\n")

  cat("Parâmetros estimados:\n")
  cat(sprintf("  a_1 = %.4f (fixo do Stage 1)\n", theta[1]))
  cat(sprintf("  a_2 = %.4f (fixo do Stage 1)\n", theta[2]))
  cat(sprintf("  a_r = %.6f (sensibilidade à taxa real)\n", theta[3]))
  cat(sprintf("  b_π1 = %.4f (fixo)\n", theta[4]))
  cat(sprintf("  b_π2 = %.4f (fixo)\n", theta[5]))
  cat(sprintf("  b_y = %.4f (fixo)\n", theta[6]))
  cat(sprintf("  g = %.4f (tendência trimestral)\n", theta[7]))
  cat(sprintf("  σ_IS = %.4f\n", theta[8]))
  cat(sprintf("  σ_π = %.4f\n", theta[9]))
  cat(sprintf("  σ_y* = %.4f\n", theta[10]))
  cat(sprintf("  σ_g = %.4f (via lambda_g)\n", lambda.g * theta[10]))
  cat(sprintf("  lambda_g = %.4f\n", lambda.g))

  cat(sprintf("\nLog-likelihood = %.4f\n", kf.final$ll))

  cat(sprintf("\nr* preliminar (Stage 2):\n"))
  cat(sprintf("  Média: %.2f%% aa\n", mean(rstar.preliminary) * 100))
  cat(sprintf("  Último: %.2f%% aa\n", tail(rstar.preliminary, 1) * 100))
  cat(sprintf("  Range: [%.2f%%, %.2f%%] aa\n",
              min(rstar.preliminary) * 100, max(rstar.preliminary) * 100))

  cat(sprintf("\nOutput Gap:\n"))
  cat(sprintf("  Média: %.2f%%\n", mean(output.gap.smoothed)))
  cat(sprintf("  Range: [%.2f%%, %.2f%%]\n",
              min(output.gap.smoothed), max(output.gap.smoothed)))

  return(list(
    theta = theta,
    param.num = param.num,
    log.likelihood = kf.final$ll,
    potential.smoothed = potential.smoothed,
    output.gap.smoothed = output.gap.smoothed,
    g.smoothed = g.smoothed,
    rstar.preliminary = rstar.preliminary,
    xi.filtered = kf.final$xi.filtered,
    lambda.g = lambda.g,
    out.mle = nloptr.out,
    y.data = y.data,
    x.data = x.data,
    T = T,
    estimar.kappa = estimar.kappa,
    kappas = c(kappa_2020 = as.numeric(kappa_2020),
               kappa_2021 = as.numeric(kappa_2021),
               kappa_2022 = as.numeric(kappa_2022)),
    kappa.vec = kappa.vec,
    anos.est = anos.est,
    trimestres.est = trimestres.est
  ))
}

cat("\n")
cat("=================================================================\n")
cat("rstar_stage2_com_kappa.R carregado\n")
cat("=================================================================\n")
cat("\n")
cat("METODOLOGIA HLW 2023 COM KAPPA:\n")
cat("  - Kappas fixos do Stage 1\n")
cat("  - kappa^2 aplicado a matriz R (erro de medicao)\n")
cat("\n")
cat("Uso:\n")
cat("  result <- rstar.stage2.com.kappa(log.output, inflation, real.interest.rate,\n")
cat("                                    stage1.results, anos, trimestres)\n")
cat("\n")
