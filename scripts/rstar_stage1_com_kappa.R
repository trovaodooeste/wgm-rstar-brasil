##------------------------------------------------------------------------------##
## File:        rstar_stage1_com_kappa.R
##
## Description: Stage 1 HLW para Brasil - COM COVID/STRINGENCY na IS CURVE
##              E COM KAPPA (volatilidade variavel no tempo)
##              E COM BREAKS ESTRUTURAIS (mudanca de tendencia no PIB potencial)
##              Seguindo EXATAMENTE a metodologia HLW 2023
##
## METODOLOGIA HLW 2023 (CORRETA):
##   - COVID/Stringency entra na IS curve (produto), NAO na Phillips
##   - Output gap e ajustado: gap = y - y* - phi*COVID
##   - Phillips usa o gap JA AJUSTADO
##   - KAPPA: volatilidade variavel durante COVID (2020 Q2-Q4, 2021, 2022)
##            kappa^2 aplicado a matriz R (erro de medicao)
##
## BREAKS ESTRUTURAIS (mudanca de tendencia no PIB potencial):
##   - Break 2008-Q4: Crise financeira global (Lehman Brothers)
##   - Break 2020-Q2: COVID-19
##   - Breaks sao rampas: c(rep(0, t_break), 1:(T-t_break))
##   - Diferente de dummy (choque temporario no hiato)
##
## Curva IS (COM COVID):
##   y_gap = phi*d + a_1*(gap_{t-1} - phi*d_{t-1}) + a_2*(gap_{t-2} - phi*d_{t-2})
##
## Output Gap:
##   gap_t = y_t - y*_t - phi*d_t
##   onde y*_t e calculado com breaks estruturais
##
## Curva de Phillips:
##   pi_t = b_pi*pi_{t-1} + (1-b_pi)*pi_bar + b_y*(gap_{t-1} - phi*d_{t-1})
##
## Parametros (13 total com kappa, 10 sem):
##   a_1, a_2, phi, b_pi1, b_pi2, b_y, g, sigma_IS, sigma_pi, sigma_y*
##   kappa_2020, kappa_2021, kappa_2022 (se estimar.kappa=TRUE)
##
## Data: 2026-02-09 (Adicionado BREAKS estruturais, removido phi_2008)
##------------------------------------------------------------------------------##

# Pacotes
if (!require("nloptr")) install.packages("nloptr")
if (!require("numDeriv")) install.packages("numDeriv")
if (!require("mFilter")) install.packages("mFilter")

library(nloptr)
library(numDeriv)
library(mFilter)

##------------------------------------------------------------------------------##
## Funcao principal: rstar.stage1.com.kappa
##------------------------------------------------------------------------------##
rstar.stage1.com.kappa <- function(log.output,
                                        inflation,
                                        covid.dummy,  # Stringency normalizado (0-1)
                                        anos,         # Vetor de anos para cada observacao
                                        trimestres,   # Vetor de trimestres (1-4)
                                        b.y.constraint = 0.25,  # Moreira & Portugal (2019) - Nota 11: b_y >= 0.25 para Brasil
                                        sigma.ystar.constraint = 0.5,
                                        a.sum.constraint = 0.99,  # a_1 + a_2 <= 0.99 (estabilidade IS)
                                        fix.phi = NA,
                                        fix.a1 = NA,   # Fixar a_1 no valor do modelo original
                                        fix.a2 = NA,   # Fixar a_2 no valor do modelo original
                                        fix.sigma.is = NA, # Fixar sigma_is no valor do modelo original
                                        # BREAKS ESTRUTURAIS (mudanca de tendencia no PIB potencial)
                                        use.break.2008 = TRUE,  # Incluir break 2008-Q4?
                                        use.break.2020 = FALSE, # Incluir break 2020-Q2? (usamos COVID dummy em vez disso)
                                        data.start.year = 2003,
                                        data.start.quarter = 1,
                                        # KAPPA (volatilidade variavel)
                                        estimar.kappa = TRUE,  # Estimar kappas ou usar fixos
                                        kappa.2020.init = 5.0,
                                        kappa.2021.init = 2.0,
                                        kappa.2022.init = 1.5,
                                        xi.00 = NA,
                                        P.00 = NA) {

  stage <- 1

  # Dados comecam 8 trimestres antes da estimacao
  T <- length(log.output) - 8

  cat("\n")
  cat("================================================================\n")
  cat("Stage 1 - METODOLOGIA HLW 2023 (COVID na IS curve)\n")
  cat("================================================================\n\n")
  cat("Trimestres para estimacao:", T, "\n")
  cat("COVID/Stringency max:", round(max(covid.dummy) * 100, 1), "%\n")
  cat("Trimestres com COVID > 0:", sum(covid.dummy > 0), "\n\n")

  ##--------------------------------------------------------------------------##
  ## 1. Estimacao inicial por OLS COM BREAKS ESTRUTURAIS
  ##--------------------------------------------------------------------------##

  # Output gap inicial via tendencia linear + BREAKS
  # Base: intercepto + tendencia linear
  x.og <- cbind(rep(1, T+4), 1:(T+4))

  # Calcular indices dos breaks baseado nas datas reais dos dados
  # Break 1: 2008-Q4 (crise financeira Lehman Brothers)
  quarters_to_2008Q4 <- (2008 - data.start.year) * 4 + (4 - data.start.quarter + 1)

  # Break 2: 2020-Q2 (inicio COVID-19 no Brasil)
  quarters_to_2020Q2 <- (2020 - data.start.year) * 4 + (2 - data.start.quarter + 1)

  cat("BREAKS ESTRUTURAIS (mudanca de tendencia no PIB potencial):\n")

  # Adicionar break 2008 se (1) habilitado E (2) dentro do periodo
  if (use.break.2008 && quarters_to_2008Q4 > 0 && quarters_to_2008Q4 < T+4) {
    x.og <- cbind(x.og, c(rep(0, quarters_to_2008Q4),
                          1:(T+4-quarters_to_2008Q4)))
    cat("  Break 2008-Q4: ATIVO (indice", quarters_to_2008Q4, ")\n")
  } else if (!use.break.2008) {
    cat("  Break 2008-Q4: DESABILITADO\n")
  } else {
    cat("  Break 2008-Q4: fora do periodo\n")
  }

  # Adicionar break 2020 se (1) habilitado E (2) dentro do periodo
  if (use.break.2020 && quarters_to_2020Q2 > 0 && quarters_to_2020Q2 < T+4) {
    x.og <- cbind(x.og, c(rep(0, quarters_to_2020Q2),
                          1:(T+4-quarters_to_2020Q2)))
    cat("  Break 2020-Q2: ATIVO (indice", quarters_to_2020Q2, ")\n")
  } else if (!use.break.2020) {
    cat("  Break 2020-Q2: DESABILITADO (usa COVID dummy na IS curve)\n")
  } else {
    cat("  Break 2020-Q2: fora do periodo\n")
  }
  cat("\n")

  y.og <- log.output[5:(T+8)]
  output.gap <- (y.og - x.og %*% solve(t(x.og) %*% x.og, t(x.og) %*% y.og)) * 100

  # Preparar dados para IS curve
  y.is    <- output.gap[5:(T+4)]
  y.is.l1 <- output.gap[4:(T+3)]
  y.is.l2 <- output.gap[3:(T+2)]

  # COVID dummy
  d       <- covid.dummy[9:(T+8)]
  d.l1    <- covid.dummy[8:(T+7)]
  d.l2    <- covid.dummy[7:(T+6)]

  # IS curve com COVID dummy (phi) - SEM phi_2008 (usa break estrutural)
  # y.is = phi*d + a_1*(gap_l1 ajustado) + a_2*(gap_l2 ajustado)
  # gap ajustado = gap - phi*d

  phi.val <- ifelse(is.na(fix.phi), -5, fix.phi)

  # OLS simples para a_1, a_2 (com dummy fixa inicialmente)
  x.is <- cbind(y.is.l1, y.is.l2)
  b.is.ols <- solve(t(x.is) %*% x.is, t(x.is) %*% y.is)

  b.is <- c(phi = phi.val, a_1 = b.is.ols[1], a_2 = b.is.ols[2])

  cat("IS curve (com COVID dummy, SEM phi_2008 - usa break estrutural):\n")
  cat("  phi (COVID) =", round(b.is["phi"], 4), "\n")
  cat("  a_1 =", round(b.is["a_1"], 4), "\n")
  cat("  a_2 =", round(b.is["a_2"], 4), "\n")

  # Residuo IS (considerando apenas COVID dummy)
  gap.ajustado <- y.is - b.is["phi"]*d
  gap.ajustado.l1 <- y.is.l1 - b.is["phi"]*d.l1
  gap.ajustado.l2 <- y.is.l2 - b.is["phi"]*d.l2

  r.is <- gap.ajustado - b.is["a_1"]*gap.ajustado.l1 - b.is["a_2"]*gap.ajustado.l2
  s.is <- sqrt(sum(r.is^2) / (length(r.is) - 3))
  cat("  sigma_IS =", round(s.is, 4), "\n\n")

  # Phillips curve por OLS
  # pi_t = b_pi*pi_{t-1} + ... + b_y*(gap_{t-1} - phi*d_{t-1})
  y.ph <- inflation[9:(T+8)]
  gap.ajustado.l1 <- y.is.l1 - b.is["phi"]*d.l1  # Gap AJUSTADO pelo COVID!

  x.ph <- cbind(
    inflation[8:(T+7)],                                                    # pi_{t-1}
    (inflation[7:(T+6)] + inflation[6:(T+5)] + inflation[5:(T+4)]) / 3,   # pi_bar_{t-2:t-4}
    (inflation[4:(T+3)] + inflation[3:(T+2)] + inflation[2:(T+1)] + inflation[1:T]) / 4,  # pi_bar_{t-5:t-8}
    gap.ajustado.l1                                                        # gap AJUSTADO!
  )

  b.ph <- solve(t(x.ph) %*% x.ph, t(x.ph) %*% y.ph)
  r.ph <- y.ph - x.ph %*% b.ph
  s.ph <- sqrt(sum(r.ph^2) / (length(r.ph) - 4))

  cat("Phillips curve (com gap ajustado pelo COVID):\n")
  cat("  b_pi1 =", round(b.ph[1], 4), "\n")
  cat("  b_pi2 =", round(b.ph[2], 4), "\n")
  cat("  b_y =", round(b.ph[4], 4), "\n")
  cat("  sigma_pi =", round(s.ph, 4), "\n\n")

  ##--------------------------------------------------------------------------##
  ## 2. Preparar dados para Kalman Filter
  ##--------------------------------------------------------------------------##

  # Extrair anos e trimestres para o periodo de estimacao
  anos.est <- anos[9:(T+8)]
  trimestres.est <- trimestres[9:(T+8)]

  cat("Periodo da estimacao:\n")
  cat("  Anos:", min(anos.est), "-", max(anos.est), "\n")
  cat("  N trimestres COVID (kappa 2020 Q2-Q4, 2021, 2022):",
      sum((anos.est == 2020 & trimestres.est >= 2) |
          anos.est == 2021 |
          anos.est == 2022), "\n\n")

  # y.data: [100*log(y), pi]
  y.data <- cbind(
    100 * log.output[9:(T+8)],
    inflation[9:(T+8)]
  )

  # x.data: lagged variables
  # [y_{t-1}, y_{t-2}, pi_{t-1}, pi_bar_{t-2:4}, pi_bar_{t-5:8}, d_t, d_{t-1}, d_{t-2}]
  x.data <- cbind(
    100 * log.output[8:(T+7)],   # 1: y_{t-1}
    100 * log.output[7:(T+6)],   # 2: y_{t-2}
    inflation[8:(T+7)],          # 3: pi_{t-1}
    (inflation[7:(T+6)] + inflation[6:(T+5)] + inflation[5:(T+4)]) / 3,   # 4: pi_bar_{t-2:t-4}
    (inflation[4:(T+3)] + inflation[3:(T+2)] + inflation[2:(T+1)] + inflation[1:T]) / 4,  # 5: pi_bar_{t-5:t-8}
    covid.dummy[9:(T+8)],        # 6: d_t (COVID)
    covid.dummy[8:(T+7)],        # 7: d_{t-1}
    covid.dummy[7:(T+6)]         # 8: d_{t-2}
  )

  ##--------------------------------------------------------------------------##
  ## 3. Parametros iniciais (SEM phi_2008)
  ##--------------------------------------------------------------------------##

  # Numero de parametros:
  # Base: a_1, a_2, phi, b_pi1, b_pi2, b_y, g, sigma_is, sigma_pi, sigma_ystar (10)
  # Com kappa: + kappa_2020, kappa_2021, kappa_2022 (13)

  n.params.base <- 10
  if (estimar.kappa) {
    n.params <- 13
    param.num <- c("a_1"=1, "a_2"=2, "phi"=3, "b_pi1"=4, "b_pi2"=5,
                   "b_y"=6, "g"=7, "sigma_is"=8, "sigma_pi"=9, "sigma_ystar"=10,
                   "kappa_2020"=11, "kappa_2021"=12, "kappa_2022"=13)
  } else {
    n.params <- 10
    param.num <- c("a_1"=1, "a_2"=2, "phi"=3, "b_pi1"=4, "b_pi2"=5,
                   "b_y"=6, "g"=7, "sigma_is"=8, "sigma_pi"=9, "sigma_ystar"=10)
  }

  # Valores iniciais (respeitando bounds)
  a.ub <- ifelse(!is.na(a.sum.constraint), a.sum.constraint / 2 + 0.25, 0.99)
  initial.parameters <- c(
    min(b.is["a_1"], a.ub * 0.95),  # 1: a_1 (respeitando bound)
    min(b.is["a_2"], a.ub * 0.95),  # 2: a_2 (respeitando bound)
    b.is["phi"],          # 3: phi (COVID)
    b.ph[1],              # 4: b_pi1
    b.ph[2],              # 5: b_pi2
    max(b.ph[4], b.y.constraint),  # 6: b_y (com constraint)
    0.5,                  # 7: g (tendencia)
    s.is,                 # 8: sigma_IS
    s.ph,                 # 9: sigma_pi
    0.5                   # 10: sigma_y*
  )

  if (estimar.kappa) {
    initial.parameters <- c(initial.parameters,
                            kappa.2020.init,   # 11: kappa_2020
                            kappa.2021.init,   # 12: kappa_2021
                            kappa.2022.init)   # 13: kappa_2022
  }

  names(initial.parameters) <- names(param.num)

  cat("Parametros iniciais:\n")
  print(round(initial.parameters, 4))
  cat("\n")

  ##--------------------------------------------------------------------------##
  ## 4. Bounds e restricoes
  ##--------------------------------------------------------------------------##

  theta.lb <- rep(-Inf, n.params)
  theta.ub <- rep(Inf, n.params)

  # b_y >= 0.25 (Moreira & Portugal 2019, Nota 11)
  if (!is.na(b.y.constraint)) {
    theta.lb[param.num["b_y"]] <- b.y.constraint
    cat("Restricao: b_y >=", b.y.constraint, "\n")
  }

  # a_1 + a_2 <= 0.99 (estabilidade IS)
  if (!is.na(a.sum.constraint)) {
    # Implementar como bound individual (simplificacao)
    a.ub <- a.sum.constraint / 2 + 0.25  # ~0.745
    theta.ub[param.num["a_1"]] <- a.ub
    theta.ub[param.num["a_2"]] <- a.ub
    cat("Restricao: a_1 + a_2 <=", a.sum.constraint, "(upper bound individual:", round(a.ub, 3), ")\n")
  }

  # sigma_y* <= constraint (HLW)
  if (!is.na(sigma.ystar.constraint)) {
    theta.ub[param.num["sigma_ystar"]] <- sigma.ystar.constraint
    initial.parameters[param.num["sigma_ystar"]] <- min(initial.parameters[param.num["sigma_ystar"]],
                                                         sigma.ystar.constraint * 0.9)
    cat("Restricao: sigma_y* <=", sigma.ystar.constraint, "\n")
  }

  # phi (COVID) deve ser negativo (COVID reduz produto)
  theta.ub[param.num["phi"]] <- 0

  # Variancias positivas
  theta.lb[param.num["sigma_is"]] <- 0.001
  theta.lb[param.num["sigma_pi"]] <- 0.001
  theta.lb[param.num["sigma_ystar"]] <- 0.001

  # Kappa bounds (se estimando)
  if (estimar.kappa) {
    theta.lb[param.num["kappa_2020"]] <- 1.0  # Variancia nao pode diminuir
    theta.lb[param.num["kappa_2021"]] <- 1.0
    theta.lb[param.num["kappa_2022"]] <- 1.0
    theta.ub[param.num["kappa_2020"]] <- 20.0  # Limite superior
    theta.ub[param.num["kappa_2021"]] <- 20.0
    theta.ub[param.num["kappa_2022"]] <- 20.0
    cat("Restricao: kappa >= 1.0 (variancia nao pode diminuir)\n")
    cat("           kappa <= 20.0 (limite superior)\n")
  }

  # Fixar parametros se solicitado
  if (!is.na(fix.phi)) {
    theta.lb[param.num["phi"]] <- fix.phi
    theta.ub[param.num["phi"]] <- fix.phi
    initial.parameters[param.num["phi"]] <- fix.phi
    cat("FIXO: phi =", fix.phi, "\n")
  }

  if (!is.na(fix.a1)) {
    theta.lb[param.num["a_1"]] <- fix.a1
    theta.ub[param.num["a_1"]] <- fix.a1
    initial.parameters[param.num["a_1"]] <- fix.a1
    cat("FIXO: a_1 =", fix.a1, "\n")
  }

  if (!is.na(fix.a2)) {
    theta.lb[param.num["a_2"]] <- fix.a2
    theta.ub[param.num["a_2"]] <- fix.a2
    initial.parameters[param.num["a_2"]] <- fix.a2
    cat("FIXO: a_2 =", fix.a2, "\n")
  }

  if (!is.na(fix.sigma.is)) {
    theta.lb[param.num["sigma_is"]] <- fix.sigma.is
    theta.ub[param.num["sigma_is"]] <- fix.sigma.is
    initial.parameters[param.num["sigma_is"]] <- fix.sigma.is
    cat("FIXO: sigma_is =", fix.sigma.is, "\n")
  }

  ##--------------------------------------------------------------------------##
  ## 5. Estados iniciais
  ##--------------------------------------------------------------------------##

  if (any(is.na(xi.00))) {
    # HP filter para tendencia inicial
    g.pot <- hpfilter(log.output[5:(T+8)], freq = 36000, type = "lambda", drift = FALSE)$trend
    xi.00 <- c(100 * g.pot[4], 100 * g.pot[3], 0)
    cat("\nxi.00 (HP filter):", round(xi.00, 2), "\n")
  }

  if (any(is.na(P.00))) {
    P.00 <- diag(c(1, 1, 1))
  }

  ##--------------------------------------------------------------------------##
  ## 6. Kalman Filter com KAPPA
  ##--------------------------------------------------------------------------##

  kalman.filter.stage1.kappa <- function(theta, y.data, x.data, xi.00, P.00, param.num,
                                          anos.est, trimestres.est, estimar.kappa) {

    # Extrair parametros
    a_1 <- theta[param.num["a_1"]]
    a_2 <- theta[param.num["a_2"]]
    phi <- theta[param.num["phi"]]
    b_pi1 <- theta[param.num["b_pi1"]]
    b_pi2 <- theta[param.num["b_pi2"]]
    b_y <- theta[param.num["b_y"]]
    g <- theta[param.num["g"]]
    sigma_is <- theta[param.num["sigma_is"]]
    sigma_pi <- theta[param.num["sigma_pi"]]
    sigma_ystar <- theta[param.num["sigma_ystar"]]

    # Kappa (apenas COVID - 2008 usa break estrutural)
    if (estimar.kappa) {
      kappa_2020 <- theta[param.num["kappa_2020"]]
      kappa_2021 <- theta[param.num["kappa_2021"]]
      kappa_2022 <- theta[param.num["kappa_2022"]]
    } else {
      kappa_2020 <- 1.0
      kappa_2021 <- 1.0
      kappa_2022 <- 1.0
    }

    T.obs <- nrow(y.data)
    n.states <- 3  # [y*_t, y*_{t-1}, gap_{t-1}]

    # Matrizes do sistema
    # Estado: [y*_t, y*_{t-1}, gap_{t-1}]
    # F: transicao
    F <- matrix(0, n.states, n.states)
    F[1,1] <- 1  # y*_t = y*_{t-1} + g + e
    F[2,1] <- 1  # y*_{t-1} = y*_t (lag)
    F[3,3] <- 0  # gap_{t-1} sera calculado

    # Q: covariancia dos choques (CONSTANTE - nao varia com kappa)
    Q <- matrix(0, n.states, n.states)
    Q[1,1] <- sigma_ystar^2  # Variancia do PIB potencial

    # R base: erros de observacao (sera multiplicado por kappa^2)
    R.base <- diag(c(sigma_is^2, sigma_pi^2))

    xi.tt <- xi.00
    P.tt <- P.00

    ll <- 0
    xi.filtered <- matrix(0, T.obs, n.states)

    for (t.i in 1:T.obs) {
      # Determinar kappa para este periodo (apenas COVID)
      ano <- anos.est[t.i]
      tri <- trimestres.est[t.i]

      # 2008 NAO usa kappa, usa break estrutural
      if (ano == 2020 && tri >= 2) {
        kappa.t <- kappa_2020
      } else if (ano == 2021) {
        kappa.t <- kappa_2021
      } else if (ano == 2022) {
        kappa.t <- kappa_2022
      } else {
        kappa.t <- 1.0
      }

      # Aplicar kappa^2 ao R (conforme HLW 2023)
      R <- R.base * (kappa.t^2)

      # Dados observados
      y.t <- y.data[t.i, ]
      y.l1 <- x.data[t.i, 1]
      y.l2 <- x.data[t.i, 2]

      # COVID dummy
      d.t <- x.data[t.i, 6]
      d.t.l1 <- x.data[t.i, 7]
      d.t.l2 <- x.data[t.i, 8]

      # Predict
      ystar.pred <- xi.tt[1] + g
      xi.ttm1 <- c(ystar.pred, xi.tt[1], xi.tt[3])
      P.ttm1 <- F %*% P.tt %*% t(F) + Q

      # Observacao IS: y_t = y*_t + phi*d_t + a_1*(gap_{t-1} - phi*d_{t-1}) + a_2*(gap_{t-2} - phi*d_{t-2}) + e_IS
      ystar.l1 <- xi.ttm1[2]
      ystar.l2 <- xi.ttm1[2]  # Aproximacao

      gap.l1.adj <- y.l1 - ystar.l1 - phi*d.t.l1
      gap.l2.adj <- y.l2 - ystar.l2 - phi*d.t.l2

      y.pred <- ystar.pred + phi*d.t + a_1*gap.l1.adj + a_2*gap.l2.adj

      # Observacao Phillips: pi_t = b_pi1*pi_{t-1} + b_pi2*pi_bar_{t-2:4} + (1-b_pi1-b_pi2)*pi_bar_{t-5:8} + b_y*gap_{t-1}
      pi.pred <- b_pi1*x.data[t.i, 3] + b_pi2*x.data[t.i, 4] +
                 (1 - b_pi1 - b_pi2)*x.data[t.i, 5] + b_y*gap.l1.adj

      y.obs.pred <- c(y.pred, pi.pred)
      v.t <- y.t - y.obs.pred

      # H matrix (linearizada)
      H <- matrix(0, n.states, 2)
      H[1, 1] <- 1  # dy/dy*_t = 1
      H[2, 1] <- -(a_1 + a_2)  # dy/dy*_{t-1}
      H[3, 2] <- b_y  # dpi/dgap_{t-1}

      # Kalman gain
      S.t <- t(H) %*% P.ttm1 %*% H + R
      S.t <- S.t + diag(1e-8, 2)  # Regularizacao

      K.t <- P.ttm1 %*% H %*% solve(S.t)

      # Update
      xi.tt <- xi.ttm1 + K.t %*% v.t
      P.tt <- P.ttm1 - K.t %*% t(H) %*% P.ttm1

      # Atualizar gap_{t-1}
      xi.tt[3] <- y.t[1] - xi.tt[1] - phi*d.t

      # Log-likelihood
      ll <- ll - 0.5 * (log(det(S.t)) + t(v.t) %*% solve(S.t) %*% v.t)

      xi.filtered[t.i, ] <- xi.tt
    }

    return(list(ll = as.numeric(ll), xi.filtered = xi.filtered,
                P.tt = P.tt, xi.tt = xi.tt))
  }

  ##--------------------------------------------------------------------------##
  ## 7. Otimizacao MLE
  ##--------------------------------------------------------------------------##

  cat("\nIniciando otimizacao MLE...\n")

  f <- function(theta) {
    result <- tryCatch({
      kf <- kalman.filter.stage1.kappa(theta, y.data, x.data, xi.00, P.00, param.num,
                                        anos.est, trimestres.est, estimar.kappa)
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
  names(theta) <- names(param.num)

  cat("\nStatus:", nloptr.out$message, "\n")

  ##--------------------------------------------------------------------------##
  ## 8. Extrair resultados
  ##--------------------------------------------------------------------------##

  kf.final <- kalman.filter.stage1.kappa(theta, y.data, x.data, xi.00, P.00, param.num,
                                          anos.est, trimestres.est, estimar.kappa)

  # Criar vetor de kappa para cada observacao
  kappa.vec <- rep(1.0, T)
  for (t.i in 1:T) {
    ano <- anos.est[t.i]
    tri <- trimestres.est[t.i]
    if (ano == 2020 && tri >= 2) {
      kappa.vec[t.i] <- theta[param.num["kappa_2020"]]
    } else if (ano == 2021) {
      kappa.vec[t.i] <- theta[param.num["kappa_2021"]]
    } else if (ano == 2022) {
      kappa.vec[t.i] <- theta[param.num["kappa_2022"]]
    }
  }

  # PIB potencial
  potential.smoothed <- kf.final$xi.filtered[, 1] / 100

  # Output gap: y - y* - phi*COVID
  d <- covid.dummy[9:(T+8)]
  output.gap.smoothed <- y.data[, 1] - kf.final$xi.filtered[, 1] - theta[param.num["phi"]] * d

  ##--------------------------------------------------------------------------##
  ## 9. Resultados
  ##--------------------------------------------------------------------------##

  cat("\n")
  cat("================================================================\n")
  cat("PARAMETROS ESTIMADOS (Metodologia HLW 2023)\n")
  cat("================================================================\n\n")

  cat("IS curve:\n")
  cat("  a_1 =", round(theta[param.num["a_1"]], 4), "\n")
  cat("  a_2 =", round(theta[param.num["a_2"]], 4), "\n")
  cat("  phi (COVID) =", round(theta[param.num["phi"]], 4), "\n")
  cat("  (2008 usa BREAK estrutural, nao phi)\n\n")

  cat("Phillips curve:\n")
  cat("  b_pi1 =", round(theta[param.num["b_pi1"]], 4), "\n")
  cat("  b_pi2 =", round(theta[param.num["b_pi2"]], 4), "\n")
  cat("  b_y =", round(theta[param.num["b_y"]], 4), "\n\n")

  cat("Outros:\n")
  cat("  g =", round(theta[param.num["g"]], 4), "\n")
  cat("  sigma_IS =", round(theta[param.num["sigma_is"]], 4), "\n")
  cat("  sigma_pi =", round(theta[param.num["sigma_pi"]], 4), "\n")
  cat("  sigma_y* =", round(theta[param.num["sigma_ystar"]], 4), "\n\n")

  if (estimar.kappa) {
    cat("Kappas COVID (volatilidade variavel - HLW 2023):\n")
    cat("  kappa_2020 (Q2-Q4)       =", round(theta[param.num["kappa_2020"]], 2), "\n")
    cat("  kappa_2021               =", round(theta[param.num["kappa_2021"]], 2), "\n")
    cat("  kappa_2022               =", round(theta[param.num["kappa_2022"]], 2), "\n")
    cat("  (2008 usa break estrutural, nao kappa)\n\n")
  }

  cat("Output gap medio:", round(mean(output.gap.smoothed), 2), "%\n")
  cat("(Ajustado por: gap = y - y* - phi*COVID)\n")
  cat("(PIB potencial calculado com break estrutural 2008)\n\n")

  # Extrair kappas para retorno
  if (estimar.kappa) {
    kappas <- c(kappa_2020 = as.numeric(theta[param.num["kappa_2020"]]),
                kappa_2021 = as.numeric(theta[param.num["kappa_2021"]]),
                kappa_2022 = as.numeric(theta[param.num["kappa_2022"]]))
  } else {
    kappas <- c(kappa_2020 = 1.0, kappa_2021 = 1.0, kappa_2022 = 1.0)
  }

  return(list(
    theta = theta,
    param.num = param.num,
    log.likelihood = kf.final$ll,
    potential.smoothed = potential.smoothed,
    output.gap.smoothed = output.gap.smoothed,
    xi.filtered = kf.final$xi.filtered,
    out.mle = nloptr.out,
    y.data = y.data,
    x.data = x.data,
    T = T,
    phi = theta[param.num["phi"]],
    use.break.2008 = use.break.2008,
    use.break.2020 = use.break.2020,
    estimar.kappa = estimar.kappa,
    kappas = kappas,
    kappa.vec = kappa.vec,
    anos.est = anos.est,
    trimestres.est = trimestres.est
  ))
}

cat("\n")
cat("================================================================\n")
cat("rstar_stage1_com_kappa.R carregado\n")
cat("================================================================\n")
cat("\n")
cat("METODOLOGIA HLW 2023 COM KAPPA E BREAKS:\n")
cat("  - COVID/Stringency entra na IS curve (produto)\n")
cat("  - Output gap = y - y* - phi*COVID\n")
cat("  - Phillips usa gap AJUSTADO\n")
cat("  - BREAK 2008: mudanca de tendencia no PIB potencial\n")
cat("  - KAPPA: volatilidade variavel durante COVID (2020 Q2-Q4, 2021, 2022)\n")
cat("  - kappa^2 aplicado a matriz R (erro de medicao)\n")
cat("\n")
cat("Uso:\n")
cat("  result <- rstar.stage1.com.kappa(log.output, inflation, covid.dummy,\n")
cat("                                    anos, trimestres, use.break.2008=TRUE)\n")
cat("\n")
