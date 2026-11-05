##------------------------------------------------------------------------------##
## File:        kalman_smoother.R
##
## Description: Rauch-Tung-Striebel (RTS) smoother
##              Backward pass que usa informacao futura dentro da amostra
##              para refinar as estimativas dos estados.
##
##              P_{t|T} <= P_{t|t} sempre (smoother reduz incerteza)
##              No ultimo periodo: P_{T|T} = P_{T|T} (sem melhoria)
##
## Algoritmo:
##   Para t = T-1, T-2, ..., 1:
##     J_t = P_{t|t} * F' * ginv(P_{t+1|t})
##     xi_{t|T} = xi_{t|t} + J_t * (xi_{t+1|T} - xi_{t+1|t})
##     P_{t|T} = P_{t|t} + J_t * (P_{t+1|T} - P_{t+1|t}) * J_t'
##
## Nota: Usa ginv (pseudoinversa) em vez de solve porque P_{t+1|t} pode
##       ser singular quando estados degenerados (gap com Q=0, F=0).
##
## Data: 2026-02-16
##------------------------------------------------------------------------------##

if (!require("MASS")) install.packages("MASS")

kalman.smoother.rts <- function(xi.tt.array, P.tt.array,
                                 xi.ttm1.array, P.ttm1.array,
                                 F.mat) {
  # xi.tt.array:   [T x n] estados filtrados
  # P.tt.array:    [n x n x T] covariancia filtrada
  # xi.ttm1.array: [T x n] estados previstos (antes do update)
  # P.ttm1.array:  [n x n x T] covariancia prevista
  # F.mat:         [n x n] matriz de transicao (constante)

  T.obs <- nrow(xi.tt.array)
  n.states <- ncol(xi.tt.array)

  # Saida: estados e covariancias suavizados
  xi.tT <- matrix(0, T.obs, n.states)
  P.tT <- array(0, dim = c(n.states, n.states, T.obs))

  # Ultimo periodo: suavizado = filtrado (sem informacao futura)
  xi.tT[T.obs, ] <- xi.tt.array[T.obs, ]
  P.tT[, , T.obs] <- P.tt.array[, , T.obs]

  # Backward pass
  for (t in (T.obs - 1):1) {
    # P_{t+1|t} — covariancia prevista do proximo periodo
    P.ttm1.tp1 <- P.ttm1.array[, , t + 1]

    # J_t = P_{t|t} * F' * ginv(P_{t+1|t})
    # ginv para lidar com estados degenerados (gap: Q=0, F=0)
    J.t <- P.tt.array[, , t] %*% t(F.mat) %*% MASS::ginv(P.ttm1.tp1)

    # Estado suavizado
    xi.tT[t, ] <- xi.tt.array[t, ] +
      J.t %*% (xi.tT[t + 1, ] - xi.ttm1.array[t + 1, ])

    # Covariancia suavizada
    P.tT[, , t] <- P.tt.array[, , t] +
      J.t %*% (P.tT[, , t + 1] - P.ttm1.array[, , t + 1]) %*% t(J.t)
  }

  return(list(xi.tT = xi.tT, P.tT = P.tT))
}
