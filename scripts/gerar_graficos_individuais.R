# ==============================================================================
# GERAR GRAFICOS INDIVIDUAIS - SELIC PF+KAPPA e DI-PRE PF+KAPPA
# ==============================================================================
# Gera 7 figuras + dashboard para cada modelo
# Baseado em gerar_graficos_kappa.R com data no nome do arquivo
# ==============================================================================

library(readxl)

wgm_path <- "."
data.arquivo <- format(Sys.Date(), "%Y%m%d")

# ==============================================================================
# FUNCAO: gerar graficos para um modelo
# ==============================================================================

gerar_graficos <- function(resultado, modelo_nome, modelo_label, sufixo) {

  stage3 <- resultado$stage3
  df <- resultado$dados

  # Criar datas para o eixo x
  T <- length(stage3$rstar.smoothed)
  anos_plot <- df$ano[9:(T+8)]
  tri_plot <- df$trimestre[9:(T+8)]
  datas <- as.Date(paste0(anos_plot, "-", (tri_plot-1)*3 + 1, "-01"))

  # Taxa real observada (anualizada)
  taxa_real_aa <- df$taxa_real_exante_tri[9:(T+8)] * 4

  # Postura monetaria (stance)
  stance <- taxa_real_aa - stage3$rstar.smoothed

  # PIB potencial e observado
  pib_obs <- exp(df$ln_pib[9:(T+8)])
  pib_pot <- exp(stage3$potential.smoothed)

  # Kappa serie
  kappa_serie <- rep(1, T)
  if (!is.null(resultado$kappas)) {
    kappa_2020 <- as.numeric(resultado$kappas["kappa_2020"])
    kappa_2021 <- as.numeric(resultado$kappas["kappa_2021"])
    kappa_2022 <- as.numeric(resultado$kappas["kappa_2022"])
  } else {
    kappa_2020 <- 1; kappa_2021 <- 1; kappa_2022 <- 1
  }
  kappa_serie[anos_plot == 2020 & tri_plot >= 2] <- kappa_2020
  kappa_serie[anos_plot == 2021] <- kappa_2021
  kappa_serie[anos_plot == 2022] <- kappa_2022

  # Indices para sombreamento
  crise2008_idx <- which((anos_plot == 2008 & tri_plot == 4) | (anos_plot == 2009 & tri_plot == 1))
  covid_idx <- which(anos_plot >= 2020 & anos_plot <= 2022)

  dir_out <- file.path(wgm_path, "outputs/graficos")

  cat("\n=== Graficos:", modelo_nome, "===\n")

  # --------------------------------------------------------------------------
  # FIG 1: r*
  # --------------------------------------------------------------------------
  fname <- file.path(dir_out, paste0("fig1_rstar_", sufixo, "_", data.arquivo, ".png"))
  png(fname, width = 1000, height = 600, res = 120)
  par(mar = c(4, 4, 3, 1))

  plot(datas, stage3$rstar.smoothed, type = "l", col = "blue", lwd = 2,
       main = paste("Taxa Natural de Juros (r*) -", modelo_label),
       xlab = "", ylab = "% a.a.", ylim = c(0, 12))
  abline(h = mean(stage3$rstar.smoothed), col = "red", lty = 2)
  abline(h = 0, col = "gray80")

  if (length(crise2008_idx) > 0)
    rect(datas[min(crise2008_idx)], -5, datas[max(crise2008_idx)], 20,
         col = rgb(0, 0, 1, 0.15), border = NA)
  if (length(covid_idx) > 0)
    rect(datas[min(covid_idx)], -5, datas[max(covid_idx)], 20,
         col = rgb(1, 0, 0, 0.1), border = NA)

  lines(datas, stage3$rstar.smoothed, col = "blue", lwd = 2)

  legend("topright",
         c(paste("r* (media =", round(mean(stage3$rstar.smoothed), 2), "%)"),
           paste("r* atual =", round(tail(stage3$rstar.smoothed, 1), 2), "%"),
           "Crise 2008 (break)", "COVID (kappa)"),
         col = c("blue", "red", rgb(0, 0, 1, 0.5), rgb(1, 0, 0, 0.4)),
         lty = c(1, 2, NA, NA), lwd = c(2, 1, NA, NA),
         pch = c(NA, NA, 15, 15), pt.cex = 2, bty = "n")
  dev.off()
  cat("  ", basename(fname), "\n")

  # --------------------------------------------------------------------------
  # FIG 2: Decomposicao g e z
  # --------------------------------------------------------------------------
  fname <- file.path(dir_out, paste0("fig2_decomposicao_gz_", sufixo, "_", data.arquivo, ".png"))
  png(fname, width = 1000, height = 600, res = 120)
  par(mar = c(4, 4, 3, 1))

  g_aa <- stage3$g.smoothed * 100
  z_aa <- stage3$z.smoothed

  plot(datas, g_aa, type = "l", col = "darkgreen", lwd = 2,
       main = paste("Decomposicao r* = g + z -", modelo_label),
       xlab = "", ylab = "% a.a.", ylim = range(c(g_aa, z_aa)))
  lines(datas, z_aa, col = "purple", lwd = 2)
  abline(h = 0, col = "gray80")

  if (length(covid_idx) > 0) {
    rect(datas[min(covid_idx)], -10, datas[max(covid_idx)], 20,
         col = rgb(1, 0, 0, 0.1), border = NA)
    lines(datas, g_aa, col = "darkgreen", lwd = 2)
    lines(datas, z_aa, col = "purple", lwd = 2)
  }

  legend("topright",
         c(paste("g (crescimento) - media:", round(mean(g_aa), 2), "%"),
           paste("z (outros fatores) - media:", round(mean(z_aa), 2), "%")),
         col = c("darkgreen", "purple"), lwd = 2, bty = "n")
  dev.off()
  cat("  ", basename(fname), "\n")

  # --------------------------------------------------------------------------
  # FIG 3: r* vs Taxa Real Observada
  # --------------------------------------------------------------------------
  fname <- file.path(dir_out, paste0("fig3_rstar_vs_taxareal_", sufixo, "_", data.arquivo, ".png"))
  png(fname, width = 1000, height = 600, res = 120)
  par(mar = c(4, 4, 3, 1))

  plot(datas, taxa_real_aa, type = "l", col = "gray50", lwd = 1,
       main = paste("r* vs Taxa Real Ex-Ante -", modelo_label),
       xlab = "", ylab = "% a.a.", ylim = range(c(stage3$rstar.smoothed, taxa_real_aa), na.rm = TRUE))
  lines(datas, stage3$rstar.smoothed, col = "blue", lwd = 2)
  abline(h = 0, col = "gray80")

  if (length(covid_idx) > 0) {
    rect(datas[min(covid_idx)], -20, datas[max(covid_idx)], 30,
         col = rgb(1, 0, 0, 0.1), border = NA)
    lines(datas, taxa_real_aa, col = "gray50", lwd = 1)
    lines(datas, stage3$rstar.smoothed, col = "blue", lwd = 2)
  }

  legend("topright",
         c("Taxa real ex-ante observada", "r* estimado"),
         col = c("gray50", "blue"), lwd = c(1, 2), bty = "n")
  dev.off()
  cat("  ", basename(fname), "\n")

  # --------------------------------------------------------------------------
  # FIG 4: Output Gap
  # --------------------------------------------------------------------------
  fname <- file.path(dir_out, paste0("fig4_output_gap_", sufixo, "_", data.arquivo, ".png"))
  png(fname, width = 1000, height = 600, res = 120)
  par(mar = c(4, 4, 3, 1))

  plot(datas, stage3$output.gap.smoothed, type = "l", col = "darkred", lwd = 2,
       main = paste("Hiato do Produto -", modelo_label),
       xlab = "", ylab = "%")
  abline(h = 0, col = "gray50", lty = 2)

  if (length(crise2008_idx) > 0)
    rect(datas[min(crise2008_idx)], -20, datas[max(crise2008_idx)], 20,
         col = rgb(0, 0, 1, 0.15), border = NA)
  if (length(covid_idx) > 0)
    rect(datas[min(covid_idx)], -20, datas[max(covid_idx)], 20,
         col = rgb(1, 0, 0, 0.1), border = NA)

  lines(datas, stage3$output.gap.smoothed, col = "darkred", lwd = 2)

  legend("topright",
         c(paste("Hiato medio:", round(mean(stage3$output.gap.smoothed), 2), "%"),
           paste("Hiato atual:", round(tail(stage3$output.gap.smoothed, 1), 2), "%"),
           "Crise 2008", "COVID"),
         col = c("darkred", "darkred", rgb(0, 0, 1, 0.5), rgb(1, 0, 0, 0.4)),
         lwd = c(2, NA, NA, NA), pch = c(NA, NA, 15, 15), pt.cex = 2, bty = "n")
  dev.off()
  cat("  ", basename(fname), "\n")

  # --------------------------------------------------------------------------
  # FIG 5: Postura Monetaria
  # --------------------------------------------------------------------------
  fname <- file.path(dir_out, paste0("fig5_stance_", sufixo, "_", data.arquivo, ".png"))
  png(fname, width = 1000, height = 600, res = 120)
  par(mar = c(4, 4, 3, 1))

  plot(datas, stance, type = "l", col = "darkorange", lwd = 2,
       main = paste("Postura Monetaria (r - r*) -", modelo_label),
       xlab = "", ylab = "p.p.")
  abline(h = 0, col = "gray50", lty = 2)

  polygon(c(datas, rev(datas)),
          c(pmax(stance, 0), rep(0, length(stance))),
          col = rgb(1, 0, 0, 0.2), border = NA)
  polygon(c(datas, rev(datas)),
          c(pmin(stance, 0), rep(0, length(stance))),
          col = rgb(0, 0, 1, 0.2), border = NA)

  if (length(covid_idx) > 0)
    rect(datas[min(covid_idx)], -20, datas[max(covid_idx)], 20,
         col = rgb(0.5, 0.5, 0.5, 0.2), border = NA)

  lines(datas, stance, col = "darkorange", lwd = 2)

  legend("topright",
         c("Contracionista (r > r*)", "Expansionista (r < r*)",
           paste("Atual:", round(tail(stance, 1), 2), "p.p.")),
         fill = c(rgb(1, 0, 0, 0.3), rgb(0, 0, 1, 0.3), NA),
         border = NA, bty = "n")
  dev.off()
  cat("  ", basename(fname), "\n")

  # --------------------------------------------------------------------------
  # FIG 6: PIB Potencial
  # --------------------------------------------------------------------------
  fname <- file.path(dir_out, paste0("fig6_pib_potencial_", sufixo, "_", data.arquivo, ".png"))
  png(fname, width = 1000, height = 600, res = 120)
  par(mar = c(4, 4, 3, 1))

  pib_obs_idx <- pib_obs / pib_obs[1] * 100
  pib_pot_idx <- pib_pot / pib_pot[1] * 100

  plot(datas, pib_obs_idx, type = "l", col = "gray50", lwd = 1,
       main = paste("PIB Observado vs Potencial -", modelo_label),
       xlab = "", ylab = "Indice (inicio = 100)",
       ylim = range(c(pib_obs_idx, pib_pot_idx)))
  lines(datas, pib_pot_idx, col = "blue", lwd = 2)

  if (length(covid_idx) > 0) {
    rect(datas[min(covid_idx)], 0, datas[max(covid_idx)], 200,
         col = rgb(1, 0, 0, 0.1), border = NA)
    lines(datas, pib_obs_idx, col = "gray50", lwd = 1)
    lines(datas, pib_pot_idx, col = "blue", lwd = 2)
  }

  legend("topleft",
         c("PIB observado", "PIB potencial"),
         col = c("gray50", "blue"), lwd = c(1, 2), bty = "n")
  dev.off()
  cat("  ", basename(fname), "\n")

  # --------------------------------------------------------------------------
  # FIG 7: Kappa
  # --------------------------------------------------------------------------
  fname <- file.path(dir_out, paste0("fig7_kappa_", sufixo, "_", data.arquivo, ".png"))
  png(fname, width = 1000, height = 600, res = 120)
  par(mar = c(4, 4, 3, 1))

  plot(datas, kappa_serie, type = "l", col = "orange", lwd = 2,
       main = paste("Kappa COVID (HLW 2023) -", modelo_label),
       xlab = "", ylab = "kappa", ylim = c(0, max(kappa_serie) * 1.2))
  abline(h = 1, col = "gray50", lty = 2)

  if (length(covid_idx) > 0)
    rect(datas[min(covid_idx)], -5, datas[max(covid_idx)], 25,
         col = rgb(1, 0.5, 0, 0.1), border = NA)

  lines(datas, kappa_serie, col = "orange", lwd = 2)
  idx_covid_kappa <- which((anos_plot == 2020 & tri_plot >= 2) | anos_plot == 2021 | anos_plot == 2022)
  points(datas[idx_covid_kappa], kappa_serie[idx_covid_kappa], col = "red", pch = 19, cex = 0.8)

  legend("topright",
         c(paste("kappa 2020 (Q2-Q4) =", round(kappa_2020, 2)),
           paste("kappa 2021 =", round(kappa_2021, 2)),
           paste("kappa 2022 =", round(kappa_2022, 2)),
           "kappa = 1 (normal)"),
         col = c("red", "red", "red", "gray50"),
         pch = c(19, 19, 19, NA), lty = c(NA, NA, NA, 2), bty = "n")
  dev.off()
  cat("  ", basename(fname), "\n")

  # --------------------------------------------------------------------------
  # DASHBOARD 6 paineis
  # --------------------------------------------------------------------------
  fname <- file.path(dir_out, paste0("dashboard_", sufixo, "_", data.arquivo, ".png"))
  png(fname, width = 1400, height = 1000, res = 120)
  par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))

  # 1. r*
  plot(datas, stage3$rstar.smoothed, type = "l", col = "blue", lwd = 2,
       main = "r* (Taxa Natural)", xlab = "", ylab = "% a.a.")
  abline(h = mean(stage3$rstar.smoothed), col = "red", lty = 2)
  if (length(covid_idx) > 0) {
    rect(datas[min(covid_idx)], -5, datas[max(covid_idx)], 20, col = rgb(1, 0, 0, 0.1), border = NA)
    lines(datas, stage3$rstar.smoothed, col = "blue", lwd = 2)
  }
  legend("topright", c(paste("Media:", round(mean(stage3$rstar.smoothed), 2), "%"),
                        paste("Atual:", round(tail(stage3$rstar.smoothed, 1), 2), "%")),
         bty = "n", cex = 0.8)

  # 2. Decomposicao
  plot(datas, g_aa, type = "l", col = "darkgreen", lwd = 2,
       main = "r* = g + z", xlab = "", ylab = "% a.a.",
       ylim = range(c(g_aa, z_aa)))
  lines(datas, z_aa, col = "purple", lwd = 2)
  if (length(covid_idx) > 0) {
    rect(datas[min(covid_idx)], -10, datas[max(covid_idx)], 20, col = rgb(1, 0, 0, 0.1), border = NA)
    lines(datas, g_aa, col = "darkgreen", lwd = 2)
    lines(datas, z_aa, col = "purple", lwd = 2)
  }
  legend("topright", c("g", "z"), col = c("darkgreen", "purple"), lwd = 2, bty = "n", cex = 0.8)

  # 3. r* vs taxa real
  plot(datas, taxa_real_aa, type = "l", col = "gray50", lwd = 1,
       main = "r* vs Taxa Real", xlab = "", ylab = "% a.a.")
  lines(datas, stage3$rstar.smoothed, col = "blue", lwd = 2)
  if (length(covid_idx) > 0) {
    rect(datas[min(covid_idx)], -20, datas[max(covid_idx)], 30, col = rgb(1, 0, 0, 0.1), border = NA)
    lines(datas, taxa_real_aa, col = "gray50", lwd = 1)
    lines(datas, stage3$rstar.smoothed, col = "blue", lwd = 2)
  }
  legend("topright", c("Obs.", "r*"), col = c("gray50", "blue"), lwd = c(1, 2), bty = "n", cex = 0.8)

  # 4. Output gap
  plot(datas, stage3$output.gap.smoothed, type = "l", col = "darkred", lwd = 2,
       main = "Hiato do Produto", xlab = "", ylab = "%")
  abline(h = 0, col = "gray50", lty = 2)
  if (length(covid_idx) > 0) {
    rect(datas[min(covid_idx)], -20, datas[max(covid_idx)], 20, col = rgb(1, 0, 0, 0.1), border = NA)
    lines(datas, stage3$output.gap.smoothed, col = "darkred", lwd = 2)
  }

  # 5. Stance
  plot(datas, stance, type = "l", col = "darkorange", lwd = 2,
       main = "Postura Monetaria (r - r*)", xlab = "", ylab = "p.p.")
  abline(h = 0, col = "gray50", lty = 2)
  polygon(c(datas, rev(datas)), c(pmax(stance, 0), rep(0, length(stance))),
          col = rgb(1, 0, 0, 0.2), border = NA)
  polygon(c(datas, rev(datas)), c(pmin(stance, 0), rep(0, length(stance))),
          col = rgb(0, 0, 1, 0.2), border = NA)
  lines(datas, stance, col = "darkorange", lwd = 2)

  # 6. Kappa
  plot(datas, kappa_serie, type = "l", col = "orange", lwd = 2,
       main = "Kappa COVID", xlab = "", ylab = "kappa", ylim = c(0, max(kappa_serie) * 1.2))
  abline(h = 1, col = "gray50", lty = 2)
  if (length(covid_idx) > 0) {
    rect(datas[min(covid_idx)], -5, datas[max(covid_idx)], 25, col = rgb(1, 0.5, 0, 0.1), border = NA)
    lines(datas, kappa_serie, col = "orange", lwd = 2)
  }
  points(datas[kappa_serie > 1], kappa_serie[kappa_serie > 1], col = "red", pch = 19, cex = 0.6)

  # Titulo geral
  mtext(modelo_nome, outer = FALSE, side = 3, line = -1.5, cex = 0.9, font = 2, at = par("usr")[1])

  dev.off()
  cat("  ", basename(fname), "\n")

  cat("  -> Total: 8 arquivos\n")
}

# ==============================================================================
# CARREGAR E GERAR
# ==============================================================================

cat("================================================================\n")
cat("GERANDO GRAFICOS INDIVIDUAIS\n")
cat("================================================================\n")

# 1. SELIC PF+KAPPA
selic_file <- file.path(wgm_path, "outputs/results/wgm_kappa_params_fixos_20260212.rds")
cat("\nCarregando:", selic_file, "\n")
selic <- readRDS(selic_file)

gerar_graficos(selic,
               modelo_nome = "Selic - Params Fixos + Kappa COVID",
               modelo_label = "Selic PF+Kappa",
               sufixo = "selic_pf")

# 2. DI-PRE PF+KAPPA
dipre_file <- file.path(wgm_path, "outputs/results/wgm_dipre_params_fixos_20260212.rds")
cat("\nCarregando:", dipre_file, "\n")
dipre <- readRDS(dipre_file)

gerar_graficos(dipre,
               modelo_nome = "DI-Pre B3 M12 - Params Fixos + Kappa COVID",
               modelo_label = "DI-Pre PF+Kappa",
               sufixo = "dipre_pf")

cat("\n================================================================\n")
cat("TODOS OS GRAFICOS INDIVIDUAIS GERADOS!\n")
cat("Data:", data.arquivo, "\n")
cat("================================================================\n")
