# ==============================================================================
# GRAFICOS DE COMPARACAO - 3 MODELOS
# ==============================================================================
# 1. SELIC original (sem kappa, sem break 2008)
# 2. SELIC PF+Kappa (params fixos + break 2008 + kappa COVID)
# 3. DI-Pre PF+Kappa (idem, com DI-Pre B3 M12)
# ==============================================================================

library(readxl)

wgm_path <- "."
data.arquivo <- format(Sys.Date(), "%Y%m%d")
dir_out <- file.path(wgm_path, "outputs/graficos")

cat("================================================================\n")
cat("COMPARACAO 3 MODELOS\n")
cat("================================================================\n")

# Carregar
orig <- readRDS(file.path(wgm_path, "outputs/results/wgm_exante_FINAL_z0_3_20250122.rds"))
selic <- readRDS(file.path(wgm_path, "outputs/results/wgm_kappa_params_fixos_20260212.rds"))
dipre <- readRDS(file.path(wgm_path, "outputs/results/wgm_dipre_params_fixos_20260212.rds"))

cat("Modelos carregados.\n")

# --------------------------------------------------------------------------
# Preparar series - usar menor comprimento comum
# --------------------------------------------------------------------------

# Original
T_o <- length(orig$stage3$rstar.smoothed)
df_o <- orig$dados
anos_o <- df_o$ano[9:(T_o+8)]
tri_o <- df_o$trimestre[9:(T_o+8)]
datas_o <- as.Date(paste0(anos_o, "-", (tri_o-1)*3+1, "-01"))
rstar_o <- orig$stage3$rstar.smoothed
g_o <- orig$stage3$g.smoothed * 100
z_o <- orig$stage3$z.smoothed
gap_o <- orig$stage3$output.gap.smoothed
taxa_o <- df_o$taxa_real_exante_tri[9:(T_o+8)] * 4

# Selic PF
T_s <- length(selic$stage3$rstar.smoothed)
df_s <- selic$dados
anos_s <- df_s$ano[9:(T_s+8)]
tri_s <- df_s$trimestre[9:(T_s+8)]
datas_s <- as.Date(paste0(anos_s, "-", (tri_s-1)*3+1, "-01"))
rstar_s <- selic$stage3$rstar.smoothed
g_s <- selic$stage3$g.smoothed * 100
z_s <- selic$stage3$z.smoothed
gap_s <- selic$stage3$output.gap.smoothed
taxa_s <- df_s$taxa_real_exante_tri[9:(T_s+8)] * 4

# DI-Pre PF
T_d <- length(dipre$stage3$rstar.smoothed)
df_d <- dipre$dados
anos_d <- df_d$ano[9:(T_d+8)]
tri_d <- df_d$trimestre[9:(T_d+8)]
datas_d <- as.Date(paste0(anos_d, "-", (tri_d-1)*3+1, "-01"))
rstar_d <- dipre$stage3$rstar.smoothed
g_d <- dipre$stage3$g.smoothed * 100
z_d <- dipre$stage3$z.smoothed
gap_d <- dipre$stage3$output.gap.smoothed
taxa_d <- df_d$taxa_real_exante_tri[9:(T_d+8)] * 4

# Indices COVID e 2008 (usar datas_s como referencia)
covid_idx_s <- which(anos_s >= 2020 & anos_s <= 2022)
crise_idx_s <- which((anos_s == 2008 & tri_s == 4) | (anos_s == 2009 & tri_s == 1))

# Cores
col_orig <- "gray40"
col_selic <- "blue"
col_dipre <- "red"
lwd_orig <- 1.5
lwd_new <- 2

# Labels
lab_orig <- "Selic Original"
lab_selic <- "Selic PF+Kappa"
lab_dipre <- "DI-Pre PF+Kappa"

sombrear <- function(datas, covid_idx, crise_idx, ylims) {
  if (length(crise_idx) > 0)
    rect(datas[min(crise_idx)], ylims[1]-10, datas[max(crise_idx)], ylims[2]+10,
         col = rgb(0, 0, 1, 0.08), border = NA)
  if (length(covid_idx) > 0)
    rect(datas[min(covid_idx)], ylims[1]-10, datas[max(covid_idx)], ylims[2]+10,
         col = rgb(1, 0, 0, 0.08), border = NA)
}

cat("Gerando graficos de comparacao...\n")

# ==============================================================================
# FIG 1: r* - 3 modelos
# ==============================================================================
fname <- file.path(dir_out, paste0("comp_rstar_3modelos_", data.arquivo, ".png"))
png(fname, width = 1100, height = 650, res = 120)
par(mar = c(4, 4, 3, 1))

ylims <- range(c(rstar_o, rstar_s, rstar_d), na.rm = TRUE)
plot(datas_o, rstar_o, type = "l", col = col_orig, lwd = lwd_orig, lty = 2,
     main = "Taxa Natural de Juros (r*) - Comparacao",
     xlab = "", ylab = "% a.a.", ylim = ylims)
sombrear(datas_s, covid_idx_s, crise_idx_s, ylims)
lines(datas_o, rstar_o, col = col_orig, lwd = lwd_orig, lty = 2)
lines(datas_s, rstar_s, col = col_selic, lwd = lwd_new)
lines(datas_d, rstar_d, col = col_dipre, lwd = lwd_new, lty = 4)
abline(h = 0, col = "gray80")

legend("topright",
       c(paste(lab_orig, "- media:", round(mean(rstar_o), 2), "%"),
         paste(lab_selic, "- media:", round(mean(rstar_s), 2), "%"),
         paste(lab_dipre, "- media:", round(mean(rstar_d), 2), "%"),
         "Crise 2008", "COVID"),
       col = c(col_orig, col_selic, col_dipre, rgb(0,0,1,0.3), rgb(1,0,0,0.3)),
       lwd = c(lwd_orig, lwd_new, lwd_new, NA, NA),
       lty = c(2, 1, 4, NA, NA),
       pch = c(NA, NA, NA, 15, 15), pt.cex = 1.5, bty = "n", cex = 0.8)
dev.off()
cat("  ", basename(fname), "\n")

# ==============================================================================
# FIG 2: g (tendencia) - 3 modelos
# ==============================================================================
fname <- file.path(dir_out, paste0("comp_g_3modelos_", data.arquivo, ".png"))
png(fname, width = 1100, height = 650, res = 120)
par(mar = c(4, 4, 3, 1))

ylims <- range(c(g_o, g_s, g_d), na.rm = TRUE)
plot(datas_o, g_o, type = "l", col = col_orig, lwd = lwd_orig, lty = 2,
     main = "Tendencia de Crescimento (g) - Comparacao",
     xlab = "", ylab = "% a.a.", ylim = ylims)
sombrear(datas_s, covid_idx_s, crise_idx_s, ylims)
lines(datas_o, g_o, col = col_orig, lwd = lwd_orig, lty = 2)
lines(datas_s, g_s, col = col_selic, lwd = lwd_new)
lines(datas_d, g_d, col = col_dipre, lwd = lwd_new, lty = 4)
abline(h = 0, col = "gray80")

legend("topright",
       c(paste(lab_orig, "- media:", round(mean(g_o), 2), "%"),
         paste(lab_selic, "- media:", round(mean(g_s), 2), "%"),
         paste(lab_dipre, "- media:", round(mean(g_d), 2), "%")),
       col = c(col_orig, col_selic, col_dipre),
       lwd = c(lwd_orig, lwd_new, lwd_new), lty = c(2, 1, 4), bty = "n", cex = 0.8)
dev.off()
cat("  ", basename(fname), "\n")

# ==============================================================================
# FIG 3: z (outros fatores) - 3 modelos
# ==============================================================================
fname <- file.path(dir_out, paste0("comp_z_3modelos_", data.arquivo, ".png"))
png(fname, width = 1100, height = 650, res = 120)
par(mar = c(4, 4, 3, 1))

ylims <- range(c(z_o, z_s, z_d), na.rm = TRUE)
plot(datas_o, z_o, type = "l", col = col_orig, lwd = lwd_orig, lty = 2,
     main = "Componente z (outros fatores) - Comparacao",
     xlab = "", ylab = "% a.a.", ylim = ylims)
sombrear(datas_s, covid_idx_s, crise_idx_s, ylims)
lines(datas_o, z_o, col = col_orig, lwd = lwd_orig, lty = 2)
lines(datas_s, z_s, col = col_selic, lwd = lwd_new)
lines(datas_d, z_d, col = col_dipre, lwd = lwd_new, lty = 4)
abline(h = 0, col = "gray80")

legend("topright",
       c(paste(lab_orig, "- media:", round(mean(z_o), 2), "%"),
         paste(lab_selic, "- media:", round(mean(z_s), 2), "%"),
         paste(lab_dipre, "- media:", round(mean(z_d), 2), "%")),
       col = c(col_orig, col_selic, col_dipre),
       lwd = c(lwd_orig, lwd_new, lwd_new), lty = c(2, 1, 4), bty = "n", cex = 0.8)
dev.off()
cat("  ", basename(fname), "\n")

# ==============================================================================
# FIG 4: Output Gap - 3 modelos
# ==============================================================================
fname <- file.path(dir_out, paste0("comp_gap_3modelos_", data.arquivo, ".png"))
png(fname, width = 1100, height = 650, res = 120)
par(mar = c(4, 4, 3, 1))

ylims <- range(c(gap_o, gap_s, gap_d), na.rm = TRUE)
plot(datas_o, gap_o, type = "l", col = col_orig, lwd = lwd_orig, lty = 2,
     main = "Hiato do Produto - Comparacao",
     xlab = "", ylab = "%", ylim = ylims)
sombrear(datas_s, covid_idx_s, crise_idx_s, ylims)
lines(datas_o, gap_o, col = col_orig, lwd = lwd_orig, lty = 2)
lines(datas_s, gap_s, col = col_selic, lwd = lwd_new)
lines(datas_d, gap_d, col = col_dipre, lwd = lwd_new, lty = 4)
abline(h = 0, col = "gray50", lty = 2)

legend("topright",
       c(paste(lab_orig, "- media:", round(mean(gap_o), 2), "%"),
         paste(lab_selic, "- media:", round(mean(gap_s), 2), "%"),
         paste(lab_dipre, "- media:", round(mean(gap_d), 2), "%")),
       col = c(col_orig, col_selic, col_dipre),
       lwd = c(lwd_orig, lwd_new, lwd_new), lty = c(2, 1, 4), bty = "n", cex = 0.8)
dev.off()
cat("  ", basename(fname), "\n")

# ==============================================================================
# FIG 5: Stance monetario - 3 modelos
# ==============================================================================
fname <- file.path(dir_out, paste0("comp_stance_3modelos_", data.arquivo, ".png"))
png(fname, width = 1100, height = 650, res = 120)
par(mar = c(4, 4, 3, 1))

stance_o <- taxa_o - rstar_o
stance_s <- taxa_s - rstar_s
stance_d <- taxa_d - rstar_d

ylims <- range(c(stance_o, stance_s, stance_d), na.rm = TRUE)
plot(datas_o, stance_o, type = "l", col = col_orig, lwd = lwd_orig, lty = 2,
     main = "Postura Monetaria (r - r*) - Comparacao",
     xlab = "", ylab = "p.p.", ylim = ylims)
sombrear(datas_s, covid_idx_s, crise_idx_s, ylims)
lines(datas_o, stance_o, col = col_orig, lwd = lwd_orig, lty = 2)
lines(datas_s, stance_s, col = col_selic, lwd = lwd_new)
lines(datas_d, stance_d, col = col_dipre, lwd = lwd_new, lty = 4)
abline(h = 0, col = "gray50", lty = 2)

legend("topright",
       c(paste(lab_orig, "- atual:", round(tail(stance_o, 1), 2), "pp"),
         paste(lab_selic, "- atual:", round(tail(stance_s, 1), 2), "pp"),
         paste(lab_dipre, "- atual:", round(tail(stance_d, 1), 2), "pp")),
       col = c(col_orig, col_selic, col_dipre),
       lwd = c(lwd_orig, lwd_new, lwd_new), lty = c(2, 1, 4), bty = "n", cex = 0.8)
dev.off()
cat("  ", basename(fname), "\n")

# ==============================================================================
# DASHBOARD - 6 paineis comparativos
# ==============================================================================
fname <- file.path(dir_out, paste0("dashboard_comparacao_3modelos_", data.arquivo, ".png"))
png(fname, width = 1500, height = 1000, res = 120)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))

# 1. r*
ylims <- range(c(rstar_o, rstar_s, rstar_d))
plot(datas_o, rstar_o, type = "l", col = col_orig, lwd = lwd_orig, lty = 2,
     main = "r* (Taxa Natural)", xlab = "", ylab = "% a.a.", ylim = ylims)
sombrear(datas_s, covid_idx_s, crise_idx_s, ylims)
lines(datas_o, rstar_o, col = col_orig, lwd = lwd_orig, lty = 2)
lines(datas_s, rstar_s, col = col_selic, lwd = lwd_new)
lines(datas_d, rstar_d, col = col_dipre, lwd = lwd_new, lty = 4)
legend("topright", c("Orig", "Selic PF", "DI-Pre PF"),
       col = c(col_orig, col_selic, col_dipre), lwd = c(lwd_orig, lwd_new, lwd_new),
       lty = c(2, 1, 4), bty = "n", cex = 0.7)

# 2. g
ylims <- range(c(g_o, g_s, g_d))
plot(datas_o, g_o, type = "l", col = col_orig, lwd = lwd_orig, lty = 2,
     main = "g (tendencia crescimento)", xlab = "", ylab = "% a.a.", ylim = ylims)
sombrear(datas_s, covid_idx_s, crise_idx_s, ylims)
lines(datas_o, g_o, col = col_orig, lwd = lwd_orig, lty = 2)
lines(datas_s, g_s, col = col_selic, lwd = lwd_new)
lines(datas_d, g_d, col = col_dipre, lwd = lwd_new, lty = 4)

# 3. z
ylims <- range(c(z_o, z_s, z_d))
plot(datas_o, z_o, type = "l", col = col_orig, lwd = lwd_orig, lty = 2,
     main = "z (outros fatores)", xlab = "", ylab = "% a.a.", ylim = ylims)
sombrear(datas_s, covid_idx_s, crise_idx_s, ylims)
lines(datas_o, z_o, col = col_orig, lwd = lwd_orig, lty = 2)
lines(datas_s, z_s, col = col_selic, lwd = lwd_new)
lines(datas_d, z_d, col = col_dipre, lwd = lwd_new, lty = 4)
abline(h = 0, col = "gray80")

# 4. Gap
ylims <- range(c(gap_o, gap_s, gap_d))
plot(datas_o, gap_o, type = "l", col = col_orig, lwd = lwd_orig, lty = 2,
     main = "Hiato do Produto", xlab = "", ylab = "%", ylim = ylims)
sombrear(datas_s, covid_idx_s, crise_idx_s, ylims)
lines(datas_o, gap_o, col = col_orig, lwd = lwd_orig, lty = 2)
lines(datas_s, gap_s, col = col_selic, lwd = lwd_new)
lines(datas_d, gap_d, col = col_dipre, lwd = lwd_new, lty = 4)
abline(h = 0, col = "gray50", lty = 2)

# 5. Stance
stance_o <- taxa_o - rstar_o
stance_s <- taxa_s - rstar_s
stance_d <- taxa_d - rstar_d
ylims <- range(c(stance_o, stance_s, stance_d), na.rm = TRUE)
plot(datas_o, stance_o, type = "l", col = col_orig, lwd = lwd_orig, lty = 2,
     main = "Postura Monetaria (r - r*)", xlab = "", ylab = "p.p.", ylim = ylims)
sombrear(datas_s, covid_idx_s, crise_idx_s, ylims)
lines(datas_o, stance_o, col = col_orig, lwd = lwd_orig, lty = 2)
lines(datas_s, stance_s, col = col_selic, lwd = lwd_new)
lines(datas_d, stance_d, col = col_dipre, lwd = lwd_new, lty = 4)
abline(h = 0, col = "gray50", lty = 2)

# 6. Tabela resumo
plot.new()
plot.window(xlim = c(0, 10), ylim = c(0, 10))

# Titulo
text(5, 9.5, "Resumo Comparativo", font = 2, cex = 1.1)

# Cabecalho
text(1.5, 8.5, "Metrica", font = 2, cex = 0.8, adj = 0)
text(5.0, 8.5, "Original", font = 2, cex = 0.8)
text(7.0, 8.5, "Selic PF", font = 2, cex = 0.8, col = col_selic)
text(9.0, 8.5, "DI-Pre PF", font = 2, cex = 0.8, col = col_dipre)
segments(0.5, 8.2, 9.8, 8.2, col = "gray50")

# Dados
metricas <- c("r* medio", "r* atual", "g medio", "z medio", "Gap medio",
              "kappa 2020", "kappa 2021", "kappa 2022")
vals_o <- c(sprintf("%.2f%%", mean(rstar_o)), sprintf("%.2f%%", tail(rstar_o,1)),
            sprintf("%.2f%%", mean(g_o)), sprintf("%.2f%%", mean(z_o)),
            sprintf("%.2f%%", mean(gap_o)), "-", "-", "-")
vals_s <- c(sprintf("%.2f%%", mean(rstar_s)), sprintf("%.2f%%", tail(rstar_s,1)),
            sprintf("%.2f%%", mean(g_s)), sprintf("%.2f%%", mean(z_s)),
            sprintf("%.2f%%", mean(gap_s)),
            sprintf("%.2f", selic$kappas["kappa_2020"]),
            sprintf("%.2f", selic$kappas["kappa_2021"]),
            sprintf("%.2f", selic$kappas["kappa_2022"]))
vals_d <- c(sprintf("%.2f%%", mean(rstar_d)), sprintf("%.2f%%", tail(rstar_d,1)),
            sprintf("%.2f%%", mean(g_d)), sprintf("%.2f%%", mean(z_d)),
            sprintf("%.2f%%", mean(gap_d)),
            sprintf("%.2f", dipre$kappas["kappa_2020"]),
            sprintf("%.2f", dipre$kappas["kappa_2021"]),
            sprintf("%.2f", dipre$kappas["kappa_2022"]))

for (i in 1:length(metricas)) {
  y <- 8.0 - i * 0.85
  text(1.5, y, metricas[i], cex = 0.75, adj = 0)
  text(5.0, y, vals_o[i], cex = 0.75)
  text(7.0, y, vals_s[i], cex = 0.75, col = col_selic)
  text(9.0, y, vals_d[i], cex = 0.75, col = col_dipre)
}

mtext("Comparacao: Selic Original vs Selic PF+Kappa vs DI-Pre PF+Kappa",
      outer = TRUE, cex = 1.1, font = 2)

dev.off()
cat("  ", basename(fname), "\n")

cat("\n================================================================\n")
cat("GRAFICOS DE COMPARACAO GERADOS!\n")
cat("Data:", data.arquivo, "\n")
cat("Arquivos em:", dir_out, "\n")
cat("================================================================\n")
