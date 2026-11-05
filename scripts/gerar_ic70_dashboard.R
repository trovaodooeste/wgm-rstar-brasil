# ==============================================================================
# IC 70% - Dashboard e resultados salvos
# Carrega SEs do IC 95% (ja calculados) e aplica z70 = qnorm(0.85)
# Salva: outputs/results/ic_hamilton_parametrico_selic_70_YYYYMMDD.rds
#        outputs/graficos/dashboard_ic70_YYYYMMDD.png
#        outputs/graficos/fig_rstar_ic70_YYYYMMDD.png
#        outputs/graficos/fig_decomposicao_ic70_YYYYMMDD.png
# ==============================================================================

wgm_path <- "."

# z-scores
z70 <- qnorm(0.85)   # = 1.0364
z95 <- 1.96

cat(sprintf("z70 = %.4f  |  z95 = %.4f\n", z70, z95))

# Carregar IC 95% (contem se.rstar, se.g, se.z, rstar.base, etc.)
ic <- readRDS(file.path(wgm_path, "outputs/results/ic_hamilton_parametrico_selic_20260213.rds"))

datas      <- ic$datas
rstar.base <- ic$rstar.base
g.base     <- ic$g.base        # em fracao (dividido por 100 la dentro)
z.base     <- ic$z.base
se.rstar   <- ic$se.rstar
se.g       <- ic$se.g
se.z       <- ic$se.z
taxa.real  <- ic$taxa.real * 4  # converter de % trimestral para % a.a. (consistente com r*)
gap        <- ic$gap
n          <- length(rstar.base)
anos.est   <- ic$anos.est
tri.est    <- ic$tri.est

# g em % a.a.
g.pct <- g.base * 100

# ==============================================================================
# CALCULAR ICs 70%
# ==============================================================================

rstar.u70 <- rstar.base + z70 * se.rstar
rstar.l70 <- rstar.base - z70 * se.rstar

g.u70 <- g.pct + z70 * se.g
g.l70 <- g.pct - z70 * se.g

z.u70 <- z.base + z70 * se.z
z.l70 <- z.base - z70 * se.z

# ==============================================================================
# TABELAS NO CONSOLE
# ==============================================================================

cat("\n================================================================\n")
cat("TABELA: r* POR PERIODO - IC 70% vs IC 95%\n")
cat("================================================================\n")
cat(sprintf("%-14s  %6s  %5s  %20s  %20s\n",
    "Periodo", "r*", "SE", "IC 70%", "IC 95%"))
cat(paste0(rep("-", 72), collapse=""), "\n")

periodos <- list(
  "2005-2008"  = anos.est >= 2005 & anos.est <= 2008,
  "2009-2013"  = anos.est >= 2009 & anos.est <= 2013,
  "2014-2019"  = anos.est >= 2014 & anos.est <= 2019,
  "2020-2022"  = anos.est >= 2020 & anos.est <= 2022,
  "2023-atual" = anos.est >= 2023,
  "Amostra"    = rep(TRUE, n)
)

for (nm in names(periodos)) {
  idx <- periodos[[nm]]
  rm  <- mean(rstar.base[idx])
  sm  <- mean(se.rstar[idx])
  cat(sprintf("%-14s  %5.2f%%  %4.2f  [%5.2f%%, %5.2f%%]    [%5.2f%%, %5.2f%%]\n",
      nm, rm, sm,
      rm - z70*sm, rm + z70*sm,
      rm - z95*sm, rm + z95*sm))
}

cat("\n----------------------------------------------------------------\n")
cat(sprintf("Ultimo trimestre: %d-Q%d\n", tail(anos.est,1), tail(tri.est,1)))
cat(sprintf("  r*    = %.2f%%\n", tail(rstar.base,1)))
cat(sprintf("  SE    = %.2f pp\n", tail(se.rstar,1)))
cat(sprintf("  IC 70%%: [%.2f%%, %.2f%%]\n",
    tail(rstar.base,1) - z70*tail(se.rstar,1),
    tail(rstar.base,1) + z70*tail(se.rstar,1)))
cat(sprintf("  IC 95%%: [%.2f%%, %.2f%%]\n",
    tail(rstar.base,1) - z95*tail(se.rstar,1),
    tail(rstar.base,1) + z95*tail(se.rstar,1)))

cat("\n--- g e z ATUAIS ---\n")
cat(sprintf("  g: %.2f%%  |  70%%: [%.2f%%, %.2f%%]  |  95%%: [%.2f%%, %.2f%%]\n",
    tail(g.pct,1),
    tail(g.pct,1)-z70*tail(se.g,1), tail(g.pct,1)+z70*tail(se.g,1),
    tail(g.pct,1)-z95*tail(se.g,1), tail(g.pct,1)+z95*tail(se.g,1)))
cat(sprintf("  z: %.2f%%  |  70%%: [%.2f%%, %.2f%%]  |  95%%: [%.2f%%, %.2f%%]\n",
    tail(z.base,1),
    tail(z.base,1)-z70*tail(se.z,1), tail(z.base,1)+z70*tail(se.z,1),
    tail(z.base,1)-z95*tail(se.z,1), tail(z.base,1)+z95*tail(se.z,1)))

# ==============================================================================
# SALVAR RDS 70%
# ==============================================================================

data.arquivo <- format(Sys.Date(), "%Y%m%d")

ic70 <- list(
  # SEs (mesmos do 95%)
  se.rstar = se.rstar,
  se.g     = se.g,
  se.z     = se.z,
  # Base
  rstar.base = rstar.base,
  g.base     = g.base,
  z.base     = z.base,
  # IC 70%
  rstar.upper70 = rstar.u70,
  rstar.lower70 = rstar.l70,
  g.upper70     = g.u70,
  g.lower70     = g.l70,
  z.upper70     = z.u70,
  z.lower70     = z.l70,
  # IC 95% (para referencia)
  rstar.upper95 = rstar.base + z95 * se.rstar,
  rstar.lower95 = rstar.base - z95 * se.rstar,
  g.upper95     = g.pct + z95 * se.g,
  g.lower95     = g.pct - z95 * se.g,
  z.upper95     = z.base + z95 * se.z,
  z.lower95     = z.base - z95 * se.z,
  # Metadados
  datas    = datas,
  anos.est = anos.est,
  tri.est  = tri.est,
  taxa.real = taxa.real,
  gap      = gap,
  z70      = z70,
  z95      = z95,
  bse      = ic$bse,
  info     = ic$info,
  niter    = ic$niter,
  n.excluded = ic$n.excluded
)

rds.path <- file.path(wgm_path, paste0("outputs/results/ic_hamilton_parametrico_selic_70_", data.arquivo, ".rds"))
saveRDS(ic70, rds.path)
cat(sprintf("\nRDS 70%% salvo: outputs/results/ic_hamilton_parametrico_selic_70_%s.rds\n", data.arquivo))

# ==============================================================================
# GRAFICO 1: r* COM IC 70%
# ==============================================================================

png(file.path(wgm_path, paste0("outputs/graficos/fig_rstar_ic70_", data.arquivo, ".png")),
    width = 1200, height = 700, res = 130)

par(mar = c(4, 4.5, 3, 1))

ylim.r <- range(c(rstar.l70 * 0.95, rstar.u70 * 1.05, taxa.real), na.rm = TRUE)

plot(datas, rstar.base, type = "n",
     main = expression(paste("Taxa Natural de Juros (", r^"*", ") com IC 70% Parametrico")),
     ylab = "% a.a.", xlab = "",
     ylim = ylim.r)

# Sombreamento COVID
rect(as.Date("2020-04-01"), par("usr")[3],
     as.Date("2022-12-01"), par("usr")[4],
     col = rgb(0.5, 0.5, 0.5, 0.08), border = NA)
text(as.Date("2021-01-01"), ylim.r[2] * 0.97, "COVID", cex = 0.7, col = "gray40")

# IC 70%
polygon(c(datas, rev(datas)),
        c(rstar.u70, rev(rstar.l70)),
        col = rgb(1, 0, 0, 0.20), border = NA)
lines(datas, rstar.u70, col = "red", lty = 2, lwd = 0.7)
lines(datas, rstar.l70, col = "red", lty = 2, lwd = 0.7)

# Taxa real e r*
lines(datas, taxa.real, col = "gray50", lwd = 1, lty = 1)
lines(datas, rstar.base, col = "red", lwd = 2.5)

legend("topright",
       c(expression(r^"*"), "IC 70%", "Taxa real observada"),
       col = c("red", rgb(1, 0, 0, 0.3), "gray50"),
       lwd = c(2.5, 8, 1), lty = c(1, 1, 1), bty = "n", cex = 0.85)

dev.off()
cat(sprintf("Grafico 1 salvo: fig_rstar_ic70_%s.png\n", data.arquivo))

# ==============================================================================
# GRAFICO 2: DECOMPOSICAO g + z COM IC 70%
# ==============================================================================

png(file.path(wgm_path, paste0("outputs/graficos/fig_decomposicao_ic70_", data.arquivo, ".png")),
    width = 1400, height = 700, res = 130)

par(mfrow = c(1, 2), mar = c(4, 4.5, 3, 1))

# g com IC 70%
plot(datas, g.pct, type = "n",
     main = "g (crescimento potencial) com IC 70%",
     ylab = "% a.a.", xlab = "",
     ylim = range(c(g.l70, g.u70)))
polygon(c(datas, rev(datas)), c(g.u70, rev(g.l70)),
        col = rgb(0, 0, 1, 0.18), border = NA)
lines(datas, g.u70, col = "blue", lty = 2, lwd = 0.7)
lines(datas, g.l70, col = "blue", lty = 2, lwd = 0.7)
lines(datas, g.pct, col = "blue", lwd = 2.5)
abline(h = 0, col = "gray", lty = 3)

# z com IC 70%
plot(datas, z.base, type = "n",
     main = "z (premio de risco e outros) com IC 70%",
     ylab = "% a.a.", xlab = "",
     ylim = range(c(z.l70, z.u70)))
polygon(c(datas, rev(datas)), c(z.u70, rev(z.l70)),
        col = rgb(0, 0.5, 0, 0.18), border = NA)
lines(datas, z.u70, col = "darkgreen", lty = 2, lwd = 0.7)
lines(datas, z.l70, col = "darkgreen", lty = 2, lwd = 0.7)
lines(datas, z.base, col = "darkgreen", lwd = 2.5)
abline(h = 0, col = "gray", lty = 3)

dev.off()
cat(sprintf("Grafico 2 salvo: fig_decomposicao_ic70_%s.png\n", data.arquivo))

# ==============================================================================
# GRAFICO 3: DASHBOARD COMPLETO (5 paineis, identico ao 95% mas com 70%)
# ==============================================================================

png(file.path(wgm_path, paste0("outputs/graficos/dashboard_ic70_", data.arquivo, ".png")),
    width = 1600, height = 1000, res = 130)

layout(matrix(c(1,1,2,3,4,5), 3, 2, byrow = TRUE), heights = c(1.3, 1, 1))
par(mar = c(4, 4.5, 3, 1))

# --- Painel 1 (largo): r* com IC 70% ---
ylim.r <- range(c(rstar.l70 * 0.95, rstar.u70 * 1.05, taxa.real), na.rm = TRUE)
plot(datas, rstar.base, type = "n",
     main = expression(paste(r^"*", " com Intervalo de Confianca 70% Parametrico")),
     ylab = "% a.a.", xlab = "",
     ylim = ylim.r)

rect(as.Date("2020-04-01"), par("usr")[3],
     as.Date("2022-12-01"), par("usr")[4],
     col = rgb(0.5, 0.5, 0.5, 0.08), border = NA)
text(as.Date("2021-01-01"), ylim.r[2] * 0.96, "COVID", cex = 0.7, col = "gray40")

polygon(c(datas, rev(datas)), c(rstar.u70, rev(rstar.l70)),
        col = rgb(1, 0, 0, 0.20), border = NA)
lines(datas, rstar.u70, col = "red", lty = 2, lwd = 0.7)
lines(datas, rstar.l70, col = "red", lty = 2, lwd = 0.7)
lines(datas, taxa.real, col = "gray50", lwd = 1)
lines(datas, rstar.base, col = "red", lwd = 2.5)

legend("topright",
       c(expression(r^"*"), "IC 70%", "Taxa real"),
       col = c("red", rgb(1, 0, 0, 0.3), "gray50"),
       lwd = c(2.5, 8, 1), bty = "n", cex = 0.85)

# --- Painel 2: g com IC 70% ---
plot(datas, g.pct, type = "n",
     main = "g (crescimento potencial)", ylab = "% a.a.", xlab = "",
     ylim = range(c(g.l70, g.u70)))
polygon(c(datas, rev(datas)), c(g.u70, rev(g.l70)),
        col = rgb(0, 0, 1, 0.18), border = NA)
lines(datas, g.u70, col = "blue", lty = 2, lwd = 0.7)
lines(datas, g.l70, col = "blue", lty = 2, lwd = 0.7)
lines(datas, g.pct, col = "blue", lwd = 2)
abline(h = 0, col = "gray", lty = 3)

# --- Painel 3: z com IC 70% ---
plot(datas, z.base, type = "n",
     main = "z (outros fatores)", ylab = "% a.a.", xlab = "",
     ylim = range(c(z.l70, z.u70)))
polygon(c(datas, rev(datas)), c(z.u70, rev(z.l70)),
        col = rgb(0, 0.5, 0, 0.18), border = NA)
lines(datas, z.u70, col = "darkgreen", lty = 2, lwd = 0.7)
lines(datas, z.l70, col = "darkgreen", lty = 2, lwd = 0.7)
lines(datas, z.base, col = "darkgreen", lwd = 2)

# --- Painel 4: Hiato do produto ---
plot(datas, gap, type = "l", col = "purple", lwd = 2,
     main = "Hiato do produto", ylab = "%", xlab = "")
abline(h = 0, col = "gray", lty = 2)
rect(as.Date("2020-04-01"), par("usr")[3],
     as.Date("2022-12-01"), par("usr")[4],
     col = rgb(0.5, 0.5, 0.5, 0.08), border = NA)

# --- Painel 5: SE ao longo do tempo ---
plot(datas, se.rstar, type = "l", col = "darkred", lwd = 2,
     main = expression(paste("Erro Padrao de ", r^"*")),
     ylab = "SE (pp)", xlab = "")
abline(h = mean(se.rstar), col = "gray", lty = 2)
text(datas[5], mean(se.rstar) * 1.15,
     paste0("Media: ", round(mean(se.rstar), 2), " pp"),
     col = "gray30", cex = 0.8)

dev.off()
cat(sprintf("Grafico 3 (dashboard) salvo: dashboard_ic70_%s.png\n", data.arquivo))

# ==============================================================================
# RESUMO FINAL
# ==============================================================================

cat("\n================================================================\n")
cat("RESUMO - IC 70% (z = 1.036)\n")
cat("================================================================\n")
cat(sprintf("r* ATUAL (%d-Q%d): %.2f%% [%.2f%%, %.2f%%]\n",
    tail(anos.est,1), tail(tri.est,1),
    tail(rstar.base,1),
    tail(rstar.l70,1), tail(rstar.u70,1)))
cat(sprintf("r* MEDIO (amostra): %.2f%% [%.2f%%, %.2f%%]\n",
    mean(rstar.base),
    mean(rstar.base) - z70*mean(se.rstar),
    mean(rstar.base) + z70*mean(se.rstar)))
cat(sprintf("SE medio r*: %.2f pp\n", mean(se.rstar)))
cat(sprintf("Largura media IC 70%%: %.2f pp  (vs %.2f pp no IC 95%%)\n",
    2*z70*mean(se.rstar), 2*z95*mean(se.rstar)))
cat(sprintf("\nRDS salvo: %s\n", basename(rds.path)))
cat("3 graficos gerados (fig_rstar, fig_decomposicao, dashboard)\n")
cat("================================================================\n")
