# ==============================================================================
# GRAFICO EXPLICATIVO: Por que o break estrutural em 2008 melhora o modelo
# ==============================================================================
# Mostra:
# 1. PIB real com tendencias (com/sem break) - a motivacao visual
# 2. z inflado sem break vs z contido com break
# 3. r* com/sem break
# ==============================================================================

library(readxl)

wgm_path <- "."

# Carregar dados e resultados
orig <- readRDS(file.path(wgm_path, "outputs/results/wgm_exante_FINAL_z0_3_20250122.rds"))
atual <- readRDS(file.path(wgm_path, "outputs/results/wgm_kappa_params_fixos_20260212.rds"))

df_orig <- orig$dados
df_atual <- atual$dados

# ==============================================================================
# RECALCULAR AS TENDENCIAS OLS (como o Stage 1 faz)
# ==============================================================================

# Dados de log(PIB) - usar os do modelo atual
ln_pib <- as.numeric(df_atual$ln_pib)
anos <- as.numeric(df_atual$ano)
trimestres <- as.numeric(df_atual$trimestre)

T_total <- length(ln_pib) - 8
y.og <- ln_pib[5:(T_total+8)]
n_og <- length(y.og)
anos_og <- anos[5:(T_total+8)]
tri_og <- trimestres[5:(T_total+8)]
datas_og <- as.Date(paste0(anos_og, "-", (tri_og-1)*3+1, "-01"))

# Tendencia SEM break (apenas intercepto + tendencia linear)
x.sem <- cbind(rep(1, n_og), 1:n_og)
beta.sem <- solve(t(x.sem) %*% x.sem, t(x.sem) %*% y.og)
trend.sem <- x.sem %*% beta.sem

# Tendencia COM break 2008-Q4
data.start.year <- 2003
data.start.quarter <- 1
quarters_to_2008Q4 <- (2008 - data.start.year) * 4 + (4 - data.start.quarter + 1)
trend.break <- c(rep(0, quarters_to_2008Q4), 1:(n_og - quarters_to_2008Q4))
x.com <- cbind(rep(1, n_og), 1:n_og, trend.break)
beta.com <- solve(t(x.com) %*% x.com, t(x.com) %*% y.og)
trend.com <- x.com %*% beta.com

# Taxas de crescimento implicitas (% a.a.)
g_sem <- beta.sem[2] * 400  # trimestral -> anual, em log -> %
g_com_antes <- beta.com[2] * 400
g_com_depois <- (beta.com[2] + beta.com[3]) * 400

cat("Crescimento potencial implicito:\n")
cat(sprintf("  Sem break: %.2f%% a.a. (constante)\n", g_sem))
cat(sprintf("  Com break (antes 2008): %.2f%% a.a.\n", g_com_antes))
cat(sprintf("  Com break (depois 2008): %.2f%% a.a.\n", g_com_depois))
cat(sprintf("  Desaceleracao: %.2f pp\n", g_com_depois - g_com_antes))

# ==============================================================================
# GRAFICO
# ==============================================================================

data.arquivo <- format(Sys.Date(), "%Y%m%d")

png(file.path(wgm_path, paste0("outputs/graficos/break_2008_explicativo_", data.arquivo, ".png")),
    width = 1600, height = 1000, res = 130)

layout(matrix(c(1,1,2,3), 2, 2, byrow=TRUE), heights=c(1.3, 1))
par(mar = c(4, 4.5, 3, 1))

# --------------------------------------------------------------------------
# PAINEL 1 (grande): PIB real + tendencias com e sem break
# --------------------------------------------------------------------------

plot(datas_og, y.og, type = "l", col = "gray40", lwd = 1.5,
     main = "PIB Real e Tendencias Estimadas (Stage 1 OLS)",
     ylab = "log(PIB) x 100", xlab = "",
     ylim = range(c(y.og, trend.sem, trend.com)))

# Sombreamento pre/pos 2008
data_break <- as.Date("2008-10-01")
rect(min(datas_og), par("usr")[3], data_break, par("usr")[4],
     col = rgb(0.7, 0.85, 1, 0.2), border = NA)
rect(data_break, par("usr")[3], max(datas_og), par("usr")[4],
     col = rgb(1, 0.85, 0.7, 0.2), border = NA)

# Re-plotar PIB por cima do sombreamento
lines(datas_og, y.og, col = "gray40", lwd = 1.5)

# Tendencia sem break
lines(datas_og, trend.sem, col = "blue", lwd = 2.5, lty = 2)

# Tendencia com break
lines(datas_og, trend.com, col = "red", lwd = 2.5)

# Linha vertical do break
abline(v = data_break, col = "black", lwd = 1.5, lty = 3)
text(data_break, max(y.og) * 0.9999, " Break 2008-Q4", pos = 4, cex = 0.9, font = 2)

# Anotacoes de crescimento
text(as.Date("2006-01-01"), trend.com[which(datas_og == as.Date("2006-01-01"))] + 0.8,
     paste0("g = ", round(g_com_antes, 1), "% a.a."),
     col = "red", cex = 0.85, font = 2)
text(as.Date("2018-01-01"), trend.com[which.min(abs(datas_og - as.Date("2018-01-01")))] - 0.8,
     paste0("g = ", round(g_com_depois, 1), "% a.a."),
     col = "red", cex = 0.85, font = 2)
text(as.Date("2015-01-01"), trend.sem[which.min(abs(datas_og - as.Date("2015-01-01")))] + 0.8,
     paste0("g = ", round(g_sem, 1), "% a.a."),
     col = "blue", cex = 0.85, font = 2)

legend("bottomright",
       c("PIB Real", "Tendencia sem break (reta unica)", "Tendencia com break 2008-Q4"),
       col = c("gray40", "blue", "red"), lwd = c(1.5, 2.5, 2.5), lty = c(1, 2, 1),
       bty = "n", cex = 0.85)

# --------------------------------------------------------------------------
# PAINEL 2: z (premio de risco) - com vs sem break
# --------------------------------------------------------------------------

n <- min(length(orig$stage3$z.smoothed), length(atual$stage3$z.smoothed))
z_orig <- orig$stage3$z.smoothed[1:n]
z_atual <- atual$stage3$z.smoothed[1:n]
anos_s3 <- orig$dados$ano[9:(n+8)]
tri_s3 <- orig$dados$trimestre[9:(n+8)]
datas_s3 <- as.Date(paste0(anos_s3, "-", (tri_s3-1)*3+1, "-01"))

plot(datas_s3, z_orig, type = "l", col = "blue", lwd = 2, lty = 2,
     main = "z (premio de risco e outros fatores)",
     ylab = "% a.a.", xlab = "",
     ylim = range(c(z_orig, z_atual)))
lines(datas_s3, z_atual, col = "red", lwd = 2)
abline(v = data_break, col = "black", lwd = 1, lty = 3)

# Destacar a area inflada
polygon(c(datas_s3, rev(datas_s3)),
        c(pmax(z_orig, z_atual), rev(pmin(z_orig, z_atual))),
        col = rgb(0.5, 0, 0.5, 0.15), border = NA)

# Anotacao
idx_max_diff <- which.max(z_orig - z_atual)
arrows(datas_s3[idx_max_diff], z_orig[idx_max_diff],
       datas_s3[idx_max_diff], z_atual[idx_max_diff],
       code = 3, length = 0.08, col = "purple", lwd = 1.5)
text(datas_s3[idx_max_diff] + 200, (z_orig[idx_max_diff] + z_atual[idx_max_diff])/2,
     paste0(round(z_orig[idx_max_diff] - z_atual[idx_max_diff], 1), " pp"),
     col = "purple", cex = 0.85, font = 2)

legend("topright",
       c("z sem break (inflado)", "z com break"),
       col = c("blue", "red"), lwd = 2, lty = c(2, 1), bty = "n", cex = 0.85)

# --------------------------------------------------------------------------
# PAINEL 3: r* com vs sem break
# --------------------------------------------------------------------------

r_orig <- orig$stage3$rstar.smoothed[1:n]
r_atual <- atual$stage3$rstar.smoothed[1:n]

plot(datas_s3, r_orig, type = "l", col = "blue", lwd = 2, lty = 2,
     main = "r* (taxa natural de juros)",
     ylab = "% a.a.", xlab = "",
     ylim = range(c(r_orig, r_atual)))
lines(datas_s3, r_atual, col = "red", lwd = 2)
abline(v = data_break, col = "black", lwd = 1, lty = 3)

# Destacar periodo mais afetado
rect(as.Date("2009-01-01"), par("usr")[3],
     as.Date("2013-12-01"), par("usr")[4],
     col = rgb(1, 0.8, 0, 0.15), border = NA)
text(as.Date("2011-01-01"), max(r_orig)*0.95,
     paste0("Dif media\n", round(mean(r_atual[anos_s3>=2009 & anos_s3<=2013]) -
            mean(r_orig[anos_s3>=2009 & anos_s3<=2013]), 1), " pp"),
     col = "darkred", cex = 0.8, font = 2)

legend("topright",
       c("r* sem break", "r* com break"),
       col = c("blue", "red"), lwd = 2, lty = c(2, 1), bty = "n", cex = 0.85)

dev.off()

cat(paste0("\nGrafico salvo: outputs/graficos/break_2008_explicativo_", data.arquivo, ".png\n"))
