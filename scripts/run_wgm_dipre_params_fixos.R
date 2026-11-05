# ==============================================================================
# MODELO WGM BRASIL - DI-PRE + KAPPA + BREAK 2008 COM PARAMETROS FIXOS
# ==============================================================================
# IDENTICO ao run_wgm_kappa_params_fixos.R (Selic) exceto:
#   - Taxa nominal: DI-Pre B3 M12 em vez de Selic
#   - Taxa real: (DI_Pre_M12 - E[IPCA_12m]) / 4
#
# USA OS MESMOS SCRIPTS MAIN (nao exante_di!):
#   - scripts/rstar_stage1_com_kappa.R
#   - scripts/rstar_stage2_com_kappa.R
#   - scripts/rstar_stage3_com_kappa.R
#
# PARAMETROS FIXOS (do modelo original Selic sem kappa):
#   - a_1 = 0.5778
#   - a_2 = 0.2506
#   - phi = -9.0629
#   - sigma_IS = 1.4267
# ==============================================================================

library(readxl)

cat("\n================================================================\n")
cat("WGM DI-PRE - PARAMS FIXOS (SCRIPTS MAIN) + BREAK 2008 + KAPPA\n")
cat("================================================================\n\n")

# Carregar funcoes - MESMOS SCRIPTS MAIN do modelo Selic params fixos
# (assume que o script e rodado a partir da raiz do repositorio)
wgm_path <- "."
source(file.path(wgm_path, "scripts/rstar_stage1_com_kappa.R"))
source(file.path(wgm_path, "scripts/rstar_stage2_com_kappa.R"))
source(file.path(wgm_path, "scripts/rstar_stage3_com_kappa.R"))
source(file.path(wgm_path, "scripts/median.unbiased.estimator.stage1.R"))
source(file.path(wgm_path, "scripts/median.unbiased.estimator.stage2.R"))

# ==============================================================================
# CARREGAR DADOS (identico ao Selic params fixos)
# ==============================================================================

cat("Carregando dados...\n")

df <- read_excel(file.path(wgm_path, "data/raw/dados_WGM_brasil_2000_2025.xlsx"))

# Expectativas Focus
exp <- read.csv(file.path(wgm_path, "data/external/expectativas_inflacao_bcb.csv"))
exp$date <- as.Date(exp$date)
exp$ano <- as.numeric(format(exp$date, "%Y"))
exp$trimestre <- ceiling(as.numeric(format(exp$date, "%m")) / 3)
exp_tri <- aggregate(expectativa_ipca_12m ~ ano + trimestre, data = exp, FUN = mean)
df <- merge(df, exp_tri, by = c("ano", "trimestre"), all.x = TRUE)
df <- df[order(df$ano, df$trimestre), ]

# Stringency Index
stringency_raw <- read.csv(file.path(wgm_path, "data/external/oxford_stringency_index.csv"))
stringency_raw$date <- as.Date(stringency_raw$date)
stringency_raw$ano <- as.numeric(format(stringency_raw$date, "%Y"))
stringency_raw$trimestre <- ceiling(as.numeric(format(stringency_raw$date, "%m")) / 3)
stringency_tri <- aggregate(stringency_index ~ ano + trimestre, data = stringency_raw, FUN = max)
stringency_tri$stringency <- stringency_tri$stringency_index / 100
df <- merge(df, stringency_tri[, c("ano", "trimestre", "stringency")], by = c("ano", "trimestre"), all.x = TRUE)
df$stringency[is.na(df$stringency)] <- 0
df <- df[order(df$ano, df$trimestre), ]

# IPCA trimestral
df_val <- read_excel(file.path(wgm_path, "data/raw/dados_completos_2000_2025_VALIDADO.xlsx"))
df$ipca_tri <- df_val$ipca_tri[match(df$data, df_val$data)]

# ==============================================================================
# CARREGAR DI-PRE B3 (M12) - UNICA DIFERENCA DO SCRIPT SELIC
# ==============================================================================

cat("Carregando DI-Pre B3 (M12)...\n")

di_pre_tri <- read.csv(file.path(wgm_path, "data/external/di_pre_b3_trimestral.csv"))
df <- merge(df, di_pre_tri[, c("ano", "trimestre", "M12")],
            by = c("ano", "trimestre"), all.x = TRUE)
df <- df[order(df$ano, df$trimestre), ]
names(df)[names(df) == "M12"] <- "di_pre_m12"

# ==============================================================================
# TAXA REAL: DI-Pre em vez de Selic (UNICA MUDANCA)
# ==============================================================================

# Taxa real ex-ante TRIMESTRAL (consistente com modelo original que gerou os params fixos)
# O modelo original (run_wgm_exante_FINAL_z0_3.R) usa /4
df$taxa_real_exante_tri <- (df$di_pre_m12 - df$expectativa_ipca_12m) / 4

# Informativo
cat("\nComparacao DI-Pre M12 vs Selic:\n")
valid <- !is.na(df$di_pre_m12) & !is.na(df$selic_aa)
cat("  DI-Pre M12 medio:", round(mean(df$di_pre_m12[valid]), 2), "% aa\n")
cat("  Selic media:     ", round(mean(df$selic_aa[valid]), 2), "% aa\n")
cat("  Spread medio:    ", round(mean(df$di_pre_m12[valid] - df$selic_aa[valid]), 2), "pp\n")
cat("  Correlacao:      ", round(cor(df$di_pre_m12[valid], df$selic_aa[valid]), 4), "\n\n")

# Filtrar
df <- df[!is.na(df$pib_dessaz), ]
df <- df[!is.na(df$taxa_real_exante_tri), ]
df <- df[!is.na(df$di_pre_m12), ]
df <- df[df$ano >= 2003, ]

cat("Periodo:", min(df$ano), "-", max(df$ano), "\n")
cat("Trimestres:", nrow(df), "\n\n")

# ==============================================================================
# PARAMETROS FIXOS (IDENTICOS ao Selic params fixos)
# ==============================================================================

cat("================================================================\n")
cat("PARAMETROS FIXOS DO MODELO ORIGINAL\n")
cat("================================================================\n")

ORIG_A1 <- 0.5778
ORIG_A2 <- 0.2506
ORIG_PHI <- -9.0629
ORIG_SIGMA_IS <- 1.4267
ORIG_AR <- -0.0956

cat("  a_1 =", ORIG_A1, "(FIXO)\n")
cat("  a_2 =", ORIG_A2, "(FIXO)\n")
cat("  phi =", ORIG_PHI, "(FIXO)\n")
cat("  sigma_IS =", ORIG_SIGMA_IS, "(FIXO)\n")
cat("  a_r =", ORIG_AR, "(FIXO)\n")
cat("  Taxa nominal: DI-Pre B3 M12 (em vez de Selic)\n")
cat("================================================================\n\n")

lambda.g <- 0.15
lambda.z <- 0.05       # usado apenas para valor inicial de sigma_z
z.initial <- 2.2       # Moreira & Portugal (2019)
sigma.z.upper <- 2.2   # sigma_z livre por MLE com teto

# ==============================================================================
# STAGE 1 - COM PARAMETROS FIXOS (IDENTICO ao Selic params fixos)
# ==============================================================================

cat("=== STAGE 1 (PARAMS FIXOS + BREAK 2008 + KAPPA COVID) ===\n")

stage1_kappa <- rstar.stage1.com.kappa(
  log.output = as.numeric(df$ln_pib),
  inflation = as.numeric(df$ipca_tri),
  covid.dummy = as.numeric(df$stringency),
  anos = as.numeric(df$ano),
  trimestres = as.numeric(df$trimestre),
  b.y.constraint = 0.25,
  sigma.ystar.constraint = 0.5,
  a.sum.constraint = NA,
  fix.a1 = ORIG_A1,
  fix.a2 = ORIG_A2,
  fix.phi = ORIG_PHI,
  fix.sigma.is = ORIG_SIGMA_IS,
  use.break.2008 = TRUE,
  use.break.2020 = FALSE,
  data.start.year = 2003,
  data.start.quarter = 1,
  estimar.kappa = TRUE,
  kappa.2020.init = 5.0,
  kappa.2021.init = 2.0,
  kappa.2022.init = 1.5
)

# Compatibilizar estrutura para Stage 2
stage1_compat <- list(
  theta = c(
    stage1_kappa$theta[stage1_kappa$param.num["a_1"]],
    stage1_kappa$theta[stage1_kappa$param.num["a_2"]],
    stage1_kappa$theta[stage1_kappa$param.num["b_pi1"]],
    stage1_kappa$theta[stage1_kappa$param.num["b_pi2"]],
    stage1_kappa$theta[stage1_kappa$param.num["b_y"]],
    stage1_kappa$theta[stage1_kappa$param.num["g"]],
    stage1_kappa$theta[stage1_kappa$param.num["sigma_is"]],
    stage1_kappa$theta[stage1_kappa$param.num["sigma_pi"]],
    stage1_kappa$theta[stage1_kappa$param.num["sigma_ystar"]]
  ),
  output.gap.smoothed = stage1_kappa$output.gap.smoothed,
  potential.smoothed = stage1_kappa$potential.smoothed,
  log.likelihood = stage1_kappa$log.likelihood,
  phi = stage1_kappa$theta[stage1_kappa$param.num["phi"]],
  estimar.kappa = stage1_kappa$estimar.kappa,
  kappas = stage1_kappa$kappas
)

# ==============================================================================
# STAGE 2 (MESMOS SCRIPTS MAIN, mesmos constraints)
# ==============================================================================

cat("\n=== STAGE 2 ===\n")
stage2 <- rstar.stage2.com.kappa(
  log.output = as.numeric(df$ln_pib),
  inflation = as.numeric(df$ipca_tri),
  real.interest.rate = as.numeric(df$taxa_real_exante_tri),
  stage1.results = stage1_compat,
  anos = as.numeric(df$ano),
  trimestres = as.numeric(df$trimestre),
  b.y.constraint = 0.25,
  a.r.constraint = ORIG_AR,     # FIXAR a_r no valor original
  a.r.lower = ORIG_AR,          # FIXAR a_r no valor original (-0.0956)
  sigma.ystar.constraint = 0.5,
  sigma.ystar.lower = 0.05,    # Floor sigma_y* = 0.05 (evitar colapso)
  lambda.g = lambda.g
)

# ==============================================================================
# STAGE 3
# ==============================================================================

cat("\n=== STAGE 3 ===\n")
stage3 <- rstar.stage3.com.kappa(
  log.output = as.numeric(df$ln_pib),
  inflation = as.numeric(df$ipca_tri),
  real.interest.rate = as.numeric(df$taxa_real_exante_tri),
  stage2.results = stage2,
  anos = as.numeric(df$ano),
  trimestres = as.numeric(df$trimestre),
  c.param = 1,
  lambda.g = lambda.g,
  lambda.z = lambda.z,
  z.initial = z.initial,
  sigma.ystar.constraint = 0.5,
  sigma.ystar.lower = 0.05,    # Floor sigma_y* = 0.05
  sigma.z.upper = sigma.z.upper
)

# ==============================================================================
# RESULTADOS E COMPARACAO
# ==============================================================================

cat("\n================================================================\n")
cat("RESULTADOS - DI-PRE (SCRIPTS MAIN, PARAMS FIXOS)\n")
cat("================================================================\n\n")

# Carregar original Selic para comparar
original <- readRDS(file.path(wgm_path, "outputs/results/wgm_exante_FINAL_z0_3_20250122.rds"))

cat("Parametros Stage 1:\n")
cat("  a_1 =", round(stage1_kappa$theta[stage1_kappa$param.num["a_1"]], 4), "(fixo:", ORIG_A1, ")\n")
cat("  a_2 =", round(stage1_kappa$theta[stage1_kappa$param.num["a_2"]], 4), "(fixo:", ORIG_A2, ")\n")
cat("  phi =", round(stage1_kappa$theta[stage1_kappa$param.num["phi"]], 4), "(fixo:", ORIG_PHI, ")\n")
cat("  sigma_IS =", round(stage1_kappa$theta[stage1_kappa$param.num["sigma_is"]], 4), "(fixo:", ORIG_SIGMA_IS, ")\n\n")

kappas_est <- stage1_kappa$kappas
cat("Kappas estimados:\n")
cat("  kappa_2020 =", round(kappas_est["kappa_2020"], 2), "\n")
cat("  kappa_2021 =", round(kappas_est["kappa_2021"], 2), "\n")
cat("  kappa_2022 =", round(kappas_est["kappa_2022"], 2), "\n\n")

cat("Stage 2:\n")
cat("  a_r =", round(stage2$theta[stage2$param.num["a_r"]], 4), "\n\n")

cat("Stage 3 (r* = g + z):\n")
cat("  g medio =", round(mean(stage3$g.smoothed) * 100, 2), "% aa\n")
cat("  z medio =", round(mean(stage3$z.smoothed), 2), "% aa\n")
cat("  r* MEDIO =", round(mean(stage3$rstar.smoothed), 2), "% aa\n")
cat("  r* ATUAL =", round(tail(stage3$rstar.smoothed, 1), 2), "% aa\n")
cat("  Gap medio =", round(mean(stage3$output.gap.smoothed), 2), "%\n")
cat("  Gap atual =", round(tail(stage3$output.gap.smoothed, 1), 2), "%\n\n")

cat("================================================================\n")
cat("COMPARACAO: SELIC orig | SELIC pf+kappa | DI-PRE pf+kappa\n")
cat("================================================================\n")

# Carregar Selic params fixos
pf_files <- list.files(file.path(wgm_path, "outputs/results"), pattern = "^wgm_kappa_params_fixos_", full.names = TRUE)
if (length(pf_files) > 0) {
  pf_selic <- readRDS(sort(pf_files, decreasing = TRUE)[1])
  cat("                    |  SELIC orig | SELIC pf+k | DI-PRE pf+k\n")
  cat("--------------------|-------------|------------|------------\n")
  cat(sprintf("r* medio            | %10.2f%% | %9.2f%% | %10.2f%%\n",
      mean(original$stage3$rstar.smoothed),
      mean(pf_selic$stage3$rstar.smoothed),
      mean(stage3$rstar.smoothed)))
  cat(sprintf("r* atual            | %10.2f%% | %9.2f%% | %10.2f%%\n",
      tail(original$stage3$rstar.smoothed, 1),
      tail(pf_selic$stage3$rstar.smoothed, 1),
      tail(stage3$rstar.smoothed, 1)))
  cat(sprintf("g medio             | %10.2f%% | %9.2f%% | %10.2f%%\n",
      mean(original$stage3$g.smoothed) * 100,
      mean(pf_selic$stage3$g.smoothed) * 100,
      mean(stage3$g.smoothed) * 100))
  cat(sprintf("z medio             | %10.2f%% | %9.2f%% | %10.2f%%\n",
      mean(original$stage3$z.smoothed),
      mean(pf_selic$stage3$z.smoothed),
      mean(stage3$z.smoothed)))
  cat(sprintf("Gap medio           | %10.2f%% | %9.2f%% | %10.2f%%\n",
      mean(original$stage3$output.gap.smoothed),
      mean(pf_selic$stage3$output.gap.smoothed),
      mean(stage3$output.gap.smoothed)))
  cat(sprintf("kappa_2020          |          - | %9.2f | %10.2f\n",
      pf_selic$kappas["kappa_2020"], kappas_est["kappa_2020"]))
  cat(sprintf("kappa_2021          |          - | %9.2f | %10.2f\n",
      pf_selic$kappas["kappa_2021"], kappas_est["kappa_2021"]))
  cat(sprintf("kappa_2022          |          - | %9.2f | %10.2f\n",
      pf_selic$kappas["kappa_2022"], kappas_est["kappa_2022"]))
} else {
  cat("Arquivo Selic params fixos nao encontrado.\n")
}

# ==============================================================================
# SALVAR
# ==============================================================================

cat("\n================================================================\n")
cat("SALVANDO RESULTADOS\n")
cat("================================================================\n")

data.arquivo <- format(Sys.Date(), "%Y%m%d")

resultado <- list(
  stage1 = stage1_kappa,
  stage2 = stage2,
  stage3 = stage3,
  lambda.g = lambda.g,
  lambda.z = lambda.z,
  z.initial = z.initial,
  kappas = stage1_kappa$kappas,
  params.fixos = list(
    a1 = ORIG_A1,
    a2 = ORIG_A2,
    phi = ORIG_PHI,
    sigma_is = ORIG_SIGMA_IS
  ),
  dados = df,
  modelo = list(
    taxa.nominal = "DI-Pre B3 M12",
    taxa.real = "r = (DI_Pre_M12 - E[IPCA_12m]) / 4",
    scripts = "MAIN (nao exante_di)",
    break.2008 = TRUE,
    kappa.covid = TRUE,
    params.fixos = TRUE
  ),
  resumo = list(
    rstar.medio = mean(stage3$rstar.smoothed),
    rstar.atual = tail(stage3$rstar.smoothed, 1),
    g.medio = mean(stage3$g.smoothed) * 100,
    z.medio = mean(stage3$z.smoothed),
    gap.medio = mean(stage3$output.gap.smoothed),
    a_r = stage2$theta[stage2$param.num["a_r"]],
    kappa_2020 = kappas_est["kappa_2020"],
    kappa_2021 = kappas_est["kappa_2021"],
    kappa_2022 = kappas_est["kappa_2022"]
  )
)

if (!dir.exists(file.path(wgm_path, "outputs/results"))) dir.create(file.path(wgm_path, "outputs/results"), recursive = TRUE)
if (!dir.exists(file.path(wgm_path, "outputs/graficos"))) dir.create(file.path(wgm_path, "outputs/graficos"), recursive = TRUE)

saveRDS(resultado, file.path(wgm_path, paste0("outputs/results/wgm_dipre_params_fixos_", data.arquivo, ".rds")))
cat(paste0("Salvo: outputs/results/wgm_dipre_params_fixos_", data.arquivo, ".rds\n"))

# Grafico
png(file.path(wgm_path, paste0("outputs/graficos/wgm_dipre_params_fixos_", data.arquivo, ".png")),
    width = 1400, height = 800, res = 120)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))

# r*
plot(original$stage3$rstar.smoothed, type = "l", col = "blue", lwd = 2,
     main = "r* Comparacao", ylab = "% aa", xlab = "Trimestre",
     ylim = range(c(original$stage3$rstar.smoothed, stage3$rstar.smoothed)))
lines(stage3$rstar.smoothed, col = "red", lwd = 2)
if (length(pf_files) > 0) lines(pf_selic$stage3$rstar.smoothed, col = "darkgreen", lwd = 1, lty = 2)
legend("topright", c("Selic (orig)", "DI-Pre (pf)", "Selic (pf)"),
       col = c("blue", "red", "darkgreen"), lwd = c(2, 2, 1), lty = c(1, 1, 2), bty = "n", cex = 0.7)

# Gap
plot(original$stage3$output.gap.smoothed, type = "l", col = "blue", lwd = 2,
     main = "Output Gap", ylab = "%", xlab = "Trimestre",
     ylim = range(c(original$stage3$output.gap.smoothed, stage3$output.gap.smoothed)))
lines(stage3$output.gap.smoothed, col = "red", lwd = 2)
if (length(pf_files) > 0) lines(pf_selic$stage3$output.gap.smoothed, col = "darkgreen", lwd = 1, lty = 2)
abline(h = 0, col = "gray", lty = 2)
legend("topright", c("Selic (orig)", "DI-Pre (pf)", "Selic (pf)"),
       col = c("blue", "red", "darkgreen"), lwd = c(2, 2, 1), lty = c(1, 1, 2), bty = "n", cex = 0.7)

# g
plot(original$stage3$g.smoothed * 100, type = "l", col = "blue", lwd = 2,
     main = "g (tendencia)", ylab = "% aa", xlab = "Trimestre",
     ylim = range(c(original$stage3$g.smoothed * 100, stage3$g.smoothed * 100)))
lines(stage3$g.smoothed * 100, col = "red", lwd = 2)
legend("topright", c("Selic (orig)", "DI-Pre (pf)"),
       col = c("blue", "red"), lwd = 2, bty = "n", cex = 0.7)

# z
plot(original$stage3$z.smoothed, type = "l", col = "blue", lwd = 2,
     main = "z (outros)", ylab = "% aa", xlab = "Trimestre",
     ylim = range(c(original$stage3$z.smoothed, stage3$z.smoothed)))
lines(stage3$z.smoothed, col = "red", lwd = 2)
legend("topright", c("Selic (orig)", "DI-Pre (pf)"),
       col = c("blue", "red"), lwd = 2, bty = "n", cex = 0.7)

# Kappa
kappa.serie <- stage1_kappa$kappa.vec
plot(kappa.serie, type = "l", col = "orange", lwd = 2,
     main = "Kappa COVID", ylab = "kappa", xlab = "Trimestre")
abline(h = 1, col = "gray", lty = 2)
legend("topright",
       c(paste("k_2020:", round(kappas_est["kappa_2020"], 2)),
         paste("k_2021:", round(kappas_est["kappa_2021"], 2)),
         paste("k_2022:", round(kappas_est["kappa_2022"], 2))),
       col = "orange", lwd = 2, bty = "n", cex = 0.7)

# Diferenca r*
diff_rstar <- stage3$rstar.smoothed - original$stage3$rstar.smoothed
plot(diff_rstar, type = "l", col = "purple", lwd = 2,
     main = "Dif r* (DI-Pre - Selic orig)", ylab = "pp", xlab = "Trimestre")
abline(h = 0, col = "gray", lty = 2)
abline(h = mean(diff_rstar), col = "red", lty = 2)

dev.off()
cat(paste0("Salvo: outputs/graficos/wgm_dipre_params_fixos_", data.arquivo, ".png\n"))

cat("\n================================================================\n")
cat("MODELO DI-PRE (SCRIPTS MAIN) PARAMS FIXOS FINALIZADO!\n")
cat("================================================================\n")
