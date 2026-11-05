# ==============================================================================
# Download da Curva DI-Pre (taxa prefixada) da B3
# Adaptado de: Werley Cordeiro (https://github.com/werleycordeiro/Brazil_Yield_Curve)
#
# Fonte: B3 (Brasil, Bolsa, Balcao) - Taxas Referenciais BM&F
# URL: https://www2.bmf.com.br/pages/portal/bmfbovespa/lumis/lum-taxas-referenciais-bmf-ptBR.asp
# Dados disponiveis desde: 08/08/2003
#
# Metodologia: A B3 calcula diariamente taxas referenciais para a curva PRE
# (prefixada) com base nos precos de ajuste dos contratos futuros de DI1.
# Os vertices sao interpolados via cubic spline para maturidades padrao.
#
# Maturidade de interesse: M12 (12 meses ~ 252 du) = taxa zero-cupom
# prefixada de 1 ano. Esta taxa eh FORWARD-LOOKING e reflete as
# expectativas do mercado para a trajetoria futura da Selic/CDI.
#
# Frequencia de download: ultimo dia util de cada mes
# Data do download: 2026-02-10
# ==============================================================================

cat("=== Download Curva DI-Pre B3 (mensal) ===\n\n")

# Pacotes necessarios
packages <- c("rvest", "httr", "functional")
new.packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) install.packages(new.packages, repos = "https://cloud.r-project.org")
suppressMessages(library(rvest))
suppressMessages(library(httr))
suppressMessages(library(functional))

# ---------------------------------------------------------------
# Gerar ultimo dia util de cada mes (aprox - pula fds)
# ---------------------------------------------------------------
gerar_fim_de_mes <- function(start_year, start_month, end_year, end_month) {
    dates <- c()
    y <- start_year
    m <- start_month
    while (y < end_year || (y == end_year && m <= end_month)) {
        # Ultimo dia do mes
        if (m == 12) {
            last_day <- as.Date(paste0(y + 1, "-01-01")) - 1
        } else {
            last_day <- as.Date(paste0(y, "-", sprintf("%02d", m + 1), "-01")) - 1
        }
        # Voltar para dia util (seg-sex)
        while (weekdays(last_day) %in% c("Saturday", "Sunday",
                                           "sábado", "domingo")) {
            last_day <- last_day - 1
        }
        dates <- c(dates, as.character(last_day))
        m <- m + 1
        if (m > 12) { m <- 1; y <- y + 1 }
    }
    return(as.Date(dates))
}

# ---------------------------------------------------------------
# Funcao para baixar curva de uma data
# ---------------------------------------------------------------
get_di_pre <- function(date_obj, maturities = c(1, 3, 6, 9, 12, 15, 18, 21, 24)) {
    date_str <- format(date_obj, "%d-%m-%Y")

    tryCatch({
        di <- GET(
            url = "https://www2.bmf.com.br/pages/portal/bmfbovespa/lumis/lum-taxas-referenciais-bmf-ptBR.asp",
            query = list(Data = date_str),
            timeout(30)
        )

        data <- read_html(di) %>%
            html_nodes("table") %>%
            html_nodes("td") %>%
            html_text()

        if (length(data) == 0) return(NULL)

        data <- data.frame(matrix(data, ncol = 3, byrow = TRUE),
                          stringsAsFactors = FALSE)
        data[, 2] <- as.numeric(gsub(",", ".", gsub("\\.", "", data[, 2])))
        data[, 3] <- as.numeric(gsub(",", ".", gsub("\\.", "", data[, 3])))

        # Converter dias corridos para meses (~21 du por mes)
        t <- as.integer(as.matrix(data[, 1])) / 21
        y <- as.numeric(as.matrix(data[, 2]))

        # Interpolar via spline
        spl <- smooth.spline(y ~ t)
        new <- predict(spl, maturities)

        result <- new$y
        names(result) <- paste0("M", maturities)
        return(result)

    }, error = function(e) {
        return(NULL)
    })
}

# ---------------------------------------------------------------
# MAIN: Baixar serie historica mensal
# ---------------------------------------------------------------

# Maturidades em meses (12 = 1 ano, que eh nosso alvo principal)
maturities <- c(1, 3, 6, 9, 12, 15, 18, 21, 24)

# Periodo: ago/2003 (inicio dos dados) a jan/2026
dates <- gerar_fim_de_mes(2003, 8, 2026, 1)
cat("Total de datas a consultar:", length(dates), "\n")
cat("Periodo:", format(min(dates), "%d/%m/%Y"), "a", format(max(dates), "%d/%m/%Y"), "\n\n")

# Matriz de resultados
results <- data.frame(
    data = as.Date(character()),
    ano = integer(),
    mes = integer(),
    trimestre = integer(),
    stringsAsFactors = FALSE
)
for (m in paste0("M", maturities)) {
    results[[m]] <- numeric()
}

# Download
n_ok <- 0
n_fail <- 0

for (i in seq_along(dates)) {
    d <- dates[i]
    date_str <- format(d, "%d/%m/%Y")

    curve <- get_di_pre(d, maturities)

    if (!is.null(curve)) {
        row <- data.frame(
            data = d,
            ano = as.integer(format(d, "%Y")),
            mes = as.integer(format(d, "%m")),
            trimestre = ceiling(as.integer(format(d, "%m")) / 3),
            stringsAsFactors = FALSE
        )
        for (m in names(curve)) {
            row[[m]] <- curve[m]
        }
        results <- rbind(results, row)
        n_ok <- n_ok + 1
        cat(sprintf("  [%3d/%d] %s: M12=%.2f%%\n", i, length(dates), date_str, curve["M12"]))
    } else {
        n_fail <- n_fail + 1
        cat(sprintf("  [%3d/%d] %s: SEM DADOS\n", i, length(dates), date_str))

        # Tentar dias anteriores (feriado?)
        for (offset in 1:5) {
            d2 <- d - offset
            if (weekdays(d2) %in% c("Saturday", "Sunday", "sábado", "domingo")) next
            curve2 <- get_di_pre(d2, maturities)
            if (!is.null(curve2)) {
                row <- data.frame(
                    data = d2,
                    ano = as.integer(format(d, "%Y")),  # manter ano/mes original
                    mes = as.integer(format(d, "%m")),
                    trimestre = ceiling(as.integer(format(d, "%m")) / 3),
                    stringsAsFactors = FALSE
                )
                for (m in names(curve2)) {
                    row[[m]] <- curve2[m]
                }
                results <- rbind(results, row)
                n_ok <- n_ok + 1
                n_fail <- n_fail - 1
                cat(sprintf("         -> fallback %s: M12=%.2f%%\n",
                            format(d2, "%d/%m/%Y"), curve2["M12"]))
                break
            }
        }
    }

    # Rate limiting
    Sys.sleep(0.5)
}

cat(sprintf("\nDownload concluido: %d OK, %d falhas\n\n", n_ok, n_fail))

# ---------------------------------------------------------------
# Salvar resultados
# ---------------------------------------------------------------

# CSV mensal (todas as maturidades)
output_dir <- "data/external"
write.csv(results, file.path(output_dir, "di_pre_b3_mensal.csv"), row.names = FALSE)
cat("Salvo:", file.path(output_dir, "di_pre_b3_mensal.csv"), "\n")

# CSV trimestral (media dos meses do trimestre, foco em M12)
tri <- aggregate(cbind(M1, M3, M6, M9, M12, M15, M18, M21, M24) ~ ano + trimestre,
                 data = results, FUN = mean)
tri <- tri[order(tri$ano, tri$trimestre), ]
write.csv(tri, file.path(output_dir, "di_pre_b3_trimestral.csv"), row.names = FALSE)
cat("Salvo:", file.path(output_dir, "di_pre_b3_trimestral.csv"), "\n")

# Mostrar serie trimestral M12
cat("\n=== Serie Trimestral DI-Pre M12 (% a.a., base 252) ===\n")
for (i in 1:nrow(tri)) {
    cat(sprintf("  %d-Q%d: %.2f%%\n", tri$ano[i], tri$trimestre[i], tri$M12[i]))
}

# Metadados
meta <- list(
    fonte = "B3 (Brasil, Bolsa, Balcao) - Taxas Referenciais BM&F PRE",
    url = "https://www2.bmf.com.br/pages/portal/bmfbovespa/lumis/lum-taxas-referenciais-bmf-ptBR.asp",
    script_original = "https://github.com/werleycordeiro/Brazil_Yield_Curve/blob/master/DI-Pre.R",
    autor_original = "Werley Cordeiro",
    descricao = "Curva de juros prefixada (DI x Pre) interpolada via cubic spline",
    vertice_principal = "M12 (12 meses ~ 252 du) - taxa zero-cupom prefixada de 1 ano",
    unidade = "% a.a. (base 252 dias uteis)",
    frequencia_download = "Ultimo dia util de cada mes",
    data_download = Sys.Date(),
    periodo = paste(format(min(results$data), "%d/%m/%Y"), "a", format(max(results$data), "%d/%m/%Y")),
    n_meses = nrow(results),
    n_trimestres = nrow(tri),
    maturidades_meses = maturities
)
saveRDS(meta, file.path(output_dir, "di_pre_b3_meta.rds"))
cat("\nSalvo:", file.path(output_dir, "di_pre_b3_meta.rds"), "\n")

cat("\n=== Concluido ===\n")
