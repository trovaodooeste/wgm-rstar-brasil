# Taxa Natural de Juros (r*) no Brasil — Modelo WGM com Volatilidade Variável (κ)

Código e dados de reprodução para a estimação da taxa natural de juros (r*) do Brasil
usando uma adaptação do modelo semiestrutural de Laubach-Williams / Holston-Laubach-Williams
(HLW) para a economia brasileira, com duas extensões:

- **Break estrutural em 2008-Q4**, para capturar a mudança de tendência no PIB potencial
  após a crise financeira global;
- **κ (kappa) de volatilidade variável no tempo** durante a pandemia de COVID-19
  (2020-2022), seguindo a metodologia de Holston, Laubach & Williams (2023).

O modelo é estimado em duas variantes de taxa nominal: **Selic** e **DI-Pre B3 (M12)**.

## Metodologia (resumo)

O modelo é um sistema de espaço de estados estimado por máxima verossimilhança via
filtro de Kalman, em 3 estágios (Stage 1/2/3), seguindo Laubach & Williams (2003) e
Holston, Laubach & Williams (2017, 2023), com as adaptações de Moreira & Portugal (2019)
para a economia brasileira (restrições em `b_y` e `a_r`, prêmio de risco `z₀`).

- **Curva IS**: hiato do produto em função de seus próprios lags e do gap de juros reais
  (taxa real menos taxa natural), com efeito do índice de rigor (*stringency*) da Oxford
  (OxCGRT) durante a pandemia.
- **Curva de Phillips**: inflação em função de seus lags, médias móveis de inflação
  passada e do hiato do produto defasado.
- **r\* = c·g + z**: a taxa natural é decomposta em crescimento potencial (`g`) e
  outros fatores / prêmio de risco (`z`).
- **κₜ²** entra multiplicando a matriz de covariância dos erros de medição (R) durante
  2020-2022, permitindo que a variância dos choques aumente temporariamente sem
  contaminar a tendência de longo prazo.
- **Intervalos de confiança**: procedimento de Hamilton (1986), com duas variantes
  incluídas — apenas incerteza paramétrica (perturbação da matriz de informação) e a
  versão completa (paramétrica + incerteza do estado via RTS smoother).

## Estrutura do repositório

```
├── data/
│   ├── raw/                         # Série principal e validação (PIB, Selic, IPCA)
│   └── external/                    # Expectativas Focus, Oxford Stringency Index, DI-Pre B3
├── scripts/
│   ├── median.unbiased.estimator.stage1.R   # Stock & Watson (1998) - correção de viés de sigma_ystar
│   ├── median.unbiased.estimator.stage2.R   # Stock & Watson (1998) - correção de viés de lambda_g
│   ├── kalman.*.R, leitao_matrix.R,
│   │   log.likelihood.wrapper.R             # Infraestrutura do filtro de Kalman (base HLW)
│   ├── kalman_smoother.R                    # RTS smoother (para o IC completo de Hamilton)
│   ├── rstar_stage1_com_kappa.R             # Stage 1: hiato do produto, break 2008, kappa
│   ├── rstar_stage2_com_kappa.R             # Stage 2: taxa real/output gap, a_r
│   ├── rstar_stage3_com_kappa.R             # Stage 3: r* = c*g + z, Kalman completo
│   ├── run_wgm_kappa_params_fixos.R         # Pipeline principal — modelo Selic
│   ├── run_wgm_dipre_params_fixos.R         # Pipeline principal — modelo DI-Pre B3
│   ├── gerar_resultados_orientador_ic.R     # IC 95% paramétrico (Hamilton) — Selic
│   ├── gerar_ic_dipre.R                     # IC 95% paramétrico (Hamilton) — DI-Pre
│   ├── gerar_ic_hamilton_completo.R         # IC 95% completo (paramétrico + estado) — Selic
│   ├── gerar_ic70_dashboard.R               # Reescala o IC para 70% (usado no texto final)
│   ├── gerar_graficos_individuais.R         # Figuras por modelo (r*, g, z, hiato, stance, kappa)
│   ├── gerar_graficos_comparacao_3modelos.R # Comparação: original x Selic+κ x DI-Pre+κ
│   ├── grafico_break_2008_explicativo.R     # Figura explicativa do break estrutural de 2008
│   └── download_di_pre_mensal.R             # Scraping da curva DI-Pre na B3 (opcional/provenance)
└── outputs/
    ├── results/   # .rds pré-computados (estimação completa e ICs), para inspeção sem re-rodar
    └── graficos/  # Saída das figuras (gerada ao rodar os scripts)
```

## Requisitos

- **R ≥ 4.1** (desenvolvido e testado em R 4.5.1)
- Pacotes: `readxl`, `nloptr`, `numDeriv`, `mFilter`, `MASS`
- Apenas para `download_di_pre_mensal.R` (opcional, requer internet): `rvest`, `httr`, `functional`

```r
install.packages(c("readxl", "nloptr", "numDeriv", "mFilter",
                    "rvest", "httr", "functional"))
```

## Como reproduzir

Todos os scripts devem ser executados **a partir da raiz do repositório** (eles usam
caminhos relativos). Os `.rds` finais já estão incluídos em `outputs/results/`, então os
passos 3 a 9 funcionam imediatamente sem precisar re-rodar os passos 1 e 2.

```bash
Rscript scripts/run_wgm_kappa_params_fixos.R      # 1. Modelo Selic (break 2008 + kappa)
Rscript scripts/run_wgm_dipre_params_fixos.R       # 2. Modelo DI-Pre B3 (idem)
Rscript scripts/gerar_resultados_orientador_ic.R   # 3. IC 95% parametrico - Selic
Rscript scripts/gerar_ic_dipre.R                   # 4. IC 95% parametrico - DI-Pre
Rscript scripts/gerar_ic_hamilton_completo.R       # 5. IC 95% completo (parametrico + estado)
Rscript scripts/gerar_ic70_dashboard.R             # 6. Reescala para IC 70%
Rscript scripts/gerar_graficos_individuais.R       # 7. Figuras por modelo
Rscript scripts/gerar_graficos_comparacao_3modelos.R # 8. Comparacao dos 3 modelos
Rscript scripts/grafico_break_2008_explicativo.R   # 9. Figura do break 2008
```

**Nota sobre nomes de arquivo:** os scripts 3-9 referenciam os `.rds` dos passos 1-2 por
nome de arquivo com data fixa (ex.: `wgm_kappa_params_fixos_20260212.rds`), que é o
arquivo incluído neste repositório. Se você rodar os passos 1-2 novamente, o resultado é
salvo com a data do dia (`outputs/results/wgm_kappa_params_fixos_YYYYMMDD.rds`); para
encadear com os passos seguintes, copie/renomeie esse arquivo para o nome esperado ou
ajuste a data dentro do script correspondente.

Todos os 9 scripts foram validados rodando do zero neste repositório e reproduzem
exatamente os números reportados abaixo.

## Resultados principais

| Métrica | Selic (break 2008 + κ) | DI-Pre B3 (break 2008 + κ) |
|---|:---:|:---:|
| r* atual (2025-T2) | 4.61% [IC 95%: 1.77%, 7.46%] | 4.83% [IC 95%: 1.78%, 7.88%] |
| r* médio (amostra) | 5.03% [IC 95%: 3.21%, 6.84%] | 5.12% [IC 95%: 3.03%, 7.22%] |
| g (crescimento potencial) | 2.10% | 1.96% |
| z (outros fatores) | 2.93% | 3.16% |
| κ_2020 (2020-Q2 a Q4) | 3.28 | 3.29 |
| κ_2021 | 1.13 | 1.11 |
| κ_2022 | 1.50 | 1.56 |

Os intervalos de confiança acima usam apenas a incerteza paramétrica (procedimento de
Hamilton, 1986, caminho A). A versão completa, que soma a incerteza do estado via RTS
smoother, produz bandas bem mais largas (ver `gerar_ic_hamilton_completo.R`) porque o
componente `z` é fracamente identificado no modelo (entra apenas via `a_r`, de magnitude
pequena, na curva IS) — um caso do *pile-up problem* discutido em Stock & Watson (1998).

## Fontes de dados

- **PIB real trimestral, Selic, IPCA**: Banco Central do Brasil (SGS) / IBGE.
- **Expectativas de inflação (Focus)**: Banco Central do Brasil.
- **Oxford COVID-19 Government Response Tracker (stringency index)**: Oxford/OxCGRT.
- **DI-Pre B3 (taxa prefixada, vértice M12)**: B3 (Brasil, Bolsa, Balcão), obtida via
  `scripts/download_di_pre_mensal.R` (adaptado de
  [Werley Cordeiro, Brazil_Yield_Curve](https://github.com/werleycordeiro/Brazil_Yield_Curve)).

## Referências

- Laubach, T.; Williams, J. C. (2003). *Measuring the Natural Rate of Interest*.
  Review of Economics and Statistics.
- Holston, K.; Laubach, T.; Williams, J. C. (2017). *Measuring the Natural Rate of
  Interest: International Trends and Determinants*. Journal of International Economics.
- Holston, K.; Laubach, T.; Williams, J. C. (2023). *Measuring the Natural Rate of
  Interest after COVID-19*. FRBNY Staff Reports.
- Moreira, R. R.; Portugal, M. S. (2019). *Estimating the Natural Rate of Interest for
  Brazil*.
- Stock, J. H.; Watson, M. W. (1998). *Median Unbiased Estimation of Coefficient
  Variance in a Time-Varying Parameter Model*. Journal of the American Statistical
  Association.
- Santos, F.; Muinhos, M. (2024). Estimação da taxa natural de juros brasileira via
  HLW (2023), com swaps DI-Pre do B3.

## Licença

Este repositório é distribuído sob a licença MIT (ver [LICENSE](LICENSE)). Os dados de
terceiros (BCB, Oxford OxCGRT, B3) seguem os termos de uso de suas respectivas fontes.
