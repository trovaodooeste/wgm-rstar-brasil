png("outputs/graficos/leitao_matrix.png", width = 700, height = 600, res = 120)
par(mar = c(5, 5, 4, 2))

# Eixos
plot(NULL, xlim = c(0, 10), ylim = c(0, 10),
     xlab = "Doida", ylab = "Gostosa",
     main = "Matriz Gostosa-Doida: Escolhas do Leitao",
     cex.lab = 1.3, cex.main = 1.4)

# Zonas
rect(0, 0, 5, 5, col = rgb(1, 0, 0, 0.15), border = NA)  # Zona de perigo
rect(5, 5, 10, 10, col = rgb(0, 1, 0, 0.15), border = NA)  # Zona segura
rect(5, 0, 10, 5, col = rgb(1, 0.5, 0, 0.15), border = NA)  # Zona cuidado

# Linha diagonal (pessoa normal)
lines(c(4, 10), c(4, 10), lty = 2, col = "gray40", lwd = 2)
text(8.5, 7.8, "Pessoa\nnormal", col = "gray40", cex = 0.8)

# Linha do Leitao (quase vertical - aceita qualquer nivel de doida)
lines(c(1, 2), c(0, 10), col = "blue", lwd = 5)

# Seta e texto LEITAO
text(4, 8.5, "LEITAO", col = "blue", font = 2, cex = 1.8)
arrows(3.5, 8, 2.3, 7, col = "blue", lwd = 2, length = 0.15)

# Texto explicativo
text(1.5, 1, "Aceita\nqualquer\ncoisa", col = "blue", cex = 0.7, font = 3)

# Texto zonas
text(2.5, 3, "Perigo", col = "red", font = 2, cex = 1.1)
text(7.5, 7.5, "Zona\nSegura", col = "darkgreen", font = 2, cex = 1.1)
text(7.5, 2.5, "Cuidado", col = "darkorange", font = 2, cex = 1.1)

# Linha minima de gostosura (5)
abline(h = 5, lty = 3, col = "gray60")
text(0.8, 5.3, "min", col = "gray50", cex = 0.7)

dev.off()
cat("Salvo: outputs/graficos/leitao_matrix.png\n")
