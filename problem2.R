library(ggplot2)
library(gridExtra)

set.seed(42)

# --- Parametry ---
N_sim   <- 10000
n_sizes <- c(1, 5, 30, 100)
lambda  <- 1

mu_X   <- 1 / lambda
sig_X  <- 1 / lambda

# E[X^2] = 2/lambda^2, Var(X^2) = 20/lambda^4
mu_X2  <- 2 / lambda^2
sig_X2 <- sqrt(20) / lambda^2

# --- Generowanie danych ---
generate_means <- function(n) {
  samples <- matrix(rexp(N_sim * n, rate = lambda), nrow = N_sim, ncol = n)
  xbar    <- rowMeans(samples)
  x2bar   <- rowMeans(samples^2)
  z_xbar  <- sqrt(n) * (xbar  - mu_X)  / sig_X
  z_x2bar <- sqrt(n) * (x2bar - mu_X2) / sig_X2
  data.frame(n = n, xbar = xbar, x2bar = x2bar,
             z_xbar = z_xbar, z_x2bar = z_x2bar)
}

all_data <- do.call(rbind, lapply(n_sizes, generate_means))
all_data$n_label <- factor(paste0("n = ", all_data$n),
                            levels = paste0("n = ", n_sizes))

# --- Histogramy ---
plot_histogram <- function(data, zvar, title_stat) {
  ggplot(data, aes_string(x = zvar)) +
    geom_histogram(aes(y = after_stat(density)),
                   bins = 50, fill = "#4C9BE8", color = "white", alpha = 0.7) +
    stat_function(fun = dnorm, args = list(mean = 0, sd = 1),
                  color = "#E84C4C", linewidth = 1) +
    facet_wrap(~n_label, nrow = 1) +
    labs(title = paste("Histogram standaryzowanej", title_stat),
         subtitle = "Czerwona krzywa = N(0,1)",
         x = "Z", y = "Gestosc") +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))
}

p_hist_xbar  <- plot_histogram(all_data, "z_xbar",  "X_n")
p_hist_x2bar <- plot_histogram(all_data, "z_x2bar", "X^2_n")

# --- KDE ---
plot_kde <- function(data, zvar, title_stat) {
  ggplot(data, aes_string(x = zvar, color = "n_label", fill = "n_label")) +
    geom_density(alpha = 0.15, linewidth = 0.9) +
    stat_function(fun = dnorm, args = list(mean = 0, sd = 1),
                  color = "black", linewidth = 1.2, linetype = "dashed") +
    scale_color_brewer(palette = "RdYlBu", name = "n") +
    scale_fill_brewer(palette  = "RdYlBu", name = "n") +
    coord_cartesian(xlim = c(-4, 4)) +
    labs(title = paste("Estymatory jadrowe (KDE) dla", title_stat),
         subtitle = "Czarna przerywana = N(0,1)",
         x = "Z", y = "Gestosc") +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "bottom")
}

p_kde_xbar  <- plot_kde(all_data, "z_xbar",  "X_n")
p_kde_x2bar <- plot_kde(all_data, "z_x2bar", "X^2_n")

# --- ECDF ---
plot_ecdf <- function(data, zvar, title_stat) {
  ggplot(data, aes_string(x = zvar, color = "n_label")) +
    stat_ecdf(linewidth = 0.8, alpha = 0.9) +
    stat_function(fun = pnorm, args = list(mean = 0, sd = 1),
                  color = "black", linewidth = 1.2, linetype = "dashed") +
    scale_color_brewer(palette = "RdYlBu", name = "n") +
    coord_cartesian(xlim = c(-4, 4)) +
    labs(title = paste("Dystrybuanta empiryczna (ECDF) dla", title_stat),
         subtitle = "Czarna przerywana = N(0,1)",
         x = "Z", y = "F(z)") +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "bottom")
}

p_ecdf_xbar  <- plot_ecdf(all_data, "z_xbar",  "X_n")
p_ecdf_x2bar <- plot_ecdf(all_data, "z_x2bar", "X^2_n")

# --- QQ ---
plot_qq <- function(data, zvar, title_stat) {
  ggplot(data, aes_string(sample = zvar)) +
    stat_qq(aes(color = n_label), alpha = 0.3, size = 0.4) +
    stat_qq_line(color = "#E84C4C", linewidth = 0.9) +
    scale_color_brewer(palette = "RdYlBu", name = "n") +
    facet_wrap(~n_label, nrow = 1) +
    labs(title = paste("Wykres QQ dla", title_stat),
         subtitle = "Czerwona prosta = N(0,1)",
         x = "Kwantyle teoretyczne N(0,1)", y = "Kwantyle empiryczne") +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "none")
}

p_qq_xbar  <- plot_qq(all_data, "z_xbar",  "X_n")
p_qq_x2bar <- plot_qq(all_data, "z_x2bar", "X^2_n")

# --- Zapis ---
png("ctg_histogramy.png", width = 1400, height = 700, res = 120)
grid.arrange(p_hist_xbar, p_hist_x2bar, nrow = 2)
dev.off()

png("ctg_kde.png", width = 900, height = 700, res = 120)
grid.arrange(p_kde_xbar, p_kde_x2bar, nrow = 2)
dev.off()

png("ctg_ecdf.png", width = 900, height = 700, res = 120)
grid.arrange(p_ecdf_xbar, p_ecdf_x2bar, nrow = 2)
dev.off()

png("ctg_qqplot.png", width = 1400, height = 700, res = 120)
grid.arrange(p_qq_xbar, p_qq_x2bar, nrow = 2)
dev.off()

cat("Zapisano: ctg_histogramy.png, ctg_kde.png, ctg_ecdf.png, ctg_qqplot.png\n")
