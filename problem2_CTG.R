# =============================================================================
# PROBLEM 2: Ilustracja Centralnego Twierdzenia Granicznego (CTG)
# =============================================================================
# Rozkład wyjściowy: Wykładniczy exp(lambda=1)
#   - Wartość oczekiwana: mu = 1
#   - Wariancja: sigma^2 = 1
#   - Odchylenie standardowe: sigma = 1
#
# CTG mówi, że dla próby rozmiaru n, standaryzowana średnia:
#   Z_n = sqrt(n) * (X_bar_n - mu) / sigma
# zbiega do N(0,1) gdy n -> infinity
#
# Badamy dwie statystyki:
#   1. X_bar_n     = średnia z próby
#   2. X2_bar_n    = średnia z X^2 (funkcja próby)
#      Dla exp(1): E[X^2] = Var(X) + E[X]^2 = 1 + 1 = 2
#                 Var(X^2) = E[X^4] - (E[X^2])^2 = 9 - 4 = 5
# =============================================================================

library(ggplot2)
library(gridExtra)  # do łączenia wykresów obok siebie

set.seed(42)  # dla powtarzalności wyników

# =============================================================================
# PARAMETRY SYMULACJI
# =============================================================================

N_sim   <- 10000   # liczba powtórzeń symulacji (im więcej, tym gładsze wykresy)
n_sizes <- c(1, 5, 30, 100)  # rozmiary próby do zbadania
lambda  <- 1       # parametr rozkładu wykładniczego

# Parametry rozkładu exp(1)
mu_X    <- 1 / lambda          # E[X] = 1
sig_X   <- 1 / lambda          # sd(X) = 1

# Parametry dla X^2 gdy X ~ exp(1)
# E[X^2] = 2/lambda^2 = 2
# Var(X^2) = E[X^4] - (E[X^2])^2 = 24/lambda^4 - 4/lambda^4 = 20... 
# Dla lambda=1: E[X^2]=2, Var(X^2)=E[X^4]-(E[X^2])^2 = 24-4 = 20, sd=sqrt(20)
mu_X2   <- 2 / lambda^2
sig_X2  <- sqrt(20) / lambda^2  # odchylenie standardowe X^2

# =============================================================================
# GENEROWANIE DANYCH
# =============================================================================
# Dla każdego rozmiaru próby n:
#   - losujemy N_sim prób, każda rozmiaru n
#   - liczymy średnią X_bar i X2_bar dla każdej próby
#   - standaryzujemy: Z = sqrt(n) * (srednia - wartosc_oczekiwana) / odchylenie

generate_means <- function(n) {
  # Generujemy macierz N_sim x n z rozkładu exp(1)
  samples <- matrix(rexp(N_sim * n, rate = lambda), nrow = N_sim, ncol = n)
  
  # Statystyka 1: średnia z próby
  xbar   <- rowMeans(samples)
  
  # Statystyka 2: średnia z X^2
  x2bar  <- rowMeans(samples^2)
  
  # Standaryzacja (żeby porównać z N(0,1))
  z_xbar  <- sqrt(n) * (xbar  - mu_X)  / sig_X
  z_x2bar <- sqrt(n) * (x2bar - mu_X2) / sig_X2
  
  data.frame(
    n       = n,
    xbar    = xbar,
    x2bar   = x2bar,
    z_xbar  = z_xbar,
    z_x2bar = z_x2bar
  )
}

# Generujemy dane dla wszystkich rozmiarów próby
all_data <- do.call(rbind, lapply(n_sizes, generate_means))
all_data$n_label <- factor(paste0("n = ", all_data$n), 
                            levels = paste0("n = ", n_sizes))

# =============================================================================
# WIZUALIZACJA 1: HISTOGRAMY ze krzywą normalną
# =============================================================================
# Pokazuje: jak kształt rozkładu standaryzowanej średniej
#           przybliża się do dzwonu N(0,1) wraz ze wzrostem n

plot_histogram <- function(data, zvar, title_stat) {
  ggplot(data, aes_string(x = zvar)) +
    # Histogram (freq=FALSE -> oś Y to gęstość, nie liczebność)
    geom_histogram(aes(y = after_stat(density)),
                   bins = 50,
                   fill = "#4C9BE8", color = "white", alpha = 0.7) +
    # Krzywa N(0,1) do porównania
    stat_function(fun = dnorm, args = list(mean = 0, sd = 1),
                  color = "#E84C4C", linewidth = 1) +
    # Osobny panel dla każdego n
    facet_wrap(~n_label, nrow = 1) +
    labs(
      title = paste("Histogram standaryzowanej", title_stat),
      subtitle = "Czerwona krzywa = N(0,1) | Niebieski histogram = symulacja",
      x = "Wartość standaryzowana Z",
      y = "Gęstość"
    ) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))
}

p_hist_xbar  <- plot_histogram(all_data, "z_xbar",  "X̄ₙ")
p_hist_x2bar <- plot_histogram(all_data, "z_x2bar", "X̄²ₙ")

# =============================================================================
# WIZUALIZACJA 2: ESTYMATORY JĄDROWE (KDE)
# =============================================================================
# Pokazuje: wygładzoną krzywą gęstości naszych danych vs. N(0,1)
# KDE = każdy punkt danych "rozmywa się" w małą krzywą, suma daje gładką linię

plot_kde <- function(data, zvar, title_stat) {
  ggplot(data, aes_string(x = zvar, color = "n_label", fill = "n_label")) +
    # KDE - wygładzona gęstość empiryczna
    geom_density(alpha = 0.15, linewidth = 0.9) +
    # Krzywa N(0,1) do porównania (czarna, gruba)
    stat_function(fun = dnorm, args = list(mean = 0, sd = 1),
                  color = "black", linewidth = 1.2, linetype = "dashed") +
    scale_color_brewer(palette = "RdYlBu", name = "Rozmiar próby") +
    scale_fill_brewer(palette  = "RdYlBu", name = "Rozmiar próby") +
    coord_cartesian(xlim = c(-4, 4)) +
    labs(
      title = paste("Estymatory jądrowe (KDE) dla", title_stat),
      subtitle = "Czarna przerywana = N(0,1) | Kolorowe krzywe = kolejne n",
      x = "Wartość standaryzowana Z",
      y = "Gęstość"
    ) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "bottom")
}

p_kde_xbar  <- plot_kde(all_data, "z_xbar",  "X̄ₙ")
p_kde_x2bar <- plot_kde(all_data, "z_x2bar", "X̄²ₙ")

# =============================================================================
# WIZUALIZACJA 3: DYSTRYBUANTA EMPIRYCZNA (ECDF)
# =============================================================================
# Pokazuje: jak empiryczna dystrybuanta (schodkowa krzywa) 
#           zbiega do dystrybuanty N(0,1) (gładka S-krzywa)
# ECDF(x) = "ile procent moich danych jest mniejsze niż x?"

plot_ecdf <- function(data, zvar, title_stat) {
  ggplot(data, aes_string(x = zvar, color = "n_label")) +
    # Empiryczna dystrybuanta
    stat_ecdf(linewidth = 0.8, alpha = 0.9) +
    # Teoretyczna dystrybuanta N(0,1)
    stat_function(fun = pnorm, args = list(mean = 0, sd = 1),
                  color = "black", linewidth = 1.2, linetype = "dashed") +
    scale_color_brewer(palette = "RdYlBu", name = "Rozmiar próby") +
    coord_cartesian(xlim = c(-4, 4)) +
    labs(
      title = paste("Dystrybuanta empiryczna (ECDF) dla", title_stat),
      subtitle = "Czarna przerywana = Φ(z) rozkładu N(0,1) | Kolorowe = ECDF z symulacji",
      x = "Wartość standaryzowana Z",
      y = "F(z)"
    ) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "bottom")
}

p_ecdf_xbar  <- plot_ecdf(all_data, "z_xbar",  "X̄ₙ")
p_ecdf_x2bar <- plot_ecdf(all_data, "z_x2bar", "X̄²ₙ")

# =============================================================================
# WIZUALIZACJA 4: WYKRES KWANTYL-KWANTYL (QQ-PLOT)
# =============================================================================
# Pokazuje: czy kwantyle naszych danych odpowiadają kwantylom N(0,1)
# Jeśli punkty leżą na prostej y=x -> rozkład jest normalny
# Im bardziej odchylają się od prostej -> tym dalej od normalności

plot_qq <- function(data, zvar, title_stat) {
  ggplot(data, aes_string(sample = zvar)) +
    # Punkty QQ
    stat_qq(aes(color = n_label), alpha = 0.3, size = 0.4) +
    # Prosta referencyjna y=x (gdyby rozkład był idealnie normalny)
    stat_qq_line(color = "#E84C4C", linewidth = 0.9) +
    scale_color_brewer(palette = "RdYlBu", name = "Rozmiar próby") +
    facet_wrap(~n_label, nrow = 1) +
    labs(
      title = paste("Wykres QQ dla", title_stat),
      subtitle = "Czerwona prosta = linia N(0,1) | Punkty = dane z symulacji",
      x = "Kwantyle teoretyczne N(0,1)",
      y = "Kwantyle empiryczne"
    ) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "none")
}

p_qq_xbar  <- plot_qq(all_data, "z_xbar",  "X̄ₙ")
p_qq_x2bar <- plot_qq(all_data, "z_x2bar", "X̄²ₙ")

# =============================================================================
# ZAPIS DO PLIKÓW PNG
# =============================================================================
# Każda figura jako osobny plik (łatwiej wkleić do raportu)

# Figura 1: Histogramy
png("ctg_histogramy.png", width = 1400, height = 700, res = 120)
grid.arrange(p_hist_xbar, p_hist_x2bar, nrow = 2)
dev.off()

# Figura 2: KDE
png("ctg_kde.png", width = 900, height = 700, res = 120)
grid.arrange(p_kde_xbar, p_kde_x2bar, nrow = 2)
dev.off()

# Figura 3: ECDF
png("ctg_ecdf.png", width = 900, height = 700, res = 120)
grid.arrange(p_ecdf_xbar, p_ecdf_x2bar, nrow = 2)
dev.off()

# Figura 4: QQ-plot
png("ctg_qqplot.png", width = 1400, height = 700, res = 120)
grid.arrange(p_qq_xbar, p_qq_x2bar, nrow = 2)
dev.off()

cat("Zapisano 4 pliki PNG w katalogu roboczym:", getwd(), "\n")
cat("  - ctg_histogramy.png\n")
cat("  - ctg_kde.png\n")
cat("  - ctg_ecdf.png\n")
cat("  - ctg_qqplot.png\n")
