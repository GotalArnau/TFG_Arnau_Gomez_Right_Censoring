################################################################################
#-------------------------------------------------------------------------------
#In this script, the supporting images for our explanations throughout the study 
#are generated.
#-------------------------------------------------------------------------------
################################################################################


################################################################################
#-------------------------------------------------------------------------------
# Kaplan-Meier and data obtenied visualization
#-------------------------------------------------------------------------------
################################################################################

library(survival)
library(ggplot2)
library(dplyr)

set.seed(188)

n <- 500
dat_ex_asim <- data.frame(times = rweibull(n, 1, 20))
censt <- rexp(n, 1/13)
tobs_ex_asim <- pmin(dat_ex_asim$times, censt)
delta <- as.numeric(dat_ex_asim$times <= censt)
data <- data.frame(times = tobs_ex_asim, cens = delta)

times <- data$times
cens <- data$cens

km_fit <- survfit(Surv(times, cens) ~ 1)
t_max <- max(times[cens == 1])

S_km_tmax <- summary(km_fit, times = t_max)$surv


ln_fit <- survreg(Surv(times, cens) ~ 1, dist = "lognormal")

mu_hat    <- ln_fit$coefficients
sigma_hat <- ln_fit$scale


S_lognorm <- function(t, mu, sigma) {
  1 - plnorm(t, meanlog = mu, sdlog = sigma)
}

S0_tmax <- S_lognorm(t_max, mu_hat, sigma_hat)


Diff_tmax <- S_km_tmax - S0_tmax

km_df <- data.frame(
  time = km_fit$time,
  surv = km_fit$surv
)


t_grid <- seq(0, max(times), length.out = 500)

param_df <- data.frame(
  time = t_grid,
  surv = S_lognorm(t_grid, mu_hat, sigma_hat)
)

points_df <- data.frame(
  time = t_max,
  surv = c(S_km_tmax, S0_tmax),
  type = c("KM", "Log-normal")
)


(p1 <- ggplot() +  geom_step(data = km_df, aes(x = time, y = surv, color = "Kaplan-Meier"), linewidth = 1.1) +
   geom_line(data = param_df, aes(x = time, y = surv, color = "Lognormal"),  linewidth = 1.1, lty = 4) +
   geom_vline(xintercept = t_max, linetype = "dashed") +
   scale_color_manual(values = c("Kaplan-Meier" = "dodgerblue4", "Lognormal" = "firebrick3")) +
   geom_hline(yintercept = c(S_km_tmax, S0_tmax), linetype = "dotted", color = "darkgreen") +
   geom_point(data = points_df, aes(x = time, y = surv), color = "darkgreen", size = 2) +
   labs(x = "Time", y = "Survival Probability", color = "Color") +
   geom_segment(aes(x = 10, xend = 10, y = S_km_tmax, yend = S0_tmax),arrow = arrow(length = unit(0.5, "cm")),
                color = "orange", linewidth = 1.3) +
   geom_segment(aes(x = 10, xend = 10, y = S0_tmax, yend = S_km_tmax),arrow = arrow(length = unit(0.5, "cm")),
                color = "orange", linewidth = 1.3) +
   theme_bw() + 
   annotate("text", x = t_max, y = S_km_tmax + 0.07, label = expression(hat(S)(t[m])), hjust = -0.2, color = "darkgreen", 
            size = 5) +
   annotate("text", x = t_max, y = S0_tmax + 0.07, label = expression((S)[0](t[m], hat(theta))), hjust = -0.1, 
            color = "darkgreen", size = 5) +
   annotate("text", x = 10, y = (S_km_tmax + S0_tmax) / 2, label = expression(Diff[t[max]]), hjust = -0.2,
            color = "orange", size = 5) +
   ylim(0,1) + theme(legend.position = "bottom", 
                     legend.title = element_text(size = 18, face = "bold"),
                     legend.text = element_text(size = 18),
                     plot.subtitle = element_text(size = 18, face = "italic", hjust = 0.5),
                     axis.title = element_text(size = 18, face = "bold"),
                     axis.text = element_text(size = 15),
                     plot.title = element_text(size = 25, face = "bold")))

ggsave(filename = "04_Output/Figures/Figure01_KM_visualization.png", plot = p1, width = 2130, height = 1600, units = "px")



################################################################################
#-------------------------------------------------------------------------------
# Random vs Administrative Censoring Comparison
#-------------------------------------------------------------------------------
################################################################################

library(GofCens)
library(grid)
library(gridExtra)
library(survival)
library(ggplot2)
library(dplyr)
library(patchwork)

set.seed(28657)

dat_asim_0.7 <- data.frame(times = rweibull(100, 2, 20))
quant <- quantile(dat_asim_0.7$times,0.7)
cens <- c()
cens <- ifelse(dat_asim_0.7$times < quant, 1, 0)
dat_asim_0.7$times[dat_asim_0.7$times > quant] <- quant
dat_asim_0.7$cens <- cens

sum(dat_asim_0.7$cens==0)/length(dat_asim_0.7$cens)


km_fit <- survfit(Surv(dat_asim_0.7$times, dat_asim_0.7$cens) ~ 1)
km_df <- data.frame(time = km_fit$time, surv = km_fit$surv)
km_surv_max <- km_df$surv[km_df$time == max(dat_asim_0.7$times[dat_asim_0.7$cens == 1])]

ln_fit <- survreg(Surv(dat_asim_0.7$times, dat_asim_0.7$cens) ~ 1, dist = "lognormal")
mu_hat <- ln_fit$coefficients
sigma_hat <- ln_fit$scale

S_lognorm <- function(t, mu, sigma) 1 - plnorm(t, mu, sigma)
t_grid <- seq(0, max(dat_asim_0.7$times), length.out = 500)
param_df <- data.frame(time = t_grid, surv = S_lognorm(t_grid, mu_hat, sigma_hat))

p2 <- ggplot() +
    geom_step(data = km_df, aes(x = time, y = surv, color = "Kaplan-Meier"), linewidth = 1.1) +
    labs(
      subtitle = "Survival function with 30% Administrative censoring",
      x = "Time",
      y = "Survival Probability"
    ) +
    scale_color_manual(values = c("Kaplan-Meier" = "blue", "Log-normal" = "red")) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, face = "italic", hjust = 0.5),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      legend.position = c(0.835, 0.82),
      legend.background = element_rect(fill = "white", color = "black"),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10)
    ) + scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.15)) + 
    geom_hline(yintercept = 0.3, color = "darkolivegreen", linetype = 2, linewidth = 1) + 
    annotate("text", x = 3, y = 0.37, label = round(km_surv_max,2), size = 5, color = "darkolivegreen")

(p2 <- p2 + theme(legend.position = "none"))



set.seed(427)

dat_asim_0.7 <- data.frame(times = rweibull(100, 2, 20))
censt <- rexp(100, 1/47.5)
quant <- rep(quantile(dat_asim_0.7$times, probs = 0.9), 100)
tobs_asim_0.7 <- pmin(dat_asim_0.7$times, censt, quant)
delta <- as.numeric(dat_asim_0.7$times <= censt & dat_asim_0.7$times <= quant)
dat_asim_0.7 <- data.frame(times = tobs_asim_0.7, cens = delta)

sum(dat_asim_0.7$cens==0)/length(dat_asim_0.7$cens)


km_fit <- survfit(Surv(dat_asim_0.7$times, dat_asim_0.7$cens) ~ 1)
km_df <- data.frame(time = km_fit$time, surv = km_fit$surv)
km_surv_max <- km_df$surv[km_df$time == max(dat_asim_0.7$times[dat_asim_0.7$cens == 1])]

ln_fit <- survreg(Surv(dat_asim_0.7$times, dat_asim_0.7$cens) ~ 1, dist = "lognormal")
mu_hat <- ln_fit$coefficients
sigma_hat <- ln_fit$scale

S_lognorm <- function(t, mu, sigma) 1 - plnorm(t, mu, sigma)
t_grid <- seq(0, max(dat_asim_0.7$times), length.out = 500)
param_df <- data.frame(time = t_grid, surv = S_lognorm(t_grid, mu_hat, sigma_hat))

p3 <- ggplot() +
    geom_step(data = km_df, aes(x = time, y = surv, color = "Kaplan-Meier"), linewidth = 1.1) +
    labs(
      subtitle = "Survival function with 30% Random censoring",
      x = "Time",
      y = "Survival Probability",
      color = "Curves"
    ) +
    scale_color_manual(values = c("Kaplan-Meier" = "blue", "Log-normal" = "red")) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, face = "italic", hjust = 0.5),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      legend.position = c(0.835, 0.82),
      legend.background = element_rect(fill = "white", color = "black"),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10)
    ) + scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.15)) + 
    geom_hline(yintercept = km_surv_max, color = "darkolivegreen", linetype = 2, linewidth = 1) + 
    annotate("text", x = 3, y = km_surv_max+0.07, label = round(km_surv_max,2),
             size = 5, color = "darkolivegreen")

(p3 <- p3 + theme(legend.position = "none"))


(combined_plots <- (p2 + p3 + plot_layout(guides = "collect") & theme(plot.subtitle = element_text(size = 18, face = "italic", hjust = 0.5),
                                                                      axis.title = element_text(size = 20, face = "bold"),
                                                                      axis.text = element_text(size = 20))) +
    plot_annotation(theme =  theme(plot.title = element_text(size = 30, face = "bold", hjust = 0.5))))

ggsave(filename = "04_Output/Figures/Figure04_Censoring_Visualization_Comparison.png", plot = combined_plots, width = 2130*2, height = 2000, units = "px")



################################################################################
#-------------------------------------------------------------------------------
# Density visualization by Skewness
#-------------------------------------------------------------------------------
################################################################################

library(ggplot2)
library(patchwork)


###############################################################################
# Weibull
################################################################################

x <- seq(0, 100, length.out = 1000)

scale_param <- 20

shape_params <- c(3.5, 1.7, 1)

df <- data.frame(
  x = rep(x, times = length(shape_params)),
  shape = factor(rep(shape_params, each = length(x)))
)

df$y <- with(df, dweibull(x, shape = as.numeric(as.character(shape)), scale = scale_param))


(p4 <- ggplot(df, aes(x = x, y = y, color = shape)) +
    geom_line(linewidth = 1) +
    labs(
      x = "Time",
      y = "Density",
    color = NULL
    ) +
    scale_color_discrete(labels = c(
      "3.5" = expression("Symmetric " * alpha == 3.5),
      "1.7" = expression("Moderately skewed " * alpha == 1.7 * " "),
      "1"   = expression("Highly skewed " * alpha == 1)
    ))  +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = c(0.9, 0.85),
      legend.justification = c(1, 1), 
      legend.background = element_rect(fill = "white", color = "black"),
      legend.key = element_rect(fill = "white", color = NA),
      legend.text = element_text(size = 20)
    ))

ggsave(filename = "04_Output/Figures/Weibull_Density_Fucntions.png", plot = p4, width = 2130*2, height = 2000, units = "px")

################################################################################
# Lognormal
################################################################################

x <- seq(0, 100, length.out = 1000)

mu_param <- 2

sdlog_params <- c(0.01, 0.4, 1.5)

df <- data.frame(
  x = rep(x, times = length(sdlog_params)),
  sdlog = factor(rep(sdlog_params, each = length(x)))
)

df$y <- with(df, dlnorm(x, sdlog = as.numeric(as.character(sdlog)), meanlog = mu_param))


(p4 <- ggplot(df, aes(x = x, y = y, color = sdlog)) +
    geom_line(linewidth = 1) +
    labs(
      x = "Time",
      y = "Density",
      color = NULL
    ) +
    scale_color_discrete(labels = c(
      "0.01" = expression("Symmetric " * sdlog == 0.01),
      "0.4" = expression("Moderately skewed " * sdlog == 0.4 * " "),
      "1.5"   = expression("Highly skewed " * sdlog == 1.5)
    ))  +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = c(0.9, 0.85),
      legend.justification = c(1, 1), 
      legend.background = element_rect(fill = "white", color = "black"),
      legend.key = element_rect(fill = "white", color = NA)
    ))


(p4_zoom <- ggplot(df, aes(x = x, y = y, color = sdlog)) +
    geom_line(linewidth = 1) +
    labs(
      x = "Time",
      y = "Density",
      color = NULL
    ) +
    scale_color_discrete(labels = c(
      "0.01" = expression("Symmetric " * sdlog == 0.01),
      "0.4" = expression("Moderately skewed " * sdlog == 0.4 * " "),
      "1.5"   = expression("Highly skewed " * sdlog == 1.5)
    ))  +
    coord_cartesian(xlim = c(0, 30), ylim = c(0, 0.2)) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = c(0.9, 0.85),
      legend.justification = c(1, 1), 
      legend.background = element_rect(fill = "white", color = "black"),
      legend.key = element_rect(fill = "white", color = NA)
    ))


(combine_plots2 <- (p4 + p4_zoom + plot_layout(guides = "collect") &
                     theme(
                       legend.position = "bottom",
                       legend.direction = "horizontal",
                       legend.text = element_text(size = 20)
                       
                     ) ))


ggsave(filename = "04_Output/Figures/Lognormal_Density_Fucntions.png", plot = combine_plots2, width = 2130*2, height = 2000, units = "px")


################################################################################
#-------------------------------------------------------------------------------
# Power comparison for complete data by Null Hypothesis
#-------------------------------------------------------------------------------
################################################################################

library(ggplot2)
library(dplyr)
library(patchwork)

#----------------------------------------------------------------------------------------------------------------------------
# X ~ Weibull
#----------------------------------------------------------------------------------------------------------------------------

df_W_logistic <- get(load("02_Data/Weibull_RandomCensoring/df_power_all4_logistic/Power/_Distribution_weibull_resultados_h0_logistic_all_.RData"))
df_W_log <- get(load("02_Data/Weibull_RandomCensoring/df_power_all4_lognormal/Power/_Distribution_weibull_resultados_h0_lognormal_all_.RData"))
df_W_wei <- get(load("02_Data/Weibull_RandomCensoring/df_power_all4_weibull/Power/_Distribution_weibull_resultados_h0_weibull_all_.RData"))


df_W_logistic2 <- df_W_logistic %>%
  select(power, cens, test) %>%
  filter(cens == "0% Censura") %>%
  mutate("Hyp_nul" = rep("Logistic"), cens = NULL)

df_W_log2 <- df_W_log %>%
  select(power, cens, test) %>%
  filter(cens == "0% Censura") %>%
  mutate("Hyp_nul" = rep("Lognormal"), cens = NULL)

df_W_wei2 <- df_W_wei %>%
  select(power, cens, test) %>%
  filter(cens == "0% Censura") %>%
  mutate("Hyp_nul" = rep("Weibull"), cens = NULL)

df_complet <- rbind(df_W_logistic2, df_W_log2, df_W_wei2 )%>%
  filter(test != "KS No Bootstrap")


(p2 <- ggplot(df_complet, aes(x = test, y = power, shape = Hyp_nul, fill = Hyp_nul)) + geom_jitter(size = 3, stroke = 0.5, height = 0, width = 0.1) +
    scale_shape_manual(values = c(22,24,23)) +
    scale_fill_manual(values = c( "Logistic" = "#CAFF70", "Lognormal" = "#79CDCD","Weibull" = "#FF6A6A")) +
    ggtitle(paste0("Power of each tests for complete data by Null Hypothesis"), subtitle = paste0("X ~ Weibull")) +
    labs(x = "Test", y = "Power", fill = "Null Hypothesis:") +
    theme_bw() +
    theme(legend.position = "right", ) +
    lims(y = c(-0.01,1.01)) +
    guides(fill = guide_legend(override.aes = list(shape = c(22,24,23), color = "black")),
           shape = "none") +
    geom_hline(yintercept = 0.1, colour = "darkblue", linetype = 2, linewidth = 0.9) +
    annotate("text", x = 0.5, y = 0.07, label = "alpha = 0.1", hjust = -0.02, color = "darkblue", 
             size = 3.5, angle = 330))

ggsave(filename = "C:/Users/arnau.gomez/Desktop/GofCens_Paper_Simulations/04_Output/02_Plots/New_plots/Figure3.1_Weibull_Power_complete_data.png", 
       plot = p2, width = 2130, height = 1600, units = "px")

p2 <- p2 + labs(title = NULL)

#----------------------------------------------------------------------------------------------------------------------------
# X ~ Logistic
#----------------------------------------------------------------------------------------------------------------------------

df_L_logistic <- get(load("02_Data/Lognormal_RandomCensoring/df_power_all4_logistic/Power/_Distribution_lognormal_resultados_h0_logistic_all_.RData"))
df_L_log <- get(load("02_Data/Lognormal_RandomCensoring/df_power_all4_lognormal/Power/_Distribution_lognormal_resultados_h0_lognormal_all_.RData"))
df_L_wei <- get(load("02_Data/Lognormal_RandomCensoring/df_power_all4_weibull/Power/_Distribution_lognormal_resultados_h0_weibull_all_.RData"))


df_L_logistic2 <- df_L_logistic %>%
  select(power, cens, test) %>%
  filter(cens == "0% Censura") %>%
  mutate("Hyp_nul" = rep("Logistic"), cens = NULL)

df_L_log2 <- df_L_log %>%
  select(power, cens, test) %>%
  filter(cens == "0% Censura") %>%
  mutate("Hyp_nul" = rep("Lognormal"), cens = NULL)

df_L_wei2 <- df_L_wei %>%
  select(power, cens, test) %>%
  filter(cens == "0% Censura") %>%
  mutate("Hyp_nul" = rep("Weibull"), cens = NULL)

df_complet <- rbind(df_L_logistic2, df_L_log2, df_L_wei2) %>%
  filter(test != "KS No Bootstrap")


(p3 <- ggplot(df_complet, aes(x = test, y = power, shape = Hyp_nul, fill = Hyp_nul)) + geom_jitter(size = 3, stroke = 0.5, height = 0, width = 0.1) +
    scale_shape_manual(values = c(22,24,23)) +
    scale_fill_manual(values = c( "Logistic" = "#CAFF70", "Lognormal" = "#79CDCD","Weibull" = "#FF6A6A")) +
    ggtitle(paste0("Power of each tests for complete data by Null Hypothesis"), subtitle = paste0("X ~ Lognormal")) +
    labs(x = "Test", y = "Power", fill = "Null Hypothesis:") +
    theme_bw() +
    theme(legend.position = "right", ) +
    lims(y = c(-0.01,1.01)) +
    guides(fill = guide_legend(override.aes = list(shape = c(22,24,23), color = "black")),
           shape = "none") +
    geom_hline(yintercept = 0.1, colour = "darkblue", linetype = 2, linewidth = 0.9) +
    annotate("text", x = 0.5, y = 0.07, label = "alpha = 0.1", hjust = -0.02, color = "darkblue", 
             size = 3.5, angle = 330))

ggsave(filename = "C:/Users/arnau.gomez/Desktop/GofCens_Paper_Simulations/04_Output/02_Plots/New_plots/Figure3.1_Lognoral_Power_complete_data.png", 
       plot = p3, width = 2130, height = 1600, units = "px")

p3 <- p3 + labs(title=NULL)

(combined_complete_data <- (p2 + p3 + plot_layout(guides = "collect") & theme(legend.position = "bottom", 
                                                                              legend.title = element_text(size = 20, face = "bold"),
                                                                              legend.text = element_text(size = 20),
                                                                              plot.subtitle = element_text(size = 20, face = "italic", hjust = 0.5),
                                                                              axis.title = element_text(size = 15, face = "bold"),
                                                                              axis.text = element_text(size = 15))) +
    plot_annotation(title = NULL, 
                    theme =  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))))

ggsave(filename = "04_Output/Figures/Power_Complete_Data_Null.png", plot = combined_complete_data, width = 2130*2, height = 1800, units = "px")



################################################################################
#-------------------------------------------------------------------------------
# Power comparison for complete data by Sample Size
#-------------------------------------------------------------------------------
################################################################################

library(ggplot2)
library(dplyr)
library(patchwork)

#----------------------------------------------------------------------------------------------------------------------------
# X ~ Weibull
#----------------------------------------------------------------------------------------------------------------------------

df_W_logistic <- get(load("02_Data/Weibull_RandomCensoring/df_power_all4_logistic/Power/_Distribution_weibull_resultados_h0_logistic_all_.RData"))
df_W_log <- get(load("02_Data/Weibull_RandomCensoring/df_power_all4_lognormal/Power/_Distribution_weibull_resultados_h0_lognormal_all_.RData"))
df_W_wei <- get(load("02_Data/Weibull_RandomCensoring/df_power_all4_weibull/Power/_Distribution_weibull_resultados_h0_weibull_all_.RData"))


df_W_logistic2 <- df_W_logistic %>%
  select(power, cens, test, n) %>%
  filter(cens == "0% Censura") %>%
  mutate("Hyp_nul" = rep("Logistic"), cens = NULL)%>%
  filter(test != "KS No Bootstrap")

df_W_log2 <- df_W_log %>%
  select(power, cens, test, n) %>%
  filter(cens == "0% Censura") %>%
  mutate("Hyp_nul" = rep("Lognormal"), cens = NULL)%>%
  filter(test != "KS No Bootstrap")

df_W_wei2 <- df_W_wei %>%
  select(power, cens, test, n) %>%
  filter(cens == "0% Censura") %>%
  mutate("Hyp_nul" = rep("Weibull"), cens = NULL)%>%
  filter(test != "KS No Bootstrap")

df_complet <- rbind(df_W_logistic2, df_W_log2, df_W_wei2 )%>%
  filter(test != "KS No Bootstrap")


(p2 <- ggplot(df_W_logistic2, aes(x = test, y = power, shape = factor(n), fill = factor(n))) + geom_jitter(size = 3, stroke = 0.5, height = 0, width = 0.1) +
    scale_shape_manual(values = c(22,24,23)) +
    scale_fill_manual(values = c( "75" = "#CAFF70", "150" = "#79CDCD","300" = "#FF6A6A")) +
    ggtitle(paste0("Power of each tests for complete data by Sample Size"), subtitle = paste0("X ~ Weibull. H0: Logistic")) +
    labs(x = "Test", y = "Power", fill = "Sample Size:") +
    theme_bw() +
    theme(legend.position = "right", ) +
    lims(y = c(-0.01,1.01)) +
    guides(fill = guide_legend(override.aes = list(shape = c(22,24,23), color = "black")),
           shape = "none") +
    geom_hline(yintercept = 0.1, colour = "darkblue", linetype = 2, linewidth = 0.9) +
    annotate("text", x = 0.5, y = 0.07, label = "alpha = 0.1", hjust = -0.02, color = "darkblue", 
             size = 3.5, angle = 330))

(p3 <- ggplot(df_W_log2, aes(x = test, y = power, shape = factor(n), fill = factor(n))) + geom_jitter(size = 3, stroke = 0.5, height = 0, width = 0.1) +
    scale_shape_manual(values = c(22,24,23)) +
    scale_fill_manual(values = c( "75" = "#CAFF70", "150" = "#79CDCD","300" = "#FF6A6A")) +
    ggtitle(paste0("Power of each tests for complete data by Sample Size"), subtitle = paste0("X ~ Weibull. H0: Lognormal")) +
    labs(x = "Test", y = "Power", fill = "Sample Size:") +
    theme_bw() +
    theme(legend.position = "right", ) +
    lims(y = c(-0.01,1.01)) +
    guides(fill = guide_legend(override.aes = list(shape = c(22,24,23), color = "black")),
           shape = "none") +
    geom_hline(yintercept = 0.1, colour = "darkblue", linetype = 2, linewidth = 0.9) +
    annotate("text", x = 0.5, y = 0.07, label = "alpha = 0.1", hjust = -0.02, color = "darkblue", 
             size = 3.5, angle = 330))

(p4 <- ggplot(df_W_wei2, aes(x = test, y = power, shape = factor(n), fill = factor(n))) + geom_jitter(size = 3, stroke = 0.5, height = 0, width = 0.1) +
    scale_shape_manual(values = c(22,24,23)) +
    scale_fill_manual(values = c( "75" = "#CAFF70", "150" = "#79CDCD","300" = "#FF6A6A")) +
    ggtitle(paste0("Power of each tests for complete data by Sample Size"), subtitle = paste0("X ~ Weibull. H0: Weibull")) +
    labs(x = "Test", y = "Power", fill = "Sample Size:") +
    theme_bw() +
    theme(legend.position = "right", ) +
    lims(y = c(-0.01,1.01)) +
    guides(fill = guide_legend(override.aes = list(shape = c(22,24,23), color = "black")),
           shape = "none") +
    geom_hline(yintercept = 0.1, colour = "darkblue", linetype = 2, linewidth = 0.9) +
    annotate("text", x = 0.5, y = 0.07, label = "alpha = 0.1", hjust = -0.02, color = "darkblue", 
             size = 3.5, angle = 330))

p2 <- p2 + labs(title = NULL)
p3 <- p3 + labs(title = NULL)
p4 <- p4 + labs(title = NULL)

(combined_complete_data_W <- (p2 + p3 + p4 + plot_layout(guides = "collect") & theme(legend.position = "bottom", 
                                                                              legend.title = element_text(size = 20, face = "bold"),
                                                                              legend.text = element_text(size = 20),
                                                                              plot.subtitle = element_text(size = 20, face = "italic", hjust = 0.5),
                                                                              axis.title = element_text(size = 15, face = "bold"),
                                                                              axis.text = element_text(size = 15))) +
    plot_annotation(title = NULL, 
                    theme =  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))))

ggsave(filename = "04_Output/Figures/Power_Complete_Data_SampleSize_Weibull.png", plot = combined_complete_data_W, width = 2130*3, height = 1800, units = "px")


#----------------------------------------------------------------------------------------------------------------------------
# X ~ Logistic
#----------------------------------------------------------------------------------------------------------------------------

df_L_logistic <- get(load("C:/Users/arnau.gomez/Desktop/GofCens_Paper_Simulations/02_Data/Lognormal_RandomCensoring/df_power_all4_logistic/SEC/_Distribution_lognormal_resultados_h0_logistic_all_.RData"))
df_L_log <- get(load("C:/Users/arnau.gomez/Desktop/GofCens_Paper_Simulations/02_Data/Lognormal_RandomCensoring/df_power_all4_lognormal/SEC/_Distribution_lognormal_resultados_h0_lognormal_all_.RData"))
df_L_wei <- get(load("C:/Users/arnau.gomez/Desktop/GofCens_Paper_Simulations/02_Data/Lognormal_RandomCensoring/df_power_all4_weibull/SEC/_Distribution_lognormal_resultados_h0_weibull_all_.RData"))


df_L_logistic2 <- df_L_logistic %>%
  select(power, cens, test, n) %>%
  filter(cens == "0% Censura") %>%
  mutate("Hyp_nul" = rep("Logistic"), cens = NULL)%>%
  filter(test != "KS No Bootstrap")

df_L_log2 <- df_L_log %>%
  select(power, cens, test, n) %>%
  filter(cens == "0% Censura") %>%
  mutate("Hyp_nul" = rep("Lognormal"), cens = NULL)%>%
  filter(test != "KS No Bootstrap")

df_L_wei2 <- df_L_wei %>%
  select(power, cens, test, n) %>%
  filter(cens == "0% Censura") %>%
  mutate("Hyp_nul" = rep("Weibull"), cens = NULL)%>%
  filter(test != "KS No Bootstrap")

df_complet <- rbind(df_L_logistic2, df_L_log2, df_L_wei2) %>%
  filter(test != "KS No Bootstrap")


(p2 <- ggplot(df_L_logistic2, aes(x = test, y = power, shape = factor(n), fill = factor(n))) + geom_jitter(size = 3, stroke = 0.5, height = 0, width = 0.1) +
    scale_shape_manual(values = c(22,24,23)) +
    scale_fill_manual(values = c( "75" = "#CAFF70", "150" = "#79CDCD","300" = "#FF6A6A")) +
    ggtitle(paste0("Power of each tests for complete data by Sample Size"), subtitle = paste0("X ~ Lognormal H0: Logistic")) +
    labs(x = "Test", y = "Power", fill = "Sample Size:") +
    theme_bw() +
    theme(legend.position = "right", ) +
    lims(y = c(-0.01,1.01)) +
    guides(fill = guide_legend(override.aes = list(shape = c(22,24,23), color = "black")),
           shape = "none") +
    geom_hline(yintercept = 0.1, colour = "darkblue", linetype = 2, linewidth = 0.9) +
    annotate("text", x = 0.5, y = 0.07, label = "alpha = 0.1", hjust = -0.02, color = "darkblue", 
             size = 3.5, angle = 330))

(p3 <- ggplot(df_L_log2, aes(x = test, y = power, shape = factor(n), fill = factor(n))) + geom_jitter(size = 3, stroke = 0.5, height = 0, width = 0.1) +
    scale_shape_manual(values = c(22,24,23)) +
    scale_fill_manual(values = c( "75" = "#CAFF70", "150" = "#79CDCD","300" = "#FF6A6A")) +
    ggtitle(paste0("Power of each tests for complete data by Sample Size"), subtitle = paste0("X ~ Lognormal H0: Lognormal")) +
    labs(x = "Test", y = "Power", fill = "Sample Size:") +
    theme_bw() +
    theme(legend.position = "right", ) +
    lims(y = c(-0.01,1.01)) +
    guides(fill = guide_legend(override.aes = list(shape = c(22,24,23), color = "black")),
           shape = "none") +
    geom_hline(yintercept = 0.1, colour = "darkblue", linetype = 2, linewidth = 0.9) +
    annotate("text", x = 0.5, y = 0.07, label = "alpha = 0.1", hjust = -0.02, color = "darkblue", 
             size = 3.5, angle = 330))

(p4 <- ggplot(df_L_wei2, aes(x = test, y = power, shape = factor(n), fill = factor(n))) + geom_jitter(size = 3, stroke = 0.5, height = 0, width = 0.1) +
    scale_shape_manual(values = c(22,24,23)) +
    scale_fill_manual(values = c( "75" = "#CAFF70", "150" = "#79CDCD","300" = "#FF6A6A")) +
    ggtitle(paste0("Power of each tests for complete data by Sample Size"), subtitle = paste0("X ~ Lognormal H0: Weibull")) +
    labs(x = "Test", y = "Power", fill = "Sample Size:") +
    theme_bw() +
    theme(legend.position = "right", ) +
    lims(y = c(-0.01,1.01)) +
    guides(fill = guide_legend(override.aes = list(shape = c(22,24,23), color = "black")),
           shape = "none") +
    geom_hline(yintercept = 0.1, colour = "darkblue", linetype = 2, linewidth = 0.9) +
    annotate("text", x = 0.5, y = 0.07, label = "alpha = 0.1", hjust = -0.02, color = "darkblue", 
             size = 3.5, angle = 330))

p2 <- p2 + labs(title = NULL)
p3 <- p3 + labs(title = NULL)
p4 <- p4 + labs(title = NULL)

(combined_complete_data_L <- (p2 + p3 + p4 + plot_layout(guides = "collect") & theme(legend.position = "bottom", 
                                                                                   legend.title = element_text(size = 20, face = "bold"),
                                                                                   legend.text = element_text(size = 20),
                                                                                   plot.subtitle = element_text(size = 20, face = "italic", hjust = 0.5),
                                                                                   axis.title = element_text(size = 15, face = "bold"),
                                                                                   axis.text = element_text(size = 15))) +
    plot_annotation(title = NULL, 
                    theme =  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))))

ggsave(filename = "04_Output/Figures/Power_Complete_Data_SampleSize_Lognormal.png", plot = combined_complete_data_L, width = 2130*3, height = 1800, units = "px")


(p_final <- (combined_complete_data_W / combined_complete_data_L + plot_layout(guides = "collect") & theme(legend.position = "bottom", 
                                                                                                         legend.title = element_text(size = 30, face = "bold"),
                                                                                                         legend.text = element_text(size = 30),
                                                                                                         plot.subtitle = element_text(size = 30, face = "italic", hjust = 0.5),
                                                                                                         axis.title = element_text(size = 15, face = "bold"),
                                                                                                         axis.text = element_text(size = 15))) +
    plot_annotation(title = NULL, 
                    theme =  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))))

ggsave(filename = "04_Output/Figures/Power_Complete_Data_SampleSize.png", plot = p_final, width = 2130*3, height = 2*1800, units = "px")

################################################################################
#-------------------------------------------------------------------------------
# Power comparison for administrative censored data by Null Hypothesis
#-------------------------------------------------------------------------------
################################################################################

library(ggplot2)
library(dplyr)
library(patchwork)

#----------------------------------------------------------------------------------------------------------------------------
# X ~ Weibull
#----------------------------------------------------------------------------------------------------------------------------

df_W_logistic <- get(load("02_Data/Weibull_AdminCensoring/df_power_all4_logistic/Power/_Distribution_weibull_resultados_h0_logistic_all_.RData"))
df_W_log <- get(load("02_Data/Weibull_AdminCensoring/df_power_all4_lognormal/Power/_Distribution_weibull_resultados_h0_lognormal_all_.RData"))
df_W_wei <- get(load("02_Data/Weibull_AdminCensoring/df_power_all4_weibull/Power/_Distribution_weibull_resultados_h0_weibull_all_.RData"))


df_W_logistic2 <- df_W_logistic %>%
  select(power, cens, test) %>%
  filter(cens == "30% Censura") %>%
  mutate("Hyp_nul" = rep("Logistic"), cens = NULL)

df_W_log2 <- df_W_log %>%
  select(power, cens, test) %>%
  filter(cens == "30% Censura") %>%
  mutate("Hyp_nul" = rep("Lognormal"), cens = NULL)

df_W_wei2 <- df_W_wei %>%
  select(power, cens, test) %>%
  filter(cens == "30% Censura") %>%
  mutate("Hyp_nul" = rep("Weibull"), cens = NULL)

df_complet <- rbind(df_W_logistic2, df_W_log2, df_W_wei2 )%>%
  filter(test != "KS No Bootstrap")


(p2 <- ggplot(df_complet, aes(x = test, y = power, shape = Hyp_nul, fill = Hyp_nul)) + geom_jitter(size = 3, stroke = 0.5, height = 0, width = 0.1) +
    scale_shape_manual(values = c(22,24,23)) +
    scale_fill_manual(values = c( "Logistic" = "#CAFF70", "Lognormal" = "#79CDCD","Weibull" = "#FF6A6A")) +
    ggtitle(paste0("Power of each Tests for Administrative Censored Data by Null Hypothesis"), subtitle = paste0("X ~ Weibull. 30% Censoring")) +
    labs(x = "Test", y = "Power", fill = "Null Hypothesis:") +
    theme_bw() +
    theme(legend.position = "right", ) +
    lims(y = c(-0.01,1.01)) +
    guides(fill = guide_legend(override.aes = list(shape = c(22,24,23), color = "black")),
           shape = "none") +
    geom_hline(yintercept = 0.1, colour = "darkblue", linetype = 2, linewidth = 0.9) +
    annotate("text", x = 0.5, y = 0.07, label = "alpha = 0.1", hjust = -0.02, color = "darkblue", 
             size = 3.5, angle = 330))

ggsave(filename = "C:/Users/arnau.gomez/Desktop/GofCens_Paper_Simulations/04_Output/02_Plots/New_plots/Figure3.1_Weibull_Power_complete_data.png", 
       plot = p2, width = 2130, height = 1600, units = "px")

p2 <- p2 + labs(title = NULL)

#----------------------------------------------------------------------------------------------------------------------------
# X ~ Logistic
#----------------------------------------------------------------------------------------------------------------------------

df_L_logistic <- get(load("02_Data/Lognormal_AdminCensoring/df_power_all4_logistic/Power/_Distribution_lognormal_resultados_h0_logistic_all_.RData"))
df_L_log <- get(load("02_Data/Lognormal_AdminCensoring/df_power_all4_lognormal/Power/_Distribution_lognormal_resultados_h0_lognormal_all_.RData"))
df_L_wei <- get(load("02_Data/Lognormal_AdminCensoring/df_power_all4_weibull/Power/_Distribution_lognormal_resultados_h0_weibull_all_.RData"))


df_L_logistic2 <- df_L_logistic %>%
  select(power, cens, test) %>%
  filter(cens == "30% Censura") %>%
  mutate("Hyp_nul" = rep("Logistic"), cens = NULL)

df_L_log2 <- df_L_log %>%
  select(power, cens, test) %>%
  filter(cens == "30% Censura") %>%
  mutate("Hyp_nul" = rep("Lognormal"), cens = NULL)

df_L_wei2 <- df_L_wei %>%
  select(power, cens, test) %>%
  filter(cens == "30% Censura") %>%
  mutate("Hyp_nul" = rep("Weibull"), cens = NULL)

df_complet <- rbind(df_L_logistic2, df_L_log2, df_L_wei2) %>%
  filter(test != "KS No Bootstrap")


(p3 <- ggplot(df_complet, aes(x = test, y = power, shape = Hyp_nul, fill = Hyp_nul)) + geom_jitter(size = 3, stroke = 0.5, height = 0, width = 0.1) +
    scale_shape_manual(values = c(22,24,23)) +
    scale_fill_manual(values = c( "Logistic" = "#CAFF70", "Lognormal" = "#79CDCD","Weibull" = "#FF6A6A")) +
    ggtitle(paste0("Power of each Tests for Administrative Censored Data by Null Hypothesis"), subtitle = paste0("X ~ Lognormal. 30% Censoring")) +
    labs(x = "Test", y = "Power", fill = "Null Hypothesis:") +
    theme_bw() +
    theme(legend.position = "right", ) +
    lims(y = c(-0.01,1.01)) +
    guides(fill = guide_legend(override.aes = list(shape = c(22,24,23), color = "black")),
           shape = "none") +
    geom_hline(yintercept = 0.1, colour = "darkblue", linetype = 2, linewidth = 0.9) +
    annotate("text", x = 0.5, y = 0.07, label = "alpha = 0.1", hjust = -0.02, color = "darkblue", 
             size = 3.5, angle = 330))

ggsave(filename = "C:/Users/arnau.gomez/Desktop/GofCens_Paper_Simulations/04_Output/02_Plots/New_plots/Figure3.1_Lognoral_Power_complete_data.png", 
       plot = p3, width = 2130, height = 1600, units = "px")

p3 <- p3 + labs(title=NULL)

(combined_complete_data <- (p2 + p3 + plot_layout(guides = "collect") & theme(legend.position = "bottom", 
                                                                              legend.title = element_text(size = 20, face = "bold"),
                                                                              legend.text = element_text(size = 20),
                                                                              plot.subtitle = element_text(size = 20, face = "italic", hjust = 0.5),
                                                                              axis.title = element_text(size = 15, face = "bold"),
                                                                              axis.text = element_text(size = 15))) +
    plot_annotation(title = NULL, 
                    theme =  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))))

ggsave(filename = "04_Output/Figures/Power_AdminCensoring_Data_Null_0.7.png", plot = combined_complete_data, width = 2130*2, height = 1800, units = "px")




################################################################################
#-------------------------------------------------------------------------------
# Power comparison for random censored data by Null Hypothesis
#-------------------------------------------------------------------------------
################################################################################

library(ggplot2)
library(dplyr)
library(patchwork)

#----------------------------------------------------------------------------------------------------------------------------
# X ~ Weibull
#----------------------------------------------------------------------------------------------------------------------------

df_W_logistic <- get(load("02_Data/Weibull_RandomCensoring/df_power_all4_logistic/Power/_Distribution_weibull_resultados_h0_logistic_all_.RData"))
df_W_log <- get(load("02_Data/Weibull_RandomCensoring/df_power_all4_lognormal/Power/_Distribution_weibull_resultados_h0_lognormal_all_.RData"))
df_W_wei <- get(load("02_Data/Weibull_RandomCensoring/df_power_all4_weibull/Power/_Distribution_weibull_resultados_h0_weibull_all_.RData"))


df_W_logistic2 <- df_W_logistic %>%
  select(power, cens, test) %>%
  filter(cens == "60% Censura") %>%
  mutate("Hyp_nul" = rep("Logistic"), cens = NULL)

df_W_log2 <- df_W_log %>%
  select(power, cens, test) %>%
  filter(cens == "60% Censura") %>%
  mutate("Hyp_nul" = rep("Lognormal"), cens = NULL)

df_W_wei2 <- df_W_wei %>%
  select(power, cens, test) %>%
  filter(cens == "60% Censura") %>%
  mutate("Hyp_nul" = rep("Weibull"), cens = NULL)

df_complet <- rbind(df_W_logistic2, df_W_log2, df_W_wei2 )%>%
  filter(test != "KS No Bootstrap")


(p2 <- ggplot(df_complet, aes(x = test, y = power, shape = Hyp_nul, fill = Hyp_nul)) + geom_jitter(size = 3, stroke = 0.5, height = 0, width = 0.1) +
    scale_shape_manual(values = c(22,24,23)) +
    scale_fill_manual(values = c( "Logistic" = "#CAFF70", "Lognormal" = "#79CDCD","Weibull" = "#FF6A6A")) +
    ggtitle(paste0("Power of each Tests for Random Censored Data by Null Hypothesis"), subtitle = paste0("X ~ Weibull. 60% Censoring")) +
    labs(x = "Test", y = "Power", fill = "Null Hypothesis:") +
    theme_bw() +
    theme(legend.position = "right", ) +
    lims(y = c(-0.01,1.01)) +
    guides(fill = guide_legend(override.aes = list(shape = c(22,24,23), color = "black")),
           shape = "none") +
    geom_hline(yintercept = 0.1, colour = "darkblue", linetype = 2, linewidth = 0.9) +
    annotate("text", x = 0.5, y = 0.07, label = "alpha = 0.1", hjust = -0.02, color = "darkblue", 
             size = 3.5, angle = 330))

ggsave(filename = "C:/Users/arnau.gomez/Desktop/GofCens_Paper_Simulations/04_Output/02_Plots/New_plots/Figure3.1_Weibull_Power_complete_data.png", 
       plot = p2, width = 2130, height = 1600, units = "px")

p2 <- p2 + labs(title = NULL)

#----------------------------------------------------------------------------------------------------------------------------
# X ~ Logistic
#----------------------------------------------------------------------------------------------------------------------------

df_L_logistic <- get(load("02_Data/Lognormal_RandomCensoring/df_power_all4_logistic/Power/_Distribution_lognormal_resultados_h0_logistic_all_.RData"))
df_L_log <- get(load("02_Data/Lognormal_RandomCensoring/df_power_all4_lognormal/Power/_Distribution_lognormal_resultados_h0_lognormal_all_.RData"))
df_L_wei <- get(load("02_Data/Lognormal_RandomCensoring/df_power_all4_weibull/Power/_Distribution_lognormal_resultados_h0_weibull_all_.RData"))


df_L_logistic2 <- df_L_logistic %>%
  select(power, cens, test) %>%
  filter(cens == "60% Censura") %>%
  mutate("Hyp_nul" = rep("Logistic"), cens = NULL)

df_L_log2 <- df_L_log %>%
  select(power, cens, test) %>%
  filter(cens == "60% Censura") %>%
  mutate("Hyp_nul" = rep("Lognormal"), cens = NULL)

df_L_wei2 <- df_L_wei %>%
  select(power, cens, test) %>%
  filter(cens == "60% Censura") %>%
  mutate("Hyp_nul" = rep("Weibull"), cens = NULL)

df_complet <- rbind(df_L_logistic2, df_L_log2, df_L_wei2) %>%
  filter(test != "KS No Bootstrap")


(p3 <- ggplot(df_complet, aes(x = test, y = power, shape = Hyp_nul, fill = Hyp_nul)) + geom_jitter(size = 3, stroke = 0.5, height = 0, width = 0.1) +
    scale_shape_manual(values = c(22,24,23)) +
    scale_fill_manual(values = c( "Logistic" = "#CAFF70", "Lognormal" = "#79CDCD","Weibull" = "#FF6A6A")) +
    ggtitle(paste0("Power of each Tests for Random Censored Data by Null Hypothesis"), subtitle = paste0("X ~ Lognormal. 60% Censoring")) +
    labs(x = "Test", y = "Power", fill = "Null Hypothesis:") +
    theme_bw() +
    theme(legend.position = "right", ) +
    lims(y = c(-0.01,1.01)) +
    guides(fill = guide_legend(override.aes = list(shape = c(22,24,23), color = "black")),
           shape = "none") +
    geom_hline(yintercept = 0.1, colour = "darkblue", linetype = 2, linewidth = 0.9) +
    annotate("text", x = 0.5, y = 0.07, label = "alpha = 0.1", hjust = -0.02, color = "darkblue", 
             size = 3.5, angle = 330))

ggsave(filename = "C:/Users/arnau.gomez/Desktop/GofCens_Paper_Simulations/04_Output/02_Plots/New_plots/Figure3.1_Lognoral_Power_complete_data.png", 
       plot = p3, width = 2130, height = 1600, units = "px")

p3 <- p3 + labs(title=NULL)

(combined_complete_data <- (p2 + p3 + plot_layout(guides = "collect") & theme(legend.position = "bottom", 
                                                                              legend.title = element_text(size = 20, face = "bold"),
                                                                              legend.text = element_text(size = 20),
                                                                              plot.subtitle = element_text(size = 20, face = "italic", hjust = 0.5),
                                                                              axis.title = element_text(size = 15, face = "bold"),
                                                                              axis.text = element_text(size = 15))) +
    plot_annotation(title = NULL, 
                    theme =  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))))

ggsave(filename = "04_Output/Figures/Power_RandomCensoring_Data_Null_0.4.png", plot = combined_complete_data, width = 2130*2, height = 1800, units = "px")

################################################################################
#-------------------------------------------------------------------------------
# Power comparison for administrative censored data by Skewness fo Tests = KS
#-------------------------------------------------------------------------------
################################################################################

library(ggplot2)
library(dplyr)
library(patchwork)

#----------------------------------------------------------------------------------------------------------------------------
# X ~ Weibull
#----------------------------------------------------------------------------------------------------------------------------

df_W_logistic <- get(load("02_Data/Weibull_AdminCensoring/df_power_all4_logistic/Power/_Distribution_weibull_resultados_h0_logistic_all_.RData"))
df_W_log <- get(load("02_Data/Weibull_AdminCensoring/df_power_all4_lognormal/Power/_Distribution_weibull_resultados_h0_lognormal_all_.RData"))
df_W_wei <- get(load("02_Data/Weibull_AdminCensoring/df_power_all4_weibull/Power/_Distribution_weibull_resultados_h0_weibull_all_.RData"))

df_KS <- df_W_logistic %>%
  select(power, test, symmetry, cens) %>%
  filter(test == "Kolmogorov-Smirnov" & cens != "0% Censura") %>%
  mutate(Distribution = "Logistic") %>%
  
  rbind(
    df_W_log %>%
      select(power, test, symmetry, cens) %>%
      filter(test == "Kolmogorov-Smirnov" & cens != "0% Censura") %>%
      mutate(Distribution = "Lognormal")
  ) %>%
  
  rbind(
    df_W_wei %>%
      select(power, test, symmetry, cens) %>%
      filter(test == "Kolmogorov-Smirnov" & cens != "0% Censura") %>%
      mutate(Distribution = "Weibull")
  ) %>%
  mutate(symmetry = factor(symmetry,
                           levels = c("Symmetric",
                                      "Asymmetric",
                                      "Extremly asymmetric")))

(p_KS_W <- ggplot(df_KS, aes(x = factor(symmetry), y = power, shape = Distribution, fill = Distribution)) +
    geom_jitter(size = 3, stroke = 0.5, height = 0, width = 0.1) +
    scale_shape_manual(values = c(22,24,23)) +
    scale_fill_manual(values = c("Logistic" = "#CAFF70", "Lognormal" = "#79CDCD", "Weibull" = "#FF6A6A")) +
    ggtitle("Power of Kolmogorov-Smirnov test under different skewness levels",
            subtitle = "X ~ Weibull") +
    labs(x = "Skewness level", y = "Power", fill = "Null Hypotheses:", shape = "Null Hypotheses:") +
    theme_bw() +
    theme(legend.position = "right") +
    lims(y = c(-0.01,1.01)) +
    guides(fill = guide_legend(override.aes = list( shape = c(22,24,23), color = "black")), shape = "none") +
    geom_hline(yintercept = 0.1, colour = "darkblue", linetype = 2, linewidth = 0.9))

p_KS_W <- p_KS_W + labs(title = NULL)


#----------------------------------------------------------------------------------------------------------------------------
# X ~ Lognormal
#----------------------------------------------------------------------------------------------------------------------------

df_L_logistic <- get(load("02_Data/Lognormal_AdminCensoring/df_power_all4_logistic/Power/_Distribution_lognormal_resultados_h0_logistic_all_.RData"))
df_L_log <- get(load("02_Data/Lognormal_AdminCensoring/df_power_all4_lognormal/Power/_Distribution_lognormal_resultados_h0_lognormal_all_.RData"))
df_L_wei <- get(load("02_Data/Lognormal_AdminCensoring/df_power_all4_weibull/Power/_Distribution_lognormal_resultados_h0_weibull_all_.RData"))

df_KS <- df_L_logistic %>%
  select(power, test, symmetry, cens) %>%
  filter(test == "Kolmogorov-Smirnov" & cens != "0% Censura") %>%
  mutate(Distribution = "Logistic") %>%
  
  rbind(
    df_L_log %>%
      select(power, test, symmetry, cens) %>%
      filter(test == "Kolmogorov-Smirnov" & cens != "0% Censura") %>%
      mutate(Distribution = "Lognormal")
  ) %>%
  
  rbind(
    df_L_wei %>%
      select(power, test, symmetry, cens) %>%
      filter(test == "Kolmogorov-Smirnov" & cens != "0% Censura") %>%
      mutate(Distribution = "Weibull")
  ) %>%
  mutate(symmetry = factor(symmetry,
                           levels = c("Symmetric",
                                      "Asymmetric",
                                      "Extremly asymmetric")))

(p_KS_L <- ggplot(df_KS, aes(x = factor(symmetry), y = power, shape = Distribution, fill = Distribution)) +
    geom_jitter(size = 3, stroke = 0.5, height = 0, width = 0.1) +
    scale_shape_manual(values = c(22,24,23)) +
    scale_fill_manual(values = c("Logistic" = "#CAFF70", "Lognormal" = "#79CDCD", "Weibull" = "#FF6A6A")) +
    ggtitle("Power of Kolmogorov-Smirnov test under different skewness levels",
            subtitle = "X ~ Lognormal") +
    labs(x = "Skewness level", y = "Power", fill = "Null Hypotheses:", shape = "Null Hypotheses:") +
    theme_bw() +
    theme(legend.position = "right") +
    lims(y = c(-0.01,1.01)) +
    guides(fill = guide_legend(override.aes = list( shape = c(22,24,23), color = "black")), shape = "none") +
    geom_hline(yintercept = 0.1, colour = "darkblue", linetype = 2, linewidth = 0.9))

p_KS_L <- p_KS_L + labs(title = NULL)

(p_KS_combined <- (p_KS_W + p_KS_L + plot_layout(guides = "collect") & theme(legend.position = "bottom", 
                                                                             legend.title = element_text(size = 20, face = "bold"),
                                                                             legend.text = element_text(size = 20),
                                                                             plot.subtitle = element_text(size = 20, face = "italic", hjust = 0.5),
                                                                             axis.title = element_text(size = 15, face = "bold"),
                                                                             axis.text = element_text(size = 15))) +
    plot_annotation(title = NULL, 
                    theme =  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))))


ggsave(filename = "04_Output/Figures/Power_Kolmogorov_Skewness_AdminCensoring.png", plot = p_KS_combined, width = 2130*2, height = 1800, units = "px")


################################################################################
#-------------------------------------------------------------------------------
# Power comparison for random censored data by Skewness fo Tests = AD
#-------------------------------------------------------------------------------
################################################################################

library(ggplot2)
library(dplyr)
library(patchwork)

#----------------------------------------------------------------------------------------------------------------------------
# X ~ Weibull
#----------------------------------------------------------------------------------------------------------------------------

df_W_logistic <- get(load("02_Data/Weibull_RandomCensoring/df_power_all4_logistic/Power/_Distribution_weibull_resultados_h0_logistic_all_.RData"))
df_W_log <- get(load("02_Data/Weibull_RandomCensoring/df_power_all4_lognormal/Power/_Distribution_weibull_resultados_h0_lognormal_all_.RData"))
df_W_wei <- get(load("02_Data/Weibull_RandomCensoring/df_power_all4_weibull/Power/_Distribution_weibull_resultados_h0_weibull_all_.RData"))

df_KS <- df_W_logistic %>%
  select(power, test, symmetry, cens) %>%
  filter(test == "Anderson-Darling" & cens != "0% Censura") %>%
  mutate(Distribution = "Logistic") %>%
  
  rbind(
    df_W_log %>%
      select(power, test, symmetry, cens) %>%
      filter(test == "Anderson-Darling" & cens != "0% Censura") %>%
      mutate(Distribution = "Lognormal")
  ) %>%
  
  rbind(
    df_W_wei %>%
      select(power, test, symmetry, cens) %>%
      filter(test == "Anderson-Darling" & cens != "0% Censura") %>%
      mutate(Distribution = "Weibull")
  ) %>%
  mutate(symmetry = factor(symmetry,
                           levels = c("Symmetric",
                                      "Asymmetric",
                                      "Extremly asymmetric")))

(p_KS_W <- ggplot(df_KS, aes(x = factor(symmetry), y = power, shape = Distribution, fill = Distribution)) +
    geom_jitter(size = 3, stroke = 0.5, height = 0, width = 0.1) +
    scale_shape_manual(values = c(22,24,23)) +
    scale_fill_manual(values = c("Logistic" = "#CAFF70", "Lognormal" = "#79CDCD", "Weibull" = "#FF6A6A")) +
    ggtitle("Power of Anderson-Darling test under different skewness levels",
            subtitle = "X ~ Weibull") +
    labs(x = "Skewness level", y = "Power", fill = "Null Hypotheses:", shape = "Null Hypotheses:") +
    theme_bw() +
    theme(legend.position = "right") +
    lims(y = c(-0.01,1.01)) +
    guides(fill = guide_legend(override.aes = list( shape = c(22,24,23), color = "black")), shape = "none") +
    geom_hline(yintercept = 0.1, colour = "darkblue", linetype = 2, linewidth = 0.9))

p_KS_W <- p_KS_W + labs(title = NULL)


#----------------------------------------------------------------------------------------------------------------------------
# X ~ Lognormal
#----------------------------------------------------------------------------------------------------------------------------

df_L_logistic <- get(load("02_Data/Lognormal_RandomCensoring/df_power_all4_logistic/Power/_Distribution_lognormal_resultados_h0_logistic_all_.RData"))
df_L_log <- get(load("02_Data/Lognormal_RandomCensoring/df_power_all4_lognormal/Power/_Distribution_lognormal_resultados_h0_lognormal_all_.RData"))
df_L_wei <- get(load("02_Data/Lognormal_RandomCensoring/df_power_all4_weibull/Power/_Distribution_lognormal_resultados_h0_weibull_all_.RData"))

df_KS <- df_L_logistic %>%
  select(power, test, symmetry, cens) %>%
  filter(test == "Anderson-Darling" & cens != "0% Censura") %>%
  mutate(Distribution = "Logistic") %>%
  
  rbind(
    df_L_log %>%
      select(power, test, symmetry, cens) %>%
      filter(test == "Anderson-Darling" & cens != "0% Censura") %>%
      mutate(Distribution = "Lognormal")
  ) %>%
  
  rbind(
    df_L_wei %>%
      select(power, test, symmetry, cens) %>%
      filter(test == "Anderson-Darling" & cens != "0% Censura") %>%
      mutate(Distribution = "Weibull")
  ) %>%
  mutate(symmetry = factor(symmetry,
                           levels = c("Symmetric",
                                      "Asymmetric",
                                      "Extremly asymmetric")))

(p_KS_L <- ggplot(df_KS, aes(x = factor(symmetry), y = power, shape = Distribution, fill = Distribution)) +
    geom_jitter(size = 3, stroke = 0.5, height = 0, width = 0.1) +
    scale_shape_manual(values = c(22,24,23)) +
    scale_fill_manual(values = c("Logistic" = "#CAFF70", "Lognormal" = "#79CDCD", "Weibull" = "#FF6A6A")) +
    ggtitle("Power of Anderson-Darling test under different skewness levels",
            subtitle = "X ~ Lognormal") +
    labs(x = "Skewness level", y = "Power", fill = "Null Hypotheses:", shape = "Null Hypotheses:") +
    theme_bw() +
    theme(legend.position = "right") +
    lims(y = c(-0.01,1.01)) +
    guides(fill = guide_legend(override.aes = list( shape = c(22,24,23), color = "black")), shape = "none") +
    geom_hline(yintercept = 0.1, colour = "darkblue", linetype = 2, linewidth = 0.9))

p_KS_L <- p_KS_L + labs(title = NULL)

(p_AD_combined <- (p_KS_W + p_KS_L + plot_layout(guides = "collect") & theme(legend.position = "bottom", 
                                                                             legend.title = element_text(size = 20, face = "bold"),
                                                                             legend.text = element_text(size = 20),
                                                                             plot.subtitle = element_text(size = 20, face = "italic", hjust = 0.5),
                                                                             axis.title = element_text(size = 15, face = "bold"),
                                                                             axis.text = element_text(size = 15))) +
    plot_annotation(title = "Anderson-Darling test", 
                    theme =  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))))


ggsave(filename = "04_Output/Figures/Power_Anderson_Darling_Skewness_RandomCensoring.png", plot = p_AD_combined, width = 2130*2, height = 1800, units = "px")


################################################################################
#-------------------------------------------------------------------------------
# Power comparison for random censored data by Skewness fo Tests = CvM
#-------------------------------------------------------------------------------
################################################################################

library(ggplot2)
library(dplyr)
library(patchwork)

#----------------------------------------------------------------------------------------------------------------------------
# X ~ Weibull
#----------------------------------------------------------------------------------------------------------------------------

df_W_logistic <- get(load("02_Data/Weibull_RandomCensoring/df_power_all4_logistic/Power/_Distribution_weibull_resultados_h0_logistic_all_.RData"))
df_W_log <- get(load("02_Data/Weibull_RandomCensoring/df_power_all4_lognormal/Power/_Distribution_weibull_resultados_h0_lognormal_all_.RData"))
df_W_wei <- get(load("02_Data/Weibull_RandomCensoring/df_power_all4_weibull/Power/_Distribution_weibull_resultados_h0_weibull_all_.RData"))

df_KS <- df_W_logistic %>%
  select(power, test, symmetry, cens) %>%
  filter(test == "Cramér von Mises" & cens != "0% Censura") %>%
  mutate(Distribution = "Logistic") %>%
  
  rbind(
    df_W_log %>%
      select(power, test, symmetry, cens) %>%
      filter(test == "Cramér von Mises" & cens != "0% Censura") %>%
      mutate(Distribution = "Lognormal")
  ) %>%
  
  rbind(
    df_W_wei %>%
      select(power, test, symmetry, cens) %>%
      filter(test == "Cramér von Mises" & cens != "0% Censura") %>%
      mutate(Distribution = "Weibull")
  ) %>%
  mutate(symmetry = factor(symmetry,
                           levels = c("Symmetric",
                                      "Asymmetric",
                                      "Extremly asymmetric")))

(p_KS_W <- ggplot(df_KS, aes(x = factor(symmetry), y = power, shape = Distribution, fill = Distribution)) +
    geom_jitter(size = 3, stroke = 0.5, height = 0, width = 0.1) +
    scale_shape_manual(values = c(22,24,23)) +
    scale_fill_manual(values = c("Logistic" = "#CAFF70", "Lognormal" = "#79CDCD", "Weibull" = "#FF6A6A")) +
    ggtitle("Power of Cramér-von Mises test under different skewness levels",
            subtitle = "X ~ Weibull") +
    labs(x = "Skewness level", y = "Power", fill = "Null Hypotheses:", shape = "Null Hypotheses:") +
    theme_bw() +
    theme(legend.position = "right") +
    lims(y = c(-0.01,1.01)) +
    guides(fill = guide_legend(override.aes = list( shape = c(22,24,23), color = "black")), shape = "none") +
    geom_hline(yintercept = 0.1, colour = "darkblue", linetype = 2, linewidth = 0.9))

p_KS_W <- p_KS_W + labs(title = NULL)


#----------------------------------------------------------------------------------------------------------------------------
# X ~ Lognormal
#----------------------------------------------------------------------------------------------------------------------------

df_L_logistic <- get(load("02_Data/Lognormal_RandomCensoring/df_power_all4_logistic/Power/_Distribution_lognormal_resultados_h0_logistic_all_.RData"))
df_L_log <- get(load("02_Data/Lognormal_RandomCensoring/df_power_all4_lognormal/Power/_Distribution_lognormal_resultados_h0_lognormal_all_.RData"))
df_L_wei <- get(load("02_Data/Lognormal_RandomCensoring/df_power_all4_weibull/Power/_Distribution_lognormal_resultados_h0_weibull_all_.RData"))

df_KS <- df_L_logistic %>%
  select(power, test, symmetry, cens) %>%
  filter(test == "Cramér von Mises" & cens != "0% Censura") %>%
  mutate(Distribution = "Logistic") %>%
  
  rbind(
    df_L_log %>%
      select(power, test, symmetry, cens) %>%
      filter(test == "Cramér von Mises" & cens != "0% Censura") %>%
      mutate(Distribution = "Lognormal")
  ) %>%
  
  rbind(
    df_L_wei %>%
      select(power, test, symmetry, cens) %>%
      filter(test == "Cramér von Mises" & cens != "0% Censura") %>%
      mutate(Distribution = "Weibull")
  ) %>%
  mutate(symmetry = factor(symmetry,
                           levels = c("Symmetric",
                                      "Asymmetric",
                                      "Extremly asymmetric")))

(p_KS_L <- ggplot(df_KS, aes(x = factor(symmetry), y = power, shape = Distribution, fill = Distribution)) +
    geom_jitter(size = 3, stroke = 0.5, height = 0, width = 0.1) +
    scale_shape_manual(values = c(22,24,23)) +
    scale_fill_manual(values = c("Logistic" = "#CAFF70", "Lognormal" = "#79CDCD", "Weibull" = "#FF6A6A")) +
    ggtitle("Power of Cramér-von Mises test under different skewness levels",
            subtitle = "X ~ Lognormal") +
    labs(x = "Skewness level", y = "Power", fill = "Null Hypotheses:", shape = "Null Hypotheses:") +
    theme_bw() +
    theme(legend.position = "right") +
    lims(y = c(-0.01,1.01)) +
    guides(fill = guide_legend(override.aes = list( shape = c(22,24,23), color = "black")), shape = "none") +
    geom_hline(yintercept = 0.1, colour = "darkblue", linetype = 2, linewidth = 0.9))

p_KS_L <- p_KS_L + labs(title = NULL)

(p_CvM_combined <- (p_KS_W + p_KS_L + plot_layout(guides = "collect") & theme(legend.position = "bottom", 
                                                                             legend.title = element_text(size = 20, face = "bold"),
                                                                             legend.text = element_text(size = 20),
                                                                             plot.subtitle = element_text(size = 20, face = "italic", hjust = 0.5),
                                                                             axis.title = element_text(size = 15, face = "bold"),
                                                                             axis.text = element_text(size = 15))) +
    plot_annotation(title = "Cramér-von Mises test", 
                    theme =  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))))


ggsave(filename = "04_Output/Figures/Power_Cramér_von_Mises_Skewness_RandomCensoring.png", plot = p_CvM_combined, width = 2130*2, height = 1800, units = "px")


################################################################################
#-------------------------------------------------------------------------------
# Power comparison for random censored data by Skewness fo Tests = KS
#-------------------------------------------------------------------------------
################################################################################

library(ggplot2)
library(dplyr)
library(patchwork)

#----------------------------------------------------------------------------------------------------------------------------
# X ~ Weibull
#----------------------------------------------------------------------------------------------------------------------------

df_W_logistic <- get(load("02_Data/Weibull_RandomCensoring/df_power_all4_logistic/Power/_Distribution_weibull_resultados_h0_logistic_all_.RData"))
df_W_log <- get(load("02_Data/Weibull_RandomCensoring/df_power_all4_lognormal/Power/_Distribution_weibull_resultados_h0_lognormal_all_.RData"))
df_W_wei <- get(load("02_Data/Weibull_RandomCensoring/df_power_all4_weibull/Power/_Distribution_weibull_resultados_h0_weibull_all_.RData"))

df_KS <- df_W_logistic %>%
  select(power, test, symmetry, cens) %>%
  filter(test == "Kolmogorov-Smirnov" & cens != "0% Censura") %>%
  mutate(Distribution = "Logistic") %>%
  
  rbind(
    df_W_log %>%
      select(power, test, symmetry, cens) %>%
      filter(test == "Kolmogorov-Smirnov" & cens != "0% Censura") %>%
      mutate(Distribution = "Lognormal")
  ) %>%
  
  rbind(
    df_W_wei %>%
      select(power, test, symmetry, cens) %>%
      filter(test == "Kolmogorov-Smirnov" & cens != "0% Censura") %>%
      mutate(Distribution = "Weibull")
  ) %>%
  mutate(symmetry = factor(symmetry,
                           levels = c("Symmetric",
                                      "Asymmetric",
                                      "Extremly asymmetric")))

(p_KS_W <- ggplot(df_KS, aes(x = factor(symmetry), y = power, shape = Distribution, fill = Distribution)) +
    geom_jitter(size = 3, stroke = 0.5, height = 0, width = 0.1) +
    scale_shape_manual(values = c(22,24,23)) +
    scale_fill_manual(values = c("Logistic" = "#CAFF70", "Lognormal" = "#79CDCD", "Weibull" = "#FF6A6A")) +
    ggtitle("Power of Kolmogorov-Smirnov test under different skewness levels",
            subtitle = "X ~ Weibull") +
    labs(x = "Skewness level", y = "Power", fill = "Null Hypotheses:", shape = "Null Hypotheses:") +
    theme_bw() +
    theme(legend.position = "right") +
    lims(y = c(-0.01,1.01)) +
    guides(fill = guide_legend(override.aes = list( shape = c(22,24,23), color = "black")), shape = "none") +
    geom_hline(yintercept = 0.1, colour = "darkblue", linetype = 2, linewidth = 0.9))

p_KS_W <- p_KS_W + labs(title = NULL)


#----------------------------------------------------------------------------------------------------------------------------
# X ~ Lognormal
#----------------------------------------------------------------------------------------------------------------------------

df_L_logistic <- get(load("02_Data/Lognormal_RandomCensoring/df_power_all4_logistic/Power/_Distribution_lognormal_resultados_h0_logistic_all_.RData"))
df_L_log <- get(load("02_Data/Lognormal_RandomCensoring/df_power_all4_lognormal/Power/_Distribution_lognormal_resultados_h0_lognormal_all_.RData"))
df_L_wei <- get(load("02_Data/Lognormal_RandomCensoring/df_power_all4_weibull/Power/_Distribution_lognormal_resultados_h0_weibull_all_.RData"))

df_KS <- df_L_logistic %>%
  select(power, test, symmetry, cens) %>%
  filter(test == "Kolmogorov-Smirnov" & cens != "0% Censura") %>%
  mutate(Distribution = "Logistic") %>%
  
  rbind(
    df_L_log %>%
      select(power, test, symmetry, cens) %>%
      filter(test == "Kolmogorov-Smirnov" & cens != "0% Censura") %>%
      mutate(Distribution = "Lognormal")
  ) %>%
  
  rbind(
    df_L_wei %>%
      select(power, test, symmetry, cens) %>%
      filter(test == "Kolmogorov-Smirnov" & cens != "0% Censura") %>%
      mutate(Distribution = "Weibull")
  ) %>%
  mutate(symmetry = factor(symmetry,
                           levels = c("Symmetric",
                                      "Asymmetric",
                                      "Extremly asymmetric")))

(p_KS_L <- ggplot(df_KS, aes(x = factor(symmetry), y = power, shape = Distribution, fill = Distribution)) +
    geom_jitter(size = 3, stroke = 0.5, height = 0, width = 0.1) +
    scale_shape_manual(values = c(22,24,23)) +
    scale_fill_manual(values = c("Logistic" = "#CAFF70", "Lognormal" = "#79CDCD", "Weibull" = "#FF6A6A")) +
    ggtitle("Power of Kolmogorov-Smirnov test under different skewness levels",
            subtitle = "X ~ Lognormal") +
    labs(x = "Skewness level", y = "Power", fill = "Null Hypotheses:", shape = "Null Hypotheses:") +
    theme_bw() +
    theme(legend.position = "right") +
    lims(y = c(-0.01,1.01)) +
    guides(fill = guide_legend(override.aes = list( shape = c(22,24,23), color = "black")), shape = "none") +
    geom_hline(yintercept = 0.1, colour = "darkblue", linetype = 2, linewidth = 0.9))

p_KS_L <- p_KS_L + labs(title = NULL)

(p_KS_combined <- (p_KS_W + p_KS_L + plot_layout(guides = "collect") & theme(legend.position = "bottom", 
                                                                             legend.title = element_text(size = 20, face = "bold"),
                                                                             legend.text = element_text(size = 20),
                                                                             plot.subtitle = element_text(size = 20, face = "italic", hjust = 0.5),
                                                                             axis.title = element_text(size = 15, face = "bold"),
                                                                             axis.text = element_text(size = 15))) +
    plot_annotation(title = "Kolmogorov-Smirnov test", 
                    theme =  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))))


ggsave(filename = "04_Output/Figures/Power_Kolmogorov_Skewness_RandomCensoring.png", plot = p_KS_combined, width = 2130*2, height = 1800, units = "px")


p_AD_combined <- wrap_elements(
  full = p_AD_combined
) + plot_annotation(title = "Anderson-Darling test", theme = theme(plot.title = element_text(size = 50, face = "bold", hjust = 0.5)))

p_CvM_combined <- wrap_elements(
  full = p_CvM_combined
) + plot_annotation(title = "Cramér-von Mises test", theme = theme(plot.title = element_text(size = 50, face = "bold", hjust = 0.5)))

p_KS_combined <- wrap_elements(
  full = p_KS_combined
) + plot_annotation(title = "Kolmogorov-Smirnov test", theme = theme(plot.title = element_text(size = 50, face = "bold", hjust = 0.5)))

(p_final <- (p_AD_combined / p_CvM_combined / p_KS_combined + plot_layout(guides = "collect") & theme(legend.position = "bottom", 
                                                                                                           legend.title = element_text(size = 30, face = "bold"),
                                                                                                           legend.text = element_text(size = 30),
                                                                                                           plot.title = element_text(size = 30, face = "italic", hjust = 0.5),
                                                                                                           plot.subtitle = element_text(size = 30, face = "italic", hjust = 0.5),
                                                                                                           axis.title = element_text(size = 15, face = "bold"),
                                                                                                           axis.text = element_text(size = 15))) +
    plot_annotation(title = NULL, 
                    theme =  theme(plot.title = element_text(size = 35, face = "bold", hjust = 0.5))))



ggsave(filename = "04_Output/Figures/Power_Skewness_RandomCensoring.png", plot = p_final, width = 2130*2, height = 3*1800, units = "px")


################################################################################
#-------------------------------------------------------------------------------
# Diff against censoring level at random censoring
#-------------------------------------------------------------------------------
################################################################################

library(ggplot2)
library(patchwork)
library(dplyr)

################################################################################
# Weibull
################################################################################

df_W_logistic <- get(load("02_Data/Weibull_RandomCensoring/df_power_all4_logistic/Power/_Distribution_weibull_resultados_h0_logistic_all_.RData"))
df_W_log <- get(load("02_Data/Weibull_RandomCensoring/df_power_all4_lognormal/Power/_Distribution_weibull_resultados_h0_lognormal_all_.RData"))
df_W_wei <- get(load("02_Data/Weibull_RandomCensoring/df_power_all4_weibull/Power/_Distribution_weibull_resultados_h0_weibull_all_.RData"))

df_W_logistic2 <- df_W_logistic %>%
  mutate(
    Diff_tmax = SEC_min_median - km_min_median,
    symmetry = factor(
      symmetry,
      levels = c("Symmetric", "Asymmetric", "Extremly asymmetric"),
      labels = c("Sym", "Asym", "Ex. Asym")
    ), cens = factor(cens, levels = c("0% Censura", "30% Censura", "60% Censura"),
                     labels = c("Complete Data", "30% Censorship", "60% Censorship"))
  )

df_W_log2 <- df_W_log %>%
  mutate(
    Diff_tmax = SEC_min_median - km_min_median,
    symmetry = factor(
      symmetry,
      levels = c("Symmetric", "Asymmetric", "Extremly asymmetric"),
      labels = c("Sym", "Asym", "Ex. Asym")
    ), cens = factor(cens, levels = c("0% Censura", "30% Censura", "60% Censura"),
                     labels = c("Complete Data", "30% Censorship", "60% Censorship"))
  )

df_W_wei2 <- df_W_wei %>%
  mutate(
    Diff_tmax = SEC_min_median - km_min_median,
    symmetry = factor(
      symmetry,
      levels = c("Symmetric", "Asymmetric", "Extremly asymmetric"),
      labels = c("Sym", "Asym", "Ex. Asym")
    ), cens = factor(cens, levels = c("0% Censura", "30% Censura", "60% Censura"),
                     labels = c("Complete Data", "30% Censorship", "60% Censorship"))
  )

df_plot <- rbind(df_W_logistic2, df_W_log2, df_W_wei2) %>%
  filter(test != "KS No Bootstrap")


(p_Diff_W <- ggplot(df_plot, aes(x = symmetry, y = Diff_tmax, color = cens)) +
  geom_jitter(width = 0.15, size = 3, shape = 16) +
  facet_wrap(~ cens) +
  scale_color_manual(values = c("Complete Data" = "#CAFF70", "30% Censorship" = "#79CDCD", "60% Censorship" = "#FF6A6A")) +
  theme_minimal() +
  labs(x = "Skewness level", y = expression(Diff[t[max]]), color = "Censoring level",
       title = expression(Diff[t[max]] ~ "vs Skewness under Random Censoring"),
       subtitle = "X ~ Weibull") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1)) + ylim(c(-0.2,0.2)))

p_Diff_W <- p_Diff_W + labs(title = NULL)


################################################################################
# Lognormal
################################################################################

df_L_logistic <- get(load("02_Data/Lognormal_RandomCensoring/df_power_all4_logistic/Power/_Distribution_lognormal_resultados_h0_logistic_all_.RData"))
df_L_log <- get(load("02_Data/Lognormal_RandomCensoring/df_power_all4_lognormal/Power/_Distribution_lognormal_resultados_h0_lognormal_all_.RData"))
df_L_wei <- get(load("02_Data/Lognormal_RandomCensoring/df_power_all4_weibull/Power/_Distribution_lognormal_resultados_h0_weibull_all_.RData"))


df_L_logistic2 <- df_L_logistic %>%
  mutate(
    Diff_tmax = SEC_min_median - km_min_median,
    symmetry = factor(
      symmetry,
      levels = c("Symmetric", "Asymmetric", "Extremly asymmetric"),
      labels = c("Sym", "Asym", "Ex. Asym")
    ), cens = factor(cens, levels = c("0% Censura", "30% Censura", "60% Censura"),
                     labels = c("Complete Data", "30% Censorship", "60% Censorship"))
  )

df_L_log2 <- df_L_log %>%
  mutate(
    Diff_tmax = SEC_min_median - km_min_median,
    symmetry = factor(
      symmetry,
      levels = c("Symmetric", "Asymmetric", "Extremly asymmetric"),
      labels = c("Sym", "Asym", "Ex. Asym")
    ), cens = factor(cens, levels = c("0% Censura", "30% Censura", "60% Censura"),
                     labels = c("Complete Data", "30% Censorship", "60% Censorship"))
  )

df_L_wei2 <- df_L_wei %>%
  mutate(
    Diff_tmax = SEC_min_median - km_min_median,
    symmetry = factor(
      symmetry,
      levels = c("Symmetric", "Asymmetric", "Extremly asymmetric"),
      labels = c("Sym", "Asym", "Ex. Asym")
    ), cens = factor(cens, levels = c("0% Censura", "30% Censura", "60% Censura"),
                     labels = c("Complete Data", "30% Censorship", "60% Censorship"))
  )

df_plot <- rbind(df_L_logistic2, df_L_log2, df_L_wei2) %>%
  filter(test != "KS No Bootstrap")


(p_Diff_L <- ggplot(df_L_logistic2, aes(x = symmetry, y = Diff_tmax, color = cens)) +
    geom_jitter(width = 0.15, size = 3, shape = 16) +
    facet_wrap(~ cens) +
    scale_color_manual(values = c("Complete Data" = "#CAFF70", "30% Censorship" = "#79CDCD", "60% Censorship" = "#FF6A6A")) +
    theme_minimal() +
    labs(x = "Skewness level", y = expression(Diff[t[max]]), color = "Censoring level",
         title = expression(Diff[t[max]] ~ "vs Skewness under Random Censoring"),
         subtitle = "X ~ Lognormal") +
    theme(axis.text.x = element_text(angle = 25, hjust = 1)) + ylim(c(-0.2,0.2)))

p_Diff_L <- p_Diff_L + labs(title = NULL)


(p_Diff_combined <- (p_Diff_W + p_Diff_L + plot_layout(guides = "collect")&
                       theme_minimal()  & theme(legend.position = "bottom",
                                                legend.title = element_text(size = 20, face = "bold"),
                                                legend.text = element_text(size = 20),
                                                plot.subtitle = element_text(size = 20, face = "italic", hjust = 0.5),
                                                axis.title = element_text(size = 20, face = "bold"),
                                                axis.text = element_text(size = 15, angle = 35),
                                                strip.text = element_text(size = 15))) +
    plot_annotation(title = NULL, 
                    theme =  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))))


ggsave(filename = "04_Output/Figures/Diff_Skewness_RandomCensoring.png", plot = p_Diff_combined, width = 2130*2, height = 1800, units = "px")


################################################################################
#-------------------------------------------------------------------------------
# Diff against sample size at random censoring
#-------------------------------------------------------------------------------
################################################################################

library(ggplot2)
library(patchwork)
library(dplyr)

################################################################################
# Weibull
################################################################################

df_W_logistic <- get(load("02_Data/Weibull_RandomCensoring/df_power_all4_logistic/Power/_Distribution_weibull_resultados_h0_logistic_all_.RData"))
df_W_log <- get(load("02_Data/Weibull_RandomCensoring/df_power_all4_lognormal/Power/_Distribution_weibull_resultados_h0_lognormal_all_.RData"))
df_W_wei <- get(load("02_Data/Weibull_RandomCensoring/df_power_all4_weibull/Power/_Distribution_weibull_resultados_h0_weibull_all_.RData"))

df_W_logistic2 <- df_W_logistic %>%
  mutate(
    Diff_tmax = SEC_min_median - km_min_median,
    n = factor(n), cens = factor(cens, levels = c("0% Censura", "30% Censura", "60% Censura"),
                     labels = c("Complete Data", "30% Censorship", "60% Censorship"))
  )

df_W_log2 <- df_W_log %>%
  mutate(
    Diff_tmax = SEC_min_median - km_min_median,
    n = factor(n), cens = factor(cens, levels = c("0% Censura", "30% Censura", "60% Censura"),
                     labels = c("Complete Data", "30% Censorship", "60% Censorship"))
  )

df_W_wei2 <- df_W_wei %>%
  mutate(
    Diff_tmax = SEC_min_median - km_min_median,
    n = factor(n), cens = factor(cens, levels = c("0% Censura", "30% Censura", "60% Censura"),
                     labels = c("Complete Data", "30% Censorship", "60% Censorship"))
  )

df_plot <- rbind(df_W_logistic2, df_W_log2, df_W_wei2) %>%
  filter(test != "KS No Bootstrap")


(p_Diff_W <- ggplot(df_plot, aes(x = n, y = Diff_tmax, color = cens)) +
    geom_jitter(width = 0.15, size = 3, shape = 16) +
    facet_wrap(~ cens) +
    scale_color_manual(values = c("Complete Data" = "#CAFF70", "30% Censorship" = "#79CDCD", "60% Censorship" = "#FF6A6A")) +
    theme_minimal() +
    labs(x = "Sample Size", y = expression(Diff[t[max]]), color = "Censoring level",
         title = expression(Diff[t[max]] ~ "vs Sample Size under Random Censoring"),
         subtitle = "X ~ Weibull") +
    theme(axis.text.x = element_text(angle = 25, hjust = 1)) + ylim(c(-0.2,0.2)))

p_Diff_W <- p_Diff_W + labs(title = NULL)


################################################################################
# Lognormal
################################################################################

df_L_logistic <- get(load("02_Data/Lognormal_RandomCensoring/df_power_all4_logistic/Power/_Distribution_lognormal_resultados_h0_logistic_all_.RData"))
df_L_log <- get(load("02_Data/Lognormal_RandomCensoring/df_power_all4_lognormal/Power/_Distribution_lognormal_resultados_h0_lognormal_all_.RData"))
df_L_wei <- get(load("02_Data/Lognormal_RandomCensoring/df_power_all4_weibull/Power/_Distribution_lognormal_resultados_h0_weibull_all_.RData"))


df_L_logistic2 <- df_L_logistic %>%
  mutate(
    Diff_tmax = SEC_min_median - km_min_median,
    n = factor(n), cens = factor(cens, levels = c("0% Censura", "30% Censura", "60% Censura"),
                     labels = c("Complete Data", "30% Censorship", "60% Censorship"))
  )

df_L_log2 <- df_L_log %>%
  mutate(
    Diff_tmax = SEC_min_median - km_min_median,
    n = factor(n), cens = factor(cens, levels = c("0% Censura", "30% Censura", "60% Censura"),
                     labels = c("Complete Data", "30% Censorship", "60% Censorship"))
  )

df_L_wei2 <- df_L_wei %>%
  mutate(
    Diff_tmax = SEC_min_median - km_min_median,
    n = factor(n), cens = factor(cens, levels = c("0% Censura", "30% Censura", "60% Censura"),
                     labels = c("Complete Data", "30% Censorship", "60% Censorship"))
  )

df_plot <- rbind(df_L_logistic2, df_L_log2, df_L_wei2) %>%
  filter(test != "KS No Bootstrap")


(p_Diff_L <- ggplot(df_L_logistic2, aes(x = n, y = Diff_tmax, color = cens)) +
    geom_jitter(width = 0.15, size = 3, shape = 16) +
    facet_wrap(~ cens) +
    scale_color_manual(values = c("Complete Data" = "#CAFF70", "30% Censorship" = "#79CDCD", "60% Censorship" = "#FF6A6A")) +
    theme_minimal() +
    labs(x = "Sample Size", y = expression(Diff[t[max]]), color = "Censoring level",
         title = expression(Diff[t[max]] ~ "vs Sample Size under Random Censoring"),
         subtitle = "X ~ Lognormal") +
    theme(axis.text.x = element_text(angle = 25, hjust = 1)) + ylim(c(-0.2,0.2)))

p_Diff_L <- p_Diff_L + labs(title = NULL)


(p_Diff_combined <- (p_Diff_W + p_Diff_L + plot_layout(guides = "collect") &
                       theme_minimal()  & theme(legend.position = "bottom",
                                                legend.title = element_text(size = 20, face = "bold"),
                                                legend.text = element_text(size = 20),
                                                plot.subtitle = element_text(size = 20, face = "italic", hjust = 0.5),
                                                axis.title = element_text(size = 20, face = "bold"),
                                                axis.text = element_text(size = 15, angle = 35),
                                                strip.text = element_text(size = 15))) +
    plot_annotation(title = NULL, 
                    theme =  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))))


ggsave(filename = "04_Output/Figures/Diff_Sample_Size_RandomCensoring.png", plot = p_Diff_combined, width = 2130*2, height = 1800, units = "px")



################################################################################
#-------------------------------------------------------------------------------
# Diff vs power for each test at random censoring
#-------------------------------------------------------------------------------
################################################################################

library(ggplot2)
library(patchwork)
library(dplyr)

################################################################################
# Weibull
################################################################################

df_W_logistic <- get(load("02_Data/Weibull_RandomCensoring/df_power_all4_logistic/Power/_Distribution_weibull_resultados_h0_logistic_all_.RData"))
df_W_log <- get(load("02_Data/Weibull_RandomCensoring/df_power_all4_lognormal/Power/_Distribution_weibull_resultados_h0_lognormal_all_.RData"))
df_W_wei <- get(load("02_Data/Weibull_RandomCensoring/df_power_all4_weibull/Power/_Distribution_weibull_resultados_h0_weibull_all_.RData"))

df_W_logistic2 <- df_W_logistic %>%
  mutate(Distribution = "Logistic",
    Diff_tmax = SEC_min_median - km_min_median,
    n = factor(n), cens = factor(cens, levels = c("0% Censura", "30% Censura", "60% Censura"),
                                 labels = c("Complete Data", "30% Censorship", "60% Censorship"))
  )

df_W_log2 <- df_W_log %>%
  mutate(Distribution = "Lognormal",
    Diff_tmax = SEC_min_median - km_min_median,
    n = factor(n), cens = factor(cens, levels = c("0% Censura", "30% Censura", "60% Censura"),
                                 labels = c("Complete Data", "30% Censorship", "60% Censorship"))
  )

df_W_wei2 <- df_W_wei %>%
  mutate(Distribution = "Weibull",
    Diff_tmax = SEC_min_median - km_min_median,
    n = factor(n), cens = factor(cens, levels = c("0% Censura", "30% Censura", "60% Censura"),
                                 labels = c("Complete Data", "30% Censorship", "60% Censorship"))
  )

################################################################################
#AD

df_plot <- rbind(df_W_logistic2, df_W_log2, df_W_wei2) %>%
  filter(test == "Anderson-Darling")


(p_AD <- ggplot(df_plot, aes(x = Diff_tmax, y = power, shape = Distribution, fill = cens)) +
  geom_point(size = 3.5, color = "black", stroke = 0.7) +
  geom_vline(xintercept = 0, color = "darkblue", linetype = "dashed", size = 1.2) +
  geom_hline(yintercept = 0.1, color = "darkred", linetype = "dotted", size = 1) +
  annotate("text", x = -0.15, y = 0.06, label = "alpha = 0.1", color = "darkred", size = 4.5, hjust = 0) +
  scale_shape_manual(values = c("Logistic" = 22, "Lognormal" = 24, "Weibull" = 23)) +
  scale_fill_manual(values = c("Complete Data" = "#bfff80", "30% Censorship" = "#7cd1d6", "60% Censorship" = "#ff796c")) +
  labs(title = "Power of the Anderson-Darling test by Censoring and Null Hypothesis", subtitle = "Anderson-Darling test", 
       y = "Power", x = expression(paste("Difference:  ", Mean(S[0](t[m]: hat(theta)) - hat(S)(t[m])))),
       fill = "Censoring", shape = "Null Hypothesis") +
  theme_bw() +
  guides(fill = guide_legend(override.aes = list(shape = 21, color = "black")),
         shape = guide_legend(override.aes = list(fill = "black"))) +
  theme(
    plot.title = element_text(size = 16, face = "plain", margin = margin(b = 5)),
    plot.subtitle = element_text(size = 13, face = "plain", margin = margin(b = 15)),
    axis.title.x = element_text(size = 12, margin = margin(t = 10)),
    axis.title.y = element_text(size = 12, margin = margin(r = 10)),
    axis.text = element_text(size = 10, color = "black"),
    panel.grid.major = element_line(color = "#e6e6e6"),
    panel.grid.minor = element_line(color = "#f2f2f2"),
    legend.position = c(0.02, 0.98),
    legend.justification = c("left", "top"),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_blank(),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)) + 
  scale_y_continuous(breaks = c(0.00, 0.25, 0.50, 0.75, 1.00), limits = c(-0.05, 1.05)) +
  scale_x_continuous(breaks = c(-0.2, -0.1, 0.0, 0.1, 0.2), limits = c(-0.22, 0.22)))

p_AD <- p_AD+ labs(title = NULL)


################################################################################
#CvM

df_plot <- rbind(df_W_logistic2, df_W_log2, df_W_wei2) %>%
  filter(test == "Cramér von Mises")


(p_CvM <- ggplot(df_plot, aes(x = Diff_tmax, y = power, shape = Distribution, fill = cens)) +
  geom_point(size = 3.5, color = "black", stroke = 0.7) +
  geom_vline(xintercept = 0, color = "darkblue", linetype = "dashed", size = 1.2) +
  geom_hline(yintercept = 0.1, color = "darkred", linetype = "dotted", size = 1) +
  annotate("text", x = -0.15, y = 0.06, label = "alpha = 0.1", color = "darkred", size = 4.5, hjust = 0) +
  scale_shape_manual(values = c("Logistic" = 22, "Lognormal" = 24, "Weibull" = 23)) +
  scale_fill_manual(values = c("Complete Data" = "#bfff80", "30% Censorship" = "#7cd1d6", "60% Censorship" = "#ff796c")) +
  labs(title = "Power of the Cramér-von Mises test by Censoring and Null Hypothesis", subtitle = "Cramér-von Mises test", 
       y = "Power", x = expression(paste("Difference:  ", Mean(S[0](t[m]: hat(theta)) - hat(S)(t[m])))),
       fill = "Censoring", shape = "Null Hypothesis") +
  theme_bw() +
  guides(fill = guide_legend(override.aes = list(shape = 21, color = "black")),
         shape = guide_legend(override.aes = list(fill = "black"))) +
  theme(
    plot.title = element_text(size = 16, face = "plain", margin = margin(b = 5)),
    plot.subtitle = element_text(size = 13, face = "plain", margin = margin(b = 15)),
    axis.title.x = element_text(size = 12, margin = margin(t = 10)),
    axis.title.y = element_text(size = 12, margin = margin(r = 10)),
    axis.text = element_text(size = 10, color = "black"),
    panel.grid.major = element_line(color = "#e6e6e6"),
    panel.grid.minor = element_line(color = "#f2f2f2"),
    legend.position = c(0.02, 0.98),
    legend.justification = c("left", "top"),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_blank(),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)) + 
  scale_y_continuous(breaks = c(0.00, 0.25, 0.50, 0.75, 1.00), limits = c(-0.05, 1.05)) +
  scale_x_continuous(breaks = c(-0.2, -0.1, 0.0, 0.1, 0.2), limits = c(-0.22, 0.22)))

p_CvM <- p_CvM + labs(title = NULL)

################################################################################
# KS

df_plot <- rbind(df_W_logistic2, df_W_log2, df_W_wei2) %>%
  filter(test == "Kolmogorov-Smirnov")


(p_KS <- ggplot(df_plot, aes(x = Diff_tmax, y = power, shape = Distribution, fill = cens)) +
  geom_point(size = 3.5, color = "black", stroke = 0.7) +
  geom_vline(xintercept = 0, color = "darkblue", linetype = "dashed", size = 1.2) +
  geom_hline(yintercept = 0.1, color = "darkred", linetype = "dotted", size = 1) +
  annotate("text", x = -0.15, y = 0.06, label = "alpha = 0.1", color = "darkred", size = 4.5, hjust = 0) +
  scale_shape_manual(values = c("Logistic" = 22, "Lognormal" = 24, "Weibull" = 23)) +
  scale_fill_manual(values = c("Complete Data" = "#bfff80", "30% Censorship" = "#7cd1d6", "60% Censorship" = "#ff796c")) +
  labs(title = "Power of the Kolmogorov-Smirnov test by Censoring and Null Hypothesis", subtitle = "Kolmogorov-Smirnov test", 
       y = "Power", x = expression(paste("Difference:  ", Mean(S[0](t[m]: hat(theta)) - hat(S)(t[m])))),
       fill = "Censoring", shape = "Null Hypothesis") +
  theme_bw() +
  guides(fill = guide_legend(override.aes = list(shape = 21, color = "black")),
         shape = guide_legend(override.aes = list(fill = "black"))) +
  theme(
    plot.title = element_text(size = 16, face = "plain", margin = margin(b = 5)),
    plot.subtitle = element_text(size = 13, face = "plain", margin = margin(b = 15)),
    axis.title.x = element_text(size = 12, margin = margin(t = 10)),
    axis.title.y = element_text(size = 12, margin = margin(r = 10)),
    axis.text = element_text(size = 10, color = "black"),
    panel.grid.major = element_line(color = "#e6e6e6"),
    panel.grid.minor = element_line(color = "#f2f2f2"),
    legend.position = c(0.02, 0.98),
    legend.justification = c("left", "top"),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_blank(),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)) + 
  scale_y_continuous(breaks = c(0.00, 0.25, 0.50, 0.75, 1.00), limits = c(-0.05, 1.05)) +
  scale_x_continuous(breaks = c(-0.2, -0.1, 0.0, 0.1, 0.2), limits = c(-0.22, 0.22)))

p_KS <- p_KS+ labs(title = NULL)


(p_final <- (p_AD + p_CvM + p_KS + plot_layout(guides = "collect") &
              theme_minimal()  & theme(legend.position = "bottom",
                                       legend.title = element_text(size = 20, face = "bold"),
                                       legend.text = element_text(size = 20),
                                       plot.subtitle = element_text(size = 20, face = "italic", hjust = 0.5),
                                       axis.title = element_text(size = 20, face = "bold"),
                                       axis.text = element_text(size = 15, angle = 35),
                                       strip.text = element_text(size = 15))))

ggsave(filename = "04_Output/Figures/Diff_Power_each_test_Weibull.png", plot = p_final, width = 2130*3, height = 1800, units = "px")





















