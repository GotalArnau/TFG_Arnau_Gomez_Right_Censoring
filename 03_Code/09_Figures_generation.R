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
   annotate("text", x = t_max, y = S_km_tmax + 0.07, label = expression(hat(S)(t[max])), hjust = -0.2, color = "darkgreen", 
            size = 5) +
   annotate("text", x = t_max, y = S0_tmax + 0.07, label = expression(hat(S)[H[0]](t[max], hat(theta)[ML])), hjust = -0.1, 
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
      subtitle = "Survival Function with 30% Administrative Censoring",
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
      subtitle = "Survival Function with 30% Random Censoring",
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

ggsave(filename = "04_Output/Figures/Figure03_Censoring_Visualization_Comparison.png", plot = combined_plots, width = 2130*2, height = 2000, units = "px")


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

df_complet <- rbind(df_W_logistic2, df_W_log2, df_W_wei2)


(p2 <- ggplot(df_complet, aes(x = test, y = power, shape = Hyp_nul, fill = Hyp_nul)) + geom_jitter(size = 2.5, stroke = 0.5, height = 0, width = 0.1) +
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

df_L_logistic <- get(load("C:/Users/arnau.gomez/Desktop/GofCens_Paper_Simulations/02_Data/Lognormal_RandomCensoring/df_power_all4_logistic/SEC/_Distribution_lognormal_resultados_h0_logistic_all_.RData"))
df_L_log <- get(load("C:/Users/arnau.gomez/Desktop/GofCens_Paper_Simulations/02_Data/Lognormal_RandomCensoring/df_power_all4_lognormal/SEC/_Distribution_lognormal_resultados_h0_lognormal_all_.RData"))
df_L_wei <- get(load("C:/Users/arnau.gomez/Desktop/GofCens_Paper_Simulations/02_Data/Lognormal_RandomCensoring/df_power_all4_weibull/SEC/_Distribution_lognormal_resultados_h0_weibull_all_.RData"))


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

df_complet <- rbind(df_L_logistic2, df_L_log2, df_L_wei2)


(p3 <- ggplot(df_complet, aes(x = test, y = power, shape = Hyp_nul, fill = Hyp_nul)) + geom_jitter(size = 2.5, stroke = 0.5, height = 0, width = 0.1) +
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
                                                                              legend.title = element_text(size = 12, face = "bold"),
                                                                              legend.text = element_text(size = 12),
                                                                              plot.subtitle = element_text(size = 12, face = "italic", hjust = 0.5),
                                                                              axis.title = element_text(size = 12, face = "bold"),
                                                                              axis.text = element_text(size = 12))) +
    plot_annotation(title = "Power of each tests for complete data by Null Hypothesis", 
                    theme =  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))))

ggsave(filename = "C:/Users/arnau.gomez/Desktop/GofCens_Paper_Simulations/04_Output/02_Plots/New_plots/Figure3_Merge_Power_complete_data.png", 
       plot = combined_complete_data, width = 2130*2, height = 1800, units = "px")