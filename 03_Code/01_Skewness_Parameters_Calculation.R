################################################################################
#-------------------------------------------------------------------------------
# We must determine the parameters of the generative distributions so that they 
# are consistent with our desired level of skewness. We start by estimating them
# visually and then refine the selection using Bowley’s and Kelly’s coefficients
# of skewness. Here we present the parameters used in the final study.
#-------------------------------------------------------------------------------
################################################################################

library(GofCens)
library(dplyr)
library(ggplot2)
library(flexsurv)
library(latex2exp)


#-----------------------------------------------------------------------------------------------------------------------------
# Weibull Simetrico
#-----------------------------------------------------------------------------------------------------------------------------

# Bowley's Coefficient:

Q1 <- qweibull(0.25, shape = 3.5, scale=20)
Q2 <- qweibull(0.5, shape = 3.5, scale=20)
Q3 <- qweibull(0.75, shape = 3.5, scale=20)

(B1 <- (Q3 + Q1 - 2*Q2) / (Q3 - Q1))


# Kelly's Coefficient:
P1 <- qweibull(0.1, shape = 3.5, scale=20)
Q2 <- qweibull(0.5, shape = 3.5, scale=20)
P3 <- qweibull(0.9, shape = 3.5, scale=20)

(K1 <- (P3 + P1 - 2*Q2) / (P3 - P1))


#-----------------------------------------------------------------------------------------------------------------------------
# Weibull Asimetrico (Actual)
#-----------------------------------------------------------------------------------------------------------------------------

# Bowley's Coefficient:

Q1 <- qweibull(0.25, shape = 2, scale=20)
Q2 <- qweibull(0.5, shape = 2, scale=20)
Q3 <- qweibull(0.75, shape = 2, scale=20)

(B2 <- (Q3 + Q1 - 2*Q2) / (Q3 - Q1))


# Kelly's Coefficient:

P1 <- qweibull(0.1, shape = 2, scale=20)
Q2 <- qweibull(0.5, shape = 2, scale=20)
P3 <- qweibull(0.9, shape = 2, scale=20)

(K2 <- (P3 + P1 - 2*Q2) / (P3 - P1))


#-----------------------------------------------------------------------------------------------------------------------------
# Weibull Asimetrico (Nueva propuesta)
#-----------------------------------------------------------------------------------------------------------------------------

# Bowley's Coefficient:

Q1 <- qweibull(0.25, shape = 1.7, scale=20)
Q2 <- qweibull(0.5, shape = 1.7, scale=20)
Q3 <- qweibull(0.75, shape = 1.7, scale=20)

(B22 <- (Q3 + Q1 - 2*Q2) / (Q3 - Q1))


# Kelly's Coefficient:

P1 <- qweibull(0.1, shape = 1.7, scale=20)
Q2 <- qweibull(0.5, shape = 1.7, scale=20)
P3 <- qweibull(0.9, shape = 1.7, scale=20)

(K22 <- (P3 + P1 - 2*Q2) / (P3 - P1))


#-----------------------------------------------------------------------------------------------------------------------------
# Weibull Muy Asimetrico
#-----------------------------------------------------------------------------------------------------------------------------

# Bowley's Coefficient:

Q1 <- qweibull(0.25, shape = 1, scale=20)
Q2 <- qweibull(0.5, shape = 1, scale=20)
Q3 <- qweibull(0.75, shape = 1, scale=20)

(B3 <- (Q3 + Q1 - 2*Q2) / (Q3 - Q1))

# Kelly's Coefficient:

P1 <- qweibull(0.1, shape = 1, scale=20)
Q2 <- qweibull(0.5, shape = 1, scale=20)
P3 <- qweibull(0.9, shape = 1, scale=20)

(K3 <- (P3 + P1 - 2*Q2) / (P3 - P1))


Skewness_Weibull <- data.frame(
  "Skewness" = c("Symmetric","Asymmetric Actual", "Asymmentric New","Extremly Asymmetric"),
  "Shape" = c(3.5,2,1.7,1),
  "Bowly's" = c(B1,B2,B22,B3),
  "Kelly's" = c(K1,K2,K22,K3)
)

set.seed(28657)

data_w <- data.frame("Rand_w" = c(rweibull(n = 10000, shape = 3.5, scale=20),
                                  rweibull(n = 10000, shape = 2, scale=20),
                                  rweibull(n = 10000, shape = 1.7, scale=20),
                                  rweibull(n = 10000, shape = 1, scale=20)),
                     "Skewness" = rep(c("Symmetric","Asymmetric Actual", "Asymmetric New", "Extremely Asymmetric"), each = 10000))


(plot_weibull <- ggplot(data_w, aes(x = Rand_w, colour = Skewness)) + geom_density(linewidth = 1) +
  scale_colour_manual(values = c("Symmetric" = "darkorchid", "Asymmetric New" = "orange",
                                 "Asymmetric Actual" = "darkorange4", "Extremely Asymmetric" = "cyan2"),
                      breaks = c("Symmetric", "Asymmetric Actual", "Asymmetric New", "Extremely Asymmetric"),
                      labels = c(expression("Symmetric "(alpha == 3.5)),
                                 expression("Asymmetric Actual "(alpha == 2)),
                                 expression("Asymmetric New "(alpha == 1.7)),
                                 expression("Extremely asymmetric "(alpha == 1)))) +
  labs(title = expression("Data ~ Weibull "( alpha, beta == 20)), x = NULL, y = "Density") +
  theme_minimal() + lims(x = c(0,110)) + 
  theme(legend.position = c(0.7,0.7)))
  
  






#-----------------------------------------------------------------------------------------------------------------------------
# Lognormal Simetrico
#-----------------------------------------------------------------------------------------------------------------------------

# Bowley's Coefficient:

Q1 <- qlnorm(0.25, sdlog = 0.01, meanlog=20)
Q2 <- qlnorm(0.5, sdlog = 0.01, meanlog=20)
Q3 <- qlnorm(0.75, sdlog = 0.01, meanlog=20)

(B1 <- (Q3 + Q1 - 2*Q2) / (Q3 - Q1))


# Kelly's Coefficient:
P1 <- qlnorm(0.1, sdlog = 0.01, meanlog=20)
Q2 <- qlnorm(0.5, sdlog = 0.01, meanlog=20)
P3 <- qlnorm(0.9, sdlog = 0.01, meanlog=20)

(K1 <- (P3 + P1 - 2*Q2) / (P3 - P1))


#-----------------------------------------------------------------------------------------------------------------------------
# Lognormal Asimetrico (Actual)
#-----------------------------------------------------------------------------------------------------------------------------

# Bowley's Coefficient:

Q1 <- qlnorm(0.25, sdlog = 0.4, meanlog=20)
Q2 <- qlnorm(0.5, sdlog = 0.4, meanlog=20)
Q3 <- qlnorm(0.75, sdlog = 0.4, meanlog=20)

(B2 <- (Q3 + Q1 - 2*Q2) / (Q3 - Q1))


# Kelly's Coefficient:

P1 <- qlnorm(0.1, sdlog = 0.4, meanlog=20)
Q2 <- qlnorm(0.5, sdlog = 0.4, meanlog=20)
P3 <- qlnorm(0.9, sdlog = 0.4, meanlog=20)

(K2 <- (P3 + P1 - 2*Q2) / (P3 - P1))


#-----------------------------------------------------------------------------------------------------------------------------
# Lognormal Muy Asimetrico
#-----------------------------------------------------------------------------------------------------------------------------

# Bowley's Coefficient:

Q1 <- qlnorm(0.25, sdlog = 1.5, meanlog=20)
Q2 <- qlnorm(0.5, sdlog = 1.5, meanlog=20)
Q3 <- qlnorm(0.75, sdlog = 1.5, meanlog=20)

(B3 <- (Q3 + Q1 - 2*Q2) / (Q3 - Q1))

# Kelly's Coefficient:

P1 <- qlnorm(0.1, sdlog = 1.5, meanlog=20)
Q2 <- qlnorm(0.5, sdlog = 1.5, meanlog=20)
P3 <- qlnorm(0.9, sdlog = 1.5, meanlog=20)

(K3 <- (P3 + P1 - 2*Q2) / (P3 - P1))


Skewness_Lognormal <- data.frame(
  "Skewness" = c("Symmetric","Asymmetric", "Extremly Asymmetric"),
  "Sdlog" = c(0.01,0.4,1.5),
  "Bowly's" = c(B1,B2,B3),
  "Kelly's" = c(K1,K2,K3)
)


set.seed(28657)

data_l <- data.frame("Rand_l" = c(rlnorm(n = 10000, sdlog = 0.01, meanlog=7),
                                  rlnorm(n = 10000, sdlog = 0.4, meanlog=7),
                                  rlnorm(n = 10000, sdlog = 1.5, meanlog=7)),
                      "Skewness" = rep(c("Symmetric","Asymmetric", "Extremely Asymmetric"), each = 10000))


(plot_lognormal <- ggplot(data_l, aes(x = Rand_l, colour = Skewness)) + geom_density(linewidth = 1) +
  scale_colour_manual(values = c("Symmetric" = "darkorchid", "Asymmetric" = "orange", "Extremely Asymmetric" = "cyan2"),
                      breaks = c("Symmetric", "Asymmetric", "Extremely Asymmetric"),
                      labels = c(expression("Symmetric "(sigma[log] == 0.01)),
                                 expression("Asymmetric "(sigma[log] == 0.4)),
                                 expression("Extremely asymmetric "(sigma[log] == 1.5)))) +
  
  labs(title = expression("Data ~ Lognormal "( sigma[log], mu == 7)), x = NULL, y = "Density") +
  theme_minimal() + lims(x = c(0,5000)) + 
  theme(legend.position = c(0.7,0.7)))

(plot_lognormal_zoom <- plot_lognormal + coord_cartesian(ylim = c(0,0.003)))










#-----------------------------------------------------------------------------------------------------------------------------
# Log - Logistic Simetrico
#-----------------------------------------------------------------------------------------------------------------------------

# Bowley's Coefficient:

Q1 <- qllogis(0.25, shape = 50, scale = 20)
Q2 <- qllogis(0.5, shape = 50, scale = 20)
Q3 <- qllogis(0.75, shape = 50, scale = 20)

(B1 <- (Q3 + Q1 - 2*Q2) / (Q3 - Q1))

# Kelly's Coefficient:

P1 <- qllogis(0.1, shape = 50, scale = 20)
Q2 <- qllogis(0.5, shape = 50, scale = 20)
P3 <- qllogis(0.9, shape = 50, scale = 20)

(K1 <- (P3 + P1 - 2*Q2) / (P3 - P1))


#-----------------------------------------------------------------------------------------------------------------------------
# Log - Logistic Asimetrico
#-----------------------------------------------------------------------------------------------------------------------------

# Bowley's Coefficient:

Q1 <- qllogis(0.25, shape = 4, scale = 20)
Q2 <- qllogis(0.5, shape = 4, scale = 20)
Q3 <- qllogis(0.75, shape = 4, scale = 20)

(B2 <- (Q3 + Q1 - 2*Q2) / (Q3 - Q1))

# Kelly's Coefficient:

P1 <- qllogis(0.1, shape = 4, scale = 20)
Q2 <- qllogis(0.5, shape = 4, scale = 20)
P3 <- qllogis(0.9, shape = 4, scale = 20)

(K2 <- (P3 + P1 - 2*Q2) / (P3 - P1))


#-----------------------------------------------------------------------------------------------------------------------------
# Log - Logistic Muy Aimetrico
#-----------------------------------------------------------------------------------------------------------------------------

# Bowley's Coefficient:

Q1 <- qllogis(0.25, shape = 1.5, scale = 2)
Q2 <- qllogis(0.5, shape = 1.5, scale = 20)
Q3 <- qllogis(0.75, shape = 1.5, scale = 20)

(B3 <- (Q3 + Q1 - 2*Q2) / (Q3 - Q1))

# Kelly's Coefficient:

P1 <- qllogis(0.1, shape = 1.5, scale = 20)
Q2 <- qllogis(0.5, shape = 1.5, scale = 20)
P3 <- qllogis(0.9, shape = 1.5, scale = 20)

(K3 <- (P3 + P1 - 2*Q2) / (P3 - P1))


Skewness_Loglogistic <- data.frame(
  "Skewness" = c("Symmetric","Asymmetric", "Extremly Asymmetric"),
  "Sdlog" = c(50,4,1.5),
  "Bowly's" = c(B1,B2,B3),
  "Kelly's" = c(K1,K2,K3)
)

set.seed(28657)

data_ll <- data.frame("Rand_ll" = c(rllogis(n = 10000, shape = 10, scale = 20),
                                      rllogis(n = 10000, shape = 3, scale = 20),
                                      rllogis(n = 10000, shape = 1.5, scale = 20)),
                      "Skewness" = rep(c("Symmetric","Asymmetric", "Extremely Asymmetric"), each = 10000))


(plot_loglogistic <- ggplot(data_ll, aes(x = Rand_ll, colour = Skewness)) + geom_density(linewidth = 1) +
    scale_colour_manual(values = c("Symmetric" = "darkorchid", "Asymmetric" = "orange", "Extremely Asymmetric" = "cyan2"),
                        breaks = c("Symmetric", "Asymmetric", "Extremely Asymmetric"),
                        labels = c(expression("Symmetric "(alpha == 10)),
                                   expression("Asymmetric "(alpha == 3)),
                                   expression("Extremely asymmetric "(alpha == 1.5)))) +
    labs(title = expression("Data ~ Log-Logistic "( alpha, beta == 20)), x = NULL, y = "Density") +
    theme_minimal() + lims(x = c(0,110)) + 
    theme(legend.position = c(0.7,0.7)))




Skewness_Weibull
Skewness_Lognormal
Skewness_Loglogistic
