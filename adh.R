# Immunofluorescence colocalozation analisys
# Bobkov, Poloyanskaya, Goncharova, 2020
# mailto: bobkov@incras.ru
# https://github.com/dan609/adh

library(ggplot2)
######################################
#1# α-Actinin-4/F-actin colocalization 
data_aA4_Fact_d <- read.csv('aA4_Fact_d.csv', sep = ';')

ggplot(data_aA4_Fact_d, aes(x = probe, y = bTau)) +
  
  
  ylim(c(-.2, .6)) +
  
  xlim(c(7, 18)) +
  
  theme_classic(base_size=14) +
  
  geom_jitter(position = position_jitter(.1),
              cex = .9,
              shape = 16) +
  
  theme(legend.position = "none") +
  
  theme(axis.text.x = element_text(angle = 0, hjust = 0)) + 
  
  geom_hline(yintercept=0, linetype="dashed", color = "red") + # add zero level
  
  labs(title = "α-Actinin-4/F-actin",
       subtitle = "colocalization",
       caption = "Lm Coeff = -0.015\n R-squared = 0.11\np < 0.0001") +
  
  geom_smooth(method='lm', formula= y ~ x) + # add linear regression
  
  scale_x_continuous(name = "Cell Passage", 
                     breaks = c(1, seq(8, 17, 1)), 
                     limits = c(8, 17)) +
  
  scale_y_continuous(name = "Kendall’s Tau", 
                     breaks = c(seq(-0.1, 0.6, .1)), 
                     limits = c(-0.1, 0.6))


# Compute linear regression

fit_aA4_Fact_d <- lm(bTau ~ probe, data = data_aA4_Fact_d)

summary(fit_aA4_Fact_d)


###########################################
#2# α-Actinin-4/Hoechst33342 colocalozation
data_aA4_Nuc_d <- read.csv('aA4_Nuc_d.csv', sep = ';')

ggplot(data_aA4_Nuc_d, aes(x = probe, y = bTau)) +
  
  
  ylim(c(-.2, .4)) +
  
  xlim(c(7, 18)) +
  
  theme_classic(base_size=14) +
  
  geom_jitter(position = position_jitter(.1),
              cex = .9,
              shape = 16) +
  
  theme(legend.position = "none") +
  
  theme(axis.text.x = element_text(angle = 0, hjust = 0)) + 
  
  geom_hline(yintercept=0, linetype="dashed", color = "red") + # add zero level
  
  labs(title = "α-Actinin-4/Hoechst33342",
       subtitle = "colocalization",
       caption = "Lm Coeff = 0.021\nR-squared = 0.29\np < 0.0001") +
  
  geom_smooth(method='lm', formula= y ~ x) + # add linear regression
  
  scale_x_continuous(name = "Cell Passage", 
                     breaks = c(1, seq(8, 17, 1)), 
                     limits = c(8, 17)) +
  
  scale_y_continuous(name = "Kendall’s Tau", 
                     breaks = c(seq(-0.2, 0.6, .1)), 
                     limits = c(-0.2, 0.4))


# Compute linear regression

fit_aA4_Nuc_d <- lm(bTau ~ probe, data = data_aA4_Nuc_d)

summary(fit_aA4_Nuc_d)


###############################
#3# RhoA/F-actin colocalization 
data_rhoa_Fact_d <- read.csv('rhoa_Fact_d.csv', sep = ';')

ggplot(data_rhoa_Fact_d, aes(x = probe, y = bTau)) +
  
  
  ylim(c(-.2, .4)) +
  
  xlim(c(7, 18)) +
  
  theme_classic(base_size=14) +
  
  geom_jitter(position = position_jitter(.1),
              cex = .9,
              shape = 16) +
  
  theme(legend.position = "none") +
  
  theme(axis.text.x = element_text(angle = 0, hjust = 0)) + 
  
  geom_hline(yintercept=0, linetype="dashed", color = "red") + # add zero level
  
  labs(title = "RhoA/F-actin",
       subtitle = "colocalization",
       caption = "Lm Coeff = -0.024\nR-squared = 0.22\np < 0.0001") +
  
  geom_smooth(method='lm', formula= y ~ x) + # add linear regression
  
  scale_x_continuous(name = "Cell Passage", 
                     breaks = c(1, seq(8, 17, 1)), 
                     limits = c(8, 17)) +
  
  scale_y_continuous(name = "Kendall’s Tau", 
                     breaks = c(seq(-0.2, 0.6, .1)), 
                     limits = c(-0.2, 0.4))


# Compute linear regression

fit_rhoa_Fact_d <- lm(bTau ~ probe, data = data_rhoa_Fact_d)

summary(fit_rhoa_Fact_d)

###############################
#4# RhoA/F-actin colocalization 
data_rhoa_Nuc_d <- read.csv('rhoa_Nuc_d.csv', sep = ';')

ggplot(data_rhoa_Nuc_d, aes(x = probe, y = bTau)) +
  
  
  ylim(c(-.2, .4)) +
  
  xlim(c(7, 18)) +
  
  theme_classic(base_size=14) +
  
  geom_jitter(position = position_jitter(.1),
              cex = .9,
              shape = 16) +
  
  theme(legend.position = "none") +
  
  theme(axis.text.x = element_text(angle = 0, hjust = 0)) + 
  
  geom_hline(yintercept=0, linetype="dashed", color = "red") + # add zero level
  
  labs(title = "RhoA/Hoechst33342",
       subtitle = "colocalization",
       caption = "Lm Coeff = 0.021\nR-squared = 0.26\np < 0.0001") +
  
  geom_smooth(method='lm', formula= y ~ x) + # add linear regression
  
  scale_x_continuous(name = "Cell Passage", 
                     breaks = c(1, seq(8, 17, 1)), 
                     limits = c(8, 17)) +
  
  scale_y_continuous(name = "Kendall’s Tau", 
                     breaks = c(seq(-0.2, 0.6, .1)), 
                     limits = c(-0.2, 0.4))



# Compute linear regression

fit_rhoa_Nuc_d <- lm(bTau ~ probe, data = data_rhoa_Nuc_d)

summary(fit_rhoa_Nuc_d)


