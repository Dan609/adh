# Immunofluorescence colocalozation analisys
# Bobkov, Poloyanskaya, Goncharova, 2020
# mailto: bobkov@incras.ru, 

library(ggplot2)

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
  
  geom_hline(yintercept=0, linetype="dashed", color = "red") + # add zero
  
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


################################

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
  
  #  scale_x_discrete(labels= CellSciGuylabs) +
  
  labs(y = 'bTau',
       x = "Cell Passage") +
  
  theme(axis.text.x = element_text(angle = 0, hjust = 0)) + 
  
  #  geom_boxplot(alpha = I(0.5), outlier.shape = NA, coef = 100) +
  
  geom_hline(yintercept=0, linetype="dashed", color = "red") +
  
  labs(title = "α-Actinin-4/Hoechst33342",
       subtitle = "colocalization",
       caption = "Lm Coeff = 0.021\nR-squared = 0.29\np < 0.0001") +
  
  geom_smooth(method='lm', formula= y ~ x) +
  
  scale_x_continuous(name = "Cell Passage", 
                     breaks = c(1, seq(8, 17, 1)), 
                     limits = c(8, 17)) +
  
  scale_y_continuous(name = "Kendall’s Tau", 
                     breaks = c(seq(-0.2, 0.6, .1)), 
                     limits = c(-0.2, 0.4))


# add linear regression

fit-aA4_Nuc_d <- lm(bTau ~ probe, data = -aA4_Nuc_d)

summary(fit-aA4_Nuc_d)


################################

# 3 RhoA/F-actin
data <- read.csv('rhoa-Fact-d.csv', sep = ';')

ggplot(data, aes(x = probe, y = bTau)) +
  
  
  ylim(c(-.2, .4)) +
  
  xlim(c(7, 18)) +
  
  theme_classic(base_size=14) +
  
  geom_jitter(position = position_jitter(.1),
              cex = .9,
              shape = 16) +
  
  theme(legend.position = "none") +
  
  #  scale_x_discrete(labels= CellSciGuylabs) +
  
  labs(y = 'bTau',
       x = "Cell Passage") +
  
  theme(axis.text.x = element_text(angle = 0, hjust = 0)) + 
  
  #  geom_boxplot(alpha = I(0.5), outlier.shape = NA, coef = 100) +
  
  geom_hline(yintercept=0, linetype="dashed", color = "red") +
  
  labs(title = "RhoA/F-actin",
       subtitle = "colocalization",
       caption = "Lm Coeff = -0.02421; p =1.41e-12\n Data source: aA4_Nuc_d.csv") +
  
  geom_smooth(method='lm', formula= y ~ x) +
  
  scale_x_continuous(name = "Cell Passage", 
                     breaks = c(1, seq(8, 17, 1)), 
                     limits = c(8, 17)) +
  
  scale_y_continuous(name = "Kendall’s Tau", 
                     breaks = c(seq(-0.2, 0.6, .1)), 
                     limits = c(-0.2, 0.4))


# add linear regression

fit3 <- lm(bTau ~ probe, data = data)

summary(fit3)

################################

# 4 RhoA/F-actin
data <- read.csv('rhoa-Nuc-d.csv', sep = ';')

ggplot(data, aes(x = probe, y = bTau)) +
  
  
  ylim(c(-.2, .4)) +
  
  xlim(c(7, 18)) +
  
  theme_classic(base_size=14) +
  
  geom_jitter(position = position_jitter(.1),
              cex = .9,
              shape = 16) +
  
  theme(legend.position = "none") +
  
  #  scale_x_discrete(labels= CellSciGuylabs) +
  
  labs(y = 'bTau',
       x = "Cell Passage") +
  
  theme(axis.text.x = element_text(angle = 0, hjust = 0)) + 
  
  #  geom_boxplot(alpha = I(0.5), outlier.shape = NA, coef = 100) +
  
  geom_hline(yintercept=0, linetype="dashed", color = "red") +
  
  labs(title = "RhoA/Hoechst33342",
       subtitle = "colocalization",
       caption = "Lm Coeff = 0.030896; p =3.82e-15\n Data source: aA4_Nuc_d.csv") +
  
  geom_smooth(method='lm', formula= y ~ x)


# add linear regression

fit4 <- lm(bTau ~ probe, data = data)

summary(fit4)


