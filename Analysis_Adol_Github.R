
#CC Joe Barnby 2025
# Libraries ---------------------------------------------------------------
library(R.matlab)  # Functions for handling mat files etc.
library(Hmisc)
library(ppcor)
library(tidyverse)
library(foreach)
library(doParallel)
library(tidyquant)
source('Utilities.R')
library(patchwork)
library(ggpubr)
library(BayesianFirstAid)

theme_set(theme_bw(base_size=18))

# Read Data ---------------------------------------------------------------

check_manipulation_tidy <- read.csv('Data/check_tidy.csv')
check_manipulation_full <- read.csv('Data/check_full.csv')
adolIntent_clean        <- read.csv('Data/AdolIntent_Clean.csv') %>% 
                            dplyr::select(-X) %>%
                            filter(trial != is.na(trial))

IDinc <- check_manipulation_full %>%
  filter(Context == 'Inc') %>%
  dplyr::select(ID) %>%
  unique()

IDexc <- check_manipulation_full %>%
  filter(Context == 'Exc') %>%
  dplyr::select(ID) %>%
  unique()

# Colours for groups ------------------------------------------------------

colour_group_adol <- c('#42E2B8', '#07004D')

# Behavioural analysis ----------------------------------------------------

beh_anal_adol <- check_manipulation_tidy %>%
  mutate(distancebeta=abs(server_beta_par-server_beta_ppt),
         distancealpha=abs(server_alpha_par-server_alpha_ppt),
         across(Q6:DASSTot, ~ifelse(.x==9999, NA, .x))
         ) %>%
  dplyr::select(ID, SI:last_col()) %>%
  distinct()

## Psychometric curves and correctness --------

curve_ex <- glm(correct ~ trial,
    data = check_manipulation_full %>% filter(Context=='Exc'),
    family = "binomial")
curve_in <- glm(correct ~ trial,
    data = check_manipulation_full %>% filter(Context=='Inc'),
    family = "binomial")

glm(correct ~ trial,
    data = check_manipulation_full,
    family = "binomial") %>%
  summary()

pred_time_dat <- data.frame(trial=seq(0, 54, 0.01))

prob_in <- curve_in %>% predict(pred_time_dat, type = "response")
prob_ex <- curve_ex %>% predict(pred_time_dat, type = "response")

pred_time_dat$prob_ex = prob_ex
pred_time_dat$prob_in = prob_in

pred_time_dat <- pred_time_dat %>%
  pivot_longer(2:3, names_to = 'Group', values_to = 'p(Correct)') %>%
  mutate(Group = ifelse(Group=='prob_in', 'in', 'ex'))

log_corp1 <- ggplot(pred_time_dat,
       aes(trial, `p(Correct)`, colour = Group))+
  geom_smooth(se = T, size = 1) +  # Smoother lines without confidence intervals
  scale_color_manual(values=colour_group_adol)+
  labs(x = 'Trial', y = 'P(Correct)') +  # Proper axis labels
  theme(
    legend.position = c(0.75, 0.25),
    legend.background = element_rect(colour = 'black', fill = 'white'),
    legend.title = element_blank(),  # Improve legend title readability
    legend.text = element_text(size = 12)  # Improve legend text readability
  )+
ggplot(check_manipulation_full %>%
         filter(Phase==2) %>%
         dplyr::select(Context, ID, correctSum) %>%
         distinct(),
       aes(as.factor(Context), correctSum/54, fill = as.factor(Context)))+
  geom_jitter(alpha = 0.7, shape = 21, width = 0.1, size = 2) +
  geom_boxplot(width = 0.2, alpha = 0.7, outlier.shape = NA) +
  geom_hline(yintercept = 0.5, size = 1.2)+
  scale_fill_manual(values=colour_group_adol)+
  scale_y_continuous(expand=c(0,0))+
  stat_compare_means(label.y = 0.25, size = 5)+
  #coord_cartesian(y = c(0, 1))+
  labs(x = 'Group', y = '% Correct')+
  theme(
      legend.position = 'none',
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10)),
      plot.title = element_text(hjust = 0.5)
    )

log_corp1

check_manipulation_full %>%
         filter(Phase==2) %>%
         dplyr::select(Context, ID, correctSum) %>%
         distinct() %>%
  ungroup() %>%
  summarise(cor = mean(correctSum)/54, corsd = sd(correctSum))

## Behavioural choice ------------------------------------------------------

### Phase 1 ------
choices_ppts_1 <- adolIntent_clean %>%
  filter(Phase==1, ID %in% beh_anal_adol$ID) %>%
  mutate(group=ifelse(ID%in%IDinc$ID, 'in', 'ex'),
         relation = case_when(
              (S1 == O1 & S2 > O2 & S2 > S1 & choice==1) | (S2 == O2 & S1 > O1 & S1 > S2 & choice == 2) ~ "Pi",
              (S1 == O1 & S2 > O2 & O1 > O2 & choice==1) | (S2 == O2 & S1 > O1 & O2 > O1 & choice == 2) ~ "Pc",
              (S1 >  O1 & S2 > O2 & S1 > S2 & choice==1) | (S1 >  O1 & S2 > O2 & S2 > S1 & choice == 2) ~ 'Ic',
              TRUE ~ NA
            )) %>%
  dplyr::select(relation, group, trial, ID) %>%
  distinct()

summary_choices_ppts_1 <- choices_ppts_1 %>%
  group_by(ID, relation, group) %>%
  summarise(total=n()) %>%
  pivot_wider(id_cols = c(ID, group), names_from = 'relation', values_from = 'total') %>%
  mutate(Pc=ifelse(is.na(Pc), 0, Pc),
         Ic=ifelse(is.na(Ic), 0, Ic),
         Pi=ifelse(is.na(Pi), 0, Pi)
         )

mean(summary_choices_ppts_1$Ic); sd(summary_choices_ppts_1$Ic)
mean(summary_choices_ppts_1$Pc); sd(summary_choices_ppts_1$Pc)
mean(summary_choices_ppts_1$Pi); sd(summary_choices_ppts_1$Pi)

summary_choices_ppts_1 %>%
  group_by(group) %>%
  mutate(
         meanPc = mean(Pc),
         meanIc = mean(Ic),
         meanPi = mean(Pi),
         sdPc = sd(Pc),
         sdIc = sd(Ic),
         sdPi = sd(Pi),
         ) %>%
  dplyr::select(meanPc:sdPi, group) %>%
  distinct() %>%
  as.data.frame()

t.test(Pc ~ group, data=summary_choices_ppts_1)
t.test(Pi ~ group, data=summary_choices_ppts_1)
t.test(Ic ~ group, data=summary_choices_ppts_1)

dplyr::full_join(check_manipulation_tidy %>% dplyr::select(ID, RGPTSPers) %>% distinct(),
                 summary_choices_ppts_1,
                 by = 'ID') %>%
  lm(Pi ~ RGPTSPers, data=.) %>%
  summary()

choices_ppts_3 <- adolIntent_clean %>%
  filter(Phase==3, ID %in% beh_anal_adol$ID) %>%
  mutate(group=ifelse(ID%in%IDinc$ID, 'in', 'ex'),
         relation = case_when(
              (S1 == O1 & S2 > O2 & S2 > S1 & choice==1) | (S2 == O2 & S1 > O1 & S1 > S2 & choice == 2) ~ "Pi",
              (S1 == O1 & S2 > O2 & O1 > O2 & choice==1) | (S2 == O2 & S1 > O1 & O2 > O1 & choice == 2) ~ "Pc",
              (S1 >  O1 & S2 > O2 & S1 > S2 & choice==1) | (S1 >  O1 & S2 > O2 & S2 > S1 & choice == 2) ~ 'Ic',
              TRUE ~ NA
            )) %>%
  dplyr::select(relation, group, trial, ID) %>%
  distinct()

summary_choices_ppts_3 <- choices_ppts_3 %>%
  group_by(ID, relation, group) %>%
  summarise(total=n()) %>%
  pivot_wider(id_cols = c(ID, group), names_from = 'relation', values_from = 'total') %>%
  mutate(
         Pc=ifelse(is.na(Pc), 0, Pc),
         Ic=ifelse(is.na(Ic), 0, Ic),
         Pi=ifelse(is.na(Pi), 0, Pi)
         )

mean(summary_choices_ppts_3$Ic); sd(summary_choices_ppts_3$Ic)
mean(summary_choices_ppts_3$Pc); sd(summary_choices_ppts_3$Pc)
mean(summary_choices_ppts_3$Pi); sd(summary_choices_ppts_3$Pi)

summary_choices_ppts_3 %>%
  group_by(group) %>%
  mutate(
         meanPc = mean(Pc),
         meanIc = mean(Ic),
         meanPi = mean(Pi),
         sdPc = sd(Pc),
         sdIc = sd(Ic),
         sdPi = sd(Pi),
         ) %>%
  dplyr::select(meanPc:sdPi, group) %>%
  distinct() %>%
  as.data.frame()

t.test(Pc ~ group, data=summary_choices_ppts_3)
t.test(Pi ~ group, data=summary_choices_ppts_3)
t.test(Ic ~ group, data=summary_choices_ppts_3)

choices_ppts_diff <- adolIntent_clean %>%
  filter(Phase%in%c(1,3)) %>%
  mutate(group=ifelse(ID%in%IDinc$ID, 'in', 'ex'),
         relation = case_when(
              (S1 == O1 & S2 > O2 & S2 > S1 & choice==1) | (S2 == O2 & S1 > O1 & S1 > S2 & choice == 2) ~ "Pi",
              (S1 == O1 & S2 > O2 & O1 > O2 & choice==1) | (S2 == O2 & S1 > O1 & O2 > O1 & choice == 2) ~ "Pc",
              (S1 >  O1 & S2 > O2 & S1 > S2 & choice==1) | (S1 >  O1 & S2 > O2 & S2 > S1 & choice == 2) ~ 'Ic',
              TRUE ~ NA
            )) %>%
  dplyr::select(relation, group, trial, Phase, ID) %>%
  distinct()

summary_choices_ppts_diff <- choices_ppts_diff %>%
  filter(ID %in% beh_anal_adol$ID) %>%
  group_by(ID, relation, group, Phase) %>%
  summarise(total=n()) %>%
  pivot_wider(id_cols = c(ID, group, Phase), names_from = 'relation', values_from = 'total') %>%
  mutate(
         Pc=ifelse(is.na(Pc), 0, Pc),
         Ic=ifelse(is.na(Ic), 0, Ic),
         Pi=ifelse(is.na(Pi), 0, Pi)
         )

summary_choices_ppts_diff %>%
  group_by(group, Phase) %>%
  mutate(
         meanPc = mean(Pc),
         meanIc = mean(Ic),
         meanPi = mean(Pi),
         sdPc = sd(Pc),
         sdIc = sd(Ic),
         sdPi = sd(Pi),
         ) %>%
  dplyr::select(meanPc:sdPi, group) %>%
  distinct() %>%
  as.data.frame()

t.test(Pc ~ Phase, data=summary_choices_ppts_diff)
t.test(Pi ~ Phase, data=summary_choices_ppts_diff)
t.test(Ic ~ Phase, data=summary_choices_ppts_diff)

# Model Simulations -------------------------------------------------------

## Simplex ----------

simpLong <- simulate_simplex(adolIntent_clean, res = 0.5, v = 0.1)

#place participants on grid
simpPPTs <- ggplot(simpLong$diff, aes(alpha, beta, fill = val)) +
  geom_tile() +
  scale_fill_gradient(low = 'black', high = 'white', name = expression(paste(Delta, 'Value'))) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = expression(paste(alpha)), y = expression(paste(beta))) +
  theme_bw() +
  theme(text = element_text(size = 20),
        panel.grid = element_blank(),
        panel.border = element_blank())
simpPPTs

## Phase 3 simulations ------
s_uncertainty <- seq(1, 15, 1)
o_uncertainty <- seq(1, 15, 1)

loop_values <- expand.grid(s_uncertainty, o_uncertainty)
shift_beta  <- matrix(NA,
                      nrow = nrow(loop_values),
                      ncol = (ncol(loop_values) + ncol(loop_values) + 2)
                      )

IDn <- sample(IDinc$ID,1)

for(i in 1:nrow(loop_values)){

  datuse <- adolIntent_clean[adolIntent_clean$ID==IDn,]

  par_alpha <- datuse$server_alpha_par[1]
  par_beta  <- datuse$server_beta_par[2]

  print(loop_values[i,])
  parms <- c(par_alpha, -20, 1, loop_values[i,1], 5, 10, 1, loop_values[i,2])
  run1  <- ABA_shift_Gen(parms,
                        datuse,
                        sim=1,
                        plot = 0,
                        proj = 0)
  ppt_act   <- simulate_phase_decisions(parms[1:2], datuse)
  par_act   <- simulate_phase_decisions(c(par_alpha, par_beta), datuse)

  match     <- ifelse(ppt_act==par_act, 1, 0); matched = sum(match[,5])

  shift_beta[i,1] <- run1$Shift[2]
  shift_beta[i,2] <- run1$Shift[4]
  shift_beta[i,3] <- matched
  shift_beta[i,4] <- loop_values[i,1]
  shift_beta[i,5] <- loop_values[i,2]

}

colnames(shift_beta)      <- c('m_shift', 'sd_shift', 'match', 's_u', 'o_u', 'disp')

shift_beta_l <- shift_beta %>%
  as.data.frame() %>%
  mutate(s_u = factor(s_u),
         o_u = factor(o_u)) %>%
  pivot_longer(1:2, names_to = 'parameter', values_to = 'shift')

sd_sim <- ggplot(shift_beta_l %>%
         filter(parameter == 'sd_shift'),
       aes(s_u, o_u, fill = shift))+
  scale_fill_gradient(low = 'black',high = 'white',
                      name = expression(paste(Delta, theta[ppt]^sigma)))+
  geom_tile()+
  labs(y = expression(paste('Other Unc.')),
       x = expression(paste('Self Unc.'))
       ) +
  scale_y_discrete(breaks = c(1, 5, 10, 15)) +
  scale_x_discrete(breaks = c(1, 5, 10, 15)) +
  theme_bw(base_size=18)

m_sim <- ggplot(shift_beta_l %>%
         filter(parameter == 'm_shift'),
       aes(s_u, o_u, fill = shift))+
  scale_fill_gradient(low = 'white', high = 'black',
                      name = expression(paste(Delta, theta[ppt]^mu)))+
  geom_tile() +
  labs(y = expression(paste('Other Unc.')),
       x = expression(paste('Self Unc.'))
       ) +
  scale_y_discrete(breaks = c(1, 5, 10, 15)) +
  scale_x_discrete(breaks = c(1, 5, 10, 15)) +
  theme_bw(base_size=18)

sd_sim/m_sim

# Group Matlab Analysis ---------------------------------------------------

# ORDER FOR ANALYSIS:

# - Fit each individual to every model (laplace)
# - Fit both groups hierarchically (HBI)
# - Simulate data from each model
# - Refit the simulated data from each individual to every model (laplace)
# - Fit both groups hierarchically (HBI)

#Fitted hierarchical Parameters
ADOL_ex      <- readMat('HBI/hbi_adolescence_excluded.mat')
ADOL_in      <- readMat('HBI/hbi_adolescence_included.mat')

#Simulated choices with best fitting model
ADOL_ex_sim       <- readMat('data_sim/data_sim_excluded.mat')
ADOL_in_sim       <- readMat('data_sim/data_sim_included.mat')

ADOL_ex_full_sim  <- readMat('data_sim/full_sim_excluded.mat')
ADOL_in_full_sim  <- readMat('data_sim/full_sim_included.mat')

#Recovered files for group
ADOL_ex_recov <- readMat('HBI/hbi_adolescence_exc_recovery.mat')
ADOL_in_recov <- readMat('HBI/hbi_adolescence_inc_recovery.mat')

## Individual parameters --------
ID_ex <- beh_anal_adol$ID[beh_anal_adol$Context=='Exc']
ID_in <- beh_anal_adol$ID[beh_anal_adol$Context=='Inc']

ADOL_parms_ex <- extract_parameters_a(ADOL_ex, unique(ID_ex))
ADOL_parms_in <- extract_parameters_a(ADOL_in, unique(ID_in))

joint_parms_adol <- rbind(ADOL_parms_ex, ADOL_parms_in) %>%
               plyr::join(., adolIntent_clean %>%
                            filter(Phase ==2) %>%
                            dplyr::select(ID, SI:correctSum, -RT),
                          by = 'ID') %>%
  distinct()

plot_parms_adol <- rbind(
                    ADOL_parms_ex,
                    ADOL_parms_in) %>%
         pivot_longer(2:9, names_to = 'parms', values_to = 'est') %>%
  plyr::join(check_manipulation_full %>%
               dplyr::select(correctSum, ID),
             by = 'ID')

ggplot(plot_parms_adol %>% filter(correctSum>30),
       aes(est, fill = group))+
  geom_density(alpha = 0.5)+
  scale_fill_manual(values=colour_group_adol)+
  facet_wrap(~parms, scale = 'free', nrow = 4,
             labeller = labeller(
               parms = c(
                      alpha=expression(paste(alpha['ppt']^m)),
                      beta=expression(paste(beta['ppt']^m)),
                      alpha_v=expression(paste(alpha['ppt']^sigma)),
                      beta_v=expression(paste(beta['ppt']^sigma)),
                      alpha_ref=expression(paste(alpha['ref']^sigma)),
                      beta_ref=expression(paste(beta['ref']^sigma))
                      )
               )
             )+
  labs(y = 'Density')+
  theme(legend.position = 'none',
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        legend.background = element_blank(),
        plot.margin = margin(t=1,r=1,b=1,l=1),
        panel.grid = element_blank())

## Responsibility ----------------------------------------------------------

mod_col <- c('#CCDBDC', '#CCDBDC', '#020100', '#CCDBDC')

responADOL <- ADOL_ex$cbm[,,1]$output[,,1]$responsibility %>%
  as.data.frame()%>%
  dplyr::select(1:4)%>%
  mutate(ID = as.vector(ADOL_parms_ex$ID), group = 'Excl.', ID_n = 1:length(ID)) %>%
  rbind(.,ADOL_in$cbm[,,1]$output[,,1]$responsibility %>%
          as.data.frame()%>%
          dplyr::select(1:4)%>%
          mutate(ID = as.vector(ADOL_parms_in$ID), group = 'Inc.', ID_n = 1:length(ID))) %>%
  mutate(ID_n = c(1:251, 1:251)) %>%
  dplyr::rename('M2'=2, 'M3'=3, 'M1'=1, 'M4'=4) %>%
  pivot_longer(1:4, names_to = 'Model', values_to = 'Respon')

freqADOL <- ADOL_ex$cbm[,,1]$output[,,1]$model.frequency %>%
  as.data.frame()%>%
  dplyr::select(1:4)%>%
  mutate(group = 'Excl.') %>%
  rbind(.,ADOL_in$cbm[,,1]$output[,,1]$model.frequency %>%
          as.data.frame()%>%
          dplyr::select(1:4)%>%
          mutate(group = 'Inc.')) %>%
  dplyr::rename('M2'=2, 'M3'=3, 'M1'=1, 'M4'=4) %>%
  pivot_longer(1:4, names_to = 'Model', values_to = 'Freq')

exprobADOL <-rbind(ADOL_ex$cbm[,,1]$output[,,1]$exceedance.prob,
               ADOL_in$cbm[,,1]$output[,,1]$exceedance.prob) %>%
          as.data.frame()%>%
  dplyr::select(1:4)%>%
  mutate(group = c('Excl.', 'Inc.')) %>%
  dplyr::rename('M2'=2, 'M3'=3, 'M1'=1, 'M4'=4) %>%
  pivot_longer(1:4, names_to = 'Model', values_to = 'ExProb')

model_comparison_outcomesADOL <- responADOL %>%
  plyr::join(., freqADOL, by = c('group', 'Model')) %>%
  plyr::join(., exprobADOL, by = c('group', 'Model'))

library(patchwork)

rplotADOL <- ggplot(responADOL,
       aes(ID_n, Respon, fill = Model))+
  geom_col(colour = 'black')+
  facet_wrap(~group, scales = 'free_x')+
  scale_fill_manual(values = mod_col)+
  labs(x = 'ID', y = 'Respon.')+
  coord_cartesian(ylim = c(0,1))+
  scale_x_continuous(breaks = c(1, 10, 20, 30, 40, 50), expand = c(0,0))+
  scale_y_continuous(breaks = c(0.01, 0.5, 1), labels = c('0','0.5','1'), expand = c(0,0))+
  theme(legend.position = 'none',
        axis.text.x = element_blank())

fplotADOL <- ggplot(freqADOL,
       aes(Model, ifelse(Freq<0.01, 0+0.01, Freq), fill = Model))+
  geom_col(position = 'dodge', width = 0.5, colour = 'black')+
  labs(x = 'Model', y = 'Frequency')+
  scale_fill_manual(values=mod_col)+
  coord_cartesian(ylim = c(0,1))+
  facet_wrap(~group)+
  scale_y_continuous(breaks = c(0.01, 0.5, 1), labels = c('0','0.5','1'), expand = c(0,0))+
  theme(
        legend.position = 'none')

explotADOL <- ggplot(exprobADOL,
       aes(Model, ifelse(ExProb<0.01, 0+0.01, ExProb), fill = Model))+
  geom_col(position = 'dodge', width = 0.5, colour = 'black')+
  scale_fill_manual(values=mod_col,
                    name='Model')+
  facet_wrap(~group)+
  coord_cartesian(ylim = c(0,1))+
  scale_y_continuous(breaks = c(0.01, 0.5, 1), labels = c('0','0.5','1'), expand = c(0,0))+
  labs(y = 'Ex. Prob')+
  theme(legend.position = 'none')

comp_plotADOL <- explotADOL/fplotADOL &
  theme(plot.margin = margin(rep(1, 4), unit = 'mm'),
        text = element_text(size = 20),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid = element_blank()
  )
comp_plotADOL

summary(lm(Respon ~ group, data = responADOL %>% filter(Model=='M1')))
summary(lm(Respon ~ group, data = responADOL %>% filter(Model=='M4')))

## Recovery ------------------------------------------------

### Model Recovery ------------------------------------------------

respon_r <- ADOL_ex_recov$cbm[,,1]$output[,,1]$responsibility %>%
  as.data.frame()%>%
  dplyr::select(1:4)%>%
  mutate(ID = as.vector(IDexc$ID), Group = 'exc', ID_n = 1:length(IDexc$ID)) %>%
  rbind(.,ADOL_in_recov$cbm[,,1]$output[,,1]$responsibility %>%
          as.data.frame()%>%
          dplyr::select(1:4)%>%
          mutate(ID = as.vector(IDinc$ID), Group = 'inc', ID_n = 1:length(IDinc$ID))) %>%
  rename('M1'=1, 'M2'=2, 'M3'=3, 'M4'=4) %>%
  pivot_longer(1:4, names_to = 'Model', values_to = 'Respon')

freq_r <- ADOL_ex_recov$cbm[,,1]$output[,,1]$model.frequency %>%
  as.data.frame()%>%
  dplyr::select(1:4)%>%
  mutate(Group = 'EXC') %>%
  rbind(.,ADOL_in_recov$cbm[,,1]$output[,,1]$model.frequency %>%
          as.data.frame()%>%
          dplyr::select(1:4)%>%
          mutate(Group = 'INC')) %>%
  rename('M1'=1, 'M2'=2, 'M3'=3, 'M4'=4) %>%
  pivot_longer(1:4, names_to = 'Model', values_to = 'Freq')

exprob_r <-rbind(ADOL_ex_recov$cbm[,,1]$output[,,1]$exceedance.prob,
               ADOL_in_recov$cbm[,,1]$output[,,1]$exceedance.prob) %>%
          as.data.frame()%>%
  dplyr::select(1:4)%>%
  mutate(Group = c('EXC', 'INC')) %>%
  rename('M1'=1, 'M2'=2, 'M3'=3, 'M4'=4) %>%
  pivot_longer(1:4, names_to = 'Model', values_to = 'ExProb')

library(patchwork)

rplot_r <- ggplot(respon_r,
       aes(ID_n, Respon, fill = Model))+
  geom_col(colour = 'black')+
  facet_wrap(~Group, scales = 'free_x')+
  scale_fill_manual(values = mod_col)+
  labs(x = 'ID', y = 'Respon')+
  scale_x_continuous(breaks = c(1, 10, 20, 30, 40, 50), expand = c(0,0))+
  scale_y_continuous(breaks = c(0, 0.5, 1), expand = c(0,0))+
  theme_bw()+
  theme(text= element_text(size=18),
        legend.position = 'none',
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

fplot_r <- ggplot(freq_r,
       aes(Model, Freq, fill = Model))+
  geom_col(position = 'dodge', width = 0.5, colour = 'black')+
  labs(x = 'Model', y = 'Frequency')+
  scale_fill_manual(values = mod_col)+
  facet_wrap(~Group)+
  scale_y_continuous(breaks = c(0, 0.5, 1), expand = c(0,0), limits=c(0,1))+
  theme_bw()+
  theme(text= element_text(size=18),
        legend.position = 'none')

explot_r <- ggplot(exprob_r,
       aes(Model, ifelse(ExProb<0.01, 0+0.01, ExProb), fill = Model))+
  geom_col(position = 'dodge', width = 0.5, colour = 'black')+
  scale_fill_manual(values = mod_col)+
  scale_y_continuous(breaks = c(0, 0.5, 1), expand = c(0,0), limits=c(0,1))+
  facet_wrap(~Group)+
  labs(y = 'Ex. Prob')+
  theme_bw()+
  theme(text= element_text(size=18),
        legend.position = 'top')

model_rec <- explot_r/fplot_r &
  theme(plot.margin = margin(rep(1, 4), unit = 'mm'),
        text = element_text(size = 20),
        axis.title.x = element_blank()
  )

model_rec

### Parameter Recovery --------------------------

ADOL_ex_recov$cbm[,,1]$output[,,1]$group.mean[,3][[1]][[1]]
ADOL_in_recov$cbm[,,1]$output[,,1]$group.mean[,3][[1]][[1]]

rec_ex_parms <- as.data.frame(matrix(NA, nrow = nrow(ADOL_ex_recov$cbm[,,1]$output[,,1]$responsibility), 9))
rec_in_parms <- as.data.frame(matrix(NA, nrow = nrow(ADOL_in_recov$cbm[,,1]$output[,,1]$responsibility), 9))

rec_ex_parms[,1:8] <- ADOL_ex_recov$cbm[,,1]$output[,,1]$parameters[3,][[1]][[1]]
rec_ex_parms[,9]   <- 'ex_r'
rec_in_parms[,1:8] <- ADOL_in_recov$cbm[,,1]$output[,,1]$parameters[3,][[1]][[1]]
rec_in_parms[,9]   <- 'in_r'

rec_joint_parms_adol <- data.frame(ID      = unique(responADOL$ID),
                           rec_alpha   = c((1/(1+exp(-rec_ex_parms[,1])))*30,
                                           (1/(1+exp(-rec_in_parms[,1])))*30),
                           rec_beta    = c(rec_ex_parms[,2],
                                           rec_in_parms[,2]),
                           rec_alpha_v= c(exp(rec_ex_parms[,3]),
                                          exp(rec_in_parms[,3])),
                           rec_beta_v = c(exp(rec_ex_parms[,4]),
                                          exp(rec_in_parms[,4])),
                           rec_alpha_par= c((1/(1+exp(-rec_ex_parms[,5])))*30,
                                            (1/(1+exp(-rec_in_parms[,5])))*30),
                           rec_beta_par = c(rec_ex_parms[,6],
                                            rec_in_parms[,6]),
                           rec_alpha_ref= c(exp(rec_ex_parms[,7]),
                                            exp(rec_in_parms[,7])),
                           rec_beta_ref = c(exp(rec_ex_parms[,8]),
                                        exp(rec_in_parms[,8])),
                           group   = c(rec_ex_parms[,9],
                                       rec_in_parms[,9])) %>%
  distinct()


recovered_ps_adol <- plyr::join(rec_joint_parms_adol %>% dplyr::select(-group),
                           joint_parms_adol %>% dplyr::select(ID, alpha:group),
                           by = 'ID')

rec_cor_parms_adol  <- rcorr(as.matrix(recovered_ps_adol %>%
                                       dplyr::select(c(rec_alpha:beta_ref,
                                                       -group))
                                   ))

ggcorrplot::ggcorrplot(rec_cor_parms_adol$r[1:8, 9:16],
                       p.mat = rec_cor_parms_adol$P[1:8, 9:16],
                       lab = T,
                       colors = c( '#291720','white','#007A5A'),
                       legend.title = 'Pearson\nR', sig.level = 0.01,insig='blank',
                       pch.cex = 15,
                       pch.col = 1) +
  scale_x_discrete(labels = c(
    expression(paste(alpha['ppt']^m)),
    expression(paste(beta['ppt']^m)),
    expression(paste(alpha['ppt']^sigma)),
    expression(paste(beta['ppt']^sigma)),
    expression(paste(alpha['par']^m)),
    expression(paste(beta['par']^m)),
    expression(paste(alpha['par']^ref)),
    expression(paste(beta['par']^ref))
  ),
  expand=c(0,0)) +
    scale_y_discrete(labels = c(
    expression(paste(alpha['ppt']^m)),
    expression(paste(beta['ppt']^m)),
    expression(paste(alpha['ppt']^sigma)),
    expression(paste(beta['ppt']^sigma)),
    expression(paste(alpha['par']^m)),
    expression(paste(beta['par']^m)),
    expression(paste(alpha['par']^ref)),
    expression(paste(beta['par']^ref))
  ),
  expand=c(0,0)) +
  labs(x = 'Recovered', y = 'Real')+
  theme_bw()+
  theme(axis.title = element_text(color = 'black', size = 16),
        axis.text  = element_text(size=15),
        legend.position = 'none')

#### Scatter plot recovery ------

correlation_plot(joint_parms_adol[,2:10], rec_joint_parms_adol[,2:10], colour_group_adol)

## Group parameters --------
ADOL_list <- list(ADOL_in, ADOL_ex)
group_names <- c('in', 'ex')
group_list <- list()

# Loop through each base object and create the data frames
for (i in seq_along(ADOL_list)) {
  base_object <- ADOL_list[[i]]
  group_name <- group_names[i]
  group_df <- data.frame(
    alpha   = c((1 / (1 + exp(-base_object$cbm[,,1]$output[,,1]$group.mean[,3][[1]][[1]][,1])))*30,
                (1 / (1 + exp(-base_object$cbm[,,1]$output[,,1]$group.hierarchical.errorbar[,3][[1]][[1]][,1])))),
    beta    = c(base_object$cbm[,,1]$output[,,1]$group.mean[,3][[1]][[1]][,2],
                base_object$cbm[,,1]$output[,,1]$group.hierarchical.errorbar[,3][[1]][[1]][,2]),
    alpha_v = c(exp(base_object$cbm[,,1]$output[,,1]$group.mean[,3][[1]][[1]][,3]),
                exp(base_object$cbm[,,1]$output[,,1]$group.hierarchical.errorbar[,3][[1]][[1]][,3])),
    beta_v  = c(exp(base_object$cbm[,,1]$output[,,1]$group.mean[,3][[1]][[1]][,4]),
                exp(base_object$cbm[,,1]$output[,,1]$group.hierarchical.errorbar[,3][[1]][[1]][,4])),
    alpha_par = c((1 / (1 + exp(-base_object$cbm[,,1]$output[,,1]$group.mean[,3][[1]][[1]][,5])))*30,
                (1 / (1 + exp(-base_object$cbm[,,1]$output[,,1]$group.hierarchical.errorbar[,3][[1]][[1]][,5])))),
    beta_par  = c(base_object$cbm[,,1]$output[,,1]$group.mean[,3][[1]][[1]][,6],
                  base_object$cbm[,,1]$output[,,1]$group.hierarchical.errorbar[,3][[1]][[1]][,6]),
    alpha_r = c(exp(base_object$cbm[,,1]$output[,,1]$group.mean[,3][[1]][[1]][,7]),
                exp(base_object$cbm[,,1]$output[,,1]$group.hierarchical.errorbar[,3][[1]][[1]][,7])),
    beta_r  = c(exp(base_object$cbm[,,1]$output[,,1]$group.mean[,3][[1]][[1]][,8]),
                exp(base_object$cbm[,,1]$output[,,1]$group.hierarchical.errorbar[,3][[1]][[1]][,8])),
    group   = group_name,
    type    = c('mean', 'error')
  )

  group_list[[group_name]] <- group_df
}


parm_group_adol <- rbind(group_list[[1]], group_list[[2]]) %>%
  pivot_longer(1:8, names_to = 'Parm', values_to = 'Val') %>%
  mutate(category = ifelse(str_detect(Parm, "beta"), "beta", "alpha"))

parm_group_adol$Parm <- factor(parm_group_adol$Parm,
                          levels = c("alpha", "beta", "alpha_v", "beta_v",
                                     'alpha_par', 'beta_par', "alpha_r", "beta_r"))

# Separate the data into mean and error datasets
means_a  <- parm_group_adol %>% filter(type == "mean")
errors_a <- parm_group_adol %>% filter(type == "error")

# Merge the datasets on group, Parm, and category
parm_group_full_adol <- means_a %>%
  left_join(errors_a, by = c("group", "Parm", "category"), suffix = c("_mean", "_error"))

ggplot(parm_group_full_adol,
       aes(Parm, Val_mean, fill = group)) +
  geom_col(position = 'dodge')+
  geom_errorbar(aes(ymin = Val_mean - Val_error, ymax = Val_mean + Val_error),
                width = 0.2, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values=colour_group_adol)+
  facet_wrap(~category, scale = 'free')+
  scale_y_continuous(expand=c(0,0.2))+
  scale_x_discrete(labels = function(x) {
    ifelse(str_detect(x, "beta"),
           c(expression(paste(beta['ppt']^m)),
             expression(paste(beta['ppt']^sigma)),
             expression(paste(beta['par']^m)),
             expression(paste(beta['par']^sigma))
             ),
           c(expression(paste(alpha['ppt']^m)),
             expression(paste(alpha['ppt']^sigma)),
             expression(paste(alpha['par']^m)),
             expression(paste(alpha['par']^sigma))
             ))
  },
  expand = c(0,0.5)) +
  theme(legend.position = c(0.9, 0.23),
        legend.background = element_rect(colour = 'black'),
        strip.text.x = element_blank(),
        legend.title = element_blank(),
        axis.title.x = element_blank())

# Distributions -----------------------------------------------------------
# Parameters
lim <- 30
res <- 0.25

# Process BPD and CON data
ADOLdist_in <- process_dist_data(ADOL_in_full_sim, unique(ID_in), 'in', 'M3', lim, res)
ADOLdist_ex <- process_dist_data(ADOL_ex_full_sim, unique(ID_ex), 'ex', 'M3', lim, res)

# Combine the results
int_dist_adol <- rbind(ADOLdist_in, ADOLdist_ex)

# Check the result
head(int_dist_adol)

## Individual Distribution Examples For Each Model -----------

ggplot(int_dist_adol %>% filter(ID == unique(int_dist_adol$ID)[2]))+
  geom_line(aes(beta, bP1, colour = '1 - P1'), size = 1.2)+
  geom_line(aes(beta, bP2a, colour = '2A - Prior'), size = 1.2)+
  geom_line(aes(beta, bP2b, colour = '2B - Post.'), size = 1.2)+
  geom_line(aes(beta, bP3, colour = '3 - Post.'), size = 1.2)+
  labs(title = '', x = expression(paste(beta)), y = expression(paste('p(',beta,')')))+
  scale_colour_manual(values = c('black','#FBD19D', '#F7941D', 'grey'))+
  theme(legend.position = 'none') &
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = 'white'),
        plot.title = element_text(hjust = 0.5),
        text=element_text(size=16),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin=margin(1,1,1,1),
        axis.title.y =element_text(vjust=-10),
        panel.border = element_blank())

## Entropy & KL -----------

entropy_group_adol <- rbind(
  int_dist_adol %>%
  group_by(ID) %>%
  mutate(P1           = calculate_entropy(bP1),
         P2a          = calculate_entropy(bP2a),
         P2b          = calculate_entropy(bP2b),
         P3           = calculate_entropy(bP3),
         type = 'Beta') %>%
  dplyr::select(ID, group, Model,
                P1, P2a, P2b, P3, type) %>%
  pivot_longer(P1:P3, names_to = 'Phase', values_to = 'Entropy') %>%
  mutate(Phase = factor(Phase, levels = c('P1', 'P2a', 'P2b', 'P3'))) %>%
  distinct(),
  int_dist_adol %>%
  group_by(ID) %>%
  mutate(P1           = calculate_entropy(aP1),
         P2a          = calculate_entropy(aP2a),
         P2b          = calculate_entropy(aP2b),
         P3           = calculate_entropy(aP3),
         type = 'Alpha') %>%
  dplyr::select(ID, group, Model,
                P1, P2a, P2b, P3, type) %>%
  pivot_longer(P1:P3, names_to = 'Phase', values_to = 'Entropy') %>%
  mutate(Phase = factor(Phase, levels = c('P1', 'P2a', 'P2b', 'P3'))) %>%
  distinct()
)

ggplot(entropy_group_adol %>%
         dplyr::select(ID, Group=group, Phase, Entropy, type) %>%
         distinct(),
       aes(Phase, Entropy, fill = Group))+
  geom_jitter(shape=21, width = 0.1, alpha = 0.3)+
  stat_summary(geom='line', aes(group = Group))+
  stat_summary(shape=21, size = 1)+
  stat_compare_means(label = 'p.signif', method = 't.test',
                     size = 8, label.y.npc = 0.15, hide.ns = T)+
  scale_fill_manual(values=colour_group_adol)+
  facet_wrap(~type)+
  labs(x = 'Phase', y = 'Entropy (Inf. Surprise)')+
  theme(panel.grid = element_blank(),
        legend.position = 'top',
        legend.title = element_blank())

## Mean value shift ------

mean_val_shift_adol <- int_dist_adol %>%
  group_by(ID) %>%
  mutate(B1  = sum(beta*bP1),
         B2a = sum(beta*bP2a),
         B2b = sum(beta*bP2b),
         B3  = sum(beta*bP3),
         A1  = sum(alpha*aP1),
         A2a = sum(alpha*aP2a),
         A2b = sum(alpha*aP2b),
         A3  = sum(alpha*aP3)) %>%
  dplyr::select(ID: A3) %>%
  distinct()

p1 <- ggplot(mean_val_shift_adol %>%
               mutate(BinComp = ifelse(B1 > 0, 'Competitive', 'Prosocial')) %>%
               pivot_longer(4:5, names_to = 'Phase', values_to = 'MeanBeta'),
             aes(Phase, MeanBeta, fill = BinComp)) +
  geom_hline(yintercept = 0, colour = 'grey', linetype = 2) +
  geom_jitter(aes(colour = BinComp), alpha = 0.1) +
  stat_summary(geom = 'line', colour = 'black', aes(group = BinComp), size = 1.1) +
  stat_summary(shape = 21, colour = 'black', size = 1) +
  coord_cartesian(ylim = c(-25, 20)) +
  scale_x_discrete(labels = c('Phase 1\nPreferences', 'Phase 2\nPrior'))+
  scale_fill_brewer(palette = 'Dark2') +
  scale_colour_brewer(palette = 'Dark2')+
  labs(x = 'Phase', y = expression(paste(beta^m)))

p2 <- ggplot(mean_val_shift_adol %>%
               mutate(BinComp = ifelse(B1 > 0, 'C', 'P'),
                      BDelta = B3 - B1),
             aes(BinComp, BDelta, fill = BinComp)) +
  geom_hline(yintercept = 0, colour = 'black') +
  geom_jitter(aes(colour = BinComp), alpha = 0.5, width = 0.2) +
  stat_summary(geom = 'bar', colour = 'black', size = 0.5, alpha = 0.7) +
  coord_cartesian(ylim = c(-25, 20)) +
  labs(x = '', y = expression(paste(Delta, beta[ppt]^m, ' [B1-B3]'))) +
  scale_fill_brewer(palette = 'Dark2') +
  scale_colour_brewer(palette = 'Dark2')

# Combine the plots side-by-side
p1 + p2 &
  theme_bw(base_size = 24) &
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = 'white'),
        legend.title = element_blank(),
        legend.position = 'none',
        axis.title.x = element_blank(),
        plot.margin = margin(1,1,1,1)) &
  plot_layout(widths = c(5, 1))

ggplot(mean_val_shift_adol %>%
               mutate(BinComp = ifelse(B1 > 0, 'Competitive', 'Prosocial')) %>%
               pivot_longer(4:5, names_to = 'Phase', values_to = 'MeanBeta'),
             aes(Phase, MeanBeta, group = ID, colour = MeanBeta)) +
  geom_hline(yintercept = 0, colour = 'grey', linetype = 2) +
  geom_line(colour = 'grey', alpha = 0.1)+
  geom_jitter(width = 0.1)+
  scale_colour_gradient(low = "#D95F02", high = "#1B9E77")+
  labs(x = 'Phase', y = expression(paste(beta[ppt]^m)))+
  theme_bw(base_size = 24) +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = 'white'),
        legend.title = element_blank(),
        legend.position = 'none')

mean_val_shift_adol %>%
  mutate(BinComp = ifelse(B1 > 0, 'Competitive', 'Prosocial')) %>%
  pivot_longer(B1:B2a, names_to = 'Phase', values_to = 'Mean') %>%
  filter(BinComp=='Prosocial') %>%
  dplyr::select(Phase, Mean, BinComp, ID) %>%
  distinct() %>%
  glm(Mean ~ Phase, data = .) %>%
  summary()

mean_val_shift_adol %>%
  mutate(BinComp = ifelse(B1 > 0, 'Competitive', 'Prosocial')) %>%
  pivot_longer(B1:B2a, names_to = 'Phase', values_to = 'Mean') %>%
  filter(BinComp=='Competitive') %>%
  dplyr::select(Phase, Mean, BinComp, ID) %>%
  distinct() %>%
  glm(Mean ~ Phase, data = .) %>%
  confint()

mean_val_shift_adol %>%
  pivot_longer(A1:A2a, names_to = 'Phase', values_to = 'Mean') %>%
  dplyr::select(Phase, Mean, ID) %>%
  distinct() %>%
  glm(Mean ~ Phase, data = .) %>%
  confint()

## Belief Updates ----

b_updates_adol <- check_manipulation_full %>%
  filter(Phase==2) %>%
  dplyr::select(ID, trial, correctSum, Context, SI, HI, agency, Q6, Q7) %>%
  mutate(group = ifelse(Context=='Exc', 'ex', 'in'),
         kl_div_a = 0,
         kl_div_b = 0,
         trial=1:54)

for(k in 1:3){
  if(k == 1){x = ADOL_in_full_sim$results; group = 'in'; ID = unique(ID_in)}
  if(k == 2){x = ADOL_ex_full_sim$results; group = 'ex'; ID = unique(ID_ex)}
  for(j in 1:length(ID)){
    kl_divsa = rep(NA, 54)
    kl_divsb = kl_divsa
    for(i in 1:54){
      x[j,][[1]][[1]][,,1]$alpha.cont[,36] = x[j,][[1]][[1]][,,1]$alpha.marg1
      x[j,][[1]][[1]][,,1]$beta.cont[,36]  = x[j,][[1]][[1]][,,1]$beta.marg1

      bsa            = x[j,][[1]][[1]][,,1]$alpha.cont[,36:90]
      b_t2a          = bsa[,i+1]
      b_t1a          = bsa[,i]
      bsb            = x[j,][[1]][[1]][,,1]$beta.cont[,36:90]
      b_t2b          = bsb[,i+1]
      b_t1b          = bsb[,i]
      kl_divsa[i]    = calculate_KL_divergence(b_t2a, b_t1a)
      kl_divsb[i]    = calculate_KL_divergence(b_t2b, b_t1b)
    }
    b_updates_adol[b_updates_adol$group==group & b_updates_adol$ID==ID[j],'kl_div_a'][1:54,] <- kl_divsa
    b_updates_adol[b_updates_adol$group==group & b_updates_adol$ID==ID[j],'kl_div_b'][1:54,] <- kl_divsb
  }
}

b_updates_adol <- b_updates_adol %>%
  mutate(kl_div_aroll = rollapply(kl_div_a, width = 5, FUN = mean, align = "right", fill = NA),
         kl_div_broll = rollapply(kl_div_b, width = 5, FUN = mean, align = "right", fill = NA)) %>%
  plyr::join(., adolIntent_clean %>% dplyr::select(ID, server_alpha_ppt:server_beta_par), by = 'ID') %>%
         distinct() %>%
         mutate(diffB  = abs(server_beta_ppt - server_beta_par),
                diffBQ = ifelse(diffB < 11, 'low',
                               ifelse(diffB >= 11 & diffB < 19.25, 'mid',
                                      ifelse(diffB >= 19.25 & diffB <= 27, 'high', NA))))

facet_labels <- c('alpha', 'beta')
names(facet_labels) <- c('kl_div_a', 'kl_div_b')
ggplot(b_updates_adol %>%
         filter(trial !=1) %>%
         pivot_longer(11:12, names_to = 'Parameter', values_to = 'KL_Div'),
       aes(trial, KL_Div, fill = factor(group), colour = factor(group)))+
  stat_smooth(alpha = 0.01, aes(group = ID), linewidth = 0.05, se = F)+
  stat_smooth(alpha = 0.7, colour = 'black')+
  scale_fill_manual(values=colour_group_adol)+
  scale_colour_manual(values=colour_group_adol)+
  coord_cartesian(ylim = c(0, 0.016))+
  scale_x_continuous(breaks = c(2, 25, 50), expand = c(0,0))+
  scale_y_continuous(breaks = c(0, 0.008, 0.016), labels = c('0', '0.008', '0.016'), expand = c(0.01,0))+
  facet_wrap(~Parameter, nrow = 2) +
  labs(y = expression(paste('D'[KL])), x = 'Trial')+
  theme_bw(base_size = 20)+
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        legend.title = element_blank(),
        strip.text.x = element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(10,1,1,1))

lmerTest::lmer(scale(kl_div_b) ~ scale(diffB) * trial + (1|ID), data = b_updates_adol) %>% summary()

# Parameter analysis cleaning ----------------------------------------------------------

#Calculate final distribution at end of Phase 2

par_parms_adol <- int_dist_adol %>%
  group_by(ID) %>%
  mutate(beta_par_m = sum(beta*bP2b),
         alpha_par_m = sum(alpha*aP2b),
         beta_par_pri = sum(beta*bP2a),
         alpha_par_pri = sum(alpha*aP2a),
         shift_beta_p2 = abs(beta_par_pri-beta_par_m),
         shift_alpha_p2 = abs(alpha_par_pri-alpha_par_m),
         beta_hat = sum(beta*bP3),
         alpha_hat= sum(alpha*aP3),
         beta_hat_sd = sqrt(sum((beta-beta_hat)^2 * bP3)),
         alpha_hat_sd = sqrt(sum((alpha-alpha_hat)^2 * aP3))) %>%
  dplyr::select(ID, group,
                beta_par_m, alpha_par_m,
                beta_par_pri, alpha_par_pri,
                beta_hat, alpha_hat,
                beta_hat_sd, alpha_hat_sd,
                shift_beta_p2, shift_alpha_p2) %>%
  distinct() %>%
  plyr::join(., adolIntent_clean %>%
         dplyr::select(ID, server_beta_par, server_alpha_par)) %>%
  mutate(disp_beta_par = abs(server_beta_par - beta_par_m),
         disp_alpha_par = abs(server_alpha_par - alpha_par_m)) %>%
  dplyr::select(ID,
                beta_par_m, alpha_par_m,
                beta_par_pri, alpha_par_pri,
                disp_beta_par, disp_alpha_par,
                beta_hat_sd, alpha_hat_sd,
                beta_hat, alpha_hat,
                shift_beta_p2, shift_alpha_p2) %>%
  distinct()

joint_parms_ex_adol <- plyr::join(par_parms_adol, joint_parms_adol, by = 'ID') %>%
                       plyr::join(., check_manipulation_tidy %>%
                                    dplyr::select(ID, Q6, Q7),
                                  by = 'ID') %>%
  mutate(disp_beta_ppt = abs(beta_hat - beta),
         disp_alpha_ppt = abs(alpha_hat - alpha),
         disp_beta_hat_sd = abs(beta_hat_sd - beta_v),
         disp_alpha_hat_sd = abs(alpha_hat_sd - alpha_v),
         distancebeta=abs(server_beta_par-server_beta_ppt),
         distancealpha=abs(server_alpha_par-server_alpha_ppt))


## Check for Self Change Confirmation --------------------------------

lim = 30
res = 0.25

for (i in 1:length(IDexc$ID)){
  xdist <-   data.frame(
    bP1  = as.vector(ADOL_ex_full_sim$results[i,][[1]][[1]][,,1]$beta.marg1),
    bP2a = as.vector(ADOL_ex_full_sim$results[i,][[1]][[1]][,,1]$beta.marg2a),
    bP2b = as.vector(ADOL_ex_full_sim$results[i,][[1]][[1]][,,1]$beta.marg2b),
    bP3  = as.vector(ADOL_ex_full_sim$results[i,][[1]][[1]][,,1]$beta.marg3),
    aP1  = as.vector(ADOL_ex_full_sim$results[i,][[1]][[1]][,,1]$alpha.marg1),
    aP2a  = as.vector(ADOL_in_full_sim$results[i,][[1]][[1]][,,1]$alpha.marg2a),
    aP2b  = as.vector(ADOL_in_full_sim$results[i,][[1]][[1]][,,1]$alpha.marg2b),
    aP3  = as.vector(ADOL_ex_full_sim$results[i,][[1]][[1]][,,1]$alpha.marg3),
    beta = seq(-lim, lim, res),
    alpha = seq(0, lim, res/2),
    ID = IDexc$ID[i],
    group = 'ex',
    Model = 'M3'
  )
  if(i>1){EXCdist <- rbind(EXCdist, xdist)} else {EXCdist <- xdist}
}
for (i in 1:length(IDinc$ID)){
  zdist <-   data.frame(
    bP1  = as.vector(ADOL_in_full_sim$results[i,][[1]][[1]][,,1]$beta.marg1),
    bP2a = as.vector(ADOL_in_full_sim$results[i,][[1]][[1]][,,1]$beta.marg2a),
    bP2b = as.vector(ADOL_in_full_sim$results[i,][[1]][[1]][,,1]$beta.marg2b),
    bP3  = as.vector(ADOL_in_full_sim$results[i,][[1]][[1]][,,1]$beta.marg3),
    aP1  = as.vector(ADOL_in_full_sim$results[i,][[1]][[1]][,,1]$alpha.marg1),
    aP2a  = as.vector(ADOL_in_full_sim$results[i,][[1]][[1]][,,1]$alpha.marg2a),
    aP2b  = as.vector(ADOL_in_full_sim$results[i,][[1]][[1]][,,1]$alpha.marg2b),
    aP3  = as.vector(ADOL_in_full_sim$results[i,][[1]][[1]][,,1]$alpha.marg3),
    beta  = seq(-lim, lim, res),
    alpha = seq(0, lim, res/2),
    ID = IDinc$ID[i],
    group = 'in',
    Model = 'M3'
  )
  if(i>1){INCdist <- rbind(INCdist, zdist)} else {INCdist <- zdist}
}

int_dist_conv <- rbind(EXCdist, INCdist)

mean_val_shift_adol <- int_dist_conv %>%
  group_by(ID) %>%
  mutate(B1  = sum(beta*bP1),
         B2a = sum(beta*bP2a),
         B2b = sum(beta*bP2b),
         B3  = sum(beta*bP3),
         A1  = sum(alpha*aP1),
         A2a = sum(alpha*aP2a),
         A2b = sum(alpha*aP2b),
         A3  = sum(alpha*aP3),
         B1_sd = sqrt(sum((beta-B1)^2 * bP1)),
         B2a_sd= sqrt(sum((beta-B2a)^2 * bP2a)),
         B2b_sd= sqrt(sum((beta-B2b)^2 * bP2b)),
         B3_sd = sqrt(sum((beta-B3)^2 * bP3)),
         A1_sd = sqrt(sum((alpha-A1)^2* aP1)),
         A2a_sd = sqrt(sum((alpha-A2a)^2* aP2a)),
         A2b_sd = sqrt(sum((alpha-A2b)^2* aP2b)),
         A3_sd = sqrt(sum((alpha-A3)^2* aP3)),
         Delta_B    = B3-B1,
         Delta_A    = A3-A1,
         Delta_B_sd = B3_sd - B1_sd,
         Delta_A_sd = A3_sd - A1_sd) %>%
  dplyr::select(ID: Delta_A_sd) %>%
  pivot_longer(Delta_B:Delta_A_sd, names_to = 'Delta', values_to = 'Shift')  %>%
  plyr::join(., PPT_PARcong_adol %>% dplyr::select(ID, matched_tot), by = 'ID') %>%
  distinct()

ggplot(mean_val_shift_adol %>% filter(Delta%in%c('Delta_B', 'Delta_A')),
       aes(Delta, Shift, fill = group))+
  geom_point(shape=21, alpha=0.3, size=3,position = position_dodge(width=0.5))+
  geom_boxplot(width=0.5, outliers = F, alpha = 0.7)+
  labs(x='', y=expression(paste(Delta, theta['ppt']^m)))+
  scale_x_discrete(labels = c(expression(paste(alpha['ppt']^m)), expression(paste(beta['ppt']^m))))+
  scale_fill_manual(values=colour_group_adol)+
  stat_compare_means(label = 'p.signif', size = 8, label.y = 15)+
  theme_bw(base_size=22)+
ggplot(mean_val_shift_adol %>% filter(Delta%in%c('Delta_B_sd', 'Delta_A_sd')),
       aes(Delta, Shift, fill = group))+
  geom_point(shape=21, alpha=0.3, size=3,position = position_dodge(width=0.5))+
  geom_boxplot(width=0.5, outliers = F, alpha = 0.7)+
  labs(x='', y=expression(paste(Delta, theta['ppt']^sigma)))+
  scale_x_discrete(labels = c(expression(paste(alpha['ppt']^sigma)), expression(paste(beta['ppt']^sigma))))+
  scale_fill_manual(values=colour_group_adol)+
  stat_compare_means(label = 'p.signif', size = 8, label.y = -15)+
  theme_bw(base_size=22)&
  theme(axis.title.x = element_blank(),
        legend.position = 'none',
        panel.grid = element_blank(),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=22),
        plot.margin = margin(t = 1, r = 1, b = 1, l = 1),
        axis.title.y = element_text(vjust=-1))

# mean shift contribution

ex_change   <- rstanarm::stan_glm(abs(Shift) ~ B1 + B1_sd + B2a_sd + B2b + B2b_sd, data = mean_val_shift_adol %>% filter(Delta%in%c('Delta_B'), group == 'ex'))
in_change   <- rstanarm::stan_glm(abs(Shift) ~ B1 + B1_sd + B2a_sd + B2b + B2b_sd, data = mean_val_shift_adol %>% filter(Delta%in%c('Delta_B'), group == 'in'))

mean_changeB <- rstanarm::stan_glm(abs(Shift) ~ B1 + B1_sd + B2a_sd + B2b + B2b_sd, data = mean_val_shift_adol %>% filter(Delta%in%c('Delta_B')))
mean_changeA <- rstanarm::stan_glm(abs(Shift) ~ A1 + A1_sd + A2a_sd + A2b + A2b_sd, data = mean_val_shift_adol %>% filter(Delta%in%c('Delta_B')))

summary(mean_changeB, prob = c(0.025, 0.975))
summary(mean_changeA, prob = c(0.025, 0.975))

summary(glm(B2b_sd ~ group, data = mean_val_shift_adol %>% filter(Delta%in%c('Delta_B'))))
summary(glm(B1_sd ~ group, data = mean_val_shift_adol %>% filter(Delta%in%c('Delta_B'))))

m_shift_plot <- rbind(
  as.data.frame(ex_change) %>% mutate(Group = 'Exc.'),
  as.data.frame(in_change) %>% mutate(Group = 'Inc.')
) %>%
  pivot_longer(c(B1, B1_sd, B2a_sd, B2b, B2b_sd), names_to = 'Par', values_to = 'Shift') %>%

  ggplot(aes(Shift, fct_reorder(Par, Shift), colour = Group))+
    tidybayes::geom_eyeh(fill = NA, size = 10)+
    geom_vline(xintercept = 0)+
    scale_colour_manual(values = colour_group_adol)+
    scale_y_discrete(labels = c(expression(paste(beta[par]^ref)),
                                expression(paste(beta[par]^sigma)),
                                expression(paste(beta[ppt]^m)),
                                expression(paste(beta[par]^m)),
                                expression(paste(beta[ppt]^sigma))))+
    labs(y = 'Parameter', x = expression(paste('|',Delta,beta[ppt]^m, '|')))+
    scale_x_continuous(breaks = c(-5, 0, 5))+
    tidybayes::theme_tidybayes()+
    theme(legend.position = 'none',
          text = element_text(size=24),
          legend.title = element_blank())
m_shift_plot

#ratio test
ratio_shift <- mean_val_shift_adol %>%
  mutate(B_sd_ratio = B1_sd/(B1_sd+B2a_sd))

ex_change_rat <- rstanarm::stan_glm(abs(Shift) ~ B_sd_ratio, data = ratio_shift %>% filter(Delta%in%c('Delta_B'), group == 'ex'))
in_change_rat <- rstanarm::stan_glm(abs(Shift) ~ B_sd_ratio, data = ratio_shift %>% filter(Delta%in%c('Delta_B'), group == 'in'))

mean_changeBrat <- rstanarm::stan_glm(abs(Shift) ~ B_sd_ratio, data = ratio_shift %>% filter(Delta%in%c('Delta_B')))
mean_changeArat <- rstanarm::stan_glm(abs(Shift) ~ B_sd_ratio, data = ratio_shift %>% filter(Delta%in%c('Delta_B')))

summary(mean_changeBrat, prob = c(0.025, 0.975))
summary(mean_changeArat, prob = c(0.025, 0.975))

rbind(
  as.data.frame(ex_change_rat) %>% mutate(Group = 'Shift'),
  as.data.frame(in_change_rat) %>% mutate(Group = 'Shift')
) %>%

  ggplot(aes(Group, B_sd_ratio))+
    ggdist::stat_slabinterval(alpha = 1, fill= NA, size = 10)+
    #ggdist::stat_dotsinterval(fill = 'grey', colour = NA, quantiles = 1000)+
    #stat_summary(colour = 'black')+
    #scale_colour_manual(values = colour_group_adol)+
    labs(y = expression(paste(beta[ppt]^sigma, ':', beta[par]^ref)), x = expression(paste(Delta, beta[ppt]^m)))+
    ggpubr::stat_compare_means(paired = F, label = 'p.signif', size = 8, label.x.npc = 0.5, label.y = 1)+
    scale_y_continuous(expand = c(0,0))+
    coord_cartesian(ylim = c(0, 20))+
    tidybayes::theme_tidybayes()+
    theme_bw()+
    theme(legend.position = 'none',
          text = element_text(size=24),
          legend.title = element_blank(),
          panel.grid = element_blank(),
          axis.text.x = element_blank())

#sd shift contribution

ex_change_sd <- rstanarm::stan_glm(abs(Shift) ~ B1 + B1_sd + B2a_sd + B2b + B2b_sd, data = mean_val_shift_adol %>% filter(Delta%in%c('Delta_B_sd'), group == 'ex'))
in_change_sd <- rstanarm::stan_glm(abs(Shift) ~ B1 + B1_sd + B2a_sd + B2b + B2b_sd, data = mean_val_shift_adol %>% filter(Delta%in%c('Delta_B_sd'), group == 'in'))

summary(glm(B2a_sd ~ group, data = mean_val_shift_adol %>% filter(Delta%in%c('Delta_B_sd'))))
summary(glm(B1_sd ~ group, data = mean_val_shift_adol %>% filter(Delta%in%c('Delta_B_sd'))))

sd_shift_plot <- rbind(
  as.data.frame(ex_change_sd) %>% mutate(Group = 'Exc.'),
  as.data.frame(in_change_sd) %>% mutate(Group = 'Inc.')
) %>%
  pivot_longer(c(B1, B1_sd, B2a_sd, B2b, B2b_sd), names_to = 'Par', values_to = 'Shift') %>%

  ggplot(aes(Shift, fct_reorder(Par, Shift), colour = Group))+
    tidybayes::geom_eyeh(fill = NA, size = 10)+
    geom_vline(xintercept = 0)+
    scale_colour_manual(values = colour_group_adol)+
    scale_y_discrete(labels = c(expression(paste(beta[par]^ref)),
                                expression(paste(beta[par]^sigma)),
                                expression(paste(beta[ppt]^m)),
                                expression(paste(beta[par]^m)),
                                expression(paste(beta[ppt]^sigma))))+
    labs(y = 'Parameter', x = expression(paste('|',Delta,beta[ppt]^sigma, '|')))+
    scale_x_continuous(breaks = c(-0.5, 0, 0.5))+
    tidybayes::theme_tidybayes()+
    theme(legend.position = c(0.75, 0.5),
          text = element_text(size=24),
          legend.title = element_blank(),
          legend.background = element_rect(colour = 'black'),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          axis.text.y = element_blank(),
          plot.margin = margin(7, 25, 7, 25))

ggarrange(m_shift_plot, sd_shift_plot)

#ratio test sd
ex_change_rat_sd <- rstanarm::stan_glm(abs(Shift) ~ B_sd_ratio, data = ratio_shift %>% filter(Delta%in%c('Delta_B_sd'), group == 'ex'))
in_change_rat_sd <- rstanarm::stan_glm(abs(Shift) ~ B_sd_ratio, data = ratio_shift %>% filter(Delta%in%c('Delta_B_sd'), group == 'in'))

rbind(
  as.data.frame(ex_change_rat_sd) %>% mutate(Group = 'Shift'),
  as.data.frame(in_change_rat_sd) %>% mutate(Group = 'Shift')
) %>%

  ggplot(aes(Group, B_sd_ratio))+
    ggdist::stat_slabinterval(alpha = 1, fill= NA, size = 10)+
    #ggdist::stat_dotsinterval(fill = 'grey', colour = NA, quantiles = 1000)+
    #stat_summary(colour = 'black')+
    #scale_colour_manual(values = colour_group_adol)+
    labs(y = expression(paste(beta[ppt]^sigma, ':', beta[par]^ref)), x = expression(paste(Delta, beta[ppt]^sigma)))+
    ggpubr::stat_compare_means(paired = F, label = 'p.signif', size = 8, label.x.npc = 0.5, label.y = 1)+
    scale_y_continuous(expand = c(0,0))+
    coord_cartesian(ylim = c(0, 7))+
    tidybayes::theme_tidybayes()+
    theme_bw()+
    theme(legend.position = 'none',
          text = element_text(size=24),
          legend.title = element_blank(),
          panel.grid = element_blank(),
          axis.text.x = element_blank())

#overall shift by context
rstanarm::stan_glm(abs(Shift) ~ group, data = ratio_shift %>% filter(Delta%in%c('Delta_B'))) %>%
  as.data.frame() %>%
  ggplot(aes(groupin))+
  geom_histogram(colour = 'black', fill = 'darkred')+
  geom_vline(xintercept = 0, size = 2) +
  labs(x = expression(paste('|',Delta,beta[ppt]^m, '|  Inc. vs Exc.')))+
  theme_bw(base_size = 18)+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank())

rstanarm::stan_glm(-abs(Shift) ~ group, data = ratio_shift %>% filter(Delta%in%c('Delta_B_sd'))) %>%
  as.data.frame() %>%
  ggplot(aes(groupin))+
  geom_histogram(colour = 'black', fill = 'darkred')+
  geom_vline(xintercept = 0, size = 2) +
  labs(x = expression(paste('|',Delta,beta[ppt]^sigma, '|  Inc. vs. Exc.')))+
  theme_bw(base_size = 18)+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank())

# Check against model predictions

shift_indiv  <- matrix(NA,
                      nrow = nrow(joint_parms_ex_adol),
                      ncol = 9)

for(i in 1:nrow(joint_parms_ex_adol)){

  datuse <- adolIntent_clean[adolIntent_clean$ID==joint_parms_ex_adol$ID[i],]
  parms  <- joint_parms_ex_adol[i,14:21]

  print(round(i/nrow(joint_parms_ex_adol), 3)*100)
  run1  <- ABA_shift_Gen(as.numeric(parms),
                        na.omit(datuse),
                        sim=1,
                        plot = 0,
                        proj=0)

  shift_indiv[i,1] <- as.numeric(run1$Shift[1])
  shift_indiv[i,2] <- as.numeric(run1$Shift[2])
  shift_indiv[i,3] <- as.numeric(run1$Shift[3])
  shift_indiv[i,4] <- as.numeric(run1$Shift[4])
  shift_indiv[i,5] <- as.numeric(parms[,3])
  shift_indiv[i,6] <- as.numeric(parms[,4])
  shift_indiv[i,7] <- as.numeric(parms[,7])
  shift_indiv[i,8] <- as.numeric(parms[,8])
  shift_indiv[i,9] <- joint_parms_ex_adol$ID[i]

}

colnames(shift_indiv) <- c('m_shift_alpha', 'm_shift_beta', 'sd_shift_alpha', 'sd_shift_beta', 's_u_a', 's_u_b','o_u_a', 'o_u_b','ID')

shift_l_indiv <- shift_indiv %>%
  as.data.frame() %>%
  pivot_longer(1:4, names_to = 'parameter', values_to = 'shift')

delta_shift_real <- joint_parms_ex_adol %>%
  mutate(
         delta_alpha_ppt_m = abs(alpha - alpha_hat),
         delta_beta_ppt_m = abs(beta - beta_hat),
         delta_alpha_ppt_sd= -abs(alpha_v - alpha_hat_sd),
         delta_beta_ppt_sd= -abs(beta_v - beta_hat_sd),
         beta_v = round(beta_v),
         beta_ref = round(beta_ref)) %>%
  dplyr::select(delta_alpha_ppt_m, delta_beta_ppt_m, delta_alpha_ppt_sd, delta_beta_ppt_sd, beta_v, beta_ref, ID)

delta_shift_model <- shift_l_indiv %>%
  pivot_wider(id_cols = c(s_u_a, s_u_b, o_u_a, o_u_b, ID), names_from = parameter, values_from = shift) %>%
  rename(alpha_v = s_u_a,
         beta_v = s_u_b,
         alpha_ref = o_u_a,
         beta_ref = o_u_b,
         delta_alpha_ppt_m = m_shift_alpha,
         delta_beta_ppt_m = m_shift_beta,
         delta_alpha_ppt_sd = sd_shift_alpha,
         delta_beta_ppt_sd = sd_shift_beta) %>%
  mutate(alpha_v = as.numeric(alpha_v),
         beta_v = as.numeric(beta_v),
         alpha_ref = as.numeric(alpha_ref),
         beta_ref = as.numeric(beta_ref),
         delta_alpha_ppt_m  = as.numeric(delta_alpha_ppt_m),
         delta_beta_ppt_m   = as.numeric(delta_beta_ppt_m),
         delta_alpha_ppt_sd = as.numeric(delta_alpha_ppt_sd),
         delta_beta_ppt_sd  = as.numeric(delta_beta_ppt_sd))

delta_shift <- left_join(delta_shift_real, delta_shift_model,
                         by = c('ID'),
                         suffix = c("_real", "_simulated"))

comparison_plot <- ggplot(delta_shift,
       aes(delta_beta_ppt_sd_real, delta_beta_ppt_sd_simulated)) +
  geom_point(shape=21, size = 3, fill = 'grey')+
  geom_smooth(method='lm', colour = 'black')+
  scale_x_continuous(breaks = c(0, -5, -10))+
  scale_y_continuous(breaks = c(0, -5, -10))+
  scale_alpha_continuous(range = c(0.2,1))+
  theme_bw(base_size=18)+
  labs(x=expression(paste('Real ',Delta, beta[ppt]^sigma)),
       y=expression(paste('Sim. ',Delta, beta[ppt]^sigma)))+
ggplot(delta_shift,
       aes(delta_beta_ppt_m_real, abs(delta_beta_ppt_m_simulated))) +
  geom_point(shape=21, size= 3, fill = 'grey')+
  geom_smooth(method='lm', colour = 'black')+
  scale_x_continuous(breaks = c(0, 10, 20, 30))+
  scale_y_continuous(breaks = c(0, 10, 20, 30))+
  theme_bw(base_size=18)+
  labs(x=expression(paste('Real ',Delta, beta[ppt]^m)),
       y=expression(paste('Sim. ',Delta, beta[ppt]^m)))+
  ggplot(delta_shift,
       aes(delta_alpha_ppt_sd_real, delta_alpha_ppt_sd_simulated)) +
  geom_point(shape=21, size = 3, fill = 'grey')+
  geom_smooth(method='lm', colour = 'black')+
  scale_x_continuous(breaks = c(0, -0.5, -1))+
  scale_y_continuous(breaks = c(0, -0.5, -1))+
  scale_alpha_continuous(range = c(0.2,1))+
  theme_bw(base_size=18)+
  labs(x=expression(paste('Real ',Delta, alpha[ppt]^sigma)),
       y=expression(paste('Sim. ',Delta, alpha[ppt]^sigma)))+
ggplot(delta_shift,
       aes(delta_alpha_ppt_m_real, abs(delta_alpha_ppt_m_simulated))) +
  geom_point(shape=21, size= 3, fill = 'grey')+
  geom_smooth(method='lm', colour = 'black')+
  scale_x_continuous(breaks = c(0, 10, 20, 30))+
  scale_y_continuous(breaks = c(0, 10, 20, 30))+
  theme_bw(base_size=18)+
  labs(x=expression(paste('Real ',Delta, alpha[ppt]^m)),
       y=expression(paste('Sim. ',Delta, alpha[ppt]^m)))&
  stat_cor(label.y.npc = 1, colour = 'red')&
  theme(legend.position = 'none',
        panel.grid = element_blank())

comparison_plot

# Preregistered Analyses --------------------------------------------------

confirm_clean <- beh_anal_adol %>%
  ungroup() %>%
  mutate(bully_total   = rowSums(select(., BVIgnored:BVSexual)),
         discrim_total = rowSums(select(., EDS1:EDS5)),
         life_total    = rowSums(select(., sum_LE_1:sum_LE_25)),
         state_para    = rowSums(select(., SP1:SP5)),
         neg_self      = rowSums(select(., IMPSelfn1:IMPSelfn6)),
         pos_self      = rowSums(select(., IMPSelfp1:IMPSelfp6)),
         neg_other     = rowSums(select(., IMPOthern1:IMPOthern6)),
         pos_other     = rowSums(select(., IMPOtherp1:IMPOtherp6)),
         ) %>%
  plyr::join(joint_parms_ex_adol %>% dplyr::select(ID, beta_par_m:beta_ref), by = 'ID') %>% 
  dplyr::select(-RT) %>%
  distinct()

t.test(confirm_clean[confirm_clean$Context=='Exc',]$life_total,
       confirm_clean[confirm_clean$Context=='Inc',]$life_total)

#Descriptives
confirm_clean %>%
  group_by(Context) %>%
  #dplyr::select(discrim_total, bully_total, life_total) %>%
  summarise(meanD = mean(discrim_total),
            meanB = mean(bully_total),
            meanL = mean(life_total))

#Hypothesis A
summary(lm(state_para ~ bully_total + attention, data = confirm_clean))
summary(lm(state_para ~ discrim_total + attention, data = confirm_clean))
summary(lm(state_para ~ life_total + attention, data = confirm_clean))

performance::check_collinearity(lm(state_para ~ bully_total + discrim_total + life_total, data = confirm_clean))
summary(lm(state_para ~ bully_total + discrim_total + life_total, data = confirm_clean))

#Hypothesis B
summary(lm(state_para ~ Context + attention, data = confirm_clean))

#Hypothesis C
confint(lm(state_para ~ bully_total:Context+attention, data = confirm_clean))
confint(lm(state_para ~ discrim_total:Context+attention, data = confirm_clean))
confint(lm(state_para ~ life_total:Context+attention, data = confirm_clean))

confirm_clean_hyp3 <- confirm_clean
y1 <- lm(state_para ~ attention, data = confirm_clean_hyp3)
confirm_clean_hyp3$state_para <- y1$residuals

ggplot(confirm_clean_hyp3 %>%
         pivot_longer(bully_total:life_total, names_to = 'Measure', values_to = 'Score') %>%
         mutate(Measure = ifelse(Measure=='bully_total', 'Bullied',
                                 ifelse(Measure=='discrim_total', 'Discrim.', 'Life Adv.'))),
       aes(Score, state_para, fill = Context))+
  geom_jitter(shape=21, alpha = 0.1)+
  geom_smooth(method='lm', colour = 'black')+
  scale_fill_manual(values=colour_group_adol)+
  scale_colour_manual(values=colour_group_adol)+
  stat_cor(method='spearman', aes(colour = Context), show.legend = F, cor.coef.name = 'rho')+
  scale_x_continuous(expand=c(0, 0.1))+
  scale_y_continuous(expand=c(0, 0.1))+
  labs(x='Score', y='State Paranoia\n(Residuals)')+
  facet_wrap(~Measure, scales = 'free')+
  theme(legend.position = 'top',
        panel.grid = element_blank(),
        strip.background = element_blank())

#Hypothesis D
summary(lm(state_para ~ neg_self + attention, data = confirm_clean))
summary(lm(state_para ~ neg_other + attention, data = confirm_clean))
summary(lm(state_para ~ neg_self:neg_other, data = confirm_clean))

ggplot(confirm_clean %>%
         pivot_longer(neg_self:pos_other, names_to = 'type', values_to = 'val'),
         aes(val, state_para, fill = type))+
  geom_jitter(shape = 21, alpha = 0.1, colour = 'black')+
  geom_smooth(method = lm, aes(colour = type))+
  scale_fill_brewer(palette = 'Paired')+
  scale_colour_brewer(palette = 'Paired')+
  stat_cor(method = 'spearman', aes(colour = type))

#Hypothesis E
summary(lm(beta_ref ~ state_para, data = confirm_clean))
summary(lm(state_para ~ alpha_ref, data = confirm_clean))

summary(lm(beta_ref ~ Context, data = confirm_clean))

ggplot(confirm_clean, aes(state_para, beta_ref))+
  geom_jitter(colour = 'black', fill = 'lightblue', shape = 21,  alpha = 0.3, width = 0.2)+
  geom_jitter(aes(y = alpha_ref),
              colour = 'black', fill = 'grey', shape = 21, alpha = 0.3, width = 0.2)+
  geom_smooth(method = 'lm', colour = 'lightblue', se = F)+
  geom_smooth(aes(y = alpha_ref),
              method = 'lm', colour = 'black', se = F)+
  coord_cartesian(ylim = c(0,15))+
  scale_x_continuous(expand=c(0,0.2))+
  scale_y_continuous(expand=c(0,0))+
  stat_cor(label.y = 10, label.x.npc = 0.4, colour = 'blue')+
  stat_cor(aes(y = alpha_ref), colour = 'black', label.y = 9, label.x.npc = 0.4)+
  labs(y = expression(paste(theta[par]^ref)), x = 'State Par.')+
  theme_bw(base_size=18)+
  theme(panel.grid = element_blank())

#Hypothesis F
hyp_f <- confirm_clean %>%
            mutate(neg_self_ord  = ifelse(neg_self  > median(neg_self), 'High', 'Low'),
                   neg_other_ord = ifelse(neg_other > median(neg_other), 'High', 'Low')
                   )

confint(lm(state_para ~ bully_total:neg_self, data = confirm_clean))
confint(lm(state_para ~ discrim_total:neg_self, data = confirm_clean))
confint(lm(state_para ~ life_total:neg_self, data = confirm_clean))

confint(lm(state_para ~ bully_total:neg_other, data = confirm_clean))
confint(lm(state_para ~ discrim_total:neg_other, data = confirm_clean))
confint(lm(state_para ~ life_total:neg_other, data = confirm_clean))

ggplot(hyp_f %>%
         pivot_longer(bully_total:life_total, names_to = 'Measure', values_to = 'Score'),
         aes(Score, state_para, fill = neg_self_ord))+
  geom_jitter(shape = 21, alpha = 0.2, colour = 'black')+
  geom_smooth(method = lm, aes(colour = neg_self_ord))+
  scale_fill_brewer(palette = 'Paired', name = 'Negative\nSelf\nBeliefs')+
  scale_colour_brewer(palette = 'Paired', name = 'Negative\nSelf\nBeliefs')+
  stat_cor(method = 'spearman', aes(colour = neg_self_ord), show.legend = F)+
  facet_wrap(~Measure, scales = 'free_x')+
  theme(legend.position = c(0.93, 0.2))
ggplot(hyp_f %>%
         pivot_longer(bully_total:life_total, names_to = 'Measure', values_to = 'Score'),
         aes(Score, state_para, fill = neg_other_ord))+
  geom_jitter(shape = 21, alpha = 0.2, colour = 'black')+
  geom_smooth(method = lm, aes(colour = neg_other_ord))+
  scale_fill_brewer(palette = 'Paired', name = 'Negative\nOther\nBeliefs')+
  scale_colour_brewer(palette = 'Paired', name = 'Negative\nOther\nBeliefs')+
  stat_cor(method = 'spearman', aes(colour = neg_other_ord), show.legend = F)+
  facet_wrap(~Measure, scales = 'free_x')+
  theme(legend.position = c(0.93, 0.2))

#Hypothesis G

int1 <- lm(state_para ~ neg_other:beta_ref, data = confirm_clean)
int2 <- lm(state_para ~ neg_other:alpha_ref, data = confirm_clean)

summary(int1)
summary(int2)

confint(int1)
confint(int2)

library(rockchalk)
ps1  <- plotSlopes(int1, plotx="neg_other", modx="beta_ref", xlab = "Neg. Other Beliefs", ylab = "Paranoia", modxVals = "std.dev")
ps2  <- plotSlopes(int2, plotx="neg_other", modx="alpha_ref", xlab = "Neg. Other Beliefs", ylab = "Paranoia", modxVals = "std.dev")

int3 <- lm(state_para ~ neg_self:beta_ref, data = confirm_clean)
int4 <- lm(state_para ~ neg_self:alpha_ref, data = confirm_clean)

summary(int3)
summary(int4)

confint(int3)
confint(int4)

library(rockchalk)
ps3  <- plotSlopes(int3, plotx="neg_self", modx="beta_ref",  xlab = "Neg. Self Beliefs", ylab = "Paranoia", modxVals = "std.dev")
ps4  <- plotSlopes(int4, plotx="neg_self", modx="alpha_ref", xlab = "Neg. Self Beliefs", ylab = "Paranoia", modxVals = "std.dev")

ps1p<- ggplot(ps1$newdata %>%
         mutate(beta_ref=ifelse(beta_ref==0.22, '-1 SD', ifelse(beta_ref==1.03, 'Mean', '+1 SD')),
                beta_ref=factor(beta_ref, levels=c('-1 SD', 'Mean', '+1 SD'))),
       aes(neg_other, fit, fill = factor(beta_ref), group=beta_ref))+
  geom_line()+
  geom_point(shape=21)+
  labs(title=expression(beta[ppt]^sigma), y='Paranoia', x='Negative Other Beliefs')+
  scale_fill_manual(values=c('#DFCAF6', '#9F61E5', '#4E178C'), name='Flexibility')+
  theme(legend.position = c(0.25, 0.75),
        panel.grid = element_blank(),
        legend.background = element_blank())
ps2p <- ggplot(ps2$newdata %>%
         mutate(alpha_v=ifelse(alpha_v==0.67, '-1 SD', ifelse(alpha_v==0.78, 'Mean', '+1 SD')),
                alpha_v=factor(alpha_v, levels=c('-1 SD', 'Mean', '+1 SD'))),
       aes(neg_other, fit, fill = factor(alpha_v), group=alpha_v))+
  geom_line()+
  geom_point(shape=21)+
  labs(title=expression(alpha[ppt]^sigma), y='Paranoia', x='Negative Other Beliefs')+
  scale_fill_manual(values=c('#DFCAF6', '#9F61E5', '#4E178C'))+
  theme(legend.position = 'none',
        panel.grid = element_blank())
ps3p <- ggplot(ps3$newdata %>%
         mutate(beta_v=ifelse(beta_v==0.88, '-1 SD', ifelse(beta_v==1.75, 'Mean', '+1 SD')),
                beta_v=factor(beta_v, levels=c('-1 SD', 'Mean', '+1 SD'))),
       aes(neg_self, fit, fill = factor(beta_v), group=beta_v))+
  geom_line()+
  geom_point(shape=21)+
  labs(title=expression(beta[ppt]^sigma), y='Paranoia', x='Negative Self Beliefs')+
  scale_fill_manual(values=c('#DFCAF6', '#9F61E5', '#4E178C'), name='Flexibility')+
  theme(legend.position = c(0.25, 0.75),
        panel.grid = element_blank(),
        legend.background = element_blank())
ps4p <- ggplot(ps4$newdata %>%
         mutate(alpha_v=ifelse(alpha_v==0.67, '-1 SD', ifelse(alpha_v==0.78, 'Mean', '+1 SD')),
                alpha_v=factor(alpha_v, levels=c('-1 SD', 'Mean', '+1 SD'))),
       aes(neg_self, fit, fill = factor(alpha_v), group=alpha_v))+
  geom_line()+
  geom_point(shape=21)+
  labs(title=expression(alpha[ppt]^sigma), y='Paranoia', x='Negative Self Beliefs')+
  scale_fill_manual(values=c('#DFCAF6', '#9F61E5', '#4E178C'))+
  theme(legend.position = 'none',
        panel.grid = element_blank())

(ps1p+ps2p)/(ps3p+ps4p) & scale_y_continuous(breaks = c(0, 5, 10, 15), limits = c(0, 15))

# Group Diffs ----

bayes_list_adol <- c('alpha',
                     'beta',
                     'alpha_v',
                     'beta_v',
                     'alpha_par',
                     'beta_par',
                     'alpha_ref',
                     'beta_ref',
                     #'alpha_par_m',
                     #'beta_par_m',
                     #'alpha_hat',
                     #'beta_hat',
                     #'disp_alpha_par',
                     #'disp_beta_par',
                     'disp_alpha_ppt',
                     'disp_beta_ppt',
                     #'disp_alpha_hat_sd',
                     #'disp_beta_hat_sd',
                     'shift_beta_p2',
                     'shift_alpha_p2'
                     )

library(BayesianFirstAid)
for(i in bayes_list_adol){

  x_parms <- joint_parms_ex_adol %>%
    distinct()

  print(paste('Now running ', i, sep = ''))

  c1 <- rstanarm::stan_glm(as.formula(paste(i, "~ group")), data = x_parms)

  c1summ  <- summary(c1)
  c1post  <- rstanarm::posterior_interval(c1, prob = 0.95)
  c1samp  <- as.data.frame(c1)

  x1_adol <- as.data.frame(t(c1summ[2,1:3])) %>%
        mutate(parameter = i,
               HDIlo = c1post[2,1],
               HDIup = c1post[2,2],
               sig = ifelse((HDIlo > 0 & HDIup > 0) | (HDIlo < 0 & HDIup < 0), 'Yes', 'No'),
               sigma = c1summ[3,1],
               sdHDIlo = c1post[3,1],
               sdHDIup = c1post[3,2],
               sigSD = ifelse((sdHDIlo > 0 & sdHDIup > 0) | (sdHDIlo < 0 & sdHDIup < 0), 'Yes', 'No'))
  x3_adol <- c1samp %>%
        rename(weight=groupADOL_in) %>%
        mutate(parameter = i)

  if(i == bayes_list_adol[1]){
    x2_adol <- x1_adol
    xmcmc_adol <- x3_adol
  } else {
    x2_adol <- rbind(x2_adol, x1_adol)
    xmcmc_adol <- rbind(xmcmc_adol, x3_adol)
  }
}

xmcmc1_adol <- plyr::join(xmcmc_adol, x2_adol %>% dplyr::select(sig, sigSD, parameter), by = 'parameter')
ggmcmc_adol <- xmcmc1_adol %>% mutate(parameter = factor(parameter, levels = bayes_list_adol))
ggx2_adol   <- x2_adol %>% mutate(var = bayes_list_adol) %>% dplyr::select(var, everything())

ggplot(ggmcmc_adol %>%
  mutate(category = ifelse(str_detect(parameter, "beta"), "beta", "alpha")) %>%
  filter(parameter %in% c('alpha', 'beta', 'alpha_v', 'beta_v')),
  aes(parameter, weight, fill = sig)) +
  geom_hline(yintercept = 0, alpha = 0.5) +
  tidybayes::stat_dist_pointinterval(aes(colour = sig)) +
  scale_colour_manual(values = c('grey', 'darkred'), name = 'HDI Crosses Zero') +
  labs(y = expression(paste(Delta, mu, ' [Inc - Exc]'))) +
  scale_x_discrete(
    labels = c(
      'beta' = expression(beta[ppt]^m),
      'beta_v' = expression(beta[ppt]^sigma),
      'alpha' = expression(alpha[ppt]^m),
      'alpha_v' = expression(alpha[ppt]^sigma)
    )
  ) +
  facet_wrap(category~parameter, scales = "free", nrow = 2)  +
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        strip.text = element_blank(),
        axis.text.x = element_text(size=24),
        axis.title.y = element_text(size=24),
        legend.background = element_blank(),
        plot.margin = margin(t = 4, r = 1, b = 1, l = 1),
        panel.grid = element_blank())

ggplot(ggmcmc_adol %>%
  mutate(category = ifelse(str_detect(parameter, "beta"), "beta", "alpha")) %>%
  filter(parameter %in% c('alpha_par', 'beta_par', 'alpha_ref', 'beta_ref')),
  aes(parameter, weight, fill = sig)) +
  geom_hline(yintercept = 0, alpha = 0.5) +
  tidybayes::stat_dist_pointinterval(aes(colour = sig)) +
  scale_colour_manual(values = c('grey', 'darkred'), name = 'HDI Crosses Zero') +
  labs(y = expression(paste(Delta, mu, ' [Inc - Exc]'))) +
  scale_x_discrete(
    labels = c(
      'beta_par' = expression(beta[par]^m),
      'beta_ref' = expression(beta[par]^ref),
      'alpha_par' = expression(alpha[par]^m),
      'alpha_ref' = expression(alpha[par]^ref)
    )
  ) +
  facet_wrap(category~parameter, scales = "free", nrow = 2)  +
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        strip.text = element_blank(),
        axis.text.x = element_text(size=24),
        axis.title.y = element_text(size=24),
        legend.background = element_blank(),
        plot.margin = margin(t = 4, r = 1, b = 1, l = 1),
        panel.grid = element_blank())

## Figure S1: Recovery and Accuracy ----------------------------------------------------------

### Check server matching ------------------------------------------------

parm_vals_adol <- adolIntent_clean %>%
  dplyr::select(ID, server_alpha_ppt, server_alpha_par) %>%
  distinct() %>%
  pivot_longer(c(server_alpha_ppt, server_alpha_par), names_to = 'player', values_to = 'alpha_val') %>%
  mutate(player = ifelse(player == 'server_alpha_ppt', 'PPT', 'PAR')) %>%
  plyr::join(.,
        adolIntent_clean %>%
  dplyr::select(ID, server_beta_ppt, server_beta_par) %>%
  distinct() %>%
  pivot_longer(c(server_beta_ppt, server_beta_par), names_to = 'player', values_to = 'beta_val') %>%
  mutate(player = ifelse(player == 'server_beta_ppt', 'PPT', 'PAR')) %>%
  dplyr::select(ID, beta_val, player),
  by = c('ID', 'player')
  )

density_plot_adol <- ggplot(parm_vals_adol %>% filter(player == 'PPT'),
                            aes(alpha_val, beta_val))+
  geom_density_2d(contour_var = 'ndensity',
                  colour = "black",
                  bins = 20)+
  labs(x = expression(alpha), y = expression(beta))+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
ggplot(parm_vals_adol %>% filter(player == 'PAR'), aes(alpha_val, beta_val))+
  geom_density_2d(contour_var = 'ndensity',
                  colour = "black",
                  bins = 20)+
  labs(x = expression(alpha), y = expression(beta))+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))&
theme(legend.position = 'none',
        text = element_text(size=12),
        panel.grid = element_blank())
density_plot_adol

### Congruency of matching between participant and partner (SimAFix/SimA) ------------------------------------------------

# Process both BPD and CON groups
ADOL_ex_results <- process_simulation_results(ADOL_ex_full_sim)
ADOL_in_results <- process_simulation_results(ADOL_in_full_sim)

FixCong_adol <- cbind(ADOL_in_results$FIX, ADOL_ex_results$FIX)
ActCong_adol <- cbind(ADOL_in_results$ACT, ADOL_ex_results$ACT)

PPT_PARcong_adol <- matrix(NA, nrow = 54, ncol = ncol(FixCong_adol))

for (i in 1:ncol(FixCong_adol)){
  match_ppt_par_adol   <- ifelse(FixCong_adol[,i] == ActCong_adol[,i], 1, 0)
  PPT_PARcong_adol[,i] <- match_ppt_par_adol
}

PPT_PARcong_adol <- as.data.frame(PPT_PARcong_adol)
names(PPT_PARcong_adol)[1:251]     <- unique(ID_in)
names(PPT_PARcong_adol)[252:502]   <- unique(ID_ex)
PPT_PARcong_adol$Trial <- 1:54

PPT_PARcong_adol <- PPT_PARcong_adol %>%
  pivot_longer(1:502, names_to = 'ID', values_to = 'match_ppt_par_adol') %>%
  group_by(ID) %>%
  mutate(matched_tot = sum(match_ppt_par_adol) / 54,
         Group = ifelse(ID %in% ID_in, 'in','ex')) %>%
  dplyr::select(ID, matched_tot, Group) %>%
  distinct()

ppt_par_congp1_adol <- ggplot(PPT_PARcong_adol, aes(matched_tot, Group, fill = Group)) +
  geom_jitter(alpha = 0.7, shape = 21, width = 0.1, height = 0.2) +
  geom_boxplot(width = 0.2, alpha = 0.7) +
  coord_cartesian(x = c(0, 1)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, '.25', '.5', '.75', '1'),
                     expand = c(0, 0)) +
  scale_fill_manual(values = colour_group_adol) +
  labs(x = 'PPT-PAR Similarity') +
  theme(legend.position = 'none',
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

ppt_par_congp1_adol

### LL ------------------------------------------------

LL_adol <- matrix(NA, nrow = 502, ncol = 6)
LL_adol[,1] <- rbind(ADOL_in_full_sim$F, ADOL_ex_full_sim$F)
LL_adol[,2] <- c(extract_and_combine(ADOL_in_full_sim$results, 'lik1', 1), extract_and_combine(ADOL_ex_full_sim$results, 'lik1', 1))
LL_adol[,3] <- c(extract_and_combine(ADOL_in_full_sim$results, 'lik2', 1), extract_and_combine(ADOL_ex_full_sim$results, 'lik2', 1))
LL_adol[,4] <- c(extract_and_combine(ADOL_in_full_sim$results, 'lik3', 1), extract_and_combine(ADOL_ex_full_sim$results, 'lik3', 1))
LL_adol[,5] <- c(unique(ID_in), unique(ID_ex))
LL_adol[,6] <- ifelse(LL_adol[,5] %in% ID_in, 'in', 'ex')

LL_adol  <- as.data.frame(LL_adol) %>% mutate(across(1:4, as.numeric))
names(LL_adol) <- c('LL', 'lik1', 'lik2', 'lik3', 'ID', 'group')

### Probabilities ------------------------------------------------

prob1_adol <- rbind(extract_and_combine(ADOL_in_full_sim$results, 'prob1', 36), extract_and_combine(ADOL_ex_full_sim$results, 'prob1', 36))
prob2_adol <- rbind(extract_and_combine(ADOL_in_full_sim$results, 'prob2', 54), extract_and_combine(ADOL_ex_full_sim$results, 'prob2', 54))
prob3_adol <- rbind(extract_and_combine(ADOL_in_full_sim$results, 'prob3', 36), extract_and_combine(ADOL_ex_full_sim$results, 'prob3', 36))

probabilities_adol <- matrix(NA, nrow = 126, ncol = nrow(prob1_adol))

for (i in 1:nrow(prob1_adol)){
  probabilities_adol[1:36, i]   <- prob1_adol[i,]
  probabilities_adol[37:90, i]  <- prob2_adol[i,]
  probabilities_adol[91:126, i] <- prob3_adol[i,]
}

probabilities_adol <- as.data.frame(probabilities_adol)
names(probabilities_adol)[1:251]     <- unique(ID_in)
names(probabilities_adol)[252:502]   <- unique(ID_ex)

probabilities_adol <- probabilities_adol %>%
  pivot_longer(1:502, names_to = 'ID', values_to = 'Probs') %>%
  group_by(ID) %>%
  mutate(Trial = 1:126,
         Group = ifelse(ID %in% ID_in, 'in','ex'),
         Phase = ifelse(Trial %in% 1:36, 1, ifelse(Trial %in% 37:90, 2, 3))) %>%
  arrange(ID, Trial)

### Check action congruency ------------------------------------------------
sim_actions_adol      <- ADOL_in_sim$data.sim[,1][[1]][[1]]
sim_actions_adol      <- as_tibble(sim_actions_adol)

for (i in 2:length(unique(beh_anal_adol$ID))){
  if(i %in% 2:251){
  sim_actions_adol        <- rbind(sim_actions_adol, as.data.frame(ADOL_in_sim$data.sim[,i][[1]][[1]]))
  }
  if(i %in% 252:502){
  sim_actions_adol        <- rbind(sim_actions_adol, as.data.frame(ADOL_ex_sim$data.sim[,i-251][[1]][[1]]))
  }
}

sim_actions_adol[,1] <- rbind(check_manipulation_full %>% filter(Context=='Inc') %>% dplyr::select(ID),
                              check_manipulation_full %>% filter(Context=='Exc') %>% dplyr::select(ID))

names(sim_actions_adol)[c(1:2,7)] <- c('ID', 'trial', 'sim_choice')
sim_actions_adol <- sim_actions_adol %>% dplyr::select(1,2,7)
sim_actions_adol <- sim_actions_adol %>% group_by(ID) %>% mutate(trial = 1:126)
check_manipulation_full <- check_manipulation_full %>% group_by(ID) %>% mutate(trial=1:126)

action_cong_adol <- plyr::join(sim_actions_adol,
                               check_manipulation_full,
                               by = c('ID', 'trial')) %>%
  mutate(cong = ifelse(sim_choice == choice, 1, 0)) %>%
  group_by(ID) %>%
  mutate(sumCong = sum(cong)/126,
         group = ifelse(ID %in% ID_in, 'in','ex')) %>%
  group_by(ID, Phase) %>%
  mutate(sumPhaseCong = ifelse(Phase %in% c(1,3), sum(cong)/36, sum(cong)/54))

### Check behavioural comparison --------------------------------------------

ADOL_in_results <- process_simulation_data(ADOL_in_sim, 'in')
ADOL_ex_results <- process_simulation_data(ADOL_ex_sim, 'ex')

# Combine prediction data
predCor_adol <- rbind(ADOL_in_results$predCor, ADOL_ex_results$predCor)

# Fit GLM models
curveSim_in <- glm(cor ~ V2, data = ADOL_in_results$predCor, family = "binomial")
curveSim_ex <- glm(cor ~ V2, data = ADOL_ex_results$predCor, family = "binomial")

# Generate prediction data
pred_time_datSim_adol <- data.frame(V2 = seq(0, 54, 0.01))
probSim_in <- predict(curveSim_in, pred_time_datSim_adol, type = "response")
probSim_ex <- predict(curveSim_ex, pred_time_datSim_adol, type = "response")

pred_time_datSim_adol <- pred_time_datSim_adol %>%
  mutate(prob_in=as.vector(probSim_in),
         prob_ex=as.vector(probSim_ex)) %>%
  pivot_longer(cols = c(prob_in, prob_ex),
               names_to = 'Group',
               values_to = 'p(Correct)') %>%
  mutate(Group = ifelse(Group=='prob_in', 'in', 'ex'))

ggplot(rbind(pred_time_datSim_adol %>% rename(trial=V2) %>% mutate(type='sim'),
             pred_time_dat %>% mutate(type='real')),
       aes(trial, `p(Correct)`, colour = Group))+
  geom_smooth(se = T, size = 1, aes(linetype=type)) +  # Smoother lines without confidence intervals
  scale_color_manual(values=colour_group_adol)+
  labs(x = 'Trial', y = 'P(Correct)') +  # Proper axis labels
  theme(
    legend.position = 'none',
    legend.background = element_rect(colour = 'black', fill = 'white'),
    legend.title = element_blank(),  # Improve legend title readability
    legend.text = element_text(size = 12)  # Improve legend text readability
  )

model_obs_adol <- rbind(
        data.frame(sumCor=c(ADOL_in_results$sumC,ADOL_ex_results$sumC)/54,
                   group=c(rep('Inc', 251), rep('Exc', 251)),
                   type='sim',
                   ID = rbind(check_manipulation_full %>%
                              filter(Context == 'Inc') %>%
                              dplyr::select(ID) %>%
                              unique(),
                              check_manipulation_full %>%
                              filter(Context == 'Exc') %>%
                              dplyr::select(ID) %>%
                              unique())),
        check_manipulation_full %>%
         filter(Phase==2) %>%
         dplyr::select(group=Context, sumCor=correctSum) %>%
         distinct() %>%
         mutate(type='real', sumCor=sumCor/54))

top_rec <- ggplot(model_obs_adol,
       aes(group, sumCor, fill = group, pattern = type, shape=type))+
  geom_jitter(shape = 21, width = 0.1, height = 0.2, size = 2, alpha = 0.1) +
  geom_boxplot(width = 0.2, alpha = 1, outlier.shape = NA) +
  geom_hline(yintercept = 0.5, size = 1.2)+
  scale_shape_manual(values=c(21, 22))+
  scale_fill_manual(values=colour_group_adol)+
  scale_y_continuous(expand=c(0,0), limits=c(0.25,1), breaks=c(0, 0.5, 1))+
  stat_compare_means(label.y = 0.1, size = 5, paired=T, label = 'p.signif')+
  coord_cartesian(y = c(0, 1))+
  labs(x = '', y = 'p(correct)')

bottom_rec <- ggplot(model_obs_adol %>%
         pivot_wider(id_cols = c(ID, group), names_from = 'type', values_from = 'sumCor'),
       aes(sim, real, fill = group))+
  geom_point(shape=21, alpha = 0.5)+
  geom_smooth(method='lm', colour = 'black')+
  stat_cor(aes(label = ..p.label.., colour = group),
           label.x = c(0.7, 0.7),
           label.y = c(0.25, 0.35))+
  stat_cor(aes(label = ..r.label.., colour = group),
           label.x = c(0.35, 0.35),
           label.y = c(0.8, 0.9))+
  scale_fill_manual(values=colour_group_adol)+
  scale_colour_manual(values=colour_group_adol)+
  scale_shape_manual(values=c(21, 22))+
  scale_y_continuous(expand=c(0,0), limits=c(0.25,1), breaks=c(0, 0.5, 1))+
  scale_x_continuous(expand=c(0,0), limits=c(0.25,1), breaks=c(0, 0.5, 1))+
  coord_cartesian(y = c(0, 1))+
  labs(x = 'Obs.', y = 'Model')+
  theme_bw()+
  theme(axis.line.x = element_line(colour = 'black'))

top_rec/bottom_rec &
  theme(plot.margin = margin(3,1,1,1),
        text = element_text(size=20),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line.y = element_line(colour = 'black'),
        legend.position = 'none',
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10)))

### Plot all ------------------------------------------------

probp1_adol <- ggplot(probabilities_adol, aes(Probs, fill = Group))+
  geom_density(alpha = 0.5)+
  scale_fill_manual(values=colour_group_adol)+
  labs(x = expression(paste('p(Choice | ', theta, ')')))+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  theme(legend.position = 'none',
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

congphase_adol <- ggplot(action_cong_adol %>%
         dplyr::select(ID, sumPhaseCong, group, Phase) %>% distinct(),
       aes(sumPhaseCong, group, fill = group))+
  geom_jitter(alpha = 0.7, shape = 21, width = 0.1, height = 0.2)+
  geom_boxplot(width = 0.2, alpha = 0.7)+
  coord_cartesian(x = c(0.5, 1))+
  scale_x_continuous(breaks = c(0.5, 0.75, 1), labels = c('.5', '.75', '1'), expand = c(0,0))+
  facet_wrap(~Phase, scales = 'free_x')+
  labs(x = 'Model Accuracy by Phase')+
  theme(legend.position = 'none')

congall_adol <- ggplot(action_cong_adol %>%
                  ungroup() %>%
                  dplyr::select(ID, sumCong, group) %>%
                  distinct(),
       aes(sumCong, group, fill = group))+
  geom_jitter(alpha = 0.7, shape = 21, width = 0.1, height = 0.2)+
  geom_boxplot(width = 0.2, alpha = 0.7)+
  geom_vline(aes(xintercept = mean(sumCong)), size = 1.2)+
  coord_cartesian(x = c(0.5, 1))+
  geom_vline(xintercept = c(0.5, 1), colour = c('grey', 'darkgreen'), size = 3)+
  coord_cartesian(xlim = c(0.5,1))+
  scale_x_continuous(breaks = c(seq(0.5, 1, 0.1)), labels = c('.5', '.6', '.7', '.8', '.9', '1'), expand = c(0,0))+
  labs(x = 'Overall Model Accuracy')+
  theme(legend.position = c(0.15, 0.65),
        legend.title = element_blank(),
        legend.background = element_rect(colour = 'black'))

LLall_adol <- ggplot(LL_adol, aes(LL, fill = group))+
  geom_density(alpha = 0.5)+
  geom_vline(aes(xintercept = mean(LL)), size = 1.2)+
  geom_vline(xintercept = c(log(0.5)*126, 0), colour = c('grey', 'darkgreen'), size = 3)+
  scale_fill_manual(values=colour_group_adol)+
  scale_x_continuous(expand = c(0,0), limits = c(-90, 1))+
  scale_y_continuous(expand = c(0,0))+
  labs(x = 'Total LL Overall')+
  theme(legend.position = 'none',
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

Accp1_adol <- ((ppt_par_congp1_adol | congall_adol)/congphase_adol) &
  scale_fill_manual(values=colour_group_adol) &
  theme(text = element_text(size = 16),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

density_plot_adol
Accp1_adol/(LLall_adol | probp1_adol)
## Figure S3: Reaction times ----------------------------------------------------------

RT_anal_adol <- adolIntent_clean %>%
  arrange(ID) %>%
  filter(ID %in% check_manipulation_tidy$ID) %>%
  mutate(group=ifelse(ID%in%IDinc$ID, 'in', 'ex'),
         distancebeta=abs(server_beta_par-server_beta_ppt),
         distancealpha=abs(server_alpha_par-server_alpha_ppt),
         total_distance = distancebeta+distancealpha,
         total_distance_q = ifelse(total_distance < 22.4, 'low', 'high'),
         total_distance_q = factor(total_distance_q, levels = c('low', 'high'), ordered = T)
         ) %>%
  distinct()

ggplot(RT_anal_adol %>%
         filter(RT<10000, Phase ==2),
       aes(trial, RT, colour= total_distance_q))+
  stat_summary(geom='line', alpha = 0.2, size=2)+
  stat_smooth(method='lm', formula = y ~ x + I(x^2))+
  scale_color_manual(name='Distance', values = c('#9A031E', 'black'))+
  coord_cartesian(ylim = c(1500, 4000))+
  stat_cor(label.y = c(4000, 3800), label.x = c(20), size=4, method = 'spearman')+
  labs(x='Trial',y='RT (ms)')+
  theme_bw(base_size=18)+
  theme(legend.position = 'none',
        panel.grid = element_blank())+

ggplot(RT_anal_adol %>%
         filter(RT<10000, Phase ==2),
       aes(trial, RT, colour= group))+
  stat_summary(geom='line', alpha = 0.2, size=2)+
  stat_smooth(method='lm', formula = y ~ x + I(x^2))+
  scale_color_manual(name='Distance', values = colour_group_adol)+
  coord_cartesian(ylim = c(1500, 4000))+
  stat_cor(label.y = c(4000, 3800), label.x = c(20), size=4, method = 'spearman')+
  labs(x='Trial',y='RT (ms)')+
  theme_bw(base_size=18)+
  theme(legend.position = 'none',
        panel.grid = element_blank())

#Phase 1
RT_model4 <- lme4::lmer(RT ~ group + (1|ID), data = RT_anal_adol %>% filter(Phase%in%c(1)))
summary(RT_model4)
confint(RT_model4)
anova(RT_model4)

#Phase 2
RT_model1a <- lmerTest::lmer(RT ~ trial +  (1|ID), data = RT_anal_adol %>% filter(Phase==2))
summary(RT_model1a)
confint(RT_model1a)

RT_model1b <- lmerTest::lmer(RT ~ total_distance +(1|ID), data = RT_anal_adol %>% filter(Phase==2))
summary(RT_model1b)
confint(RT_model1b)

RT_model1c <- lmerTest::lmer(RT ~ trial:total_distance + (1|ID), data = RT_anal_adol %>% filter(Phase==2))
summary(RT_model1c)
confint(RT_model1c)

RT_model2 <- lme4::lmer(RT ~ group + (1|ID), data = RT_anal_adol %>% filter(Phase==2))
summary(RT_model2)
confint(RT_model2)

#Phase 1-3
RT_model3a <- lme4::lmer(RT ~ Phase + group + (1|ID), data = RT_anal_adol %>% filter(Phase%in%c(1,3)))
summary(RT_model3a)
confint(RT_model3a)

RT_model3b <- lme4::lmer(RT ~ group + (1|ID), data = RT_anal_adol %>% filter(Phase%in%c(3)))
summary(RT_model3b)
confint(RT_model3b)
anova(RT_model3b)
