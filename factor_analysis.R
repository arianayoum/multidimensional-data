
pacman::p_load(dplyr, lme4, afex, ez, ggplot2, tidyverse, RColorBrewer, 
               data.table, psycho, formattable, ggpubr, effsize, 
               robustHD, DescTools, retimes, car)

# Set working directory
Dir <- "setpath"
OA_Dir<-"setpath"

#############################################################
################# ALL BEHAVIOURAL MEASURES ##################
#############################################################
Study2_slopes <- Study2_slopes %>% rename(subject=i)

measures <- merge(shipley_raw, all.submeasures, by="subject") %>% 
  merge(SAT.tau, by="subject") %>% 
  merge(Study2_PC, by="subject") %>% 
  merge(Study2_slopes, by="subject")

measures.long <- measures %>% 
  gather(Measure, Score, Shipley_Raw:ind_slope, factor_key=TRUE)

# Check distribution for each measure

pdf("histograms.pdf")
histograms <- measures.long %>%
  ggplot(aes(x=Score)) +
  geom_histogram(alpha=0.6) +
  theme_classic() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)
  ) +
  xlab("") +
  ylab("Frequency") +
  facet_wrap(~Measure, scales='free')

print(histograms)
dev.off()

# Check normality
# Shapiro-Wilk normality test for OA
# with(SAT.tau, shapiro.test(tau))# p < 0.01

# Check outliers 
measures.long %>%
  filter(Measure=='Digits_forward') %>% 
  ggplot(aes(x=Measure, y=Score, label = subject)) +
  geom_boxplot(width=.5)+
  # jittered text with geom_text
  geom_text(check_overlap = TRUE,
            position=position_jitter(width=0.15)) +
  facet_grid(Measure ~ .) +
  theme_classic() +
  theme(text=element_text(size=20))

# Winsorize 
library(DescTools)
measures.w <- measures %>%
  dplyr::select(-subject) %>% 
  mutate(across(everything(), ~ Winsorize(., na.rm = TRUE))) %>% 
  rename_with(~ paste0(., "_w"))

measures.all <- cbind(measures, measures.w)

winsorized <- measures.all %>% dplyr::select(c(Shipley_Raw_w:ind_slope_w))

# Make data long for plot
winsorized.long <- winsorized %>% 
  gather(Measure, Score, Shipley_Raw_w:ind_slope_w, factor_key=TRUE)


# #############################################################
# ####################### CORR MATRIX #########################
# #############################################################
# 
library(psych)
library(corrplot)
library(psych)
library(ggplot2)
library(car)

# ## Look at all measures - used completeobs which is case-wise deletion
cormatrix <- cor(winsorized, method = 'spearman', use = 'complete.obs')

# Clustering
library(cluster)
library(maptree)
dissimilarity = 1 - cormatrix
distance = as.dist(dissimilarity)
clus <- hclust(distance)
op_k <- kgs(clus, distance, maxclus = 20)
plot (names (op_k), op_k, xlab="# clusters", ylab="Measures")
min(op_k)

# Plot

png(file="corr.png", res=300, width=6000, height=6000)
corrplot(as.matrix(cormatrix),
         #type = 'lower', diag = FALSE,
         tl.cex = 1, tl.col = "black", method = "color",
         outline = T, order = 'hclust',
         addCoef.col = "black", number.digits = 2, number.cex = 1,
         cl.pos = 'b', cl.cex = 1, addrect = 9, rect.lwd = 3,
         col = colorRampPalette(c("midnightblue", "white","darkred"))(100))
dev.off()

## Descriptive stats
describe(measures)

library(missForest)
# seed 10% missing values
measures.mis <- prodNA(measures, noNA = 0.1)
summary(measures.mis)

# impute missing values, using all parameters as default values
measures.imp <- missForest(measures.mis, ntree = 2000, mtry = 5)

# check imputed values
measures.imp$ximp

# check imputation error - NRMSE is normalized mean squared error. It represents
# error derived from imputing continuous values
measures.imp$OOBerror

measures.dat <- data.frame(measures.imp$ximp)
datamatrix <- cor(measures.dat)
corrplot(datamatrix, method='number')


# LANGUAGE MODEL
lang <- winsorized %>% 
  dplyr::select(BNT_SponSemantic_w, Shipley_Raw_w, Animals_TOTAL_w)

# EF MODEL (Processing and Control)
ef <- winsorized %>% 
  dplyr::select(Hayling2_ConvertedScore_w, Digits_total_w,
                FAS_TotalCorr_w, TMT_res_w, Stroop_res_w)

# PS MODEL
ps <- winsorized %>% 
  dplyr::select(ind_slope_w, Door_total_w, Name_total_w)

# PC MODEL
pc <- winsorized %>% 
  dplyr::select(automatic_w, VPAII_recall_score_w)

#############################################################
###################### FACTOR ANALYSIS ######################
#############################################################

# Can we conduct a factor analysis on our data?
KMO(r=cor(lang)) #should be over 0.6
cortest.bartlett(measures) #p.val should be less than 0.05
det(cor(measures.dat)) #should be a positive determinant

# Number of factors to extract
library(ggplot2)
# Define data
data<-lang

fafitfree <- fa(data, nfactors = ncol(data), rotate = "none")
n_factors <- length(fafitfree$e.values)
scree     <- data.frame(
  Factor_n =  as.factor(1:n_factors), 
  Eigenvalue = fafitfree$e.values)
ggplot(scree, aes(x = Factor_n, y = Eigenvalue, group = 1)) + 
  geom_point() + geom_line() +
  xlab("Number of factors") +
  ylab("Initial eigenvalue") +
  labs( title = "Scree Plot", 
        subtitle = "(Based on the unreduced correlation matrix)")

parallel <- fa.parallel(data)

# Language 
fa.none <- fa(r=lang, 
              nfactors = 1, 
              # covar = FALSE, SMC = TRUE,
              fm='pa', # type of factor analysis we want to use (“pa” is principal axis factoring)
              max.iter=100, # (50 is the default, but we have changed it to 100
              rotate='oblimin') # none rotation
print(fa.none)

fa.diagram(fa.none)

lang_fs <- factor.scores(lang, fa.none)
lang_fs <- lang_fs$scores
measures.all <- cbind(measures.all, lang_fs)
names(measures.all)[names(measures.all) == 'PA1'] <- 'language'

# EF 
fa.none <- fa(r=ef, 
              nfactors = 1, 
              # covar = FALSE, SMC = TRUE,
              fm='pa', # type of factor analysis we want to use (“pa” is principal axis factoring)
              max.iter=100, # (50 is the default, but we have changed it to 100
              rotate='oblimin') # none rotation
print(fa.none)

fa.diagram(fa.none)

ef_fs <- factor.scores(ef, fa.none)
ef_fs <- ef_fs$scores
measures.all <- cbind(measures.all, ef_fs)
names(measures.all)[names(measures.all) == 'PA1'] <- 'executive'

# Behavioural Pattern Separation
fa.none <- fa(r=ps, 
              nfactors = 1, 
              # covar = FALSE, SMC = TRUE,
              fm='pa', # type of factor analysis we want to use (“pa” is principal axis factoring)
              max.iter=100, # (50 is the default, but we have changed it to 100
              rotate='oblimin') # none rotation
print(fa.none)

fa.diagram(fa.none)

ps_fs <- factor.scores(ps, fa.none)
ps_fs <- ps_fs$scores
measures.all <- cbind(measures.all, ps_fs)
names(measures.all)[names(measures.all) == 'PA1'] <- 'ps'

# Behavioural Pattern Completion
fa.none <- fa(r=pc, 
              nfactors = 1, 
              # covar = FALSE, SMC = TRUE,
              fm='pa', # type of factor analysis we want to use (“pa” is principal axis factoring)
              max.iter=100, # (50 is the default, but we have changed it to 100
              rotate='oblimin') # none rotation
print(fa.none)

fa.diagram(fa.none)

pc_fs <- factor.scores(pc, fa.none)
pc_fs <- pc_fs$scores
measures.all <- cbind(measures.all, pc_fs)
names(measures.all)[names(measures.all) == 'PA1'] <- 'pc'


# Set up Interaction Model Data
interaction.m <- measures.all %>% 
  dplyr::select(subject, tau_boxcox_transformed, executive, language_log_transformed, pc, ps) %>% 
  gather(testtype, score, pc:ps, factor_key=TRUE)

# Effect code Test Type
interaction.m$testtype.ec <- ifelse(interaction.m$testtype=='pc', -1, 1)

int.m <- lmer(score ~ testtype.ec*tau_boxcox_transformed + testtype.ec*executive + testtype.ec*language_log_transformed + 
                (1|subject), data=interaction.m)
summary(int.m)

library(ggeffects)
ggpredict(int.m, terms=c('tau_boxcox_transformed', 'testtype.ec'))
mydf <- ggpredict(int.m, se=TRUE)
mydf <- ggpredict(int.m, terms = c('tau_boxcox_transformed', 'testtype.ec'))

# Invert tau on graph
mydf$tau.inv <- 1-mydf$x

ggplot(mydf, aes(x = tau.inv, y = predicted, colour = group)) +
  geom_line(size=4) +
  geom_ribbon( aes(ymin = conf.low, ymax = conf.high, fill = group, color = NULL), 
               alpha = .25) +
  scale_color_manual(values=c('mistyrose4','grey'))+
  scale_fill_manual(values=c('mistyrose2','grey')) +
  labs(
    x= "1-Tau"
  ) +
  theme_classic() +
  theme(legend.text = element_text(margin = margin(r = 10, unit = "pt"), size=20),
        axis.text=element_text(size=35),
        axis.title=element_text(size=35),
        plot.title = element_text(hjust = 0.5, size=20))

##############################################
########### HYPOTHESIS TESTING  #############
#############################################

# Check distribution of all major variables

hist(measures.all$automatic_w)
hist(measures.all$inclusion_w)
hist(measures.all$ind_slope_w)
hist(measures.all$executive)
hist(measures.all$tau_w) # positive skew
hist(measures.all$language) # negative skew

# Transform skewed variables

# Transform Tau
library(MASS)
bc <- boxcox(measures.all$tau_w ~ 1, lambda = seq(-5, 5, 0.1))
best_lambda <- bc$x[which.max(bc$y)]
measures.all$tau_boxcox_transformed <- (measures.all$tau_w^best_lambda - 1) / best_lambda

# Transform Language
constant <- max(measures.all$language, na.rm=TRUE) + 1
reflected <- constant - measures.all$language
measures.all$language_log_transformed <- log(reflected)

# Lower tau is better pattern sep (so we hypothesize a negative correlation)
ps.m <- lm(ind_slope_w ~ tau_boxcox_transformed + language_log_transformed, data=measures.all)
summary(ps.m)

ps.m <- lm(ind_slope_w ~ tau_boxcox_transformed + language_log_transformed + executive, data=measures.all)
summary(ps.m)

library(ISLR)
par(mfrow=c(2,2))
plot(ps.m)
cooksD <- cooks.distance(ps.m)
influential <- cooksD[(cooksD > (3 * mean(cooksD, na.rm = TRUE)))]
influential_indices <- which(cooksD > (3 * mean(cooksD, na.rm = TRUE)))
influential_subjects <- measures.all$subject[influential_indices]
# Remove influential subjects
measures.all_clean <- measures.all[!measures.all$subject %in% influential_subjects, ]
# Run clean dataset
ps2.m <- lm(ind_slope_w ~ tau_boxcox_transformed + language_log_transformed + executive, data=measures.all_clean)
summary(ps2.m)

library(faraway)
pairs(measures.all %>% dplyr::select(ind_slope_w, tau_boxcox_transformed, language_log_transformed, executive), col = "dodgerblue")
round(cor(measures.all %>% dplyr::select(ind_slope_w, tau_boxcox_transformed, language, executive)), 2)
# Does our executive factor score correlate with controlled processes?
cor.test(measures.all$executive, measures.all$controlled_w)
# Does our executive factor score correlate with automatic processes?
cor.test(measures.all$executive, measures.all$automatic_w)

# Higher tau is better pattern comp (so we hypothesize a positive correlation)
# First covary out language
pc.m <- lm(automatic_w ~ tau_boxcox_transformed + language_log_transformed + executive, data=measures.all)
summary(pc.m)
pairs(measures.all %>% dplyr::select(automatic_w, tau_boxcox_transformed, language_log_transformed, executive), col = "dodgerblue")
round(cor(measures.all %>% dplyr::select(automatic_w, tau_boxcox_transformed, language_log_transformed, executive)), 2)

par(mfrow=c(2,2))
plot(pc.m)
cooksD <- cooks.distance(pc.m)
pc.influential <- cooksD[(cooksD > (3 * mean(cooksD, na.rm = TRUE)))]
pc.influential_indices <- which(cooksD > (3 * mean(cooksD, na.rm = TRUE)))
pc.influential_subjects <- measures.all$subject[pc.influential_indices]
# Remove influential subjects
pc.measures.all_clean <- measures.all[!measures.all$subject %in% pc.influential_subjects, ]
# Run clean dataset
pc2.m <- lm(ind_slope_w ~ tau_boxcox_transformed + language_log_transformed + executive, data=pc.measures.all_clean)
summary(pc2.m)

# Interaction model

# Remove influential subjects
all_influential_subjects <- c(306, 322, 325, 329, 354, 365, 357)
measures.all_clean <- measures.all[!measures.all$subject %in% all_influential_subjects, ]

# scale measures
measures.all_clean <- measures.all_clean %>% 
  mutate(
    automatic_w.s=scale(automatic_w, scale=TRUE),
    ind_slope_w.s=scale(ind_slope_w, scale=TRUE)
  )

# Set up Interaction Model Data
interaction.m <- measures.all_clean %>% 
  dplyr::select(subject, tau_boxcox_transformed, executive, language_log_transformed, automatic_w.s, ind_slope_w.s) %>% 
  gather(testtype, score, automatic_w.s:ind_slope_w.s, factor_key=TRUE)

# Effect code Test Type
interaction.m$testtype.ec <- ifelse(interaction.m$testtype=='automatic_w.s', -1, 1)

int.m <- lmer(score ~ testtype.ec*tau_boxcox_transformed + testtype.ec*executive + testtype.ec*language_log_transformed + 
                (1|subject), data=interaction.m)
summary(int.m)

# Visualization
library(ggeffects)
ggpredict(int.m, terms=c('tau_boxcox_transformed', 'testtype.ec'))
mydf <- ggpredict(int.m, se=TRUE)
mydf <- ggpredict(int.m, terms = c('tau_boxcox_transformed', 'testtype.ec'))

# Invert tau on graph
mydf$tau.inv <- 1-mydf$x

int.plot <- ggplot(mydf, aes(x = tau.inv, y = predicted, colour = group)) +
  geom_line(size=4) +
  geom_ribbon( aes(ymin = conf.low, ymax = conf.high, fill = group, color = NULL), 
               alpha = .25) +
  scale_color_manual(values=c('mistyrose4','grey'))+
  scale_fill_manual(values=c('mistyrose2','grey')) +
  labs(
    x= "1-Tau"
  ) +
  theme_classic() +
  theme(legend.text = element_text(margin = margin(r = 10, unit = "pt"), size=20),
        axis.text=element_text(size=35),
        axis.title=element_text(size=35),
        plot.title = element_text(hjust = 0.5, size=20))

ggsave("Int.png", plot=int.plot)

# BPS and auto cor after covariates are factored out?
slopes.pure <- resid(ps.m)
auto.pure <- resid(pc.m)
cor.test(slopes.pure, auto.pure, method='spearman')

############ TRY WITH INCLUSION INSTEAD OF AUTO #############

# scale measures
measures.all_clean <- measures.all_clean %>% 
  mutate(
    inclusion_w.s=scale(inclusion_w, scale=TRUE),
    ind_slope_w.s=scale(ind_slope_w, scale=TRUE)
  )

# Set up Interaction Model Data
interaction.m <- measures.all_clean %>% 
  dplyr::select(subject, tau_boxcox_transformed, executive, language_log_transformed, inclusion_w.s, ind_slope_w.s) %>% 
  gather(testtype, score, inclusion_w.s:ind_slope_w.s, factor_key=TRUE)

# Effect code Test Type
interaction.m$testtype.ec <- ifelse(interaction.m$testtype=='inclusion_w.s', -1, 1)

int.m <- lmer(score ~ testtype.ec*tau_boxcox_transformed + testtype.ec*executive + testtype.ec*language_log_transformed + 
                (1|subject), data=interaction.m)
summary(int.m)

# Visualization
library(ggeffects)
ggpredict(int.m, terms=c('tau_boxcox_transformed', 'testtype.ec'))
mydf <- ggpredict(int.m, se=TRUE)
mydf <- ggpredict(int.m, terms = c('tau_boxcox_transformed', 'testtype.ec'))

# Invert tau on graph
mydf$tau.inv <- 1-mydf$x

int.plot <- ggplot(mydf, aes(x = tau.inv, y = predicted, colour = group)) +
  geom_line(size=4) +
  geom_ribbon( aes(ymin = conf.low, ymax = conf.high, fill = group, color = NULL), 
               alpha = .25) +
  scale_color_manual(values=c('mistyrose4','grey'))+
  scale_fill_manual(values=c('mistyrose2','grey')) +
  labs(
    x= "1-Tau"
  ) +
  theme_classic() +
  theme(legend.text = element_text(margin = margin(r = 10, unit = "pt"), size=20),
        axis.text=element_text(size=35),
        axis.title=element_text(size=35),
        plot.title = element_text(hjust = 0.5, size=20))

ggsave("Int.png", plot=int.plot)


#############################################################
######### HYPOTHESIS TESTING WITH -1*d'[L,F] MEASURE #########
#############################################################

# Lower tau is better pattern sep (so we hypothesize a negative correlation)
ps.m <- lm(LF_dprime_inv_w ~ tau_w + language + executive, data=measures.w)
summary(ps.m)

# Higher tau is better pattern comp (so we hypothesize a positive correlation)
# First covary out language
pc.m <- lm(automatic_w ~ tau_w + language + executive, data=measures.w)
summary(pc.m)

# Interaction model

# scale measures
measures.w <- measures.w %>% 
  mutate(
    automatic_w.s=scale(automatic_w, scale=TRUE),
    LF_dprime_inv_w.s=scale(LF_dprime_inv_w, scale=TRUE)
  )

interaction.m <- measures.w %>% 
  dplyr::select(subject, tau_w, executive, language, automatic_w.s, LF_dprime_inv_w.s) %>% 
  gather(testtype, score, automatic_w.s:LF_dprime_inv_w.s, factor_key=TRUE)

# Effect code Test Type
interaction.m$testtype.ec <- ifelse(interaction.m$testtype=='automatic_w.s', -1, 1)

int.m <- lmer(score ~ testtype.ec*tau_w + testtype.ec*executive + testtype.ec*language + 
                (1|subject), data=interaction.m)
summary(int.m)

# Visualization
library(ggeffects)
ggpredict(int.m, terms=c('tau_w', 'testtype.ec'))
mydf <- ggpredict(int.m, se=TRUE)
mydf <- ggpredict(int.m, terms = c('tau_w', 'testtype.ec'))

# Invert tau on graph
mydf$tau.inv <- 1-mydf$x

int.plot <- ggplot(mydf, aes(x = tau.inv, y = predicted, colour = group)) +
  geom_line(size=4) +
  geom_ribbon( aes(ymin = conf.low, ymax = conf.high, fill = group, color = NULL), 
               alpha = .25) +
  scale_color_manual(values=c('mistyrose4','grey'))+
  scale_fill_manual(values=c('mistyrose2','grey')) +
  labs(
    x= "1-Tau"
  ) +
  theme_classic() +
  theme(legend.text = element_text(margin = margin(r = 10, unit = "pt"), size=20),
        axis.text=element_text(size=35),
        axis.title=element_text(size=35),
        plot.title = element_text(hjust = 0.5, size=20))

ggsave("Int.png", plot=int.plot)

# BPS and auto cor after covariates are factored out?
BPS.pure <- resid(ps.m)
auto.pure <- resid(pc.m)
cor.test(BPS.pure, auto.pure, method='spearman')

# Factor analysis with new "pure" pattern separation scores
measures.w$BPS.pure <- BPS.pure

door.m <- lm(measures.w$Door_total_w ~ measures.w$executive)
door.pure <- resid(door.m)
measures.w$door.pure <- door.pure

name.m <- lm(measures.w$Name_total_w ~ measures.w$executive)
name.pure <- resid(name.m)
measures.w$name.pure <- name.pure

ps_final <- measures.w %>% select(BPS.pure, door.pure, name.pure)

fa.none <- fa(r=ps_final, 
              nfactors = 1, 
              # covar = FALSE, SMC = TRUE,
              fm='pa', # type of factor analysis we want to use (“pa” is principal axis factoring)
              max.iter=100, # (50 is the default, but we have changed it to 100
              rotate='oblimin') # none rotation
print(fa.none)

fa.diagram(fa.none)

## Visualization
ggplotRegression <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}
ggplotRegression(ps.m)

ggplotRegression(pc.m)

##########################################
############ EXPLORATORY #################
#########################################

# Fit a quadratic model
quad_model <- lm(ind_slope_w ~ poly(tau_boxcox_transformed, 2, raw=TRUE) + language_log_transformed, data = measures.all)
summary(quad_model)

library(ggplot2)

# Predict values
measures.all$fitted <- predict(quad_model, newdata = measures.all)

# Plot the data and the fitted curve
ggplot(measures.all, aes(x = tau_boxcox_transformed, y = ind_slope_w)) +
  geom_point() +
  geom_line(aes(y = fitted), color = "blue") +
  labs(title = "Quadratic Fit (Inverted U-Shape)", x = "x", y = "y") +
  theme_minimal()
