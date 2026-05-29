source("./Rscripts/FigwStats.R")
library(readxl)
library(dplyr)
library(purrr)
library(ggplot2)


# run FigwStats2 on each set of input parameters (ie, each kind of bleo day, and area or cfe, and by sex or not)
runAnalysis <- function(spec, pmth, adjp, addptb) {
  do.call(FigwStats2, c(spec, list(
    term="*",
    pval_method=pmth,
    adjust_p=adjp,
    add_emm_table=addptb
  )))
}



# ---- setup ----
pvalmethod <- "t.test"  # emmeans
adjust_p <- FALSE
add_p_table <- FALSE  # if want to use the included table outside of the plot boundary, keep this set to FALSE

orgarea <- read_xlsx('./2025 2026 organoid area CFU.xlsx')
orgarea$Day <- factor(orgarea$Day, levels=as.character(sort(unique(orgarea$Day))))
orgarea$Genotype <- factor(orgarea$Genotype, levels = c('WT','HZ'))
# add calc thing
orgarea <- orgarea %>% 
  group_by(Genotype, Sex, Day, `Bleo day`, `Experiment date`, `Mouse ID`) %>% 
  mutate(max_num_organoid = max(`No. of organoids`),
         CFE = (max_num_organoid / 15000) * 100)

#baseline only timecourse
blonly <- orgarea %>%
  filter(Day==10)
# combine days 10 and 14
blonly$`Bleo day` <- as.character(blonly$`Bleo day`)
blonly$`Bleo day`[which(blonly$`Bleo day` %in% c(10, 14))] <- "10-14"
# change names from 0 to 'baseline' and 10-14 to 'post-bleo'
blonly$`Bleo day` <- ifelse(blonly$`Bleo day` == "0", "Baseline",
                            ifelse(blonly$`Bleo day` == "10-14", "Post-Bleo", NA_character_))

blonly$log_area <- log(blonly$`Area (^inch)`)

blonly2 <- orgarea %>% 
  filter(`Bleo day` %in% c(10, 14))
blonly2$log_area <- log(blonly2$`Area (^inch)`)

# check for normalicy (p < 0.05 if NOT normally distributed)
blonly2 %>% 
  mutate(log_area = log(`Area (^inch)`)) %>% 
  group_by(Day, Genotype) %>% 
  # have to take a subset of the data points to use for testing, cause can't be over 5000????
  sample_n(size=500) %>% 
  rstatix::shapiro_test(log_area)
# none of these are normally distributed


# subset columns for CFE
CFEsub <- unique(orgarea[ , c("Day", "Genotype", "Bleo day", "Mouse ID", "max_num_organoid", "CFE")])
CFEsub_bl <- unique(blonly[ , c("Sex", "Day", "Genotype", "Bleo day", "Mouse ID", "max_num_organoid", "CFE")])
CFEsub_bl2 <- unique(blonly2[ , c("Sex", "Day", "Genotype", "Bleo day", "Mouse ID", "max_num_organoid", "CFE")])



# ---- make ALL figures at once ----

# define parameters for each
specs <- list(
  area_bleo = list(
    df=blonly, indep_var="log_area", dep_vars=c("Genotype", "`Bleo day`"),
    jitter_color="Genotype"),
  
  area_bleo_sex = list(
    df=blonly, indep_var="log_area", dep_vars=c("Genotype", "`Bleo day`","Sex"),
    jitter_color="Genotype"),
  
  area_postbleo = list(
    df=blonly2, indep_var="log_area", dep_vars=c("Genotype", "Day"),
    jitter_color="Genotype", outlier_shape=NA, table_place="left"),
  
  area_postbleo_sex = list(
    df=blonly2, indep_var="log_area", dep_vars=c("Genotype","Day","Sex"),
    jitter_color="Genotype", outlier_shape=NA, table_place="left"),
  
  cfe_bleo = list(
    df=CFEsub_bl, indep_var="CFE", dep_vars=c("Genotype","`Bleo day`"),
    jitter_color=NULL),
  
  cfe_bleo_sex = list(
    df=CFEsub_bl, indep_var="CFE", dep_vars=c("Genotype","`Bleo day`","Sex"),
    jitter_color=NULL),
  
  cfe_postbleo = list(
    df=CFEsub_bl2, indep_var="CFE", dep_vars=c("Genotype","Day"),
    jitter_color=NULL),
  
  cfe_postbleo_sex = list(
    df=CFEsub_bl2, indep_var="CFE", dep_vars=c("Genotype","Day","Sex"),
    jitter_color=NULL)
)

# make figures/stats on everything
results <- imap(specs, ~runAnalysis(.x, pmth=pvalmethod, adjp=adjust_p, addptb=add_p_table))

# and save them all
iwalk(figs, ~ggsave(
  filename = paste0("./figures/", .y, ".png"),
  plot = .x,
  width = 12.21, height = 12.375, units = "in"
))





# ---- Area Plots - OLD w/ repeated function call ----
# see if this works. Yes, this works.
# Just Bleo Day(s) v Genotype
blday_area <- FigwStats2(df=blonly, indep_var="log_area", dep_vars=c("Genotype", "`Bleo day`"), term="*", jitter_color="Genotype",
                        pval_method=pvalmethod, adjust_p=adjust_p, add_emm_table=add_p_table)
blday_area$fig
blday_area$anova
blday_area$es
blday_area$emm

# Just Bleo Day(s) v Genotype and Sex
blday_area_sex <- FigwStats2(df=blonly, indep_var="log_area", dep_vars=c("Genotype", "`Bleo day`", "Sex"), term="*", jitter_color="Genotype",
                            pval_method=pvalmethod, adjust_p=adjust_p, add_emm_table=add_p_table)
blday_area_sex$fig
blday_area_sex$anova
blday_area_sex$es
blday_area_sex$emm

# Experiment days v Genotype 
experday_area <- FigwStats2(df=blonly2, indep_var="log_area", dep_vars=c("Genotype", "Day"), term="*", jitter_color="Genotype", outlier_shape = NA,
                           pval_method=pvalmethod, adjust_p=adjust_p, table_place="left", add_emm_table=add_p_table)
experday_area$fig
experday_area$anova
experday_area$es
experday_area$emm

# Experiment days v Genotype  and Sex
experday_area_sex <- FigwStats2(df=blonly2, indep_var="log_area", dep_vars=c("Genotype", "Day", "Sex"), term="*", jitter_color="Genotype", outlier_shape = NA,
                               pval_method=pvalmethod, adjust_p=adjust_p, table_place="left", add_emm_table=add_p_table)
experday_area_sex$fig
experday_area_sex$anova
experday_area_sex$es
experday_area_sex$emm



# ---- CFE plots - OLD w/ repeated function call ----

# Bleo days v genotype
blday_cfe <- FigwStats2(df=CFEsub_bl, indep_var="CFE", dep_vars=c("Genotype", "`Bleo day`"), term="*", jitter_color=NULL,
                       pval_method=pvalmethod, adjust_p=adjust_p, add_emm_table=add_p_table)
blday_cfe$fig
blday_cfe$anova
blday_cfe$es
blday_cfe$emm

# Bleo days v genotype and sex
blday_cfe_sex <- FigwStats2(df=CFEsub_bl, indep_var="CFE", dep_vars=c("Genotype", "`Bleo day`", "Sex"), term="*", jitter_color=NULL,
                           pval_method=pvalmethod, adjust_p=adjust_p, add_emm_table=add_p_table)
blday_cfe_sex$fig
blday_cfe_sex$anova
blday_cfe_sex$es
blday_cfe_sex$emm


# Experiment day v genotype
experday_cfe <- FigwStats2(df=CFEsub_bl2, indep_var="CFE", dep_vars=c("Genotype", "Day"), term="*", jitter_color=NULL,
                          pval_method=pvalmethod, adjust_p=adjust_p, add_emm_table=add_p_table)
experday_cfe$fig
experday_cfe$anova
experday_cfe$es
experday_cfe$emm


# Experiment day v genotype and sex
# NOTE: if wish to check different post-hoc comparisons (e.g., CFE between sex FOR a given genotype on a given day),
# can just change the order of the dep_vars like so:
# experday_cfe_sex <- FigwStats(df=CFEsub_bl2, indep_var="CFE", dep_vars=c("Sex", "Genotype", "Day"), term="*", jitter_color=NULL)
experday_cfe_sex <- FigwStats2(df=CFEsub_bl2, indep_var="CFE", dep_vars=c("Genotype", "Day", "Sex"), term="*", jitter_color=NULL,
                              pval_method=pvalmethod, adjust_p=adjust_p, add_emm_table=add_p_table)
experday_cfe_sex$fig
experday_cfe_sex$anova
experday_cfe_sex$es
experday_cfe_sex$emm



# ---- save figures ----
areafigs <- list(area_bleo = blday_area$fig,
                 area_bleo_sex = blday_area_sex$fig,
                 area_postbleo = experday_area$fig,
                 area_postbleo_sex = experday_area_sex$fig)

cfefigs <- list(cfe_bleo = blday_cfe$fig,
                cfe_bleo_sex = blday_cfe_sex$fig,
                cfe_postbleo = experday_cfe$fig,
                cfe_postbleo_sex = experday_cfe_sex$fig)


lapply(names(areafigs), function(x) {
  ggsave(filename=paste0("./figures/", x, ".png"), plot=areafigs[[x]],
         width = 12.21, height = 12.375, units = "in")
         # width = 10.175, height = 10.3125, units = "in")
         # width = 8.14, height = 8.25, units = "in")  # the default sizes
})

lapply(names(cfefigs), function(x) {
  ggsave(filename=paste0("./figures/", x, ".png"), plot=cfefigs[[x]],
         width = 12.21, height = 12.375, units = "in")
         # width = 8.14, height = 8.25, units = "in")  # the default sizes
})

# ggsave(filename="./figures/area_bleo.png", plot=areafigs$area_bleo)





# ----  Area Plots (not using function) ----

# plots use:

# Bleo Day
ggplot(blonly, aes(x=`Bleo day`, y=log(`Area (^inch)`), fill=Genotype)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(aes(color=Genotype),
              position=position_jitterdodge(jitter.width=0.5),
              alpha=0.1) +
  stat_compare_means(aes(group = Genotype), method="t.test", label="p.signif")
# anova
# area.mod.bl <- aov(log(`Area (^inch)`) ~ `Bleo day` * Genotype, data=blonly)
# BIG NOTE:: do NOT want Type I ANOVA here, since data is unbalanced and clearly driven by Day
# run type II or III instead
area.mod.bl <- car::Anova(lm(log(`Area (^inch)`) ~ `Bleo day` * Genotype, data=blonly), type=3)
# summary(area.mod.bl)  # difference between bleo day, but not genotype
plot(density(residuals(area.mod.bl)))  # it's not *terribly* skewed...
plot(area.mod.bl, which=2)
# post-hoc test
# print(agricolae::LSD.test(area.mod.bl, "interaction('Bleo day', Genotype)"))
emm.bl <- emmeans::emmeans(area.mod.bl, ~ Genotype | `Bleo day`)
# pairs(emm.bl, adjust="none")  # similar to fisher LSD
pairs(emm.bl, adjust="tukey")  # with such a large N and strong enough effetcs, can use parametric tukey.

# calculate effect size (eta-squared)
# note: partial eta squared is relative to the residual
effectsize::eta_squared(area.mod.bl, partial=T)
# very small effect of genotype. but day has a moderate effect.




# Bleo Day - split by sex
ggplot(blonly, aes(x=`Bleo day`, y=log(`Area (^inch)`), fill=Genotype)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(aes(color=Genotype),
              position=position_jitterdodge(jitter.width=0.5),
              alpha=0.1) +
  facet_wrap(~ Sex) +
  stat_compare_means(aes(group = Genotype), method="t.test", label="p.signif")
# anova
# area.mod.bl.sex <- aov(log(`Area (^inch)`) ~ `Bleo day` * Genotype * Sex, data=blonly)
# area.mod.bl.sex <- aov(log(`Area (^inch)`) ~ Genotype * `Bleo day` * Sex, data=blonly)
area.mod.bl.sex.lm <- lm(log(`Area (^inch)`) ~ `Bleo day` * Genotype * Sex, data=blonly)
area.mod.bl.sex= car::Anova(area.mod.bl.sex.lm, type=3)
summary(area.mod.bl.sex)  # difference between bleo day, and sex, but not genotype
plot(density(residuals(area.mod.bl.sex)))
plot(area.mod.bl.sex, which=2)
emm.bl.sex <- emmeans::emmeans(area.mod.bl.sex.lm, ~ Genotype | `Bleo day` * Sex)
# pairs(emm.bl, adjust="none")  # similar to fisher LSD
pairs(emm.bl.sex, adjust="tukey")  # with such a large N and strong enough effetcs, can use parametric tukey.
effectsize::eta_squared(area.mod.bl.sex, partial=T)
# bleo days matters the most here too


# see if this works
blday_area_sex <- FigwStats(df=blonly, indep_var="log_area", dep_vars=c("Genotype", "`Bleo day`", "Sex"), term="*", jitter_color="Genotype")







# Experiment Day
ggplot(blonly2, aes(x=Day, y=log(`Area (^inch)`), fill=Genotype)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(aes(color=Genotype),
              position=position_jitterdodge(jitter.width=0.5),
              alpha=0.03) +
  stat_compare_means(aes(group = Genotype), method="t.test", label="p.signif")
# anova
area.mod.day <- aov(log(`Area (^inch)`) ~ Genotype * Day, data=blonly2)
summary(area.mod.day)  # difference across day and genotype
plot(density(residuals(area.mod.day)))  # kind of right-skewed..
plot(area.mod.day, which=2)
emm.day <- emmeans::emmeans(area.mod.day, ~ Genotype | Day)
# pairs(emm.bl, adjust="none")  # similar to fisher LSD
pairs(emm.day, adjust="tukey")  # with such a large N and strong enough effetcs, can use parametric tukey.
effectsize::eta_squared(area.mod.day, partial=T)
# again, day matters most here

# compare results w/ nonparametric, kruskal-wallis test
kruskal.test(log(`Area (^inch)`) ~ interaction(Day, Genotype), data = blonly2)


# Experiment Day - split by sex
emm.day.sex <- emmeans::emmeans(area.mod.day.sex, ~ Genotype | Day * Sex)
# emm.day.sex.df <- as.data.frame(pairs(emm.day.sex))
# emm.day.sex.df$label <- ifelse(emm.day.sex.df$p.value < 0.001, "***",
#                                ifelse(emm.day.sex.df$p.value < 0.01, "**",
#                                       ifelse(emm.day.sex.df$p.value < 0.05, "*", "ns")))

ggplot(blonly2, aes(x=Day, y=log(`Area (^inch)`), fill=Genotype)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(aes(color=Genotype),
              position=position_jitterdodge(jitter.width=0.5),
              alpha=0.03) +
  facet_wrap(~ Sex) +
  # stat_pvalue_manual(emm.day.sex.df, label="label", y.position=unique(blonly2$Day), xmin="Day", xmax="Day")
  stat_compare_means(aes(group = Genotype), method="t.test", label="p.signif")
# anova
area.mod.day.sex <- aov(log(`Area (^inch)`) ~ Genotype * Day * Sex, data=blonly2)
summary(area.mod.day.sex)  # difference across all day, genotype, and sex
plot(density(residuals(area.mod.day.sex)))  # kind of right-skewed..
plot(area.mod.day.sex, which=2)
emm.day.sex <- emmeans::emmeans(area.mod.day.sex, ~ Genotype | Day * Sex)
# pairs(emm.bl, adjust="none")  # similar to fisher LSD
pairs(emm.day.sex, adjust="tukey")  # with such a large N and strong enough effetcs, can use parametric tukey.
effectsize::eta_squared(area.mod.day.sex, partial=T)
# day still matters most



# ANOVA
# bleo day


# experiment days







# ggplot(orgarea, aes(x=Day, y=`Area (^inch)`, color=Genotype)) +
#   geom_boxplot() +
#   geom_jitter(position=position_jitter(width=0.2, seed=10))
# 
# 
# # The 0 and 10-14 Bleo Treated Area Plot
# ggplot(blonly, aes(x=`Bleo day`, y=log(`Area (^inch)`), fill=Genotype)) +
#   geom_boxplot(outlier.shape=NA) +
#   geom_jitter(aes(color=Genotype),
#               position=position_jitterdodge(jitter.width=0.5),
#               alpha=0.1)
#   # geom_jitter(width=0.1, alpha=0.3)
#   # geom_point(aes(group=Genotype))
#   # geom_jitter(position=position_jitter(seed=10, width=0.05), alpha=0.3)
# 
# ggplot(blonly, aes(x=`Bleo day`, y=log10(`Area (^inch)`), fill=Genotype)) +
#   geom_boxplot(outlier.shape=NA) +
#   geom_jitter(aes(color=Genotype),
#               position=position_jitterdodge(jitter.width=0.5),
#               alpha=0.1)













# ggplot(blonly, aes(x=`Bleo day`, y=log2(`Area (^inch)`), color=Genotype)) +
#   geom_boxplot() #+
# #geom_jitter()
# 
# ggplot(blonly, aes(x=`Bleo day`, y=log(`Area (^inch)`), color=Genotype)) +
#   geom_boxplot() +
#   facet_wrap(~Sex)
# #geom_jitter()


# ---- CFE Plots (not using function) ----


# check for normalicy (p < 0.05 if NOT normally distributed)
CFEsub_bl2 %>% 
  group_by(Day, Genotype) %>% 
  rstatix::shapiro_test(CFE)
# hmm... 

# CFE per bleo day
expday_cfe <- FigwStats(df=CFEsub_bl, indep_var="CFE", dep_vars=c("Genotype", "`Bleo day`"), term="*")
expdaysex_cfe <- FigwStats(df=CFEsub_bl, indep_var="CFE", dep_vars=c("Genotype", "`Bleo day`", "Sex"), term="*")

# TODO add the post-hoc to these CFE ones

# plots to use:
# Bleo Day
# anova
cfe.mod.bl <- aov(CFE ~ Genotype * `Bleo day`, data = CFEsub_bl)
summary(cfe.mod.bl)  # difference between bleo day but NOT genotype
plot(density(residuals(cfe.mod.bl)))  # much more normal
plot(cfe.mod.bl, which=2)

emm.cfe.bl <- emmeans::emmeans(cfe.mod.bl, ~ Genotype | `Bleo day`)
# pairs(emm.bl, adjust="none")  # similar to fisher LSD
pairs(emm.cfe.bl, adjust="tukey")  # with such a large N and strong enough effetcs, can use parametric tukey.
emm.cfe.bl.df <- as.data.frame(pairs(emm.cfe.bl, adjust="tukey"))
emm.cfe.bl.df$label <- 
  ifelse(emm.cfe.bl.df$p.value < 0.001, "***",
         ifelse(emm.cfe.bl.df$p.value < 0.01, "**",
                ifelse(emm.cfe.bl.df$p.value < 0.05, "*", "ns")))
# add position
emm.cfe.bl.df$x <- emm.cfe.bl.df$`Bleo day`
y_pos.cfe.bl <- CFEsub_bl %>%
  group_by(`Bleo day`) %>%
  summarise(y = max(CFE, na.rm = TRUE) * 1.05)
# add y pos
emm.cfe.bl.df <- emm.cfe.bl.df %>% 
  left_join(y_pos.cfe.bl, by="Bleo day")



ggplot(CFEsub_bl, aes(x=`Bleo day`, y=CFE, fill=Genotype)) +
  geom_boxplot() +
  geom_jitter(position=position_jitterdodge(jitter.width=0.5),
              alpha=0.9) +
  geom_text(data = emm.cfe.bl.df,
            aes(x = x, y = y, label = label),
            inherit.aes = FALSE,
            size = 5)

  # stat_compare_means(aes(group = Genotype), method="t.test", label="p.signif")

# This should return same. it does.
expday_cfe <- FigwStats(df=CFEsub_bl, indep_var="CFE", dep_vars=c("Genotype", "`Bleo day`"), term="*")
print(expday_cfe$fig)


effectsize::eta_squared(cfe.mod.bl, partial=T)
# here, bleo day has a huge effect, genotype alone is a moderate effect,
# and day*genotype has a huge effect


# Bleo Day by Sex
ggplot(CFEsub_bl, aes(x=`Bleo day`, y=CFE, fill=Genotype)) +
  geom_boxplot() +
  geom_jitter(position=position_jitterdodge(jitter.width=0.5),
              alpha=0.9) +
  facet_wrap(~ Sex)
# anova
cfe.mod.bl.sex <- aov(CFE ~ `Bleo day` + Genotype * Sex, data = CFEsub_bl)
summary(cfe.mod.bl.sex)  # difference between bleo day and MARGINALLY sex, but NOT genotype..
plot(density(residuals(cfe.mod.bl.sex)))  # tiny bit left skewed, but otherwise fine
plot(cfe.mod.bl.sex, which=2)
effectsize::eta_squared(cfe.mod.bl.sex, partial=T)
# huge effect of bleo day, medium effect of sex alone,
# but small effect of genotype alone, and even smaller genotype*sex

# have to modify other formula to not use all "|"?
expdaysex_cfe <- FigwStats(df=CFEsub_bl, indep_var="CFE", dep_vars=c("Genotype", "`Bleo day`", "Sex"), term="*")



# Experiment Day
ggplot(CFEsub_bl2, aes(x=Day, y=CFE, fill=Genotype)) +
geom_boxplot() +
  geom_jitter(position=position_jitterdodge(jitter.width=0.5),
              alpha=0.9)
# anova
cfe.mod.day <- aov(CFE ~ Day * Genotype, data = CFEsub_bl2)
summary(cfe.mod.day)  # ok, NOW a difference between day AND genotype
plot(density(residuals(cfe.mod.day)))  # little bit right skewed..
plot(cfe.mod.day, which=2)
effectsize::eta_squared(cfe.mod.day, partial=T)
# large effect of day alone and genotype alone,
# but small effect of day * genotype



# Experiment Day by Sex
ggplot(CFEsub_bl2, aes(x=Day, y=CFE, fill=Genotype)) +
  geom_boxplot() +
  geom_jitter(position=position_jitterdodge(jitter.width=0.5),
              alpha=0.9) +
  facet_wrap(~ Sex)
# anova
cfe.mod.day.sex <- aov(CFE ~ Day + Genotype * Sex, data = CFEsub_bl2)
summary(cfe.mod.day.sex)  # difference between day, genotype, and sex
plot(density(residuals(cfe.mod.day.sex)))  # mostly fine--there's a little bump in the right tail, though
plot(cfe.mod.day.sex, which=2)
effectsize::eta_squared(cfe.mod.day.sex, partial=T)
# large effect of day alone, genotype alone, and sex alone,
# but no real effect of genotype * sex



# overall
ggplot(CFEsub, aes(x=Day, y=CFE, fill=Genotype)) +
  geom_boxplot() +
  geom_jitter()






# as repeated-measures
# ............ok, it says a column that is clearly in the table doesn't exist. And I'm not going to fight with this.
# cfe.mod.day.RM <- rstatix::anova_test(data=CFEsub_bl2, dv=CFE, wid=max_num_organoid, within=Day)






# ggplot(CFEsub_bl, aes(x=`Bleo day`, y=CFE, color=Genotype)) +
#   geom_boxplot() +
#   geom_jitter() + 
#   facet_wrap(~ Sex)


# ggplot(CFEsub_bl2, aes(x=Day, y=CFE, fill=Genotype)) +

