# df=CFEsub_bl
# formula=as.formula("CFE ~ `Bleo day` * Genotype")

FigwStats <- function(df, indep_var="CFE", dep_vars=c("Genotype", "`Bleo day`"), 
                      term="*", jitter_color=NULL, outlier_shape=19) { # formula) {
  formula <- as.formula(paste0(indep_var, " ~ ", paste0(dep_vars, collapse=term)))
  message("ANOVA formula: ", formula)
  
  # mod <- aov(formula, data = df)
  # Use Type III ANOVA instead
  mod.lm <- lm(formula, data = df)
  mod <- car::Anova(mod.lm, type=3)
  # summary(mod)  # difference between bleo day but NOT genotype
  densplt <- plot(density(residuals(mod.lm)))  # much more normal
  qqplt <- plot(mod.lm, which=2)
  
  # if testing interactions, get all these PER day (e.g., the first dep var PER the second dep var.)
  # emm <- emmeans::emmeans(mod, ~ Genotype | `Bleo day`)
  
  # NOTE: the fact that I can add "|" between every term just so happens to work out for emmeans,
  # since the order of my inputs is what I want anyway.
  # basically, if I always pass my dependent vars as: Geno, day, then sex,
  # this will perform the | between geno and day, and then wind up doing just the ineraction
  # between day and sex anyway.
  separator <- ifelse(term == "*", "|", "+")
  emm_form <- as.formula(paste0("~", paste(dep_vars, collapse=separator)))
  message("emmeans formula: ", emm_form)
  
  # This is my post-hoc test
  emm <- emmeans::emmeans(mod.lm, emm_form)
  # pairs(emm.bl, adjust="none")  # similar to fisher LSD
  emm_res <- pairs(emm, adjust="tukey")  # with such a large N and strong enough effetcs, can use parametric tukey.
  emm.df <- as.data.frame(emm_res)
  emm.df$label <- 
    ifelse(emm.df$p.value < 0.001, "***",
           ifelse(emm.df$p.value < 0.01, "**",
                  ifelse(emm.df$p.value < 0.05, "*", "ns")))
  # add position
  # check for which day param used
  day_param <- grep("day", dep_vars, value=T, ignore.case=T)
  day_param2 <- gsub("`", "", day_param)
  # emm.df$x <- emm.df$`Bleo day`
  emm.df$x <- emm.df[,day_param2]
  y_pos <- df %>%
    group_by(!!sym(day_param2)) %>%
    summarise(y = max(!!sym(indep_var), na.rm = TRUE) * 1.05)
  # add y pos
  emm.df <- emm.df %>% 
    left_join(y_pos, by=day_param2)
  
  not_day_dep_vars <- dep_vars[grep("day", dep_vars, ignore.case=T, invert=T)]
  
  
  if (!is.null(jitter_color)) {
    jitter_layer <- geom_jitter(
      aes(color = .data[[jitter_color]]),
      position = position_jitterdodge(jitter.width = 0.5),
      alpha = 0.1
    )
  } else {
    jitter_layer <- geom_jitter(
      position = position_jitterdodge(jitter.width = 0.5),
      alpha = 0.9
    )
  }

  
  mainfig <- ggplot(df, aes(x=.data[[day_param2]], y=.data[[indep_var]], fill=.data[[not_day_dep_vars[1]]])) +
    geom_boxplot(outlier.shape=outlier_shape) +
    jitter_layer +
    geom_text(data = emm.df,
              aes(x = x, y = y, label = label),
              inherit.aes = FALSE,
              size = 5)
  
  if (length(not_day_dep_vars) == 2) {
    mainfig <- mainfig +
      facet_wrap(as.formula(paste0("~", not_day_dep_vars[2])))
  }
  
  return(list(fig=mainfig, density=densplt, qq=qqplt, anova=mod, emm=emm_res))
}
