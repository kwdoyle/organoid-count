library(gridExtra)
library(patchwork)

# df=CFEsub_bl
# formula=as.formula("CFE ~ `Bleo day` * Genotype")

FigwStats <- function(df, indep_var="CFE", dep_vars=c("Genotype", "`Bleo day`"), 
                      term="*", jitter_color=NULL, outlier_shape=19, 
                      add_emm_table=TRUE, pval_method="emmeans", adjust_p=FALSE,
                      show_density_and_qq=FALSE, table_place="right") { # formula) {
  
  if (table_place == "right") {
    x_place <- Inf
    h_just <- 1.1
  } else if (table_place == "left") {
    x_place <- -Inf
    h_just <- -0.05
  } else {
    stop("table_place must be 'right' or 'left'")
  }
  # check for which day param used
  day_param <- grep("day", dep_vars, value=T, ignore.case=T)
  day_param2 <- gsub("`", "", day_param)
  not_day_dep_vars <- dep_vars[grep("day", dep_vars, ignore.case=T, invert=T)]
  
  formula <- as.formula(paste0(indep_var, " ~ ", paste0(dep_vars, collapse=term)))
  message("ANOVA formula: ", formula)
  
  # mod <- aov(formula, data = df)
  # Use Type III ANOVA instead
  mod.lm <- lm(formula, data = df)
  mod <- car::Anova(mod.lm, type=3)
  # summary(mod)  # difference between bleo day but NOT genotype
  if (show_density_and_qq) {
    densplt <- plot(density(residuals(mod.lm)))  # much more normal
    qqplt <- plot(mod.lm, which=2)
  } else {
    densplt <- NULL
    qqplt <- NULL
  }

  # effect sizes
  es <- effectsize::eta_squared(mod, partial=T)
  
  # if testing interactions, get all these PER day (e.g., the first dep var PER the second dep var.)
  # emm <- emmeans::emmeans(mod, ~ Genotype | `Bleo day`)
  
  # NOTE: the fact that I can add "|" between every term just so happens to work out for emmeans,
  # since the order of my inputs is what I want anyway.
  # basically, if I always pass my dependent vars as: Geno, day, then sex,
  # this will perform the | between geno and day, and then wind up doing just the ineraction
  # between day and sex anyway.
  separator <- ifelse(term == "*", "|", "+")
  # TODO if want to check comparisons BETWEEN sex FOR each Genotype and day, formula can be:
  # "~ Sex | Genotype | Day" (i.e., "~ Sex | Genotype * Day")
  # ....can this be done just by changing the order of dep_vars I pass? Yes, yes it can.
  emm_form <- as.formula(paste0("~", paste(dep_vars, collapse=separator)))
  message("emmeans formula: ", emm_form)
  
  # This is my post-hoc test
  emm <- emmeans::emmeans(mod.lm, emm_form)
  # pairs(emm.bl, adjust="none")  # similar to fisher LSD
  if (pval_method == "emmeans") {
    emm_res <- pairs(emm, adjust="tukey")  # with such a large N and strong enough effetcs, can use parametric tukey.
    emm.df <- as.data.frame(emm_res)
    
  } else if (pval_method == "t.test") {
    emm_res <- NULL
    
    # Identify grouping variables
    day_var <- day_param2
    group_var <- not_day_dep_vars[1]
    facet_var <- ifelse(length(not_day_dep_vars) > 1, not_day_dep_vars[2], NA)
    
    # Perform per-group t-tests
    if (!is.na(facet_var)) {
      emm.df <- df %>%
        group_by(.data[[day_var]], .data[[facet_var]]) %>%
        summarise(
          p.value = tryCatch(
            t.test(
              .data[[indep_var]] ~ .data[[group_var]]
            )$p.value,
            error = function(e) NA_real_
          ),
          .groups = "drop"
        )
    } else {
      emm.df <- df %>%
        group_by(.data[[day_var]]) %>%
        summarise(
          p.value = tryCatch(
            t.test(
              .data[[indep_var]] ~ .data[[group_var]]
            )$p.value,
            error = function(e) NA_real_
          ),
          .groups = "drop"
        )
    }
    
    if (adjust_p) {
      # emm.df$p.value <- p.adjust(emm.df$p.value, method="BH")
      # group by sex, if it is used
      if (!is.na(facet_var)) {
        emm.df <- emm.df %>%
          group_by(.data[[facet_var]])
      }
      emm.df <- emm.df %>%
        mutate(p.value = p.adjust(p.value, method = "BH")) %>%
        ungroup()
      
    }
    
  } else {
    stop("pval_method must be 'emmeans' or 't.test'")
  }
  
  # ensure consistent naming
  emm.df$p.value <- as.numeric(emm.df$p.value)
  # format nicely
  # emm.df$p.value <- sprintf("%.3g", emm.df$p.value)
  

  emm.df$label <- 
    ifelse(emm.df$p.value < 0.001, "***",
           ifelse(emm.df$p.value < 0.01, "**",
                  ifelse(emm.df$p.value < 0.05, "*", "ns")))
  # add position

  # emm.df$x <- emm.df$`Bleo day`
  emm.df$x <- emm.df[,day_param2]
  y_pos <- df %>%
    group_by(!!sym(day_param2)) %>%
    summarise(y = max(!!sym(indep_var), na.rm = TRUE) * 1.05)
  # add y pos
  emm.df <- emm.df %>% 
    left_join(y_pos, by=day_param2)
  
  
  
  
  
  
  # table to add on the figure
  emm_table <- emm.df %>%
    # select(!!sym(day_param2), Sex, p.value) %>%   # adjust if column names differ
    select(any_of(c(day_param2, not_day_dep_vars, "p.value"))) %>% 
    mutate(
      p.value = signif(p.value, 2)
      # p.value = sprintf("%.2g", p.value)
    )
  if (any(grepl("sex", not_day_dep_vars, ignore.case=T))) {
    emm_table <- emm_table %>% 
      arrange(!!sym(day_param2), !!sym(not_day_dep_vars[grep("sex", not_day_dep_vars, ignore.case=T)]))
    
    # split table by sex (if it's being analyzed) and add to each facet plot separately
    facet_var <- not_day_dep_vars[grep("sex", not_day_dep_vars, ignore.case=T)]
    
    emm_labels <- emm_table %>%
      group_by(!!sym(facet_var)) %>%
      summarise(
        label = paste0(
          sprintf("%-6s %12s\n", "Day", "p"),
          paste0(
            sprintf(
              "%-6s %12.2e",
              .data[[day_param2]],
              as.numeric(p.value)
            )
            ,
            collapse = "\n"
          )
        ),
        .groups = "drop"
      )
    
    
    
    coords <- df %>%
      group_by(!!sym(facet_var)) %>%
      summarise(
        x = x_place,
        y = Inf,
        .groups = "drop"
      )
    
    emm_labels <- left_join(emm_labels, coords, by = facet_var)
    
    
    
    
    
  } else {
    emm_table <- emm_table %>% 
      arrange(!!sym(day_param2))
    
    
    emm_labels <- emm_table %>%
      # group_by(!!sym(facet_var)) %>%
      summarise(
        label = paste0(
          sprintf("%-6s %12s\n", "Day", "p"),
          paste0(
            sprintf(
              "%-6s %12.2e",
              .data[[day_param2]],
              as.numeric(p.value)
            ),
            collapse = "\n"
          )
        ),
        .groups = "drop"
      )%>% 
      mutate(x=x_place, y=Inf)
    
    # 
    # coords <- df %>%
    #   group_by(!!sym(facet_var)) %>%
    #   summarise(
    #     x = Inf,
    #     y = Inf,
    #     .groups = "drop"
    #   )
    
    # emm_labels <- left_join(emm_labels, coords, by = facet_var)
  }
    
  
  
  
  
  
  
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
    jitter_layer #+
    # geom_text(data = emm.df,
    #           aes(x = x, y = y, label = label),
    #           inherit.aes = FALSE,
    #           size = 5)
  
  if (length(not_day_dep_vars) == 2) {
    mainfig <- mainfig +
      facet_wrap(as.formula(paste0("~", not_day_dep_vars[2])))
  }
  
  
  if (add_emm_table) {
    mainfig <- mainfig +
      geom_text(
        data = emm_labels,
        aes(x = x, y = y, label = label),
        inherit.aes = FALSE,
        hjust = h_just, #-0.05, #1.1,
        vjust = 1.1,
        size = 3.5,
        family = "mono"
      )
  }
  
  mainfig <- mainfig +
    theme_bw()
  
  
  return(list(fig=mainfig, density=densplt, qq=qqplt, anova=mod, emm=emm_res, es=es))
}


# This version will use table_grob to put the p value table outside of the plot instead.
# Technically can just merge it with the above version, but I made it as a test and now
# don't feel like modifying the first one.
FigwStats2 <- function(df, indep_var="CFE", dep_vars=c("Genotype", "`Bleo day`"), 
                      term="*", jitter_color=NULL, outlier_shape=19, 
                      add_emm_table=TRUE, pval_method="emmeans", adjust_p=FALSE,
                      show_density_and_qq=FALSE, table_place="right") { # formula) {
  
  if (table_place == "right") {
    x_place <- Inf
    h_just <- 1.1
  } else if (table_place == "left") {
    x_place <- -Inf
    h_just <- -0.05
  } else {
    stop("table_place must be 'right' or 'left'")
  }
  # check for which day param used
  day_param <- grep("day", dep_vars, value=T, ignore.case=T)
  day_param2 <- gsub("`", "", day_param)
  not_day_dep_vars <- dep_vars[grep("day", dep_vars, ignore.case=T, invert=T)]
  
  formula <- as.formula(paste0(indep_var, " ~ ", paste0(dep_vars, collapse=term)))
  message("ANOVA formula: ", formula)
  
  # mod <- aov(formula, data = df)
  # Use Type III ANOVA instead
  mod.lm <- lm(formula, data = df)
  mod <- car::Anova(mod.lm, type=3)
  # summary(mod)  # difference between bleo day but NOT genotype
  if (show_density_and_qq) {
    densplt <- plot(density(residuals(mod.lm)))  # much more normal
    qqplt <- plot(mod.lm, which=2)
  } else {
    densplt <- NULL
    qqplt <- NULL
  }
  
  # effect sizes
  es <- effectsize::eta_squared(mod, partial=T)
  
  # if testing interactions, get all these PER day (e.g., the first dep var PER the second dep var.)
  # emm <- emmeans::emmeans(mod, ~ Genotype | `Bleo day`)
  
  # NOTE: the fact that I can add "|" between every term just so happens to work out for emmeans,
  # since the order of my inputs is what I want anyway.
  # basically, if I always pass my dependent vars as: Geno, day, then sex,
  # this will perform the | between geno and day, and then wind up doing just the ineraction
  # between day and sex anyway.
  separator <- ifelse(term == "*", "|", "+")
  # TODO if want to check comparisons BETWEEN sex FOR each Genotype and day, formula can be:
  # "~ Sex | Genotype | Day" (i.e., "~ Sex | Genotype * Day")
  # ....can this be done just by changing the order of dep_vars I pass? Yes, yes it can.
  emm_form <- as.formula(paste0("~", paste(dep_vars, collapse=separator)))
  message("emmeans formula: ", emm_form)
  
  # This is my post-hoc test
  emm <- emmeans::emmeans(mod.lm, emm_form)
  # pairs(emm.bl, adjust="none")  # similar to fisher LSD
  if (pval_method == "emmeans") {
    emm_res <- pairs(emm, adjust="tukey")  # with such a large N and strong enough effetcs, can use parametric tukey.
    emm.df <- as.data.frame(emm_res)
    
  } else if (pval_method == "t.test") {
    emm_res <- NULL
    
    # Identify grouping variables
    day_var <- day_param2
    group_var <- not_day_dep_vars[1]
    facet_var <- ifelse(length(not_day_dep_vars) > 1, not_day_dep_vars[2], NA)
    
    # Perform per-group t-tests
    if (!is.na(facet_var)) {
      emm.df <- df %>%
        group_by(.data[[day_var]], .data[[facet_var]]) %>%
        summarise(
          p.value = tryCatch(
            t.test(
              .data[[indep_var]] ~ .data[[group_var]]
            )$p.value,
            error = function(e) NA_real_
          ),
          .groups = "drop"
        )
    } else {
      emm.df <- df %>%
        group_by(.data[[day_var]]) %>%
        summarise(
          p.value = tryCatch(
            t.test(
              .data[[indep_var]] ~ .data[[group_var]]
            )$p.value,
            error = function(e) NA_real_
          ),
          .groups = "drop"
        )
    }
    
    if (adjust_p) {
      # emm.df$p.value <- p.adjust(emm.df$p.value, method="BH")
      # group by sex, if it is used
      if (!is.na(facet_var)) {
        emm.df <- emm.df %>%
          group_by(.data[[facet_var]])
      }
      emm.df <- emm.df %>%
        mutate(p.value = p.adjust(p.value, method = "BH")) %>%
        ungroup()
      
    }
    
  } else {
    stop("pval_method must be 'emmeans' or 't.test'")
  }
  
  # ensure consistent naming
  emm.df$p.value <- as.numeric(emm.df$p.value)
  # format nicely
  # emm.df$p.value <- sprintf("%.3g", emm.df$p.value)
  
  
  emm.df$label <- 
    ifelse(emm.df$p.value < 0.001, "***",
           ifelse(emm.df$p.value < 0.01, "**",
                  ifelse(emm.df$p.value < 0.05, "*", "ns")))
  # add position
  
  # emm.df$x <- emm.df$`Bleo day`
  emm.df$x <- emm.df[,day_param2]
  y_pos <- df %>%
    group_by(!!sym(day_param2)) %>%
    summarise(y = max(!!sym(indep_var), na.rm = TRUE) * 1.05)
  # add y pos
  emm.df <- emm.df %>% 
    left_join(y_pos, by=day_param2)
  
  
  
  
  
  
  # table to add on the figure
  emm_table <- emm.df %>%
    # select(!!sym(day_param2), Sex, p.value) %>%   # adjust if column names differ
    select(any_of(c(day_param2, not_day_dep_vars, "p.value"))) %>% 
    mutate(
      p.value = signif(p.value, 2)
      # p.value = sprintf("%.2g", p.value)
    )
  
  table_grob <- tableGrob(
    emm_table,
    rows = NULL
  )
  
  
  if (any(grepl("sex", not_day_dep_vars, ignore.case=T))) {
    emm_table <- emm_table %>% 
      arrange(!!sym(day_param2), !!sym(not_day_dep_vars[grep("sex", not_day_dep_vars, ignore.case=T)]))
    
    # split table by sex (if it's being analyzed) and add to each facet plot separately
    facet_var <- not_day_dep_vars[grep("sex", not_day_dep_vars, ignore.case=T)]
    
    emm_labels <- emm_table %>%
      group_by(!!sym(facet_var)) %>%
      summarise(
        label = paste0(
          sprintf("%-6s %12s\n", "Day", "p"),
          paste0(
            sprintf(
              "%-6s %12.2e",
              .data[[day_param2]],
              as.numeric(p.value)
            )
            ,
            collapse = "\n"
          )
        ),
        .groups = "drop"
      )
    
    
    
    coords <- df %>%
      group_by(!!sym(facet_var)) %>%
      summarise(
        x = x_place,
        y = Inf,
        .groups = "drop"
      )
    
    emm_labels <- left_join(emm_labels, coords, by = facet_var)
    
    
    
    
    
  } else {
    emm_table <- emm_table %>% 
      arrange(!!sym(day_param2))
    
    
    emm_labels <- emm_table %>%
      # group_by(!!sym(facet_var)) %>%
      summarise(
        label = paste0(
          sprintf("%-6s %12s\n", "Day", "p"),
          paste0(
            sprintf(
              "%-6s %12.2e",
              .data[[day_param2]],
              as.numeric(p.value)
            ),
            collapse = "\n"
          )
        ),
        .groups = "drop"
      )%>% 
      mutate(x=x_place, y=Inf)
    
    # 
    # coords <- df %>%
    #   group_by(!!sym(facet_var)) %>%
    #   summarise(
    #     x = Inf,
    #     y = Inf,
    #     .groups = "drop"
    #   )
    
    # emm_labels <- left_join(emm_labels, coords, by = facet_var)
  }
  
  
  
  
  
  
  
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
    jitter_layer #+
  # geom_text(data = emm.df,
  #           aes(x = x, y = y, label = label),
  #           inherit.aes = FALSE,
  #           size = 5)
  
  if (length(not_day_dep_vars) == 2) {
    mainfig <- mainfig +
      facet_wrap(as.formula(paste0("~", not_day_dep_vars[2])))
  }
  
  
  if (add_emm_table) {
    mainfig <- mainfig +
      geom_text(
        data = emm_labels,
        aes(x = x, y = y, label = label),
        inherit.aes = FALSE,
        hjust = h_just, #-0.05, #1.1,
        vjust = 1.1,
        size = 3.5,
        family = "mono"
      )
  }
  
  mainfig <- mainfig +
    theme_bw() +
    table_grob +
    plot_layout(widths = c(4, 1))  # plot wider than table
  
  
  return(list(fig=mainfig, density=densplt, qq=qqplt, anova=mod, emm=emm_res, es=es))
}
