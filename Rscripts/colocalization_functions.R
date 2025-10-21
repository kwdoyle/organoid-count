runCountSum <- function(x, formnms, funcrun="split", dapithresh=1000) {
  id <- unique(x$id)
  type <- unique(x$type)
  
  formnms_use <- formnms[[type]]
  
  # this is set up to perform whichever calculations specified via the 'formnms' list object
  # instead of needing a separate function for each
  if (funcrun == "split") {
    out <- autoCountSumNamed(df=x, id=id, formnms=formnms_use, dapi_thresh=dapithresh, tdt_thresh=10)
  } else if (funcrun == "full") {
    out <- autoFullSumNamed(df=x, formnms=formnms_use, dapi_thresh=dapithresh, tdt_thresh=10)
  }
}


coreSummaryFunc <- function(df2, formnms, id) {
  keynmidx_fulltbl <- which(names(df2) %in% c("section", "id", "type", "dapi_colocalize"))
  
  # need to remove keynmidx after having removed the "express" names
  names_not_express <- names(df2)[-grep("express", names(df2))]
  keynmidx <- which(names_not_express %in% c("section", "id", "type", "dapi_colocalize"))
  regnms <- names_not_express[-keynmidx]
  usenms <- c(regnms, formnms)
  
  dfuse <- df2[,-keynmidx_fulltbl]
  
  if (ncol(dfuse) != length(usenms)) {
    stop("unequal number of columns to plot vs column names to use")
  }

  par(mfrow = c(3, ceiling(ncol(dfuse) / 4)))
  lapply(1:ncol(dfuse), function(x) {
    if (!all(is.na(dfuse[,x,drop=T]))) {
      hist(dfuse[,x,drop=T], main=paste0(id, ": ", usenms[x]), xlab=names(dfuse)[x])
    }
  })
  par(mfrow = c(1, 1))
  
  sumstats <- lapply(dfuse, summary)
  names(sumstats) <- usenms
  
  sumstatsdf <- lapply(names(sumstats), function(x) {
    tmp <- sumstats[[x]]
    tmpdf <- data.frame(unclass(tmp))
    tmpdf$name <- rownames(tmpdf)
    names(tmpdf)[1] <- "value"
    tmpdf2 <- tmpdf %>% pivot_wider(names_from=name, values_from=value)
    tmpdf2$metric <- x
    tmpdf3 <- tmpdf2[,c(ncol(tmpdf2), 1:(ncol(tmpdf2)-1))]
    return(tmpdf3)
  })
  
  sumstatsdf2 <- bind_rows(sumstatsdf)
  # add sample id
  sumstatsdf2$id <- id
  sumstatsdf2 <- sumstatsdf2[,c(ncol(sumstatsdf2), c(1:ncol(sumstatsdf2)-1))]
  
  return(list(df=df2, stats=sumstatsdf2))
}


parseNm <- function(formnm) {
  parts <- str_split(formnm, "/", simplify = TRUE)
  num <- paste0("n_", gsub(" ", "_", str_trim(parts[1])))
  denom <- paste0("n_", gsub(" ", "_", str_trim(parts[2])))
  list(num = num, denom = denom)
}


autoCountSumNamed <- function(df, id, formnms, dapi_thresh = 1000, tdt_thresh = 10) {
  df1 <- df %>%
    filter(n_dapi >= dapi_thresh, n_tdt > tdt_thresh)
  
  names(df1) <- gsub("670", "", names(df1))
  
  for (nm in names(formnms)) {
    formnm <- formnms[[nm]]
    cols <- parseNm(formnm)
    df1 <- df1 %>%
      mutate(!!nm := !!sym(cols$num) / !!sym(cols$denom) * 100)
  }
  
  out <- coreSummaryFunc(df2 = df1, formnms = unname(formnms), id = id)
  return(out)
}


autoFullSumNamed <- function(df, formnms, rm_names = c("section", "id", "type"), dapi_thresh = 1000, tdt_thresh = 10) {
  # Remove low DAPI count and bad sections
  df1 <- df %>%
    filter(n_dapi >= dapi_thresh, n_tdt > tdt_thresh)
  
  names(df1) <- gsub("670", "", names(df1))
  # Drop rm_names columns, sum the rest
  col_keep <- setdiff(names(df1), rm_names)
  tmp <- as.data.frame(lapply(df1[, col_keep, drop=FALSE], sum, na.rm = TRUE))
  
  # Convert to tibble for mutate
  tmp2 <- as_tibble(tmp)
  
  # Add ratio columns from formnms (named vector/list)
  for (nm in names(formnms)) {
    formnm <- formnms[[nm]]
    cols <- parseNm(formnm)
    tmp2 <- tmp2 %>%
      mutate(!!nm := !!sym(cols$num) / !!sym(cols$denom) * 100)
  }
  
  return(tmp2)
}


dayGTcompare <- function(day, d, yvar, xvar="Day") {
  GTvals <- unique(d$GT)
  d$GT_id <- paste(d$GT, d$id, sep=".")
  d$GT_id <- as.factor(d$GT_id)
  if (length(GTvals) != 2) {
    warning("2 unique GTs not found when running t test")
  }
  
  tstdat <- d[d[,xvar,drop=T] == day & d$Sex != "Combined", ]
  
  xdat <- tstdat[tstdat$GT == GTvals[1], yvar, drop=T]
  ydat <- tstdat[tstdat$GT == GTvals[2], yvar, drop=T]
  
  # check number of unique IDs
  u.id.chk <- tstdat %>% ungroup() %>% dplyr::select(c(id, Sex, GT, !!sym(xvar))) %>% distinct()
  
  # have 2 slices per animal in this analysis
  if (any(duplicated(u.id.chk$id))) {
    stop("Duplicate data found for wilcox test! On day ", day, " for param ", yvar )
  }
  
  # test, check plots for tha current day
  d$GT <- factor(d$GT, levels=c("Hz","WT")) # to set color correctly, just swap factor order
  tstplt1 <- ggplot(d[d[,xvar,drop=T] == day & d$Sex != "Combined", ], aes(x=.data[[xvar]], y=.data[[yvar]], fill=GT, group=GT_id)) + geom_boxplot()
  tstplt2 <- ggplot(d[d[,xvar,drop=T] == day & d$Sex != "Combined", ], aes(x=.data[[xvar]], y=.data[[yvar]], fill=GT)) + geom_boxplot()
  tstplts <- list(plt1=tstplt1, plt2=tstplt2)
  
  if ((all(xdat == 0) | all(is.na(xdat))) | (all(ydat == 0) | all(is.na(ydat)))) {
    wout <- NULL
  } else {
    wout <- wilcox.test(x=xdat, y=ydat, alternative="two.sided", conf.int=T)
    wout[[as.character(GTvals[1])]] <- median(xdat)
    wout[[as.character(GTvals[2])]] <- median(ydat)
    
  }
  
  return(list(wilcox=wout, plots=tstplts))
}


plotInLoop_multi <- function(df, y_var, x_var="Day", col_var=NULL, calcKey, rm_outliers=TRUE, show_p_correct=TRUE, plotfull=FALSE, 
                             plotallagg=FALSE, pvalsize=5, labelsize=12, ptsize=0.9, make_manual_plot=FALSE, ids_rm=NULL, makemainfig=FALSE) {
  if (show_p_correct) {
    pcoluse <- "p_adj"
  } else {
    pcoluse <- "p_value"
  }
  # reorder GT
  df$GT <- factor(df$GT, levels=c("WT", "Hz"))
  
  if (all(range(na.omit(df[,y_var])) == 0)) {
    sub_dist <- 0.01
  } else if (diff(range(na.omit(df[,y_var]))) <= 5) {
    sub_dist <- 0.1
  } else {
    sub_dist <- 5
  }
  
  df <- df %>%
    filter(!is.na(.data[[y_var]])) %>% 
  group_by(as.factor(Day))
  
  if (!plotfull & rm_outliers) {
    df <- df %>% 
    mutate(is_outlier = !between(.data[[y_var]], 
                                 quantile(.data[[y_var]], 0.25) - 1.5 * IQR(.data[[y_var]]), 
                                 quantile(.data[[y_var]], 0.75) + 1.5 * IQR(.data[[y_var]]))) %>%
    filter(!is_outlier) %>%
      distinct()
  }
  
  if (!is.null(ids_rm)) {
    df <- df %>% 
      filter(!id %in% ids_rm)
  }
  
  df$Day <- factor(df$Day, levels=sort(unique(df$Day)))
  # perform tests
  testdays <- sort(unique(df[,x_var,drop=T]))
  GTdaytests <- lapply(testdays, dayGTcompare, d=df, yvar=y_var, xvar=x_var)
  names(GTdaytests) <- testdays
  
  # extract to data frame
  wdf.1 <- t(sapply(GTdaytests, function(x) {
    wil <- x$wilcox
    
    if (is.null(wil)) {
      outt <- c(p_value = NA, WT_median=NA, Hz_median=NA, diff_in_location=NA, lower_bound = NA, upper_bound = NA)
      return(outt)
    } else {
      p_value <- wil$p.value
      conf_int <- wil$conf.int
      est <- unname(wil$estimate)
      wtmed <- wil$WT
      hzmed <- wil$Hz
      c(p_value = p_value, WT_median=wtmed, Hz_median=hzmed, diff_in_location=est, lower_bound = conf_int[1], upper_bound = conf_int[2])
    }
    
  }))
  wdf <- as.data.frame(wdf.1)
  wdf$Day <- rownames(wdf)
  wdf$param <- y_var
  rownames(wdf) <- NULL
  # multiple comparisons correct here instead
  wdf$p_adj <- p.adjust(wdf$p_value, method="bonf")
  
  wiltestplts <- lapply(GTdaytests, function(x) {
    pts <- x$plots
    return(pts)
  })
  names(wiltestplts) <- names(GTdaytests)
  
  # summary stats per 'box' in the boxplots
  df_smry <- df %>%
    filter(Sex != "Combined") %>% 
    group_by(across(any_of(c('Virus', 'Day', 'id', 'Sex', 'GT')))) %>%
    summarise(param=y_var,
              mean = mean(.data[[y_var]]),
              sd = sd(.data[[y_var]]),
              median = median(.data[[y_var]]),
              qt_25 = quantile(.data[[y_var]])[2],
              qt_75 = quantile(.data[[y_var]])[4])
  
  wilcoxreturn <- list(wilcoxdf=wdf, wilcoxplots=wiltestplts, dfsmry=df_smry)
  
  if (make_manual_plot) {
    if (!("M" %in% df$Sex & "F" %in% df$Sex)) {
      df <- df[df$Sex != "Combined", ]
    }
    
    colors2 <- c("WT" = "#00BFC4", "Hz" = "#F8766D")
    colors3 <- c("No AAV" = "grey", "WT-AAV" = "green", "Mut-AAV" = "red")
    df$GT_id <- paste(df$GT, df$id, sep=".")
    ordertbl2 = unique(df[,c("GT", "sex_orig", "id", "GT_id")]) %>%
      arrange(GT, sex_orig, id)
    factor_order2 <- ordertbl2$GT_id
    
    df$GT_id <- factor(df$GT_id, levels=factor_order2)
    # new way to make same plot as above
    df$slice_id <- paste(df$slice, df$id, sep=".")
    
    if ("Virus" %in% names(df)) {
      ordertblnew = unique(df[,c("slice",  "Virus", "sex_orig",  "GT", "id", "slice_id")]) %>% 
        arrange(Virus,  sex_orig, id, GT, id)
    } else {
      ordertblnew = unique(df[,c("GT", "slice", "sex_orig", "id", "slice_id")]) %>% 
        arrange(GT, sex_orig, id)
    }
    factor_order3 <- ordertblnew$slice_id
    
    df$slice_id <- factor(df$slice_id, levels=factor_order3)
    
    
    idlabdat <- df %>% 
      group_by(id, Day, GT, slice, Sex, sex_orig, slice_id) %>%
      group_by(across(c(id, !!sym(x_var),!!sym(col_var), GT, slice, Sex, sex_orig, slice_id))) %>% 
      summarise(ypos = (min(!!sym(y_var)) * 0.9) - sub_dist) %>%
      ungroup() %>% 
      dplyr::select(-sex_orig)
    
    
    # fill using Virus instead
    plt1 <- ggplot(df, aes(x=.data[[x_var]], y=.data[[y_var]], fill=Virus, group=slice_id, color=Virus))
      
    plt1 <- plt1 +
      geom_boxplot(outlier.shape=NA, color="black") +
      facet_wrap(~ Sex, ncol=1, labeller = as_labeller(c(M="M only", F="F only", Combined="M+F")), scales="free_x") +
      scale_fill_manual(values = colors3) +
      # this is to color the IDs with the same virus color 
      scale_color_manual(values = colors3) +
      scale_x_discrete(limits = levels(df[,x_var, drop=T])) +
      labs(x="Days", y="Percent", title=calcKey[[y_var]]) +
      theme_bw() +
      theme(axis.text=element_text(size=labelsize),
            axis.text.x=element_text(size=labelsize),
            axis.text.y=element_text(size=labelsize),
            axis.title=element_text(size=labelsize),
            legend.key.size=unit(1.5, 'cm'),
            legend.title=element_text(size=14),
            legend.text=element_text(size=12),
            plot.title=element_text(size=labelsize + 5),
            strip.text=element_text(size=labelsize)) +
      geom_text(data=idlabdat, aes(label=as.character(id), y=ypos),
                position=position_dodge(width=0.75),
                show.legend=F,
                size=5, angle=90, fontface="bold")
    
    print(plt1)
    
    pltlst <- list(plt1)
  } else if ((length(unique(df$slice)) == 1 & 'section' %in% names(df)) & !makemainfig) {
    if (!("M" %in% df$Sex & "F" %in% df$Sex)) {
      df <- df[df$Sex != "Combined", ]
    }
    colors2 <- c("WT" = "#00BFC4", "Hz" = "#F8766D")
    df$GT_id <- paste(df$GT, df$id, sep=".")
    ordertbl2 = unique(df[,c("GT", "sex_orig", "id", "GT_id")]) %>%
      arrange(GT, sex_orig, id)
    factor_order2 <- ordertbl2$GT_id

    df$GT_id <- factor(df$GT_id, levels=factor_order2)
    # new way to make same plot as above
    df$slice_id <- paste(df$slice, df$id, sep=".")
    
    if ("Virus" %in% names(df)) {
      ordertblnew = unique(df[,c("slice", "GT", "Virus", "sex_orig", "id", "slice_id")]) %>% 
        arrange(GT, Virus, sex_orig, id)
    } else {
      ordertblnew = unique(df[,c("GT", "slice", "sex_orig", "id", "slice_id")]) %>% 
        arrange(GT, sex_orig, id)
    }
    factor_order3 <- ordertblnew$slice_id
    
    df$slice_id <- factor(df$slice_id, levels=factor_order3)


    idlabdat <- df %>% 
      group_by(id, Day, GT, slice, Sex, sex_orig, slice_id) %>%
      group_by(across(c(id, !!sym(x_var), GT, slice, Sex, sex_orig, slice_id))) %>% 
      summarise(ypos = (min(!!sym(y_var)) * 0.9) - sub_dist) %>%  #0.05)
      ungroup() %>% 
      dplyr::select(-sex_orig)
    
     # M and F separate via facet_wrap
    if (!is.null(col_var)) {
       plt1 <- ggplot(df, aes(x=.data[[x_var]], y=.data[[y_var]], fill=GT, group=slice_id, color=.data[[col_var]]))
      
    } else if (is.null(col_var)) {
       plt1 <- ggplot(df, aes(x=.data[[x_var]], y=.data[[y_var]], fill=GT, group=slice_id))
      
    } else {
      stop("col_var not set correctly")
    }
    plt1 <- plt1 +
      geom_boxplot(outlier.shape=NA) +
       facet_wrap(~ Sex, ncol=1, labeller = as_labeller(c(M="M only", F="F only", Combined="M+F")), scales="free_x") +
      scale_fill_manual(values = colors2) +
         scale_x_discrete(limits = levels(df[,x_var, drop=T])) +
      geom_quasirandom(dodge.width=0.75, width=0.05, varwidth=T, groupOnX=T) +
      labs(x="Days", y="Percent", title=calcKey[[y_var]]) +
      theme_bw() +
      stat_compare_means(label="p.format", method="t.test", size=pvalsize) +
       theme(axis.text=element_text(size=labelsize),
             axis.text.x=element_text(size=labelsize),
             axis.text.y=element_text(size=labelsize),
             axis.title=element_text(size=labelsize),
             legend.key.size=unit(1.5, 'cm'),
             legend.title=element_text(size=14),
             legend.text=element_text(size=12),
             plot.title=element_text(size=labelsize + 5),
             strip.text=element_text(size=labelsize)) +
       geom_text(data=idlabdat, aes(label=as.character(id), y=ypos),
                 position=position_dodge(width=0.75),
                 show.legend=F,
                 size=5, angle=90, fontface="bold")
    
     pltlst <- list(plt1)
    
  } else if ((length(unique(df$slice)) == 2 & 'section' %in% names(df)) | makemainfig) {
    
    colors2 <- c("WT" = "#00BFC4", "Hz" = "#F8766D") 
    # make new variable for the 'slices combined' plot
    df$GT_id <- paste(df$GT, df$id, sep=".")
    ordertbl2 = unique(df[,c("GT", "sex_orig", "id", "GT_id")]) %>%
      arrange(GT, sex_orig, id)
    factor_order2 <- ordertbl2$GT_id

    df$GT_id <- factor(df$GT_id, levels=factor_order2)
    # new way to make same plot as above
    df$slice_id <- paste(df$slice, df$id, sep=".")
    
    ordertblnew = unique(df[,c("GT", "slice", "sex_orig", "id", "slice_id")]) %>% 
      arrange(GT, sex_orig, id)
    factor_order3 <- ordertblnew$slice_id
    
    df$slice_id <- factor(df$slice_id, levels=factor_order3)


    idlabdat <- df %>% 
      group_by(across(c(id, !!sym(x_var), GT, slice, Sex, sex_orig, slice_id))) %>% 
      summarise(ypos = (min(!!sym(y_var)) * 0.9) - sub_dist) %>%
      ungroup() %>% 
      dplyr::select(-sex_orig)
    
   plt1 <- ggplot(df, aes(x=.data[[x_var]], y=.data[[y_var]], fill=GT, group=slice_id)) +
        geom_boxplot(outlier.shape=NA) +
         facet_wrap(~ Sex, ncol=1, labeller = as_labeller(c(M="M only", F="F only", Combined="M+F")), scales="free_x") +
        scale_fill_manual(values = colors2) +
         scale_x_discrete(limits = levels(df[,x_var, drop=T])) +
        geom_quasirandom(dodge.width=0.75, width=0.05, varwidth=T, groupOnX=T) + #dodge.width=0.95 #width=0.2
        labs(x="Days", y="Percent", fill="GT", title=paste(calcKey[[y_var]], "-- separate slices")) +
        theme_bw() +
         theme(axis.text=element_text(size=labelsize),
               axis.text.x=element_text(size=labelsize),
               axis.text.y=element_text(size=labelsize),
               axis.title=element_text(size=labelsize),
               legend.key.size=unit(1.5, 'cm'),
               legend.title=element_text(size=14),
               legend.text=element_text(size=12),
               plot.title=element_text(size=17),
               strip.text=element_text(size=labelsize)) +
         geom_text(data=idlabdat, aes(label=as.character(id), y=ypos),
                                position=position_dodge(width=0.75),
                                show.legend=F,
                   size=5, angle=90, fontface="bold")
    
    
    # label data
    idlabdat2 <- df %>% 
      group_by(across(c(id, !!sym(x_var), GT, Sex, sex_orig, GT_id))) %>% 
      summarise(ypos = (min(!!sym(y_var)) * 0.9) - sub_dist) %>%
      ungroup() %>% 
      dplyr::select(-sex_orig)
    if (plotallagg) {
      idlabdat2_use <- subset(idlabdat2, Sex != "Combined")
    } else {
      idlabdat2_use <- idlabdat2
    }
    
    # label with GTday wilcox test p values too for the 'combined' plot
    gttestlabdat <- df %>% 
      ungroup() %>% 
      left_join(wdf, by="Day") %>% 
      mutate(ypos = max(!!sym(y_var)) + 5) %>% 
      dplyr::select(c(!!sym(x_var), GT, GT_id, Sex, p_value, p_adj, ypos)) %>% 
      distinct() %>% 
      filter(Sex == "Combined")
        
    if (plotallagg) {
      plt2 <- ggplot(df, aes(x=.data[[x_var]], y=.data[[y_var]], fill=GT))
      plttitle <- "slices and mice combined"
    } else {
      plt2 <- ggplot(df, aes(x=.data[[x_var]], y=.data[[y_var]], fill=GT, group=GT_id))
      plttitle <- "slices combined"
      
    }
    
    plt2 <- plt2 +
      geom_boxplot(outlier.shape=NA, position=position_dodge(width=0.75)) +
       facet_wrap(~ Sex, ncol=1, labeller = as_labeller(c(M="M only", F="F only", Combined="M+F")), scales="free_x") +
      scale_fill_manual(values = colors2) +
     scale_color_manual(values = colors2) +
     scale_x_discrete(limits = levels(df[,x_var, drop=T])) +
      geom_quasirandom(dodge.width=0.75, width=0.05, varwidth=T, groupOnX=T) + 
      labs(x="Days", y="Percent", title=paste(calcKey[[y_var]], "--", plttitle)) +
      theme_bw() +
       theme(axis.text=element_text(size=labelsize),
             axis.text.x=element_text(size=labelsize),
             axis.text.y=element_text(size=labelsize),
             axis.title=element_text(size=labelsize),
             legend.key.size=unit(1.5, 'cm'),
             legend.title=element_text(size=14),
             legend.text=element_text(size=12),
             plot.title=element_text(size=17),
             strip.text=element_text(size=labelsize)) +
       geom_text(data=idlabdat2_use, aes(label=as.character(id), y=ypos),
                              position=position_dodge(width=0.75),
                              show.legend=F,
                 size=5, angle=90, fontface="bold") +
     geom_text(data=gttestlabdat, aes(label=paste("p =", as.character(format(get(pcoluse), scientific=T, digits=3))), y=ypos),
               size=5)
     
    # add comparisons separately -- ignore the 'combined' group
    plt2 <- plt2 +
        stat_compare_means(label = "p.format", method = "wilcox.test", size = pvalsize,
                     data = subset(df, Sex != "Combined"))
    
    
    plt3 <- ggplot(filter(df, Sex != "Combined"), aes(x=.data[[x_var]], y = .data[[y_var]], color = GT)) +
      geom_jitter(position=position_jitter(seed=42, width=0.2), alpha=0.6, size=ptsize) +
      scale_color_manual(values = colors2) +
      stat_summary(fun = "median", fun.min = "median", fun.max= "median", size= 0.5, width = 0.3, geom = "crossbar")+
      stat_summary(fun = 'median', geom = 'line', aes(group = GT, color = GT), size = 1.5) +
      theme_bw() +
      theme(axis.text=element_text(size=labelsize),
            axis.text.x=element_text(size=labelsize),
            axis.text.y=element_text(size=labelsize),
            axis.title=element_text(size=labelsize),
            legend.key.size=unit(1.5, 'cm'),
            legend.title=element_text(size=14),
            legend.text=element_text(size=12),
            plot.title=element_text(size=labelsize + 5),
            strip.text=element_text(size=labelsize)) +
      labs(x="Days after Bleomycin", title=paste(calcKey[[y_var]])) +
      scale_y_continuous(labels = function(x) paste0(x, "%")) +
      theme(axis.title.y = element_blank()) +
      geom_text(data=gttestlabdat, aes(label=paste("p =", as.character(format(get(pcoluse), scientific=T, digits=3))), y=ypos),
                size=5, color='black')
    
    
    pltlst <- list(plt1, plt2, plt3)
    
    
  } else if (length(unique(df$slice)) == 3) {
    # plot 1 with all slices separate
    colors1 <- c("WT.1" = "#00BFC4", "WT.2" = "#7CAE00", "WT.3" = "#CD9600",
              "Hz.1" = "#F8766D", "Hz.2" = "#E76BF3", "Hz.3" = "#00A9FF")
    
    plt1 <- ggplot(df, aes(x=as.factor(Day), y=.data[[y_var]], fill=interaction(GT, slice))) +
      geom_boxplot(outlier.shape=NA) +
      scale_fill_manual(values = colors1) +
      geom_jitter(position=ggbeeswarm::position_beeswarm(cex=1.5, dodge.width=0.75)) +
      labs(x="Days", y="Percent", fill="GT/slice", title=calcKey[[y_var]]) +
      theme_bw() +
       theme(axis.text=element_text(size=labelsize))
    
    # plot 2 with slice 2 and 3 together
    df$slice_combine <- ifelse(df$slice %in% c(2,3), "2&3", df$slice)
    df$slice_combine <- factor(df$slice_combine, levels=c("1", "2&3"))
    colors2 <- c("WT.1" = "#00BFC4", "WT.2&3" = "#CCCCCC",
              "Hz.1" = "#F8766D", "Hz.2&3" = "#6B6B6B") 
    
    plt2 <- ggplot(df, aes(x=as.factor(Day), y=.data[[y_var]], fill=interaction(GT, slice_combine))) +
      geom_boxplot(outlier.shape=NA) +
      scale_fill_manual(values = colors2) +
      geom_jitter(position=ggbeeswarm::position_beeswarm(cex=1.5, dodge.width=0.75)) +
      labs(x="Days", y="Percent", fill="GT/slice", title=calcKey[[y_var]]) +
      theme_bw() +
       theme(axis.text=element_text(size=labelsize))#+
    
    # plot 3 with all slices combined
    colors3 <- c("WT"="#DCDCDC", "Hz"="#778899")
    plt3 <- ggplot(df, aes(x=as.factor(Day), y=.data[[y_var]], fill=GT)) +
      geom_boxplot(outlier.shape=NA) +
      scale_fill_manual(values = colors3) +
      geom_jitter(position=ggbeeswarm::position_beeswarm(cex=1.5, dodge.width=0.75)) +
      labs(x="Days", y="Percent", title=calcKey[[y_var]]) +
      theme_bw() + 
      stat_compare_means(label="p.format", method="t.test", size=pvalsize) +
       theme(axis.text=element_text(size=labelsize))
    
    pltlst <- list(plt1, plt2, plt3)
    
    
  } else if (!'section' %in% names(df)) {
    # make new variable for the 'slices combined' plot
    df$GT_id <- paste(df$GT, df$id, sep=".")
    ordertbl2 = unique(df[,c("GT", "sex_orig", "id", "GT_id")]) %>% 
      arrange(GT, sex_orig, id)
    factor_order2 <- ordertbl2$GT_id
    
    df$GT_id <- factor(df$GT_id, levels=factor_order2)
    
    colors2 <- c("WT" = "#00BFC4", "Hz" = "#F8766D")
  
    plt1 <- NULL
    
    # label data
    idlabdat2 <- df %>% 
      group_by(id, Day, GT, Sex, GT_id) %>% 
      summarise(ypos = (min(!!sym(y_var)) * 0.9) - sub_dist)
    
    # label with GTday wilcox test p values too for the 'combined' plot
    gttestlabdat <- df %>% 
      ungroup() %>% 
      left_join(wdf, by="Day") %>% 
      mutate(ypos = max(!!sym(y_var)) + 5) %>% 
      dplyr::select(Day, GT, GT_id, Sex, p_value, p_adj, ypos) %>% 
      distinct() %>% 
      filter(Sex == "Combined")
    
   plt2 <- ggplot(df, aes(x=as.factor(Day), y=.data[[y_var]], fill=GT)) +
      geom_boxplot(outlier.shape=NA, position=position_dodge(width=0.75)) +
       facet_wrap(~ Sex, ncol=1, labeller = as_labeller(c(M="M only", F="F only", Combined="M+F"))) +
      scale_fill_manual(values = colors2) +
     scale_color_manual(values = colors2) +
       scale_x_discrete(limits = levels(df$Day)) +
      geom_quasirandom(aes(color=NULL), dodge.width=0.75, width=0.02, varwidth=T, groupOnX=T, size=3) +
     # add trendline
     geom_smooth(aes(group=GT, color=GT), method="loess", se=TRUE, position=position_dodge(width=0.75)) +
      labs(x="Days", y="Percent", title=paste(calcKey[[y_var]], "-- aggregate calculation over sums across both replicate slices")) +
      theme_bw() +
       theme(axis.text=element_text(size=labelsize),
             axis.text.x=element_text(size=labelsize),
             axis.text.y=element_text(size=labelsize),
             axis.title=element_text(size=labelsize),
             legend.key.size=unit(1.5, 'cm'),
             legend.title=element_text(size=14),
             legend.text=element_text(size=12),
             plot.title=element_text(size=17),
             strip.text=element_text(size=labelsize)) +
 
     geom_text(data=gttestlabdat, aes(label=paste("p =", as.character(format(.data[["p_value"]], scientific=T, digits=3))), 
                                      y=ypos), 
               size=5, color="#000000")
     
    # add comparisons separately -- ignore the 'combined' group
    plt2 <- plt2 +
      stat_compare_means(label="p.format", method="wilcox.test", size=pvalsize,
                         data=subset(df, Sex != "Combined"),
                         aes(x=as.factor(Day), y=.data[[y_var]], group=GT),
                         show.legend=F)
    
    
    
    pltlst <- list(plt1, plt2)
    
    # also return the wilcox test
    return(list(plots=pltlst, dayGTtest=wilcoxreturn))
  } else {
    stop("Slice length is ", length(unique(df$slice)))
  }
  
  return(list(plots=pltlst, dayGTtest=wilcoxreturn))
  
}


createOuts <- function(datfls, formnms, dapithresh=1000) {
  outputdat <- lapply(datfls, read_xlsx)
  outnms1 <- lapply(strsplit(datfls, "/"), function(x) {
    # protein/stain name is now a parent directory. extract that name and add onto the ids with a "_"
    paste(x[c((length(x)-1), (length(x)-2))], collapse="_")
  })
  names(outputdat) <- unlist(outnms1)
  # add sample id to data frames
  outputdat2 <- lapply(names(outputdat), function(x) {
    tmp <- outputdat[[x]]
    id <- unlist(lapply(strsplit(x, "_"), `[[`, 1))
    type <- unlist(lapply(strsplit(x, "_"), `[[`, 2))
    tmp$id <- id
    tmp$type <- type
    return(tmp)
  })
  
  names(outputdat2) <- unlist(outnms1)
  
  sumOuts <- lapply(outputdat2, runCountSum, funcrun="split", formnms=formnms, dapithresh=dapithresh)
  names(sumOuts) <- unlist(outnms1)
  # full lung counts
  fullOuts <- lapply(outputdat2, runCountSum, funcrun="full", formnms=formnms, dapithresh=dapithresh)
  names(fullOuts) <- unlist(outnms1)
  
  return(list(split_section_output=sumOuts, full_slide_output=fullOuts))
  
}


makePlotdf <- function(sumout, type, minfo, sex_analyze=NULL) {
  SumOutUse <- sumout[grep(type, names(sumout))]
  
  stopifnot("No data found specified type"=length(SumOutUse) > 0)
  
  # check if sep section or full lung analysis was used
  # (separate section data would have had the summaries calculated, which were stored as another list element)
  thechk <- unique(unlist(lapply(SumOutUse, class)))
  if ("list" %in% thechk) {
    alldfs <- bind_rows(lapply(SumOutUse, `[[`, 1))
  } else if ("data.frame" %in% thechk) {
    alldfs <- bind_rows(SumOutUse)
  }
  
  if (class(minfo$ID) == "numeric") {
    alldfs$id <- as.numeric(alldfs$id)
  } else if (class(minfo$ID) == "character") {
    alldfs$id <- as.character(alldfs$id)
  }
  
  alldfs2 <- alldfs %>% 
    left_join(minfo, by=c("id"="ID")) %>% 
    filter(Sex %in% sex_analyze)
  
  return(alldfs2)
  
}


# proccess each via a function
processSection <- function(sumout, slice, type, minfo, sex_analyze) {
  if (length(sumout) == 0) {
    return(NULL)
  }
  # check if data exists for this stain
  if (!any(grepl(type, names(sumout)))) {
    return(NULL)
  }
  df <- makePlotdf(sumout = sumout, type = type, minfo = mouseinfo, sex_analyze = sex_analyze)
  df$slice <- slice
  return(df)
}


processDF <- function(df, calc_key, start_idx, x_var="Day", col_var=NULL, show_p_correct=TRUE, plotallagg, ptsize=0.9, pvalsize=5, labelsize=12, rm_outliers=TRUE,
                      make_manual_plot=FALSE, ids_rm=NULL, makemainfig=FALSE) {
  if (is.null(df)) {
    return(NULL)
  }
  
  # make 'orig sex' column so I can refer back to the actual sex for the new 'combined' ones
  df$sex_orig <- df$Sex
  df$sex_orig <- factor(df$sex_orig, levels=c("M", "F"))
  
  df_copy <- df
  df_copy$Sex <- "Combined"
  df_new <- rbind(df, df_copy)
  
  # account for new sex_orig param
  plot_vars <- names(df_new)[start_idx:(ncol(df_new) - 1)]
  # remove where the main stain is in the denominator
  plot_vars <- plot_vars[!grepl("^lcn2_express|^krt8_express", plot_vars, perl=T)]
  # use only variables that aren't just all NAs or 0s
  chkvars <- unlist(lapply(plot_vars, function(x) {all(is.na(df_new[,x]) | df_new[,x] == 0)}))
  names(chkvars) <- plot_vars
  plot_vars <- names(chkvars[which(!chkvars)])
  plot_lst <- plot_vars %>% 
    map(~ plotInLoop_multi(df_new, .x, x_var=x_var, col_var=col_var, calcKey = calc_key, show_p_correct=show_p_correct, rm_outliers=rm_outliers, plotfull = FALSE, ptsize=ptsize,
                           plotallagg=plotallagg, pvalsize=pvalsize, labelsize=labelsize, make_manual_plot=make_manual_plot, ids_rm=ids_rm, makemainfig=makemainfig)) %>% 
    try(., silent = TRUE)
  
  return(list(df=df_new, plot_lst=plot_lst, plot_vars=plot_vars))
}
