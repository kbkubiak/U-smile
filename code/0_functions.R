# candidate predictors - pretty names ####
train_dataset <- read_excel('data/training_dataset.xlsx')

all_vars <- colnames(train_dataset)

all_vars_pretty <- c("Disease", "Age", "Sex", "Chest pain", "Blood pressure", "Cholesterol", "Glucose", "ECG", "Heart rate", "Exercise angina", "ST depression", 
                     "Location",
                     "Rnd normal", 
                     "Rnd uniform", 
                     "Rnd exponential", 
                     "Rnd Bernoulli",  
                     "Rnd binomial",
                     "Rnd Poisson", 
                     "Str Rnd normal", 
                     "Str Rnd uniform", 
                     "Str Rnd exponential", 
                     "Str Rnd Bernoulli" ,
                     "Str Rnd binomial",
                     "Str Rnd Poisson",
                     "Rnd normal 0.1", 
                     "Rnd normal 0.2",  
                     "Rnd normal 0.3",
                     "Rnd normal 0.4", 
                     "Rnd normal 0.5",  
                     "Rnd normal 0.6",
                     "Rnd normal 0.7", 
                     "Rnd normal 0.8", 
                     "Rnd normal 0.9" ,
                     "Str Rnd normal 0.1", 
                     "Str Rnd normal 0.2",  
                     "Str Rnd normal 0.3",
                     "Str Rnd normal 0.4", 
                     "Str Rnd normal 0.5",  
                     "Str Rnd normal 0.6",
                     "Str Rnd normal 0.7", 
                     "Str Rnd normal 0.8", 
                     "Str Rnd normal 0.9" 
                     )
usmile_all_vars_pretty <- c("Disease", "Age", "Sex", "Chest~pain", "Blood~pressure", "Cholesterol", "Glucose", "ECG", "Heart~rate", "Exercise~angina", "ST~depression",
                            "Location",
                            "Rnd~normal", 
                            "Rnd~uniform", 
                            "Rnd~exponential", 
                            "Rnd~Bernoulli",
                            "Rnd~binomial",
                            "Rnd~Poisson", 
                            "Str~Rnd~normal", 
                            "Str~Rnd~uniform", 
                            "Str~Rnd~exponential", 
                            "Str~Rnd~Bernoulli" ,
                            "Str~Rnd~binomial",
                            "Str~Rnd~Poisson",
                            "Rnd~normal~0.1", 
                            "Rnd~normal~0.2",
                            "Rnd~normal~0.3",
                            "Rnd~normal~0.4", 
                            "Rnd~normal~0.5",
                            "Rnd~normal~0.6",
                            "Rnd~normal~0.7", 
                            "Rnd~normal~0.8", 
                            "Rnd~normal~0.9",
                            "Str~Rnd~normal~0.1", 
                            "Str~Rnd~normal~0.2", 
                            "Str~Rnd~normal~0.3",
                            "Str~Rnd~normal~0.4",
                            "Str~Rnd~normal~0.5",
                            "Str~Rnd~normal~0.6",
                            "Str~Rnd~normal~0.7",
                            "Str~Rnd~normal~0.8",
                            "Str~Rnd~normal~0.9"
)

# extract names and indices of candidate variables
new.vars <- function(ref_vars){
  new_vars <- all_vars[!all_vars %in% c('id', 'location', 'disease', ref_vars)]
  return(new_vars)
}

new.vars.index <- function(new_vars){
  index <- all_vars %in% new_vars
  return(index)
}


# usmile facet titles
# comps - comparisons
usmile.pretty.models <- function(new_vars_index, comps){
  result <- paste0(usmile_all_vars_pretty[new_vars_index],
                   ifelse(comps$LRT_p.value < 0.05, "~'*'", ""), 
                   ifelse(comps$AUC_p.value < 0.05, "^{', #'}", ""))
  return(result)
}


# ROC & PID plots titles
# copmps - comparisons
plots.pretty.models <- function(new_vars_index, comps){
  result <- paste0(all_vars_pretty[new_vars_index],  
                   ifelse(comps$LRT_p.value < 0.05, " *", ""), 
                   ifelse(comps$AUC_p.value < 0.05, "<sup>, #</sup>", ""))
  return(result)
}

# models ####
build.ref.model <- function(ref_vars, dataset){
  ref_model <- glm(as.formula(paste('disease ~ ', paste(ref_vars, sep="", collapse=" + ") )), data = dataset, family = binomial) 
  return(ref_model)
}

build.new.models <- function(ref_vars, new_vars, dataset){
  new_models <- list()
  for(i in seq_along(new_vars)){
    
    new_models[[i]] <- glm(as.formula(paste0(paste('disease ~ ', paste(ref_vars, sep="", collapse=" + ")), '+', new_vars[i])),
                           data = dataset, family = binomial)
  }
  return(new_models)
}

# ROCs ####
build.ROC.ref <- function(ref_model){
  roc_ref <- roc(ref_model$y, ref_model$fitted.values, direction = '<')
  return(roc_ref)
}

# 
build.ROC.new <- function(new_models){
  roc_new <- lapply(new_models, function(x) roc(x$y, x$fitted.values, direction = '<'))
  return(roc_new)
}

build.ROC.ref.test <- function(ref_model, test_dataset){
  roc_ref_test <- roc(test_dataset$disease, predict(ref_model, test_dataset, type = 'r'), direction = '<')
  return(roc_ref_test)
}

build.ROC.new.test <- function(new_models, test_dataset){
  roc_new_test <- lapply(new_models, function(x) roc(test_dataset$disease, predict(x, test_dataset, type = 'r'), direction = '<' ) )
  return(roc_new_test)
}

# comparisons ####

compare.ref.with.new.models <- function(ref_model, new_models, roc_ref, roc_new, new_vars){
  
  auc_ref <- roc_ref$auc
  LRT_p.value <- sapply(new_models, function(x) lrtest(ref_model, x)$`Pr(>Chisq)`[2] )
  AUC_p.value <- sapply(roc_new, function(x) roc.test(roc_ref, x)$p.value  )
  AUC <- sapply(roc_new, function(x) x$auc)
  deltaAUC <- AUC - auc_ref
  
  result <- data.frame(model = new_vars, 
                       LRT_p.value = round(LRT_p.value, 6),
                       AUC_p.value = round(AUC_p.value, 6),
                       AUC         = round(AUC, 6),
                       deltaAUC    = round(deltaAUC, 6)
  )

  rownames(result) <- NULL 
  result$LRT_p.value <- as.numeric(result$LRT_p.value)
  result$AUC_p.value <- as.numeric(result$AUC_p.value)
  result$AUC <- as.numeric(result$AUC)
  result$deltaAUC <- as.numeric(result$deltaAUC)
  
  return(result)
}


# probabilities ####
calculate.probs <- function(ref_model, new_models){
  
  y <- ref_model$y
  p_ref <- ref_model$fitted.values
  res_ref <- abs(ref_model$y - ref_model$fitted.values)
  res_sq_ref <- (ref_model$y - ref_model$fitted.values)^2
  n <- length(y)
  
  probs <- list()
  for(i in seq_along(new_models)){
    p <- new_models[[i]]$fitted.values
    delta <- p - p_ref
    flag <- ifelse(delta > 0, 'up',
                   ifelse(delta < 0, 'dw', 'c'))
    subclass <- ifelse(y == 0 & flag == 'dw',                      'nonev_be',
                       ifelse(y == 0 & flag == 'up',               'nonev_wo',
                              ifelse(y == 1 & flag == 'dw',        'event_wo',
                                     ifelse(y == 1 & flag == 'up', 'event_be', 'unknown'))))
    res <- abs(new_models[[i]]$fitted.values - y)
    res_sq <- (new_models[[i]]$fitted.values - y)^2
    
    temp <- data.frame(y, 
                       p_ref,
                       p, 
                       subclass,
                       res_sq_ref,
                       res_sq)
    probs[[i]] <- temp
    
  }
  names(probs) <- names(new_models)
  return(probs)
}

calculate.probs.test <- function(ref_model, new_models, test_dataset){ 
  
  y <- test_dataset$disease
  p_ref <- predict(ref_model, test_dataset, type = 'response')
  res_ref <- abs(y - p_ref)
  res_sq_ref <- (y - p_ref)^2
  
  probs_test <- list()
  for(i in seq_along(new_models)){
    
    p <- predict(new_models[[i]], test_dataset, type = 'response')
    delta <- p - p_ref
    flag <- ifelse(delta > 0, 'up',
                   ifelse(delta < 0, 'dw', 'c'))
    subclass <- ifelse(y == 0 & flag == 'dw',                      'nonev_be',
                       ifelse(y == 0 & flag == 'up',               'nonev_wo',
                              ifelse(y == 1 & flag == 'dw',        'event_wo',
                                     ifelse(y == 1 & flag == 'up', 'event_be', 'unknown'))))
    res <- abs(p - y)
    res_sq <- (p - y)^2
    
    temp <- data.frame(y,
                       p_ref,
                       p,
                       subclass,
                       res_sq_ref,
                       res_sq)
    probs_test[[i]] <- temp
  }
  
  names(probs_test) <- names(new_models)
  return(probs_test)
}

#
# RB & I ####
# col1 must be model, cols 2:9 must be RB, I - needed for prepare.rbi.for.usmile()
calculate.rbi <- function(probs){
  
  n0 <- sum(probs[[1]]$y == 0)
  n1 <- sum(probs[[1]]$y == 1)
  SS_nonev_ref <- sum(probs[[1]]$res_sq_ref[probs[[1]]$y == 0])
  SS_event_ref <- sum(probs[[1]]$res_sq_ref[probs[[1]]$y == 1])
  
  SS_nonev_dw_ref <- sapply(probs, function(x) sum(x$res_sq_ref[x$subclass == 'nonev_be']) )
  SS_nonev_up_ref <- sapply(probs, function(x) sum(x$res_sq_ref[x$subclass == 'nonev_wo']) )
  SS_event_dw_ref <- sapply(probs, function(x) sum(x$res_sq_ref[x$subclass == 'event_wo']) )
  SS_event_up_ref <- sapply(probs, function(x) sum(x$res_sq_ref[x$subclass == 'event_be']) )
  SS_nonev_dw     <- sapply(probs, function(x) sum(x$res_sq[x$subclass == 'nonev_be']) )
  SS_nonev_up     <- sapply(probs, function(x) sum(x$res_sq[x$subclass == 'nonev_wo']) )
  SS_event_dw     <- sapply(probs, function(x) sum(x$res_sq[x$subclass == 'event_wo']) )
  SS_event_up     <- sapply(probs, function(x) sum(x$res_sq[x$subclass == 'event_be']) )
  delta_SS_nonev_dw <- SS_nonev_dw_ref - SS_nonev_dw
  delta_SS_nonev_up <- SS_nonev_up_ref - SS_nonev_up
  delta_SS_event_dw <- SS_event_dw_ref - SS_event_dw
  delta_SS_event_up <- SS_event_up_ref - SS_event_up
  
  RB_nonev_be <-  (delta_SS_nonev_dw) / SS_nonev_ref
  RB_nonev_wo <- -(delta_SS_nonev_up) / SS_nonev_ref
  RB_event_wo <- -(delta_SS_event_dw) / SS_event_ref
  RB_event_be <-  (delta_SS_event_up) / SS_event_ref
  I_nonev_be  <- sapply(probs, function(x) sum(x$subclass == 'nonev_be') / n0 )
  I_nonev_wo  <- sapply(probs, function(x) sum(x$subclass == 'nonev_wo') / n0 )
  I_event_wo  <- sapply(probs, function(x) sum(x$subclass == 'event_wo') / n1 )
  I_event_be  <- sapply(probs, function(x) sum(x$subclass == 'event_be') / n1 )
  
  result <- data.frame(model = names(SS_nonev_dw_ref), 
                       RB_nonev_be,
                       RB_nonev_wo,
                       RB_event_wo,
                       RB_event_be,
                       I_nonev_be,
                       I_nonev_wo,
                       I_event_wo,
                       I_event_be)
  rownames(result) <- NULL
  return(result)
}

# x - result of calculate.rbi
# comps - comparisons
prepare.rbi.for.usmile <- function(x, new_vars, new_vars_index, comps){
  results_rbi_long <- pivot_longer(x[, 1:9], -model)
  colnames(results_rbi_long)[2] <- 'coefficient'
  coeff <- ifelse(str_starts(results_rbi_long$coefficient, 'RB'), 'RB', 'I')
  subclass <- str_sub(results_rbi_long$coefficient, -8, -1)
  
  result <- data.frame(model = factor(results_rbi_long$model, levels = new_vars),
                       coefficient = coeff,
                       subclass = subclass,
                       value = results_rbi_long$value)
  
  result$model_pretty <- factor(result$model, levels = new_vars, labels = usmile.pretty.models(new_vars_index, comps))
  
  return(result)
}

# settings for plots ####
subclass_order <- c('nonev_be',
                    'nonev_wo',
                    'event_wo',
                    'event_be')
usmile_colors <- c('nonev_be' = '#0F3C78',
                   'nonev_wo' = '#0F3C78',
                   'event_wo' = '#D51424',
                   'event_be' = '#D51424')
usmile_fills  <- c('nonev_be' = '#0F3C78',
                   'nonev_wo' = '#BED2FA',
                   'event_wo' = '#FBCDB9',
                   'event_be' = '#D51424')
usmile_labels <- c('non-events with better prediction',
                   'non-events with worse prediction',
                   'events with worse prediction',
                   'events with better prediction')
#
# plot U-smiles  ####
plot.usmile <- function(rbi_for_usmile, results_rbi_train, results_rbi_test, title){
  usmile_caption <- "*, P value < 0.05 of the likelihood-ratio test; <sup>#</sup>, P value < 0.05 of the DeLong's test for two correlated ROC curves"

  RB_ymax <- max(with(results_rbi_train, max(RB_nonev_be,
                                             RB_nonev_wo,
                                             RB_event_wo,
                                             RB_event_be)),
                 with(results_rbi_test , max(RB_nonev_be,
                                             RB_nonev_wo,
                                             RB_event_wo,
                                             RB_event_be))) #* 1.1
  
  rbi_for_usmile %>% 
    filter(coefficient == 'RB') %>% 
    ggplot(aes(x = subclass, y = value, group = model)) +
    geom_line() + 
    geom_point(aes(col = subclass, fill = subclass), shape = 21, size = 3, stroke = 1) + 
    scale_x_discrete(limits = subclass_order) +
    scale_color_manual(values = usmile_colors,
                       breaks = subclass_order,
                       labels = usmile_labels) +
    scale_fill_manual(values = usmile_fills,
                      breaks = subclass_order,
                      labels = usmile_labels) + 
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x  = element_blank(),
          axis.text.y  = element_text(size = 9),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),
          strip.text.x = element_text(angle = 90, hjust = 0, size = 10),
          strip.text.y = element_text(size = 13),
          legend.position = 'none',
          plot.title = element_text(size = 14),
          plot.title.position = "plot") +
    ggtitle(title) +
    facet_grid(coefficient ~ model_pretty, labeller = "label_parsed") +
    ylim(c(0, 1.14*RB_ymax)) -> plot_RB1 
  
  rbi_for_usmile %>% 
    filter(coefficient == 'I') %>% 
    ggplot(aes(x = subclass, y = value, group = model)) +
    geom_line() + 
    geom_point(aes(col = subclass, fill = subclass), shape = 21, size = 3, stroke = 1) + 
    scale_x_discrete(limits = subclass_order) +
    scale_color_manual(values = usmile_colors,
                       breaks = subclass_order,
                       labels = usmile_labels) +
    scale_fill_manual(values = usmile_fills,
                      breaks = subclass_order,
                      labels = usmile_labels) + 
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y  = element_text(size = 9),
          axis.title.y = element_blank(),
          strip.text.x = element_blank(),
          strip.text.y = element_text(size = 13),
          legend.position = "bottom",
          plot.caption = element_markdown(size = 9)) +
    facet_grid(coefficient ~ model) +
    ylim(c(0, 1)) +
    labs(caption = usmile_caption) -> plot_I1
  
  g1 <- rbind(ggplotGrob(plot_RB1), 
              ggplotGrob(plot_I1),
              size = 'first')
  return(g1)
}

# plot ROCs ####
plot.roc <- function(roc_ref, roc_new, new_vars_index, comps, title){
  # panel titles
  panel_titles <- plots.pretty.models(new_vars_index, comps)
  # extract legend
  ROC_gg <- ggroc(list(`Reference model` = roc_ref, `New model` = roc_new[[1]])) +
    scale_colour_manual(values = c('black', '#D51424')) +
    theme(legend.title = element_blank()) +
    guides(colour = guide_legend(nrow = 1))
  
  tmp <- ggplot_gtable(ggplot_build(ROC_gg))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  ROC_legend <- tmp$grobs[[leg]]
  
  ROC_plots <- list()
  for(i in seq_along(ROC_new)){
    ROC_plots[[i]] <- ggroc(list(roc_ref, roc_new[[i]]),
                            legacy.axes = TRUE) +
      scale_colour_manual(values = c('black', '#D51424')) +
      geom_abline(intercept = 0, slope = 1, color = 'grey45') +
      scale_x_continuous(breaks = c(0, 0.5, 1)) +
      scale_y_continuous(breaks = c(0, 0.5, 1)) +
      theme(aspect.ratio = 1,
            legend.position = 'none',
            plot.title = element_markdown(size = 8),
            axis.title.x = element_markdown(size = 7), 
            axis.title.y = element_markdown(size = 7), 
            axis.text = element_text(size = 6),
            plot.margin = unit(c(0.2, 0.1, 0, 0.1), 'cm')
      ) +
      ggtitle(panel_titles[i])
  }
  
  roc_caption <- expression(paste("*, P value < 0.05 of the likelihood-ratio test; " ^"#" *", P value < 0.05 of the DeLong's test for two correlated ROC curves"))
  
  grid.arrange(ROC_plots[[1]],  ROC_plots[[7]],  ROC_plots[[13]],  ROC_plots[[2]],  ROC_plots[[8]],  ROC_plots[[14]], ROC_plots[[3]], 
               ROC_plots[[9]],  ROC_plots[[15]],  ROC_plots[[4]], ROC_plots[[10]], ROC_plots[[16]], 
               ROC_plots[[5]], ROC_plots[[11]], ROC_plots[[17]], ROC_plots[[6]], 
               ROC_plots[[12]], ROC_plots[[18]], 
               ROC_legend, 
               layout_matrix = rbind(c(1:3),
                                     c(4:6),
                                     c(7:9),
                                     c(10:12),
                                     c(13:15),
                                     c(16:18),
                                     c(19,19, 19)
               ),
               widths  = rep(1, 3),
               heights = c(1, 1, 1, 1, 1, 1, 0.1), 
               top = textGrob(title,
                              gp = gpar(fontsize = 14),
                              hjust = 0, x = 0),
               bottom = textGrob(roc_caption,
                                 gp = gpar(fontsize = 9),
                                 hjust = 1, x = 1)
  )
}

# plot PIWs ####
plot.piw <- function(probs, new_vars_index, comps, title){
  # panel titles
  panel_titles <- plots.pretty.models(new_vars_index, comps)
  # extract legend 
  PIW_gg <- probs[[1]] %>% 
    ggplot(aes(p_ref, p)) +
    geom_point(aes(col = subclass, fill = subclass), shape = 21, size = 1, stroke = 1) +
    scale_color_manual(values = usmile_colors,
                       breaks = subclass_order,
                       labels = usmile_labels) +
    scale_fill_manual(values = usmile_fills,
                      breaks = subclass_order,
                      labels = usmile_labels) +
    theme(legend.position = "bottom")
  #guides(colour = guide_legend(nrow = 1),
  #       fill = guide_legend(nrow = 1))
  
  tmp <- ggplot_gtable(ggplot_build(PIW_gg))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  PIW_legend <- tmp$grobs[[leg]]
  
  PIW_plots <- list()
  for(i in seq_along(probs)){
    PIW_plots[[i]] <- 
      ggplot(probs[[i]], aes(p_ref, p)) +
      geom_abline(intercept = 0, slope = 1, color = 'grey45') +
      geom_point(aes(col = subclass, fill = subclass), shape = 21, size = 1) +
      scale_color_manual(values = usmile_colors,
                         breaks = subclass_order,
                         labels = usmile_labels) +
      scale_fill_manual(values = usmile_fills,
                        breaks = subclass_order,
                        labels = usmile_labels) +
      labs(x = "*p*<sub>(ref)</sub>", y = "*p*") +
      scale_x_continuous(breaks = c(0, 0.5, 1)) +
      scale_y_continuous(breaks = c(0, 0.5, 1)) +
      theme(aspect.ratio = 1,
            legend.position = 'none',
            plot.title = element_markdown(size = 8),
            axis.title.x = element_markdown(size = 8), 
            axis.title.y = element_markdown(size = 8), 
            axis.text = element_text(size = 6),
            plot.margin = unit(c(0.2, 0.1, 0, 0.1), 'cm')
      ) +
      ggtitle(panel_titles[i])
  }
  
  piw_caption <- expression(paste("*, P value < 0.05 of the likelihood-ratio test; " ^"#" *", P value < 0.05 of the DeLong's test for two correlated ROC curves"))
  
  grid.arrange(PIW_plots[[1]],  PIW_plots[[7]],   PIW_plots[[13]],  PIW_plots[[2]],  PIW_plots[[8]],  PIW_plots[[14]], PIW_plots[[3]], 
               PIW_plots[[9]],  PIW_plots[[15]],  PIW_plots[[4]],   PIW_plots[[10]], PIW_plots[[16]], 
               PIW_plots[[5]],  PIW_plots[[11]],  PIW_plots[[17]],  PIW_plots[[6]], 
               PIW_plots[[12]], PIW_plots[[18]],
               PIW_legend,
               layout_matrix = rbind(c(1:3),
                                     c(4:6),
                                     c(7:9),
                                     c(10:12),
                                     c(13:15),
                                     c(16:18),
                                     c(19,19, 19)
               ),
               widths  = rep(1, 3),
               heights = c(1, 1, 1, 1, 1, 1, 0.1), 
               top = textGrob(title, 
                              gp = gpar(fontsize = 14),
                              hjust = 0, x = 0),
               bottom = textGrob(piw_caption,
                                 gp = gpar(fontsize = 9),
                                 hjust = 1, x = 1)
  )
}