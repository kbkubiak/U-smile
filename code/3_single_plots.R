# for plotting one U-smile plot and one PIW plot #

# U-smile plot of RB coeffs for one new model ####
# arguments: ref_model, new_model

plot.single.usmile.rb <- function(ref_model, new_model){
  # calculate RB
  y <- ref_model$y
  p_ref <- ref_model$fitted.values
  res_sq_ref <- (ref_model$y - ref_model$fitted.values)^2
  p <- new_model$fitted.values
  delta <- p - p_ref
  flag <- ifelse(delta > 0, 'up',
                 ifelse(delta < 0, 'dw', 'c'))
  subclass <- ifelse(y == 0 & flag == 'dw',                      'nonev_be',
                     ifelse(y == 0 & flag == 'up',               'nonev_wo',
                            ifelse(y == 1 & flag == 'dw',        'event_wo',
                                   ifelse(y == 1 & flag == 'up', 'event_be', 'unknown'))))
  res_sq <- (new_model$fitted.values - y)^2
  
  probs <- data.frame(y, p_ref, p, subclass, res_sq_ref, res_sq)
  
  SS_nonev_ref <- sum(probs$res_sq_ref[probs$y == 0])
  SS_event_ref <- sum(probs$res_sq_ref[probs$y == 1])
  
  SS_nonev_dw_ref <- sum(probs$res_sq_ref[probs$subclass == 'nonev_be']) 
  SS_nonev_up_ref <- sum(probs$res_sq_ref[probs$subclass == 'nonev_wo']) 
  SS_event_dw_ref <- sum(probs$res_sq_ref[probs$subclass == 'event_wo']) 
  SS_event_up_ref <- sum(probs$res_sq_ref[probs$subclass == 'event_be']) 
  SS_nonev_dw     <- sum(probs$res_sq[probs$subclass == 'nonev_be']) 
  SS_nonev_up     <- sum(probs$res_sq[probs$subclass == 'nonev_wo']) 
  SS_event_dw     <- sum(probs$res_sq[probs$subclass == 'event_wo']) 
  SS_event_up     <- sum(probs$res_sq[probs$subclass == 'event_be']) 
  delta_SS_nonev_dw <- SS_nonev_dw_ref - SS_nonev_dw
  delta_SS_nonev_up <- SS_nonev_up_ref - SS_nonev_up
  delta_SS_event_dw <- SS_event_dw_ref - SS_event_dw
  delta_SS_event_up <- SS_event_up_ref - SS_event_up
  
  RB_nonev_be <-  (delta_SS_nonev_dw) / SS_nonev_ref
  RB_nonev_wo <- -(delta_SS_nonev_up) / SS_nonev_ref
  RB_event_wo <- -(delta_SS_event_dw) / SS_event_ref
  RB_event_be <-  (delta_SS_event_up) / SS_event_ref
  
  data.rb <- data.frame(subclass = c('nonev_be',
                                     'nonev_wo',
                                     'event_wo',
                                     'event_be'),
                        value = c(RB_nonev_be,
                                  RB_nonev_wo,
                                  RB_event_wo,
                                  RB_event_be),
                        coefficient = rep('RB', 4))
  
  # plot settings
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
  
  # U-smile plot
  data.rb %>% 
    ggplot(aes(x = subclass, y = value, group = coefficient)) +
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
          axis.title.y = element_blank() ) -> single_usmile
  
  return(single_usmile)
}

# PIW plot for one new model ####
# arguments: ref_model, new_model
plot.single.piw <- function(ref_model, new_model){
  y <- ref_model$y
  p_ref <- ref_model$fitted.values
  res_sq_ref <- (ref_model$y - ref_model$fitted.values)^2
  p <- new_model$fitted.values
  delta <- p - p_ref
  flag <- ifelse(delta > 0, 'up',
                 ifelse(delta < 0, 'dw', 'c'))
  subclass <- ifelse(y == 0 & flag == 'dw',                      'nonev_be',
                     ifelse(y == 0 & flag == 'up',               'nonev_wo',
                            ifelse(y == 1 & flag == 'dw',        'event_wo',
                                   ifelse(y == 1 & flag == 'up', 'event_be', 'unknown'))))
  res_sq <- (new_model$fitted.values - y)^2
  
  probs <- data.frame(y, p_ref, p, subclass, res_sq_ref, res_sq)
  
  # plot settings
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
  
  # PIW plot
  probs %>% 
    ggplot(aes(p_ref, p)) +
    geom_abline(intercept = 0, slope = 1, color = 'grey45') +
    geom_point(aes(col = subclass, fill = subclass), shape = 21) +
    scale_color_manual(values = usmile_colors,
                       breaks = subclass_order,
                       labels = usmile_labels) +
    scale_fill_manual(values = usmile_fills,
                      breaks = subclass_order,
                      labels = usmile_labels) +
    scale_x_continuous(breaks = c(0, 0.5, 1)) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    theme(aspect.ratio = 1) -> single_piw
  
  return(single_piw)
}