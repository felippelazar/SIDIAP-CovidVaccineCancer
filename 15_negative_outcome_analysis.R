# ============================================================================ #
# 15. REM Analysis - Negative Outcomes Analysis                                #
# Author: Felippe Lazar, IDIAP Jordi Gol, 2023 #
# ============================================================================ #

library(EmpiricalCalibration)
library(here)
library(glue)
library(tidyverse)
library(openxlsx)
library(grid)
library(gridExtra)

# Setting Temporary WD
# Loading Table of Negative Outcome Results
getNegativeOutcomesResults <- function(concept_id, file_location){
      
      file_full_path = glue(file_location, 'outcome_three_periods', concept_id, '.csv') 
      
      return(read.table(file_full_path, sep=';', header = T) %>%
            select(term, estimate, std.error, conf.low, conf.high, n_event) %>%
            mutate(est_hr.conf.interval = sprintf('%.2f (%.2f - %.2f)', estimate, conf.low, conf.high)) %>%
            mutate(ConceptId = concept_id))
}

NCO <- read.csv(file=here("NCO.csv"), sep = ";") # File with the conditions chosen to assess residual confounding
# Choosing Only Validated NCO
NCO <- NCO[c(1:43), ] # Original
# NCO <- NCO[c(1:54), ] # Expanded Negative Outcomes

negOutcomes12 <- lapply(NCO$ConceptId, getNegativeOutcomesResults, glue('Results/', 'dose_12/negative outcomes/'))
negOutcomes12 <- do.call(bind_rows, negOutcomes12)
negOutcomes12 <- negOutcomes12 %>% left_join(NCO)

negOutcomes3 <- lapply(NCO$ConceptId, getNegativeOutcomesResults, glue('Results/', 'dose_3/negative outcomes/'))
negOutcomes3 <- do.call(bind_rows, negOutcomes3)
negOutcomes3 <- negOutcomes3 %>% left_join(NCO)

# Plot of Negative Control Outcome Estimates
plot_neg_outcomes <- bind_rows(negOutcomes12, negOutcomes3) %>% filter(n_event != 0 & term != 'periodno-vax') %>% 
      filter((term != 'periodV1 0-14D') & (term != 'periodV3 0-14D')) %>%
      ggplot(aes(y=reorder(OutcomeName, n_event), x=estimate, size=n_event)) +
      geom_errorbarh(aes(xmax = conf.high, xmin = conf.low), size = .5, height = .2, color = "gray50") +
      geom_point(aes(group=term), alpha = 0.8) + scale_size(range = c(0, 3)) +
      geom_vline(xintercept=1, size = .25, linetype = "dashed") +
      facet_grid(.~term) + scale_x_log10() + coord_cartesian(xlim=c(0.3,3)) + 
      theme_minimal() +
      theme(panel.grid.minor = element_blank(), text = element_text(family='Helvetica'),
            legend.position = "bottom", strip.background =element_rect(fill="lightgrey")) +
      labs(x = 'Hazard Ratio', y = 'Outcomes', size = 'Number of Events')

ggsave('Figures/NegOutcomesPlot_Figure_21022.png', height = 150, width = 250, units = "mm", dpi = "print")

plot_neg_outcomes_text <- bind_rows(negOutcomes12, negOutcomes3) %>% filter(n_event != 0 & term != 'periodno-vax' & estimate < 10000) %>% 
      filter((term != 'periodV1 0-14D') & (term != 'periodV3 0-14D')) %>%
      ggplot(aes(y=reorder(OutcomeName, n_event), x=term, label = est_hr.conf.interval,
                 fill = ifelse(estimate >= 1, TRUE, FALSE))) +
      geom_tile(alpha = 1, color = 'black', fill = 'white') + 
      geom_text() + theme_minimal() + 
      #scale_fill_manual(values = c('darkgray', 'white')) + 
      theme(panel.grid.minor = element_blank(), text = element_text(family='Helvetica'),
            legend.position = "none", strip.background =element_rect(fill="lightgrey")) +
      labs(x = 'Hazard Ratio (95% Confidence Interval)', y = 'Outcomes')

ggsave('Figures/NegOutcomesPlot_Figure_2005_Text.png', height = 150, width = 250, units = "mm", dpi = "print")

neg_outcomes_combined_figure <- cowplot::plot_grid(plot_neg_outcomes, 
                                         plot_neg_outcomes_text, 
                                         ncol = 1, 
                                         labels = c("A", "B"))

ggsave('Figures/negative_outcomes_combined.pdf', height = 29, width = 21, units = 'cm')

# Plotting Calibration
negOutcomes12_toPlot <- negOutcomes12 %>% filter(term == 'periodV1 14D+') %>%
      select(estimate, std.error) %>% drop_na()

v1_14 <- plotCiCalibrationEffect(log(negOutcomes12_toPlot$estimate), negOutcomes12_toPlot$std.error, rep(0,nrow(negOutcomes12_toPlot)))

negOutcomes12_toPlot <- negOutcomes12 %>% filter(term == 'periodV2 7D+') %>%
      select(estimate, std.error) %>% drop_na()

v2_7 <- plotCiCalibrationEffect(log(negOutcomes12_toPlot$estimate), negOutcomes12_toPlot$std.error, rep(0,nrow(negOutcomes12_toPlot)))

negOutcomes3_toPlot <- negOutcomes3 %>% filter(term == 'periodV3 14-60D') %>%
      select(estimate, std.error) %>% drop_na()

v3_14 <- plotCiCalibrationEffect(log(negOutcomes3_toPlot$estimate), negOutcomes3_toPlot$std.error, rep(0,nrow(negOutcomes3_toPlot)))

negOutcomes3_toPlot <- negOutcomes3 %>% filter(term == 'periodV3 60+') %>%
      select(estimate, std.error) %>% drop_na()

v3_60 <- plotCiCalibrationEffect(log(negOutcomes3_toPlot$estimate), negOutcomes3_toPlot$std.error, rep(0,nrow(negOutcomes3_toPlot)))

lay = rbind(c(1, 3), c(2, 4))

tt <- grid.arrange(v1_14, v2_7, v3_14, v3_60, 
                   layout_matrix = lay)

tt <- egg::ggarrange(v1_14, v2_7, v3_14, v3_60, ncol = 2,
                     labels = c('A', 'B', 'C', 'D'),
                     label.args = list(gp=gpar(fontface='bold', fontsize=30), x=unit(2,"line"), hjust=-0.5, vjust=1.5))

ggsave("Figures/NegOutcomes_CalibrationPlot_210224.png", plot = tt, height = 260, width = 260, units = "mm", dpi = "print")

# Loading Main Outcome Results
getCalibratedResults <- function(neg_outcomes_dataframe, main_results_dataframe, term_to_filter){
      neg_outcomes_dataframe <- neg_outcomes_dataframe %>%
            filter(term == term_to_filter) %>% drop_na()
      main_results_dataframe <- main_results_dataframe %>%
            filter(term == term_to_filter) %>% drop_na()
      
      model <- fitSystematicErrorModel(log(neg_outcomes_dataframe$estimate), neg_outcomes_dataframe$std.error, rep(0,nrow(neg_outcomes_dataframe))) # use NCO to fit model
      result_full <- calibrateConfidenceInterval(log(main_results_dataframe$estimate), main_results_dataframe$std.error, model, ciWidth = 0.95) # use model to calibrate treatment estimates
      
      result_full <- result_full %>%
            mutate(term = term_to_filter) %>%
            mutate(estimate = exp(logRr), 
                   conf.low = exp(logLb95Rr),
                   conf.high = exp(logUb95Rr),
                   std.error = seLogRr,
                   type = 'calibrated_estimated') %>%
            select(term, estimate, std.error, conf.low, conf.high, type)
      
      return(result_full)
}

outcomeHR12 <- read.table(glue('Results/',  'dose_12/rem_main_analysis/outcome_hosp_period_three.csv'),  sep=';', header = T) %>%
      mutate(OutcomeName = 'COVID-19 Hosp.')

outcomeHR3 <- read.table(glue('Results/',  'dose_3/rem_main_analysis/outcome_hosp_period_three.csv'),  sep=';', header = T) %>%
      mutate(OutcomeName = 'COVID-19 Hosp.')

tidy_results <- bind_rows(getCalibratedResults(negOutcomes12, outcomeHR12, 'periodV1 14D+'),
          getCalibratedResults(negOutcomes12, outcomeHR12, 'periodV2 7D+'),
          getCalibratedResults(negOutcomes3, outcomeHR3, 'periodV3 14-60D'),
          getCalibratedResults(negOutcomes3, outcomeHR3, 'periodV3 60+')) %>%
      mutate(estimate_ve = if_else(estimate>1, -(1-(1/estimate))*100, (1-estimate)*100),
             conf.high_ve = if_else(conf.low>1, -(1-(1/conf.low))*100, (1-conf.low)*100),
             conf.low_ve= if_else(conf.high>1, -(1-(1/conf.high))*100, (1-conf.high)*100)) %>%
      mutate(est_ve.conf.interval = sprintf('%.2f (%.2f - %.2f)', estimate_ve, conf.low_ve, conf.high_ve))

bind_rows(negOutcomes12 %>% mutate(data = '1-2 dose'), negOutcomes3 %>% mutate(data = '3 dose'))

# Exporting To Excel
wb <- createWorkbook()
addWorksheet(wb, "calibrated_results")
writeData(wb, sheet = 1, x = tidy_results)
addWorksheet(wb, "negative_outcomes_1-2nd")
writeData(wb, sheet = 2, x = negOutcomes12 %>% mutate(data = '1-2 dose'))
addWorksheet(wb, "negative_outcomes_3rd")
writeData(wb, sheet = 3, x = negOutcomes3 %>% mutate(data = '3 dose'))
saveWorkbook(wb, glue('Results Tidy/',  'neg_outcomes_tidy.xlsx'), overwrite = TRUE)

             