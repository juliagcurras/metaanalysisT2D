##############################################################################-

#################       METANALYSIS - MAIN STRUCTURE       ###################-

##############################################################################-


# Julia G Currás - 2025/05/28
rm(list=ls())
graphics.off()
library(dplyr)
library(Biostatech)
library(meta)
library(metafor)

# SMD and SE estimation
getSMD_escalc <- function(df, proteina){
  # Get data
  df <- df[proteina,]
  dm <- df %>% dplyr::select(starts_with("T2DM_")) %>% as.vector %>% unlist
  nDM <- length(dm)
  control <- df %>% dplyr::select(starts_with("Control_")) %>% as.vector %>% unlist
  nControl <- length(control)
  
  # get SMD and SE pooled
  result <-  metafor::escalc(measure = "SMD", 
                    m1i = mean(dm), m2i = mean(control), 
                    sd1i = sd(dm), sd2i = sd(control), 
                    n1i = nDM, n2i = nControl, 
                    add = 0, 
                    vtype = "LS")
  
  return(list(SMD = result$yi, Var = result$vi, 
              mean.e =  mean(dm), mean.c = mean(control),
              sd.e = sd(dm), sd.c = sd(control),
              n.e = nDM, n.c = nControl) )
}


# Loading curated, raw data ####
lista <- readRDS(file = "2_allRawDataProcessed.rds")

# Loading interesting proteins ####
hojas <- readxl::excel_sheets("3_proteinsForMetanalisis.xlsx")
allDatasets <- list()
for (i in hojas){
  allDatasets[[i]] <- readxl::read_xlsx(path = "3_proteinsForMetanalisis.xlsx", 
                  sheet = i, col_names = T)
  
}

# Loading extra data ####
basal <- readxl::read_xlsx(path = "../3_Extraction/Database.xlsx", 
                                         sheet = "1. Basal data", col_names = T)
sampleInfo <- basal %>% 
  dplyr::select(ID, Reference, `Type of sample (grouped)`, `Data acquisition mode`) %>%
  rename(Sample = `Type of sample (grouped)`)
sampleInfo$Sample <- as.character(sampleInfo$Sample)
sampleInfo$Sample <- ifelse(sampleInfo$Sample == "Pancreas biopsy (pancreatic island cells)", 
                                    "Pancreas", 
                                    ifelse(sampleInfo$Sample == "Extracellular vesicles (blood)", "Blood (EV)", 
                                    ifelse(sampleInfo$Sample == "Serum", "Blood (Serum)", 
                                    ifelse(sampleInfo$Sample == "Plasma", "Blood (Plasma)", 
                                           sampleInfo$Sample))))


# Meta-analysis ####
allMeta_data_esc <- list()

for(j in hojas){
  data <- as.data.frame(allDatasets[[j]])
  proteinas <- unique(data$ProteinID)
  nombre <- gsub(j, pattern = "[^A-Za-z]", replacement = "")
  for (proteina in proteinas){
    # Select info
    ids <- data[which(data$ProteinID == proteina), "ID"]
    selectedList <- lista[ids]
    
    # Get SDM and SE
    result <- lapply(selectedList, getSMD_escalc, proteina)
    df <- data.frame(
      ID = names(result),
      SMD = sapply(result, "[[", 1),
      Var = sapply(result, "[[", 2),
      mean.e = sapply(result, "[[", 3),
      mean.c = sapply(result, "[[", 4),
      sd.e = sapply(result, "[[", 5), 
      sd.c = sapply(result, "[[", 6),
      n.e = sapply(result, "[[", 7), 
      n.c = sapply(result, "[[", 8),
      row.names = NULL
    )
    
    df <- merge(df, sampleInfo, by = "ID", all.x = T)
    
    finalName <- paste0(nombre, "_", proteina)
    allMeta_data_esc[[finalName]] <- df
    
  } 
}


allMeta_data_esc_def <- allMeta_data_esc[1:9]

#######============== Results of metaanalysys
library(dmetar)
library(ggrepel)


#####================= Model meta ==============================================

models_results_meta <- function(data) {
  
  model_meta <- meta::metacont(n.e, mean.e, sd.e,
                               n.c, mean.c, sd.c,
                               studlab = Reference, data = data,
                               sm = "SMD", MH.exact = TRUE,
                               fixed = FALSE, random = TRUE, method.tau = "REML", hakn = TRUE, 
                               method.random.ci = "KR", method.predict = "KR", prediction = TRUE)
  
  model <- model_meta
  
  res <- data.frame(Estimate = model$TE.random, ci.lb = model$lower.random, ci.ub = model$upper.random,
                    pi.lb = model$lower.predict, pi.ub = model$upper.predict)
  
  model_metafor <- metafor::rma(measure = "SMD", 
                                yi = SMD, 
                                vi = Var, 
                                method="REML", 
                                test="hksj",
                                data = data)
  
  #----- OUTLIERS ------------------------------------------------------------------------
  out_find <- find.outliers(model)
  out_find_studies <- out_find$out.study.random
  
  res_out <- data.frame(Estimate = out_find$m.random$TE.random, ci.lb = out_find$m.random$lower.random, ci.ub = out_find$m.random$upper.random,
                        pi.lb = out_find$m.random$lower.predict, pi.ub = out_find$m.random$upper.predict)
  
  heter_mod <- data.frame(tau2 = model$tau2, tau2.lb = model$lower.tau2, tau2.ub = model$upper.tau2,
                          I2 = model$I2, I2.ci.lb = model$lower.I2, I2.ci.ub = model$upper.I2,
                          H2 = model$H, H2.ci.lb = model$lower.H, H2.ci.ub = model$upper.H)%>%
    dplyr::mutate_if(is.numeric, round, 3)
  heter_out <- data.frame(tau2 = out_find$m.random$tau2, tau2.lb = out_find$m.random$lower.tau2, tau2.ub = out_find$m.random$upper.tau2,
                          I2 = out_find$m.random$I2, I2.ci.lb = out_find$m.random$lower.I2, I2.ci.ub = out_find$m.random$upper.I2,
                          H2 = out_find$m.random$H, H2.ci.lb = out_find$m.random$lower.H, H2.ci.ub = out_find$m.random$upper.H)%>%
    dplyr::mutate_if(is.numeric, round, 3)
  heter <- data.frame(Model = c("with outliers", "without outliers"), rbind(heter_mod, heter_out))
  outliers <- list(out_find, out_find$out.study.random, heter)
  names(outliers) <- c("model_out", "out", "heterogeneity")
  
  #----- BIAS ------------------------------------------------------------------------
  bias <- metabias(model, method.bias = "Egger", k.min = 4)
  bias_test <- paste0("t(df) = ", round(bias$statistic, 3), "(", round(bias$df, 0), ")", ", p-value = ", round(bias$pval, 4))
  bias_table <- data.frame(bias$estimate) 
  fail_safe <- fsn(x = SMD, vi = Var, data = data, type = "Rosenberg")
  fail_safe <- paste0("Observed Significance Level:",  round(fail_safe$pval, 4),", ", 
                      "Target Significance Level:",  fail_safe$alpha,", ", 
                      "Fail safe N: ", fail_safe$fsnum)
  excess_signif <- tes(x = SMD, vi = Var, data = data, test = "chi2")
  excess_signif_res <- paste0("Test of Excess Significance: ", round(excess_signif$pval, 4), " (chi-square: ", 
                              round(excess_signif$X2, 4), ", Observed Number / Expected Number: ", 
                              round(excess_signif$OEratio, 3), ")")
  bias_res <- list(bias, bias_test, bias_table, fail_safe, excess_signif, excess_signif_res)
  names(bias_res) <- c("model_bias", "test", "estimate_bias", "fail_safe", "excess_signif", "excess_signif_res")
  #----- TRIM and FILL ------------------------------------------------------------------------
  trim_fill <- trimfill(model)
  res_tf <- data.frame(Estimate = trim_fill$TE.random, ci.lb = trim_fill$lower.random, ci.ub = trim_fill$upper.random,
                       pi.lb = trim_fill$lower.predict, pi.ub = trim_fill$upper.predict)%>%
    dplyr::mutate_if(is.numeric, round, 3)
  
  heter_tf <- data.frame(
    tau2 = trim_fill$tau2, tau2.ci.lb = trim_fill$lower.tau2, tau2.ci.ub = trim_fill$upper.tau2,
    I2 = trim_fill$I2, I2.ci.lb = trim_fill$lower.I2, I2.ci.ub = trim_fill$upper.I2,
    H2 = trim_fill$H, H2.ci.lb = trim_fill$lower.H, H2.ci.ub = trim_fill$upper.H)%>%
    dplyr::mutate_if(is.numeric, round, 3)
  names(heter_tf) <- names(heter)[-c(1)]
  
  res_trim_fill <- data.frame(Model = c("tf with outliers", "model with outliers", "model without outliers"), 
                              rbind(res_tf,  res, res_out))
  heter_trim_fill <- data.frame(Model = c("tf with outliers", "model with outliers", "model without outliers"), 
                                rbind(heter_tf, heter[, -c(1)]))
  
  trim_fill <- list(trim_fill, res_trim_fill, heter_trim_fill)
  names(trim_fill) <- c("Model_tf",  "Results_tf", "Results_tf_out")
  #----- SELECTION MODEL  ------------------------------------------------------------------------
  set.seed(666)
  selection_model <- selmodel(model_metafor, type = "stepfun", steps = 0.025)
  test_sel <- paste0("Test for Selection Model Parameters: LRT(df = ", selection_model$LRTdf, ") = ", 
                     round(selection_model$LRT, 3), ", p-value = ", round(selection_model$LRTp, 4))
  res_sel <- data.frame(Estimate = selection_model$beta, ci.lb = selection_model$ci.lb, ci.ub = selection_model$ci.ub)%>%
    dplyr::mutate_if(is.numeric, exp)%>%
    dplyr::mutate_if(is.numeric, round, 4)
  res_sel <- paste0(res_sel$Estimate, " (95%CI = ", res_sel$ci.lb, " - ", res_sel$ci.ub, ")")
  p_sel <- cbind(selection_model$ptable, 
                 estimated_delta = selection_model$delta, 
                 ci.lb = selection_model$ci.lb.delta, ci.ub = selection_model$ci.ub.delta, 
                 zval_delta = selection_model$zval.delta, pval_delta = selection_model$pval.delta)%>%
    dplyr::mutate_if(is.numeric, round, 3)
  
  sel_model <- list(selection_model, test_sel, res_sel, p_sel)
  names(sel_model) <- c("model", "test", "estimates", "sel_model_results")
  #----- INFLUENCE  ------------------------------------------------------------------------
  # infl <- metainf(model, pooled = "random")
  # inf.an <- dmetar::InfluenceAnalysis(model, return.separate.plots = TRUE, random = TRUE)
  # tabla_influencia <- inf.an$Data%>%
  #   dplyr::mutate_if(is.numeric, round, 3)
  # names(tabla_influencia)[1] <- "Omitted"
  # bjt <- tabla_influencia%>%
  #   dplyr::mutate(studlab = gsub('[Omitting ]','', tabla_influencia$Omitted))
  # baujat_plot <- ggplot(bjt, aes(x = HetContrib, y = InfluenceEffectSize)) + geom_point(aes(size = (1/se)), color = "blue", alpha = 0.75) +
  #   geom_rug(color = "lightgray", alpha = 0.5) + theme(legend.position = "none") + xlab("Overall heterogeneity contribution") +
  #   ylab("Influence on pooled result") + 
  #   geom_label_repel(label = bjt$studlab, color = "black", size = 4, alpha = 0.7, 
  #                    box.padding = 0.5, label.padding = 0.5, alpha = 0.7, max.overlaps = 30)
  # 
  # res_influence <- list(infl, inf.an, inf.metafor, tabla_influencia, baujat_plot)
  # names(res_influence) <- c("mod_infl", "mod_infl_dmetar", "mod_infl_metafor", "tabla", "baujat")
  
  res <- list(model, outliers, bias_res, trim_fill, sel_model) #, res_influence)
  names(res) <- c("model", "outliers", "bias", "trim_fill", "selection_model") # "res_influence")
  return(res)
}

results <- lapply(allMeta_data_esc_def, models_results_meta)
saveRDS(results, file = "4_resultsRandomEffectMetanalysis.rds")



















