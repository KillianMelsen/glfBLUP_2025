#` Adjusted values of a phenotype based on genotypic effects and residuals
#`
#` @description `BLUE_adjust` returns a vector the same length as the phenotype containing 
#` BLUP adjusted values, with BLUPs coming from a full random-effects model.
#` @param model A linear-mixed model with fixed genotypic effect of class lmerMod
#` @param gid the variable name (character) of the genotypic effect
#` @return a vector of the same length as the dependent variable in the model
#` @examples adjust_values(lattice.lmer, 'GID')
BLUE_adjust <- function(model, gid, alpha_scaling = TRUE){
  # older style, when model without intercept is used
  # temp.residuals <- residuals(model)
  # fix_effects <- lme4::fixef(model)
  # names(fix_effects) <- gsub(gid, '', names(fix_effects))
  # matched.gid.blues <- fix_effects[match(model@frame[[gid]], names(fix_effects))]
  # adjusted_output <- matched.gid.blues + temp.residuals
  
  temp.residuals <- residuals(model)
  fix_effects <- lme4::fixef(model)
  names(fix_effects) <- levels(model@frame[[gid]])
  fix_effects <- c(fix_effects[1],(fix_effects[2:length(fix_effects)] + fix_effects[1]))
  matched.gid.blues <- fix_effects[match(model@frame[[gid]], names(fix_effects))]
  adjusted_output <- matched.gid.blues + temp.residuals
  
  if(alpha_scaling){
    ve <- as.data.frame(lme4::VarCorr(model))$vcov[length(as.data.frame(lme4::VarCorr(model))$vcov)]
    n.rep.vector <- as.integer(table(model@frame[[gid]]))
    replicates <- (sum(n.rep.vector) - sum(n.rep.vector^2)/sum(n.rep.vector))/(length(n.rep.vector) - 1)
    vg_anova <- (anova(model)$`Mean Sq` - ve) / replicates
    constant <- vg_anova / ve
    pseudo_df <- data.frame(
      gid = names(adjusted_output),
      vals = as.numeric(adjusted_output)
    )
    pseudo.aov <- anova(lm(vals ~ gid, pseudo_df))
    msg <- pseudo.aov$`Mean Sq`[1]
    mse <- pseudo.aov$`Mean Sq`[2]
    alpha <- (msg / mse) * ( 1 / (replicates * constant + 1))
    adjusted_output <- matched.gid.blues + (temp.residuals * sqrt(alpha))
    
    alpha_repeat <- heritability::repeatability(adjusted_output, pseudo_df$gid)
    
    if(round(alpha_repeat$repeatability, 4) != round(vg_anova / (vg_anova + ve), 4)){
      warning('Heritability alpha-scaled and estimated design model unequal')
    }
  }
  return(adjusted_output)
}
