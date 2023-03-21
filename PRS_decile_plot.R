# Required packages
library(dplyr)
library(ggplot2)
library(ggthemes)

# Convert variables to factors and get quantiles
prs_table <- train %>%
  mutate(
    decile1 = ntile(train$PRS1, 10),
    decile2 = ntile(train$PRS2, 10),
    phenotype = as.factor(train$PHENOTYPE)
  )

prs_table$decile1 <- as.factor(prs_table$decile1)
prs_table$decile2 <- as.factor(prs_table$decile2)

# Fit regression model
prs_glm1 <- glm(phenotype ~ decile1,
               data = prs_table,
               family = 'binomial')
prs_glm2 <- glm(phenotype ~ decile2,
                data = prs_table,
                family = 'binomial')

# Put results in data.frame
summs1 <- prs_glm1 %>% summary()
summs2 <- prs_glm2 %>% summary()

# Get point estimates and SEs
results1 <- bind_cols(coef(prs_glm1),
                     summs1$coefficients[, 2]) %>%
  setNames(c("estimate", "se"))  %>%
  mutate(decile = 1:10)

results2 <- bind_cols(coef(prs_glm2),
                      summs2$coefficients[, 2]) %>%
  setNames(c("estimate", "se"))  %>%
  mutate(decile = 1:10)
# Your coefficients are on the log odds scale, such that a coefficient is
# log(odds_y==1 / odds_y == 0). We can exponentiate to get odds instead.
results_odds1 <- results1 %>% mutate(across(.cols = -decile,
                                          ~ exp(.x)))
results_odds2 <- results2 %>% mutate(across(.cols = -decile,
                                            ~ exp(.x)))

# Need SEs on the odds scale too
results_odds1 <- results_odds1 %>%
  mutate(var_diag = diag(vcov(prs_glm1)),
         se = sqrt(estimate ^ 2 * var_diag))

results_odds2 <- results_odds2 %>%
  mutate(var_diag = diag(vcov(prs_glm2)),
         se = sqrt(estimate ^ 2 * var_diag))

results_odds1$PRS <- 1
results_odds2$PRS <- 2

results_odds <- rbind(results_odds1,results_odds2)
# Plot with +/- 1 SE
p <- ggplot(results_odds, aes(x = as.factor(decile), y = estimate, color=factor(PRS),group=factor(PRS))) +
  geom_point(stat = "identity", shape = 15, position = position_dodge(width = 0.40)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey") +
  geom_errorbar(aes(ymin = estimate - se, ymax = estimate + se), width = 0.4,
                position="dodge") +
  xlab("PRS Decile") +
  ylab("Odds Ratio (95%CI)")

p + theme_hc() +
  scale_color_hue(name = "Polygenic Risk Score",labels = c("PRS1", "PRS2"))
