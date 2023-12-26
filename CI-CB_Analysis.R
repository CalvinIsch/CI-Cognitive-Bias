library(dplyr); library(tidyr); library(DescTools)

# Read the CSV file
df <- read.csv('ci_cb_2023-06-26.csv')

# Data manipulation
df <- df %>%
  mutate(
    Anchor = ifelse(Condition == 'Low', 118, 353),
    Mag_Rev1 = abs(response_2 - response_1),
    Mag_Rev2 = abs(response_3 - response_2),
    Mag_Rev = abs(response_3 - response_1),
    Initial_Error = abs(response_1 - true_answer),
    Second_Error = abs(response_2 - true_answer),
    Final_Error = abs(response_3 - true_answer),
    Dist_Anchor_1 = abs(response_1 - Anchor),
    Dist_Anchor_2 = abs(response_2 - Anchor),
    Dist_Anchor_3 = abs(response_3 - Anchor),
    Difference_Initial_Anchor = response_1 - Anchor,
    Difference_Initial_Truth = response_1 - 246,
    Matched = ifelse((Difference_Initial_Anchor * Difference_Initial_Truth) > 0, 1, 0),
    Change_accuracy = Initial_Error - Final_Error
  )


#### JonckheereTerpstra
df_network <- df %>% filter(Net_type != 'Control') %>% drop_na(c('Mag_Rev'))
df_1 <- df_network %>% 
  select(group, node_id, Initial_Error, Mag_Rev) %>% drop_na()
df_1 <- df_1 %>% mutate(quantile = ntile(Initial_Error, 10))
q <- quantile(df_1$Initial_Error, probs = seq(0, 1, 1/10),na.rm=TRUE)
df_1$quantile <- factor(df_1$quantile)

JonckheereTerpstraTest(df_1$Mag_Rev,df_1$quantile)
