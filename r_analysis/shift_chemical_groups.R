library("dplyr")
library("car")
library(MASS)
library("ggplot2")




df <- read.csv("../data/tp53_synthetic_mutants_embedding_metrics_annotated.csv")



metric <- "k5_cosine_distance_to_WT"


physicochemical_features <- c(
  "charge_change",
  "class_change",
  "size_change",
  "introduces_proline",
  "introduces_glycine",
  "introduces_cysteine",
  "removes_proline",
  "removes_glycine",
  "removes_cysteine",
  "aromatic_change",
  "hydro_group_change",
  "is_conservative_like"
)


df <- df %>%
  mutate(
    across(
      all_of(physicochemical_features),
      ~ . == "True"
    )
  )


results <- lapply(physicochemical_features, function(feature){
  
  x <- df %>%
    filter(.data[[feature]] == TRUE) %>%
    pull(all_of(metric))
  
  y <- df %>%
    filter(.data[[feature]] == FALSE) %>%
    pull(all_of(metric))
  
  mw <- wilcox.test(
    x,
    y,
    exact = FALSE
  )
  
  cliff <- cliff.delta(x, y)
  
  data.frame(
    feature = feature,
    n_true = length(x),
    n_false = length(y),
    median_true = median(x),
    median_false = median(y),
    W = as.numeric(mw$statistic),
    p_value = mw$p.value,
    cliff_delta = cliff$estimate,
    cliff_magnitude = cliff$magnitude
  )
  
}) %>%
  bind_rows()


results <- results %>%
  mutate(
    p_adjusted = p.adjust(p_value, method = "holm")
  )

results <- results %>%
  mutate(
    direction = case_when(
      cliff_delta > 0 ~ "TRUE > FALSE",
      cliff_delta < 0 ~ "TRUE < FALSE",
      TRUE ~ "equal"
    )
  )

results %>%
  arrange(p_adjusted)

