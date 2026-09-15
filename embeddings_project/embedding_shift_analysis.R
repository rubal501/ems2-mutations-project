library("dplyr")
library("car")
library(MASS)
library("ggplot2")




df <- read.csv("tp53_synthetic_mutants_embedding_metrics_annotated.csv")

group_tags<- unique(df$group)


## Checamos la normalidad de los grupos 
df %>%
  group_by(group) %>%
  summarise(
    p_value = shapiro.test(k5_cosine_distance_to_WT)$p.value
  )


ggplot(df, aes(sample = k5_cosine_distance_to_WT)) +
  stat_qq() +
  stat_qq_line() +
  facet_wrap(~ group)

ggplot(df, aes(x = k5_cosine_distance_to_WT)) +
  geom_histogram()+
  facet_wrap(~ group)

ggplot(df, aes(x = k5_cosine_distance_to_WT,y=group)) +
  geom_boxplot()+
  labs(x="Distancia coseno K5 a WT", y="grupo de mutación")


## Pruebas de hipotesis


#Hacemos anova sin transformacion 

# primero hacemos anova
res.aov <- aov(lm(k5_cosine_distance_to_WT ~ group, data = df))

res.aov %>% summary()

leveneTest(k5_cosine_distance_to_WT ~ group, data=df)

aov_residuals <- residuals(object = res.aov)
shapiro.test(x = aov_residuals)

plot(res.aov,1)

plot(res.aov,2)



df <- df %>%
  mutate(log_k5 = log(k5_cosine_distance_to_WT))


df %>%
  group_by(group) %>%
  summarise(
    p_value = shapiro.test(log_k5)$p.value
  )



anova(lm(log_k5 ~ group, data = df))


kruskal.test(k5_cosine_distance_to_WT ~ group, data = df)

ggplot(df, aes(x = group, y = k5_cosine_distance_to_WT)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.15) +
  theme_bw() +
  labs(
    x = "Grupo",
    y = "Distancia coseno local k5"
  )

df %>%
  group_by(group) %>%
  summarise(
    n = n(),
    mediana = median(k5_cosine_distance_to_WT),
    media = mean(k5_cosine_distance_to_WT),
    sd = sd(k5_cosine_distance_to_WT),
    iqr = IQR(k5_cosine_distance_to_WT)
  )

kruskal.test(k5_cosine_distance_to_WT ~ group, data = df)

pairwise.wilcox.test(
  df$k5_cosine_distance_to_WT,
  df$group,
  p.adjust.method = "holm"
)

library(effsize)

grupos <- unique(df$group)

combn(grupos, 2, function(g) {
  sub_df <- df %>% filter(group %in% g)
  
  test <- cliff.delta(
    k5_cosine_distance_to_WT ~ group,
    data = sub_df
  )
  
  data.frame(
    grupo_1 = g[1],
    grupo_2 = g[2],
    cliff_delta = test$estimate,
    magnitude = test$magnitude
  )
}, simplify = FALSE) %>%
  bind_rows() -> cliff_df



comparaciones <- combn(unique(df$group), 2, simplify = FALSE)


resultados <- lapply(comparaciones, function(gr) {
  
  x <- df %>%
    filter(group == gr[1]) %>%
    pull(k5_cosine_distance_to_WT)
  
  y <- df %>%
    filter(group == gr[2]) %>%
    pull(k5_cosine_distance_to_WT)
  
  test <- wilcox.test(x, y)
  
  data.frame(
    grupo_1 = gr[1],
    grupo_2 = gr[2],
    W = test$statistic,
    p_value = test$p.value
  )
})

bind_rows(resultados)

resultados

metric <- "k5_cosine_distance_to_WT"


# Kruskal-Wallis global
kruskal_result <- kruskal.test(
  formula = as.formula(paste(metric, "~ group")),
  data = df
)

kruskal_result

# Mann-Whitney / Wilcoxon por pares
pairwise_result <- pairwise.wilcox.test(
  x = df[[metric]],
  g = df$group,
  p.adjust.method = "holm",
  exact = FALSE
)

pairwise_result

groups <- unique(df$group)

## MW test por grupo de mutacion

mann_whitney_table <- combn(groups, 2, simplify = FALSE) %>%
  lapply(function(pair) {
    
    x <- df %>%
      filter(group == pair[1]) %>%
      pull(all_of(metric))
    
    y <- df %>%
      filter(group == pair[2]) %>%
      pull(all_of(metric))
    
    test <- wilcox.test(x, y, exact = FALSE)
    
    data.frame(
      group_1 = pair[1],
      group_2 = pair[2],
      W = as.numeric(test$statistic),
      p_value = test$p.value
    )
  }) %>%
  bind_rows() %>%
  mutate(
    p_adjusted_holm = p.adjust(p_value, method = "holm")
  )

mann_whitney_table

mann_whitney_table <- mann_whitney_table %>%
  mutate(
    significance = case_when(
      p_adjusted_holm < 0.001 ~ "***",
      p_adjusted_holm < 0.01  ~ "**",
      p_adjusted_holm < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )

mann_whitney_table


## Temporal 



metric <- "k5_cosine_distance_to_WT"

# Pruebas para uno contra todos
results <- lapply(unique(df$group), function(target_group) {
  
  x <- df %>%
    filter(group == target_group) %>%
    pull(all_of(metric))
  
  y <- df %>%
    filter(group != target_group) %>%
    pull(all_of(metric))
  
  mw <- wilcox.test(
    x,
    y,
    exact = FALSE
  )
  
  cliff <- cliff.delta(x, y)
  
  data.frame(
    group = target_group,
    n_group = length(x),
    n_rest = length(y),
    median_group = median(x),
    median_rest = median(y),
    W = as.numeric(mw$statistic),
    p_value = mw$p.value,
    cliff_delta = cliff$estimate,
    cliff_magnitude = cliff$magnitude
  )
})

results <- bind_rows(results) %>%
  mutate(
    p_adjusted = p.adjust(p_value, method = "holm")
  )

results



metric <- "k5_cosine_distance_to_WT"

group_pairs <- combn(unique(df$group), 2, simplify = FALSE)

pairwise_results <- lapply(group_pairs, function(pair) {
  
  g1 <- pair[1]
  g2 <- pair[2]
  
  x <- df %>%
    filter(group == g1) %>%
    pull(all_of(metric))
  
  y <- df %>%
    filter(group == g2) %>%
    pull(all_of(metric))
  
  mw <- wilcox.test(
    x,
    y,
    exact = FALSE
  )
  
  cliff <- cliff.delta(x, y)
  
  data.frame(
    group_1 = g1,
    group_2 = g2,
    n_1 = length(x),
    n_2 = length(y),
    median_1 = median(x),
    median_2 = median(y),
    W = as.numeric(mw$statistic),
    p_value = mw$p.value,
    cliff_delta = cliff$estimate,
    cliff_magnitude = cliff$magnitude
  )
})

## Se revisa el efecto de la 
pairwise_results <- bind_rows(pairwise_results) %>%
  mutate(
    p_adjusted = p.adjust(p_value, method = "holm"),
    direction = case_when(
      cliff_delta > 0 ~ paste(group_1, ">", group_2),
      cliff_delta < 0 ~ paste(group_1, "<", group_2),
      TRUE ~ "equal"
    )
  )

print("Comparaciones a pares con MW test (grupos de mutaciones)")
pairwise_results







