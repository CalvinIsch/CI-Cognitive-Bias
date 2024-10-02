library(dplyr); library(tidyr); library(DescTools); library(igraph)

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



# Cor guess and neighborhood
df <- df %>% filter(Net_type == "Network")
df$Dist_Neighbor = abs(as.numeric(df$response_1) - as.numeric(df$neighbors_1))
df$neighbor_confidence = factor(df$neighbor_confidence, 
                         levels = c("Not at all confident", "Not too confident", 
                                    "Somewhat confident", "Very confident"))

cor.test(df$Dist_Neighbor, as.numeric(df$neighbor_confidence), method = "spearman")

df$Accuracy_1 <- abs(df$response_1 - 246)
df$initial_confidence = factor(df$initial_confidence, 
                                levels = c("Not at all confident", "Not too confident", 
                                           "Somewhat confident", "Very confident"))
cor.test(df$Accuracy_1, as.numeric(df$initial_confidence), method = "spearman")



############### Degroot Simulations
library(ggplot2); library(ggpubr); library(igraph); library(viridis);
library(tidyverse); library(sandwich); library(ggpattern); library(sf); 
library(DescTools); library(network); library(reshape2); library(boot); library(bootnet))

set.seed(42)

# Generate a normal distribution based on mean and sd of our data
df = read.csv('ci_cb_2023-06-26.csv')
# Winsorize the response_1
df$response_1 <- Winsorize(df$response_1, probs = c(0.02, 0.98), na.rm=T)
mean_value <- mean(df %>% pull(response_1), na.rm = TRUE)
sd_value <- sd(df %>% pull(response_1), na.rm = TRUE)

get_value <- function(){
  value <- rnorm(1, mean = mean_value, sd = sd_value)
  if(value < 0){
    value <- 0
  }
  return(value)
}

# Simulate degroot model - Egalitarian network
df <- data.frame()
sim_degroot <- function(n) {
  # 30 initial guesses
  guess_1 <- sapply(1:30, function(x) get_value())
  
  # Generate a correlation between self-weights and errors ~.6
  errors_1 <- abs(guess_1 - 35)
  self_weights = -errors_1*.5 + rnorm(30, 0, 1000) + rnorm(30, 0, 1000) 
  
  # Standardize self_weights between 0 and 1
  self_weights <- (self_weights - min(self_weights)) / (max(self_weights) - min(self_weights))
  
  neighborhood_means <- c()
  for (i in 1:30){
    node <- egal$Node[i]
    neighbors <- c(egal$N1[i], egal$N2[i], egal$N3[i], egal$N4[i])
    neighborhood_means <- c(neighborhood_means, mean(guess_1[neighbors]))
  }
  guess_2 <- self_weights * guess_1 + (1 - self_weights) * neighborhood_means
  
  neighborhood_means_2 <- c()
  for (i in 1:30){
    node <- egal$Node[i]
    neighbors <- c(egal$N1[i], egal$N2[i], egal$N3[i], egal$N4[i])
    neighborhood_means_2 <- c(neighborhood_means_2, mean(guess_2[neighbors]))
  }
  guess_3 <- self_weights * guess_2 + (1 - self_weights) * neighborhood_means_2
  df_t <- data.frame(X = 1:30, Change_1 = (guess_2 - 246) - (guess_1 - 246), Change_2 = (guess_3 - 246) - (guess_1 - 246))
  df_t$Instance <- n
  return(df_t)
}
egal <- read.csv('neighbors.csv')
for (i in 1:1000){
  d_t <- sim_degroot(i)
  df <- rbind(df, d_t)
}

df <- df %>% group_by(X) %>% summarise(Change_1 =mean(Change_1), Change_2 = mean(Change_2))


### Make nice layout

egal_g <- make_empty_graph(n = nrow(egal), directed = TRUE)

for(i in 1:nrow(egal)){
  node <- egal$Node[i]
  neighbors <- c(egal$N1[i], egal$N2[i], egal$N3[i], egal$N4[i])
  for(neigh in neighbors){
    egal_g <- add_edges(egal_g, c(node, neigh))
  }
}
egal_g <- as.undirected(egal_g, mode="collapse")
N <- vcount(egal_g)
half_N <- N/2

# Define radii for inner and outer circles
r_inner <- 0.8
r_outer <- 1
angles <- seq(0, 2*pi, length.out = half_N + 1)
angles <- angles[-(half_N + 1)]  
x_inner <- r_inner * cos(angles)
y_inner <- r_inner * sin(angles)
x_outer <- r_outer * cos(angles)
y_outer <- r_outer * sin(angles)
x_coords <- c(x_inner, x_outer)
y_coords <- c(y_inner, y_outer)

label_size = 4
vertex_size = 8

df_network = df
# Get nice color gradient
color_gradient <- colorRampPalette(c("darkgreen",'green',"lightgreen","yellow", "orange","darkorange","red"))
i <- 60
values <- seq(-i, i, length.out = i*2)
colors <- color_gradient(length(values))
get_color <- function(value) {
  if (value < -i){
    value <- -i
  } else if (value > i-1){
    value <- i-1
  }
  scaled_value <- as.integer(round(value,0)) + i+ 1
  colors[scaled_value]
}


png(filename="Figures/degroot-sims-egal.png", width=2400, height=600)
par(mfrow=c(1,3), mar=c(5, 5, 5, 5), mgp=c(3, 1, 0))

# Round 1
example_colors <- rep("gray", 30)
plot(egal_g,  layout=matrix(c(x_coords, y_coords), ncol=2), vertex.label=NA, vertex.size=vertex_size, vertex.color=example_colors)
title(main="Round 1", cex.main=label_size, font.main=1.5)

# Round 2
guesses <- c(df_network %>% pull(Change_1), rep(NA, 30 - length(df_network$X)))
example_colors <- sapply(guesses, get_color)
plot(egal_g,  layout=matrix(c(x_coords, y_coords), ncol=2), vertex.label=NA, vertex.size=vertex_size, vertex.color=example_colors)
title(main="Round 2", cex.main=label_size, font.main=1.5)

# Round 3
guesses <- c(df_network %>% pull(Change_2), rep(NA, 30 - length(df_network$X)))
example_colors <- sapply(guesses, get_color)
plot(egal_g,  layout=matrix(c(x_coords, y_coords), ncol=2), vertex.label=NA, vertex.size=vertex_size, vertex.color=example_colors)
title(main="Round 3", cex.main=label_size, font.main=1.5)

dev.off()

df_ind <- df


##### Individuals ####
df <- data.frame()
sim_degroot_individual <- function(n) {
  # 30 initial guesses
  guess_1 <- sapply(1:30, function(x) get_value())
  
  # No neighbors, just add random noise
  guess_2 <- guess_1 + rnorm(30, 0, 400)
  guess_3 <- guess_2 + rnorm(30, 0, 400)
  df_t <- data.frame(X = 1:30, Change_1 = (guess_2 - 246) - (guess_1 - 246), Change_2 = (guess_3 - 246) - (guess_1 - 246))
  df_t$Instance <- n
  return(df_t)
}

for (i in 1:1000){
  d_t <- sim_degroot_individual(i)
  df <- rbind(df, d_t)
}

df <- df %>% group_by(X) %>% summarise(Change_1 =mean(Change_1), Change_2 = mean(Change_2))
df_network <- df
nodes_only_graph <- make_empty_graph(n = vcount(egal_g))

color_gradient <- colorRampPalette(c("darkgreen",'green',"lightgreen","yellow", "orange","darkorange","red"))
i <- 60
values <- seq(-i, i, length.out = i*2)
colors <- color_gradient(length(values))
get_color <- function(value) {
  if (value < -i){
    value <- -i
  } else if (value > i-1){
    value <- i-1
  }
  scaled_value <- as.integer(round(value,0)) + i+ 1
  colors[scaled_value]
}




png(filename="Figures/degroot-sims-individuals.png", width=2400, height=600)
par(mfrow=c(1,3), mar=c(5, 5, 5, 5), mgp=c(3, 1, 0))

# Round 1
example_colors <- rep("gray", 30)
plot(nodes_only_graph,  layout=matrix(c(x_coords, y_coords), ncol=2), vertex.label=NA, vertex.size=vertex_size, vertex.color=example_colors)
title(main="Round 1", cex.main=label_size, font.main=1.5)

# Round 2
guesses <- c(df_network %>% pull(Change_1), rep(NA, 30 - length(df_network$X)))
example_colors <- sapply(guesses, get_color)
plot(nodes_only_graph,  layout=matrix(c(x_coords, y_coords), ncol=2), vertex.label=NA, vertex.size=vertex_size, vertex.color=example_colors)
title(main="Round 2", cex.main=label_size, font.main=1.5)

# Round 3
guesses <- c(df_network %>% pull(Change_2), rep(NA, 30 - length(df_network$X)))
example_colors <- sapply(guesses, get_color)
plot(nodes_only_graph,  layout=matrix(c(x_coords, y_coords), ncol=2), vertex.label=NA, vertex.size=vertex_size, vertex.color=example_colors)
title(main="Round 3", cex.main=label_size, font.main=1.5)

dev.off()








# More simulations of size 50 per Damon's request

# Fully connected
df <- data.frame()
sim_fc_degroot <- function(n) {
  guess_1 <- rnorm(50, mean = mean_value, sd = sd_value)
  
  # Generate self-weights that are anticorrelated with error and between 0 and 1
  errors_1 <- abs(guess_1 - 246)
  self_weights = -errors_1*.5 + rnorm(50, 0, 1000)  
  self_weights <- (self_weights - min(self_weights)) / (max(self_weights) - min(self_weights))
  
  guess_2 <- self_weights * guess_1 + (1 - self_weights) * mean(guess_1)
  guess_3 <- self_weights * guess_2 + (1 - self_weights) * mean(guess_2)
  
  df_t <- data.frame(X = 1:50, Change_1 = abs((guess_2 - 246)) - errors_1, Change_2 = abs((guess_3 - 246)) - errors_1)
  df_t$Instance <- n
  return(df_t)
}

for (i in 1:1000){
  d_t <- sim_fc_degroot(i)
  df <- rbind(df, d_t)
}
# To show nodes
#df_fc <- df %>% group_by(X) %>% summarise(Change_1 =mean(Change_1), Change_2 = mean(Change_2))
# To show networks
df_fc <- df %>% group_by(Instance) %>% summarize(Change_1 = mean(Change_1), Change_2 = mean(Change_2))


# Random network - Degree 4
df <- data.frame()
g <- sample_k_regular(no.of.nodes = 50, k = 4)
sim_random_degroot <- function(n) {
  guess_1 <- rnorm(50, mean = mean_value, sd = sd_value)
  
  # Generate self-weights that are anticorrelated with error and between 0 and 1
  errors_1 <- abs(guess_1 - 246)
  self_weights = -errors_1*.5 + rnorm(50, 0, 1000)  
  self_weights <- (self_weights - min(self_weights)) / (max(self_weights) - min(self_weights))
  
  neighborhood_means <- c()
  for (i in 1:50){
    neighborhood_means <- c(neighborhood_means, mean(guess_1[g[[i]][[1]]]))
  }
  guess_2 <- self_weights * guess_1 + (1 - self_weights) * neighborhood_means
  
  neighborhood_means_2 <- c()
  for (i in 1:50){
    neighborhood_means_2 <- c(neighborhood_means_2, mean(guess_2[g[[i]][[1]]]))
  }
  guess_3 <- self_weights * guess_2 + (1 - self_weights) * neighborhood_means_2
  
  df_t <- data.frame(X = 1:50, Change_1 = abs((guess_2 - 246)) - errors_1, Change_2 = abs((guess_3 - 246)) - errors_1)
  df_t$Instance <- n
  return(df_t)
}

for (i in 1:1000){
  d_t <- sim_random_degroot(i)
  df <- rbind(df, d_t)
}
# To show mean change for each node 
#df_ran4 <- df %>% group_by(X) %>% summarise(Change_1 =mean(Change_1), Change_2 = mean(Change_2))
# To show mean change for each network
df_ran4 <- df %>% group_by(Instance) %>% summarize(Change_1 = mean(Change_1), Change_2 = mean(Change_2))

g <- sample_k_regular(no.of.nodes = 50, k = 6)
for (i in 1:1000){
  d_t <- sim_random_degroot(i)
  df <- rbind(df, d_t)
}
df_ran6 <- df %>% group_by(Instance) %>% summarize(Change_1 = mean(Change_1), Change_2 = mean(Change_2))




create_hexagonal_grid <- function(rows, cols) {
  g <- make_empty_graph(n = 0, directed = FALSE)
  for (i in 0:(rows - 1)) {
    for (j in 0:(cols - 1)) {
      node_id <- i * cols + j + 1
      g <- add_vertices(g, 1)
      if (j > 0) { # connect to the left neighbor
        g <- add_edges(g, c(node_id, node_id - 1))
      }
      if (i > 0) { # connect to the top neighbor
        g <- add_edges(g, c(node_id, node_id - cols))
        if (j > 0 && (i %% 2 == 0)) { # connect to the top-left neighbor on even rows
          g <- add_edges(g, c(node_id, node_id - cols - 1))
        }
        if (j < cols - 1 && (i %% 2 == 1)) { # connect to the top-right neighbor on odd rows
          g <- add_edges(g, c(node_id, node_id - cols + 1))
        }
      }
    }
  }
  
  # Add extra connections to the edge
  
  # From left to right
  for (i in 0:(rows - 1)) {
    g <- add_edges(g, c(i*(cols)+1, (i+1)*(cols)))
    print(paste(i*(cols)+1, (i+1)*(cols)))
  }
  # From top to bottom
  for (i in 0:(cols - 1)) {
    g <- add_edges(g, c(i+1, (rows-1)*(cols)+i+1))
  }
  # Diagnol from top to bottom
  for (i in 1:(cols - 1)) {
    g <- add_edges(g, c(i+1, (rows-1)*(cols)+i))
  }
  g <- add_edges(g, c(1, 14))
  g <- add_edges(g, c(1, 49))
  g <- add_edges(g, c(42, 43))
  g <- add_edges(g, c(42, 29))
  g <- add_edges(g, c(28, 29))
  g <- add_edges(g, c(28, 15))
  g <- add_edges(g, c(15, 14))
  
  return(g)
}

# Create the hexagonal grid graph
rows <- 7
cols <- 7
g <- create_hexagonal_grid(rows, cols)



sim_random_hex<- function(n) {
  # 49 bc needs square
  guess_1 <- rnorm(49, mean = mean_value, sd = sd_value)
  
  # Generate self-weights that are anticorrelated with error and between 0 and 1
  errors_1 <- abs(guess_1 - 246)
  self_weights = -errors_1*.5 + rnorm(49, 0, 1000)  
  self_weights <- (self_weights - min(self_weights)) / (max(self_weights) - min(self_weights))
  
  neighborhood_means <- c()
  for (i in 1:49){
    neighborhood_means <- c(neighborhood_means, mean(guess_1[g[[i]][[1]]]))
  }
  guess_2 <- self_weights * guess_1 + (1 - self_weights) * neighborhood_means
  
  neighborhood_means_2 <- c()
  for (i in 1:49){
    neighborhood_means_2 <- c(neighborhood_means_2, mean(guess_2[g[[i]][[1]]]))
  }
  guess_3 <- self_weights * guess_2 + (1 - self_weights) * neighborhood_means_2
  
  df_t <- data.frame(X = 1:49, Change_1 = abs((guess_2 - 246)) - errors_1, Change_2 = abs((guess_3 - 246)) - errors_1)
  df_t$Instance <- n
  return(df_t)
}



df <- data.frame()
for (i in 1:1000){
  d_t <- sim_random_hex(i)
  df <- rbind(df, d_t)
}
df_hex6 <- df %>% group_by(Instance) %>% summarize(Change_1 = mean(Change_1), Change_2 = mean(Change_2))

mean(df_hex6$Change_2)
sd(df_hex6$Change_2) / sqrt(1000)





# Plot individual for comparison
df <- data.frame()
sim_degroot_individual <- function(n) {
  guess_1 <- sapply(1:50, function(x) get_value())
  guess_2 <- guess_1 + rnorm(50, 0, 100)
  guess_3 <- guess_2 + rnorm(50, 0, 100)
  df_t <- data.frame(X = 1:50, Change_1 = (guess_2 - 246) - (guess_1 - 246), Change_2 = (guess_3 - 246) - (guess_1 - 246))
  df_t$Instance <- n
  return(df_t)
}

for (i in 1:1000){
  d_t <- sim_degroot_individual(i)
  df <- rbind(df, d_t)
}
df_ind <- df %>% group_by(X) %>% summarise(Change_1 =mean(Change_1), Change_2 = mean(Change_2))





df <-data.frame()
sim_degroot_centralized <- function(n) {
  guess_1 <- sapply(1:50, function(x) get_value())
  # Generate self-weights that are anticorrelated with error and between 0 and 1
  errors_1 <- abs(guess_1 - 246)
  self_weights = -errors_1*.5 + rnorm(50, 0, 1000)  
  self_weights <- (self_weights - min(self_weights)) / (max(self_weights) - min(self_weights))
  
  # Only central agent matters
  guess_2 <- self_weights * guess_1 + (1 - self_weights) * guess_1[1]
  guess_2[1] <- self_weights[1] * guess_1[1] + (1 - self_weights[1]) * mean(guess_1[-1])
  
  guess_3 <- self_weights * guess_2 + (1 - self_weights) * guess_2[1]
  guess_3[1] <- self_weights[1] * guess_2[1] + (1 - self_weights[1]) * mean(guess_2[-1])
  df_t <- data.frame(X = 1:50, Change_1 = (guess_2 - 246) - (guess_1 - 246), Change_2 = (guess_3 - 246) - (guess_1 - 246))
  df_t$Instance <- n
  return(df_t)
}

for (i in 1:1000){
  d_t <- sim_degroot_centralized(i)
  df <- rbind(df, d_t)
}
df_central <- df %>% group_by(Instance) %>% summarize(Change_1 = mean(Change_1), Change_2 = mean(Change_2))
mean(df_central$Change_2)

# Get mean and se for each tyoe
networks <- c("FC", "Random 4", "Random 6", "Hexagonal 6", "Individuals", "Centralized")
means <- c(mean(df_fc$Change_2), mean(df_ran4$Change_2), mean(df_ran6$Change_2), 
           mean(df_hex6$Change_2), mean(df_ind$Change_2), mean(df_central$Change_2))
ses <- c(sd(df_fc$Change_2) / sqrt(1000), sd(df_ran4$Change_2) / sqrt(1000), 
         sd(df_ran6$Change_2) / sqrt(1000), sd(df_hex6$Change_2) / sqrt(1000),
         sd(df_ind$Change_2) / sqrt(1000), sd(df_central$Change_2) / sqrt(1000))
df_plot <- data.frame(networks, means, ses)
colnames(df_plot) <- c("Network", "means", "ses")

df_plot$Network <- factor(df_plot$Network, levels = c("FC", "Hexagonal 6", "Random 6", "Random 4", "Individuals", "Centralized"))





ggplot(df_plot, aes(x = Network, y = means, fill = Network)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = means - ses, ymax = means + ses), width = 0.2) +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +
  labs(title = "1,000 Degroot simulations w/Rev. Coef", x = "Network Type", y = "Mean Change in Error")






