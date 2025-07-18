# Article title: "Biodiversity mediates the effects of socio-environmental factors on coastal fish biomass: 
#                 a structural equation modeling approach in three ecoregions" 
#
# Script authors: David Mouillot, Jean-Baptiste Juhel
# date: 2025-01-14

# This script provides Structural Equation Model (SEM) 

#### Libraries  & workspace ####
library(piecewiseSEM)
library(multcompView)
library(MuMIn)
library(semEff)
library(MASS)
library(car)
library(visreg)
library(DiagrammeR)
library(lavaan)
library(semPlot)
library(dplyr)
library(ggplot2)
library(patchwork)


setwd("C:/Users/jeanb/Documents/Project_MACOBIOS")

#### SEM function ####
# SEM
plot.psem <- function(x, return=FALSE,
                      node_attrs = data.frame(shape = "rectangle", color = "black",
                                              fillcolor = "white"),
                      edge_attrs = data.frame(style = "solid", color="black"),
                      ns_dashed = T, alpha=0.05,
                      show = "std", digits = 3, 
                      add_edge_label_spaces = TRUE, ...
){
  
  #get the coefficients table
  ctab <- coefs(x)
  ctab$Response <- as.character(ctab$Response)
  ctab$Predictor <- as.character(ctab$Predictor)
  
  #make a nodes DF
  unique_nodes <- unique(c(ctab$Response, ctab$Predictor))
  nodes <- create_node_df(n = length(unique_nodes),
                          nodes = unique_nodes,
                          type = "lower",
                          label = unique_nodes)
  nodes <- cbind(nodes, node_attrs)
  nodes[] <- lapply(nodes, as.character)
  nodes$id <- as.numeric(nodes$id)
  
  #make an edges DF
  edges <- create_edge_df(
    from = match(ctab$Predictor, unique_nodes),
    to = match(ctab$Response, unique_nodes))
  
  edges <- data.frame(edges, edge_attrs)
  edges[] <- lapply(edges, as.character)
  edges$id <- as.numeric(edges$id)
  edges$from <- as.numeric(edges$from)
  edges$to <- as.numeric(edges$to)
  if(ns_dashed) edges$style[which(ctab$P.Value>alpha)] <- "dashed"
  if(show == "std") edges$label = round(ctab$`Std.Estimate`, digits)
  if(show == "unstd") edges$label = round(ctab$Estimate, digits)
  if(add_edge_label_spaces) edges$label = paste0(" ", edges$label, " ")
  
  #turn into a graph
  sem_graph <- create_graph(nodes, edges, directed=TRUE)
  
  if(return) return(sem_graph)
  
  render_graph(sem_graph, ...)
  
}

#### France ####
data <- read.csv("data/France_data.csv",header=T,row.names=1)

# Transform data
hist(data$ES_totbiom_fished)
cor(log10(data$ES_totbiom+1),log10(data$ES_totbiom_fished+1))
logBiom <- log10(data$ES_totbiom+1)
data <- cbind(data,logBiom)
data <- data[data$BIO_qual>0.70,]
hist(data$logBiom,xlab ="log(Total Fish Biomass)", main=NULL)

# Variable selection
mod.Biom <- lm(logBiom ~ ENV_Depth+ENV_Hab_block+ENV_Hab_scramble+ENV_Hab_rock+ENV_Hab_sand+
                 BIO_Taxo_Ric+BIO_Taxo_Ent+BIO_Func_Ric+BIO_Func_Ent+BIO_Phyl_Ric+BIO_Phyl_Ent, data)
stepAIC(mod.Biom)

tax_ric <- lm(BIO_Taxo_Ric ~ENV_Depth+ENV_Hab_block+ENV_Hab_scramble+ENV_Hab_rock+ENV_Hab_sand, data)
stepAIC(tax_ric)

phy_ent <- lm(BIO_Phyl_Ent ~ ENV_Depth+ENV_Hab_block+ENV_Hab_scramble+ENV_Hab_rock+ENV_Hab_sand, data)
stepAIC(phy_ent)


# Parsimonious model
model.simp <- psem(
  lm(BIO_Taxo_Ric ~ ENV_Hab_sand, 
     data = data),
  lm(BIO_Phyl_Ent ~ ENV_Hab_sand, data = data),
  lm(logBiom ~ BIO_Taxo_Ric + BIO_Phyl_Ent, data = data)
)

anova(model.simp)
coefs(model.simp)
rsquared(model.simp)


plot(model.simp)
plot.psem(model.simp)
x <- data.frame(shape = "rectangle", color = "black",fillcolor = "orange",y=1:2)
y <- data.frame(shape = "rectangle", color = "black",fillcolor = "black",y=3)
z <- data.frame(shape = "rectangle", color = "black",fillcolor = "blue",y=4)
col <- rbind(x,y,z)
plot.psem(model.simp,node_attrs = col)
line <- data.frame(style = "solid", color = "black")
plot.psem(model.simp,node_attrs = col,edge_attrs=line) 

# Model residuals
# Objects of class psem do not have global residuals, they have several imbricated models
# Normality test (Shapiro-Wilk)
model.simp <- psem(
  lm(BIO_Taxo_Ric ~ ENV_Hab_sand, 
     data = data),
  lm(BIO_Phyl_Ent ~ ENV_Hab_sand, data = data),
  lm(logBiom ~ BIO_Taxo_Ric + BIO_Phyl_Ent, data = data)
)

mod_list <- list(
  TaxR_model =  lm(BIO_Taxo_Ric ~ ENV_Hab_sand, 
                   data = data),
  PhyE_model = lm(BIO_Phyl_Ent ~ ENV_Hab_sand, data = data),
  logBiom_model = lm(logBiom ~ BIO_Taxo_Ric + BIO_Phyl_Ent, data = data)
)
res_list <- lapply(mod_list, residuals)

resid_data <- bind_rows(
  lapply(names(mod_list), function(name) {
    res <- residuals(mod_list[[name]])
    pval <- shapiro.test(res)$p.value
    data.frame(
      residuals = res,
      model = name,
      shapiro_p = sprintf("p = %.3f", pval)
    )
  })
)

p_dens_france <- ggplot(resid_data, aes(x = residuals)) +
  geom_density(fill = "steelblue", alpha = 0.6) +
  facet_wrap(~model, ncol = 1, scales = "free") +
  theme_minimal() +
  labs(
    x = "Residuals",
    y = "Density",
    title = "France"
  )

p_combo_france <- ggplot(resid_data, aes(x = residuals)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 30, fill = "lightblue", color = "black", alpha = 0.6) +
  geom_density(color = "darkblue", linewidth = 1) +
  facet_wrap(~model, ncol = 1, scales = "free") +
  geom_text(data = resid_data %>%
              group_by(model, shapiro_p) %>%
              summarise(.groups = "drop", x = Inf, y = Inf),
            aes(x = x, y = y, label = shapiro_p),
            hjust = 1.1, vjust = 1.5, size = 3.5) +
  theme_minimal() +
  labs(
    x = "Residuals",
    y = "Density",
    title = "France"
  )

#### Bonaire ####
data <- read.csv("data/Bonaire_data.csv",sep=";",header=T,row.names=1)

# Transform data
hist(data$ES_totbiom)
logBiom <- log10(data$ES_totbiom+1)
data <- cbind(data,logBiom)
data <- data[data$BIO_qual>0.70,]
hist(data$logBiom,xlab ="log(Total Fish Biomass)", main=NULL)

# Variable selection
mod.Biom <- lm(logBiom ~ ENV_Depth+No.Take+ENV_Wave_energy+ENV_Terrace_width+ENV_Substratum_hard+ENV_Rugosity+
                 ENV_Height_vertical+ENV_Coral_live+ ENV_Coral_reef_building+
                 BIO_Taxo_Ric+BIO_Taxo_Ent+BIO_Func_Ric+BIO_Func_Ent+BIO_Phyl_Ric+BIO_Phyl_Ent, data)
stepAIC(mod.Biom)

tax_ric <- lm(BIO_Taxo_Ric ~ENV_Depth+No.Take+ENV_Wave_energy+ENV_Terrace_width+ENV_Substratum_hard+ENV_Rugosity+
                ENV_Height_vertical+ENV_Coral_live+ ENV_Coral_reef_building, data)
stepAIC(tax_ric)

tax_ent <- lm(BIO_Taxo_Ent ~ENV_Depth+No.Take+ENV_Wave_energy+ENV_Terrace_width+ENV_Substratum_hard+ENV_Rugosity+
                ENV_Height_vertical+ENV_Coral_live+ ENV_Coral_reef_building, data)
stepAIC(tax_ent)

fun_ent <- lm(BIO_Func_Ent ~ ENV_Depth+No.Take+ENV_Wave_energy+ENV_Terrace_width+ENV_Substratum_hard+ENV_Rugosity+
                ENV_Height_vertical+ENV_Coral_live+ ENV_Coral_reef_building, data)
stepAIC(fun_ent)


# Parsimonious model
model.simp <- psem(
  lm(BIO_Taxo_Ric ~ ENV_Depth + No.Take + ENV_Wave_energy + 
       ENV_Rugosity + ENV_Coral_live + ENV_Coral_reef_building, 
     data = data),
  lm(BIO_Taxo_Ent ~ ENV_Depth + ENV_Rugosity, data = data),
  lm(BIO_Func_Ent ~ ENV_Depth + No.Take + ENV_Height_vertical + 
       ENV_Coral_live, data = data),
  lm(logBiom ~ ENV_Depth + ENV_Terrace_width + ENV_Substratum_hard + 
       ENV_Coral_live + ENV_Coral_reef_building + BIO_Taxo_Ric + 
       BIO_Taxo_Ent + BIO_Func_Ent, data = data)
)

anova(model.simp)
coefs(model.simp)
rsquared(model.simp)

plot(model.simp)
plot.psem(model.simp)
x <- data.frame(shape = "rectangle", color = "black",fillcolor = "orange",y=1:3)
y <- data.frame(shape = "rectangle", color = "black",fillcolor = "black",y=4)
z <- data.frame(shape = "rectangle", color = "black",fillcolor = "blue",y=5:13)
col <- rbind(x,y,z)
plot.psem(model.simp,node_attrs = col)
line <- data.frame(style = "solid", color = "black")
plot.psem(model.simp,node_attrs = col,edge_attrs=line) 

# Model residuals
# Objects of class psem do not have global residuals, they have several imbricated models
model.simp <- psem(
  lm(BIO_Taxo_Ric ~ ENV_Depth + No.Take + ENV_Wave_energy + 
       ENV_Rugosity + ENV_Coral_live + ENV_Coral_reef_building, 
     data = data),
  lm(BIO_Taxo_Ent ~ ENV_Depth + ENV_Rugosity, data = data),
  lm(BIO_Func_Ent ~ ENV_Depth + No.Take + ENV_Height_vertical + 
       ENV_Coral_live, data = data),
  lm(logBiom ~ ENV_Depth + ENV_Terrace_width + ENV_Substratum_hard + 
       ENV_Coral_live + ENV_Coral_reef_building + BIO_Taxo_Ric + 
       BIO_Taxo_Ent + BIO_Func_Ent, data = data)
)


mod_list <- list(
  TaxR_model =  lm(BIO_Taxo_Ric ~ ENV_Depth + No.Take + ENV_Wave_energy + 
                     ENV_Rugosity + ENV_Coral_live + ENV_Coral_reef_building, 
                   data = data),
  TaxE_model = lm(BIO_Taxo_Ent ~ ENV_Depth + ENV_Rugosity, data = data),
  FunE_model = lm(BIO_Func_Ent ~ ENV_Depth + No.Take + ENV_Height_vertical + 
                    ENV_Coral_live, data = data),
  logBiom_model = lm(logBiom ~ ENV_Depth + ENV_Terrace_width + ENV_Substratum_hard + 
                       ENV_Coral_live + ENV_Coral_reef_building + BIO_Taxo_Ric + 
                       BIO_Taxo_Ent + BIO_Func_Ent, data = data)
)
res_list <- lapply(mod_list, residuals)

resid_data <- bind_rows(
  lapply(names(mod_list), function(name) {
    res <- residuals(mod_list[[name]])
    pval <- shapiro.test(res)$p.value
    data.frame(
      residuals = res,
      model = name,
      shapiro_p = sprintf("p = %.3f", pval)
    )
  })
)

p_dens_bonaire <- ggplot(resid_data, aes(x = residuals)) +
  geom_density(fill = "steelblue", alpha = 0.6) +
  facet_wrap(~model, ncol = 1, scales = "free") +
  theme_minimal() +
  labs(
    x = "Residuals",
    y = "Density",
    title = "Bonaire"
  )

p_combo_bonaire <- ggplot(resid_data, aes(x = residuals)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 30, fill = "lightblue", color = "black", alpha = 0.6) +
  geom_density(color = "darkblue", linewidth = 1) +
  facet_wrap(~model, ncol = 1, scales = "free") +
  geom_text(data = resid_data %>%
              group_by(model, shapiro_p) %>%
              summarise(.groups = "drop", x = Inf, y = Inf),
            aes(x = x, y = y, label = shapiro_p),
            hjust = 1.1, vjust = 1.5, size = 3.5) +
  theme_minimal() +
  labs(
    x = "Residuals",
    y = "Density",
    title = "Bonaire"
  )


#### Martinique ####
data_mart_var <- read.csv("data/Martinique_data.csv",sep=";",header=T)
data_mart_env <- read.csv("data/Martinique_env.csv",sep=";",header=T)
data_mart_merge <- merge(data_mart_var, data_mart_env, by = "spl_id")

# Transform data
data_mart <- data_mart_merge %>% 
  mutate_at(c("coral_cover", "macroalgae_cover", "turf_cover", "other_cover"), as.numeric) %>% 
  mutate(spl_management = recode(spl_management, "Fishing ?" = "Fishing")) %>%
  mutate(spl_management = recode(spl_management, "chlordecone_no-take_zone" = "No-take area")) %>% 
  mutate(spl_management = as.factor(spl_management)) %>% 
  mutate(spl_management = as.numeric(spl_management))

hist(data_mart$ES_totbiom)
logBiom=log10(data_mart$ES_totbiom)
hist(logBiom)
data_mart=cbind(data_mart,logBiom)

# Variable selection
model <- lm(logBiom ~ BIO_Taxo_Ric + BIO_Taxo_Ent + BIO_Func_Ric + BIO_Func_Ent + BIO_Phyl_Ric + BIO_Phyl_Ent +
              spl_management + spl_depth + coral_cover + macroalgae_cover + turf_cover + other_cover, 
            data_mart)

stepAIC(model)

tax_ric <- lm(BIO_Taxo_Ric ~ spl_management + spl_depth + coral_cover + macroalgae_cover + turf_cover + other_cover, 
              data_mart)
stepAIC(tax_ric)

tax_ent <- lm(BIO_Taxo_Ent ~ spl_management + spl_depth + coral_cover + macroalgae_cover + turf_cover + other_cover, 
              data_mart)
stepAIC(tax_ent)


model.final <- lm(logBiom ~ BIO_Taxo_Ent + BIO_Phyl_Ent, data = data_mart)

summary(model.final)
visreg(model.final,scale="response")

# Parsimonious model
model.simp <- psem(
  lm(BIO_Taxo_Ric ~ spl_depth + coral_cover + macroalgae_cover + turf_cover + other_cover, 
     data = data_mart),
  lm(BIO_Taxo_Ent ~ spl_depth  + macroalgae_cover + turf_cover + other_cover, 
     data = data_mart),
  lm(logBiom ~ BIO_Taxo_Ric + BIO_Taxo_Ent + spl_management + spl_depth + 
       macroalgae_cover, data = data_mart)
)

anova(model.simp)
coefs(model.simp)
rsquared(model.simp)


plot(model.simp)
plot.psem(model.simp)
x <- data.frame(shape = "rectangle", color = "black",fillcolor = "orange",y=1:2)
y <- data.frame(shape = "rectangle", color = "black",fillcolor = "black",y=3)
z <- data.frame(shape = "rectangle", color = "black",fillcolor = "blue",y=4:7)
col <- rbind(x,y,z)
line <- data.frame(style = "solid", color = "red")
plot.psem(model.simp) 

# pdf("results/sem_Martinique.pdf",height = 10, width = 10, pointsize = 5, useDingbats = FALSE)
plot(model.simp)
# dev.off()


# Model residuals
# Objects of class psem do not have global residuals, they have several imbricated models
mod_list <- list(
  TaxR_model =  lm(BIO_Taxo_Ric ~ spl_depth + coral_cover + macroalgae_cover + turf_cover + other_cover, 
                   data = data_mart),
  TaxE_model = lm(BIO_Taxo_Ent ~ spl_depth  + macroalgae_cover + turf_cover + other_cover, 
                  data = data_mart),
  logBiom_model = lm(logBiom ~ BIO_Taxo_Ric + BIO_Taxo_Ent + spl_management + spl_depth + 
                       macroalgae_cover, data = data_mart)
)
res_list <- lapply(mod_list, residuals)

resid_data <- bind_rows(
  lapply(names(mod_list), function(name) {
    res <- residuals(mod_list[[name]])
    pval <- shapiro.test(res)$p.value
    data.frame(
      residuals = res,
      model = name,
      shapiro_p = sprintf("p = %.3f", pval)
    )
  })
)

p_dens_martinique <- ggplot(resid_data, aes(x = residuals)) +
  geom_density(fill = "steelblue", alpha = 0.6) +
  facet_wrap(~model, ncol = 1, scales = "free") +
  theme_minimal() +
  labs(
    x = "Residuals",
    y = "Density",
    title = "Martinique"
  )

p_combo_martinique <- ggplot(resid_data, aes(x = residuals)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 30, fill = "lightblue", color = "black", alpha = 0.6) +
  geom_density(color = "darkblue", linewidth = 1) +
  facet_wrap(~model, ncol = 1, scales = "free") +
  geom_text(data = resid_data %>%
              group_by(model, shapiro_p) %>%
              summarise(.groups = "drop", x = Inf, y = Inf),
            aes(x = x, y = y, label = shapiro_p),
            hjust = 1.1, vjust = 1.5, size = 3.5) +
  theme_minimal() +
  labs(
    x = "Residuals",
    y = "Density",
    title = "Martinique"
  )

#### Ireland ####
data <- read.csv("data/Ireland_data.csv",header=T,row.names=1)

# Transform data
data <- na.omit(data)
hist(data$FBiom,xlab = "Total Fished Biomass", main=NULL)
logFB <- log10(data$FBiom)
data <- cbind(data,logFB)
hist(data$logFB,xlab = "log(Total Fished Biomass)", main=NULL)

# Variable selection
mod.FB <- lm(logFB ~ Depth+SST+Chla+Eup+Sal+NPP+Sub+
               TaxR+TaxE+FunR+FunE+PhyR+PhyE, data)
stepAIC(mod.FB)

tax_ric <- lm(TaxR ~ Depth+SST+Chla+Eup+Sal+NPP+Sub, data)
stepAIC(tax_ric)

tax_ent <- lm(TaxE ~ Depth+SST+Chla+Eup+Sal+NPP+Sub, data)
stepAIC(tax_ent)

fun_ric <- lm(FunR ~ Depth+SST+Chla+Eup+Sal+NPP+Sub, data)
stepAIC(fun_ric)

phy_ric <- lm(PhyR ~ Depth+SST+Chla+Eup+Sal+NPP+Sub, data)
stepAIC(phy_ric)

phy_ent <- lm(PhyE ~ Depth+SST+Chla+Eup+Sal+NPP+Sub, data)
stepAIC(phy_ent)


# Parsimonious model
model.simp <- psem(
  lm(TaxR ~ Depth + SST + Chla + Eup, data),
  lm(TaxE ~ SST + Chla + NPP, data),
  lm(FunR ~ Depth + Eup + Sal + NPP, data),
  lm(PhyR ~ Depth + SST + Eup + Sal, data),
  lm(PhyE ~ SST + NPP, data),
  lm(logFB ~ SST + Chla + Eup + Sal + NPP + TaxR + 
       TaxE + FunR + PhyR + PhyE, data)
)

anova(model.simp)

coefs(model.simp)

rsquared(model.simp)

plot(model.simp)

plot.psem(model.simp)

x <- data.frame(shape = "rectangle", color = "black",fillcolor = "orange",y=1:5)
y <- data.frame(shape = "rectangle", color = "black",fillcolor = "black",y=6)
z <- data.frame(shape = "rectangle", color = "black",fillcolor = "blue",y=7:12)

col <- rbind(x,y,z)

plot.psem(model.simp,node_attrs = col)

line <- data.frame(style = "solid", color = "black")

plot.psem(model.simp,node_attrs = col,edge_attrs=line) 

# Model residuals
# Objects of class psem do not have global residuals, they have several imbricated models
mod_list <- list(
  TaxR_model = lm(TaxR ~ Depth + SST + Chla + Eup, data),
  TaxE_model = lm(TaxE ~ SST + Chla + NPP, data),
  FunR_model = lm(FunR ~ Depth + Eup + Sal + NPP, data),
  PhyR_model = lm(PhyR ~ Depth + SST + Eup + Sal, data),
  PhyE_model = lm(PhyE ~ SST + NPP, data),
  logFB_model = lm(logFB ~ SST + Chla + Eup + Sal + NPP + TaxR + 
                     TaxE + FunR + PhyR + PhyE, data)
)
res_list <- lapply(mod_list, residuals)

resid_data <- bind_rows(
  lapply(names(mod_list), function(name) {
    res <- residuals(mod_list[[name]])
    pval <- shapiro.test(res)$p.value
    data.frame(
      residuals = res,
      model = name,
      shapiro_p = sprintf("p = %.3f", pval)
    )
  })
)

p_dens_ireland <- ggplot(resid_data, aes(x = residuals)) +
  geom_density(fill = "steelblue", alpha = 0.6) +
  facet_wrap(~model, ncol = 1, scales = "free") +
  theme_minimal() +
  labs(
    x = "Residuals",
    y = "Density",
    title = "Ireland"
  )

p_combo_ireland <- ggplot(resid_data, aes(x = residuals)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 30, fill = "lightblue", color = "black", alpha = 0.6) +
  geom_density(color = "darkblue", linewidth = 1) +
  facet_wrap(~model, ncol = 1, scales = "free") +
  geom_text(data = resid_data %>%
              group_by(model, shapiro_p) %>%
              summarise(.groups = "drop", x = Inf, y = Inf),
            aes(x = x, y = y, label = shapiro_p),
            hjust = 1.1, vjust = 1.5, size = 3.5) +
  theme_minimal() +
  labs(
    x = "Residuals",
    y = "Density",
    title = "Ireland"
  )


#### Med Spain ####
data_spain <- read.csv("data/Spain_data.csv",sep=";",header=T,row.names=1)

# Transform data
hist(data_spain$densf)
logESF <- log10(data_spain$densf+1)
data_spain <- cbind(data_spain,logESF)
data_spain <- data_spain[data_spain$BIOq>0.70,]
hist(data_spain$logESF,xlab ="log(Total number of fish individuals)", main=NULL)

# Variable selection
mod.ESF <- lm(logESF ~ depth+maerl+reserve+fishing_effort+
             TaxR+TaxE+FunR+FunE+PhyR+PhyE, data_spain)
stepAIC(mod.ESF)

tax_ric <- lm(TaxR ~ depth+maerl+reserve+fishing_effort, data_spain)
stepAIC(tax_ric)

fun_ric <- lm(FunR ~ depth+maerl+reserve+fishing_effort, data_spain)
stepAIC(fun_ric)


# Parsimonious model
model.simp <- psem(
  lm(TaxR ~ maerl + fishing_effort, data_spain),
  lm(FunR ~ depth + reserve + fishing_effort, data_spain),
  lm(logESF ~ TaxR + FunR, data_spain)
)

anova(model.simp)
coefs(model.simp)
rsquared(model.simp)

plot.psem(model.simp)
x <- data.frame(shape = "rectangle", color = "black",fillcolor = "orange",y=1:2)
y <- data.frame(shape = "rectangle", color = "black",fillcolor = "black",y=3)
z <- data.frame(shape = "rectangle", color = "black",fillcolor = "blue",y=4:7)
col <- rbind(x,y,z)

line <- data.frame(style = "solid", color = "black")

plot.psem(model.simp,node_attrs = col,edge_attrs=line) 


# Model residuals
# Objects of class psem do not have global residuals, they have several imbricated models
mod_list <- list(
  TaxR_model = lm(TaxR ~ maerl + fishing_effort, data_spain),
  FunR_model = lm(FunR ~ depth + reserve + fishing_effort, data_spain),
  logESF_model= lm(logESF ~ TaxR + FunR, data_spain)
)
res_list <- lapply(mod_list, residuals)


resid_data <- bind_rows(
  lapply(names(mod_list), function(name) {
    res <- residuals(mod_list[[name]])
    pval <- shapiro.test(res)$p.value
    data.frame(
      residuals = res,
      model = name,
      shapiro_p = sprintf("p = %.3f", pval)
    )
  })
)

p_dens_spain <- ggplot(resid_data, aes(x = residuals)) +
  geom_density(fill = "steelblue", alpha = 0.6) +
  facet_wrap(~model, ncol = 1, scales = "free") +
  theme_minimal() +
  labs(
    x = "Residuals",
    y = "Density",
    title = "Spain"
  )

p_combo_spain <- ggplot(resid_data, aes(x = residuals)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 30, fill = "lightblue", color = "black", alpha = 0.6) +
  geom_density(color = "darkblue", linewidth = 1) +
  facet_wrap(~model, ncol = 1, scales = "free") +
  geom_text(data = resid_data %>%
              group_by(model, shapiro_p) %>%
              summarise(.groups = "drop", x = Inf, y = Inf),
            aes(x = x, y = y, label = shapiro_p),
            hjust = 1.1, vjust = 1.5, size = 3.5) +
  theme_minimal() +
  labs(
    x = "Residuals",
    y = "Density",
    title = "Spain"
  )


#### Residuals plots ####
ls_plots <- list(
  p_combo_france,
  p_combo_bonaire,
  p_combo_martinique,
  p_combo_ireland,
  p_combo_spain
)

Fig_res <- (
  (p_combo_france / p_combo_bonaire / p_combo_martinique ) | 
    (p_combo_spain / p_combo_ireland + plot_layout(heights = c(1, 2)))
) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A')

# save the figure
ggsave("results/Models_residuals2.jpeg",
       plot = Fig_res,
       width = 15,
       height = 25,
       units = "cm")
pdf("results/Models_residuals2.pdf",height = 15, width = 10, pointsize = 5, useDingbats = FALSE)
Fig_res
dev.off()


