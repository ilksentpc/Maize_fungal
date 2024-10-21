##Maize_fungal###

###Alpha_diversity####

library(phyloseq)

otu_leafbac <- read.csv("otu_fungal2024.csv", row.names = 1)
tax_leaf <- read.csv("tax_fungal_2024.csv", row.names = 1)
otumat <- as.matrix(otu_leafbac)
taxmat <- as.matrix(tax_leaf)
metadata <- read.csv("metadata_ITS.csv", row.names = 1)
sampledata = sample_data(data.frame(metadata))
OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
OTU

physeq = phyloseq(OTU, TAX, sampledata)
physeq

library(microeco)
library(file2meco)
library(magrittr)
meco_dataset <- phyloseq2meco(physeq)
meco_dataset

meco_dataset$sample_table$Plant_Group %<>% factor(levels = c("Teosinte", "Mexican maize", "US maize"))
meco_dataset$sample_table$Genotype %<>% factor(levels = c("Perennial teosinte", "Balsas teosinte", "Mexican landrace", "Mexican inbred", "US landrace", "US inbred"))

meco_dataset$sample_table$Zone %<>% factor(levels = c("BS", "RS", "PH"))

meco_dataset$sample_table$Genotype

###shannon
t1 <- trans_alpha$new(dataset = meco_dataset, group = "Genotype", by_group = "Zone")
t1$plot_alpha(measure = "Shannon")

t2 <- trans_alpha$new(dataset = meco_dataset, group = "Plant_Group", by_group = "Zone")
t2$data_alpha$Genotype

Plant_group_alpha_plot_shannon <- t2$plot_alpha(measure = "Shannon")

png("Plant_group_alpha_plot_shannon.png", width=13, height=9, units="in", res=300)
Plant_group_alpha_plot_shannon + 
  theme(
    text = element_text(family = "Arial"),  # Sets all text to Arial
    plot.title = element_text(size = 16, face = "bold", colour = "black"),  # Title styling
    axis.title.x = element_text(size = 16, face = "bold", colour = "black"),  # X-axis label styling
    axis.title.y = element_text(size = 16, face = "bold", colour = "black"),  # Y-axis label styling
    axis.text.x = element_text(size = 14, face = "bold", colour = "black"),  # X-axis tick labels
    axis.text.y = element_text(size = 14, face = "bold", colour = "black"),  # Y-axis tick labels
    legend.text = element_text(size = 12, face = "bold", colour = "black"),   # Make legend text bold and black
    legend.title = element_text(size = 14, face = "bold", colour = "black")  # Make legend title bold and black
  )
dev.off()

#Genotype nested within Plant group
meco_dataset$sample_table$Genotype_PlantGroup <- with(meco_dataset$sample_table, paste(Genotype, Plant_Group, sep = " : "))

# For Perennial teosinte
for (i in 1:nrow(meco_dataset$sample_table)) {
  if (meco_dataset$sample_table$Genotype[i] == "Perennial teosinte" && meco_dataset$sample_table$Plant_Group[i] != "Teosinte") {
    meco_dataset$sample_table$Genotype_PlantGroup[i] <- NA
  }
}

# For Balsas teosinte
for (i in 1:nrow(meco_dataset$sample_table)) {
  if (meco_dataset$sample_table$Genotype[i] == "Balsas teosinte" && meco_dataset$sample_table$Plant_Group[i] != "Teosinte") {
    meco_dataset$sample_table$Genotype_PlantGroup[i] <- NA
  }
}

# For US inbred
for (i in 1:nrow(meco_dataset$sample_table)) {
  if (meco_dataset$sample_table$Genotype[i] == "US inbred" && meco_dataset$sample_table$Plant_Group[i] != "US maize") {
    meco_dataset$sample_table$Genotype_PlantGroup[i] <- NA
  }
}

# For US landrace
for (i in 1:nrow(meco_dataset$sample_table)) {
  if (meco_dataset$sample_table$Genotype[i] == "US landrace" && meco_dataset$sample_table$Plant_Group[i] != "US maize") {
    meco_dataset$sample_table$Genotype_PlantGroup[i] <- NA
  }
}

# For Mexican landrace
for (i in 1:nrow(meco_dataset$sample_table)) {
  if (meco_dataset$sample_table$Genotype[i] == "Mexican landrace" && meco_dataset$sample_table$Plant_Group[i] != "Mexican maize") {
    meco_dataset$sample_table$Genotype_PlantGroup[i] <- NA
  }
}

# For Mexican inbred
for (i in 1:nrow(meco_dataset$sample_table)) {
  if (meco_dataset$sample_table$Genotype[i] == "Mexican inbred" && meco_dataset$sample_table$Plant_Group[i] != "Mexican maize") {
    meco_dataset$sample_table$Genotype_PlantGroup[i] <- NA
  }
}

#remove rows with NA 
valid_rows <- !is.na(meco_dataset$sample_table$Genotype_PlantGroup)
meco_dataset$sample_table <- meco_dataset$sample_table[valid_rows, ]

meco_dataset$sample_table$Genotype_PlantGroup <- factor(meco_dataset$sample_table$Genotype_PlantGroup, 
                                                        levels = c("Perennial teosinte : Teosinte", "Balsas teosinte : Teosinte", 
                                                                   "Mexican landrace : Mexican maize", "Mexican inbred : Mexican maize", 
                                                                   "US landrace : US maize", "US inbred : US maize"))



t3 <- trans_alpha$new(dataset = meco_dataset, group = "Genotype_PlantGroup", by_group = "Zone")
t3$plot_alpha(measure = "Shannon")

library(ggplot2)
alpha_plot_shannon <- t3$plot_alpha(measure = "Shannon")
png("shannon_nestedPlantGroup_styled.png", width=13, height=9, units="in", res=300)
alpha_plot_shannon + 
  theme(
    text = element_text(family = "Arial"),  # Sets all text to Arial
    plot.title = element_text(size = 16, face = "bold", colour = "black"),  # Title styling
    axis.title.x = element_text(size = 16, face = "bold", colour = "black"),  # X-axis label styling
    axis.title.y = element_text(size = 16, face = "bold", colour = "black"),  # Y-axis label styling
    axis.text.x = element_text(size = 14, face = "bold", colour = "black"),  # X-axis tick labels
    axis.text.y = element_text(size = 14, face = "bold", colour = "black"),  # Y-axis tick labels
    legend.text = element_text(size = 12, face = "bold", colour = "black"),   # Make legend text bold and black
    legend.title = element_text(size = 14, face = "bold", colour = "black")  # Make legend title bold and black
  )
dev.off()

#####
meco_dataset$alpha_diversity
meco_dataset$sample_table$Shannon <- meco_dataset$alpha_diversity$Shannon

#Nested ANOVA on Shannon diversity
anova_result <- aov(Shannon ~ Plant_Group / Genotype, data = meco_dataset$sample_table)
summary(anova_result)

# Tukey's Honest
tukey_result <- TukeyHSD(anova_result, which = "Plant_Group:Genotype")
print(tukey_result)

print(rownames(tukey_result$`Plant_Group:Genotype`))

comparison_labels <- c(
  "Teosinte:Balsas teosinte-Teosinte:Perennial teosinte", 
  "Mexican maize:Mexican landrace-Teosinte:Balsas teosinte",
  "Mexican maize:Mexican inbred-Mexican maize:Mexican landrace",
  "US maize:US landrace-Mexican maize:Mexican landrace",
  "US maize:US inbred-US maize:US landrace"
)


filtered_comparisons <- tukey_result$`Plant_Group:Genotype`[comparison_labels, ]
print(filtered_comparisons)


###FOR EACH ZONE

meco_dataset$sample_table$Shannon <- meco_dataset$alpha_diversity$Shannon

unique_zones <- unique(meco_dataset$sample_table$Zone)  # Assuming 'Zone' is the variable name


#comparisons labels for Tukey's HSD
comparison_labels <- c(
  "Teosinte:Balsas teosinte-Teosinte:Perennial teosinte", 
  "Mexican maize:Mexican landrace-Teosinte:Balsas teosinte",
  "Mexican maize:Mexican inbred-Mexican maize:Mexican landrace",
  "US maize:US landrace-Mexican maize:Mexican landrace",
  "US maize:US inbred-US maize:US landrace"
)

# Loop over each zone (BS, RS, PH) 
for (zone in unique_zones) {
  # Filter the dataset for the current zone
  zone_data <- subset(meco_dataset$sample_table, Zone == zone)
  
  #Nested ANOVA on Shannon diversity
  anova_result <- aov(Shannon ~ Plant_Group / Genotype, data = zone_data)
  
  print(paste("ANOVA results for zone:", zone))
  print(summary(anova_result))
  
  # Perform Tukey's 
  tukey_result <- TukeyHSD(anova_result, which = "Plant_Group:Genotype")
  
  # the row names (comparison labels) in the TukeyHSD result 
  print(paste("Tukey's HSD comparison labels for zone:", zone))
  print(rownames(tukey_result$`Plant_Group:Genotype`))
  
  filtered_comparisons <- tukey_result$`Plant_Group:Genotype`[comparison_labels, ]
  
  print(paste("Filtered comparisons for zone:", zone))
  print(filtered_comparisons)
}
