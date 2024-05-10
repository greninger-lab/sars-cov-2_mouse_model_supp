library(tidyverse)
library(cowplot)
library(ggrepel)

# Not in function
'%!in%' <- function(x,y)!('%in%'(x,y))

####################################################################################
# Set the directory based on command line argument or hard-coded location
####################################################################################
# rava.parent.dir should be a directory that contains all the RAVA output directories
# that you want to analyze
rava.parent.dir <- "/Volumes/lizso_backup_drive/GSU_CL_VOCBeta_Swift/sars-cov-2_mouse_model_supp"

args <- commandArgs(trailingOnly=TRUE)

if(!is.na(args[1])){
  rava.parent.dir <- args[1]
}

####################################################################################
# make a directory to put the figures in the rava parent directory
####################################################################################
figure.location <- paste0(rava.parent.dir, "/figures")

if(!dir.exists(figure.location)){
  dir.create(figure.location)
}
####################################################################################
# read in visualization information
####################################################################################

# set the working directory to wherever you downloaded the data :)
setwd(rava.parent.dir)

# Get all the visualization.csv files
# Note, we will read the inoculum RAVA visualization file as well, but this data is filtered 
# out later on in the script. 
visualization.files.simple <- list.files(path = ".", pattern = "_visualization.csv", recursive = TRUE, include.dirs = TRUE)

visualization.table <- visualization.files.simple %>% 
  map_df(~read_csv(.))

# make columns for treatment, timepoint, and animal from the "Passage" column
visualization.table$treatment <- str_split_i(visualization.table$Passage, "_",1)
visualization.table$timepoint <- str_split_i(visualization.table$Passage, "_",2)
visualization.table$animal <- str_split_i(visualization.table$Passage, "_",3)
visualization.table$experiment <- str_split_i(visualization.table$Passage, "_",4)

####################################################################################
# Combine the data frames and put figures in the new folder
####################################################################################
setwd(figure.location)

####################################################################################
# Filter the variant list, keep variants with:
# mean allele frequency >= 5%, 
# mean depth >= 100, 
# and appears in both technical replicates (count = 2)
####################################################################################
# Get the relevant columns, use the "mature peptide" to uniquely identify variants ... they are listed twice for orf1a and orf1ab based on rava annotations :-/
# Get the mean allele frequency and depth across technical animals (prep1 and prep2), count the number of technical animals where the variant appears (1 or 2)
visualization.table <- visualization.table %>% select(treatment,timepoint,animal,Syn,MatPeptide,`Amino Acid Change`, Position, LetterChange,AF,Depth,experiment) %>%
  mutate(MatPeptide = case_when(MatPeptide == "-" ~ `Amino Acid Change`, MatPeptide != "-" ~ MatPeptide)) %>%
  select(-c(`Amino Acid Change`)) %>% unique() %>% group_by(treatment,timepoint,animal,Syn,Position,LetterChange, MatPeptide) %>% summarise(AF=mean(AF),Depth=mean(Depth), count = n())

# Make a Sample id column by [treatment]_[timepoint]_[animal]
visualization.table$Sample <- paste(visualization.table$treatment, visualization.table$timepoint, visualization.table$animal, sep="_")

# Filter the variant list, keep variants with mean allele frequency >= 5%, mean depth >= 100, and appears in both technical replicates (count = 2)
visualization.table.filtered <- visualization.table %>% filter(AF >= 5, Depth >= 100, count > 1)

####################################################################################
# Summarize the number of variants per animal so we can display this on top of 
# the scatter plot of all variants. 
# Since there are different numbers of animals in each group, it's not valid to 
# just display the scatter and intuit the number of variants per group :-(
####################################################################################
num.variants.per.animal <- visualization.table.filtered %>% ungroup() %>% group_by(treatment,timepoint,animal) %>% summarise(num.variants = n())

####################################################################################
# One-way ANOVA of mean variants per animal
####################################################################################
# One-way ANOVA to compare mean mutations per animal across treatment groups at D21
test.D21 <- num.variants.per.animal %>% ungroup() %>% filter(timepoint == "D21")
res.aov.D21 <- aov(num.variants ~ treatment, data = test.D21)
tukey_hsd.D21 <- TukeyHSD(res.aov.D21)

tukey_hsd.D21
write.csv(x = tukey_hsd.D21$treatment, file = "stats_mean-var-per-animal_D21.csv")

# One-way ANOVA to compare mean mutations per animal across timepoints for vehicle-treated animals
test.veh <- num.variants.per.animal %>% ungroup() %>% filter(treatment == "Vehicle")
res.aov.veh <- aov(num.variants ~ timepoint, data = test.veh)
tukey_hsd.veh <- TukeyHSD(res.aov.veh)

tukey_hsd.veh
write.csv(x = tukey_hsd.veh$timepoint, file = "stats_mean-var-per-animal_Vehicle.csv")

# One-way ANOVA to compare mean mutations per animal across timepoints for Nirmatrelvir/Ritonavir-treated animals
test.pax <- num.variants.per.animal %>% ungroup() %>% filter(treatment == "Nirmatrelvir")
res.aov.pax <- aov(num.variants ~ timepoint, data = test.pax)
tukey_hsd.pax <- TukeyHSD(res.aov.pax)

tukey_hsd.pax
write.csv(x = tukey_hsd.pax$timepoint, file = "stats_mean-var-per-animal_Nirmatrelvir.csv")

####################################################################################
# Prepare the variants per animal table for plotting
####################################################################################
# Factor the table of mutations per animal by treatment group
num.variants.per.animal$treatment <- factor(num.variants.per.animal$treatment, levels = c("Vehicle", "Nirmatrelvir","Molnupiravir","4'-FlU"))

# Get mean and sd variants per group and number of animals in each group 
num.variants.per.animal.mean <- num.variants.per.animal %>% ungroup() %>% group_by(treatment, timepoint) %>% summarise(mean.variants = mean(num.variants),
                                                                                                                         sd.variants = sd(num.variants),
                                                                                                                         num.animals = n())

# Make a column of labels for the venn diagrams including mean, sd, and number of animals (n) for mean mutations per animal
num.variants.per.animal.mean <- num.variants.per.animal.mean %>% mutate(labels = paste0("mean=", round(mean.variants, digits = 1), "\n sd=", round(sd.variants, digits = 1), "\n n=", num.animals))

# Make a column for the circle size ... ultimately did not scale, so mean.variants = circle.size
num.variants.per.animal.mean$circle.size <- num.variants.per.animal.mean$mean.variants

# Separate the dataframe of mean variants per animal by timepoint
num.variants.per.animal.mean.D14 <- num.variants.per.animal.mean %>% filter(timepoint == "D14")
num.variants.per.animal.mean.D21 <- num.variants.per.animal.mean %>% filter(timepoint == "D21")
num.variants.per.animal.mean.D28 <- num.variants.per.animal.mean %>% filter(timepoint == "D28")

####################################################################################
# Prepare the table of relevant variants for plotting
####################################################################################

# Clean up the names of the mature peptides so they're easier to display
# Get the list of all variants to be plotted
MatPeptides.unique <- unique(visualization.table.filtered$MatPeptide)
# Make a copy to update the names
MatPeptides.unique.new <- MatPeptides.unique
# Remove all the colons and semicolons from the name
MatPeptides.unique.new <- str_replace_all(MatPeptides.unique.new,c(";" = "", ":" = ""))
# The protein name is the first element in the list separated by spaces
MatPeptides.unique.new.protein <- str_split_i(MatPeptides.unique.new, " ", 1)

# Shorten the protein names:
  # RNA-dependent_RNA_polymerase_rib_26 -> nsp12
  # 2'-O-ribose_methyltransferase -> nsp16
  # remove "_protein" and "_glycoprotein"
  # nucleocapsid_phosphoprotein -> N
  # membrane -> M
  # envelope -> E
  # leader -> nsp1
  # helicase -> nsp13
  # endoRNAse -> nsp15
  # 3'-to-5'_exonuclease -> nsp14
  # 3C-likease -> nsp5
  # surface -> S
MatPeptides.unique.new.protein <- str_replace_all(
  str_replace_all(
    str_replace_all(
      str_replace_all(
        MatPeptides.unique.new.protein, "RNA-dependent_RNA_polymerase_rib_26", "nsp12"), "2'-O-ribose_methyltransferase", "nsp16"), "_protein", ""), "_glycoprotein", "")
MatPeptides.unique.new.protein <- str_replace_all(
  str_replace_all(
    str_replace_all(
      str_replace_all(
        str_replace_all(
          str_replace_all(
            str_replace_all(
              MatPeptides.unique.new.protein, "nucleocapsid_phosphoprotein", "N"), "membrane", "M"), "envelope", "E"), "leader", "nsp1"), "helicase", "nsp13"), "endoRNAse", "nsp15"), "3'-to-5'_exonuclease", "nsp14")
MatPeptides.unique.new.protein <- str_replace_all(
  str_replace_all(
    MatPeptides.unique.new.protein, "3C-likease", "nsp5"), "surface", "S")

# The amino acid change is the second element in the list separated by spaces
MatPeptides.unique.new.aachange <- str_split_i(MatPeptides.unique.new, " ", 2)

# combine the shortened protein names with the amino acid change, separated by a colon
MatPeptides.unique.new <- paste(MatPeptides.unique.new.protein, MatPeptides.unique.new.aachange, sep=":")

# make a dictionary with keys = old annotations and values = new annotations
names(MatPeptides.unique.new) <- MatPeptides.unique

# Use the dictionary to change the names in the dataframe of all relevant variants
visualization.table.filtered$MatPeptide <- str_replace_all(visualization.table.filtered$MatPeptide, MatPeptides.unique.new)

# Save the table of all relevant variants
write.csv(visualization.table.filtered, "figure3H_NGS_table.csv")

# Make a list of variants per treatment group, calculate mean of mean allele frequency and depth for variants that appear in >1 animal in a group. 
# group = treatment + timepoint
variants.per.group <- visualization.table.filtered %>% ungroup() %>% group_by(treatment,timepoint, Syn, MatPeptide) %>% summarise(AF.mean = mean(AF),
                                                                                                                                   Depth.mean = mean(Depth),
                                                                                                                                   num.animals = n())

# Factor the table by treatment
variants.per.group$treatment <- factor(variants.per.group$treatment, levels = c("Vehicle", "Nirmatrelvir","Molnupiravir","4'-FlU"))

# Change the integer number of animals in a group with a given variant to a character so we can 
# factor that too
variants.per.group$num.animals <- as.character(variants.per.group$num.animals)
variants.per.group$num.animals <- factor(variants.per.group$num.animals, levels = c("1","2","3","4","5","6"))

# Keep the annotation for nonsynonymous variants with allele frequency >=20%, or any variants that appears in >1 animal per group
variants.per.group <- variants.per.group %>% mutate(MatPeptide = case_when(Syn == "synonymous SNV" & num.animals != "1" ~ MatPeptide,
                                                                             Syn == "synonymous SNV" & num.animals == "1" ~ "",
                                                                             Syn != "synonymous SNV" & AF.mean >= 20 ~ MatPeptide,
                                                                             Syn != "synonymous SNV" & AF.mean < 20 ~ ""))
# Custom colors for treatment groups
palette.treatment <- list("Vehicle" = "#000000",
                          "Nirmatrelvir" = "#ff2600",
                          "Molnupiravir" = "#0433FF",
                          "4'-FlU" = "#00F900")

# Custom shapes for the number of animals per group with a given variant
shape.animals <- list("1" = 16,
                      "2" = 17,
                      "3" = 15,
                      "4" = 8,
                      "5" = 18,
                      "6" = 9)

# Split the dataframe by timepoint
variants.per.group.D14 <- variants.per.group %>% filter(timepoint == "D14")
variants.per.group.D21 <- variants.per.group %>% filter(timepoint == "D21")
variants.per.group.D28 <- variants.per.group %>% filter(timepoint == "D28")

####################################################################################
# Plot the table of variants, separated by treatment group and timepoint
# Note, if saving in pdf format, you can include the alpha (transparency on points)
# eds format doesn't accept this transparency change, so use plot commands without
# alpha specification
####################################################################################
# Plot the D14 variants per treatment group
plot.variants.per.group.D14 <- ggplot(variants.per.group.D14, aes(x = treatment, y = AF.mean, color = treatment, label = MatPeptide, shape=num.animals)) + 
  geom_jitter(size = 2, alpha = 0.6, position = position_jitter(seed = 1,width = 0.2, height = 0)) +
  # geom_jitter(size = 2, position = position_jitter(seed = 1,width = 0.2, height = 0)) +
  geom_text_repel(size = 3,  max.overlaps = 50, position = position_jitter(seed = 1,width = 0.2, height = 0), min.segment.length = 0)+
  theme_classic(base_size = 12) + 
  theme(legend.position = c(0.4,0.8))+
  scale_shape_manual(name = "# of animals", values=shape.animals,drop = FALSE) +
  guides(colour = "none") +
  scale_y_continuous(limits = c(0,101), breaks = seq(0,100, by=10)) +
  scale_color_manual(values=palette.treatment) +
  xlab("") +
  ylab("allele frequency (%)") +
  labs(title="") + 
  theme(legend.text = element_text(size = 7), legend.title = element_text(size = 7))

# Look at the plot
# plot.variants.per.group.D14 

# Plot the D21 variants per treatment group
plot.variants.per.group.D21 <- ggplot(variants.per.group.D21, aes(x = treatment, y = AF.mean, color = treatment, label = MatPeptide, shape=num.animals)) + 
  geom_jitter(size = 2, alpha = 0.6, position = position_jitter(seed = 1,width = 0.2, height = 0)) +
  # geom_jitter(size = 2, position = position_jitter(seed = 1,width = 0.2, height = 0)) +
  geom_text_repel(size = 3,  max.overlaps = 50, position = position_jitter(seed = 1,width = 0.2, height = 0), min.segment.length = 0)+
  theme_classic(base_size = 12) + 
  theme(legend.position = "none",axis.text.y = element_blank()) +
  scale_y_continuous(limits = c(0,101), breaks = seq(0,100, by=10)) +
  scale_color_manual(name = "treatment", values=palette.treatment) +
  scale_shape_manual(values=shape.animals,drop = FALSE) +
  xlab("") +
  ylab("") +
  labs(title="")

# Look at the plot
# plot.variants.per.group.D21

# Plot the D28 variants per treatment group
plot.variants.per.group.D28 <- ggplot(variants.per.group.D28, aes(x = treatment, y = AF.mean, color = treatment, label = MatPeptide, shape=num.animals)) + 
  geom_jitter(size = 2, alpha = 0.6, position = position_jitter(seed = 1,width = 0.2, height = 0)) +
  # geom_jitter(size = 2, position = position_jitter(seed = 1,width = 0.2, height = 0)) +
  geom_text_repel(size = 3,  max.overlaps = 50, position = position_jitter(seed = 1,width = 0.2, height = 0), min.segment.length = 0)+
  theme_classic(base_size = 12) + 
  theme(legend.position = "none", axis.text.y = element_blank()) +
  scale_y_continuous(limits = c(0,101), breaks = seq(0,100, by=10)) +
  scale_color_manual(values=palette.treatment) +
  scale_shape_manual(values=shape.animals) +
  xlab("") +
  ylab("") +
  labs(title="") 

# Look at the plot
# plot.variants.per.group.D28 

# Combine the variant scatter plots into a panel including all timepoints
panel <- plot_grid(plot.variants.per.group.D14, plot.variants.per.group.D21, plot.variants.per.group.D28, nrow = 1, rel_widths = c(1.3,3.5,2.2)) + theme(text=element_text(family="Arial"))

# Look at the panel
# panel 
# ggsave("panel_v1.pdf", plot=panel, width = 15, height = 7)

####################################################################################
# Plot venn diagrams of mean variants per animal, where size of circle corresponds
# to mean variants per animal
####################################################################################
# Plot the mean variants per animal within a treatment group at D14
circles.D14 <- ggplot(data = num.variants.per.animal.mean.D14, aes(x=treatment, y=-1, size = circle.size, fill = treatment, color = treatment, label = labels)) +
  geom_point(alpha = 0.3) +
  # geom_point() +
  scale_color_manual(values=palette.treatment) +
  geom_text( size = 3, color = "black", nudge_y = 27) +
  scale_y_continuous(limits = c(-12,40)) +
  scale_size(range = c(3.7,14.5), limits = c(min(num.variants.per.animal.mean$circle.size),max(num.variants.per.animal.mean$circle.size))) +
  theme_void() +
  xlab("") +
  labs(title = "Day 14") +
  ylab("mut per animal") +
  theme(legend.position = "NA",axis.text.y = element_blank(), axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5), axis.title.y=element_text(angle=90,size=10), axis.line.y = element_line(color = "white"), axis.text.y.left = element_text(color = "white")) 

# Look at the plot
# circles.D14

# Plot the mean variants per animal within a treatment group at D21
circles.D21 <- ggplot(data = num.variants.per.animal.mean.D21, aes(x=treatment, y=-1, size = circle.size, fill = treatment, color = treatment, label = labels)) +
  geom_point(alpha = 0.3) +
  # geom_point() +
  scale_color_manual(values=palette.treatment) +
  geom_text( size = 3, color = "black", nudge_y = 27) +
  scale_y_continuous(limits = c(-12,40)) +
  scale_size(range = c(3.7,14.5), limits = c(min(num.variants.per.animal.mean$circle.size),max(num.variants.per.animal.mean$circle.size))) +
  theme_classic() +
  xlab("") +
  ylab("") +
  labs(title = "Day 21") +
  theme_void()+
  theme(legend.position = "NA",axis.text.y = element_blank(), axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5),axis.line.y = element_line(color = "white"), axis.text.y.left = element_text(color = "white")) 

# look at the plot
# circles.D21

# Plot the mean variants per animal within a treatment group at D21
circles.D28 <- ggplot(data = num.variants.per.animal.mean.D28, aes(x=treatment, y=-1, size = circle.size, fill = treatment, color = treatment, label = labels)) +
  geom_point(alpha = 0.3) +
  # geom_point() +
  scale_color_manual(values=palette.treatment) +
  geom_text( size = 3, color = "black", nudge_y =27) +
  scale_y_continuous(limits = c(-12,40)) +
  scale_size(range = c(3.7,14.5), limits = c(min(num.variants.per.animal.mean$circle.size),max(num.variants.per.animal.mean$circle.size))) +
  theme_classic() +
  xlab("") +
  ylab("") +
  labs(title = "Day 28") +
  theme_void()+
  theme(legend.position = "NA",axis.text.y = element_blank(), axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5),axis.line.y = element_line(color = "white"), axis.text.y.left = element_text(color = "white")) 

# look at the plot
# circles.D28

# Combine the mean variants per animal venn diagrams into a panel including all timepoints
panel2 <- plot_grid(circles.D14, circles.D21, circles.D28, nrow = 1, rel_widths = c(1.3,3.5,2.2)) + theme(text=element_text(family="Arial"))

# look at the panel
# panel2

####################################################################################
# Combine the mean variants per animal venn diagrams with the 
# variant scatter plots into a single figure :) 
####################################################################################
# Combine both panels 
panel.final <- plot_grid(panel2, panel, nrow = 2, ncol = 1, rel_heights = c(1,4))

# Look at the final panel
# panel.final

# Save the combined panel
# ggsave(filename = "figure3H_NGS.eps", plot = panel.final, width = 10, height = 6)
ggsave(filename = "figure3H_NGS.pdf", plot = panel.final, width = 10, height = 6)
