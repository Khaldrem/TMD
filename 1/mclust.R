library("ggpubr")
library("corrplot")
library("psych")
library("mclust")
library("corrplot")

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

GetMeasuresDisp <- function(x) {
  print(sprintf("Var: %f", var(x)))
  print(sprintf("Sd: %f", sd(x)))
  print(sprintf("Mean: %f", mean(x)))
  print(sprintf("Median: %i", median(x)))
  print(sprintf("Mode: %i", Mode(x)))
}

# --------- First approach -------------
#Working dir
setwd("~/Desktop/code/TMD")

# Dataset path
dataset_path <- "./dataset/breast-cancer-wisconsin.data"
data <- read.csv(dataset_path, header = FALSE)

# Rename variables
names(data) <- c("id", "clump_thickness", "uniformity_cell_size", "uniformity_cell_shape", "marginal_adhesion", "single_epithelial_cell_size", "bare_nuclei", "bland_chromatin", "normal_nucleoli", "mitoses", "class")

# Get first 3 samples
head(data, n = 3)

# Summary data
summary(data)

# --------- Cleanning data -------------

# Dropping ID variable
data$id <- NULL

# Removing rows that contains '?' in 'bare_nuclei' column
data <- data[-which(data$bare_nuclei == "?"),]

# Converting bare_nuclei into int
data$bare_nuclei <- as.integer(data$bare_nuclei)


# --------- Measures of dispersion -------------
GetMeasuresDisp(data$clump_thickness)
GetMeasuresDisp(data$uniformity_cell_size)
GetMeasuresDisp(data$uniformity_cell_shape)
GetMeasuresDisp(data$marginal_adhesion)
GetMeasuresDisp(data$single_epithelial_cell_size)
GetMeasuresDisp(data$bare_nuclei)
GetMeasuresDisp(data$bland_chromatin)
GetMeasuresDisp(data$normal_nucleoli)
GetMeasuresDisp(data$mitoses)

# --------- Graphs -----------------
# class: (2) benign, (4) malignant

type_table <- table(data$class)
print(type_table)

type_percentage <- 100*prop.table(table(data$class))
print(type_percentage)

type_df <- data.frame(class=c("benign", "malignant"), n_class=c(type_table[1], type_table[2]))

# Barplot of type of 
ggbarplot(type_df, x = "class", y = "n_class", xlab=c("Class"), ylab="# of samples", 
          fill="class", color="class", palette = c("#f6f578", "#06623b"), width = c(0.4, 0.4), 
          title = "Breast cancer type", label = TRUE, label.pos = "out")

#clump_thickness: 1 - 10
gghistogram(data, x = "clump_thickness", y = "..count..", ylab = "# of samples", fill = "#E7B800", size = 0.2,
            add = "mean", rug = TRUE, add_density = TRUE, bins = 10, title = "Clump Thickness")

#uniformity_cell_size: 1 - 10
gghistogram(data, x = "uniformity_cell_size", y = "..count..", ylab = "# of samples", fill = "#e74d00", size = 0.2,
            add = "mean", rug = TRUE, add_density = TRUE, bins = 10, title = "Uniformity of Cell Size")

#uniformity_cell_shape: 1 - 10
gghistogram(data, x = "uniformity_cell_shape", y = "..count..", ylab = "# of samples", fill = "#e70f00", size = 0.2,
            add = "mean", rug = TRUE, add_density = TRUE, bins = 10, title = "Uniformity of Cell Shape")

#marginal_adhesion: 1 - 10
gghistogram(data, x = "marginal_adhesion", y = "..count..", ylab = "# of samples", fill = "#00e79a", size = 0.2,
            add = "mean", rug = TRUE, add_density = TRUE, bins = 10, title = "Marginal Adhesion")

#single_epithelial_cell_size: 1 - 10
gghistogram(data, x = "single_epithelial_cell_size", y = "..count..", ylab = "# of samples", fill = "#3128b5", size = 0.2,
            add = "mean", rug = TRUE, add_density = TRUE, bins = 10, title = "Single Epithelial Cell Size")

#bare_nuclei: 1 - 10
gghistogram(data, x = "bare_nuclei", y = "..count..", ylab = "# of samples", fill = "#9d28b5", size = 0.2,
            add = "mean", rug = TRUE, add_density = TRUE, bins = 10, title = "Bare Nuclei")

#bland_chromatin: 1 - 10
gghistogram(data, x = "bland_chromatin", y = "..count..", ylab = "# of samples", fill = "#e3d536", size = 0.2,
            add = "mean", rug = TRUE, add_density = TRUE, bins = 10, title = "Bland Chromatin")

#normal_nucleoli: 1 - 10
gghistogram(data, x = "normal_nucleoli", y = "..count..", ylab = "# of samples", fill = "#9536e3", size = 0.2,
            add = "mean", rug = TRUE, add_density = TRUE, bins = 10, title = "Normal Nucleoli")

#mitoses: 1 - 10
gghistogram(data, x = "mitoses", y = "..count..", ylab = "# of samples", fill = "#e33667", size = 0.2,
            add = "mean", rug = TRUE, add_density = TRUE, bins = 10, title = "Mitoses")


# scatter matrix
pairs.panels(data[,-10], method = "pearson", density = TRUE, ellipses = FALSE)

clPairs(data[,1:9], data$class)

# Correlation matrix
cor_mat <- cor(data[,1:9])
corrplot(cor_mat, method="color")

# Covariance matrix
cov_mat <- cov(data[,1:9])
print(cov_mat)

#Test de normalidad
shapiro.test(data$clump_thickness)
shapiro.test(data$uniformity_cell_size)
shapiro.test(data$uniformity_cell_shape)
shapiro.test(data$marginal_adhesion)
shapiro.test(data$single_epithelial_cell_size)
shapiro.test(data$bare_nuclei)
shapiro.test(data$bland_chromatin)
shapiro.test(data$normal_nucleoli)
shapiro.test(data$mitoses)
shapiro.test(data$class)

# Point-biseral correlation
cor.test(data$clump_thickness, data$class)
cor.test(data$uniformity_cell_size, data$class)
cor.test(data$uniformity_cell_shape, data$class)
cor.test(data$marginal_adhesion, data$class)
cor.test(data$single_epithelial_cell_size, data$class)
cor.test(data$bare_nuclei, data$class)
cor.test(data$bland_chromatin, data$class)
cor.test(data$normal_nucleoli, data$class)
cor.test(data$mitoses, data$class)

#--------------------------------------------------------
# mclust & BIC
#default
mclust_default <- Mclust(data[,1:9])
summary(mclust_default, parameters = "TRUE")

plot(mclust_default, what = "uncertainty")


sort(mclust_default$uncertainty, decreasing = TRUE) %>% head()

legend_args <- list(x = "bottomright", ncol = 5)
plot(mclust_default, what = 'BIC', legendArgs = legend_args)
plot(mclust_default, what = "classification")

legend("topright", legend = 1:2,
       col = mclust.options("classPlotColors"),
       pch = mclust.options("classPlotSymbols"),title = "Class labels:")


mclust_default_density <- densityMclust(data[,1:9])
summary(mclust_default_density)

plot(mclust_default_density, what = "density")

classes <- as.character(data$class)
classes[classes == "2"] <- "benign"
classes[classes == "4"] <- "malignant"

table(classes, mclust_default$classification)

#EXPERIMENTO 1
#BIC
BIC = mclustBIC(data[, 1:9],  prior = priorControl(functionName="defaultPrior", shrinkage=0.1))
summary(BIC)
plot(BIC)

mclust_bic_based <- Mclust(data[,1:9], x=BIC)
summary(mclust_bic_based)
plot(mclust_bic_based, what = "classification")

legend("topright", legend = 1:5,
       col = mclust.options("classPlotColors"),
       pch = mclust.options("classPlotSymbols"),title = "Class labels:")

table(classes, mclust_bic_based$classification)


#ICL
ICL <- mclust::mclustICL(data[,1:9])
summary(ICL)
plot(ICL)

mclust_icl_based = Mclust(data[,1:9], G=2, prior = priorControl(functionName="defaultPrior", shrinkage=0.1), modelNames ="VEV")
plot(mclust_icl_based, what = "classification")
legend("topright", legend = 1:2,
       col = mclust.options("classPlotColors"),
       pch = mclust.options("classPlotSymbols"),title = "Class labels:")


table(classes, mclust_icl_based$classification)

#EXPERIMENTO 2
#Removiendo la variable menos importante mitoses
#BIC
data_wo_mitoses <- data
data_wo_mitoses$mitoses <- NULL

BIC_wo_mitoses = mclustBIC(data_wo_mitoses[, 1:8],  prior = priorControl(functionName="defaultPrior", shrinkage=0.1))
summary(BIC_wo_mitoses)
plot(BIC_wo_mitoses)

mclust_bic_wo_mitoses <- Mclust(data_wo_mitoses[,1:8], x=BIC_wo_mitoses)
summary(mclust_bic_wo_mitoses)
plot(mclust_bic_wo_mitoses, what = "classification")

legend("topright", legend = 1:4,
       col = mclust.options("classPlotColors"),
       pch = mclust.options("classPlotSymbols"),title = "Class labels:")

table(classes, mclust_bic_wo_mitoses$classification)

#ICL
ICL_wo_mitoses <- mclust::mclustICL(data_wo_mitoses[,1:8])
summary(ICL_wo_mitoses)
plot(ICL_wo_mitoses)

mclust_icl_wo_mitoses = Mclust(data_wo_mitoses[,1:8], G=3, prior = priorControl(functionName="defaultPrior", shrinkage=0.1), modelNames ="VEV")
plot(mclust_icl_wo_mitoses, what = "classification")
legend("topright", legend = 1:3,
       col = mclust.options("classPlotColors"),
       pch = mclust.options("classPlotSymbols"),title = "Class labels:")


table(classes, mclust_icl_wo_mitoses$classification)

#EXPERIMENTO 3
#Eliminando la variable mas significativa Bare Nucleoi
data_wo_barenucleoi <- data
data_wo_barenucleoi$bare_nuclei <- NULL

#BIC
BIC_wo_barenucleoi = mclustBIC(data_wo_barenucleoi[, 1:8],  prior = priorControl(functionName="defaultPrior", shrinkage=0.1))
summary(BIC_wo_barenucleoi)
plot(BIC_wo_barenucleoi)

mclust_bic_wo_barenucleoi <- Mclust(data_wo_barenucleoi[,1:8], x=BIC_wo_barenucleoi)
summary(mclust_bic_wo_barenucleoi)
plot(mclust_bic_wo_barenucleoi, what = "classification")

legend("topright", legend = 1:3,
       col = mclust.options("classPlotColors"),
       pch = mclust.options("classPlotSymbols"),title = "Class labels:")

table(classes, mclust_bic_wo_barenucleoi$classification)

#ICL
ICL_wo_barenucleoi <- mclust::mclustICL(data_wo_barenucleoi[,1:8])
summary(ICL_wo_barenucleoi)
plot(ICL_wo_barenucleoi)

mclust_icl_wo_barenucleoi = Mclust(data_wo_barenucleoi[,1:8], G=5, prior = priorControl(functionName="defaultPrior", shrinkage=0.1), modelNames ="VEI")
plot(mclust_icl_wo_barenucleoi, what = "classification")
legend("topright", legend = 1:5,
       col = mclust.options("classPlotColors"),
       pch = mclust.options("classPlotSymbols"),title = "Class labels:")


table(classes, mclust_icl_wo_barenucleoi$classification)