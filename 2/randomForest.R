library("randomForest")
library("dplyr")
library("ROCR")
library("reshape2")
library("ggplot2")
library("imbalance")
library("caret")
library("pROC")
library("ROSE")

#----------- Functions --------------
generateRandomForestModels <- function(data, ntrees) {
  models <- list()
  for (i in 1:(ncol(data)-1)) {
    for(j in 1:length(ntrees)) {
      model <- randomForest(class ~., data = data, importance = TRUE, proximity = TRUE, ntree=ntrees[j], mtry=i)
      models <- c(models, list(model))
    }
  }
  return(models)
}

getNTreeByModel <- function(models) {
  data <- list()
  for (i in 1:length(models)) {
    data <- c(data, models[[i]]$ntree)
  }
  return(data)
}

getMtryByModel <- function(models) {
  data <- list()
  for (i in 1:length(models)) {
    data <- c(data, models[[i]]$mtry)
  }
  return(data)
}

getOOBByModel <- function(models) {
  data <- list()
  for (i in 1:length(models)) {
    data <- c(data, models[[i]]$err.rate[models[[i]]$ntree, 1])
  }
  return(data)
}

saveOOBPlot <- function(path, models, exp) {
  for (i in 1:length(models)) {
    ntree = models[[i]]$ntree
    mtry = models[[i]]$mtry
    
    #Create dataframe to plot
    oob <- models[[i]]$err.rate[,1]
    benign <- models[[i]]$err.rate[,2]
    malignant <- models[[i]]$err.rate[,3]
    
    model_df <- data.frame(index=(1:ntree), oob, benign, malignant)
    df2plot <- melt(model_df ,  id.vars = 'index', variable.name = 'series')
    
    ggplot(df2plot, aes(index, value)) +
      geom_line(linetype = "dashed", aes(colour = series), size=0.3) +
      ggtitle(sprintf("Gráfico Error OOB Exp. %s - (%i, %i)", exp, ntree, mtry)) +
      xlab("Número de árboles") + 
      ylab("%") +
      theme_minimal()
    
    ggsave(sprintf("%s/oob_plot_%i_%i.png", path, ntree, mtry), width=1920, height=1080, units = "px")
  }
}

saveVarImpPlot <- function(path, models, exp) {
  for (i in 1:length(models)) {
    ntree = models[[i]]$ntree
    mtry = models[[i]]$mtry
    
    png(filename=sprintf("%s/varimpplot_%i_%i.png", path, ntree, mtry), width=1024, height=700,units = "px")
    varImpPlot(models[[i]], main=sprintf("Exp. %s - (%i, %i)", exp, ntree, mtry))
    dev.off()
  }
}

saveROCPlot <- function(path, models, dataset, exp) {
  for (i in 1:length(models)) {
    ntree = models[[i]]$ntree
    mtry = models[[i]]$mtry
    
    predict_df <- as.data.frame(predict(models[[i]], dataset, type = "prob"))
    predict_df$predict <- names(predict_df)[1:2][apply(predict_df[,1:2], 1, which.max)]
    predict_df$observed <- dataset$class
    
    roc.benign <- roc(ifelse(predict_df$observed=="benign", "benign", "non-benign"), as.numeric(predict_df$benign))
    roc.malignant <- roc(ifelse(predict_df$observed=="malignant", "malignant", "non-malignant"), as.numeric(predict_df$malignant))
    
    png(filename=sprintf("%s/roc_%i_%i.png", path, ntree, mtry), width=1024, height=700,units = "px")
    plot(roc.benign, col = "green")
    lines(roc.malignant, col = "red")
    legend("bottomright", inset = 0.1, legend=c("benign", "malignant"),
           col=c("green", "red"), lty=1:2, cex=0.7)
    dev.off()
  }
  
  plot_list = list()
  for(i in 1:length(models)) {
    predict_df_1 <- as.data.frame(predict(models[[i]], dataset, type = "prob"))
    pred_test_roc <- roc(dataset$class, predict_df_1$benign)
    
    p =  ggroc(pred_test_roc, lwd=0.5, col="blue")+
         geom_abline(intercept = 1, slope = 1, color = "red", linetype = "dashed", lwd=0.5)
    
    plot_list[[i]] = p
  }
  
  for (i in 1:length(models)) {
    ntree = models[[i]]$ntree
    mtry = models[[i]]$mtry
    
    file_name = sprintf("%s/roc_all_%i_%i.png", path, ntree, mtry)
    png(file_name)
    print(plot_list[[i]])
    dev.off()
  }
}

saveConfusionMatrix <- function(models, dataset, path) {
  cm_df <- data_frame()
  for (i in 1:length(models)) {
    pred_test_class <- predict(models[[i]], dataset, type="class", norm.votes=TRUE, predict.all=FALSE, proximity=FALSE, nodes=FALSE)
    results <- confusionMatrix(pred_test_class, dataset$class, positive="benign", mode="everything")
    
    cm_df_append <- data.frame(cbind(t(round(results$overall, 4)),t(round(results$byClass, 4))))
    cm_df <- rbind(cm_df, cm_df_append)
  }
  
  cm_df["nro_model"] <- c(1:length(models))
  write.csv(cm_df, file=path)
}


# --------- Loading data -------------
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

# Converting class to string
data[which(data$class == 2),10] <- "benign"
data[which(data$class == 4),10] <- "malignant"
data$class <- as.factor(data$class)

#----------- Preparing data -------------
set.seed(505)
list_of_ntrees <- c(500, 1000, 1500)
#Para Exp. 3


#----------------------------------------
#Experimento 1: randomForest con clases desbalanceadas
#Para Exp. 1
ind <- sample(2, nrow(data), replace = TRUE,prob = c(0.7, 0.3))
data_train_1 <- data[ind==1,]
data_test_1 <- data[ind==2,]

#Imbalance ratio (Size minority class)/(Size majority class)
imbalanceRatio(data, classAttr = "class")

#Experimento 1.1: Con todas las variables

#Observamos que la clase esta desbalanceada
#0.65% para bening y 0.35% para malignant
table(data_train_1$class)
prop.table(table(data_train_1$class))

#Generamos la combinatoria de modelos para 500, 1000 y 1500 arboles
models <- generateRandomForestModels(data_train_1, list_of_ntrees)

ntree_by_model <- getNTreeByModel(models)
mtry_by_model <- getMtryByModel(models)
oob_by_model <- getOOBByModel(models)

#Creamos dataframe con los resultados
exp_1_1_df <- data_frame(nroModel=1:length(models), ntree=as.numeric(ntree_by_model), mtry=as.numeric(mtry_by_model), oob=round(as.double(oob_by_model), 4))

write.csv(exp_1_1_df, file="~/Desktop/code/TMD/2/results/exp_1_1.csv")

exp_1_1_oob_plot_path <- "~/Desktop/code/TMD/2/images/1/1/oob"
saveOOBPlot(exp_1_1_oob_plot_path, models, "1.1")

#dotchart variable importance
exp_1_1_varimpplot_path <- "~/Desktop/code/TMD/2/images/1/1/varimpplot"
saveVarImpPlot(exp_1_1_varimpplot_path, models, "1.1")

#ROC
exp_1_1_roc_path <- "~/Desktop/code/TMD/2/images/1/1/roc"
saveROCPlot(exp_1_1_roc_path, models, data_test_1, "1.1")

#ConfusionMatrix results
cm_path <- "~/Desktop/code/TMD/2/results/cm_exp_1_1.csv"
saveConfusionMatrix(models, data_test_1, cm_path)


#Experimento 1.2: Removiendo variable mas importante

#Removiendo Variable: bare nuclei
data_train_1$bare_nuclei <- NULL
data_test_1$bare_nuclei <- NULL

models <- generateRandomForestModels(data_train_1, list_of_ntrees)
ntree_by_model <- getNTreeByModel(models)
mtry_by_model <- getMtryByModel(models)
oob_by_model <- getOOBByModel(models)

exp_1_2_df <- data_frame(nroModel=1:length(models), ntree=as.numeric(ntree_by_model), mtry=as.numeric(mtry_by_model), oob=round(as.double(oob_by_model), 4))
write.csv(exp_1_2_df, file="~/Desktop/code/TMD/2/results/exp_1_2.csv")

exp_1_2_oob_plot_path <- "~/Desktop/code/TMD/2/images/1/2/oob"
saveOOBPlot(exp_1_2_oob_plot_path, models, "1.2")

#dotchart variable importance
exp_1_2_varimpplot_path <- "~/Desktop/code/TMD/2/images/1/2/varimpplot"
saveVarImpPlot(exp_1_2_varimpplot_path, models, "1.2")

#ROC
exp_1_2_roc_path <- "~/Desktop/code/TMD/2/images/1/2/roc"
saveROCPlot(exp_1_2_roc_path, models, data_test_1, "1.2")

#ConfusionMatrix results
cm_path <- "~/Desktop/code/TMD/2/results/cm_exp_1_2.csv"
saveConfusionMatrix(models, data_test_1, cm_path)

#Experimento 1.3: Removiendo variable menos importante
#Removiendo Variable: mitoses
data_train_1 <- data[ind==1,]
data_test_1 <- data[ind==2,]
data_train_1$mitoses <- NULL
data_test_1$mitoses <- NULL

models <- generateRandomForestModels(data_train_1, list_of_ntrees)
ntree_by_model <- getNTreeByModel(models)
mtry_by_model <- getMtryByModel(models)
oob_by_model <- getOOBByModel(models)

exp_1_3_df <- data_frame(nroModel=1:length(models), ntree=as.numeric(ntree_by_model), mtry=as.numeric(mtry_by_model), oob=round(as.double(oob_by_model), 4))
write.csv(exp_1_3_df, file="~/Desktop/code/TMD/2/results/exp_1_3.csv")

exp_1_3_oob_plot_path <- "~/Desktop/code/TMD/2/images/1/3/oob"
saveOOBPlot(exp_1_3_oob_plot_path, models, "1.3")

#dotchart variable importance
exp_1_3_varimpplot_path <- "~/Desktop/code/TMD/2/images/1/3/varimpplot"
saveVarImpPlot(exp_1_3_varimpplot_path, models, "1.3")

#ROC
exp_1_3_roc_path <- "~/Desktop/code/TMD/2/images/1/3/roc"
saveROCPlot(exp_1_3_roc_path, models, data_test_1, "1.3")

#ConfusionMatrix results
cm_path <- "~/Desktop/code/TMD/2/results/cm_exp_1_3.csv"
saveConfusionMatrix(models, data_test_1, cm_path)


#Experimento 2: randomForest con oversampling
#Para Exp. 2
#MWMOTE - Oversample

numBening <- length(which(data$class == "benign"))
numMalignant <- length(which(data$class == "malignant"))
nInstances <- numBening - numMalignant
mwmote_data <- mwmote(data, nInstances, classAttr = "class")

data_oversampled_mwmote <- rbind(data, mwmote_data)

# Plot a visual comparison between new and old dataset
plotComparison(data, data_oversampled_mwmote, attrs = names(data)[1:3], classAttr = "class")
imbalanceRatio(data_oversampled_mwmote, classAttr = "class")

ind_oversampling <- sample(2, nrow(data_oversampled_mwmote), replace = TRUE,prob = c(0.7, 0.3))
data_train_2 <- data_oversampled_mwmote[ind_oversampling==1,]
data_test_2 <- data_oversampled_mwmote[ind_oversampling==2,]

#Experimento 2.1: Con todas las variables
models_2_1 <- generateRandomForestModels(data_train_2, list_of_ntrees)

ntree_by_model <- getNTreeByModel(models_2_1)
mtry_by_model <- getMtryByModel(models_2_1)
oob_by_model <- getOOBByModel(models_2_1)

#Creamos dataframe con los resultados
exp_2_1_df <- data_frame(nroModel=1:length(models_2_1), ntree=as.numeric(ntree_by_model), 
                        mtry=as.numeric(mtry_by_model), oob=round(as.double(oob_by_model), 4))

write.csv(exp_2_1_df, file="~/Desktop/code/TMD/2/results/exp_2_1.csv")

exp_2_1_oob_plot_path <- "~/Desktop/code/TMD/2/images/2/1/oob"
saveOOBPlot(exp_2_1_oob_plot_path, models_2_1, "2.1")

#dotchart variable importance
exp_2_1_varimpplot_path <- "~/Desktop/code/TMD/2/images/2/1/varimpplot"
saveVarImpPlot(exp_2_1_varimpplot_path, models_2_1, "2.1")

#ROC
exp_2_1_roc_path <- "~/Desktop/code/TMD/2/images/2/1/roc"
saveROCPlot(exp_2_1_roc_path, models_2_1, data_test_2, "2.1")

#ConfusionMatrix results
cm_path <- "~/Desktop/code/TMD/2/results/cm_exp_2_1.csv"
saveConfusionMatrix(models_2_1, data_test_2, cm_path)

#Experimento 2.2: Removiendo variable mas importante
data_train_2$bare_nuclei <- NULL
data_test_2$bare_nuclei <- NULL

models_2_2 <- generateRandomForestModels(data_train_2, list_of_ntrees)
ntree_by_model <- getNTreeByModel(models_2_2)
mtry_by_model <- getMtryByModel(models_2_2)
oob_by_model <- getOOBByModel(models_2_2)

exp_2_2_df <- data_frame(nroModel=1:length(models_2_2), ntree=as.numeric(ntree_by_model), mtry=as.numeric(mtry_by_model), oob=round(as.double(oob_by_model), 4))
write.csv(exp_2_2_df, file="~/Desktop/code/TMD/2/results/exp_2_2.csv")

exp_2_2_oob_plot_path <- "~/Desktop/code/TMD/2/images/2/2/oob"
saveOOBPlot(exp_2_2_oob_plot_path, models_2_2, "2.2")

#dotchart variable importance
exp_2_2_varimpplot_path <- "~/Desktop/code/TMD/2/images/2/2/varimpplot"
saveVarImpPlot(exp_2_2_varimpplot_path, models_2_2, "2.2")

#ROC
exp_2_2_roc_path <- "~/Desktop/code/TMD/2/images/2/2/roc"
saveROCPlot(exp_2_2_roc_path, models_2_2, data_test_2, "2.2")

#ConfusionMatrix results
cm_path <- "~/Desktop/code/TMD/2/results/cm_exp_2_2.csv"
saveConfusionMatrix(models_2_2, data_test_2, cm_path)

#Experimento 2.3: Removiendo variable menos importante
#Removiendo Variable: mitoses
data_train_2 <- data_oversampled_mwmote[ind_oversampling==1,]
data_test_2 <- data_oversampled_mwmote[ind_oversampling==2,]
data_train_2$mitoses <- NULL
data_test_2$mitoses <- NULL

models_2_3 <- generateRandomForestModels(data_train_2, list_of_ntrees)
ntree_by_model <- getNTreeByModel(models_2_3)
mtry_by_model <- getMtryByModel(models_2_3)
oob_by_model <- getOOBByModel(models_2_3)

exp_2_3_df <- data_frame(nroModel=1:length(models_2_3), ntree=as.numeric(ntree_by_model), mtry=as.numeric(mtry_by_model), oob=round(as.double(oob_by_model), 4))
write.csv(exp_2_3_df, file="~/Desktop/code/TMD/2/results/exp_2_3.csv")

exp_2_3_oob_plot_path <- "~/Desktop/code/TMD/2/images/2/3/oob"
saveOOBPlot(exp_2_3_oob_plot_path, models_2_3, "2.3")

#dotchart variable importance
exp_2_3_varimpplot_path <- "~/Desktop/code/TMD/2/images/2/3/varimpplot"
saveVarImpPlot(exp_2_3_varimpplot_path, models_2_3, "2.3")

#ROC
exp_2_3_roc_path <- "~/Desktop/code/TMD/2/images/2/3/roc"
saveROCPlot(exp_2_3_roc_path, models_2_3, data_test_2, "2.3")

#ConfusionMatrix results
cm_path <- "~/Desktop/code/TMD/2/results/cm_exp_2_3.csv"
saveConfusionMatrix(models_2_3, data_test_2, cm_path)


#Experimento 3: randomForest con undersampling
#data_undersampled <- ovun.sample(class ~., data = data, method="under", p=0.5, seed=1443)$data
#ind_under <- sample(2, nrow(data_undersampled), replace = TRUE,prob = c(0.7, 0.3))
#data_train_3 <- data_undersampled[ind_under==1,]
#data_test_3 <- data_undersampled[ind_under==2,]

#imbalanceRatio(data_train_3, classAttr = "class")

#Experimento 3.1: Con todas las variables
#Experimento 3.2: Removiendo variable mas importante
#Experimento 3.3: Removiendo variable menos importante














#---------------------------------------
#General
rf_fullDataset <- randomForest(data_train$class ~., data = data_train[,1:9], importance = TRUE, proximity = TRUE, ntree=1000)

#Results
print(rf_fullDataset)

#Plot
plot(rf_fullDataset)

#Importance
importance(rf_fullDataset)
varImpPlot(rf_fullDataset)

#Proximity
rf_fullDataset.mds <- cmdscale(1 - rf_fullDataset$proximity, eig = TRUE)
op <- par(pty="s")
pairs(cbind(data_train[, 1:9], rf_fullDataset.mds$points), 
            cex = 0.6, gap = 0, col=c("red", "blue")[as.numeric(data_train$class)],
            main = "Breast Cancer Wisconsin: Predictors and MDS of Proximity Based on RandomForest")
par(op)

#MDSplot
print(rf_fullDataset.mds$GOF)
MDSplot(rf_fullDataset, data_train$class)


#OOB por cada config del modelo
print(sprintf("Nro. | ntree | mtry | OBB "))
for (i in 1:length(models)) {
  print(sprintf("%i | %i | %i | %f ", i, models[[i]]$ntree, models[[i]]$mtry, models[[i]]$err.rate[models[[i]]$ntree, 1]))
}


model.mds <- cmdscale(1 - models[[1]]$proximity, eig = TRUE)
print(model.mds)

op <- par(pty="s")
pairs(cbind(data_train_1[, 1:9], model.mds$points), 
      cex = 0.6, gap = 0, col=c("red", "blue")[as.numeric(data_train_1$class)],
      main = "Breast Cancer Wisconsin: Predictors and MDS of Proximity Based on RandomForest")
plot(plotMDS)
legend(1, 95, legend=c("Line 1", "Line 2"),
       col=c("red", "blue"), lty=1:2, cex=0.8)

proximity.plot(models[[1]])

#MDSplot
print(rf_fullDataset.mds$GOF)
MDSplot(rf_fullDataset, data_train$class)
