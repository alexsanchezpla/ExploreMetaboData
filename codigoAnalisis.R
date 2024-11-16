## ----setup, echo=FALSE-----------------------------------------------------
library(knitr)

## Global options
options(max.print="75")
opts_chunk$set(echo=TRUE,
	             cache=FALSE,
               prompt=FALSE,
               tidy=TRUE,
               comment=NA,
               message=FALSE,
               warning=FALSE)

opts_knit$set(width=75)


## ----selectRep-------------------------------------------------------------
set.seed(123)
sample(1:6,1)


## ----readFromURL-----------------------------------------------------------
metadatosURL <- "https://raw.githubusercontent.com/nutrimetabolomics/metaboData/refs/heads/main/Datasets/2018-MetabotypingPaper/DataInfo_S013.csv"
datosURL <- "https://raw.githubusercontent.com/nutrimetabolomics/metaboData/refs/heads/main/Datasets/2018-MetabotypingPaper/DataValues_S013.csv"


## --------------------------------------------------------------------------
metaDatos <- read.csv(metadatosURL, row.names = 1)
datos <- read.csv(datosURL, row.names = 1)
class(datos)


## ----splitData-------------------------------------------------------------
infoSamples <- datos[,1:5]
if (ncol(datos) ==695)
  datos <- datos[,-c(1:5)]
if (nrow(metaDatos) ==695)
  metaDatos <- metaDatos[-c(1:5),]


## ----fig.align='center', out.width="80%"-----------------------------------
knitr::include_graphics("images/clipboard-388301461.png")


## ----installSE-------------------------------------------------------------
if(!require(BiocManager)) 
  install.packages("BiocManager")
if(!require(SummarizedExperiment)) 
  BiocManager.install("SummarizedExperiment")


## ----dataprepare-----------------------------------------------------------
library(dplyr)

# Convertir las columnas 2 (SURGERY), 4 (GENDER) y 5 (Group) a factores
infoSamples <- infoSamples %>%
  mutate(
    SURGERY = as.factor(SURGERY),
    GENDER = as.factor(GENDER),
    Group = as.factor(Group)
  )


# Crear el nuevo nombre para las filas
infoSamples <- infoSamples %>%
  mutate(
    RowName = paste0(toupper(substr(SURGERY, 1, 1))
                     , "_", Group, "_", SUBJECTS)
  )
rownames(infoSamples) <- infoSamples$RowName
infoSamples <- infoSamples %>% select(-RowName)
head(infoSamples)


## --------------------------------------------------------------------------
matDatos <- as.matrix(datos)
if (sum(rownames(matDatos)!=infoSamples$SUBJECTS)==0)
  rownames(matDatos) <- rownames(infoSamples)



## ----createSE--------------------------------------------------------------
mySE <- SummarizedExperiment(assays = 
            list (rawValues= matDatos),
            rowData=infoSamples, 
            colData=metaDatos)
show(mySE)


## ----metaDatos-------------------------------------------------------------
data <- rowData(mySE)[, c("SURGERY", "AGE", "GENDER", "Group")]

# Convertir a data.frame para facilidad de manipulación
data <- as.data.frame(data)


## ----describbMetaDatos-----------------------------------------------------
# Calcular tablas de frecuencias para las columnas 1 (SURGERY), 3 (GENDER) y 4 (Group)
(freq_surgery <- table(data$SURGERY))
(freq_gender <- table(data$GENDER))
(freq_group <- table(data$Group))
(freq_surgery_Group <- table(data$SURGERY,data$Group))
(freq_surgery_Group_Gender <- table(data$SURGERY,data$Group, data$GENDER))
# Calcular el resumen numérico para la columna 2 (AGE)
(summary_age <- summary(data$AGE))



## ----naniarPackage---------------------------------------------------------
if (!require(naniar)) install.packages("naniar")
library(naniar)


## ----missings--------------------------------------------------------------
# Visualización de valores faltantes
dataf <- as.data.frame(t(assay(mySE, "rawValues")))
# rownames(dataf) <- rownames(rawVals)
# colnames(dataf) <- colnames(rawVals)


## ----plotMissings----------------------------------------------------------
library(ggplot2)

# Visualización con ajuste del tamaño de etiquetas
vis_miss(dataf) +
  theme(axis.text.x = element_text(size = 6, angle = 90, hjust = 1))


## --------------------------------------------------------------------------
n_miss(dataf)
n_complete(dataf)
prop_miss(dataf)
prop_complete(dataf)
pct_miss(dataf)
pct_complete(dataf)


## ----imputeMissings--------------------------------------------------------
rawValues <- assay(mySE, "rawValues")
rawNoMissings <- rawValues
rawNoMissings[is.na(rawNoMissings)] <- 1
sum(is.na(rawNoMissings))
assays(mySE)$rawNoMissings <- rawNoMissings
show(mySE)


## ----boxplot1--------------------------------------------------------------
boxplot(assay(mySE, "rawNoMissings"),
        xlab = "Concentraciones",
        cex.lab=0.8,
        horizontal=TRUE,
        cex.axis=0.4, 
        las=2, 
        main="Distribucion de los valores por metabolito", cex.main=.8)


## ----boxplot2--------------------------------------------------------------
boxplot(t(assay(mySE, "rawNoMissings")),
        ylab="Muestra",
        xlab = "Concentraciones",
        cex.lab=0.8,
        horizontal=TRUE,
        cex.axis=.6, las=2, 
        main="Distribucion de los valores por muestra", cex.main=.8)


## ----logNoMissings---------------------------------------------------------
logNoMissings <- log(assay(mySE, "rawNoMissings"))
assays(mySE)$logNoMissings <- logNoMissings
show(mySE)


## --------------------------------------------------------------------------
save(mySE, file="metabodatSE.Rda")


## ----boxplot1log-----------------------------------------------------------
boxplot(assay(mySE, "logNoMissings"),
        xlab = "logConcentraciones",
        cex.lab=0.8,
        horizontal=TRUE,
        cex.axis=.4, las=2, 
        main="Distribucion de los logvalores por metabolito", cex.main=.8)



## ----boxplot2b-------------------------------------------------------------
boxplot(t(assay(mySE, "logNoMissings")),
        ylab="Muestra",
        xlab = "log-Concentraciones",
        cex.lab=0.8,
        horizontal=TRUE,
        cex.axis=.6, las=2, 
        main="Distribucion de los log-valores por muestra", cex.main=.8)


## ----pheatmap--------------------------------------------------------------
if(!require(pheatmap)) 
  install.packages("pheatmap")
library(pheatmap)


## ----plotHeatmap1----------------------------------------------------------

# Heatmap sin clustering
pheatmap(t(assay(mySE, "rawNoMissings")),
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         main = "Heatmap sin Clustering",
    fontsize_row = 3,         
         fontsize_col = 6,
         color = colorRampPalette(c("blue", "white", "red"))(100))


## ----plotHeatmap2----------------------------------------------------------
# Heatmap sin clustering
pheatmap(t(assay(mySE, "rawValues")),
         cluster_rows = FALSE, 
         cluster_cols = TRUE, 
         main = "Heatmap sin Clustering",
         fontsize_row = 3,  
         fontsize_col = 6,
         color = colorRampPalette(c("blue", "white", "red"))(100))


## ----installPCAtools-------------------------------------------------------
if (!require(PCAtools)) BiocManager::install("PCAtools", dep=TRUE)


## ----runPCA----------------------------------------------------------------
library(PCAtools)
x <- t(assay(mySE, "rawNoMissings"))
p <- pca(x, 
         metadata = rowData(mySE), 
         removeVar = 0.1)
class(p)


## ----plotScreeplot---------------------------------------------------------
screeplot(p, axisLabSize = 10, titleLabSize = 22)


## ----plotBiplotxSURGERY----------------------------------------------------
biplot(p,
       showLoadings = FALSE,
       labSize = 3, 
       pointSize = 3,
       axisLabSize = 8,
       title = "PCA de los datos",
       sizeLoadingsNames = 3,
       colby = 'SURGERY',
       hline = 0, vline = 0,
       legendPosition = 'right')


## ----plotBiplotxSURGERYandEllipse------------------------------------------
biplot(p,
       showLoadings = FALSE,
       labSize = 3, 
       pointSize = 3,
       axisLabSize = 8,
       ellipse=TRUE,
       title = "PCA de los datos",
       sizeLoadingsNames = 3,
       colby = 'SURGERY',
       hline = 0, vline = 0,
       legendPosition = 'right')


## ----plotBiplotxGroup------------------------------------------------------
biplot(p,
       showLoadings = FALSE,
       labSize = 3, 
       pointSize = 3,
       axisLabSize = 8,
       title = "PCA de los datos",
       sizeLoadingsNames = 3,
       colby = 'Group',
       hline = 0, 
       vline = 0,
       legendPosition = 'right')


## ----plotBiplotLoadings----------------------------------------------------
biplot(p,
       showLoadings = TRUE,
       labSize = 3, 
       pointSize = 3,
       axisLabSize = 8,
       title = "PCA de los datos",
       sizeLoadingsNames = 3,
       colby = 'SURGERY',
       lab= NULL,
       hline = 0, vline = 0,
       legendPosition = 'right')


## --------------------------------------------------------------------------
res.hc <- assay(mySE, "rawNoMissings") %>%
  scale() %>%                    # Scale the data
  dist(method = "euclidean") %>% # Compute dissimilarity matrix
  hclust(method = "ward.D2")     # Compute hierachical clustering

library(factoextra)
fviz_dend(res.hc, cex=0.6,
          k = 2,
          k_colors = c("#2E9FDF",  "#FC4E07"),
          color_labels_by_k = TRUE, # color labels by groups
          rect = TRUE) # Add rectangle around groups)

distMat <- dist(assay(mySE, "rawNoMissings"))
hc <-hclust(distMat) 
plot(hclust(distMat))


## ----preparaRepo-----------------------------------------------------------
# Definir nombres de archivos
informe_origen <- "UOC-MU-AD0-2024-25-S2-PEC1-Solucion.pdf"
informe_destino <- "InformeAnalisis.pdf"

codigo_origen <- "UOC-MU-AD0-2024-25-S2-PEC1-Solucion.R"
codigo_destino <- "codigoAnalisis.R"

objeto_binario <- "metaboDatSE.Rda"
metadatos_destino <- "metaDatos.md"

# Realizar las copias de los archivos
file.copy(informe_origen, informe_destino, overwrite = TRUE)
file.copy(codigo_origen, codigo_destino, overwrite = TRUE)
# Supongamos que ya has generado los metadatos como un archivo "metaDatos.md"

# Crear un nuevo directorio para el repositorio
repo_name <- "ExploreMetaboData"
repo_path <- file.path(getwd(), repo_name)
dir.create(repo_path, showWarnings = FALSE)

# Mover los archivos al nuevo directorio
file.copy(informe_destino, file.path(repo_path, informe_destino), overwrite = TRUE)
file.copy(codigo_destino, file.path(repo_path, codigo_destino), overwrite = TRUE)
file.copy(objeto_binario, file.path(repo_path, objeto_binario), overwrite = TRUE)
file.copy(metadatos_destino, file.path(repo_path, metadatos_destino), overwrite = TRUE)


