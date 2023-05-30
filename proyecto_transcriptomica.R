---
title: "Proyecto Transcriptómica"
author: "César Esparza"
date: "29/5/2023"
output: html_document
---

## RNA-seq
```{r librerias_rna-seq, message=FALSE, warning=FALSE}
suppressPackageStartupMessages(library("edgeR"))
suppressPackageStartupMessages(library("ExploreModelMatrix"))
suppressPackageStartupMessages(library("limma"))
suppressPackageStartupMessages(library("recount3"))
suppressPackageStartupMessages(library("clusterProfiler"))
suppressPackageStartupMessages(library("org.Mm.eg.db"))
suppressPackageStartupMessages(library(AnnotationDbi))
```

## Rendimiento
```{r librerias_rendimiento, message=FALSE, warning=FALSE}
suppressPackageStartupMessages(library('stringr'))
suppressPackageStartupMessages(library("biocthis"))
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("postcards"))
suppressPackageStartupMessages(library("pryr"))
suppressPackageStartupMessages(library("usethis"))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library('plyr'))
```
## Tablas
```{r librerias_tablas, message=FALSE, warning=FALSE}
suppressPackageStartupMessages(library("SummarizedExperiment"))
suppressPackageStartupMessages(library("iSEE"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("patchwork"))
suppressPackageStartupMessages(library("pheatmap"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library(ggpubr))
```

## Importe del set de datos 
```{r}
# archivos de salida de featureCounts tienen la extensión .txt y están en el mismo directorio
files <- list.files(path = "/Documentos/cesar/proyrecto_transcriptomica/counts", pattern = "*.txt", full.names = TRUE)
```

```{r}
# Inicializa una lista vacía para almacenar tus datos
data_list <- list()

# Lee cada archivo y añádelo a la lista
for (i in seq_along(files)) {
  data <- read.table(files[i], header = TRUE, row.names = 1)
  # Asegúrate de que solo estás obteniendo la columna de conteo
  data <- data[, "count"]
  data_list[[i]] <- data
}

# Combina todos los datos en una sola matriz de conteo
count_matrix <- do.call(cbind, data_list)

# Asegúrate de que los nombres de las columnas son los nombres de tus muestras
colnames(count_matrix) <- c("SRR5068328", "SRR5068329", "SRR5068330", "SRR5068331", "SRR5068332", "SRR5068333", "SRR5068334", "SRR5068335", "SRR5068336", "SRR5068337", "SRR5068338", "SRR5068339")
```

```{r}
# Crea un data.frame con la información de las condiciones de tus muestras
condition <- factor(c("Mecp2 knockout", "Mecp2 knockout", "Mecp2 knockout", "wildtype", "wildtype", "wildtype", "Mecp2 knockout", "Mecp2 knockout", "Mecp2 knockout", "wildtype", "wildtype", "wildtype"))
coldata <- data.frame(row.names = colnames(count_matrix), condition)
```

```{r}
# Crea un objeto DESeqDataSet a partir de tu matriz de conteo
rse_gene <- DESeqDataSetFromMatrix(countData = count_matrix, colData = coldata, design = ~ condition)

```

```{r}
rse_gene
```

```{r revisar_datos,message=FALSE, warning=FALSE}

# Obtener los numeros de lecturas
assay(rse_gene, "counts") <- compute_read_counts(rse_gene)

# Facilitar el uso de la informacion del experimento
rse_gene <- expand_sra_attributes(rse_gene)

# Explorar los parametros del experimento y algunas de sus variables
colData(rse_gene)[,
    grepl("^sra_attribute", colnames(colData(rse_gene)))]

```
## Formateo


```{r arreglar_atributos}
## Arreglar los atributos para que no generen problemas  
# Remplazar los espacios en blanco para que no hayan errores
rse_gene$sra.sample_attributes <- gsub("week 5", "week_5", rse_gene$sra.sample_attributes)
rse_gene$sra.sample_attributes <- gsub("week 24", "week_24", rse_gene$sra.sample_attributes)


```


```{r mod_data}
# Pasar los atributos a factores para poder manipularlos
rse_gene$sra_attribute.genotype_variation  <- factor(rse_gene$"sra_attribute.genotype/variation")
rse_gene$sra_attribute.Stage <- factor(rse_gene$sra_attribute.Stage)
names(rse_gene) <- gsub("/", "_", names(rse_gene))


# Generar resumen de las variables de interes
summary(as.data.frame(colData(rse_gene)[,
    grepl("^sra_attribute.[cellular_fraction|gender|genotype|source_name|strain|age]", colnames(colData(rse_gene)))]))
```

```{r}
rse_gene$sra_attribute.genotype_variation
```

## Revisar calidad

```{r qual_check1}
## http://rna.recount.bio/docs/quality-check-fields.html
# Analisis de calidad
rse_gene$assigned_gene_prop <- rse_gene$recount_qc.gene_fc_count_all.assigned / rse_gene$recount_qc.gene_fc_count_all.total
summary(rse_gene$assigned_gene_prop)
```

```{r qual_check2}
assignedGeneProp <- data.frame(Categoria = rse_gene$assigned_gene_prop)
fig1 <-ggplot(assignedGeneProp, aes(x=Categoria)) + 
  geom_histogram(color="black", fill= 'darkblue') +
  xlab('Assigned gene prop')
fig1
```

```{r}
rse_gene$sra_attribute.Stage
rse_gene$sra_attribute.genotype_variation
```


```{r attrbt_check}
# Graficar los niveles de expresion en distintas edades
with(colData(rse_gene), plot(assigned_gene_prop, sra_attribute.Stage))
abline(v=0.3585,col = "red")
# Graficar los niveles de expresion en los ratones con y sin RTT
with(colData(rse_gene), plot(assigned_gene_prop, sra_attribute.genotype_variation))
abline(v=0.3585,col = "red")

```

```{r genotyp_check}
# Revisar el atributo de fenotipo
with(colData(rse_gene), tapply(assigned_gene_prop, sra_attribute.genotype_variation, summary))
```

```{r age_check}
# Revisar el atributo de fenotipo
with(colData(rse_gene), tapply(assigned_gene_prop, sra_attribute.Stage, summary))
```


```{r estd_gene}
# Obtener estadísticas de la expresión de genes
gene_means <- rowMeans(assay(rse_gene, "counts"))
summary(gene_means)
```

## Normalización y Expresión diferencial 

```{r }
# Construir un objeto con el cual se podran normalizar los datos
dge <- DGEList(
    counts = assay(rse_gene, "counts"),
    genes = rowData(rse_gene))
dge <- calcNormFactors(dge)
```

```{r}
# Grafica la expresion diferencial entre el control y los Mecp2_KO
ggplot(as.data.frame(colData(rse_gene)), aes(y = assigned_gene_prop, x = sra_attribute.genotype_variation)) +
    geom_boxplot() +
    theme_bw(base_size = 20) +
    ylab("Assigned Gene Prop") +
    xlab("genotype")

```

```{r}
# Graficar la expresion diferencial entre las edades
ggplot(as.data.frame(colData(rse_gene)), aes(y = assigned_gene_prop, x = sra_attribute.Stage)) +
    geom_boxplot() +
    theme_bw(base_size = 20) +
    ylab("Assigned Gene Prop") +
    xlab("Edades")
```


```{r}
# Generar el modelo linear estadistico
mod <- model.matrix(~   sra_attribute.Stage + sra_attribute.genotype_variation + assigned_gene_prop,
    data = colData(rse_gene))
# Observar las variables que componen el modelo
colnames(mod)
```

```{R }
vGene <- voom(dge, mod, plot = TRUE)
```

```{R}
eb_results <- eBayes(lmFit(vGene))
de_results <- topTable(
    eb_results,
    coef = 2,
    number = nrow(rse_gene),
    sort.by = "none"
)
dim(de_results)
```

```{R}
head(de_results)
```

```{R}
table(de_results$adj.P.Val < 0.05)
```

```{R}
## Visualicemos los resultados estadísticos
limma::plotMA(eb_results, coef = 2)
```

```{R}
volcanoplot(eb_results, coef = 2, highlight = 10, names = de_results$gene_name)
```



```{r}
## Extraer valores de los genes de interés
exprs_heatmap <- vGene$E[rank(de_results$adj.P.Val) <= 30, ]
## Creemos una tabla con información de las muestras
## y con nombres de columnas más amigables
df <- as.data.frame(colData(rse_gene)[, 
        c("sra_attribute.Stage","sra_attribute.genotype_variation" )])
colnames(df) <- c("Age", "Genotype")
## Hacer el p heatmap
pheatmap(
    exprs_heatmap,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = FALSE,
    show_colnames = FALSE,
    annotation_col = df
)
```

```{r}
# Obtén la tabla de genes diferencialmente expresados
de_results <- topTable(eb_results, coef = 2, number = nrow(rse_gene), sort.by = "none")

# Filtra los genes con un valor p ajustado < 0.05
significant_genes <- de_results[de_results$adj.P.Val < 0.05, ]

```

```{r}
# Resumen de los resultados
summary(de_results)

# Primeras filas de los resultados
head(de_results)

```


```{r}
# Extrae los nombres de los genes
gene_list <- rownames(significant_genes)

# Elimina la versión del gen de los identificadores
gene_list <- sub("\\..*", "", gene_list)

```

```{r}
# Asegúrate de que 'your_gene_ids' es el vector que contiene tus IDs de genes
ensembl_ids <- mapIds(org.Mm.eg.db,
                      keys = gene_list,
                      column = "ENSEMBL",
                      keytype = "ENSEMBL", # O "SYMBOL", "ENTREZID", etc., dependiendo de tu tipo de ID de genes
                      multiVals = "first")
```

```{r}
# Realiza el análisis de enriquecimiento de GO
ego <- enrichGO(gene         = gene_list,
                OrgDb         = org.Mm.eg.db,
                keyType       = "ENSEMBL",
                ont           = "BP", # Cambia esto a "CC" o "MF" si estás interesado en componentes celulares o funciones moleculares
                pAdjustMethod = "BH",
                qvalueCutoff  = 0.05)

# Muestra los resultados
print(ego)
```

```{r}
barplot(ego, showCategory=10)
dotplot(ego, showCategory=10)

```
