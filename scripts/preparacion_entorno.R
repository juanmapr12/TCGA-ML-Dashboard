# Descomentar estas líneas si no se tiene instalado el paquete 'BiocManager'
# ni la herramienta 'GDCRNATools:
if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
  BiocManager::install('GDCRNATools')
  BiocManager::install('impute')
  BiocManager::install('EnhancedVolcano')
  #BiocManager::install('Limma')
  #BiocManager::install("multiMiR")
  #BiocManager::install("clusterProfiler")
  #BiocManager::install("org.Hs.eg.db")       # Anotación para humano
  #BiocManager::install("enrichplot")         # Para visualizaciones
}
library(GDCRNATools)
library(impute)
library(EnhancedVolcano)
library(limma)

library(moments)
library(ggplot2)
library(dplyr)
library(tidyr)
library(vcd)

library(multiMiR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

library(pheatmap)
library(tidymodels)
library(discrim)  # Para Naive-Bayes



system("wget https://gdc.cancer.gov/files/public/file/gdc-client_v1.6.1_Ubuntu_x64.zip")
unzip("gdc-client_v1.6.1_Ubuntu_x64.zip", exdir=".")

gdcGetURL_new <- function(project.id, data.type) {
  urlAPI <- 'https://api.gdc.cancer.gov/files?'
  if (data.type=='RNAseq') {
    data.category <- 'Transcriptome Profiling'
    data.type <- 'Gene Expression Quantification'
    workflow.type <- 'STAR - Counts'
  } else if (data.type=='miRNAs') {
    data.category <- 'Transcriptome Profiling'
    data.type <- 'Isoform Expression Quantification'
    workflow.type <- 'BCGSC miRNA Profiling'
  } else if (data.type=='Clinical') {
    data.category <- 'Clinical'
    data.type <- 'Clinical Supplement'
    workflow.type <- NA
  } else if (data.type=='pre-miRNAs') {
    data.category <- 'Transcriptome Profiling'
    data.type <- 'miRNA Expression Quantification'
    workflow.type <- 'BCGSC miRNA Profiling'
  }
  
  project <- paste('{"op":"in","content":{"field":"cases.',
                   'project.project_id","value":["', 
                   project.id, '"]}}', sep='')
  dataCategory <- paste('{"op":"in","content":{"field":"files.', 
                        'data_category","value":"', data.category, '"}}', sep='')
  dataType <- paste('{"op":"in","content":{"field":"files.data_type",',
                    '"value":"', data.type, '"}}', sep='')
  workflowType <- paste('{"op":"in","content":{"field":"files.',
                        'analysis.workflow_type","value":"', workflow.type, '"}}', sep='')
  
  
  if (is.na(workflow.type)) {
    dataFormat <- paste('{"op":"in","content":{"field":"files.',
                        'data_format","value":"', 'BCR XML', '"}}', sep='')
    content <- paste(project, dataCategory, dataType, dataFormat, sep=',')
  } else {
    content <- paste(project, dataCategory, dataType, 
                     workflowType, sep=',')
  }
  
  filters <- paste('filters=',URLencode(paste('{"op":"and","content":[', 
                                              content, ']}', sep='')),sep='')
  
  expand <- paste('analysis', 'analysis.input_files', 'associated_entities',
                  'cases', 'cases.diagnoses','cases.diagnoses.treatments', 
                  'cases.demographic', 'cases.project', 'cases.samples', 
                  'cases.samples.portions', 'cases.samples.portions.analytes', 
                  'cases.samples.portions.analytes.aliquots',
                  'cases.samples.portions.slides', sep=',')
  
  expand <- paste('expand=', expand, sep='')
  
  payload <- paste(filters, 'pretty=true', 'format=JSON', 
                   'size=10000', expand, sep='&')
  url <- paste(urlAPI, payload, sep='')
  
  return (url)
}

toolenv <- environment(get("gdcGetURL", envir = asNamespace("GDCRNATools")))
# Sacamos el entorno donde está definida la función 'gdcGetURL' dentro de el
# espacio de nombres de paquete 'GDCmiRNATools'.
# Dentro del entorno de la función puede haber variables y otras funciones que
# hagan referencia a ésta, de ahí que lo necesitemos.
# ¿Qué es un espacio de nombres de paquete? Es un ente que organiza y gestiona
# las funciones, variables y objetos dentro de un paquete.

unlockBinding("gdcGetURL", toolenv)
# Desbloqueamos variable que está bloqueada (por defecto) en el entorno con 
# el fin de modificarla

assignInNamespace("gdcGetURL", gdcGetURL_new, ns="GDCRNATools", envir=toolenv)
# Asignamos en el espacio de nombres o ns 'GDCRNATools' la nueva función
# definida sustituyendo a la antigua. Todo esto dentro del entorno.

assign("gdcGetURL", gdcGetURL_new)
# Asignación global (fuera de todo namespace)

lockBinding("gdcGetURL", toolenv)
# La volvemos a bloquear una vez realizado el cambio.


#### Descarga de datos clínicos: ####

descarga_clinicos <- function(project){
  
  clinical_dir <- paste('metadatos', project, 'Clinical', sep='/')
  
  # Download clinical data
  gdcClinicalDownload(project.id = project,
                      directory = clinical_dir, 
                      write.manifest = FALSE,
                      method = "gdc-client")
}

crea_df_clinico <- function(project){
  clinical_dir <- paste('metadatos', project, 'Clinical', sep='/')
  clinical_df <- gdcClinicalMerge(clinical_dir)
  return(clinical_df)
}


#### Descarga de datos de expresiones de mirnas: ####


descarga_mirnas <- function(project){
  
  isoform_miRNAs_dir <- paste('metadatos', project, 'miRNAs', sep='/')
  
  # Download miRNA data 
  gdcRNADownload(project.id     = project, 
                 data.type      = 'miRNAs', 
                 write.manifest = FALSE,
                 method         = 'gdc-client',
                 directory      = isoform_miRNAs_dir)
}

crea_df_mirnas <- function(project){
  metaMatrix.MIR <- gdcParseMetadata(project.id = project,
                                     data.type  = 'miRNAs', 
                                     write.meta = FALSE)
  return(metaMatrix.MIR)
}



#### Eliminamos carpetas ####

elimina_carpetas_GDC <- function() {
  # Lista todas las carpetas en el WD (no recursivo)
  carpetas <- list.dirs(path = ".", recursive = FALSE, full.names = TRUE)
  
  # Expresión regular para carpetas tipo UUID: xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx
  patron_uuid <- "^[0-9a-fA-F]{8}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{12}$"
  
  # Filtra carpetas que cumplen el patrón
  carpetas_a_borrar <- carpetas[grepl(patron_uuid, basename(carpetas))]
  
  # Borra esas carpetas
  for(carpeta in carpetas_a_borrar){
    # Aseguramos que existe y luego borramos
    if(dir.exists(carpeta)){
      unlink(carpeta, recursive = TRUE, force = TRUE)
    }
  }
  message("Eliminadas ", length(carpetas_a_borrar), " carpetas")
}



