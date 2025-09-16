
# Aquí unimos todos los archivos de cuantificación de isoformas en un solo df
merge_MIR <- function(metadata, fdir){
  filelist <- list.files(fdir, pattern="*.txt$", 
                         recursive = TRUE, full.names=TRUE)
  for (i in 1:length(filelist)){
    iname <- basename(filelist[i])
    isamplename <- metadata[metadata$file_name==iname, "sample"]
    if(length(isamplename) == 0)
      next
    idf <- read.csv(filelist[i], sep="\t", header=TRUE)
    idf <- idf %>% filter(startsWith(idf$miRNA_region,"mature"), idf$cross.mapped == "N")
    subset_idf <- idf[, c("miRNA_ID", "reads_per_million_miRNA_mapped")]
    subset_idf <- subset_idf %>% group_by(miRNA_ID) %>% summarise(conteo_normalizado = sum(reads_per_million_miRNA_mapped, na.rm = TRUE))
    names(subset_idf)[2] <- isamplename
    if (i==1){
      combined_df <- subset_idf
      rm(subset_idf)
    } else {
      combined_df <- merge(combined_df, subset_idf, by='miRNA_ID', all=TRUE)
      rm(subset_idf)
    }
  }
  rownames(combined_df) <- combined_df$miRNA_ID
  combined_df <- combined_df[,-which(names(combined_df) %in% "miRNA_ID")]
  return(combined_df)
}


preprocesado_mirnas <- function(df, project, porcentaje_max_nulos, 
                                porcentaje_max_na, valor_varianza){
  
  # Filtrar duplicados
  df <- gdcFilterDuplicate(df)
  
  # Filtrar las muestras non-Primary Tumor y non-Solid Tissue Normal en 
  # los metadatos
  df <- gdcFilterSampleType(df)
  
  
  isoform_miRNAs_dir <- paste('metadatos', project, 'miRNAs', sep='/')
  
  
  # Obtenemos el dataframe filtrado y con los conteos normalizados
  MIRCounts <- merge_MIR(df, isoform_miRNAs_dir)
  
  # Ponemos en marcha los filtros
  filtro_cerca_del_cero <- (rowSums(MIRCounts < 1, na.rm = TRUE) 
                            / ncol(MIRCounts)) <= (porcentaje_max_nulos / 100)
  
  valores_na_por_fila <- apply(MIRCounts, 1, function(x) sum(is.na(x)))
  filtro_na <- valores_na_por_fila / ncol(MIRCounts) <= (porcentaje_max_na / 100)
  
  varianza <- apply(MIRCounts, 1, function(x) var(x, na.rm=TRUE))
  filtro_desviacion <- varianza >= valor_varianza
  
  
  filtros <- filtro_cerca_del_cero & filtro_na & filtro_desviacion
  expr_mirnas_tras_filtros <- MIRCounts[filtros,]
  
  # Por último, imputamos valores mediante knn
  expresiones_mirna_traspuesto <- t(expr_mirnas_tras_filtros)
  expr_mirnas_tras_filtros <- t(impute.knn(as.matrix(expresiones_mirna_traspuesto))$data)
  
  return(expr_mirnas_tras_filtros)
}

