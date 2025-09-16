# Análisis de expresión diferencial:

saca_pareadas <- function(muestras_totales){
  
  # 1.
  id_paciente <- sapply(strsplit(muestras_totales, "-"), 
                        function(x) paste(x[1:3], collapse = "-"))
  
  # 2. 
  tipo_muestra <- sapply(strsplit(muestras_totales, "-"), 
                         function(x) x[4])
  
  # 3.
  df <- data.frame(
    paciente = id_paciente,
    tipo = tipo_muestra
  )
  
  # 4. 
  muestras_pareadas <- df %>%
    group_by(paciente) %>%
    filter(all(c("01", "11") %in% tipo)) %>%
    ungroup() %>%
    arrange(paciente) %>%
    mutate(union = paste(paciente, tipo, sep="-")) %>%
    pull(union)
  
  return(muestras_pareadas)
}


analisis_expr_dif <- function(muestras_escogidas, df){
  
  grupo <- vector(mode="character")
  
  if (muestras_escogidas == "Totales"){
    muestras <- colnames(df)
  }else if (muestras_escogidas == "Pareadas"){
    muestras <- saca_pareadas(colnames(df))
  }
  
  for (i in 1:length(muestras)){
    tipo_muestra <- sapply(strsplit(muestras[i], "-"),
                           function(x) x[4])
    if (tipo_muestra == "01"){
      grupo <- c(grupo, "enfermo")
    }else{
      grupo <- c(grupo, "sano")
    }
  }
  
  grupo <- factor(grupo)
  grupo <- relevel(grupo, ref = "enfermo")
  design <- model.matrix(~ grupo)
  
  # Simplemente expresamos los CPM en log2 para poder aplicarles limma
  log_expr <- log2(df[,muestras] + 1)
  
  fit <- lmFit(log_expr, design)
  fit <- eBayes(fit)
  tt <- topTable(fit, coef = 2, number = Inf)
  
  return(tt)
}


volcano <- function(tabla_resultados, muestras_escogidas, umbral_logFC, umbral_pvalor_ajustado){
  
  if (muestras_escogidas == "Totales"){
    titulo <- "Volcano Plot para el total de muestras"
  }else if (muestras_escogidas == "Pareadas"){
    titulo <- "Volcano Plot para el conjunto de muestras pareadas"
  }
  print(EnhancedVolcano(tabla_resultados,
                        title = titulo,
                        lab = rownames(tabla_resultados),
                        x = 'logFC',
                        y = 'adj.P.Val',
                        FCcutoff = umbral_logFC,
                        pCutoff = umbral_pvalor_ajustado))
}


mirnas_que_cumplen_ambos_filtros <- function(resultado_analisis, umbral_logFC, umbral_pvalor_ajustado){
  tt_filtrado <- subset(
    resultado_analisis,
    adj.P.Val < umbral_pvalor_ajustado & abs(logFC) > umbral_logFC
  )
  return(tt_filtrado[,c("logFC", "adj.P.Val")])
}

