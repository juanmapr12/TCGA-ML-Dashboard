
clustering_dendrograma <- function(df_mirnas, df_completo, N, metodo_jerarquico){
  
  media_mirnas <- apply(df_mirnas, 1, mean)
  sd_mirnas <- apply(df_mirnas, 1, sd)
  
  # Coeficiente de variación = sd / media
  cv_mirna <- sd_mirnas / media_mirnas
  topN_cv <- sort(cv_mirna, decreasing = TRUE)[1:N]
  mirnas_topN <- t(df_completo[,names(topN_cv)])
  
  mirnas_escalados <- t(scale(t(mirnas_topN)))
  D <- as.dist(1 - cor(t(mirnas_escalados), method = "spearman"))
  conglo <- hclust(D, method=metodo_jerarquico)
  p <- plot(conglo, cex=0.8, labels=rownames(mirnas_escalados),
       main="Dendrograma", xlab=conglo$method,
       ylab="Distancias", col="blue")
  
  p$plot
  
  return(list(mirnas_escalados, conglo))
}


clustering_heatmap <- function(mirnas_escalados, conglomerado, k, umbral){
  agrupamiento_mirnas <- cutree(conglomerado, k)
  p <- pheatmap(mirnas_escalados,
                breaks = seq(-umbral, umbral, length.out = 100),
                cluster_rows = TRUE, 
                cluster_cols = TRUE,
                show_colnames = FALSE,
                annotation_row = data.frame(Cluster = factor(agrupamiento_mirnas)))
  
  p$plot
  
  return(agrupamiento_mirnas)
}


tabla_de_datos_por_cada_cluster <- function(df_mirnas, agrupamiento_mirnas, cluster, 
                                            limite, clinical_df, df_completo){
  
  mirnas_del_cluster <- names(agrupamiento_mirnas[agrupamiento_mirnas == cluster])
  
  vector_columnas <- c("tipo_muestra",colnames(clinical_df))
  
  submatriz_mirnas_cluster <- df_mirnas[mirnas_del_cluster, ]
  
  
  # Promedio de expresión por muestra (columna)
  expr_media_por_muestra <- sort(abs(colMeans(submatriz_mirnas_cluster)), decreasing=TRUE)
  expr_ordenado <- sort(expr_media_por_muestra, decreasing = TRUE)
  muestras_top_expr_alta <- names(expr_ordenado)[1:limite]
  
  resultado_mirnas <- df_completo[muestras_top_expr_alta, mirnas_del_cluster]
  resultado_mirnas_dos_decimales <- as.data.frame(
    lapply(resultado_mirnas, function(x) round(as.numeric(x), 2))
  )
  resultado_clinico <- df_completo[muestras_top_expr_alta, vector_columnas]
  
  return(list(resultado_mirnas_dos_decimales, 
              resultado_clinico, 
              muestras_top_expr_alta))
}


analisis_columna_clinica <- function(df_completo, 
                                      conj_muestras_sobreexpresadas,
                                      variable_seleccionada){
  
  
  datos_columna_clinica <- df_completo[conj_muestras_sobreexpresadas,]
  
  if(variable_seleccionada %in% lista_columnas_numericas){
    if(variable_seleccionada == "time")
      banda = 300
    else
      banda = 3
    p <- ggplot(datos_columna_clinica, aes(x=get(variable_seleccionada))) + 
      geom_histogram(binwidth = banda, fill="skyblue", color = "black") +
      xlab(variable_seleccionada)
  }else{
    p <- ggplot(datos_columna_clinica, aes(x=get(variable_seleccionada))) + 
      geom_bar(fill = "skyblue", color = "black") +
      xlab(variable_seleccionada)
  }
  
  p$plot
}

