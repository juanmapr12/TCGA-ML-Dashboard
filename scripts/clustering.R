
clustering_dendrograma <- function(df_mirnas, df_completo, tipo, dea_y_cox, N, metodo_jerarquico){
  
  if(tipo == "microARN resultantes de los análisis de Expresión y Supervivencia"){
    mirnas_seleccionados <- t(df_completo[,dea_y_cox])
  }else{
    media_mirnas <- apply(df_mirnas, 1, mean)
    sd_mirnas <- apply(df_mirnas, 1, sd)
    
    # Coeficiente de variación = sd / media
    cv_mirna <- sd_mirnas / media_mirnas
    topN_cv <- sort(cv_mirna, decreasing = TRUE)[1:N]
    mirnas_seleccionados <- t(df_completo[,names(topN_cv)])
  }
  
  mirnas_escalados <- t(scale(t(mirnas_seleccionados)))
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
                breaks = seq(0, umbral, length.out = 100),
                cluster_rows = TRUE, 
                cluster_cols = TRUE,
                show_colnames = FALSE,
                annotation_row = data.frame(Cluster = factor(agrupamiento_mirnas)))
  
  p$plot
  
  return(agrupamiento_mirnas)
}


tabla_de_datos_por_cada_cluster <- function(df_mirnas, agrupamiento_mirnas, cluster, 
                                            tipo, limite, clinical_df, df_completo){
  
  mirnas_del_cluster <- names(agrupamiento_mirnas[agrupamiento_mirnas == cluster])
  
  vector_columnas <- c("tipo_muestra",colnames(clinical_df))
  
  submatriz_mirnas_cluster <- df_mirnas[mirnas_del_cluster,]
  
  
  # Promedio de expresión por muestra (columna)
  expr_media_por_muestra <- sort(abs(colMeans(submatriz_mirnas_cluster)), decreasing=TRUE)
  if(tipo == "up-regulated (color rojo)"){
    expresiones_mirna_ordenadas <- sort(expr_media_por_muestra, decreasing = TRUE)
  }else if(tipo == "down-regulated (color azul)"){
    expresiones_mirna_ordenadas <- sort(expr_media_por_muestra, decreasing = FALSE)
  }
  muestras_top_expresion <- names(expresiones_mirna_ordenadas)[1:limite]
  
  resultado_mirnas <- df_completo[muestras_top_expresion, mirnas_del_cluster]
  resultado_clinico <- df_completo[muestras_top_expresion, vector_columnas]
  
  return(list(resultado_mirnas, 
              resultado_clinico, 
              muestras_top_expresion))
}


analisis_columna_clinica <- function(df_completo, 
                                      conj_muestras_sobreexpresadas,
                                      variable_seleccionada){
  
  
  datos_columna_clinica <- df_completo[conj_muestras_sobreexpresadas,]
  
  if(variable_seleccionada %in% lista_columnas_numericas){
    if(variable_seleccionada == "time")
      banda = 10
    else
      banda = 1
    p <- ggplot(datos_columna_clinica, aes(x=get(variable_seleccionada))) + 
      geom_histogram(binwidth = banda, fill="skyblue", color = "black") +
      xlab(variable_seleccionada) +
      theme(
        axis.text.x = element_text(size = 14),   # números/labels en eje X
        axis.text.y = element_text(size = 14),   # números/labels en eje Y
        axis.title.x = element_text(size = 16),  # título eje X
        axis.title.y = element_text(size = 16),  # título eje Y
        plot.title  = element_text(size = 18, hjust = 0.5) # título principal
      )
  }else{
    p <- ggplot(datos_columna_clinica, aes(x=get(variable_seleccionada))) + 
      geom_bar(fill = "skyblue", color = "black") +
      xlab(variable_seleccionada) +
      theme(
        axis.text.x = element_text(size = 14),   # números/labels en eje X
        axis.text.y = element_text(size = 14),   # números/labels en eje Y
        axis.title.x = element_text(size = 16),  # título eje X
        axis.title.y = element_text(size = 16),  # título eje Y
        plot.title  = element_text(size = 18, hjust = 0.5) # título principal
      )
  }
  
  print(p)
}

