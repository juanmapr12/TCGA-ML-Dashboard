##### Genes diana #####


resultados_genes_diana <- function(mirnas_seleccionados, tipo_datos){
  
  valor <- switch(
    tipo_datos,
    "Interacciones validadas experimentalmente" = "validated",
    "Predicciones computacionales" = "predicted",
    "Ambas" = "all"
  )

  multimir_results <- get_multimir(org     = 'hsa',
                                   mirna   = mirnas_seleccionados,
                                   table   = valor,
                                   summary = TRUE)
  
  return(multimir_results@data)
}


##### Análisis de enriquecimiento #####


analisis_enriquecimiento <- function(genes_diana, tipo_analisis, p_valor, q_valor, ontologia){
  
  entrez_ids <- genes_diana$target_entrez
  
  # Biological ID translation, convierte identificadores de genes entre distintos sistemas.
  # En nuestro caso, nos interesa pasar de un entrez_id a su símbolo correspondiente
  bitr(entrez_ids, 
       fromType = "ENTREZID", 
       toType = "SYMBOL", 
       OrgDb = org.Hs.eg.db)
  
  valor <- switch(
    ontologia,
    "Proceso biológico" = "BP",
    "Función molecular" = "MF",
    "Componente celular" = "CC"
  )
  
  if(tipo_analisis == "GO"){
    enrich_object <- enrichGO(gene          = entrez_ids,
                              OrgDb         = org.Hs.eg.db,
                              keyType       = "ENTREZID",
                              ont           = valor,
                              pAdjustMethod = "BH",
                              pvalueCutoff  = p_valor,
                              qvalueCutoff  = q_valor)
  }else if(tipo_analisis == "KEGG"){
    enrich_object <- enrichKEGG(gene         = entrez_ids,
                                organism     = "hsa",
                                pvalueCutoff = p_valor,
                                qvalueCutoff = q_valor)
  }
  return(enrich_object)
}


diagrama_barras <- function(enrich_object, max_categorias, titulo){
  barplot(enrich_object,
          showCategory = max_categorias,
          title = titulo)
}
