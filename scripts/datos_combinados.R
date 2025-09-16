# Función para crear df que contenga tanto datos clínicos como 
# de expresiones de mirnas
procesar_datos_mirna_y_clinicos <- function(x, y) {
  
  # Extraer nombres de muestra
  nombres_muestras <- colnames(x)

  # Clasificar según los últimos dos dígitos
  tipo_muestra <- ifelse(
    grepl("-01$", nombres_muestras), "enferma",
    ifelse(grepl("-11$", nombres_muestras), "sana", NA)
  )

  # Crear nueva fila con clasificación
  x <- rbind("tipo_muestra" = tipo_muestra, x)
  
  columnas_clinicas <- colnames(y)
  for(i in 1:length(columnas_clinicas)){
    nueva_obs <- rep(NA, ncol(x))   # Columna de NA's que se repite el número de columnas de x
    names(nueva_obs) <- colnames(x)
    x <- rbind(x, nueva_obs)        # La añadimos al df            
    rownames(x)[nrow(x)] <- columnas_clinicas[i]   
    # Le asociamos el nombre de la columna clínica correspondiente
  }
  
  id_pacientes <- rownames(y)
  for(i in 1:length(id_pacientes)){
    # Comprobamos si hay match entre un id del paciente de los datos clínicos
    # y las muestras provenientes de los datos de mirnas
    if(any(sapply(id_pacientes[i], function(p) grepl(p, colnames(x))))){
      # En caso afirmativo, vemos en qué indices coinciden
      vector_coincidencias <- which(grepl(id_pacientes[i],colnames(x)))
      for(j in 1:length(vector_coincidencias)){
        pos <- vector_coincidencias[j]
        x[,pos][(nrow(x)-(length(columnas_clinicas)-1)):nrow(x)] <- t(y[id_pacientes[i],])
      }
    }
  }
  
  return(as.data.frame(t(x)))
}


transformacion_variables <- function(df, clinical_df){
  
  df$tipo_muestra <- as.factor(df$tipo_muestra)
  col_clinicas <- colnames(clinical_df)
  for (i in 2:(ncol(df) - length(col_clinicas))){
    df[,i] <- as.numeric(df[,i])
  }
  for (i in (ncol(df) - length(col_clinicas)):ncol(df)){
    columna_iesima <- colnames(df)[i]
    if (columna_iesima %in% lista_columnas_numericas){
      df[,i] <- as.numeric(df[,i])
    }else if(columna_iesima %in% lista_columnas_categoricas){
      df[,i] <- as.factor(df[,i])
    }
  }
  
  return(df)
}


limpieza_final <- function(df1, df2){
  muestras_a_eliminar <- rownames(df1)[is.na(df1$time) | is.na(df1$vital_status)]
  
  df1_filtrado <- df1[!(rownames(df1) %in% muestras_a_eliminar), ]
  df2_filtrado <- df2[,!(colnames(df2) %in% muestras_a_eliminar)]
  
  return(list(df1_filtrado, df2_filtrado))
}
