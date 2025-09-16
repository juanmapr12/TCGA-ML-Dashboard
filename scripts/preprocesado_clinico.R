
lista_columnas_numericas <- c("age",
                               "time")

lista_columnas_categoricas_nominales <- c("ethnicity",
                                           "gender",
                                           "race",
                                           "vital_status")

lista_columnas_categoricas_ordinales <- 
  c("clinical_stage",
    "clinical_T",
    "clinical_N",
    "clinical_M",
    "pathologic_stage",
    "pathologic_T",
    "pathologic_M",
    "pathologic_N",
    "tumor_event_after_treatment_no",
    "tumor_event_after_treatment_yes",
    "additional_pharmaceutical_therapy_yes",
    "additional_pharmaceutical_therapy_no",
    "additional_radiation_therapy_yes",
    "additional_radiation_therapy_no")


lista_columnas_categoricas <- c(lista_columnas_categoricas_nominales,
                                 lista_columnas_categoricas_ordinales)


preprocesado_clinico <- function(df){
  
  rownames_orig <- rownames(df)
  
  names(df)[names(df) == "age_at_initial_pathologic_diagnosis"] <- "age"
  names(df)[names(df) == "new_tumor_event_after_initial_treatment_no"] <- "tumor_event_after_treatment_no"
  names(df)[names(df) == "new_tumor_event_after_initial_treatment_yes"] <- "tumor_event_after_treatment_yes"
  
  
  # Columna 'time'
  
  # Se creará la columna 'time', que representará el tiempo de un paciente 
  # desde el diagnóstico inicial del cáncer hasta que se produce el evento 
  # de interés (la muerte del paciente) o hasta el último chequeo. 
  # Para ello utilizamos las variables del dataset 'days_to_death' y 
  # 'days_to_last_follow_up', eliminándolas posteriormente debido a que 
  # contienen información oculta sobre nuestra variable objetivo 
  # 'vital_status'.
  df$time <- ifelse(!is.na(df$days_to_death),
                             df$days_to_death,
                             df$days_to_last_followup)
  df$days_to_death <- NULL
  df$days_to_last_followup <- NULL
  
  
  # Observamos que:
  # - En muchas de las variables, los valores 'NA' están representados como cadenas
  # de texto en lugar de como valores faltantes.
  # - Hay varias variables con todos sus valores nulos.
  # - La variable 'pathologic_categories' parece un resumen de la terna ('pathologic_T',
  # 'pathologic_N','pathologic_M'). 
  
  
  # 1.
  df[df == "NA"] <- NA
  
  # 2.
  porcentaje_nulos <- 60
  umbral_valores_nulos <- trunc(nrow(df) * (porcentaje_nulos / 100))
  datos_nulos <- colSums(is.na(df))
  for (i in 1:length(datos_nulos)){
    if (datos_nulos[i] >= umbral_valores_nulos)
      df <- df[, -which(names(df) == names(datos_nulos[i]))]
  }
  
  # 3.
  df$pathologic_categories <- NULL
  
  
  # 4.
  df <- df %>% filter(!is.na(time)) %>% filter(!is.na(vital_status))
  
  
  lista_variables <- colnames(df)
  
  
  # 5. 
  valores_faltantes <- function(df){
    for (i in 1:length(colnames(df))){
      col <- colnames(df)[[i]]
      a <- df[,col][is.na(df[,col])]
      if (length(a) != 0){
        if (col %in% lista_columnas_numericas){
          df[,col][is.na(df[,col])] <- 0
        }else{
          df[,col][is.na(df[,col])] <- "is_missing"
        }
      }
    }
    return(df)
  }
  df <- valores_faltantes(df)
  
  
  
  ############### ANÁLISIS UNIVARIADO Y BIVARIADO DE CADA VARIABLE  ##################
  
  
  # Agilizamos las referencias a las columnas del dataframe:
  attach(df)
  
  # Primero las numéricas:
  
  # Función que las convierte de caracter a numéricas:
  convierte_aNumerica <- function(nombre_var){
    df[[nombre_var]] <- as.numeric(df[[nombre_var]])
    return(df)
  }
  
  # Función para comprobar normalidad 
  sigueDistribucionNormal <- function(nombre_var, variable){
    bool <- FALSE
    muestra_alive <- variable[vital_status == "Alive"]
    muestra_dead <- variable[vital_status == "Dead"]
    
    pvalor_alive <- shapiro.test(muestra_alive)$p.value
    pvalor_dead <- shapiro.test(muestra_dead)$p.value
    if (pvalor_alive > 0.05 & pvalor_dead > 0.05)
      bool <- TRUE
    
    return(list(bool = bool, muestra_alive = muestra_alive,
                muestra_dead = muestra_dead))
  }
  
  lista_var_signif <- list("vital_status", "time")
  
  testNumericas_VariablesAsociadas <- function(nombre_var){
    lista <- list()
    resultado <- sigueDistribucionNormal(nombre_var, df[[nombre_var]])
    if (resultado$bool){
      p_valor_medias <- 0  # Inicializo
      p_valor_varianza <- var.test(resultado$muestra_alive, resultado$muestra_dead)
      if(p_valor_varianza < 0.05){
        # CASO 1: Las varianzas son iguales
        p_valor_medias <- t.test(resultado$muestra_alive, resultado$muestra_dead, 
                                 alternative = "two.sided",
                                 mu          = 0,
                                 var.equal   = TRUE,
                                 conf.level  = 0.95)$p.value
      }else{
        # CASO 2: Las varianzas no son iguales
        p_valor_medias <- t.test(resultado$muestra_alive, resultado$muestra_dead, 
                                 alternative = "two.sided",
                                 mu          = 0,
                                 var.equal   = FALSE,
                                 conf.level  = 0.95)$p.value
      }
      if (p_valor_medias < 0.05){
        lista <- list(nombre_var)
      }
    }else{
      # Prueba U de Mann-Whitney
      p_valor <- wilcox.test(resultado$muestra_alive, resultado$muestra_dead, paired = FALSE)$p.value
      if(p_valor < 0.05){
        lista <- list(nombre_var)
      }
    }
    return(lista)
  }
  
  # Y ahora las categóricas:
  
  # Función que las convierte de caracter a categóricas:
  convierte_aCategorica <- function(nombre_var){
    df[[nombre_var]] <- as.factor(df[[nombre_var]])
    return(df)
  }
  
  testCategoricas_VariablesAsociadas <- function(nombre_var){
    lista <- list()
    variable <- df[[nombre_var]]
    df <- data.frame(vital_status, variable)
    tabla <- table(df$vital_status, df$variable)
    # Aplicamos test Chi-cuadrado considerando las condiciones de Cochran
    chi = chisq.test(tabla)
    condicion1 <- mean(chi$expected < 5) <= 0.20
    condicion2 <- all(chi$expected >= 1)
    if(condicion1 & condicion2){
      if(chi$p.value < 0.05){
        lista <- list(nombre_var)
      }
    }else{
      v <- assocstats(tabla)$cramer
      if (v >= 0.5) {
        lista <- list(nombre_var)
      }
    }
    return(lista)
  }
  
  
  ####### Bucle variables numéricas ########
  
  for(i in 1:length(lista_columnas_numericas)){
    if(lista_columnas_numericas[i] %in% lista_variables){
      nombre <- lista_columnas_numericas[i]
      df <- convierte_aNumerica(nombre)
      if(!lista_columnas_numericas[i] %in% c("vital_status","time"))
        lista_var_signif <- append(lista_var_signif,
                                   testNumericas_VariablesAsociadas(nombre))
    }
  }
  
  
  ####### Bucle variables categóricas #######
  
  for(i in 1:length(lista_columnas_categoricas)){
    if(lista_columnas_categoricas[i] %in% lista_variables){
      nombre <- lista_columnas_categoricas[i]
      df <- convierte_aCategorica(nombre)
      if(!lista_columnas_categoricas[i] %in% c("vital_status","time"))
        lista_var_signif <- append(lista_var_signif,
                                   testCategoricas_VariablesAsociadas(nombre))
    }
  }
  
  
  comprueba_relaciones <- function(lista_pares_variables, lista){
    for(i in 1:(length(lista_pares_variables)/2)){
      var1 <- lista_pares_variables[,i][[1]]
      var2 <- lista_pares_variables[,i][[2]]
      res <- assocstats(table(get(var1), get(var2)))$cramer
      if(!is.nan(res) & res >= 0.5){
        lista <- append(lista, var1)
      }
    }
    return(lista)
  }
  
  significativas_sin_tiempo_ni_vital_status <- lista_var_signif[-c(1,2)]
  pares <- combn(significativas_sin_tiempo_ni_vital_status, 2)
  lista_negra <- list()
  lista_negra <- unique(comprueba_relaciones(pares, lista_negra))
  
  lista_negra <- lista_negra[lista_negra != "vital_status"]
  
  # Nos quedamos con la diferencia de las significativas menos las que ya
  # explican a otras (la "lista negra")
  lista_var_definitiva <- setdiff(lista_var_signif, lista_negra)
  
  df <- df[,unlist(lista_var_definitiva)]
  
  return(list(df, unlist(lista_var_definitiva)))
}


# Función para el análisis univariado
analisisUnivariado_Numericas <- function(df, nombre_var){
  
  variable <- df[[nombre_var]]
  
  # Calculamos estadísticas
  media     <- mean(variable, na.rm = TRUE)
  mediana   <- median(variable, na.rm = TRUE)
  sd_val    <- sd(variable, na.rm = TRUE)
  rango     <- range(variable, na.rm = TRUE)
  cuartiles <- quantile(variable, probs = c(0.25, 0.5, 0.75, 0.9), na.rm = TRUE)
  iqr_val   <- IQR(variable, na.rm = TRUE)
  asim      <- skewness(variable, na.rm = TRUE)
  apunt     <- kurtosis(variable, na.rm = TRUE)
  
  # Imprimimos resumen
  cat(sprintf("Resumen descriptivo de la variable %s:\n", nombre_var))
  cat("\n")
  cat(sprintf("Media:                %.2f\n", media))
  cat(sprintf("Mediana:              %.2f\n", mediana))
  cat(sprintf("Desviación estándar:  %.2f\n", sd_val))
  cat(sprintf("Rango:                %.2f – %.2f\n", rango[1], rango[2]))
  cat(sprintf("IQR (Q3 - Q1):        %.2f\n", iqr_val))
  cat(sprintf("Asimetría (skewness): %.2f\n", asim))
  cat(sprintf("Apuntamiento (kurtosis): %.2f\n", apunt))
  cat("\n Cuartiles:\n")
  cat(sprintf("  Q1 (25%%):  %.2f\n", cuartiles[1]))
  cat(sprintf("  Q2 (50%%):  %.2f\n", cuartiles[2]))
  cat(sprintf("  Q3 (75%%):  %.2f\n", cuartiles[3]))
  cat(sprintf("  90%%:       %.2f\n", cuartiles[4]))
  
}

resultadoBivariado_Numericas <- function(df, col){
  return(ggplot(df, aes(x=get(col),
                        fill=vital_status,
                        colour=vital_status)) + 
           geom_histogram(position="identity", alpha=0.5) +
           xlab(col))
}


# Función para el análisis univariado
analisisUnivariado_Categoricas <- function(df, nombre_var){
  cat(sprintf("Resumen de la variable %s, donde indicamos\nel porcentaje sobre 1 en el que aparece cada valor: \n", nombre_var))
  print(prop.table(table(df[[nombre_var]])))
}

resultadoBivariado_Categoricas <- function(df, col){
  return(ggplot(df, aes(x=get(col),
                 fill=vital_status,
                 color=vital_status)) + 
    geom_bar(position=position_dodge()) +
    xlab(col))
}


grafica_relacion_vitalstatus <- function(df, col){
  if(col %in% lista_columnas_categoricas){
    resultadoBivariado_Categoricas(df, col)
  }else if(col %in% lista_columnas_numericas){
    resultadoBivariado_Numericas(df, col)
  }
}

datos_analisis_univariado <- function(df, col){
  if(col %in% lista_columnas_categoricas){
    analisisUnivariado_Categoricas(df, col)
  }else if(col %in% lista_columnas_numericas){
    analisisUnivariado_Numericas(df, col)
  }
}




