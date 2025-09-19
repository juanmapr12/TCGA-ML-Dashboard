

df_escogido <- function(datos_a_usar, df_clinico, df_completo, mirnas_dea, mirnas_cox){
  if(datos_a_usar == "Variables clínicas"){
    return(df_clinico)
    
  }else if(datos_a_usar == "Variables clínicas + microARN tras Expresión Diferencial"){
    return(df_completo[,union(colnames(df_clinico), mirnas_dea)])
    
  }else if(datos_a_usar == "Variables clínicas + microARN tras Supervivencia"){
    return(df_completo[,union(colnames(df_clinico), mirnas_cox)])
    
  }else if(datos_a_usar == "Variables clínicas + microARN tras Expresión Diferencial y Supervivencia"){
    return(df_completo[,union(colnames(df_clinico), union(mirnas_dea, mirnas_cox))])
    
  }
}

pipeline_tidyverse_1 <- function(df, p, folds, modelo_a_ejecutar){
  
  set.seed(123)
  # Separamos en entrenamiento y prueba de manera estratificada:
  data_split <- initial_split(df, 
                             prop = p, 
                             strata = vital_status)
  data_train <- training(data_split)
  data_test  <- testing(data_split)
  
  
  # Creamos conjuntos para la validación cruzada:
  data_folds <- vfold_cv(data_train, 
                              v = folds, 
                              strata = vital_status)
  
  # Receta del preprocesado:
  data_recipe <- 
    recipe(vital_status ~ ., 
           data = data_train) %>%
    step_range(all_numeric_predictors(), ranges = c(0,1)) %>%
    step_dummy(all_nominal(), -all_outcomes()) %>%
    step_lincomb(all_predictors())
  
  # Especificación del modelo:
  modelo <- switch(modelo_a_ejecutar, 
             "Random Forest" = rand_forest(
               trees = tune(),
               min_n = tune(),
             ),
             "Gradient Boosting" = boost_tree(
               trees = tune(),
               tree_depth = tune(),
               learn_rate = tune(),
               loss_reduction = tune(),
               sample_size = tune(),
               min_n = tune(),
             ),
             "Regresión logística" = logistic_reg(
               penalty = tune(), 
               mixture = tune()
             ),
             "SVM (Support Vector Machines)" = svm_rbf(
               cost = tune(), 
               rbf_sigma = tune()
              )
          )
         
  motor <- switch(modelo_a_ejecutar, 
                  "Random Forest" = "ranger",
                  "Gradient Boosting" = "xgboost",
                  "Regresión logística" = "glmnet",
                  "SVM (Support Vector Machines)" = "kernlab"
          )
  
  modelo_escogido <- modelo %>%
    set_engine(motor) %>%
    set_mode("classification")
  
  # Establecimiento del flujo de trabajo:
  model_tune_wf <-
    workflow() %>%
    add_recipe(data_recipe) %>%
    add_model(modelo_escogido)
  
  
  # Definimos métricas:
  metricas <- metric_set(accuracy, roc_auc, brier_class)
  
  
  # Resultados del modelo con los mejores hiperparámetros:
  model_results_after_tuning <- 
    model_tune_wf %>% 
    tune_grid(resamples = data_folds,
              metrics = metricas)
  
  
  return(list(model_results_after_tuning, model_tune_wf, data_split))
}

pipeline_tidyverse_2 <- function(model_results_after_tuning, model_tune_wf, data_split){
  
  set.seed(123)
  # Especificación de la mejor métrica a seguir:
  # metrica_escogida <- switch(metrica_a_seguir, 
  #                            "Accuracy" = "accuracy",
  #                            "Área bajo la curva ROC" = "roc_auc",
  #                            "Puntuación Brier" = "brier_class")
  
  # Seleccionamos la mejor según la métrica
  data_best <-
    model_results_after_tuning %>% 
    select_best(metric = "accuracy")
  
  # Creamos el workflow final
  data_wf_definitivo <- 
    model_tune_wf %>%
    finalize_workflow(data_best)
  
  # Entrenamos el workflow final en el conjunto de entrenamiento y evaluamos
  # el modelo en el conjunto de prueba
  data_test_results <-
    data_wf_definitivo %>% 
    last_fit(split = data_split)
  
  # Por último, obtenemos las predicciones
  predicciones <- data_test_results %>% 
    collect_predictions()
  
  metrica_accuracy <- accuracy(predicciones, 
                               truth = vital_status, 
                               estimate = .pred_class)
  
  metrica_precision <- precision(predicciones, 
                                 truth = vital_status, 
                                 estimate = .pred_class)
  
  metrica_recall <- recall(predicciones, 
                           truth = vital_status, 
                           estimate = .pred_class)
  
  metrica_f1_score <- f_meas(predicciones, 
                             truth = vital_status, 
                             estimate = .pred_class)
  
  metrica_roc_auc <- 
    predicciones %>%
    roc_auc(truth = vital_status, .pred_Alive)
  
  valores_curva_roc <- 
    predicciones %>% 
    roc_curve(truth = vital_status, .pred_Alive) %>%
    autoplot() +
    geom_abline(linetype = "dashed", color = "orange") +
    labs(
      title = "Curva ROC",
      subtitle = "Clasificación de vital_status",
      x = "Tasa de falsos positivos (1 - Especificidad)",
      y = "Tasa de verdaderos positivos (Sensibilidad)",
      color = "Clase"
    ) +
    theme_minimal(base_size = 17) +
    theme(
      plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
      plot.subtitle = element_text(size = 13, hjust = 0.5),
      axis.title = element_text(size = 13, face = "bold"),
      axis.text = element_text(size = 11),
      legend.position = "bottom",
      legend.title = element_text(face = "bold")
    )
  
  valores_matriz_conf <- 
    predicciones %>% 
    conf_mat(truth = vital_status, estimate = .pred_class) %>% 
    autoplot(type = "heatmap") +
    scale_fill_gradient(low = "white", high = "#2c7bb6") +
    labs(
      title = "Matriz de confusión",
      subtitle = "Clasificación de vital_status",
      x = "Valor real",
      y = "Predicción"
    ) +
    theme_minimal(base_size = 17) +
    theme(
      plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
      plot.subtitle = element_text(size = 13, hjust = 0.5),
      axis.text = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 13, face = "bold"),
      legend.position = "right",
      legend.title = element_text(face = "bold")
    )
  
  
  return(list(metrica_accuracy %>% pull(.estimate),
              metrica_precision %>% pull(.estimate),
              metrica_recall %>% pull(.estimate),
              metrica_f1_score %>% pull(.estimate),
              metrica_roc_auc %>% pull(.estimate),
              valores_curva_roc,
              valores_matriz_conf)
         )
}