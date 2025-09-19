options(warn=-1)

library(shiny)
library(shinydashboard)
library(DT)
library(shinyWidgets)
library(shinyalert)

# Cargamos scripts necesarios
source("scripts/preparacion_entorno.R")
source("scripts/preprocesado_clinico.R")
source("scripts/preprocesado_mirnas.R")
source("scripts/analisis_expr_dif.R")
source("scripts/datos_combinados.R")
source("scripts/analisis_supervivencia.R")
source("scripts/genes_diana.R")
source("scripts/clustering.R")
source("scripts/clasificacion_ML.R")


ui <- dashboardPage(
  
  dashboardHeader(
    title = "TCGA Dashboard",
    tags$li(class = "dropdown",
            tags$b(textOutput("proyecto_utilizado"), style = "color:white; font-size:17px;"))
  ),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Introducción", tabName = "intro", icon = icon("book-open")),
      menuItem("Preprocesado",
               menuSubItem("microARN", tabName = "pre_mirna", icon = icon("dna")),
               menuSubItem("Clínico", tabName = "pre_clinico", icon = icon("stethoscope"))
      ),
      menuItem("Expresión diferencial", tabName = "expr_dif", icon = icon("volcano")),
      menuItem("Análisis supervivencia", tabName = "surv", icon = icon("line-chart")),
      menuItem("Enriquecimiento", tabName = "enrich", icon = icon("chart-bar")),
      menuItem("Clustering", tabName = "cluster", icon = icon("circle-nodes")),
      menuItem("Clasificación ML", tabName = "clasif_ml", icon = icon("robot"))
    )
  ),
  
  dashboardBody(
    tags$head(
      tags$style(HTML("
      .slider-custom .js-irs-1 .irs-bar{
        background: #f5f5f5;
        height: 9.5px;
      }
      
      .slider-custom .js-irs-1 .irs-line{
        background: #3c8dbc;
        height: 8px;
      }
      
      .inputs-espaciados .shiny-input-container, 
      .inputs-espaciados .form-group,
      .inputs-espaciados .help-block{
        margin-bottom: 22px !important;
      }
      
      body, label, input, select, textarea, button, .dropdown span, .sidebar-menu .treeview-menu > li > a {
        font-size: 18px !important;
      }
      
      .center-col { 
        text-align: center; 
        margin-bottom: 30px;
      }                       
      
      .center-col .shiny-input-container {
        display: inline-block;                                
        text-align: left;                                     
      }
      
      .center-col .control-label {
        width: 100%;                                          
        text-align: center;                                 
        margin-bottom: 6px;
      }
      "))
    ),
    
    tabItems(
      
      #### Introducción: Preparación del entorno (ui) ####
      tabItem(tabName = "intro",
              fluidRow(
                box(
                  width = 8,
                  title = "Datos del TCGA",
                  status = "primary",
                  solidHeader = TRUE,
                  selectInput("project", "Seleccione un proyecto TCGA:",
                              choices = c(
                                "TCGA-ACC", "TCGA-BLCA", "TCGA-CESC", "TCGA-CHOL", "TCGA-COAD",
                                "TCGA-DLBC", "TCGA-GBM", "TCGA-HNSC", "TCGA-KICH", "TCGA-KIRC",
                                "TCGA-KIRP", "TCGA-LAML", "TCGA-LGG", "TCGA-LIHC", "TCGA-LUAD",
                                "TCGA-LUSC", "TCGA-MESO", "TCGA-OV", "TCGA-PAAD", "TCGA-PCPG",
                                "TCGA-PRAD", "TCGA-READ", "TCGA-SARC", "TCGA-SKCM", "TCGA-STAD",
                                "TCGA-TGCT", "TCGA-THCA", "TCGA-THYM", "TCGA-UCS", "TCGA-UCEC",
                                "TCGA-UVM"
                              ), select = "TCGA-HNSC", selectize = TRUE
                  ),
                  tags$div(style = "text-align: center;",
                    helpText("Recuerde que puede buscar el proyecto que desee borrando y escribiendo su nombre dentro de la propia barra de búsqueda."),
                  ),
                  selectInput("data_type", "Seleccione un tipo de datos:",
                              choices = c("miRNAs", "Clinical")),
                  tags$div(style = "text-align: center;",
                           actionButton("download", 
                                        "Descargar/Obtener datos", 
                                        style = "background-color: #4CAF50; color: white; font-weight: bold; border: none;"))
                ),
                
                box(
                  width = 4,
                  title = "Aviso",
                  status = "warning",
                  solidHeader = TRUE,
                  helpText(tags$strong("Los archivos del proyecto se almacenan en local (carpeta metadatos). Por tanto, mientras más archivos contenga el proyecto más grande será la espera, por tanto se ruega que sea paciente. Si éste ya viene descargado (como el del caso de prueba que se da, el TCGA-HNSC) no será necesario descargar nada adicional."))
                ),
                
                box(
                  width = 12,
                  title = "Estado de descarga",
                  status = "info",
                  solidHeader = TRUE,
                  textOutput("status")
                ), 
                
                box(
                  width = 12,
                  title = "Proyecto que se va a usar",
                  status = "success",
                  solidHeader = TRUE,
                  uiOutput("proyecto_a_usar")
                )
              )
              
      ),
      
      #### Preprocesado mirnas (ui) ####
      tabItem(tabName = "pre_mirna",
              
              fluidRow(
                box(
                  width = 12, title = "Parámetros a modificar en el preprocesado", status = "primary", solidHeader = TRUE,
                  tags$div(class = "inputs-espaciados",
                    helpText("Una vez normalizados los microARN a lecturas mapeadas por millón (lo que quiere decir que los conteos de expresiones de todos los microARN para una sola muestra no suman uno, sino un millón), se le pide al usuario que:"),
                    numericInput(
                      inputId = "porcentaje_max_nulos",
                      label = "Establezca un porcentaje máximo de pacientes permitidos para las que las distintas expresiones del microARN superen el valor 1 (un valor cercano al 0)",
                      value = 50, min = 0, max = 100, step = 1
                    ),
                    helpText("Es decir, si se escoge el 50% (valor establecido por defecto) quiere decir que se rechazarán los microARN cuya expresión en pacientes es inferior a 1 en, al menos, el 50% del total."),
                    numericInput(
                      inputId = "porcentaje_max_na",
                      label = "Establezca un porcentaje máximo de pacientes a los que se les permite que las expresiones del microARN sean un valor nulo",
                      value = 10, min = 0, max = 20, step = 1
                    ),
                    helpText("Es decir, si se escoge el 10% (valor establecido por defecto) quiere decir que se rechazarán los microARN cuyas expresiones o conteos en la totalidad de pacientes superan el 10% de valores nulos. Se recomienda no superar el 20% para evitar introducir sesgos en la posterior imputación de los valores NA."),
                    numericInput(
                      inputId = "valor_varianza",
                      label = "Establezca un valor mínimo de la varianza a superar por cada microARN",
                      value = 0.5, min = 0, max = 500, step = 0.1
                    )
                  ),
                  helpText("Por último, se imputan los valores de expresión nulos resultantes mediante la función impute.knn del proyecto Bioconductor"),
                  tags$div(style = "text-align: center;",
                           actionButton(
                             inputId = "preprocesado_mirnas",
                             label = "Ejecutar preprocesado con estos valores",
                             style = "background-color: #4CAF50; color: white; font-weight: bold; border: none;"
                           )
                  )
                ),
                
                box(
                  width = 12,
                  title = "Estado del preprocesado",
                  status = "info",
                  solidHeader = TRUE,
                  uiOutput("estado_pre_mirnas")
                )
              )
      ),
      
      
      #### Preprocesado clínico (ui) ####
      
      tabItem(tabName = "pre_clinico",
      
              fluidRow(
                box(
                  width = 5,
                  title = "Selección de características o columnas",
                  status = "primary",
                  solidHeader = TRUE,
                  uiOutput("pre_clinico_seleccion"),
                  tags$div(style = "text-align: center;",
                           actionButton(
                             inputId = "preprocesado_clinico",
                             label = "Ejecutar selección",
                             style = "background-color: #4CAF50; color: white; font-weight: bold; border: none;"
                           )
                  )
                ),
                box(
                  width = 7, title = "Columnas filtradas en la selección", status = "info", solidHeader = TRUE,
                  uiOutput("columnas_signif")
                ),
                tabBox(
                  width = 7, side = "left",
                  tabPanel(
                    "Gráfica resultante de la relación",
                    plotOutput("grafica_pre_clinico")
                  ),
                  tabPanel(
                    "Análisis univariado",
                    tags$style("#univariado {font-size: 17px; color: darkblue;}"),
                    verbatimTextOutput("univariado")
                  )
                )
              )
      ),
      
      
      
      #### Análisis de expresión diferencial (ui) ####
      tabItem(tabName = "expr_dif",
              
              fluidRow(
                box(
                  width = 4, title = "Opciones para el Volcano Plot", status = "primary", solidHeader = TRUE,
                  radioButtons(
                    inputId = "comparacion",
                    label = "Seleccione el análisis que le gustaría realizar:",
                    choices = list("Vivo vs muerto", "Sano vs enfermo")
                  ),
                  radioButtons(
                    inputId = "muestras",
                    label = "¿Qué conjunto de muestras desea utilizar para el análisis? (sólo para Sano vs enfermo)",
                    choices = list("Totales", "Pareadas")
                  ),
                  numericInput(
                    inputId = "logFC",
                    label = "Establezca el valor mínimo del log Fold Change en valor absoluto",
                    value = 0.5, min = 0, max = 5, step = 0.1
                  ),
                  numericInput(
                    inputId = "p_valor_adj",
                    label = "Establezca el menos logaritmo en base 10 del p-valor ajustado",
                    value = 3, min = 0, max = 1000, step = 1
                  ),
                  helpText("Para ayudarle en los cálculos, x correspondería con un p-valor ajustado de 1e-x"),
                  tags$div(style = "text-align: center;",
                           actionButton(
                             inputId = "analisis_dea",
                             label = "Ejecutar análisis",
                             style = "background-color: #4CAF50; color: white; font-weight: bold; border: none;"
                           )
                  )
                ),
                
                tabBox(
                  width = 8, side = "left",
                  tabPanel("Gráfica resultante", 
                           plotOutput("volcano_plot",
                                      height = "500px",
                                      width = "100%")
                  ),
                  tabPanel("Tabla de microARNs que cumplen ambos filtros",
                          DT::dataTableOutput("mirnas_filtrados")
                  )
                )
              )
      ),
      
      #### Análisis de supervivencia (ui) ####
      tabItem(tabName = "surv",
              
              fluidRow(
                box(
                  width = 4, title = "Parámetros modelo de Cox", status = "primary", solidHeader = TRUE,
                  numericInput(
                    inputId = "p_valor",
                    label = "Establezca el p-valor a usar",
                    value = 0.05, min = 0, max = 1, step = 0.01
                  ),
                  sliderInput(
                    inputId = "minHR",
                    label = "HR mínimo",
                    value = 0.85, min = 0, max = 1, step = 0.01
                  ),
                  tags$div(class = "slider-custom",
                           sliderInput(
                             inputId = "maxHR",
                             label = "HR máximo",
                             value = 1.15, min = 1, max = 2, step = 0.01
                           )
                  ),
                  tags$div(
                    style = "text-align: center;",
                    helpText("Se escogerán los microARN cuyo HR (Hazard Ratio) quede dentro de los intervalos azules")
                  ),
                  tags$div(style = "text-align: center;",
                           actionButton(
                             inputId = "analisis_cox",
                             label = "Ejecutar análisis",
                             style = "background-color: #4CAF50; color: white; font-weight: bold; border: none;"
                           )
                  )
                ),
                
                box(
                  width = 8, title = "Lista de miARN filtrados", status = "info", solidHeader = TRUE,
                  uiOutput("mirnas")
                ),
                
                box(
                  width = 8, title = "Resultados", status = "success", solidHeader = TRUE,
                  plotOutput("kaplan_meier",
                             height = "500px",
                             width = "100%")
                )
              )
      ),
      
      
      #### Análisis de enriquecimiento [genes diana] (ui) ####
      
      tabItem(tabName = "enrich",
              
              fluidRow(
                tabBox(
                  width = 5, side = "left",
                  tabPanel(
                    "Genes diana",
                    radioButtons(
                      inputId = "mirna_genes",
                      label = "Seleccione los microARN a los que desea obtener sus genes diana:",
                      choices = list("Los resultantes tras Expresión diferencial", "Los resultantes tras Análisis de supervivencia", "La unión de ambos"),
                      selected = "La unión de ambos"
                    ),
                    radioButtons(
                      inputId = "tipo_info",
                      label = "Escoja el tipo de información que le gustaría obtener de la base de datos:",
                      choices = list("Interacciones validadas experimentalmente", "Predicciones computacionales", "Ambas")
                    ),
                    uiOutput("siguiente_seccion")
                  ),
                  tabPanel(
                    "Análisis de enriquecimiento",
                    radioButtons(
                      inputId = "tipo_enrich",
                      label = "Escoja el tipo de análisis:",
                      choices = list("GO","KEGG")
                    ),
                    numericInput(
                      inputId = "p_valor_enrich",
                      label = "Seleccione el p-valor",
                      value = 0.05, min = 0, max = 1, step = 0.01
                    ),
                    numericInput(
                      inputId = "q_valor_enrich",
                      label = "Seleccione el q-valor o p-valor ajustado",
                      value = 0.2, min = 0, max = 1, step = 0.01
                    ),
                    selectInput(
                      inputId = "ont_enrich",
                      label = "Seleccione la ontología (sólo para GO)",
                      choices = list("Proceso biológico", 
                                     "Función molecular",
                                     "Componente celular"), 
                      selected = "Proceso biológico"
                    ),
                    numericInput(
                      inputId = "max_categorias",
                      label = "¿Cuántas categorías desea visualizar?",
                      value = 10, min = 1, max = 15, step = 1
                    ),
                    tags$div(style = "text-align: center;",
                             actionButton(
                               inputId = "enrich_analysis",
                               label = "Ejecutar análisis",
                               style = "background-color: #4CAF50; color: white; font-weight: bold; border: none;"
                             )
                    )
                  )
                ),
                box(
                  width = 7, title = "Resultados",
                  status = "success", solidHeader = TRUE,
                  plotOutput("enrich_barras",
                             height = "500px",
                             width = "100%")
                )
              )
      ),
      
      
      #### Clustering (ui) ####
      
      tabItem(tabName = "cluster", 
              
              fluidRow(
                tabBox(
                  width = 12, side = "left",
                  tabPanel(
                    title = "Parámetros dendrograma",
                    fluidRow(
                      column(4, 
                             div(class = "center-col",
                                 radioButtons(
                                   inputId = "mirnas_escogidos",
                                   label = "Seleccione el conjunto de microARN a usar para su agrupación",
                                   choices = list("microARN resultantes de los análisis de Expresión y Supervivencia", "microARN con más variabilidad (ordenados de mayor a menor según coeficiente de variación)")
                                 )
                             )
                      ),
                      column(4, 
                             div(class = "center-col",
                                 numericInput(
                                   inputId = "num_mirnas",
                                   label = "En el segundo caso, seleccione los N primeros microARN que le gustaría visualizar",
                                   value = 50, min = 1, max = 100, step = 1
                                 )
                             )
                      ),
                      column(4, 
                             div(class = "center-col",
                                 selectInput(
                                   inputId = "jerarquia",
                                   label = "Seleccione el método de enlace que permitirá medir las distancias entre los microARN",
                                   choices = list("ward.D2", "complete", "average")
                                 )
                             )
                      )
                    ),
                    fluidRow(
                      column(4, offset = 4,
                             div(style = "text-align: center;",
                               actionButton(
                                 inputId = "dendrograma",
                                 label = "Ver dendrograma",
                                 style = "background-color: #4CAF50; color: white; font-weight: bold; border: none;"
                               )
                             )
                      )
                    )
                  ),
                  tabPanel(
                    title = "Parámetros heatmap",
                    fluidRow(
                      column(6, 
                             div(class = "center-col",
                                 sliderInput(
                                   inputId = "k_clustering",
                                   label = "Establezca el número de clústeres que desea formar",
                                   value = 4, min = 2, max = 6, step = 1
                                 )
                             )
                      ),
                      column(6, 
                             div(class = "center-col",
                                 numericInput(
                                   inputId = "umbral_heatmap",
                                   label = "Seleccione el valor que limitará por arriba los valores representados en el mapa",
                                   value = 5, min = 1, max = 15, step = 1
                                 )
                             )
                      )
                    ),
                    fluidRow(
                      column(4, offset = 4,
                             div(style = "text-align: center;",
                                 actionButton(
                                   inputId = "heatmap",
                                   label = "Ver mapa de calor",
                                   style = "background-color: #4CAF50; color: white; font-weight: bold; border: none;"
                                 )
                             )
                      )
                    )
                  )
                ),
                tabBox(
                  width = 12, side = "left",
                  tabPanel(
                    "Gráfico resultante",
                    plotOutput("resultado_clustering", 
                               height = "800px", 
                               width = "100%")
                  ),
                  tabPanel(
                    "Agrupamiento por clúster",
                    uiOutput("parametro_clinicos_cluster"),
                    uiOutput("resultados_clinicos_cluster"),
                    DT::dataTableOutput("tabla_clinicos_cluster")
                  )
                )
              )
      ),
      
      
      #### Clasificación ML (ui) ####
      
      tabItem(tabName = "clasif_ml", 
              
              fluidRow(
                box(
                  width = 4, title = "Opciones para el modelado", status = "primary", solidHeader = TRUE,
                  numericInput(
                    inputId = "proporcion",
                    label = "Seleccione el porcentaje (entre 0 y 1) en el que le gustaría dividir los conjuntos de entrenamiento y test",
                    value = 0.6, min = 0, max = 1, step = 0.01
                  ),
                  numericInput(
                    inputId = "particiones",
                    label = "Seleccione el número de subconjuntos en el que desea separar el conjunto de entrenamiento dentro de la validación cruzada",
                    value = 10, min = 0, max = 20, step = 1
                  ),
                  selectInput(
                    inputId = "conjunto_datos_a_usar",
                    label = "Seleccione el conjunto de datos",
                    choices = list("Variables clínicas", "Variables clínicas + microARN tras Expresión Diferencial", "Variables clínicas + microARN tras Supervivencia", "Variables clínicas + microARN tras Expresión Diferencial y Supervivencia")
                  ),
                  radioButtons(
                    inputId = "modelo_a_usar",
                    label = "Seleccione el modelo de clasificación a usar",
                    choices = list("Random Forest", "Regresión logística", "Gradient Boosting", "SVM (Support Vector Machines)")
                  ),
                  tags$div(style = "text-align: center;",
                           actionButton(
                             inputId = "ejecucion_tidyverse",
                             label = "Ejecutar con estos valores",
                             style = "background-color: #4CAF50; color: white; font-weight: bold; border: none;"
                           )
                  )
                ),
                
                tabBox(
                  width = 8, side = "left",
                  tabPanel("Métricas obtenidas", 
                           uiOutput("metricas_obtenidas")
                  ),
                  tabPanel("Curva ROC",
                           plotOutput("curva_roc",
                                      height = "600px", 
                                      width = "100%"
                           ),
                           uiOutput("valor_auc")
                  ),
                  tabPanel("Matriz de confusión",
                           plotOutput("matr_conf",
                                      height = "600px", 
                                      width = "100%"
                           )
                  )
                )
              )
      )
      
    )
  )
)


server <- function(input, output, session) {
  
  #### Introducción: Preparación del entorno (server) ####
  
  proyecto <- eventReactive(input$download, {
    input$project
  })
  
  tipo_datos <- eventReactive(input$download, {
    input$data_type
  })
  
  output$proyecto_utilizado <- renderText(
    paste0("→ El tipo de cáncer a analizar es el ", proyecto(), ". Si selecciona otro distinto, es necesario ejecutar todos los análisis que se deseen desde cero.")
  )
  
  
  observeEvent(input$download, {
    
    # Ruta donde se guardaría el archivo simulado
    file_name <- paste0("metadatos/", proyecto(), "/", tipo_datos())
    
    if (file.exists(file_name)) {
      # Caso: ya descargado
      output$status <- renderText(
          paste("Los datos de tipo", tipo_datos() ,"del proyecto", proyecto(), "ya están listos. No hace falta descargar nada.")
      )
    } else {
      # Caso: no descargado → lanzamos descarga de lo que se nos pida
      progressSweetAlert(
        session = session,
        id = "myprogress",
        title = "Descargando datos del GDC...",
        display_pct = TRUE, value = 0
      )
      
      for (i in 1:4) {
        updateProgressBar(
          session = session,
          id = "myprogress",
          value = i * 25
        )
        if(input$data_type == "miRNAs") {
          descarga_mirnas(input$project)
        }else if(input$data_type == "Clinical"){
          descarga_clinicos(input$project)
        }
      }
      
      closeSweetAlert(session = session)
      
      elimina_carpetas_GDC()
      
      output$status <- renderText(
        paste("Datos del proyecto", proyecto(), "descargados correctamente.")
      )
    } 
  })
  
  
  output$proyecto_a_usar <- renderUI({
    titulo <- "Se usará el último proyecto que se haya descargado u obtenido y que contenga ambos tipos de datos. En nuestro caso:"
    file_name_cli <- paste0("metadatos/", proyecto(), "/Clinical")
    file_name_mir <- paste0("metadatos/", proyecto(), "/miRNAs")
    
    if(file.exists(file_name_cli) & file.exists(file_name_cli)){
      items <- c(paste("Clínicos: ", proyecto()),
                 paste("Expresiones de miRNA o miARN: ", proyecto())
      )
    }else if(!file.exists(file_name_cli)){
      shinyalert(
        title = "¡Atención!",
        text = paste("Faltan los datos clínicos del proyecto", proyecto(), ". Si quiere usar todas las funcionalidades de la aplicación, se recomienda que los descargue."),
        type = "warning"
      )
      items <- c("Clínicos: Faltan", paste("Expresiones de miRNA o miARN: ", proyecto()))
    }else{
      shinyalert(
        title = "¡Atención!",
        text = paste("Faltan los datos de expresiones de microARN del proyecto", proyecto(), ". Si quiere usar todas las funcionalidades de la aplicación, se recomienda que los descargue."),
        type = "warning"
      )
      items <- c(paste("Clínicos: ", proyecto()), "Expresiones de miRNA o miARN: Faltan")
    }
    HTML(paste0(
      "<strong>", titulo, "</strong>",  # título en negrita
      "<ul>",
      paste("<li>", items, "</li>", collapse = ""),
      "</ul>"
    ))
  })
  
  
  #### Preprocesado mirnas (server) ####
  
  
  df_resultado_pre_mirnas <- eventReactive(input$preprocesado_mirnas, {
    preprocesado_mirnas(
      df = crea_df_mirnas(proyecto()),
      project = proyecto(),
      porcentaje_max_nulos = input$porcentaje_max_nulos,
      porcentaje_max_na = input$porcentaje_max_na, 
      valor_varianza = input$valor_varianza
    )
  })
  
  num_mirnas_selecc <- reactive({
    nrow(df_resultado_pre_mirnas())
  })
  
  output$estado_pre_mirnas <- renderUI(
    "Por ahora nada. Recuerde descargar u obtener los datos del proyecto deseado para su análisis en la ventana de 'Introducción' antes de continuar"
  )
  
  
  observeEvent(input$preprocesado_mirnas, {
    
    validate(
      need(input$download, "⚠️ Recuerde descargar u obtener los datos del proyecto deseado para su análisis antes de continuar")
    )
    
    sendSweetAlert(
      session = session,
      title = "Espere un poco...",
      text = "Se está ejecutando el preprocesado con los parámetros elegidos. Suele tardar entre 2 y 4 minutos.",
      type = "info"
    )
    
    notif_mirnas <- showNotification("Espere un poco...", duration = NULL, type = "warning")
    
    if(nrow(df_resultado_pre_mirnas()) > 0){
      output$estado_pre_mirnas <- renderUI(
        HTML(paste0("<strong>Datos preprocesados correctamente. Hay ", 
                    num_mirnas_selecc(), " microARN.</strong>"))
      )
      closeSweetAlert(session)
      removeNotification(notif_mirnas)
    }
  })
  
  
  #### Preprocesado clínico (server) ####
  
  observeEvent(input$preprocesado_clinico, {
    
    validate(
      need(input$download, "⚠️ Debe primero descargar u obtener los datos del proyecto deseado desde la ventana de 'Introducción'")
    )
    
    sendSweetAlert(
      session = session,
      title = "Espere un poco...",
      text = "Se está ejecutando el filtrado de características. Suele tardar entre 1 y 2 minutos.",
      type = "info"
    )
    
    notif_clinica <- showNotification("Espere un poco...", duration = NULL, type = "warning")
    
    df_resultado_pre_clinico()
    
    if(ncol(df_resultado_pre_clinico()) == length(lista_signif())){
      closeSweetAlert(session)
      removeNotification(notif_clinica)
    }
  })
  
  output$pre_clinico_seleccion <- renderUI({
    texto_principal <- "Se ha optado por una selección de características manual. En ella se consideran los siguientes aspectos:"
    elementos <-  c("De todas las columnas clínicas se eliminan las que tengan un porcentaje elevado de valores faltantes. Se ha escogido el 60% como umbral inferior.",
                    "Después se filtran aquellas que tengan relación con la que se considera la variable objetivo, el estado vital del paciente. La relación se comprueba a partir de test paramétricos, no paramétricos o de independencia. Si ésta es estadísticamente significativa (p-valor menor que 0.05) se selecciona la columna.",
                    "Por último, de las columnas resultantes se comprueba si alguna está fuertemente relacionada con otra mediante la prueba V de Cramer. Si fuese el caso se procede a eliminar una de ellas, buscando evitar así la multicolinealidad.")
      
    HTML(paste0(
      "<p> <strong>", texto_principal, "</strong> </p>",
      "<ul>",
      paste("<p> <li>", elementos, "</li> </p>", collapse = ""),
      "</ul>"
    ))
  })
  
  df_resultado_pre_clinico <- eventReactive(input$preprocesado_clinico, {
    preprocesado_clinico(
      df = crea_df_clinico(proyecto())
    )[[1]]
  })
  
  lista_signif <- eventReactive(input$preprocesado_clinico, {
    validate(
      need(input$download, "⚠️ Debe primero descargar u obtener los datos del proyecto deseado desde la ventana de 'Introducción'")
    )
    preprocesado_clinico(
      df = crea_df_clinico(proyecto())
    )[[2]]
  })
  
  output$columnas_signif <- renderUI({
    selectInput("caracteristica",
                "Elija una de ellas para ver su análisis y su relación con la variable objetivo 'vital_status'",
                choices = lista_signif()
    )
  })
  
  grafica_caracteristica_selecc <- reactive({
    grafica_relacion_vitalstatus(
      df = df_resultado_pre_clinico(),
      col = input$caracteristica
    )
  })
  
  output$grafica_pre_clinico <- renderPlot({
    req(input$caracteristica)
    grafica_caracteristica_selecc()
  })
  
  datos_univariados_caracteristica_selecc <- reactive({
    datos_analisis_univariado(
      df = df_resultado_pre_clinico(),
      col = input$caracteristica
    )
  })
  
  output$univariado <- renderPrint({
    req(input$caracteristica)
    datos_univariados_caracteristica_selecc()
  })
  
  
  #### Combinación de datos clínicos y de expresiones ####
  
  df_conjunto <- reactive({
    req(input$preprocesado_mirnas)
    req(input$preprocesado_clinico)
    limpieza_final(
      transformacion_variables(
        procesar_datos_mirna_y_clinicos(
          df_resultado_pre_mirnas(),
          df_resultado_pre_clinico()
        ),
        df_resultado_pre_clinico()
      ), 
      df_resultado_pre_mirnas()
    )[[1]]
  })
  
  df_resultado_pre_mirnas_filtrado <- reactive({
    req(input$preprocesado_mirnas)
    req(input$preprocesado_clinico)
    limpieza_final(
      transformacion_variables(
        procesar_datos_mirna_y_clinicos(
          df_resultado_pre_mirnas(),
          df_resultado_pre_clinico()
        ),
        df_resultado_pre_clinico()
      ), 
      df_resultado_pre_mirnas()
    )[[2]]
  })
  
  
  #### Análisis de expresión diferencial (server) ####
  
  resultados_dea <- eventReactive(input$analisis_dea, {
    validate(
      need(input$preprocesado_mirnas, "⚠️ Debe ejecutar primero el preprocesado de los microARN. Una vez lo haga vuelva a ejecutar el análisis.")
    )
    validate(
      need(input$preprocesado_clinico, "⚠️ Debe ejecutar primero el preprocesado clínico. Una vez lo haga vuelva a ejecutar el análisis.")
    )
    analisis_expr_dif(
      muestras_escogidas = input$muestras,
      comparacion = input$comparacion,
      df = df_resultado_pre_mirnas(),
      df_completo = df_conjunto()
    )
  })
  
  grafica_dea <- eventReactive(input$analisis_dea, {
    volcano(
      tabla_resultados = resultados_dea(),
      umbral_logFC = input$logFC, 
      umbral_pvalor_ajustado = input$p_valor_adj
    )
  })
  
  output$volcano_plot <- renderPlot({
    req(input$analisis_dea)
    grafica_dea()
  })
  
  resultados_dea_tras_filtro <- eventReactive(input$analisis_dea, {
    validate(
      need(input$preprocesado_mirnas, "⚠️ Debe ejecutar primero el preprocesado de los microARN. Una vez lo haga vuelva a ejecutar el análisis.")
    )
    validate(
      need(input$preprocesado_clinico, "⚠️ Debe ejecutar primero el preprocesado clínico. Una vez lo haga vuelva a ejecutar el análisis.")
    )
    mirnas_que_cumplen_ambos_filtros(
      resultado = resultados_dea(),
      umbral_logFC = input$logFC, 
      umbral_pvalor_ajustado = input$p_valor_adj
    )
  })
  
  output$mirnas_filtrados <- DT::renderDataTable({
    req(input$analisis_dea)
    resultados_dea_tras_filtro()
  }, options = list(
      pageLength = 10,        # número de filas visibles
      scrollY = "400px",      # altura fija
      scrollCollapse = TRUE,  # que no genere espacio extra si hay pocas filas
      scrollX = TRUE          # scroll horizontal
  ))
  
  
  
  #### Análisis de supervivencia (server) ####
  
  
  resultados_cox <- eventReactive(input$analisis_cox, {
    validate(
      need(input$preprocesado_mirnas, "⚠️ Debe ejecutar primero el preprocesado de los microARN. Una vez lo haga vuelva a ejecutar el análisis.")
    )
    validate(
      need(input$preprocesado_clinico, "⚠️ Debe ejecutar primero el preprocesado clínico. Una vez lo haga vuelva a ejecutar el análisis.")
    )
    mirnas_modelo_cox(
      p_value = input$p_valor,
      hazard_ratio_inferior = input$minHR,
      hazard_ratio_superior = input$maxHR,
      df_mirnas = df_resultado_pre_mirnas_filtrado(),
      df_completo = df_conjunto()
    )[[1]]
  })
  
  hr_cox <- eventReactive(input$analisis_cox, {
    validate(
      need(input$preprocesado_mirnas, "⚠️ Debe ejecutar primero el preprocesado de los microARN. Una vez lo haga vuelva a ejecutar el análisis.")
    )
    validate(
      need(input$preprocesado_clinico, "⚠️ Debe ejecutar primero el preprocesado clínico. Una vez lo haga vuelva a ejecutar el análisis.")
    )
    mirnas_modelo_cox(
      p_value = input$p_valor,
      hazard_ratio_inferior = input$minHR,
      hazard_ratio_superior = input$maxHR,
      df_mirnas = df_resultado_pre_mirnas_filtrado(),
      df_completo = df_conjunto()
    )[[2]]
  })
  
  observeEvent(input$analisis_cox, {
    
    validate(
      need(input$preprocesado_mirnas, "⚠️ Debe ejecutar primero el preprocesado de los microARN. Una vez lo haga vuelva a ejecutar el análisis.")
    )
    validate(
      need(input$preprocesado_clinico, "⚠️ Debe ejecutar primero el preprocesado clínico. Una vez lo haga vuelva a ejecutar el análisis.")
    )
    
    notif_supervivencia <- showNotification("Espere un poco...", duration = NULL, type = "warning")
    
    resultados_cox()
    
    if(length(resultados_cox()) > 0){
      removeNotification(notif_supervivencia)
    }
  })
  
  output$mirnas <- renderUI({
    selectInput("columna",
                "Seleccione uno de los microARN",
                choices = resultados_cox()
    )
  })
  
  grafica_mirna <- reactive({
    curva_kaplan_meier(
      lista_mirnas = resultados_cox(),
      lista_hr = hr_cox(),
      mirna_escogido = input$columna,
      df_completo = df_conjunto()
    )
  })
  
  output$kaplan_meier <- renderPlot({
    req(input$columna)
    grafica_mirna()
  })
  
  
  #### Análisis de enriquecimiento [genes diana] (server) ####
  
  mirnas_seleccionados <- reactive({
    validate(
      need(input$analisis_dea, "⚠️ Debe ejecutar primero el análisis de Expresión diferencial. Una vez lo haga vuelva a ejecutar este análisis.")
    )
    validate(
      need(input$analisis_cox, "⚠️ Debe ejecutar primero el análisis de supervivencia. Una vez lo haga vuelva a ejecutar este análisis.")
    )
    if(input$mirna_genes == "Los resultantes tras Análisis de supervivencia"){
      resultados_cox()
    }else if(input$mirna_genes == "Los resultantes tras Expresión diferencial"){
      rownames(resultados_dea_tras_filtro())
    }else if(input$mirna_genes == "La unión de ambos"){
      union(resultados_cox(), rownames(resultados_dea_tras_filtro()))
    }
  })
  
  genes_diana <- reactive({
    resultados_genes_diana(
      mirnas_seleccionados = mirnas_seleccionados(),
      tipo_datos = input$tipo_info
    )
  })
  
  output$siguiente_seccion <- renderUI({
    HTML("<strong>Una vez configurada esta parte, pase a la siguiente pestaña de esta misma ventana ↑</strong>")
  })
  
  datos_enriquecimiento <- eventReactive(input$enrich_analysis, {
    validate(
      need(input$analisis_dea, "⚠️ Debe ejecutar primero el análisis de Expresión diferencial. Una vez lo haga vuelva a ejecutar este análisis.")
    )
    validate(
      need(input$analisis_cox, "⚠️ Debe ejecutar primero el análisis de supervivencia. Una vez lo haga vuelva a ejecutar este análisis.")
    )
    analisis_enriquecimiento(
      genes_diana = genes_diana(),
      tipo_analisis = input$tipo_enrich,
      p_valor = input$p_valor_enrich,
      q_valor = input$q_valor_enrich,
      ontologia = input$ont_enrich
    )
  })
  
  observeEvent(input$enrich_analysis, {
    validate(
      need(input$analisis_dea, "⚠️ Debe ejecutar primero el análisis de Expresión diferencial. Una vez lo haga vuelva a ejecutar este análisis.")
    )
    validate(
      need(input$analisis_cox, "⚠️ Debe ejecutar primero el análisis de supervivencia. Una vez lo haga vuelva a ejecutar este análisis.")
    )
    
    notif_enriquecimiento <- showNotification("Espere un poco...", duration = NULL, type = "warning")
    
    datos_enriquecimiento()
    
    if(nrow(datos_enriquecimiento()) > 0){
      removeNotification(notif_enriquecimiento)
    }
  })
  
  tipo <- eventReactive(input$enrich_analysis, {
    input$tipo_enrich
  })
  
  categorias <- eventReactive(input$enrich_analysis, {
    input$max_categorias
  })
  
  output$enrich_barras <- renderPlot({
    req(input$enrich_analysis)
    diagrama_barras(
      enrich_object = datos_enriquecimiento(),
      max_categorias = categorias(),
      titulo = paste("Diagrama de barras para el análisis", tipo())
    )
  })
  
  
  
  #### Clustering (server) ####
  
  
  # guardamos qué botón fue pulsado
  grafico_a_mostrar <- reactiveVal(NULL)
  
  observeEvent(input$dendrograma, {
    grafico_a_mostrar("dendrograma")
  })
  
  observeEvent(input$heatmap, {
    grafico_a_mostrar("heatmap")
  })
  
  union_mirnas <- reactive({
    validate(
      need(input$analisis_dea, "⚠️ Debe ejecutar primero el análisis de Expresión diferencial. Una vez lo haga vuelva a ejecutar este análisis.")
    )
    validate(
      need(input$analisis_cox, "⚠️ Debe ejecutar primero el análisis de supervivencia. Una vez lo haga vuelva a ejecutar este análisis.")
    )
    union(resultados_cox(), rownames(resultados_dea_tras_filtro()))
  })
  
  resultado_dendrograma <- eventReactive(input$dendrograma, {
    validate(
      need(input$preprocesado_mirnas, "⚠️ Debe ejecutar primero el preprocesado de los microARN. Una vez lo haga vuelva a ejecutar el análisis.")
    )
    validate(
      need(input$preprocesado_clinico, "⚠️ Debe ejecutar primero el preprocesado clínico. Una vez lo haga vuelva a ejecutar el análisis.")
    )
    clustering_dendrograma(
      df_mirnas = df_resultado_pre_mirnas_filtrado(),
      df_completo = df_conjunto(),
      tipo = input$mirnas_escogidos,
      dea_y_cox = union_mirnas(),
      N = input$num_mirnas,
      metodo_jerarquico = input$jerarquia
    )
  })
    
  resultado_heatmap <- eventReactive(input$heatmap, {
    validate(
      need(input$preprocesado_mirnas, "⚠️ Debe ejecutar primero el preprocesado de los microARN. Una vez lo haga vuelva a ejecutar el análisis.")
    )
    validate(
      need(input$preprocesado_clinico, "⚠️ Debe ejecutar primero el preprocesado clínico. Una vez lo haga vuelva a ejecutar el análisis.")
    )
    clustering_heatmap(
      mirnas_escalados = resultado_dendrograma()[[1]],
      conglomerado = resultado_dendrograma()[[2]],
      k = input$k_clustering,
      umbral = input$umbral_heatmap
    )
  })
  
  output$resultado_clustering <- renderPlot({
    req(grafico_a_mostrar())
    
    if (grafico_a_mostrar() == "dendrograma") {
      resultado_dendrograma()
    } else if (grafico_a_mostrar() == "heatmap") {
      resultado_heatmap()
    }
  })
  
  output$parametro_clinicos_cluster <- renderUI({
    fluidRow(
      column(3, 
             div(class = "center-col",
                 selectInput(
                   inputId = "num_cluster_selecc",
                   label = "Seleccione el cluster con el que desee trabajar",
                   choices = setNames(1:input$k_clustering, paste("Cluster", 1:input$k_clustering)),
                   selected = 1
                 )
             )
      ),
      column(3,
             div(class = "center-col",
                 radioButtons(
                   inputId = "tipo_regulated",
                   label = "Indique con qué tipo de microARNs del clúster escogido le gustaría quedarse",
                   choices = list("up-regulated (color rojo)", "down-regulated (color azul)")
                 )
             )
      ),
      column(3,
             div(class = "center-col",
                 numericInput(
                   inputId = "limite_selecc",
                   label = "Establezca el número de muestras donde los microARN del clúster escogido se encuentran diferencialmente expresados según el criterio de la izquierda",
                   value = 5, min = 1, max = 50, step = 1
                 )
             )
      ),
      column(3, 
             div(class = "center-col",
                 selectInput(
                   inputId = "columna_selecc",
                   label = "Seleccione la columna clínica que desee visualizar",
                   choices = c("tipo_muestra", lista_signif())
                 )
             )
      )
    )
  })
  
  output$resultados_clinicos_cluster <- renderUI({
    fluidRow(
      tabBox(
        width = 12, side = "left",
        tabPanel(
          "Tabla microARNs",
          div(style = "margin-top: 30px;",
              DT::dataTableOutput("mirnas_cluster")
          )
        ),
        tabPanel(
          "Tabla variables clínicas",
          div(style = "margin-top: 30px;",
              DT::dataTableOutput("variables_clinicas_cluster")
          )
        ),
        tabPanel(
          "Análisis de la variable clínica seleccionada para el clúster",
          plotOutput("analisis_univariado_cluster")
        )
      )
    )
  })
  
  tabla_mirnas_cluster <- reactive({
    tabla_de_datos_por_cada_cluster(
      df_mirnas = df_resultado_pre_mirnas_filtrado(),
      agrupamiento_mirnas = resultado_heatmap(),
      cluster = input$num_cluster_selecc,
      tipo = input$tipo_regulated,
      limite = input$limite_selecc,
      clinical_df = df_resultado_pre_clinico(),
      df_completo = df_conjunto()
    )[[1]]
  })
  
  tabla_clinicas_cluster <- reactive({
    tabla_de_datos_por_cada_cluster(
      df_mirnas = df_resultado_pre_mirnas_filtrado(),
      agrupamiento_mirnas = resultado_heatmap(),
      cluster = input$num_cluster_selecc,
      limite = input$limite_selecc,
      tipo = input$tipo_regulated,
      clinical_df = df_resultado_pre_clinico(),
      df_completo = df_conjunto()
    )[[2]]
  })
  
  muestras_sobreexpresadas <- reactive({
    tabla_de_datos_por_cada_cluster(
      df_mirnas = df_resultado_pre_mirnas_filtrado(),
      agrupamiento_mirnas = resultado_heatmap(),
      cluster = input$num_cluster_selecc,
      limite = input$limite_selecc,
      tipo = input$tipo_regulated,
      clinical_df = df_resultado_pre_clinico(),
      df_completo = df_conjunto()
    )[[3]]
  })
  
  
  output$mirnas_cluster <- renderDataTable({
    tabla_mirnas_cluster()
  }, options = list(
    pageLength = 10,        # número de filas visibles
    scrollY = "800px",      # altura fija
    scrollCollapse = TRUE,  # que no genere espacio extra si hay pocas filas
    scrollX = TRUE,          # scroll horizontal
    fixedColumns = list(leftColumns = 1)
  ))
  
  output$variables_clinicas_cluster <- renderDataTable({
    tabla_clinicas_cluster()
  }, options = list(
    pageLength = 5,        # número de filas visibles
    scrollY = "600px",      # altura fija
    scrollCollapse = TRUE,  # que no genere espacio extra si hay pocas filas
    scrollX = TRUE,          # scroll horizontal
    fixedColumns = list(leftColumns = 1)
  ))
  
  plot_variable_clinica_selecc <- reactive({
    analisis_columna_clinica(
      df_completo = df_conjunto(),
      conj_muestras_sobreexpresadas = muestras_sobreexpresadas(),
      variable_seleccionada = input$columna_selecc
    )
  })
  
  output$analisis_univariado_cluster <- renderPlot({
    plot_variable_clinica_selecc()
  })
  
  
  
  #### Clasificación ML (server) ####
  
  
  modelo_a_usar <- eventReactive(input$ejecucion_tidyverse, {
    input$modelo_a_usar
  })
  
  df_a_usar <- eventReactive(input$ejecucion_tidyverse, {
    df_escogido(
      datos_a_usar = input$conjunto_datos_a_usar,
      df_clinico = df_resultado_pre_clinico(),
      df_completo = df_conjunto(),
      mirnas_dea = rownames(resultados_dea_tras_filtro()),
      mirnas_cox = resultados_cox()
    )
  })
  
  lista_elementos <- reactive({
    validate(
      need(input$analisis_dea, "⚠️ Debe ejecutar primero el análisis de Expresión diferencial. Una vez lo haga vuelva a ejecutar este análisis.")
    )
    validate(
      need(input$analisis_cox, "⚠️ Debe ejecutar primero el análisis de supervivencia. Una vez lo haga vuelva a ejecutar este análisis.")
    )
    pipeline_tidyverse_1(
      df = df_a_usar(),
      p = input$proporcion,
      folds = input$particiones,
      modelo_a_ejecutar = input$modelo_a_usar
    )
  })
  
  resultado_final <- eventReactive(input$ejecucion_tidyverse, {
    validate(
      need(input$analisis_dea, "⚠️ Debe ejecutar primero el análisis de Expresión diferencial. Una vez lo haga vuelva a ejecutar este análisis.")
    )
    validate(
      need(input$analisis_cox, "⚠️ Debe ejecutar primero el análisis de supervivencia. Una vez lo haga vuelva a ejecutar este análisis.")
    )
    pipeline_tidyverse_2(
      model_results_after_tuning = lista_elementos()[[1]],
      model_tune_wf = lista_elementos()[[2]], 
      data_split = lista_elementos()[[3]]
    )
  })
  
  observeEvent(input$ejecucion_tidyverse, {
    validate(
      need(input$analisis_dea, "⚠️ Debe ejecutar primero el análisis de Expresión diferencial. Una vez lo haga vuelva a ejecutar este análisis.")
    )
    validate(
      need(input$analisis_cox, "⚠️ Debe ejecutar primero el análisis de supervivencia. Una vez lo haga vuelva a ejecutar este análisis.")
    )
    
    sendSweetAlert(
      session = session,
      title = "Espere un poco...",
      text = "Se están ejecutando los pasos necesarios para el modelo escogido.",
      type = "info"
    )
    
    notif_modelo <- showNotification("Espere un poco...", duration = NULL, type = "warning")
    
    resultado_final()
    
    if(length(resultado_final()) > 0){
      closeSweetAlert(session)
      removeNotification(notif_modelo)
    }
  })
  
  output$metricas_obtenidas <- renderUI({
    tagList(
      h4(paste("Resultados del modelo",modelo_a_usar(), "para la clasificación de la variable 'vital_status'")),
      tags$ul(
        tags$li(HTML(paste0("<strong>Accuracy:</strong> ", round(resultado_final()[[1]], 3)))),
        tags$li(HTML(paste0("<strong>Precision:</strong> ", round(resultado_final()[[2]], 3)))),
        tags$li(HTML(paste0("<strong>Recall:</strong> ", round(resultado_final()[[3]], 3)))),
        tags$li(HTML(paste0("<strong>F1 Score:</strong> ", round(resultado_final()[[4]], 3)))),
        tags$li(HTML(paste0("<strong>AUC (Area under curve):</strong> ", round(resultado_final()[[5]], 3))))
      )
    )
  })
  
  output$valor_auc <- renderUI(
    HTML(paste("El valor AUC (Area under curve) es de: <strong>",
               round(resultado_final()[[5]],3),
               "</strong>"))
  )
  
  output$curva_roc <- renderPlot(
    resultado_final()[[6]]
  )
  
  output$matr_conf <- renderPlot(
    resultado_final()[[7]]
  )
  
}

shinyApp(ui, server)
