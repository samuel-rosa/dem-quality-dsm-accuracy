# USER DEFINED FUNCTIONS ######################################################################################

# Validação cruzada dos modelos ajustados
crossValidation <- 
  function (model, digits = 2) {
    
    # Criar vetor para armazenar os resultados
    cv <- vector()
    
    # Definir indice sobre o qual o loop atuará
    idx <- nrow(model.frame(model))
    
    # Iniciar loop
    for (i in 1:idx) {
      
      # Ajustar modelo de regressão logística multinomial múltipla deixando de fora a i-ésima observação
      tmp <- nnet::multinom(formula = formula(model), data = model.frame(model)[-i, ])
      
      # Predizer a classe taxonômica na i-ésima observação usando o modelo ajustado acima
      cv[i] <- predict(object = tmp, newdata = model.frame(model)[i, ])
    }
    
    # Construir a matriz de erro usando o valor predito e o valor observado
    obs <- model.frame(model)[, 1]
    cm <- table(cv, obs)
    
    # Calcular a acurácia global, ou seja, o percentual de valores na diagonal da matriz de error
    res <- round(sum(diag(cm)) / nrow(model.frame(model)) * 100, digits = digits)
    
    return (res)
  }

# Computar a acurácia da calibração
fitAccuracy <- 
  function (model, digits = 2) {
    
    res <- predict(model, model.frame(model))
    cm <- table(res, model.frame(model)[, 1])
    res <- round(sum(diag(cm)) / nrow(model.frame(model)) * 100, digits = digits)
    return (res)
  }

# Identificar as covariáveis com correlação superior a determinado limite
getCorrelation <- 
  function (mde, max = 0.85, digits = 2) {
    
    # Identificar as covariáveis que são do tipo 'factor'
    id <- colnames(mde)[!sapply(mde, is.factor)]
    
    # Identificar as colunas contendo a identificação e coordenadas das observações
    id <- id[!id %in% c("ID", "X", "Y")]
    
    # Calcular a correlação entre as covariáveis
    correl <- abs(cor(mde[, id]))
    
    # Identificar os pares de covariáveis com correlação superior à 'max' 
    # ou igual à 1, ou seja, da covariável com ela mesma
    idx <- which(correl > max & correl < 1, arr.ind = TRUE)
    covar <- data.frame(x = rownames(correl)[idx[, 1]], y = colnames(correl)[idx[, 2]])
    
    # Adicionar coluna com a correlação e ordenar da maior para a menor
    res <- cbind(covar, correl = correl[idx])
    res <- res[rev(order(res$correl)), ]
    
    # Remover comparações repetidas e arredondar
    res <- res[seq(1, nrow(res), 2), ]
    res <- data.frame(res, row.names = NULL)
    res$correl <- round(res$correl, digits = digits)
    
    # Selecionar covariáveis para entrar no modelo (correlação < 'max')
    id <- id[-match(res$y, id)]
    res <- list(correlation = res, selected = id)
    
    return (res)
  }
