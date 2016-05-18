###Script para an?lise dos MDEs utilizando a de regress?o logistica multinominal e classes de solo
#Predicting Multiple Discrete Values with Multinomials, Neural Networks and the {nnet} Package
#link:http://amunategui.github.io/multinomial-neuralnetworks-walkthrough/

rm(list = ls())

# MDEs avaliados:
# IBGE05
# MDE05
# MDE Pontos (RJ05)
# IBGE20
# MDE20
# MDE-RJ25
# IBGE30
# MDE30
# Topodata

# Buscar os objetos do R feitos e salvos nos processamentos anteriores

# Funções definidas pelo usuário ##############################################################################

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
    correl <- cor(mde[, id])
    
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

### leitura dos dados de calibração da bacia ------------------------------------------------------------------
MDE5a <- read.csv("data/mde05/dados.csv", sep = ";", dec = ",")     # MDE05             5a
MDE5b <- read.csv("data/IBGE05/dados.csv", sep = ";", dec = ",")    # IBGE05            5b
MDE5c <- read.csv("data/mde5/dados.csv", sep = ";", dec = ",")      # MDE Pontos (RJ05) 5c
MDE20a <- read.csv("data/mde20/dados.csv", sep = ";", dec = ",")    # MDE20             20a
MDE20b <- read.csv("data/IBGE20/dados.csv", sep = ";", dec = ",")   # IBGE20            20b
MDE20c <- read.csv("data/RJ/dados.csv", sep = ";", dec = ",")       # MDE-RJ25          20c
MDE30a <- read.csv("data/mde30/dados.csv", sep = ";", dec = ",")    # MDE30             30a
MDE30b <- read.csv("data/IBGE30/dados.csv", sep = ";", dec = ",")   # IBGE30            30b
MDE30c <- read.csv("data/Topodata/dados.csv", sep = ";", dec = ",") # Topodata          30c

# trasformar para fator as variáveis categóricas --------------------------------------------------------------
id <- c("geologia", "Solos", "Flow.Direct")

MDE5a[, id] <- lapply(MDE5a[, id], as.factor)
MDE5b[, id] <- lapply(MDE5b[, id], as.factor)
MDE5c[, id] <- lapply(MDE5c[, id], as.factor)

MDE20a[, id] <- lapply(MDE20a[, id], as.factor)
MDE20b[, id] <- lapply(MDE20b[, id], as.factor)
MDE20c[, id] <- lapply(MDE20c[, id], as.factor)

MDE30a[, id] <- lapply(MDE30a[, id], as.factor)
MDE30b[, id] <- lapply(MDE30b[, id], as.factor)
MDE30c[, id] <- lapply(MDE30c[, id], as.factor)

rm(id)

# Correlação dos dados para cada MDE --------------------------------------------------------------------------
# Calcular a correlação para todos os MDEs a fim de selecionar as covariáveis para entrar no modelo reduzido.
# O valor de corte para a correlação é de 0.85. A lista de covariáveis selecionadas mais curta é usada para
# calibrar o modelo reduzido.
correl <- lapply(
  list(MDE5a, MDE5b, MDE5c, MDE20a, MDE20b, MDE20c, MDE30a, MDE30b, MDE30c), getCorrelation, max = 0.85)
lapply(correl, function (x) x[["correlation"]])
selected <- lapply(correl, function (x) x[["selected"]])
selected <- c(selected[[which.min(sapply(selected, length))]], "Flow.Direct")
selected

# reclassificação para agrupamento das classes ################################################################

reclasse <- as.character(MDE5a$CLASSE)
reclasse <- gsub(pattern = "^CY$", replacement = "CY.RY", x = reclasse)
reclasse <- gsub(pattern = "^RY$", replacement = "CY.RY", x = reclasse)

reclasse <- gsub(pattern = "^LVA$", replacement = "LA.LVA", x = reclasse)
reclasse <- gsub(pattern = "^LA$", replacement = "LA.LVA", x = reclasse)

reclasse <- gsub(pattern = "^MT$", replacement = "MT.NX.RL", x = reclasse)
reclasse <- gsub(pattern = "^NX$", replacement = "MT.NX.RL", x = reclasse)
reclasse <- gsub(pattern = "^RL$", replacement = "MT.NX.RL", x = reclasse)

reclasse <- gsub(pattern = "^PA$", replacement = "PA.PV.PVA", x = reclasse)
reclasse <- gsub(pattern = "^PV$", replacement = "PA.PV.PVA", x = reclasse)
reclasse <- gsub(pattern = "^PVA$", replacement = "PA.PV.PVA", x = reclasse)

reclasse <- gsub(pattern = "^SX$", replacement = "SX.GX", x = reclasse)
reclasse <- gsub(pattern = "^GX$", replacement = "SX.GX", x = reclasse)

MDE5a$reclasse <- as.factor(reclasse)
MDE5b$reclasse <- as.factor(reclasse)
MDE5c$reclasse <- as.factor(reclasse)

MDE20a$reclasse <- as.factor(reclasse)
MDE20b$reclasse <- as.factor(reclasse)
MDE20c$reclasse <- as.factor(reclasse)

MDE30a$reclasse <- as.factor(reclasse)
MDE30b$reclasse <- as.factor(reclasse)
MDE30c$reclasse <- as.factor(reclasse)

par(mfrow = c(2, 1))
plot(MDE5a$CLASSE)
plot(MDE5a$reclasse)
dev.off()

# Ajustar modelos + validação cruzada #########################################################################

# Definir modelos 
base_full <- c("banda1", "banda2", "banda3", "banda4", "banda5", "ndvi", "savi")
base_reduz <- base_full[which(base_full %in% selected)]
full_full <- colnames(MDE5a)
full_full <- full_full[-match(c("reclasse", "CLASSE", "Solos", "geologia", "ID", "X", "Y"), full_full)]
full_reduz <- selected
sapply(list(base_full, base_reduz, full_full, full_reduz), length)

base_full <- formula(paste("reclasse ~ ", paste(base_full, collapse = " + "), sep = ""))
base_reduz <- formula(paste("reclasse ~ ", paste(base_reduz, collapse = " + "), sep = ""))
full_full <- formula(paste("reclasse ~ ", paste(full_full, collapse = " + "), sep = ""))
full_reduz <- formula(paste("reclasse ~ ", paste(full_reduz, collapse = " + "), sep = ""))

# Ajustar modelos
modelo0 <- lapply(list(base_full, base_reduz), nnet::multinom, MDE5a, model = TRUE)

modelo1 <- lapply(list(full_full, full_reduz), nnet::multinom, MDE5a, model = TRUE)
modelo2 <- lapply(list(full_full, full_reduz), nnet::multinom, MDE5b, model = TRUE)
modelo3 <- lapply(list(full_full, full_reduz), nnet::multinom, MDE5c, model = TRUE)

modelo4 <- lapply(list(full_full, full_reduz), nnet::multinom, MDE20a, model = TRUE)
modelo5 <- lapply(list(full_full, full_reduz), nnet::multinom, MDE20b, model = TRUE)
modelo6 <- lapply(list(full_full, full_reduz), nnet::multinom, MDE20c, model = TRUE)

modelo7 <- lapply(list(full_full, full_reduz), nnet::multinom, MDE30a, model = TRUE)
modelo8 <- lapply(list(full_full, full_reduz), nnet::multinom, MDE30b, model = TRUE)
modelo9 <- lapply(list(full_full, full_reduz), nnet::multinom, MDE30c, model = TRUE)

# Acurácia do ajuste
accuracy <- lapply(
  list(modelo0, modelo1, modelo2, modelo3, modelo4, modelo5, modelo6, modelo7, modelo8, modelo9),
  function (x) sapply(x, fitAccuracy))
accuracy <- do.call(rbind, accuracy)
colnames(accuracy) <- c("full", "reduz")
rownames(accuracy) <- paste("modelo", 0:9, sep = "")

# Validação cruzada
validation <- lapply(
  list(modelo0, modelo1, modelo2, modelo3, modelo4, modelo5, modelo6, modelo7, modelo8, modelo9),
  function (x) sapply(x, crossValidation))
validation <- do.call(rbind, validation)
colnames(validation) <- c("full", "reduz")

# Resultado final
res <- cbind(accuracy, validation)
res <- round(res[, c(1, 3, 2, 4)], 1)
res











# carregando os grides de cada MDE ###############################

grid_MDE1 <- read.csv("C:/MDE/dados2/IBGE05/grid.csv", sep=";", dec=",") #IBGE05
grid_MDE2 <- read.csv("C:/MDE/dados2/mde05/grid.csv", sep=";", dec=",") #MDE05 
grid_MDE3 <- read.csv("C:/MDE/dados2/mde5/grid.csv", sep=";", dec=",") #MDE Pontos (RJ05)
grid_MDE4 <- read.csv("C:/MDE/dados2/IBGE20/grid.csv", sep=";", dec=",") #IBGE20
grid_MDE5 <- read.csv("C:/MDE/dados2/mde20/grid.csv", sep=";", dec=",") #MDE20
grid_MDE6 <- read.csv("C:/MDE/dados2/RJ/grid.csv", sep=";", dec=",")    #MDE-RJ25
grid_MDE7 <- read.csv("C:/MDE/dados2/IBGE30/grid.csv", sep=";", dec=",") #IBGE30
grid_MDE8 <- read.csv("C:/MDE/dados2/mde30/grid.csv", sep=";", dec=",")  #MDE30
grid_MDE9 <- read.csv("C:/MDE/dados2/Topodata/grid.csv", sep=";", dec=",") #Topodata

#retirando Nas dos grids
grid_MDE1 <- na.omit(grid_MDE1)
grid_MDE2 <- na.omit(grid_MDE2)
grid_MDE3 <- na.omit(grid_MDE3)
grid_MDE4 <- na.omit(grid_MDE4)
grid_MDE5 <- na.omit(grid_MDE5)
grid_MDE6 <- na.omit(grid_MDE6)
grid_MDE7 <- na.omit(grid_MDE7)
grid_MDE8 <- na.omit(grid_MDE8)
grid_MDE9 <- na.omit(grid_MDE9)

#trasformando as covariaveis categoricas em fator para da MDE
#MDE1
grid_MDE1$uso <- as.factor(grid_MDE1$uso)
grid_MDE1$geologia <- as.factor(grid_MDE1$geologia)
grid_MDE1$Solos <- as.factor(grid_MDE1$Solos)
grid_MDE1$Flow.Direct <- as.factor(grid_MDE1$Flow.Direct)

#MDE2
grid_MDE2$uso <- as.factor(grid_MDE2$uso)
grid_MDE2$geologia <- as.factor(grid_MDE2$geologia)
grid_MDE2$Solos <- as.factor(grid_MDE2$Solos)
grid_MDE2$Flow.Direct <- as.factor(grid_MDE2$Flow.Direct)

#MDE3
grid_MDE3$uso <- as.factor(grid_MDE3$uso)
grid_MDE3$geologia <- as.factor(grid_MDE3$geologia)
grid_MDE3$Solos <- as.factor(grid_MDE3$Solos)
grid_MDE3$Flow.Direct <- as.factor(grid_MDE3$Flow.Direct)

#MDE4
grid_MDE4$uso <- as.factor(grid_MDE4$uso)
grid_MDE4$geologia <- as.factor(grid_MDE4$geologia)
grid_MDE4$Solos <- as.factor(grid_MDE4$Solos)
grid_MDE4$Flow.Direct <- as.factor(grid_MDE4$Flow.Direct)

#MDE5
grid_MDE5$uso <- as.factor(grid_MDE5$uso)
grid_MDE5$geologia <- as.factor(grid_MDE5$geologia)
grid_MDE5$Solos <- as.factor(grid_MDE5$Solos)
grid_MDE5$Flow.Direct <- as.factor(grid_MDE5$Flow.Direct)

#MDE6
grid_MDE6$uso <- as.factor(grid_MDE6$uso)
grid_MDE6$geologia <- as.factor(grid_MDE6$geologia)
grid_MDE6$Solos <- as.factor(grid_MDE6$Solos)
grid_MDE6$Flow.Direct <- as.factor(grid_MDE6$Flow.Direct)

#MDE7
grid_MDE7$uso <- as.factor(grid_MDE7$uso)
grid_MDE7$geologia <- as.factor(grid_MDE7$geologia)
grid_MDE7$Solos <- as.factor(grid_MDE7$Solos)
grid_MDE7$Flow.Direct <- as.factor(grid_MDE7$Flow.Direct)

#MDE8
grid_MDE8$uso <- as.factor(grid_MDE8$uso)
grid_MDE8$geologia <- as.factor(grid_MDE8$geologia)
grid_MDE8$Solos <- as.factor(grid_MDE8$Solos)
grid_MDE8$Flow.Direct <- as.factor(grid_MDE8$Flow.Direct)

#MDE9
grid_MDE9$uso <- as.factor(grid_MDE9$uso)
grid_MDE9$geologia <- as.factor(grid_MDE9$geologia)
grid_MDE9$Solos <- as.factor(grid_MDE9$Solos)
grid_MDE9$Flow.Direct <- as.factor(grid_MDE9$Flow.Direct)



#fazer as predições e salvar em formato .asc
require(sp)
#Testemunha
pred1 <- predict(modelo1, newdata = grid_MDE1, na.action = na.omit, type = "raw")
pred1 <-data.frame(x=grid_MDE1$X, y=grid_MDE1$Y, pred1=pred1)
str(pred1)
coordinates(pred1) <- ~x+y
str(pred1)
gridded(pred1) <- TRUE# transformar "grid_clas" de DataFrame para SpatialPixelsDataFrame
str(pred1)
spplot(pred1)
proj4string(obj = pred1)<-CRS("+init=epsg:32723")
str(pred1)
pred1 <-as(object = pred1, "SpatialGridDataFrame")
writeRaster(raster(pred1), "Testemunha.asc", format="ascii", NAflag=-9999) # esportar em formato asc

#MDE1
pred2 <- predict(modelo2, newdata = grid_MDE1, na.action = na.omit, type = "raw")
pred2 <-data.frame(x=grid_MDE1$X, y=grid_MDE1$Y, pred2=pred2)
str(pred2)
coordinates(pred2) <- ~x+y
str(pred2)
gridded(pred2) <- TRUE# transformar "grid_clas" de DataFrame para SpatialPixelsDataFrame
str(pred2)
spplot(pred2)
proj4string(obj = pred2)<-CRS("+init=epsg:32723")
str(pred2)
pred2 <-as(object = pred2, "SpatialGridDataFrame")
writeRaster(raster(pred2), "MDE1.asc", format="ascii", NAflag=-9999) # esportar em formato asc

#MDE2
pred3 <- predict(modelo3, newdata = grid_MDE2, na.action = na.omit, type = "raw")
pred3 <-data.frame(x=grid_MDE2$X, y=grid_MDE2$Y, pred3=pred3)
str(pred3)
coordinates(pred3) <- ~x+y
str(pred3)
gridded(pred3) <- TRUE# transformar "grid_clas" de DataFrame para SpatialPixelsDataFrame
str(pred3)
spplot(pred3)
proj4string(obj = pred3)<-CRS("+init=epsg:32723")
str(pred3)
pred3 <-as(object = pred3, "SpatialGridDataFrame")
writeRaster(raster(pred3), "MDE2.asc", format="ascii", NAflag=-9999) # esportar em formato asc

#MDE3
pred4 <- predict(modelo4, newdata = grid_MDE3, na.action = na.omit, type = "raw")
pred4 <-data.frame(x=grid_MDE3$X, y=grid_MDE3$Y, pred4=pred4)
str(pred4)
coordinates(pred4) <- ~x+y
str(pred4)
gridded(pred4) <- TRUE# transformar "grid_clas" de DataFrame para SpatialPixelsDataFrame
str(pred4)
spplot(pred4)
proj4string(obj = pred4)<-CRS("+init=epsg:32723")
str(pred4)
pred4 <-as(object = pred4, "SpatialGridDataFrame")
writeRaster(raster(pred4), "MDE3.asc", format="ascii", NAflag=-9999) # esportar em formato asc

#MDE4
pred5 <- predict(modelo5, newdata = grid_MDE4, na.action = na.omit, type = "raw")
pred5 <-data.frame(x=grid_MDE4$X, y=grid_MDE4$Y, pred5=pred5)
str(pred5)
coordinates(pred5) <- ~x+y
str(pred5)
gridded(pred5) <- TRUE# transformar "grid_clas" de DataFrame para SpatialPixelsDataFrame
str(pred5)
spplot(pred5)
proj4string(obj = pred5)<-CRS("+init=epsg:32723")
str(pred5)
pred5 <-as(object = pred5, "SpatialGridDataFrame")
writeRaster(raster(pred5), "MDE4.asc", format="ascii", NAflag=-9999) # esportar em formato asc

#MDE5
pred6 <- predict(modelo6, newdata = grid_MDE5, na.action = na.omit, type = "raw")
pred6 <-data.frame(x=grid_MDE5$X, y=grid_MDE5$Y, pred6=pred6)
str(pred6)
coordinates(pred6) <- ~x+y
str(pred6)
gridded(pred6) <- TRUE# transformar "grid_clas" de DataFrame para SpatialPixelsDataFrame
str(pred6)
spplot(pred6)
proj4string(obj = pred6)<-CRS("+init=epsg:32723")
str(pred6)
pred6 <-as(object = pred6, "SpatialGridDataFrame")
writeRaster(raster(pred6), "MDE5.asc", format="ascii", NAflag=-9999) # esportar em formato asc

#MDE6
pred7 <- predict(modelo7, newdata = grid_MDE6, na.action = na.omit, type = "raw")
pred7 <-data.frame(x=grid_MDE6$X, y=grid_MDE6$Y, pred7=pred7)
str(pred7)
coordinates(pred7) <- ~x+y
str(pred7)
gridded(pred7) <- TRUE# transformar "grid_clas" de DataFrame para SpatialPixelsDataFrame
str(pred7)
spplot(pred7)
proj4string(obj = pred7)<-CRS("+init=epsg:32723")
str(pred7)
pred7 <-as(object = pred7, "SpatialGridDataFrame")
writeRaster(raster(pred7), "MDE6.asc", format="ascii", NAflag=-9999) # esportar em formato asc

#MDE7
pred8 <- predict(modelo8, newdata = grid_MDE7, na.action = na.omit, type = "raw")
pred8 <-data.frame(x=grid_MDE7$X, y=grid_MDE7$Y, pred8=pred8)
str(pred8)
coordinates(pred8) <- ~x+y
str(pred8)
gridded(pred8) <- TRUE# transformar "grid_clas" de DataFrame para SpatialPixelsDataFrame
str(pred8)
spplot(pred8)
proj4string(obj = pred8)<-CRS("+init=epsg:32723")
str(pred8)
pred8 <-as(object = pred8, "SpatialGridDataFrame")
writeRaster(raster(pred8), "MDE7.asc", format="ascii", NAflag=-9999) # esportar em formato asc

#MDE8
pred9 <- predict(modelo9, newdata = grid_MDE8, na.action = na.omit, type = "raw")
pred9 <-data.frame(x=grid_MDE8$X, y=grid_MDE8$Y, pred9=pred9)
str(pred9)
coordinates(pred9) <- ~x+y
str(pred9)
gridded(pred9) <- TRUE# transformar "grid_clas" de DataFrame para SpatialPixelsDataFrame
str(pred9)
spplot(pred9)
proj4string(obj = pred9)<-CRS("+init=epsg:32723")
str(pred9)
pred9 <-as(object = pred9, "SpatialGridDataFrame")
writeRaster(raster(pred9), "MDE8.asc", format="ascii", NAflag=-9999) # esportar em formato asc


#MDE9
pred10 <- predict(modelo10, newdata = grid_MDE9, na.action = na.omit, type = "raw")
pred10 <-data.frame(x=grid_MDE9$X, y=grid_MDE9$Y, pred10=pred10)
str(pred10)
coordinates(pred10) <- ~x+y
str(pred10)
gridded(pred10) <- TRUE# transformar "grid_clas" de DataFrame para SpatialPixelsDataFrame
str(pred10)
spplot(pred10)
proj4string(obj = pred10)<-CRS("+init=epsg:32723")
str(pred10)
pred10 <-as(object = pred10, "SpatialGridDataFrame")
writeRaster(raster(pred10), "MDE9.asc", format="ascii", NAflag=-9999) # esportar em formato asc

                          #The end
                  
