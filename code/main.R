### Script para análise dos MDEs utilizando a de regressão logistica multinominal e classes de solo

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

# Source user defined functions
source("code/helper.R")

# Buscar os objetos do R feitos e salvos nos processamentos anteriores

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
  list(MDE5a, MDE5b, MDE5c, MDE20a, MDE20b, MDE20c, MDE30a, MDE30b, MDE30c), getCorrelation, max = 0.90)
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

# Save calibrated models
save(modelo0, modelo1, modelo2, modelo3, modelo4, modelo5, modelo6, modelo7, modelo8, modelo9,
     file = "data/R/models.rda")

# Acurácia do ajuste
# Aqui usamos a matriz de confusão para calcular a acurácia geral.
accuracy <- lapply(
  list(modelo0, modelo1, modelo2, modelo3, modelo4, modelo5, modelo6, modelo7, modelo8, modelo9),
  function (x) sapply(x, fitAccuracy))
accuracy <- do.call(rbind, accuracy)
colnames(accuracy) <- c("full", "reduz")
rownames(accuracy) <- paste("modelo", 0:9, sep = "")

# Validação cruzada 
# Usar uma observação como observação de validação de cada vez.
# Aqui usamos a matriz de erro para calcular a acurácia geral.
validation <- lapply(
  list(modelo0, modelo1, modelo2, modelo3, modelo4, modelo5, modelo6, modelo7, modelo8, modelo9),
  function (x) sapply(x, crossValidation))
validation <- do.call(rbind, validation)
colnames(validation) <- c("full", "reduz")

# Resultado final
# Agregar resultados da calibração e da validação.
# Salvar os dados para construir tabela para o artigo.
accuracy <- cbind(accuracy, validation)
accuracy <- round(accuracy[, c(1, 3, 2, 4)], 1)
accuracy
write.csv(accuracy, file = "res/tab/accuracy.csv")

# Predição espacial ###########################################################################################

# Carregar o pacote nnet para usar a função 'predict' para objetos 'multinom'.
require(nnet)

# Rampa de cores para a entropia de Shannon
traffic.light <- colorRampPalette(c("olivedrab", "khaki", "maroon1"))

# MDE5a
grid_MDE5a <- read.csv("data/mde05/grid.csv", sep = ";", dec = ",")
grid_MDE5a$Flow.Direct <- as.factor(grid_MDE5a$Flow.Direct)
pred_MDE5a <- data.frame(predict(modelo1[[2]], newdata = grid_MDE5a, type = "probs"))
pred_MDE5a$class <- unlist(predict(modelo1[[2]], newdata = grid_MDE5a, type = "class"))
pred_MDE5a$entropy <- -rowSums(pred_MDE5a[, 1:6] * log(pred_MDE5a[, 1:6], base = 6), na.rm = TRUE)
pred_MDE5a <- cbind(grid_MDE5a[, c("X", "Y")], pred_MDE5a)
sp::gridded(pred_MDE5a) <- ~ X + Y
sp::proj4string(obj = pred_MDE5a) <- sp::CRS("+init=epsg:32723")
# sp::spplot(pred_MDE5a, "class")
# sp::spplot(pred_MDE5a, "entropy", at = seq(0, 1, 0.1), col.regions = traffic.light(10))
save(pred_MDE5a, file = "data/R/pred_MDE5a.rda")
rm(pred_MDE5a)
gc()

# Testemunha
# Usar o mesmo grid de predição carregado acima.
pred_teste <- data.frame(predict(modelo0[[2]], newdata = grid_MDE5a, type = "probs"))
pred_teste$class <- unlist(predict(modelo0[[2]], newdata = grid_MDE5a, type = "class"))
pred_teste$entropy <- -rowSums(pred_teste[, 1:6] * log(pred_teste[, 1:6], base = 6), na.rm = TRUE)
pred_teste <- cbind(grid_MDE5a[, c("X", "Y")], pred_teste)
sp::gridded(pred_teste) <- ~ X + Y
sp::proj4string(pred_teste) <- sp::CRS("+init=epsg:32723")
# sp::spplot(pred_teste, "class")
# sp::spplot(pred_teste, "entropy", at = seq(0, 1, 0.1), col.regions = traffic.light(10))
save(pred_teste, file = "data/R/pred_teste.rda")
rm(pred_teste, grid_MDE5a)
gc()

# MDE5b
grid_MDE5b <- read.csv("data/IBGE05/grid.csv", sep = ";", dec = ",")
grid_MDE5b$Flow.Direct <- as.factor(grid_MDE5b$Flow.Direct)
pred_MDE5b <- data.frame(predict(modelo2[[2]], newdata = grid_MDE5b, type = "probs"))
pred_MDE5b$class <- unlist(predict(modelo2[[2]], newdata = grid_MDE5b, type = "class"))
pred_MDE5b$entropy <- -rowSums(pred_MDE5b[, 1:6] * log(pred_MDE5b[, 1:6], base = 6), na.rm = TRUE)
pred_MDE5b <- cbind(grid_MDE5b[, c("X", "Y")], pred_MDE5b)
sp::gridded(pred_MDE5b) <- ~ X + Y
sp::proj4string(pred_MDE5b) <- sp::CRS("+init=epsg:32723")
# sp::spplot(pred_MDE5b, "class")
# sp::spplot(pred_MDE5b, "entropy", at = seq(0, 1, 0.1), col.regions = traffic.light(10))
save(pred_MDE5b, file = "data/R/pred_MDE5b.rda")
rm(pred_MDE5b, grid_MDE5b)
gc()

# MDE5c
grid_MDE5c <- read.csv("data/mde5/grid.csv", sep = ";", dec = ",")
grid_MDE5c$Flow.Direct <- as.factor(grid_MDE5c$Flow.Direct)
pred_MDE5c <- data.frame(predict(modelo3[[2]], newdata = grid_MDE5c, type = "probs"))
pred_MDE5c$class <- unlist(predict(modelo3[[2]], newdata = grid_MDE5c, type = "class"))
pred_MDE5c$entropy <- -rowSums(pred_MDE5c[, 1:6] * log(pred_MDE5c[, 1:6], base = 6), na.rm = TRUE)
pred_MDE5c <- cbind(grid_MDE5c[, c("X", "Y")], pred_MDE5c)
sp::gridded(pred_MDE5c) <- ~ X + Y
sp::proj4string(pred_MDE5c) <- sp::CRS("+init=epsg:32723")
# sp::spplot(pred_MDE5c, "class")
# sp::spplot(pred_MDE5c, "entropy", at = seq(0, 1, 0.1), col.regions = traffic.light(10))
save(pred_MDE5c, file = "data/R/pred_MDE5c.rda")
rm(pred_MDE5c, grid_MDE5c)
gc()

# MDE20a
grid_MDE20a <- read.csv("data/mde20/grid.csv", sep = ";", dec = ",")
grid_MDE20a$Flow.Direct <- as.factor(grid_MDE20a$Flow.Direct)
pred_MDE20a <- data.frame(predict(modelo4[[2]], newdata = grid_MDE20a, type = "probs"))
pred_MDE20a$class <- unlist(predict(modelo4[[2]], newdata = grid_MDE20a, type = "class"))
pred_MDE20a$entropy <- -rowSums(pred_MDE20a[, 1:6] * log(pred_MDE20a[, 1:6], base = 6), na.rm = TRUE)
pred_MDE20a <- cbind(grid_MDE20a[, c("X", "Y")], pred_MDE20a)
sp::gridded(pred_MDE20a) <- ~ X + Y
sp::proj4string(pred_MDE20a) <- sp::CRS("+init=epsg:32723")
# sp::spplot(pred_MDE20a, "class")
# sp::spplot(pred_MDE20a, "entropy", at = seq(0, 1, 0.1), col.regions = traffic.light(10))
save(pred_MDE20a, file = "data/R/pred_MDE20a.rda")
rm(pred_MDE20a, grid_MDE20a)
gc()

# MDE20b
grid_MDE20b <- read.csv("data/IBGE20/grid.csv", sep = ";", dec = ",")
grid_MDE20b$Flow.Direct <- as.factor(grid_MDE20b$Flow.Direct)
pred_MDE20b <- data.frame(predict(modelo5[[2]], newdata = grid_MDE20b, type = "probs"))
pred_MDE20b$class <- unlist(predict(modelo5[[2]], newdata = grid_MDE20b, type = "class"))
pred_MDE20b$entropy <- -rowSums(pred_MDE20b[, 1:6] * log(pred_MDE20b[, 1:6], base = 6), na.rm = TRUE)
pred_MDE20b <- cbind(grid_MDE20b[, c("X", "Y")], pred_MDE20b)
sp::gridded(pred_MDE20b) <- ~ X + Y
sp::proj4string(pred_MDE20b) <- sp::CRS("+init=epsg:32723")
# sp::spplot(pred_MDE20b, "class")
# sp::spplot(pred_MDE20b, "entropy", at = seq(0, 1, 0.1), col.regions = traffic.light(10))
save(pred_MDE20b, file = "data/R/pred_MDE20b.rda")
rm(pred_MDE20b, grid_MDE20b)
gc()

# MDE20c
grid_MDE20c <- read.csv("data/RJ/grid.csv", sep = ";", dec = ",")
grid_MDE20c$Flow.Direct <- as.factor(grid_MDE20c$Flow.Direct)
pred_MDE20c <- data.frame(predict(modelo6[[2]], newdata = grid_MDE20c, type = "probs"))
pred_MDE20c$class <- unlist(predict(modelo6[[2]], newdata = grid_MDE20c, type = "class"))
pred_MDE20c$entropy <- -rowSums(pred_MDE20c[, 1:6] * log(pred_MDE20c[, 1:6], base = 6), na.rm = TRUE)
pred_MDE20c <- cbind(grid_MDE20c[, c("X", "Y")], pred_MDE20c)
sp::gridded(pred_MDE20c) <- ~ X + Y
sp::proj4string(pred_MDE20c) <- sp::CRS("+init=epsg:32723")
# sp::spplot(pred_MDE20c, "class")
# sp::spplot(pred_MDE20c, "entropy", at = seq(0, 1, 0.1), col.regions = traffic.light(10))
save(pred_MDE20c, file = "data/R/pred_MDE20c.rda")
rm(pred_MDE20c, grid_MDE20c)
gc()

# MDE30a
grid_MDE30a <- read.csv("data/mde30/grid.csv", sep = ";", dec = ",")
grid_MDE30a$Flow.Direct <- as.factor(grid_MDE30a$Flow.Direct)
pred_MDE30a <- data.frame(predict(modelo7[[2]], newdata = grid_MDE30a, type = "probs"))
pred_MDE30a$class <- unlist(predict(modelo7[[2]], newdata = grid_MDE30a, type = "class"))
pred_MDE30a$entropy <- -rowSums(pred_MDE30a[, 1:6] * log(pred_MDE30a[, 1:6], base = 6), na.rm = TRUE)
pred_MDE30a <- cbind(grid_MDE30a[, c("X", "Y")], pred_MDE30a)
sp::gridded(pred_MDE30a) <- ~ X + Y
sp::proj4string(pred_MDE30a) <- sp::CRS("+init=epsg:32723")
# sp::spplot(pred_MDE30a, "class")
# sp::spplot(pred_MDE30a, "entropy", at = seq(0, 1, 0.1), col.regions = traffic.light(10))
save(pred_MDE30a, file = "data/R/pred_MDE30a.rda")
rm(pred_MDE30a, grid_MDE30a)
gc()

# MDE30b
grid_MDE30b <- read.csv("data/IBGE30/grid.csv", sep = ";", dec = ",")
grid_MDE30b$Flow.Direct <- as.factor(grid_MDE30b$Flow.Direct)
pred_MDE30b <- data.frame(predict(modelo8[[2]], newdata = grid_MDE30b, type = "probs"))
pred_MDE30b$class <- unlist(predict(modelo8[[2]], newdata = grid_MDE30b, type = "class"))
pred_MDE30b$entropy <- -rowSums(pred_MDE30b[, 1:6] * log(pred_MDE30b[, 1:6], base = 6), na.rm = TRUE)
pred_MDE30b <- cbind(grid_MDE30b[, c("X", "Y")], pred_MDE30b)
sp::gridded(pred_MDE30b) <- ~ X + Y
sp::proj4string(pred_MDE30b) <- sp::CRS("+init=epsg:32723")
# sp::spplot(pred_MDE30b, "class")
# sp::spplot(pred_MDE30b, "entropy", at = seq(0, 1, 0.1), col.regions = traffic.light(10))
save(pred_MDE30b, file = "data/R/pred_MDE30b.rda")
rm(pred_MDE30b, grid_MDE30b)
gc()

# MDE30c
grid_MDE30c <- read.csv("data/Topodata/grid.csv", sep = ";", dec = ",")
grid_MDE30c$Flow.Direct <- as.factor(grid_MDE30c$Flow.Direct)
pred_MDE30c <- data.frame(predict(modelo9[[2]], newdata = grid_MDE30c, type = "probs"))
pred_MDE30c$class <- unlist(predict(modelo9[[2]], newdata = grid_MDE30c, type = "class"))
pred_MDE30c$entropy <- -rowSums(pred_MDE30c[, 1:6] * log(pred_MDE30c[, 1:6], base = 6), na.rm = TRUE)
pred_MDE30c <- cbind(grid_MDE30c[, c("X", "Y")], pred_MDE30c)
sp::gridded(pred_MDE30c) <- ~ X + Y
sp::proj4string(pred_MDE30c) <- sp::CRS("+init=epsg:32723")
# sp::spplot(pred_MDE30c, "class")
# sp::spplot(pred_MDE30c, "entropy", at = seq(0, 1, 0.1), col.regions = traffic.light(10))
save(pred_MDE30c, file = "data/R/pred_MDE30c.rda")
rm(pred_MDE30c, grid_MDE30c)
gc()

# Preparar figuras ############################################################################################

# Preparar os nomes dos arquivos das figuras
files <- expand.grid(c(5, 20, 30), letters[1:3])
files <- files[order(files$Var1), ]
files <- 
  c(sapply(1:nrow(files), function (i) paste("pred_MDE", files[i, 1], files[i, 2], sep = "")), "pred_teste")
file_names <- paste("load('data/R/", files, ".rda')", sep = "")
file_names <- sapply(1:length(file_names), function (i) parse(text = file_names[i]))
for (i in 1:length(file_names)) {
  eval(file_names[[i]])
}

# Preparar figuras com predições
maps <- list()
for (i in 1:length(files)) {
  obj <- parse(text = paste("maps[[i]] <- ", files[i], "$class", sep = ""))
  eval(obj)
}

# Número de píxeis por classe
tmp <- lapply(maps, summary)
tmp <- lapply(1:10, function (i) tmp[[i]][1:6])
tmp <- as.data.frame(tmp)
colnames(tmp) <- 
  c(paste(rep(paste("DEM", c(5, 20, 30), sep = ""), each = 3), c("a", "b", "c"), sep = ""), "Baseline")
nam <- rownames(tmp)
tmp <- stack(tmp)
tmp$um <- rep(nam, 10)
col <- sp::bpy.colors(length(unique(tmp$um)))
p <- lattice::barchart(
  values ~ um | ind, tmp, layout = c(3, 4), col = col, 
  scales = list(y = list(alternating = FALSE), x = list(draw = FALSE)),
  key = list(corner = c(0.9, 0.9), text = list(unique(tmp$um)), rectangles = list(col = col)),
  panel = function (...) {
    lattice::panel.grid(v = -1, h = -1)
    lattice::panel.barchart(...)
  })
p$index.cond[[1]] <- c(5:7, 2:4, 8:10, 1)
dev.off()
png("res/fig/pixels.png", width = 480 * 2.5, height = 480 * 3, res = 150)
p
dev.off()
rm(p)

# Continuar com os mapas
maps <- cbind(pred_MDE30a@coords, as.data.frame(maps))
col_names <- gsub("pred_", "", files)
col_names <- gsub("MDE", "DEM", col_names)
col_names[length(col_names)] <- "Baseline"
colnames(maps) <- c(colnames(maps)[1:2], col_names)
sp::gridded(maps) <- ~ X + Y
sp::proj4string(maps) <- sp::CRS("+init=epsg:32723")
p <- sp::spplot(
  maps, layout = c(3, 4),
  colorkey = FALSE,
  key =
    list(corner = c(0.9, 0.9), text = list(unique(tmp$um)), rectangles = list(col = col)),
  panel = function (...) {
    lattice::panel.grid(v = -1, h = -1)
    lattice::panel.levelplot(...)
  })
p$index.cond[[1]] <- c(7:9, 4:6, 1:3, 10)
dev.off()
png("res/fig/predictions.png", width = 480*2, height = 480*3.5, res = 150)
p
dev.off()
rm(p)
gc()

# Preparas figuras com entropia de Shannon
maps <- list()
for (i in 1:length(files)) {
  obj <- parse(text = paste("maps[[i]] <- ", files[i], "$entropy", sep = ""))
  eval(obj)
}

# Distribuição de frequência
tmp <- data.frame(maps)
colnames(tmp) <- 
  c(paste(rep(paste("DEM", c(5, 20, 30), sep = ""), each = 3), c("a", "b", "c"), sep = ""), "Baseline")
tmp <- stack(tmp)
p <- lattice::histogram(
  ~ values | ind, tmp, layout = c(3, 4), nint = 15, col = traffic.light(15), xlab = "Entropy")
p$index.cond[[1]] <- c(5:7, 2:4, 8:10, 1)
dev.off()
png("res/fig/entropy_histogram.png", width = 480 * 2.5, height = 480 * 3, res = 150)
p
dev.off()
rm(p, tmp)

# Continuar com os mapas
maps <- cbind(pred_MDE30a@coords, as.data.frame(maps))
col_names <- gsub("pred_", "", files)
col_names <- gsub("MDE", "DEM", col_names)
col_names[length(col_names)] <- "Baseline"
colnames(maps) <- c(colnames(maps)[1:2], col_names)

# Preparas mapas
sp::gridded(maps) <- ~ X + Y
sp::proj4string(maps) <- sp::CRS("+init=epsg:32723")
p <- sp::spplot(
  maps, layout = c(3, 4),
  col.regions = traffic.light,
  panel = function (...) {
    lattice::panel.grid(v = -1, h = -1)
    lattice::panel.levelplot(...)
  })
p$index.cond[[1]] <- c(7:9, 4:6, 1:3, 10)
names(p$legend) <- "inside"
p$legend$inside$x <- 0.795
p$legend$inside$y <- 0.875
p$legend$inside$args$key$height <- 0.2
# maps$index.cond[[1]] <- c(7:9, 4:6, 1:3, 10)
dev.off()
png("res/fig/entropy.png", width = 480*2, height = 480*3.5, res = 150)
p
dev.off()
