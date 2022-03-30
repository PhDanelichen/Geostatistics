# Disciplina de Geoestatistica
# Universidade de Cuiaba
# Professores: Victor Danelichen
# Aula - Geoestatistica

# instalar pacotes e verficar diretorio de trabalho
update.packages(ask=FALSE)
install.packages("raster")
install.packages("sf")
install.packages("sp")
install.packages("rgdal")
install.packages("tmap")
install.packages("PROJ")#acrescentado
library(raster) # 2º
library(sf) # 3º
library(sp) #executar 1º
library(rgdal) 
library(tmap)
library(PROJ)
## outros pacotes que vamos usar ao longo da aula
install.packages("dismo")
install.packages("stars")
install.packages("gstat")
install.packages("fields")
install.packages("rcompanion")
install.packages("bestNormalize")
install.packages("geoR")
install.packages("automap")
library(dismo)
library(stars)
library(gstat)
library(fields)
library(rcompanion)
library(bestNormalize)
library(geoR)
library(automap)

getwd()
setwd("C:/Aula_Geoestatistica/aula9")

## Importando dados
st_layers("aula9.gpkg") #comando que le o arquivo
mun <- st_read("aula9.gpkg", layer="mun_abc")
st_crs(mun)<-31983 
estacoes <- st_read("aula9.gpkg", layer="estacoes")
st_crs(estacoes)<-31983
View(estacoes)
tm_shape(estacoes) + tm_dots(col="chuva", pal="Blues") +
tm_shape(mun) + tm_borders(col="black")

# Vizinho mais proximo (instalar deldir)
install.packages("deldir")
library(dismo)
library(deldir)
library(maptools)
install.packages("gpclib")
library(gpclib)
estacoes_sp <- as(estacoes,"Spatial")
voronoi <- voronoi(estacoes_sp)
voronoi_sf <- st_as_sf(voronoi)
View(voronoi_sf)
plot(st_geometry(voronoi_sf))
plot(st_geometry(estacoes), pch=20, cex=0.4, col="blue", add=TRUE)
library("tmap") # acrescentado
tm_shape(voronoi_sf) + 
    tm_fill(col="chuva", style="quantile", pal="Blues") + tm_borders() +
 tm_shape(mun) + tm_borders(col="black")
tm_shape(voronoi_sf) + 
    tm_fill(col="temperat", style="quantile", pal="Reds") + tm_borders() +
  tm_shape(mun) + tm_borders(col="black")

# Inverso da Distancia

modelo_raster <- raster("srtm_abc.tif", values=FALSE)
modelo_raster #com 30m de resolucao espacial
modelo_raster <- aggregate(modelo_raster, fact = 10)
modelo_raster #com 300 m
library(stars)
library(gstat)
st_crs(estacoes) <- crs(modelo_raster)
gs_idw_2 <- gstat(formula=chuva~1, locations=estacoes, set=list(idp=4))
cv_idw_2 <- gstat.cv(gs_idw_2)
View(as.data.frame(cv_idw_2))
sqrt(mean(cv_idw_2$residual^2))                        # raiz do erro medio quadrado
1-(var(cv_idw_2$residual)/var(cv_idw_2$observed))      # estimativa da porcentagem da variacao explicada (R2)
chuva_idw_2 <- interpolate(object=modelo_raster, model=gs_idw_2) # Funcao do pacote raster - demora um pouco 
plot(chuva_idw_2)
plot(st_geometry(mun), add=TRUE)

#Exercicio 1 (scrip anterior) - Faca a interpolacao do inverso da distancia com pesos 1 e 6, comparando os mapas, a validacao cruzada e o R2

### Otimizacao
erro_medio_quadrado <- function(peso) {
                       gs <- gstat(formula=chuva~1, locations=estacoes, set=list(idp=peso))
                       cv <- gstat.cv(gs)
                       sqrt(mean(cv$residual^2))}
erro_medio_quadrado(2)
peso_otimo <- optimize(f=erro_medio_quadrado, interval=c(1,10)) # demora um pouco 
peso_otimo$minimum
peso_otimo$objective
gs_otimo <- gstat(formula=chuva~1, locations=estacoes, set=list(idp=peso_otimo$minimum))
cv_otimo <- gstat.cv(gs_otimo)
1-(var(cv_otimo$residual)/var(cv_otimo$observed))
chuva_idw_otimo <- interpolate(modelo_raster, gs_otimo)
plot(chuva_idw_otimo)
plot(st_geometry(mun), add=TRUE)

# Superficie de tendencia
gs_superficie_1 <- gstat(formula=chuva~1, location=estacoes, degree=1)
cv_superficie_1 <- gstat.cv(gs_superficie_1)
sqrt(mean(cv_superficie_1$residual^2))
1-(var(cv_superficie_1$residual)/var(estacoes$chuva))
superficie_1 <- interpolate(modelo_raster, gs_superficie_1)
plot(superficie_1)
plot(st_geometry(mun), add=TRUE)
persp(superficie_1, border=NA, col="blue", theta = 320, phi=40)

# Exercicio 2 (script anterior)- Fazer a superficie de tendencia de segunda ordem, avaliando a validacao cruzada, o R2 e visualizando o mapa e a perspectiva 3D

## Funcao de base radial - Thin plate spline
library(fields)
modelo_tps <- Tps(x=st_coordinates(estacoes), Y=estacoes$chuva)
sqrt(mean(modelo_tps$residual^2))                    # root mean squared error
1-(var(modelo_tps$residual)/var(estacoes$chuva))     # % de explicacao do modelo
chuva_tps <- interpolate(modelo_raster, modelo_tps)
plot(chuva_tps)
plot(st_geometry(mun), add=TRUE)
persp(chuva_tps, border=NA, col="blue", shade=0.3, theta = 320, phi=40)

# Exercicio 3 - Fazer uma interpolacao spline para os dados de temperatura, avaliando a validacao cruzada, e visualizando o mapa e a perspectiva 3D

#Krigagem
##normalizacao de variavel
library(rcompanion)
plotNormalHistogram(estacoes_sp$chuva)
library(bestNormalize)
modelo_normal <- bestNormalize(estacoes_sp$chuva)
modelo_normal
plotNormalHistogram(modelo_normal$x.t)
estacoes_sp$chuva_n <- modelo_normal$x.t
chuva_original <- predict(modelo_normal, newdata = modelo_normal$x.t, inverse=TRUE)
plotNormalHistogram(chuva_original)

##Variograma
variograma <- variogram(chuva_n~1, estacoes_sp)
plot(variograma)
variograma <- variogram(chuva_n~1, estacoes_sp, cutoff=100000)
plot(variograma, plot.numbers=TRUE)

##ajuste visual de variograma
library("geoR")
estacoes_geodata <- as.geodata(estacoes_sp[,"chuva_n"])
variograma_geor <- variog(estacoes_geodata)
plot(variograma_geor)
dev.new()
eyefit_geor <- eyefit(variograma_geor)

# Exercicio 4 - Faca o ajuste visual do modelo sobre o variograma e exporte o grafico final.

dev.off()
eyefit_geor
plot(variograma_geor)
lines(eyefit_geor)
variofit_geor <- variofit(variograma_geor, ini.cov.pars = eyefit_geor)
variofit_geor
plot(variograma_geor)
lines(variofit_geor)

## Ajuste automatico de variograma
library(automap)
variofit_geor$nugget         # efeito pepita
variofit_geor$cov.pars[1]    # patamar
variofit_geor$cov.pars[2]    # alcance estruturado (range)
variogram_auto <- autofitVariogram(chuva_n~1, estacoes_sp, start_vals=c(variofit_geor$nugget, variofit_geor$cov.pars[2], variofit_geor$cov.pars[1]))
View(variogram_auto$var_model)
plot(variogram_auto)

## Interpolacao por krigagem
library('PROJ')
writeRaster(modelo_raster,"modelo_raster.tif")
modelo_grid <- readGDAL("modelo_raster.tif")
crs(estacoes_sp) <- crs(modelo_grid)
krigagem_auto <- autoKrige(chuva_n~1, estacoes_sp, new_data=modelo_grid, start_vals=c(variofit_geor$nugget, variofit_geor$cov.pars[2], variofit_geor$cov.pars[1]))
plot(krigagem_auto)
par(mfrow=c(1,1))
plot(krigagem_auto$krige_output["var1.pred"])
plot(krigagem_auto$krige_output["var1.stdev"])
par(mfrow=c(1,1))
krigagem_auto$krige_output$original <- predict(modelo_normal, newdata=krigagem_auto$krige_output$var1.pred, inverse=TRUE)
plot(krigagem_auto$krige_output["original"])
krigagem_raster <- raster(krigagem_auto$krige_output["original"])
par(mfrow=c(1,1))
plot(krigagem_raster)
plot(st_geometry(mun), add=TRUE)
persp(krigagem_raster, border=NA, col="blue", shade=1, theta = 320, phi=40)
krigagem_cv <- autoKrige.cv(chuva_n~1, estacoes_sp, start_vals=c(variofit_geor$nugget, variofit_geor$cov.pars[2], variofit_geor$cov.pars[1]))
sqrt(mean(krigagem_cv$krige.cv_output$residual^2))
1-(var(krigagem_cv$krige.cv_output$residual)/var(krigagem_cv$krige.cv_output$observed))

#Exercicio 5 -  faca uma krigagem dos dados de temperatura, visualize o mapa de predicao e de desvio padrao, e analise a validacao cruzada

## Krigagem indicadora
chuva_2000 <- predict(modelo_normal, newdata=2000)
chuva_2000
fit_indicadora <- autofitVariogram(formula=I(chuva_n<chuva_2000)~1, estacoes_sp)
plot(fit_indicadora)
fit_indicadora$var_model
krigagem_indicadora <- krige(formula=I(chuva_n<chuva_2000)~1, locations=estacoes_sp, model=fit_indicadora$var_model, nmax=12, newdata=modelo_grid)
plot(krigagem_indicadora["var1.pred"])
plot(st_geometry(mun), add=TRUE)

# Exercicio 6 - Faca uma krigagem indicadora de uma regiao do mapa ter menos do que 18 graus centigrados

# Regression-kriging
srtm <- raster("srtm_abc.tif")
srtm <- aggregate(srtm, fact=10)
par(mfrow=c(1,1))
plot(srtm)
estacoes_sp$elevacao <- extract(srtm, estacoes_sp)
plot(estacoes_sp$temperat ~ estacoes_sp$elevacao)
lm_temperatura <- lm(data=estacoes_sp, temperat ~ elevacao)
summary(lm_temperatura)
plot(estacoes_sp$temperat ~ estacoes_sp$elevacao)
abline(lm_temperatura, col="red", lwd=2)
variograma_temp_el <- variogram(temperat~elevacao, estacoes_sp)
plot(variograma_temp_el)
variogram_fit_temp_el <- autofitVariogram(temperat ~ elevacao, estacoes_sp)
plot(variogram_fit_temp_el)
writeRaster(srtm, "srtm_agregado.tif")
srtm_grid <- readGDAL("srtm_agregado.tif")
View(srtm_grid@data)
colnames(srtm_grid@data) <- "elevacao"
View(srtm_grid@data)
krigagem_temp_el = autoKrige(temperat ~ elevacao, estacoes_sp, new_data = srtm_grid)
plot(krigagem_temp_el)
krigagem_cv_temp_el <- autoKrige.cv(temperat ~ elevacao, estacoes_sp)
sqrt(mean(krigagem_cv_temp_el$krige.cv_output$residual^2))
1-(var(krigagem_cv_temp_el$krige.cv_output$residual)/var(krigagem_cv_temp_el$krige.cv_output$observed))

# Exercicio 7 - Faca uma krigagem regressiva universal, usando elevacao, latitude e longitude para prever a temperatura.
##              Dica: formula = temperat ~ elevacao + x + y
##              O gstat entende o “x” como latitude e o “y” como longitude

# Anisotropia
## Mapa de semivariancia
mapa_variograma <- variogram(chuva_n~1, estacoes_sp, cutoff=90000, width=2000, map=TRUE)
plot(mapa_variograma)
variograma_anis <- variogram(chuva_n~1, estacoes_sp, cutoff=90000, alpha=c(55, 65, 75, 85, 95, 130, 140, 150, 160, 170))
plot(variograma_anis)
variograma_anis <- variogram(chuva_n~1, estacoes_sp, cutoff=90000, alpha=c(75, 150))
plot(variograma_anis)

# Exercicio 8 - Avalie a anisotropia para os dados de temperatura

# Gravar os dados antes de sair
save.image("D:/R_CTA/aula9/aula9.R.RData")
