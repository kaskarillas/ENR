---
title: "S1_bone"
author: "JFGG"
date: "2024-09-10"
output: html_document
---


# Sintaxis para realizar una huella metabolomica de la masa ósea considerando:
  - Mujeres
  - Menos de 51

## Definicion de las BBDD:
MO.rankNorm = BM v2 y todos los metabolitos de v2
MO1.rankNorm = BM v2 y metabolitos v2 menos acilcarnitinas
MO2.rankNorm = delta BM y todos los metabolitos de v2
MO.deltas.rN1 = deltas de BM v8-v2 y metabolitos v8 - v2
MO.RL = BBDD con el score de MO2.rankNorm, MO en cat, MO basa, edad, sex, BMI basal.

## 1. In less than 50 years old women to analyse, at baseline, associations between metabolites and bone mass using elastic net regression treating bone mass as a continuous variable.

```{r}
#Apertura de las BBDD:
library(readxl)
library(xlsx)
library(missForest)
library(RNOmni)
library(caTools)
library(glmnet)
library(pROC)
library(rio)
library(haven)
library(agricolae)
library(tidyverse)
library(caret)
library(magrittr)
library(dplyr)
library(gridExtra)
```

##Analisis transversal con todos los metabolitos (v0):

```{r}
colnames(bq)
bq1 = bq[,c("allocation_no","v2_C0","v2_C2.0","v2_C3.0","v2_C4.1","v2_C3.0.2M","v2_C4.0","v2_C4.OH",  
            "v2_C5.OH","v2_C6.OH","v2_C4.1.2M","v2_C4.0.2M","v2_C4.0.3M","v2_C6.0","v2_C7.0","v2_C4.DC",    
            "v2_C3.DC","v2_C3.DC.M","v2_C5.DC","v2_C6.DC","v2_C8.1","v2_C5.M.DC","v2_C8.0","v2_C9.0",      
            "v2_C10.2","v2_C10.1","v2_C10.0","v2_C11.0","v2_C12.1","v2_C12.0","v2_C7.DC","v2_C14.2",      
            "v2_C13.0","v2_C14.1","v2_C14.0","v2_C15.0","v2_C16.OH","v2_C16.2","v2_C16.1","v2_C18.2",       
            "v2_C16.0","v2_C18.1","v2_C18.0","v2_5-HT","v2_5-HIAA","v2_TMAO","v2_FChr","v2_Echr",        
            "v2_TChr","v2_TG_RMN","v2_PhCh","v2_LCh","v2_Sph","v2_FAC","v2_w3","v2_ARA+EPA",
            "v2_DHA","v2_LINOLEIC","v2_PUFA","v2_MUFA","v2_LPC_14.0","v2_LPC_15.0","v2_LPC_16.0","v2_LPC_16.0.e", 
            "v2_LPC_16.1","v2_LPC_16.1.e","v2_LPC_17.0","v2_LPC_18.0","v2_LPC_18.0.e","v2_LPC_18.1","v2_LPC_18.2","v2_LPC_20.0",    
            "v2_LPC_20.1","v2_LPC_20.3","v2_LPC_20.4","v2_LPC_22.6","v2_PC_30.0","v2_PC_32.0","v2_PC_32.1","v2_PC_32.1.e",   
            "v2_PC_32.2","v2_PC_33.1","v2_PC_34.0","v2_PC_34.1","v2_PC_34.1.e","v2_PC_34.2","v2_PC_34.2.e","v2_PC_34.3.e",   
            "v2_PC_34.4","v2_PC_35.1","v2_PC_35.2","v2_PC_36.1","v2_PC_36.2","v2_PC_36.2.e","v2_PC_36.3","v2_PC_36.4",  
            "v2_PC_36.4.e","v2_PC_36.5","v2_PC_36.5.e","v2_PC_37.4","v2_PC_38.3","v2_PC_38.4","v2_PC_38.4.e","v2_PC_38.5",    
            "v2_PC_38.5.e","v2_PC_38.6","v2_PC_40.4","v2_PC_40.4.e","v2_PC_40.5.e","v2_PC_40.6","v2_PC_42.5.e","v2_PE_36.5.e",   
            "v2_PE_38.5.e","v2_PE_38.6.e","v2_SM_32.1","v2_SM_32.2","v2_SM_33.1","v2_SM_34.1","v2_SM_34.2","v2_SM_35.1",    
            "v2_SM_36.0","v2_SM_36.1","v2_SM_36.2","v2_SM_38.1","v2_SM_38.2","v2_SM_40.1","v2_SM_40.2","v2_SM_41.1",    
            "v2_SM_41.2","v2_SM_42.1","v2_SM_42.2","v2_SM_42.3","v2_TG_48.1","v2_TG_48.2","v2_TG_50.1","v2_TG_50.2",    
            "v2_TG_50.3","v2_TG_50.4","v2_TG_51.2","v2_TG_52.1","v2_TG_52.2","v2_TG_52.3","v2_TG_52.4","v2_TG_52.5",    
            "v2_TG_54.2","v2_TG_54.3","v2_TG_54.4","v2_TG_54.5","v2_LA","v2_GlyA","v2_Ala","v2_Glycine",
            "v2_2-HbutA","v2_3-HbutA","v2_Val","v2_Leu","v2_Glycerol","v2_Isoleu","v2_Proline","v2_GliA",        
            "v2_Ser","v2_Thr","v2_Meth","v2_Orn","v2_GlutA","v2_Phe","v2_Lys","v2_CitA",        
            "v2_Fruc","v2_Glu","v2_Tyr","v2_PalA","v2_LinoA","v2_OleicA","v2_StearicA","v2_Tryp", 
            "v2_Sucrose","v2_alphaToco","v2_Choles")]

bq2 = bq[,c("allocation_no","v2_FChr","v2_Echr",        
            "v2_TChr","v2_TG_RMN","v2_PhCh","v2_LCh","v2_Sph","v2_FAC","v2_w3","v2_ARA+EPA",
            "v2_DHA","v2_LINOLEIC","v2_PUFA","v2_MUFA","v2_LPC_14.0","v2_LPC_15.0","v2_LPC_16.0","v2_LPC_16.0.e", 
            "v2_LPC_16.1","v2_LPC_16.1.e","v2_LPC_17.0","v2_LPC_18.0","v2_LPC_18.0.e","v2_LPC_18.1","v2_LPC_18.2","v2_LPC_20.0",    
            "v2_LPC_20.1","v2_LPC_20.3","v2_LPC_20.4","v2_LPC_22.6","v2_PC_30.0","v2_PC_32.0","v2_PC_32.1","v2_PC_32.1.e",   
            "v2_PC_32.2","v2_PC_33.1","v2_PC_34.0","v2_PC_34.1","v2_PC_34.1.e","v2_PC_34.2","v2_PC_34.2.e","v2_PC_34.3.e",   
            "v2_PC_34.4","v2_PC_35.1","v2_PC_35.2","v2_PC_36.1","v2_PC_36.2","v2_PC_36.2.e","v2_PC_36.3","v2_PC_36.4",  
            "v2_PC_36.4.e","v2_PC_36.5","v2_PC_36.5.e","v2_PC_37.4","v2_PC_38.3","v2_PC_38.4","v2_PC_38.4.e","v2_PC_38.5",    
            "v2_PC_38.5.e","v2_PC_38.6","v2_PC_40.4","v2_PC_40.4.e","v2_PC_40.5.e","v2_PC_40.6","v2_PC_42.5.e","v2_PE_36.5.e",   
            "v2_PE_38.5.e","v2_PE_38.6.e","v2_SM_32.1","v2_SM_32.2","v2_SM_33.1","v2_SM_34.1","v2_SM_34.2","v2_SM_35.1",    
            "v2_SM_36.0","v2_SM_36.1","v2_SM_36.2","v2_SM_38.1","v2_SM_38.2","v2_SM_40.1","v2_SM_40.2","v2_SM_41.1",    
            "v2_SM_41.2","v2_SM_42.1","v2_SM_42.2","v2_SM_42.3","v2_TG_48.1","v2_TG_48.2","v2_TG_50.1","v2_TG_50.2",    
            "v2_TG_50.3","v2_TG_50.4","v2_TG_51.2","v2_TG_52.1","v2_TG_52.2","v2_TG_52.3","v2_TG_52.4","v2_TG_52.5",    
            "v2_TG_54.2","v2_TG_54.3","v2_TG_54.4","v2_TG_54.5","v2_LA","v2_GlyA","v2_Ala","v2_Glycine",
            "v2_2-HbutA","v2_3-HbutA","v2_Val","v2_Leu","v2_Glycerol","v2_Isoleu","v2_Proline","v2_GliA",        
            "v2_Ser","v2_Thr","v2_Meth","v2_Orn","v2_GlutA","v2_Phe","v2_Lys","v2_CitA",        
            "v2_Fruc","v2_Glu","v2_Tyr","v2_PalA","v2_LinoA","v2_OleicA","v2_StearicA","v2_Tryp", 
            "v2_Sucrose","v2_alphaToco","v2_Choles")]
colnames(BBDD[4000:4236])
BBDD1 = BBDD[,c("allocation_no","v2_bone_mass","v8_bone_mass","v14_bone_mass","v1_gq_11", "v1_gq_1")]
BBDD2 = merge(BBDD1, bq1, by = "allocation_no")
BBDD3 = merge(BBDD1, bq2, by = "allocation_no")
```

Filtro por sexo y edad:

```{r}
table(BBDD2$v1_gq_1) #184 mujeres, 47 hombres

table(BBDD2$v1_gq_11)

BBDD2 = subset(BBDD2, (v1_gq_11 <= 50 & v1_gq_1 == 1) | v1_gq_1 == 2) #Se eliminan 74 por ser mayores de 50 --> 115 mujeres, 47 hombres
BBDD3 = subset(BBDD3, (v1_gq_11 <= 50 & v1_gq_1 == 1) | v1_gq_1 == 2) #Se eliminan 74 por ser mayores de 50 --> 115 mujeres, 47 hombres
```

Modelos continuos:

Quality control:
```{r}
BBDD2$delta = BBDD2$v8_bone_mass - BBDD2$v2_bone_mass

MO = cbind(BBDD2[2], BBDD2[7:184])
MO1 = cbind(BBDD3[2], BBDD3[7:139])
MO2 = cbind(BBDD2$delta, BBDD2[7:184])
MO3 = cbind(BBDD2[1:4], BBDD2[5:6], BBDD2$delta)
MO3.1 = data.frame(cbind(BBDD$allocation_no, BBDD$v2_4.1, BBDD$v2_4.2, BBDD$v2_4.3, BBDD$v2_4.4))
colnames(MO3.1) = c("allocation_no", "v2_4.1", "v2_4.2", "v2_4.3", "v2_4.4")
MO3 = merge(MO3, MO3.1, by = "allocation_no")

MO = subset(MO, v2_bone_mass > 0) #Se eliminan valores NA. 2
MO1 = subset(MO1, v2_bone_mass > 0) #Se eliminan valores NA. 2
MO2 = subset(MO2, BBDD2$delta > -10) #Se eliminan valores NA. 44
MO3 = subset(MO3, BBDD2$delta > -10)  #Se eliminan valores NA. 44

#NAs en cada variables
imp1 = MO[2:179]
imp2 = MO1[2:134]
imp3 = MO2[2:179]

na1 = sapply(imp3, function(imp3) sum(length(which(is.na(imp3)))))
na2 = sapply(imp3, function(imp3) (100*sum(length(which(is.na(imp3))))/sum(length((imp3)))))
na =  cbind("NA"=na1, "% NA"=na2)
na

#Metabolites con NA mayor de 20%
A_NA=rownames(na[na[,2]>20 ,])

#ReMGveMGs > 20%
drop = A_NA
imp1 = imp1[,!(names(imp1) %in% drop)] #10
sapply(imp1, function(imp1) sum(length(which(is.na(imp1)))))
sapply(imp1, function(imp1) (100*sum(length(which(is.na(imp1))))/sum(length((imp1)))))

imp2 = imp2[,!(names(imp2) %in% drop)] #10
sapply(imp2, function(imp2) sum(length(which(is.na(imp2)))))
sapply(imp2, function(imp2) (100*sum(length(which(is.na(imp2))))/sum(length((imp2)))))

imp3 = imp3[,!(names(imp3) %in% drop)] #10
sapply(imp3, function(imp3) sum(length(which(is.na(imp3)))))
sapply(imp3, function(imp3) (100*sum(length(which(is.na(imp3))))/sum(length((imp3)))))
```

Imputacion con randoforest:
```{r}
colnames(imp1)

set.seed(1)

imp1 = missForest(imp1, verbose = T)
imp2 = missForest(imp2, verbose = T)
imp3 = missForest(imp3, verbose = T)

imp1$OOBerror

imp1 = data.frame(imp1$ximp)
imp2 = data.frame(imp2$ximp)
imp3 = data.frame(imp3$ximp)

colnames(imp1)
```

Normalización Ranknorm.
```{r}
#imp1
imp1.rankNorm = apply(imp1, 2, rankNorm)

MO.rankNorm = data.frame(cbind(MO[1], imp1.rankNorm))

save(MO.rankNorm, file = "MO.rN.22102020.rda")

rm(imp1, imp1.rankNorm, na, na1, na2, drop)

#imp2
imp2.rankNorm = apply(imp2, 2, rankNorm)

MO1.rankNorm = data.frame(cbind(MO1[1], imp2.rankNorm))

save(MO1.rankNorm, file = "MO1.rN.22102020.rda")

rm(imp2, imp2.rankNorm, na, na1, na2, drop)

#imp3
imp3.rankNorm = apply(imp3, 2, rankNorm)

MO2.rankNorm = data.frame(cbind(MO2[1], imp3.rankNorm))

save(MO2.rankNorm, file = "MO2.rN.22102020.rda")

rm(imp3, imp3.rankNorm, na, na1, na2, drop)
```

Split t-v
```{r}
set.seed(1002)

train.data = c()
test.data = c()
cv = c()
bT = c()
coef.met = c()
predictions = c()
metricas = c()

for (i in 1:10) {
  training = MO1.rankNorm$v2_bone_mass %>% createDataPartition(p = 0.8, list = F)
  
  train.data[[i]] = MO1.rankNorm[training,]
  test.data[[i]] = MO1.rankNorm[-training,]

  cv[[i]] = train(v2_bone_mass ~ . , data = train.data[[i]],  method = "glmnet", trControl = trainControl("cv", number = 10), tuneLength = 10)
  
  bT[[i]] = cv[[i]]$bestTune

  x.test = model.matrix(v2_bone_mass~., test.data[[i]])[,-1] #Predictores
  y.test = test.data[[i]]$v2_bone_mass #Outcome
  
  predictions[[i]] = cv[[i]] %>% predict(x.test)
  
  metricas[[i]] = data.frame(RMSE = RMSE(predictions[[i]], test.data[[i]]$v2_bone_mass), 
                             Rsquare = R2(predictions[[i]], test.data[[i]]$v2_bone_mass),
                             Pearson.IC = cor.test(predictions[[i]], test.data[[i]]$v2_bone_mass, conf.level = 0.95)$conf.int,
                             Pearson.cor = cor.test(predictions[[i]], test.data[[i]]$v2_bone_mass)$estimate)
}

l = 0.05976 #cv[[i]]$bestTune$lambda
alp = 0.6 #cv[[i]]$bestTune$alpha #1

train_rows = c()
training = c()
validation = c()
modelito = c()
X = c()
Y = c()
X1 = c()
Y1 = c()
pred = c()
met = c()

for (i in 1:10) {

  train_rows[[i]] = MO1.rankNorm$v2_bone_mass %>% createDataPartition(p = 0.80, list = F)
  training[[i]] = MO1.rankNorm[train_rows[[i]],]
  validation[[i]] = MO1.rankNorm[-train_rows[[i]],]
  
  a=as.data.frame(training[[i]])
  X[[i]] = as.matrix(a[,2:(dim(a)[2]-1)])
  Y[[i]] = a$v2_bone_mass
  
  b=as.data.frame(validation[[i]])
  X1[[i]] = as.matrix(b[,2:(dim(a)[2]-1)])
  Y1[[i]] = b$v2_bone_mass

  modelito[[i]] = glmnet(X[[i]], Y[[i]], alpha = alp, lambda = l)

  pred[[i]] = modelito[[i]] %>% predict(X1[[i]])
  
  met[[i]] = data.frame(RMSE = RMSE(pred[[i]], Y1[[i]]), Rcuadrado = R2(pred[[i]], Y1[[i]]), 
                        Pearson.IC = cor.test(pred[[i]], Y1[[i]], conf.level = 0.95)$conf.int,
                        Pearson.cor = cor.test(pred[[i]], Y1[[i]])$estimate)
}

coef_met_min = c()

for (i in 1:10) {
  
  metabolitos.min = coef(modelito[[i]], s = l)
  
  coef.min = metabolitos.min@x
  coef.min = coef.min[-1]
  
  seleccion.min = c(metabolitos.min@i)
  nombres = as.data.frame(colnames(X[[i]]))
  seleccion.min1 = nombres[seleccion.min,]
  seleccion.min1 = as.character(seleccion.min1)
  coef_met_min[[i]] = cbind(seleccion.min1, coef.min)
}

(coef.min.ajust = print(coef_met_min))

BBDD.met.min = data.frame(colnames(MO1.rankNorm))
colnames(BBDD.met.min)[colnames(BBDD.met.min)=="colnames.MO1.rankNorm."] = "seleccion.min1"

for (i in 1:10) {
  BBDD.met.min= merge(BBDD.met.min, coef_met_min[[i]], by = "seleccion.min1", all.x = T)
}
colnames(BBDD.met.min) = c("metabolitos","v1","v2","v3","v4","v5","v6","v7","v8","v9","v10")

export(BBDD.met.min, "06102020.ValMO1v0.all.xlsx")
```

Internal Validation 

```{r}
model.met = data.frame(met$metabolitos, met$media)
str(model.met)
model.met$met.metabolitos = as.character(model.met$met.metabolitos)
model.met = subset(model.met, met.media != 0)

X2 = data.frame(MO1.rankNorm[,2:(dim(MO1.rankNorm)[2])])
Y2 = MO1.rankNorm$v2_bone_mass
b = model.met

if (length(b) != 0) {
  stdout  = vector('character')
  con = textConnection('stdout', 'wr', local = TRUE)
  sink(con)
  for (i in 1:(dim(b)[1]))
  {
    if (i == 1) {cat("X2$model=", sep="")}
    if ((i < (dim(b)[1])) & (i > 0)) {cat("(","X2$",b[[i,1]],"*", b[[i,2]],")"," +", "\n", sep="")}
    if (i == (dim(b)[1])) {cat("(","X2$",b[[i,1]],"*", b[[i,2]],")", sep="")}
  }
  sink()
  
  close(con)
  
  stdout_full = paste(unlist(stdout), collapse =" ")
  stdout_full[1]
  
  eval(parse(text=stdout_full[1]))
  
  modelo = X2$model
  
  correr = cor.test(Y2,modelo, conf.level = .95)
} 
```

Coef con la BBDD completa y cv:
```{r}
set.seed(1000)

cvfit = c()
coef_met_min = c()

X = data.frame(MO1.rankNorm[,2:(dim(MO1.rankNorm)[2])])
Y = MO1.rankNorm$v2_bone_mass

for (i in 1:10) {
  
  cvfit[[i]] = cv.glmnet(as.matrix(X), as.matrix(Y), family="gaussian", type.measure = "mse", alpha = alp)
  
  metabolitos.min = coef(cvfit[[i]], s = "lambda.min")
  
  coef.min = metabolitos.min@x
  coef.min = coef.min[-1]
  
  seleccion.min = c(metabolitos.min@i)
  nombres = as.data.frame(colnames(X))
  seleccion.min1 = nombres[seleccion.min,]
  seleccion.min1 = as.character(seleccion.min1)
  coef_met_min[[i]] = cbind(seleccion.min1, coef.min)
}

BBDD.met.min = data.frame(colnames(MO1.rankNorm))
colnames(BBDD.met.min)[colnames(BBDD.met.min)=="colnames.MO1.rankNorm."] = "seleccion.min1"

for (i in 1:10) {
  BBDD.met.min= merge(BBDD.met.min, coef_met_min[[i]], by = "seleccion.min1", all.x = T)
}
colnames(BBDD.met.min) = c("metabolitos","v1","v2","v3","v4","v5","v6","v7","v8","v9","v10")

export(BBDD.met.min, "06102020.CoefMO1v0.all.m50.xlsx")

save.image(file = "06102020.CoefMO1v0.all.m50.RData")
```

Comparacion DEXA:
```{r}
hist(BBDD3$v2_bone_mass)
hist(BBDD3$v8_bone_mass)

wilcox.test(BBDD3$v8_bone_mass, BBDD3$v2_bone_mass, paired = T,conf.int = 0.95)

mean(BBDD3$v2_bone_mass, na.rm = T)
mean(BBDD3$v8_bone_mass, na.rm = T)

BBDD4 = subset(BBDD3, v8_bone_mass >= -100)

mean(BBDD4$v2_bone_mass, na.rm = T)
mean(BBDD4$v8_bone_mass, na.rm = T)

hist(BBDD3$v14_bone_mass)

wilcox.test(BBDD3$v14_bone_mass, BBDD3$v2_bone_mass, paired = T,conf.int = 0.95)

mean(BBDD3$v2_bone_mass, na.rm = T)
mean(BBDD3$v14_bone_mass, na.rm = T)

BBDD3$delta = BBDD3$v8_bone_mass - BBDD3$v2_bone_mass
mean(BBDD1$delta, na.rm = T)
median(BBDD1$delta, na.rm = T)
sd(BBDD1$delta, na.rm = T)

BBDD3$delta = BBDD3$v8_bone_mass - BBDD3$v2_bone_mass
mean(BBDD3$delta, na.rm = T)
median(BBDD3$delta, na.rm = T)
sd(BBDD3$delta, na.rm = T)
```

##2. In less than 50 years old women, to analyse whether bone mass significantly changes during the 8-week LCD. Then, to predict bone mass change by metabolic signatures at baseline. The metabolic signature will be calculated as the weighted sum of the selected metabolites in Aim 1 with weights equal to coefficients from the elastic net regression. 10-fold cross-validation in order to avoid overfitting.

Prediccion a partir de los valores obtenidos en V0, de los valores que se pueden obtener en v8 y v14, es decir, validacion semi-externa del modelo generado

Se cargan los modelos en V0:
```{r}
modelos = read_excel("modelos.xlsx", col_types = c("text", "text", "numeric"))
mo.MO.v8 = subset(modelos, MO.v0 > -10, select = c("metabolitos...1","MO.v0"))
mo.MO.v14 = subset(modelos, MO.v0 > -10, select = c("metabolitos...2","MO.v0"))
```

BBDD de validación:
V8:
```{r}
bq3 = bq[,c("allocation_no","v8_FChr","v8_Echr",        
            "v8_TChr","v8_TG_RMN","v8_PhCh","v8_LCh","v8_Sph","v8_FAC","v8_w3","v8_ARA+EPA",
            "v8_DHA","v8_LINOLEIC","v8_PUFA","v8_MUFA","v8_LPC_14.0","v8_LPC_15.0","v8_LPC_16.0","v8_LPC_16.0.e", 
            "v8_LPC_16.1","v8_LPC_16.1.e","v8_LPC_17.0","v8_LPC_18.0","v8_LPC_18.0.e","v8_LPC_18.1","v8_LPC_18.2","v8_LPC_20.0",    
            "v8_LPC_20.1","v8_LPC_20.3","v8_LPC_20.4","v8_LPC_22.6","v8_PC_30.0","v8_PC_32.0","v8_PC_32.1","v8_PC_32.1.e",   
            "v8_PC_32.2","v8_PC_33.1","v8_PC_34.0","v8_PC_34.1","v8_PC_34.1.e","v8_PC_34.2","v8_PC_34.2.e","v8_PC_34.3.e",   
            "v8_PC_34.4","v8_PC_35.1","v8_PC_35.2","v8_PC_36.1","v8_PC_36.2","v8_PC_36.2.e","v8_PC_36.3","v8_PC_36.4",  
            "v8_PC_36.4.e","v8_PC_36.5","v8_PC_36.5.e","v8_PC_37.4","v8_PC_38.3","v8_PC_38.4","v8_PC_38.4.e","v8_PC_38.5",    
            "v8_PC_38.5.e","v8_PC_38.6","v8_PC_40.4","v8_PC_40.4.e","v8_PC_40.5.e","v8_PC_40.6","v8_PC_42.5.e","v8_PE_36.5.e",   
            "v8_PE_38.5.e","v8_PE_38.6.e","v8_SM_32.1","v8_SM_32.2","v8_SM_33.1","v8_SM_34.1","v8_SM_34.2","v8_SM_35.1",    
            "v8_SM_36.0","v8_SM_36.1","v8_SM_36.2","v8_SM_38.1","v8_SM_38.2","v8_SM_40.1","v8_SM_40.2","v8_SM_41.1",    
            "v8_SM_41.2","v8_SM_42.1","v8_SM_42.2","v8_SM_42.3","v8_TG_48.1","v8_TG_48.2","v8_TG_50.1","v8_TG_50.2",    
            "v8_TG_50.3","v8_TG_50.4","v8_TG_51.2","v8_TG_52.1","v8_TG_52.2","v8_TG_52.3","v8_TG_52.4","v8_TG_52.5",    
            "v8_TG_54.2","v8_TG_54.3","v8_TG_54.4","v8_TG_54.5","v8_LA","v8_GlyA","v8_Ala","v8_Glycine",
            "v8_2-HbutA","v8_3-HbutA","v8_Val","v8_Leu","v8_Glycerol","v8_Isoleu","v8_Proline","v8_GliA",        
            "v8_Ser","v8_Thr","v8_Meth","v8_Orn","v8_GlutA","v8_Phe","v8_Lys","v8_CitA",        
            "v8_Fruc","v8_Glu","v8_Tyr","v8_PalA","v8_LinoA","v8_OleicA","v8_StearicA","v8_Tryp", 
            "v8_Sucrose","v8_alphaToco","v8_Choles")]

BBDD1 = BBDD[,c("allocation_no","v2_bone_mass","v8_bone_mass","v14_bone_mass","v1_gq_11", "v1_gq_1")]
BBDD5 = merge(BBDD1, bq3, by = "allocation_no")
BBDD5 = subset(BBDD5, (v1_gq_11 <= 50 & v1_gq_1 == 1) | v1_gq_1 == 2)
```

Modelos continuos:
```{r}
MO1 = cbind(BBDD5[3], BBDD5[7:139])

MO1 = subset(MO1, v8_bone_mass > 0) #Se eliminan valores NA. 2

#NAs en cada variables
imp2 = MO1[2:134]

na1 = sapply(imp2, function(imp2) sum(length(which(is.na(imp2)))))
na2 = sapply(imp2, function(imp2) (100*sum(length(which(is.na(imp2))))/sum(length((imp2)))))
na =  cbind("NA"=na1, "% NA"=na2)
na

#ReMGveMGs > 20%
drop = A_NA
imp2 = imp2[,!(names(imp2) %in% drop)] #10
sapply(imp2, function(imp2) sum(length(which(is.na(imp2)))))
sapply(imp2, function(imp2) (100*sum(length(which(is.na(imp2))))/sum(length((imp2)))))
```

Imputacion con randoforest:
```{r}
set.seed(1)

imp2 = missForest(imp2, verbose = T)

imp2 = data.frame(imp2$ximp)
```

Normalización Ranknorm.
```{r}
imp2.rankNorm = apply(imp2, 2, rankNorm)

MO1.rankNorm = data.frame(cbind(MO1[1], imp2.rankNorm))

save(MO1.rankNorm, file = "MO1.rN.v8.06102020.rda")

rm(imp2, imp2.rankNorm, na, na1, na2, drop)
```

Validación V8:

```{r}
X2 = data.frame(MO1.rankNorm[,2:(dim(MO1.rankNorm)[2])])
Y2 = MO1.rankNorm$v8_bone_mass
b = mo.MO.v8

if (length(b) != 0) {
  stdout  = vector('character')
  con = textConnection('stdout', 'wr', local = TRUE)
  sink(con)
  for (i in 1:(dim(b)[1]))
  {
    if (i == 1) {cat("X2$model=", sep="")}
    if ((i < (dim(b)[1])) & (i > 0)) {cat("(","X2$",b[[i,1]],"*", b[[i,2]],")"," +", "\n", sep="")}
    if (i == (dim(b)[1])) {cat("(","X2$",b[[i,1]],"*", b[[i,2]],")", sep="")}
  }
  sink()
  
  close(con)
  
  stdout_full = paste(unlist(stdout), collapse =" ")
  stdout_full[1]
  
  eval(parse(text=stdout_full[1]))
  
  modelo = X2$model
  
  correr = cor.test(Y2,modelo, conf.level = .95)
  
} 
```

V14:
```{r}
bq3 = bq[,c("allocation_no","v14_FChr","v14_Echr",        
            "v14_TChr","v14_TG_RMN","v14_PhCh","v14_LCh","v14_Sph","v14_FAC","v14_w3","v14_ARA+EPA",
            "v14_DHA","v14_LINOLEIC","v14_PUFA","v14_MUFA","v14_LPC_14.0","v14_LPC_15.0","v14_LPC_16.0","v14_LPC_16.0.e", 
            "v14_LPC_16.1","v14_LPC_16.1.e","v14_LPC_17.0","v14_LPC_18.0","v14_LPC_18.0.e",
            "v14_LPC_18.1","v14_LPC_18.2","v14_LPC_20.0",    
            "v14_LPC_20.1","v14_LPC_20.3","v14_LPC_20.4","v14_LPC_22.6","v14_PC_30.0",
            "v14_PC_32.0","v14_PC_32.1","v14_PC_32.1.e",   
            "v14_PC_32.2","v14_PC_33.1","v14_PC_34.0","v14_PC_34.1","v14_PC_34.1.e",
            "v14_PC_34.2","v14_PC_34.2.e","v14_PC_34.3.e",   
            "v14_PC_34.4","v14_PC_35.1","v14_PC_35.2","v14_PC_36.1","v14_PC_36.2","v14_PC_36.2.e","v14_PC_36.3","v14_PC_36.4",  
            "v14_PC_36.4.e","v14_PC_36.5","v14_PC_36.5.e","v14_PC_37.4","v14_PC_38.3",
            "v14_PC_38.4","v14_PC_38.4.e","v14_PC_38.5",    
            "v14_PC_38.5.e","v14_PC_38.6","v14_PC_40.4","v14_PC_40.4.e","v14_PC_40.5.e",
            "v14_PC_40.6","v14_PC_42.5.e","v14_PE_36.5.e",   
            "v14_PE_38.5.e","v14_PE_38.6.e","v14_SM_32.1","v14_SM_32.2","v14_SM_33.1","v14_SM_34.1","v14_SM_34.2","v14_SM_35.1",    
            "v14_SM_36.0","v14_SM_36.1","v14_SM_36.2","v14_SM_38.1","v14_SM_38.2","v14_SM_40.1","v14_SM_40.2","v14_SM_41.1",    
            "v14_SM_41.2","v14_SM_42.1","v14_SM_42.2","v14_SM_42.3","v14_TG_48.1","v14_TG_48.2","v14_TG_50.1","v14_TG_50.2",    
            "v14_TG_50.3","v14_TG_50.4","v14_TG_51.2","v14_TG_52.1","v14_TG_52.2","v14_TG_52.3","v14_TG_52.4","v14_TG_52.5",    
            "v14_TG_54.2","v14_TG_54.3","v14_TG_54.4","v14_TG_54.5","v14_LA","v14_GlyA","v14_Ala","v14_Glycine",
            "v14_2-HbutA","v14_3-HbutA","v14_Val","v14_Leu","v14_Glycerol","v14_Isoleu","v14_Proline","v14_GliA",        
            "v14_Ser","v14_Thr","v14_Meth","v14_Orn","v14_GlutA","v14_Phe","v14_Lys","v14_CitA",        
            "v14_Fruc","v14_Glu","v14_Tyr","v14_PalA","v14_LinoA","v14_OleicA","v14_StearicA","v14_Tryp", 
            "v14_Sucrose","v14_alphaToco","v14_Choles")]

BBDD5 = merge(BBDD1, bq3, by = "allocation_no")
BBDD5 = subset(BBDD5, (v1_gq_11 <= 50 & v1_gq_1 == 1) | v1_gq_1 == 2)
```

Modelos continuos:
```{r}
MO1 = cbind(BBDD5[4], BBDD5[7:139])

MO1 = subset(MO1, v14_bone_mass > 0) #Se eliminan valores NA. 2

#NAs en cada variables
imp2 = MO1[2:134]

na1 = sapply(imp2, function(imp2) sum(length(which(is.na(imp2)))))
na2 = sapply(imp2, function(imp2) (100*sum(length(which(is.na(imp2))))/sum(length((imp2)))))
na =  cbind("NA"=na1, "% NA"=na2)
na

#ReMGveMGs > 20%
drop = A_NA
imp2 = imp2[,!(names(imp2) %in% drop)] #10
sapply(imp2, function(imp2) sum(length(which(is.na(imp2)))))
sapply(imp2, function(imp2) (100*sum(length(which(is.na(imp2))))/sum(length((imp2)))))
```

Imputacion con randoforest:
```{r}
set.seed(1)

imp2 = missForest(imp2, verbose = T)

imp2 = data.frame(imp2$ximp)
```

Normalización Ranknorm.
```{r}
imp2.rankNorm = apply(imp2, 2, rankNorm)

MO1.rankNorm = data.frame(cbind(MO1[1], imp2.rankNorm))

save(MO1.rankNorm, file = "MO1.rN.v14.06102020.rda")

rm(imp2, imp2.rankNorm, na, na1, na2, drop)
```

Validación V14:
```{r}
X2 = data.frame(MO1.rankNorm[,2:(dim(MO1.rankNorm)[2])])
Y2 = MO1.rankNorm$v14_bone_mass
b = mo.MO.v14

if (length(b) != 0) {
  stdout  = vector('character')
  con = textConnection('stdout', 'wr', local = TRUE)
  sink(con)
  for (i in 1:(dim(b)[1]))
  {
    if (i == 1) {cat("X2$model=", sep="")}
    if ((i < (dim(b)[1])) & (i > 0)) {cat("(","X2$",b[[i,1]],"*", b[[i,2]],")"," +", "\n", sep="")}
    if (i == (dim(b)[1])) {cat("(","X2$",b[[i,1]],"*", b[[i,2]],")", sep="")}
  }
  sink()
  
  close(con)
  
  stdout_full = paste(unlist(stdout), collapse =" ")
  stdout_full[1]
  
  eval(parse(text=stdout_full[1]))
  
  modelo = X2$model
  
  correr = cor.test(Y2,modelo, conf.level = .95)
} 
```

##3. In less than 50 years old women, to examine the relationship between changes in circulating metabolites during the LCD and changes in bone mass. In this case we need to adjust for confounders such as age, sex, weight changes, baseline bone mass, baseline value of each metabolite. Multivariable linear regression models with bonferroni correction.

```{r}
BD1 = BBDD[,c("allocation_no","v2_bone_mass","v8_bone_mass","v1_gq_11", "v1_gq_1", "v2_4.3", "v2_4.4", "v8_4.1", "v8_4.2")]

BD1$delta_bone = BD1$v8_bone_mass - BD1$v2_bone_mass
BD1$delta_peso = ((BD1$v8_4.1+BD1$v8_4.2)/2) - ((BD1$v2_4.3+BD1$v2_4.4)/2)

BD2 = subset(BD1, (v1_gq_11 <= 50 & v1_gq_1 == 1) | v1_gq_1 == 2) #Se eliminan mujeres mayores de 50 --> 74

BD2 = subset(BD2, delta_bone > -100) #Se eliminan valores NA. 44

bq4 = bq[,c("allocation_no","v2_FChr","v2_Echr",        
            "v2_TChr","v2_TG_RMN","v2_PhCh","v2_LCh","v2_Sph","v2_FAC","v2_w3","v2_ARA+EPA",
            "v2_DHA","v2_LINOLEIC","v2_PUFA","v2_MUFA","v2_LPC_14.0","v2_LPC_15.0","v2_LPC_16.0","v2_LPC_16.0.e", 
            "v2_LPC_16.1","v2_LPC_16.1.e","v2_LPC_17.0","v2_LPC_18.0","v2_LPC_18.0.e","v2_LPC_18.1","v2_LPC_18.2","v2_LPC_20.0",    
            "v2_LPC_20.1","v2_LPC_20.3","v2_LPC_20.4","v2_LPC_22.6","v2_PC_30.0","v2_PC_32.0","v2_PC_32.1","v2_PC_32.1.e",   
            "v2_PC_32.2","v2_PC_33.1","v2_PC_34.0","v2_PC_34.1","v2_PC_34.1.e","v2_PC_34.2","v2_PC_34.2.e","v2_PC_34.3.e",   
            "v2_PC_34.4","v2_PC_35.1","v2_PC_35.2","v2_PC_36.1","v2_PC_36.2","v2_PC_36.2.e","v2_PC_36.3","v2_PC_36.4",  
            "v2_PC_36.4.e","v2_PC_36.5","v2_PC_36.5.e","v2_PC_37.4","v2_PC_38.3","v2_PC_38.4","v2_PC_38.4.e","v2_PC_38.5",    
            "v2_PC_38.5.e","v2_PC_38.6","v2_PC_40.4","v2_PC_40.4.e","v2_PC_40.5.e","v2_PC_40.6","v2_PC_42.5.e","v2_PE_36.5.e",   
            "v2_PE_38.5.e","v2_PE_38.6.e","v2_SM_32.1","v2_SM_32.2","v2_SM_33.1","v2_SM_34.1","v2_SM_34.2","v2_SM_35.1",    
            "v2_SM_36.0","v2_SM_36.1","v2_SM_36.2","v2_SM_38.1","v2_SM_38.2","v2_SM_40.1","v2_SM_40.2","v2_SM_41.1",    
            "v2_SM_41.2","v2_SM_42.1","v2_SM_42.2","v2_SM_42.3","v2_TG_48.1","v2_TG_48.2","v2_TG_50.1","v2_TG_50.2",    
            "v2_TG_50.3","v2_TG_50.4","v2_TG_51.2","v2_TG_52.1","v2_TG_52.2","v2_TG_52.3","v2_TG_52.4","v2_TG_52.5",    
            "v2_TG_54.2","v2_TG_54.3","v2_TG_54.4","v2_TG_54.5","v2_LA","v2_GlyA","v2_Ala","v2_Glycine",
            "v2_2-HbutA","v2_3-HbutA","v2_Val","v2_Leu","v2_Glycerol","v2_Isoleu","v2_Proline","v2_GliA",        
            "v2_Ser","v2_Thr","v2_Meth","v2_Orn","v2_GlutA","v2_Phe","v2_Lys","v2_CitA",        
            "v2_Fruc","v2_Glu","v2_Tyr","v2_PalA","v2_LinoA","v2_OleicA","v2_StearicA","v2_Tryp", 
            "v2_Sucrose","v2_alphaToco","v2_Choles")]
```

Normalizacion de los datos:
```{r}
#NAs en cada variables
imp2 = bq4[2:134]

na1 = sapply(imp2, function(imp2) sum(length(which(is.na(imp2)))))
na2 = sapply(imp2, function(imp2) (100*sum(length(which(is.na(imp2))))/sum(length((imp2)))))
na =  cbind("NA"=na1, "% NA"=na2)
na

#ReMGveMGs > 20%
drop = c("v2_PC_34.1","v2_PC_36.2","v2_PC_36.3","v2_PC_36.4","v2_TG_48.1",
         "v2_TG_48.2","v2_TG_50.4","v2_TG_51.2","v2_TG_52.1","v2_TG_52.5")

#ReMGveMGs > 20%
imp2 = imp2[,!(names(imp2) %in% drop)] #10 (marcados previamente)
sapply(imp2, function(imp2) sum(length(which(is.na(imp2)))))
sapply(imp2, function(imp2) (100*sum(length(which(is.na(imp2))))/sum(length((imp2)))))
```

Imputacion con randoforest:
```{r}
colnames(imp2)

imp2 = data.frame(imp2)

set.seed(1)

imp2 = missForest(imp2, verbose = T)

imp2 = data.frame(imp2$ximp)

colnames(imp2)
```

Normalización Ranknorm.
```{r}
imp2.rankNorm = apply(imp2, 2, rankNorm)

bq4 = data.frame(cbind(bq4[1], imp2.rankNorm))


bq5 = bq[,c("allocation_no",
            "v8_FChr","v8_Echr","v8_TChr","v8_TG_RMN","v8_PhCh","v8_LCh","v8_Sph","v8_FAC","v8_w3","v8_ARA+EPA",
            "v8_DHA","v8_LINOLEIC","v8_PUFA","v8_MUFA","v8_LPC_14.0","v8_LPC_15.0","v8_LPC_16.0","v8_LPC_16.0.e", 
            "v8_LPC_16.1","v8_LPC_16.1.e","v8_LPC_17.0","v8_LPC_18.0","v8_LPC_18.0.e","v8_LPC_18.1","v8_LPC_18.2","v8_LPC_20.0",    
            "v8_LPC_20.1","v8_LPC_20.3","v8_LPC_20.4","v8_LPC_22.6","v8_PC_30.0","v8_PC_32.0","v8_PC_32.1","v8_PC_32.1.e",   
            "v8_PC_32.2","v8_PC_33.1","v8_PC_34.0","v8_PC_34.1","v8_PC_34.1.e","v8_PC_34.2","v8_PC_34.2.e","v8_PC_34.3.e",   
            "v8_PC_34.4","v8_PC_35.1","v8_PC_35.2","v8_PC_36.1","v8_PC_36.2","v8_PC_36.2.e","v8_PC_36.3","v8_PC_36.4",  
            "v8_PC_36.4.e","v8_PC_36.5","v8_PC_36.5.e","v8_PC_37.4","v8_PC_38.3","v8_PC_38.4","v8_PC_38.4.e","v8_PC_38.5",    
            "v8_PC_38.5.e","v8_PC_38.6","v8_PC_40.4","v8_PC_40.4.e","v8_PC_40.5.e","v8_PC_40.6","v8_PC_42.5.e","v8_PE_36.5.e",   
            "v8_PE_38.5.e","v8_PE_38.6.e","v8_SM_32.1","v8_SM_32.2","v8_SM_33.1","v8_SM_34.1","v8_SM_34.2","v8_SM_35.1",    
            "v8_SM_36.0","v8_SM_36.1","v8_SM_36.2","v8_SM_38.1","v8_SM_38.2","v8_SM_40.1","v8_SM_40.2","v8_SM_41.1",    
            "v8_SM_41.2","v8_SM_42.1","v8_SM_42.2","v8_SM_42.3","v8_TG_48.1","v8_TG_48.2","v8_TG_50.1","v8_TG_50.2",    
            "v8_TG_50.3","v8_TG_50.4","v8_TG_51.2","v8_TG_52.1","v8_TG_52.2","v8_TG_52.3","v8_TG_52.4","v8_TG_52.5",    
            "v8_TG_54.2","v8_TG_54.3","v8_TG_54.4","v8_TG_54.5","v8_LA","v8_GlyA","v8_Ala","v8_Glycine",
            "v8_2-HbutA","v8_3-HbutA","v8_Val","v8_Leu","v8_Glycerol","v8_Isoleu","v8_Proline","v8_GliA",        
            "v8_Ser","v8_Thr","v8_Meth","v8_Orn","v8_GlutA","v8_Phe","v8_Lys","v8_CitA",        
            "v8_Fruc","v8_Glu","v8_Tyr","v8_PalA","v8_LinoA","v8_OleicA","v8_StearicA","v8_Tryp", 
            "v8_Sucrose","v8_alphaToco","v8_Choles")]
```

Normalizacion de los datos:
```{r}
#NAs en cada variables
imp2 = bq5[2:134]

na1 = sapply(imp2, function(imp2) sum(length(which(is.na(imp2)))))
na2 = sapply(imp2, function(imp2) (100*sum(length(which(is.na(imp2))))/sum(length((imp2)))))
na =  cbind("NA"=na1, "% NA"=na2)
na

#ReMGveMGs > 20%
drop = c("v8_PC_34.1","v8_PC_36.2","v8_PC_36.3","v8_PC_36.4","v8_TG_48.1",
         "v8_TG_48.2","v8_TG_50.4","v8_TG_51.2","v8_TG_52.1","v8_TG_52.5")

#ReMGveMGs > 20%
imp2 = imp2[,!(names(imp2) %in% drop)] #10 (marcados previamente)
sapply(imp2, function(imp2) sum(length(which(is.na(imp2)))))
sapply(imp2, function(imp2) (100*sum(length(which(is.na(imp2))))/sum(length((imp2)))))
```

Imputacion con randoforest:
```{r}
colnames(imp2)

imp2 = data.frame(imp2)

set.seed(1)

imp2 = missForest(imp2, verbose = T)

imp2 = data.frame(imp2$ximp)

colnames(imp2)
```

Normalización Ranknorm.
```{r}
imp2.rankNorm = apply(imp2, 2, rankNorm)

bq5 = data.frame(cbind(bq5[1], imp2.rankNorm))

bq6 = bq5[2:124] - bq4[2:124]

colnames(bq6)

colnames(bq6) = c("FChr","Echr","TChr","TG_RMN","PhCh","LCh","Sph","FAC","w3",
                  "ARA.EPA","DHA","LINOLEIC","PUFA","MUFA","LPC_14.0","LPC_15.0","LPC_16.0","LPC_16.0.e",
                  "LPC_16.1","LPC_16.1.e","LPC_17.0","LPC_18.0","LPC_18.0.e","LPC_18.1","LPC_18.2","LPC_20.0","LPC_20.1",
                  "LPC_20.3","LPC_20.4","LPC_22.6","PC_30.0","PC_32.0","PC_32.1","PC_32.1.e","PC_32.2","PC_33.1",
                  "PC_34.0","PC_34.1.e","PC_34.2","PC_34.2.e","PC_34.3.e","PC_34.4","PC_35.1","PC_35.2","PC_36.1",
                  "PC_36.2.e","PC_36.4.e","PC_36.5","PC_36.5.e","PC_37.4","PC_38.3","PC_38.4","PC_38.4.e","PC_38.5",
                  "PC_38.5.e","PC_38.6","PC_40.4","PC_40.4.e","PC_40.5.e","PC_40.6","PC_42.5.e","PE_36.5.e","PE_38.5.e",
                  "PE_38.6.e","SM_32.1","SM_32.2","SM_33.1","SM_34.1","SM_34.2","SM_35.1","SM_36.0","SM_36.1",
                  "SM_36.2","SM_38.1","SM_38.2","SM_40.1","SM_40.2","SM_41.1","SM_41.2","SM_42.1","SM_42.2",
                  "SM_42.3","TG_50.1","TG_50.2","TG_50.3","TG_52.2","TG_52.3","TG_52.4","TG_54.2","TG_54.3",
                  "TG_54.4","TG_54.5","LA","GlyA","Ala","Glycine","2.HbutA","3.HbutA","Val",
                  "Leu","Glycerol","Isoleu","Proline","GliA","Ser","Thr","Meth","Orn",
                  "GlutA","Phe","Lys","CitA","Fruc","Glu","Tyr","PalA","LinoA",
                  "OleicA","StearicA","Tryp","Sucrose","alphaToco","Choles")

bq6$allocation_no = bq5$allocation_no

bq7 = merge(bq4, bq6, by = "allocation_no")

MO.deltas.rN = merge(BD2, bq7, by = "allocation_no")

colnames(MO.deltas.rN)

MO.deltas.rN1 = data.frame(cbind(MO.deltas.rN$delta_bone, MO.deltas.rN[134:256]))

save(MO.deltas.rN1, file = "MO.deltas.rN1.22102020.rda")
```

Regresiones lineales
```{r}
resumen = c()
IC = c()

for (i in 135:257) {
  resumen[[i]] = summary(lm((MO.deltas.rN[[10]]) ~ as.numeric(MO.deltas.rN[[i]]) + (MO.deltas.rN[[4]]) + (MO.deltas.rN[[5]]) + (MO.deltas.rN[[11]])
                           + (MO.deltas.rN[[2]]) + (MO.deltas.rN[[i - 123]])))[[4]][[23]]
}

resumen.df = data.frame(matrix(unlist(resumen), ncol=1, byrow=T),stringsAsFactors=FALSE)

colnames(resumen.df) = c("P-value")
resumen.df$names = colnames(MO.deltas.rN)

export(resumen.df, "res4pvalues.xlsx")
```

No sale nada con la regresion lineal y transformacion rankNorm. Procedo a repetir el proceso con log y 1SD.

Normalizacion de los datos:
```{r}
bq4 = bq[,c("allocation_no","v2_FChr","v2_Echr",        
            "v2_TChr","v2_TG_RMN","v2_PhCh","v2_LCh","v2_Sph","v2_FAC","v2_w3","v2_ARA+EPA",
            "v2_DHA","v2_LINOLEIC","v2_PUFA","v2_MUFA","v2_LPC_14.0","v2_LPC_15.0","v2_LPC_16.0","v2_LPC_16.0.e", 
            "v2_LPC_16.1","v2_LPC_16.1.e","v2_LPC_17.0","v2_LPC_18.0","v2_LPC_18.0.e","v2_LPC_18.1","v2_LPC_18.2","v2_LPC_20.0",    
            "v2_LPC_20.1","v2_LPC_20.3","v2_LPC_20.4","v2_LPC_22.6","v2_PC_30.0","v2_PC_32.0","v2_PC_32.1","v2_PC_32.1.e",   
            "v2_PC_32.2","v2_PC_33.1","v2_PC_34.0","v2_PC_34.1","v2_PC_34.1.e","v2_PC_34.2","v2_PC_34.2.e","v2_PC_34.3.e",   
            "v2_PC_34.4","v2_PC_35.1","v2_PC_35.2","v2_PC_36.1","v2_PC_36.2","v2_PC_36.2.e","v2_PC_36.3","v2_PC_36.4",  
            "v2_PC_36.4.e","v2_PC_36.5","v2_PC_36.5.e","v2_PC_37.4","v2_PC_38.3","v2_PC_38.4","v2_PC_38.4.e","v2_PC_38.5",    
            "v2_PC_38.5.e","v2_PC_38.6","v2_PC_40.4","v2_PC_40.4.e","v2_PC_40.5.e","v2_PC_40.6","v2_PC_42.5.e","v2_PE_36.5.e",   
            "v2_PE_38.5.e","v2_PE_38.6.e","v2_SM_32.1","v2_SM_32.2","v2_SM_33.1","v2_SM_34.1","v2_SM_34.2","v2_SM_35.1",    
            "v2_SM_36.0","v2_SM_36.1","v2_SM_36.2","v2_SM_38.1","v2_SM_38.2","v2_SM_40.1","v2_SM_40.2","v2_SM_41.1",    
            "v2_SM_41.2","v2_SM_42.1","v2_SM_42.2","v2_SM_42.3","v2_TG_48.1","v2_TG_48.2","v2_TG_50.1","v2_TG_50.2",    
            "v2_TG_50.3","v2_TG_50.4","v2_TG_51.2","v2_TG_52.1","v2_TG_52.2","v2_TG_52.3","v2_TG_52.4","v2_TG_52.5",    
            "v2_TG_54.2","v2_TG_54.3","v2_TG_54.4","v2_TG_54.5","v2_LA","v2_GlyA","v2_Ala","v2_Glycine",
            "v2_2-HbutA","v2_3-HbutA","v2_Val","v2_Leu","v2_Glycerol","v2_Isoleu","v2_Proline","v2_GliA",        
            "v2_Ser","v2_Thr","v2_Meth","v2_Orn","v2_GlutA","v2_Phe","v2_Lys","v2_CitA",        
            "v2_Fruc","v2_Glu","v2_Tyr","v2_PalA","v2_LinoA","v2_OleicA","v2_StearicA","v2_Tryp", 
            "v2_Sucrose","v2_alphaToco","v2_Choles")]

#NAs en cada variables
imp2 = bq4[2:134]

na1 = sapply(imp2, function(imp2) sum(length(which(is.na(imp2)))))
na2 = sapply(imp2, function(imp2) (100*sum(length(which(is.na(imp2))))/sum(length((imp2)))))
na =  cbind("NA"=na1, "% NA"=na2)
na

#Metabolites con NA mayor de 20%
A_NA=rownames(na[na[,2]>20 ,])

#ReMGveMGs > 20%
drop = c("v2_PC_34.1","v2_PC_36.2","v2_PC_36.3","v2_PC_36.4","v2_TG_48.1","v2_TG_48.2","v2_TG_50.4","v2_TG_51.2","v2_TG_52.1","v2_TG_52.5")

#ReMGveMGs > 20%
imp2 = imp2[,!(names(imp2) %in% drop)] #10 (marcados previamente)
sapply(imp2, function(imp2) sum(length(which(is.na(imp2)))))
sapply(imp2, function(imp2) (100*sum(length(which(is.na(imp2))))/sum(length((imp2)))))
```

Imputacion con randoforest:
```{r}
colnames(imp2)

imp2 = data.frame(imp2)

set.seed(1)

imp2 = missForest(imp2, verbose = T)

imp2 = data.frame(imp2$ximp)

colnames(imp2)
```

Normalización log
```{r}
imp2.rankNorm = apply(imp2, 2, log10)

sd3 = apply(imp2, 2, sd)

bq4 = data.frame(cbind(bq4[1], (imp2.rankNorm/sd3)))

bq5 = bq[,c("allocation_no",
            "v8_FChr","v8_Echr","v8_TChr","v8_TG_RMN","v8_PhCh","v8_LCh","v8_Sph","v8_FAC","v8_w3","v8_ARA+EPA",
            "v8_DHA","v8_LINOLEIC","v8_PUFA","v8_MUFA","v8_LPC_14.0","v8_LPC_15.0","v8_LPC_16.0","v8_LPC_16.0.e", 
            "v8_LPC_16.1","v8_LPC_16.1.e","v8_LPC_17.0","v8_LPC_18.0","v8_LPC_18.0.e","v8_LPC_18.1","v8_LPC_18.2","v8_LPC_20.0",    
            "v8_LPC_20.1","v8_LPC_20.3","v8_LPC_20.4","v8_LPC_22.6","v8_PC_30.0","v8_PC_32.0","v8_PC_32.1","v8_PC_32.1.e",   
            "v8_PC_32.2","v8_PC_33.1","v8_PC_34.0","v8_PC_34.1","v8_PC_34.1.e","v8_PC_34.2","v8_PC_34.2.e","v8_PC_34.3.e",   
            "v8_PC_34.4","v8_PC_35.1","v8_PC_35.2","v8_PC_36.1","v8_PC_36.2","v8_PC_36.2.e","v8_PC_36.3","v8_PC_36.4",  
            "v8_PC_36.4.e","v8_PC_36.5","v8_PC_36.5.e","v8_PC_37.4","v8_PC_38.3","v8_PC_38.4","v8_PC_38.4.e","v8_PC_38.5",    
            "v8_PC_38.5.e","v8_PC_38.6","v8_PC_40.4","v8_PC_40.4.e","v8_PC_40.5.e","v8_PC_40.6","v8_PC_42.5.e","v8_PE_36.5.e",   
            "v8_PE_38.5.e","v8_PE_38.6.e","v8_SM_32.1","v8_SM_32.2","v8_SM_33.1","v8_SM_34.1","v8_SM_34.2","v8_SM_35.1",    
            "v8_SM_36.0","v8_SM_36.1","v8_SM_36.2","v8_SM_38.1","v8_SM_38.2","v8_SM_40.1","v8_SM_40.2","v8_SM_41.1",    
            "v8_SM_41.2","v8_SM_42.1","v8_SM_42.2","v8_SM_42.3","v8_TG_48.1","v8_TG_48.2","v8_TG_50.1","v8_TG_50.2",    
            "v8_TG_50.3","v8_TG_50.4","v8_TG_51.2","v8_TG_52.1","v8_TG_52.2","v8_TG_52.3","v8_TG_52.4","v8_TG_52.5",    
            "v8_TG_54.2","v8_TG_54.3","v8_TG_54.4","v8_TG_54.5","v8_LA","v8_GlyA","v8_Ala","v8_Glycine",
            "v8_2-HbutA","v8_3-HbutA","v8_Val","v8_Leu","v8_Glycerol","v8_Isoleu","v8_Proline","v8_GliA",        
            "v8_Ser","v8_Thr","v8_Meth","v8_Orn","v8_GlutA","v8_Phe","v8_Lys","v8_CitA",        
            "v8_Fruc","v8_Glu","v8_Tyr","v8_PalA","v8_LinoA","v8_OleicA","v8_StearicA","v8_Tryp", 
            "v8_Sucrose","v8_alphaToco","v8_Choles")]
```

Normalizacion de los datos:
```{r}
#NAs en cada variables
imp2 = bq5[2:134]

na1 = sapply(imp2, function(imp2) sum(length(which(is.na(imp2)))))
na2 = sapply(imp2, function(imp2) (100*sum(length(which(is.na(imp2))))/sum(length((imp2)))))
na =  cbind("NA"=na1, "% NA"=na2)
na

#Metabolites con NA mayor de 20%
A_NA=rownames(na[na[,2]>20 ,])

#ReMGveMGs > 20%
drop = c("v8_PC_34.1","v8_PC_36.2","v8_PC_36.3","v8_PC_36.4","v8_TG_48.1",
         "v8_TG_48.2","v8_TG_50.4","v8_TG_51.2","v8_TG_52.1","v8_TG_52.5")

#ReMGveMGs > 20%
imp2 = imp2[,!(names(imp2) %in% drop)] #10 (marcados previamente)
sapply(imp2, function(imp2) sum(length(which(is.na(imp2)))))
sapply(imp2, function(imp2) (100*sum(length(which(is.na(imp2))))/sum(length((imp2)))))
```

Imputacion con randoforest:
```{r}
colnames(imp2)

imp2 = data.frame(imp2)

set.seed(1)

imp2 = missForest(imp2, verbose = T)

imp2 = data.frame(imp2$ximp)

colnames(imp2)
```

Normalización log10
```{r}
imp2.rankNorm = apply(imp2, 2, log10)

sd2 = apply(imp2, 2, sd)

bq5 = data.frame(cbind(bq5[1], (imp2.rankNorm/sd2)))

bq6 = bq5[2:124] - bq4[2:124]

colnames(bq6)

colnames(bq6) = c("FChr","Echr","TChr","TG_RMN","PhCh","LCh","Sph","FAC","w3",
                  "ARA.EPA","DHA","LINOLEIC","PUFA","MUFA","LPC_14.0","LPC_15.0","LPC_16.0","LPC_16.0.e",
                  "LPC_16.1","LPC_16.1.e","LPC_17.0","LPC_18.0","LPC_18.0.e","LPC_18.1","LPC_18.2","LPC_20.0","LPC_20.1",
                  "LPC_20.3","LPC_20.4","LPC_22.6","PC_30.0","PC_32.0","PC_32.1","PC_32.1.e","PC_32.2","PC_33.1",
                  "PC_34.0","PC_34.1.e","PC_34.2","PC_34.2.e","PC_34.3.e","PC_34.4","PC_35.1","PC_35.2","PC_36.1",
                  "PC_36.2.e","PC_36.4.e","PC_36.5","PC_36.5.e","PC_37.4","PC_38.3","PC_38.4","PC_38.4.e","PC_38.5",
                  "PC_38.5.e","PC_38.6","PC_40.4","PC_40.4.e","PC_40.5.e","PC_40.6","PC_42.5.e","PE_36.5.e","PE_38.5.e",
                  "PE_38.6.e","SM_32.1","SM_32.2","SM_33.1","SM_34.1","SM_34.2","SM_35.1","SM_36.0","SM_36.1",
                  "SM_36.2","SM_38.1","SM_38.2","SM_40.1","SM_40.2","SM_41.1","SM_41.2","SM_42.1","SM_42.2",
                  "SM_42.3","TG_50.1","TG_50.2","TG_50.3","TG_52.2","TG_52.3","TG_52.4","TG_54.2","TG_54.3",
                  "TG_54.4","TG_54.5","LA","GlyA","Ala","Glycine","2.HbutA","3.HbutA","Val",
                  "Leu","Glycerol","Isoleu","Proline","GliA","Ser","Thr","Meth","Orn",
                  "GlutA","Phe","Lys","CitA","Fruc","Glu","Tyr","PalA","LinoA",
                  "OleicA","StearicA","Tryp","Sucrose","alphaToco","Choles")

bq6$allocation_no = bq5$allocation_no

bq7 = merge(bq4, bq6, by = "allocation_no")

MO.deltas.rN1 = merge(BD2, bq7, by = "allocation_no")

colnames(MO.deltas.rN1)
```

Regresiones lineales
```{r}
resumen = c()
IC = c()

for (i in 135:257) {
  resumen[[i]] = summary(lm((MO.deltas.rN1[[10]]) ~ as.numeric(MO.deltas.rN1[[i]]) + (MO.deltas.rN1[[4]]) + (MO.deltas.rN1[[5]]) + (MO.deltas.rN1[[11]])
                            + (MO.deltas.rN1[[2]]) + (MO.deltas.rN1[[i - 123]])))[[4]][[23]]
}

resumen.df = data.frame(matrix(unlist(resumen), ncol=1, byrow=T),stringsAsFactors=FALSE)

colnames(resumen.df) = c("P-value")
resumen.df$names = colnames(MO.deltas.rN1)

export(resumen.df, "res5pvalues.xlsx")
```

Ajuste de la p
```{r}
p = read_excel("res5pvalues.xlsx")

p1 = p$`P-value`

p.adj.m = p.adjust(p1, "bonferroni")
```

Prueba: comprobación cambios v2 --> v8/v14 con t-test:
```{r}
bq8 = merge(bq4, bq5, by = "allocation_no")

t.test(bq8$v8_ARA.EPA, bq8$v2_ARA.EPA, paired = T) #
t.test(bq8$v8_Glycerol, bq8$v2_Glycerol, paired = T)
t.test(bq8$v8_LinoA, bq8$v2_LinoA, paired = T)
t.test(bq8$v8_LPC_22.6, bq8$v2_LPC_22.6, paired = T) #
t.test(bq8$v8_SM_32.2, bq8$v2_SM_32.2, paired = T) #
t.test(bq8$v8_SM_40.2, bq8$v2_SM_40.2, paired = T) #
t.test(bq8$v8_Tryp, bq8$v2_Tryp, paired = T) #
```

## 4. ENR con valores basales y delta de BMD (MO2.ranknorm):

Split t-v
```{r}
set.seed(1002)

train.data = c()
test.data = c()
cv = c()
bT = c()
coef.met = c()
predictions = c()
metricas = c()

for (i in 1:10) {
  
  training = MO2.rankNorm$BBDD2.delta %>% createDataPartition(p = 0.8, list = F)
  
  train.data[[i]] = MO2.rankNorm[training,]
  test.data[[i]] = MO2.rankNorm[-training,]
  
  cv[[i]] = train(BBDD2.delta ~ . , data = train.data[[i]],  method = "glmnet", trControl = trainControl("cv", number = 10), tuneLength = 10)
  
  bT[[i]] = cv[[i]]$bestTune
  
  x.test = model.matrix(BBDD2.delta~., test.data[[i]])[,-1] #Predictores
  y.test = test.data[[i]]$BBDD2.delta #Outcome
  
  predictions[[i]] = cv[[i]] %>% predict(x.test)
  
  metricas[[i]] = data.frame(RMSE = RMSE(predictions[[i]], test.data[[i]]$BBDD2.delta), Rsquare = R2(predictions[[i]], test.data[[i]]$BBDD2.delta), Pearson.IC = cor.test(predictions[[i]], test.data[[i]]$BBDD2.delta, conf.level = 0.95)$conf.int,
                             Pearson.cor = cor.test(predictions[[i]], test.data[[i]]$BBDD2.delta)$estimate)
}

l = 0.0116801054 #cv[[i]]$bestTune$lambda
alp = 0.8 #cv[[i]]$bestTune$alpha #1

train_rows = c()
training = c()
validation = c()
modelito = c()
X = c()
Y = c()
X1 = c()
Y1 = c()
pred = c()
met = c()

for (i in 1:10) {
  
  train_rows[[i]] = MO2.rankNorm$BBDD2.delta %>% createDataPartition(p = 0.80, list = F)
  training[[i]] = MO2.rankNorm[train_rows[[i]],]
  validation[[i]] = MO2.rankNorm[-train_rows[[i]],]
  
  a=as.data.frame(training[[i]])
  X[[i]] = as.matrix(a[,2:(dim(a)[2]-1)])
  Y[[i]] = a$BBDD2.delta
  
  b=as.data.frame(validation[[i]])
  X1[[i]] = as.matrix(b[,2:(dim(a)[2]-1)])
  Y1[[i]] = b$BBDD2.delta
  
  modelito[[i]] = glmnet(X[[i]], Y[[i]], alpha = alp, lambda = l)
  
  pred[[i]] = modelito[[i]] %>% predict(X1[[i]])
  
  met[[i]] = data.frame(RMSE = RMSE(pred[[i]], Y1[[i]]), Rcuadrado = R2(pred[[i]], Y1[[i]]), Pearson.IC = cor.test(pred[[i]], Y1[[i]], conf.level = 0.95)$conf.int,
                        Pearson.cor = cor.test(pred[[i]], Y1[[i]])$estimate)
}

coef_met_min = c()

for (i in 1:10) {
  
  metabolitos.min = coef(modelito[[i]], s = l)
  
  coef.min = metabolitos.min@x
  coef.min = coef.min[-1]
  
  seleccion.min = c(metabolitos.min@i)
  nombres = as.data.frame(colnames(X[[i]]))
  seleccion.min1 = nombres[seleccion.min,]
  seleccion.min1 = as.character(seleccion.min1)
  coef_met_min[[i]] = cbind(seleccion.min1, coef.min)
}

(coef.min.ajust = print(coef_met_min))

BBDD.met.min = data.frame(colnames(MO2.rankNorm))
colnames(BBDD.met.min)[colnames(BBDD.met.min)=="colnames.MO2.rankNorm."] = "seleccion.min1"

for (i in 1:10) {
  BBDD.met.min= merge(BBDD.met.min, coef_met_min[[i]], by = "seleccion.min1", all.x = T)
}
colnames(BBDD.met.min) = c("metabolitos","v1","v2","v3","v4","v5","v6","v7","v8","v9","v10")

export(BBDD.met.min, "10102020.ValMO2v0.3.all.xlsx")
```

Internal Validation 
```{r}
met = read_excel("10102020.ValMO2v0.3.all.xlsx")

model.met = data.frame(met$metabolitos, met$media)
str(model.met)
model.met$met.metabolitos = as.character(model.met$met.metabolitos)
model.met = subset(model.met, met.media != 0)

X2 = data.frame(MO2.rankNorm[,2:(dim(MO2.rankNorm)[2])])
Y2 = MO2.rankNorm$BBDD2.delta
b = model.met

if (length(b) != 0) {
  stdout  = vector('character')
  con = textConnection('stdout', 'wr', local = TRUE)
  sink(con)
  for (i in 1:(dim(b)[1]))
  {
    if (i == 1) {cat("X2$model=", sep="")}
    if ((i < (dim(b)[1])) & (i > 0)) {cat("(","X2$",b[[i,1]],"*", b[[i,2]],")"," +", "\n", sep="")}
    if (i == (dim(b)[1])) {cat("(","X2$",b[[i,1]],"*", b[[i,2]],")", sep="")}
  }
  sink()
  
  close(con)
  
  stdout_full = paste(unlist(stdout), collapse =" ")
  stdout_full[1]
  
  eval(parse(text=stdout_full[1]))
  
  modelo = X2$model
  
  correr = cor.test(Y2,modelo, conf.level = .95)
} 
```

Coef con la BBDD completa y cv:
```{r}
set.seed(102)

cvfit = c()
coef_met_min = c()

X = data.frame(MO2.rankNorm[,2:(dim(MO2.rankNorm)[2])])
Y = MO2.rankNorm$BBDD2.delta

for (i in 1:10) {
  
  cvfit[[i]] = cv.glmnet(as.matrix(X), as.matrix(Y), family="gaussian", type.measure = "mse", alpha = alp)
  
  metabolitos.min = coef(cvfit[[i]], s = "lambda.min")
  
  coef.min = metabolitos.min@x
  coef.min = coef.min[-1]
  
  seleccion.min = c(metabolitos.min@i)
  nombres = as.data.frame(colnames(X))
  seleccion.min1 = nombres[seleccion.min,]
  seleccion.min1 = as.character(seleccion.min1)
  coef_met_min[[i]] = cbind(seleccion.min1, coef.min)
}

BBDD.met.min = data.frame(colnames(MO2.rankNorm))
colnames(BBDD.met.min)[colnames(BBDD.met.min)=="colnames.MO2.rankNorm."] = "seleccion.min1"

for (i in 1:10) {
  BBDD.met.min= merge(BBDD.met.min, coef_met_min[[i]], by = "seleccion.min1", all.x = T)
}
colnames(BBDD.met.min) = c("metabolitos","v1","v2","v3","v4","v5","v6","v7","v8","v9","v10")

export(BBDD.met.min, "10102020.CoefMO2v0.all.m50.5.xlsx")

save.image(file = "10102020.CoefMO2v0.all.m50.RData")
```

## 5. ENR binomial aumento+cambio vs disminución
ENR por grupo con el delta de BM categorizado: MO2.rankNorm
```{r}
MO2.rankNorm$cat = car::recode(MO2.rankNorm$BBDD2.delta, "0:100 = '0'; else = '1'")

table(MO2.rankNorm$cat)

MO2.rankNorm1 = MO2.rankNorm
MO2.rankNorm1$BBDD2.delta = NULL
MO2.rankNorm1$cat = as.factor(MO2.rankNorm1$cat)

set.seed(221020202)

train.data = c()
test.data = c()
cv = c()
bT = c()
predictions = c()
roc.min = c()

for (i in 1:10) {
  
  training = MO2.rankNorm1$cat %>% createDataPartition(p = 0.8, list = F)
  
  train.data[[i]] = MO2.rankNorm1[training,]
  test.data[[i]] = MO2.rankNorm1[-training,]
  
  cv[[i]] = train(cat ~ . , data = train.data[[i]],  method = "glmnet", metric = "Accuracy", trControl = trainControl("cv", number = 10), tuneLength = 10)
  
  bT[[i]] = cv[[i]]$bestTune
  
  x.test = model.matrix(cat~., test.data[[i]])[,-1] #Predictores
  y.test = test.data[[i]]$cat #Outcome
  
  predictions[[i]] = cv[[i]] %>% predict(x.test, type = "raw", s = "lambda.min")
  
  roc.min[[i]] = roc(as.numeric(y.test), as.numeric(predictions[[i]]), ci=TRUE)
}

l = 0.003312363 #cv[[i]]$bestTune$lambda
alp = 0.5 #cv[[i]]$bestTune$alpha

train_rows = c()
training = c()
validation = c()
modelito = c()
X = c()
Y = c()
X1 = c()
Y1 = c()
pred = c()
roc.min1 = c()

for (i in 1:10) {
  
  train_rows[[i]] = MO2.rankNorm1$cat %>% createDataPartition(p = 0.80, list = F)
  training[[i]] = MO2.rankNorm1[train_rows[[i]],]
  validation[[i]] = MO2.rankNorm1[-train_rows[[i]],]
  
  a=as.data.frame(training[[i]])
  X = as.matrix(a[,2:(dim(a)[2]-1)])
  Y = as.numeric(a$cat)
  
  b=as.data.frame(validation[[i]])
  X1 = as.matrix(b[,2:(dim(a)[2]-1)])
  Y1 = as.numeric(b$cat)
  
  modelito[[i]] = glmnet(X, as.matrix(Y), alpha = alp, lambda = l, family = "binomial")
  
  pred = predict(modelito[[i]], newx = X1, type = "class", s = l)
  
  roc.min1[[i]] = roc((Y1), as.numeric(pred), ci=TRUE)
}

coef_met_min = c()

for (i in 1:10) {
  
  metabolitos.min = coef(modelito[[i]], s = l)
  
  coef.min = metabolitos.min@x
  coef.min = coef.min[-1]
  
  seleccion.min = c(metabolitos.min@i)
  nombres = as.data.frame(colnames(X))
  seleccion.min1 = nombres[seleccion.min,]
  seleccion.min1 = as.character(seleccion.min1)
  coef_met_min[[i]] = cbind(seleccion.min1, coef.min)
  
}

(coef.met.min.AUC.log = print(coef_met_min))

BBDD.met.min = data.frame(colnames(MO2.rankNorm1))
colnames(BBDD.met.min)[colnames(BBDD.met.min)=="colnames.MO2.rankNorm1."] = "seleccion.min1"

for (i in 1:10) {
  BBDD.met.min= merge(BBDD.met.min, coef_met_min[[i]], by = "seleccion.min1", all.x = T)
}
colnames(BBDD.met.min) = c("metabolitos","v1","v2","v3","v4","v5","v6","v7","v8","v9","v10")

export(BBDD.met.min, "23102020.cat.MO2.rankNorm1.xlsx")
```

Validación interna:
```{r}
met = read_excel("22102020.cat.MO2.rankNorm1.xlsx")

model.met = data.frame(met$metabolitos, met$media)
str(model.met)
model.met$met.metabolitos = as.character(model.met$met.metabolitos)
model.met = subset(model.met, met.media != 0)

X2 = data.frame(MO2.rankNorm1[,2:(dim(MO2.rankNorm1)[2])])[-168]
Y2 = MO2.rankNorm1$cat
b = model.met
roc1 = c()
roc.m = matrix(NA,1,3)

if (length(b) != 0) {
  stdout  = vector('character')
  con = textConnection('stdout', 'wr', local = TRUE)
  sink(con)
  for (i in 1:(dim(b)[1]))
  {
    if (i == 1) {cat("X2$model=", sep="")}
    if ((i < (dim(b)[1])) & (i > 0)) {cat("(","X2$",b[[i,1]],"*", b[[i,2]],")"," +", "\n", sep="")}
    if (i == (dim(b)[1])) {cat("(","X2$",b[[i,1]],"*", b[[i,2]],")", sep="")}
  }
  sink()
  
  close(con)
  
  stdout_full = paste(unlist(stdout), collapse =" ")
  stdout_full[1]
  
  eval(parse(text=stdout_full[1]))
  
  modelo = X2$model
  
  roc1  =  roc(Y2, modelo, ci=TRUE)
  roc.m[,1] = roc1$auc
  roc.m[,2] = roc1$ci[1]
  roc.m[,3] = roc1$ci[3]
}    
```

Coef con la BBDD completa y cv:
```{r}
set.seed(231020202)

train_rows = c()
training = c()
validation = c()
modelito = c()
coef_met_min = c()

for (i in 1:10) {

  train_rows[[i]] = MO2.rankNorm1$cat %>% createDataPartition(p = 0.80, list = F)
  training[[i]] = MO2.rankNorm1[train_rows[[i]],]
  validation[[i]] = MO2.rankNorm1[-train_rows[[i]],]
  
  a=as.data.frame(training[[i]])
  X = as.matrix(a[,2:(dim(a)[2]-1)])
  Y = as.numeric(a$cat)
  
  b=as.data.frame(validation[[i]])
  X1 = as.matrix(b[,2:(dim(a)[2]-1)])
  Y1 = as.numeric(b$cat)
  
  modelito[[i]] = glmnet(X, as.matrix(Y), alpha = alp, lambda = l, family = "binomial")
  
  metabolitos.min = coef(modelito[[i]], s = "lambda.min")
  
  coef.min = metabolitos.min@x
  coef.min = coef.min[-1]
  
  seleccion.min = c(metabolitos.min@i)
  nombres = as.data.frame(colnames(X))
  seleccion.min1 = nombres[seleccion.min,]
  seleccion.min1 = as.character(seleccion.min1)
  coef_met_min[[i]] = cbind(seleccion.min1, coef.min)
  
}

BBDD.met.min = data.frame(colnames(MO2.rankNorm1))
colnames(BBDD.met.min)[colnames(BBDD.met.min)=="colnames.MO2.rankNorm1."] = "seleccion.min1"

for (i in 1:10) {
  BBDD.met.min= merge(BBDD.met.min, coef_met_min[[i]], by = "seleccion.min1", all.x = T)
}
colnames(BBDD.met.min) = c("metabolitos","v1","v2","v3","v4","v5","v6","v7","v8","v9","v10")

export(BBDD.met.min, "23102020.catMet.MO2.rankNorm1.xlsx")

save.image(file = "23102020.catMet.MO2.rankNorm1.RData")
```

## 6. ENR binomial aumento+cambio vs disminución con cambios en los metabolitos
ENR por grupo con el delta de BM categorizado: MO.deltas.rN1

```{r}
MO.deltas.rN1$cat = car::recode(MO.deltas.rN1$MO.deltas.rN.delta_bone, "0:100 = '0'; else = '1'")

table(MO2.rankNorm$cat)

MO.deltas.rN1$MO.deltas.rN.delta_bone = NULL
MO.deltas.rN1$cat = as.factor(MO.deltas.rN1$cat)

set.seed(231020203)

train.data = c()
test.data = c()
cv = c()
bT = c()
predictions = c()
roc.min = c()

for (i in 1:10) {
  
  training = MO.deltas.rN1$cat %>% createDataPartition(p = 0.8, list = F)
  
  train.data[[i]] = MO.deltas.rN1[training,]
  test.data[[i]] = MO.deltas.rN1[-training,]
  
  cv[[i]] = train(cat ~ . , data = train.data[[i]],  method = "glmnet", 
                  metric = "Accuracy", trControl = trainControl("cv", number = 10), tuneLength = 10)
  
  bT[[i]] = cv[[i]]$bestTune
  
  x.test = model.matrix(cat~., test.data[[i]])[,-1] #Predictores
  y.test = test.data[[i]]$cat #Outcome
  
  predictions[[i]] = cv[[i]] %>% predict(x.test, type = "raw", s = "lambda.min")
  
  roc.min[[i]] = roc(as.numeric(y.test), as.numeric(predictions[[i]]), ci=TRUE)
}

l = 0.01842466 #cv[[i]]$bestTune$lambda
alp = 0.4 #cv[[i]]$bestTune$alpha

train_rows = c()
training = c()
validation = c()
modelito = c()
X = c()
Y = c()
X1 = c()
Y1 = c()
pred = c()
roc.min1 = c()

for (i in 1:10) {
  
  train_rows[[i]] = MO.deltas.rN1$cat %>% createDataPartition(p = 0.80, list = F)
  training[[i]] = MO.deltas.rN1[train_rows[[i]],]
  validation[[i]] = MO.deltas.rN1[-train_rows[[i]],]
  
  a=as.data.frame(training[[i]])
  X = as.matrix(a[,2:(dim(a)[2]-1)])
  Y = as.numeric(a$cat)
  
  b=as.data.frame(validation[[i]])
  X1 = as.matrix(b[,2:(dim(a)[2]-1)])
  Y1 = as.numeric(b$cat)
  
  modelito[[i]] = glmnet(X, as.matrix(Y), alpha = alp, lambda = l, family = "binomial")
  
  pred = predict(modelito[[i]], newx = X1, type = "class", s = l)
  
  roc.min1[[i]] = roc((Y1), as.numeric(pred), ci=TRUE)
}

coef_met_min = c()

for (i in 1:10) {
  
  metabolitos.min = coef(modelito[[i]], s = l)
  
  coef.min = metabolitos.min@x
  coef.min = coef.min[-1]
  
  seleccion.min = c(metabolitos.min@i)
  nombres = as.data.frame(colnames(X))
  seleccion.min1 = nombres[seleccion.min,]
  seleccion.min1 = as.character(seleccion.min1)
  coef_met_min[[i]] = cbind(seleccion.min1, coef.min)
  
}

(coef.met.min.AUC.log = print(coef_met_min))

BBDD.met.min = data.frame(colnames(MO.deltas.rN1))
colnames(BBDD.met.min)[colnames(BBDD.met.min)=="colnames.MO.deltas.rN1."] = "seleccion.min1"

for (i in 1:10) {
  BBDD.met.min= merge(BBDD.met.min, coef_met_min[[i]], by = "seleccion.min1", all.x = T)
}
colnames(BBDD.met.min) = c("metabolitos","v1","v2","v3","v4","v5","v6","v7","v8","v9","v10")

export(BBDD.met.min, "23102020.cat.MO.deltas.rN1.xlsx")
```

Validación interna:
```{r}
met = read_excel("23102020.cat.MO.deltas.rN1.xlsx")

model.met = data.frame(met$metabolitos, met$media)
str(model.met)
model.met$met.metabolitos = as.character(model.met$met.metabolitos)
model.met = subset(model.met, met.media != 0)

X2 = data.frame(MO.deltas.rN1[,2:(dim(MO.deltas.rN1)[2])])[-168]
Y2 = MO.deltas.rN1$cat
b = model.met
roc1 = c()
roc.m = matrix(NA,1,3)

if (length(b) != 0) {
  stdout  = vector('character')
  con = textConnection('stdout', 'wr', local = TRUE)
  sink(con)
  for (i in 1:(dim(b)[1]))
  {
    if (i == 1) {cat("X2$model=", sep="")}
    if ((i < (dim(b)[1])) & (i > 0)) {cat("(","X2$",b[[i,1]],"*", b[[i,2]],")"," +", "\n", sep="")}
    if (i == (dim(b)[1])) {cat("(","X2$",b[[i,1]],"*", b[[i,2]],")", sep="")}
  }
  sink()
  
  close(con)
  
  stdout_full = paste(unlist(stdout), collapse =" ")
  stdout_full[1]
  
  eval(parse(text=stdout_full[1]))
  
  modelo = X2$model
  
  roc1  =  roc(Y2, modelo, ci=TRUE)
  roc.m[,1] = roc1$auc
  roc.m[,2] = roc1$ci[1]
  roc.m[,3] = roc1$ci[3]
  
}    
```

Coef con la BBDD completa y cv:
```{r}
set.seed(231020204)

train_rows = c()
training = c()
validation = c()
modelito = c()
coef_met_min = c()

for (i in 1:10) {
  
  train_rows[[i]] = MO.deltas.rN1$cat %>% createDataPartition(p = 0.80, list = F)
  training[[i]] = MO.deltas.rN1[train_rows[[i]],]
  validation[[i]] = MO.deltas.rN1[-train_rows[[i]],]
  
  a=as.data.frame(training[[i]])
  X = as.matrix(a[,2:(dim(a)[2]-1)])
  Y = as.numeric(a$cat)
  
  b=as.data.frame(validation[[i]])
  X1 = as.matrix(b[,2:(dim(a)[2]-1)])
  Y1 = as.numeric(b$cat)
  
  modelito[[i]] = glmnet(X, as.matrix(Y), alpha = alp, lambda = l, family = "binomial")
  
  metabolitos.min = coef(modelito[[i]], s = "lambda.min")
  
  coef.min = metabolitos.min@x
  coef.min = coef.min[-1]
  
  seleccion.min = c(metabolitos.min@i)
  nombres = as.data.frame(colnames(X))
  seleccion.min1 = nombres[seleccion.min,]
  seleccion.min1 = as.character(seleccion.min1)
  coef_met_min[[i]] = cbind(seleccion.min1, coef.min)
  
}

BBDD.met.min = data.frame(colnames(MO.deltas.rN1))
colnames(BBDD.met.min)[colnames(BBDD.met.min)=="colnames.MO.deltas.rN1."] = "seleccion.min1"

for (i in 1:10) {
  BBDD.met.min= merge(BBDD.met.min, coef_met_min[[i]], by = "seleccion.min1", all.x = T)
}
colnames(BBDD.met.min) = c("metabolitos","v1","v2","v3","v4","v5","v6","v7","v8","v9","v10")

export(BBDD.met.min, "23102020.catMet.MO.deltas.rN1.xlsx")

save.image(file = "23102020.catMet.MO.deltas.rN1.RData")
```

## 7. Logistic regression analyses examining the associations of a 1-SD increase in predictive models with high bone mass (median or higher)

```{r}
View(acti)

acti1 = c()

acti1$id = acti$id
acti1$MVPA = acti$MVPA1

acti1 = data.frame(acti1)

colnames(acti1) = c("allocation_no", "MVPA")

RL = data.frame(cbind(MO3, modelo, MO2.rankNorm1[169]))

RL$BMI = ((RL$v2_4.3 + RL$v2_4.4)/2)/((RL$v2_4.1 + RL$v2_4.2)/2)^2

colnames(RL)

RL = merge(RL, acti1, by = "allocation_no")

MO.RL = data.frame(cbind(RL[2:6], RL[12:15]))

m1 = glm(cat ~ modelo, data = MO.RL, family = "binomial")
summary(m1)

MO.RL$scorem1 = 0.3431 + (0.7021*MO.RL$modelo)

sd(MO.RL$scorem1)

MO.RL$score1SD = (MO.RL$scorem1)/1.990999

m1.1 = glm(cat ~ score1SD, data = MO.RL, family = "binomial")
summary(m1.1)
exp(m1.1$coefficients)
exp(confint(m1.1))

m2 = glm(cat ~ v2_bone_mass + v1_gq_11 + (v1_gq_1) + BMI + MVPA + modelo, data = MO.RL, family = "binomial")
summary(m2)

MO.RL$scorem2 = 9.571455 + (-0.704990 * MO.RL$v2_bone_mass) + (0.008701 * MO.RL$v1_gq_11) + 
  (-0.694098 * MO.RL$v1_gq_1) + (-0.192870 * MO.RL$BMI) + (-0.009234 * MO.RL$MVPA) + (0.768544 * MO.RL$modelo)

sd(MO.RL$scorem2, na.rm = T)

MO.RL$score2SD = (MO.RL$scorem2)/2.198457

m2.1 = glm(cat ~ score2SD, data = MO.RL, family = "binomial")
summary(m2.1)
exp(m2.1$coefficients)
exp(confint(m2.1))

m3 = glm(cat ~ v2_bone_mass + v1_gq_11 + (v1_gq_1) + BMI + MVPA, data = MO.RL, family = "binomial")
summary(m3)

MO.RL$scorem3 = 8.7058685 + (-0.2449550 * MO.RL$v2_bone_mass) + 
  (0.0269426 * MO.RL$v1_gq_11) + (0.0914429 * MO.RL$v1_gq_1) + 
  (-0.2919638  * MO.RL$BMI) + (0.0009696   * MO.RL$MVPA)

sd(MO.RL$scorem3, na.rm = T)

MO.RL$score3SD = (MO.RL$scorem3)/0.6313127

m3.1 = glm(cat ~ score3SD, data = MO.RL, family = "binomial")
summary(m3.1)
exp(m3.1$coefficients)
exp(confint(m3.1))
```

Curva de Roc
```{r}
roc.m1.1 =  roc(MO.RL$cat, MO.RL$score1SD, ci=TRUE)
roc.m2.1 =  roc(MO.RL$cat, MO.RL$score2SD, ci=TRUE)
roc.m3.1 =  roc(MO.RL$cat, MO.RL$score3SD, ci=TRUE)

rtest1 = roc.test(roc.m1.1,roc.m2.1, paired = T)
rtest2 = roc.test(roc.m1.1,roc.m3.1, paired = T)
rtest3 = roc.test(roc.m2.1,roc.m3.1, paired = T)
```