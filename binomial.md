### 2.1.1. Entrenamiento de los modelos:

#### Estimacion del valor de alpha:

```{r}
set.seed(12)

rows = sample(nrow(diabetes))
diabetes = diabetes [rows,]
folds = cut(seq(1,nrow(diabetes)),breaks=698,labels=FALSE)

train.data = c()
test.data = c()
cv = c()
bT = c()

for (i in 1:698) {
  training = which(folds==i, arr.ind=TRUE)
  
  train.data[[i]] = diabetes [-training, ]
  test.data[[i]] = diabetes [training, ] 
  
  cv[[i]] = train(as.factor(cens) ~ . , data = train.data[[i]],  method = "glmnet", trControl = trainControl("cv", number = 5), tuneLength = 5)
  
  bT[[i]] = cbind(cv[[i]]$bestTune,accuracy = max(cv[[i]]$results$Accuracy))
}
```

Selección del mejor Alpha y lambda.

```{r}
bT = data.frame(matrix(unlist(bT), ncol = 3, byrow= T))
names(bT) = c("alpha","lambda","accuracy")
l = bT$lambda[bT$accuracy== max(bT$accuracy)] #0.01597965
#l = 0.01597965
alp = bT$alpha[bT$accuracy== max(bT$accuracy)] #0.325
#alp = 0.325
```

#### Selección de los metabolitos

```{r}
set.seed(22)

type.lambda=c("lambda.min")
#type.lambda=c("lambda.1se")

fit_full=list()
cvfit_full=list()
variables_values_full=list()

full_data=as.data.frame(diabetes)
X.full=full_data[,2:(dim(full_data)[2])]
X.full=as.matrix(X.full)
Y.full = full_data$cens

for (i in 1:10)
{
  cat("##################################", '\n')
  cat("Iteration WHOLE DATASET ", i, '\n')
  cat("##################################", '\n')
  
  fit_full[[i]] = glmnet(X.full, as.numeric(Y.full), lambda = l, family="binomial", alpha=alp)
  plot(fit_full[[i]])
  
  cvfit_full[[i]] = cv.glmnet(X.full, as.numeric(Y.full), family="binomial", type.measure = "class", nfolds = 10, alpha=alp)
  plot(cvfit_full[[i]])
  
  list_metabo_full=coef(cvfit_full[[i]], s = type.lambda)
  
  values_metabo_full=list_metabo_full@x
  values_metabo_full=values_metabo_full[-1]
  
  listado_selected_full=c(list_metabo_full@i)
  variables_model_full=as.data.frame(colnames(X.full))
  variables_get_full=variables_model_full[listado_selected_full,]
  variables_get_full=as.character(variables_get_full)
  
  variables_values_full[[i]]=cbind(variables_get_full, values_metabo_full)
}

for (i in 1:10)
{
  list_metabo_full=coef(cvfit_full[[i]], s = type.lambda)
  
  values_metabo_full=list_metabo_full@x
  values_metabo_full=values_metabo_full[-1]
  
  listado_selected_full=c(list_metabo_full@i)
  variables_model_full=as.data.frame(colnames(X.full))
  variables_get_full=variables_model_full[listado_selected_full,]
  variables_get_full=as.character(variables_get_full)
  
  variables_values_full[[i]]=cbind(variables_get_full, values_metabo_full)
}
```

#### Extracción de los metabolitos

```{R}
inter=as.data.frame(diabetes)
variables2=as.data.frame(colnames(inter[,2:(length(inter))])) #Save variables

colnames(variables2)=c("variables_get_full")
variables2$check=c('NA')

abc_set = data.frame(colnames(diabetes))
colnames(abc_set)[colnames(abc_set)=="colnames.diabetes."] = "variables_get_full"

for (i in 1:10) {
  abc_set= merge(abc_set, variables_values_full[[i]], by = "variables_get_full", all.x = T)
}
colnames(abc_set) = c("variables_get_full","v1","v2","v3","v4","v5","v6","v7","v8","v9","v10")

metabos_NA_model=as.data.frame(10 - rowSums(is.na(abc_set)))
metabos_NA_model=cbind(abc_set$variables_get_full, metabos_NA_model)
colnames(metabos_NA_model)=c("variables_get_full", "Join")
metabos_NA_model=metabos_NA_model[order(metabos_NA_model$variables_get_full),]

abc_set=merge(abc_set, metabos_NA_model, by = "variables_get_full", all=T)
colnames(abc_set)=c("variables_get_full", "It1", "It2", "It3", "It4", "It5", "It6", "It7", "It8", "It9", "It10", "JOIN")

abc_set_values=abc_set[,2:11]
str(abc_set[,2:11], list.len=(ncol(abc_set))-3)
indx = sapply(abc_set[,2:11], is.character)
abc_set_values[indx] = lapply(abc_set_values[indx], function(x) as.numeric(as.character(x)))

abc_full=cbind("Metabolites"=abc_set[,1], abc_set_values, "JOIN"=abc_set[,12])

abc_full$mean=apply(abc_full[,2:11], 1, function(x) { mean(x, na.rm=TRUE) })
abc_full$sd=apply(abc_full[,2:11], 1, function(x) { sd(x, na.rm=TRUE) })

mean_values_metabo = abc_full[ which(abc_full$JOIN==10), c(1,13,14)]

export(mean_values_metabo, "Only coef_binomial_diabetes.xlsx")

unlist(variables_values_full[[1]])

mean_values_metabo_matrix=as.matrix(mean_values_metabo)

match_names = read_delim("G:/trabajo/ARTICULOS/BBDD/PREDIMED/names_vars_280218.csv", 
    delim = ";", escape_double = FALSE, trim_ws = TRUE)

new_abc_full=merge(abc_full, match_names, by = "Metabolites")
```

#### Figuras

```{r}
new_mean_values_metabo=merge(mean_values_metabo, match_names, by = "Metabolites")

mean_values_metabo_sorted=mean_values_metabo[order(mean_values_metabo$mean),c(1,2,3)]
dim(mean_values_metabo_sorted)[1]

new_mean_values_metabo_sorted=merge(mean_values_metabo_sorted, match_names, by = "Metabolites")
new_mean_values_metabo_sorted=data.frame("Metabolites"=new_mean_values_metabo_sorted$name, "mean"=new_mean_values_metabo_sorted$mean, "sd"=new_mean_values_metabo_sorted$sd)
new_mean_values_metabo_sorted=new_mean_values_metabo_sorted[order(new_mean_values_metabo_sorted$mean),c(1,2,3)]
dim(new_mean_values_metabo_sorted)[1]

selected_mean_values_metabo_sorted_negative=new_mean_values_metabo_sorted[(new_mean_values_metabo_sorted$mean<0), ]

selected_mean_values_metabo_sorted_positive=new_mean_values_metabo_sorted[(new_mean_values_metabo_sorted$mean>=0), ]

minimo=min(selected_mean_values_metabo_sorted_negative$mean)
maximo=max(selected_mean_values_metabo_sorted_positive$mean)

plot_negative=ggplot(selected_mean_values_metabo_sorted_negative, aes(x = mean, y = reorder(Metabolites, -mean))) +
  geom_point() + scale_y_discrete(position = "left") + xlab("Coefficient value") + ylab("Metabolites") +
  scale_x_continuous(limits = c(-0.4, 0)) + geom_errorbarh(aes(xmax = mean + sd, xmin = mean - sd, height = .2)) + 
  theme_bw()

plot_positive=ggplot(selected_mean_values_metabo_sorted_positive, aes(x = mean, y = reorder(Metabolites, +mean))) +
  geom_point() + scale_y_discrete(position = "right") + xlab("Coefficient value") + ylab("Metabolites") +
  scale_x_continuous(limits = c(0, 1)) + geom_errorbarh(aes(xmax = mean + sd, xmin = mean - sd, height = .2)) + 
  theme_bw()

require(gridExtra)
png(filename = "coef_diabetes.png", width = 30, height = 20, res=300, unit="cm")
grid.arrange(plot_negative, plot_positive, ncol=2)
dev.off()

dim(selected_mean_values_metabo_sorted_negative)[1]
dim(selected_mean_values_metabo_sorted_positive)[1]
```
### 2.1.2. Validación interna

#### Validación Training-testing

```{r}
rows = sample(nrow(diabetes))
diabetes = diabetes[rows, ]
folds = cut(seq(1,nrow(diabetes)),breaks=10,labels=FALSE)

set_train=list()
set_test=list()
X.train_saved=list()
Y.train_saved=list()
X.test_saved=list()
Y.test_saved=list()    
fit=list()
cvfit=list()
variables_values=list()

for (i in 1:10)
{
  cat("##################################", '\n')
  cat("Iteration TRAINING-testing ", i, '\n')
  cat("##################################", '\n \n')

  sample = which(folds==i,arr.ind=TRUE)
  set_train[[i]] = diabetes[-sample, ] 
  set_test[[i]]  = diabetes[sample, ]
  
  a=as.data.frame(set_train[[i]])
  train=a[,2:(dim(a)[2])]
  X.train = train
  X.train = as.matrix(X.train)
  Y.train = a$cens
  
  a=as.data.frame(set_test[[i]])
  test=a[,2:(dim(a)[2])]
  X.test = test
  X.test = as.matrix(X.test)
  Y.test = a$cens
  
  X.train_saved[[i]]=X.train
  Y.train_saved[[i]]=Y.train
  X.test_saved[[i]]=X.test
  Y.test_saved[[i]]=Y.test    
  
  fit[[i]] = glmnet(X.train, as.numeric(Y.train), family="binomial", alpha=alp)
  plot(fit[[i]])
  
  cvfit[[i]] = cv.glmnet(X.train, as.numeric(Y.train), family="binomial", type.measure = "class", nfolds = 10, alpha=alp)
  plot(cvfit[[i]])
  
  list_metabo=coef(cvfit[[i]], s = type.lambda)
  values_metabo=list_metabo@x
  values_metabo=values_metabo[-1]
  
  listado_selected=c(list_metabo@i)
  variables_model=as.data.frame(colnames(X.train))
  variables_get=variables_model[listado_selected,]
  variables_get=as.character(variables_get)
  variables_values[[i]]=cbind(variables_get, values_metabo)
}
```

#### AUC para el training-testing

```{r}
corr_test_pearson=list()
corr_test_spearman=list()
standard_roc = c()
correlaciones_tt=matrix(NA,10,4)
AUCs_tt = matrix(NA,10,3)

type.lambda=c("lambda.min")
#type.lambda=c("lambda.1se")

si_model = 0 

for (i in 1:10)
{
  list_metabo=coef(cvfit[[i]], s = type.lambda)
  
  values_metabo=list_metabo@x
  values_metabo=values_metabo[-1]
  
  listado_selected=c(list_metabo@i)
  variables_model=as.data.frame(colnames(X.train))
  variables_get=variables_model[listado_selected,]
  variables_get=as.character(variables_get)
  
  variables_values[[i]]=cbind(variables_get, values_metabo)
  
  if (dim(variables_values[[i]])[1] == 0)
  {
    cat("No metabolites found for iteration", i, '\n')
  }
  
  if (dim(variables_values[[i]])[1] != 0)
  {
    standard_roc[[i]] = roc(Y.test_saved[[i]], as.numeric(predict(cvfit[[i]], X.test_saved[[i]], type = "response")), ci=TRUE)
    corr_test_pearson[[i]]=cor.test(Y.test_saved[[i]], predict(cvfit[[i]],newx=X.test_saved[[i]], s=type.lambda), method="pearson")
    corr_test_spearman[[i]]=cor.test(Y.test_saved[[i]], predict(cvfit[[i]],newx=X.test_saved[[i]], s=type.lambda), method="spearman")
    
    correlaciones_tt[i,1]=round(corr_test_pearson[[i]]$estimate,2)
    correlaciones_tt[i,2]=round(corr_test_pearson[[i]]$conf.int[1],2)
    correlaciones_tt[i,3]=round(corr_test_pearson[[i]]$conf.int[2],2)
    correlaciones_tt[i,4]=round(corr_test_pearson[[i]]$p.value,3)
    
    AUCs_tt[i,1]=round(standard_roc[[i]]$auc[1],2)
    AUCs_tt[i,2]=round(standard_roc[[i]]$ci[1],2)
    AUCs_tt[i,3]=round(standard_roc[[i]]$ci[2],2)

    si_model = cbind(si_model, i)
  }
}

ci.mean(correlaciones_tt[,1],normal=T) # 0.42 [0.37;0.47]
ci.mean(AUCs_tt[,1],normal=T) #0.80 [0.77;0.82]
```
#### AUC del modelo

```{r}
correlacion_all = c()
AUCs_all = c()

X2 = data.frame(BBDDmet[,199:580])
Y2 = BBDDmet$cens
b = mean_values_metabo

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
  
  correlacion_all = cor.test(Y2,modelo, conf.level = .95) #0.58 (0.54, 0.63)
  AUCs_all = roc(Y2, modelo, ci=TRUE) #0.90 (0.88, 0.93)
} 

score_diab = data.frame(cbind(id = BBDDmet$id, cens = BBDDmet$cens, diab_est = modelo))
```

