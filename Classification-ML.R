library(ggplot2)
library(corrplot)
library(glmnet)
library(scales)
library(gridExtra)
library(ROCR)
library(naniar)
library(VIM)
library(dplyr)
library(factoextra)
library(FactoMineR)
library(psych)
library(cluster)
library(tidyverse)
library(pROC)
library(caret)

#TODO 
#http://juliejosse.com/wp-content/uploads/2018/06/DataAnalysisMissingR.html#4)_multilevel_(mixed)_data_with_missing_values
#https://www.guru99.com/r-replace-missing-values.html
library(readxl)
data <- read_excel("C://Users//Lenovo//Desktop//Mes fichiers//ENSIIE//S3_ZAARI//Régression régularisée//Project//Data_Cortex_Nuclear.xls")
data$MouseID <- NULL
#Description of the dataset
#L'ensemble de donn�es contient les niveaux d'expression de 77 prot�ines 
#mesur�es dans le cortex c�r�bral de 8 classes de souris t�moins et 
#trisomiques expos�es au conditionnement de la peur contextuelle, 
#une t�che utilis�e pour �valuer l'apprentissage des souris.

#Lors de ce rapport nous allons se consacrer d'une exploration des 
#donn�es les bien comprendre et aussi pour bien comprendre les �tapes � suivre 
#afin de r�aliser une bonne classification du mod�le

#head(data,2)

#representing the number of missing entries in each variable


#The plot shows the real pourcentage of the missing and the present values
#for each variable of our dataset

vis_miss(data, sort_miss = FALSE) 


#Afin de r�soudre le probl�me des variables manquantes je vais remplacer les 
#valeurs par des moyennes par chaque groupe 
data <- data %>%
  group_by(class) %>%
  mutate_each(funs(replace(., which(is.na(.)), mean(., na.rm=TRUE)))) %>%
  as.data.frame()

data_new <- sapply(data, unclass) 


mycols <- c("#E59866", "#DC7633", "#BA4A00", "#BDC3C7","#717D7E",
            "#616A6B", "#85C1E9", "#154360")

ggplot(data,aes(class))+
  geom_bar(aes(fill=class),alpha=0.8)+
  scale_fill_manual(values = mycols) 



#Ensuite je vais transformer les variables quantitatives � des variables 
#qualitatives
data$class <- as.numeric(factor(data$class)) 
data$Behavior <- as.numeric(factor(data$Behavior)) - 1 
data$Treatment <- as.numeric(factor(data$Treatment))
data$Genotype <- as.numeric(factor(data$Genotype))


#Maintenant on va visualiser la correlation de nos donn�es
#Pour avoir une id�e claire sur la d�pendance des variables

mycol <- c("#85C1E9", "#154360","#E59866")
M <- cor(data)
corrplot(M,tl.col = "black",tl.cex=0.25,type="full",method="color",
         col = colorRampPalette(mycol)(100))
#



ggplot(data, aes(x = "",y="", fill = class)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  geom_text(aes(), color = "white")+
  scale_fill_manual(values = mycols) +
  theme_void()


h1 <- ggplot(data,aes(x=class))+
  geom_histogram(col='white',fill='blue',binwidth = 1)+labs(x='class')

h2 <- ggplot(data,aes(x=Behavior))+
  geom_histogram(col='red',fill='black',binwidth = 0.5)+labs(x='Behavior')
h3 <- ggplot(data,aes(x=Treatment))+
  geom_histogram(col='blue',fill='red',binwidth = 0.5)+labs(x='Treatment')
h4 <- ggplot(data,aes(x=Genotype))+
  geom_histogram(col='black',fill='white',binwidth = 0.5)+labs(x='Genotype')

grid.arrange(h1,h2,h3,h4,ncol=2,nrow=2)

#Visualisations et calculons l'ACP

X <- data %>% 
  select(starts_with("class"))


res.pca <- prcomp(data, scale = TRUE, center = TRUE)
#fviz_eig(res.pca)

p1 <- fviz_pca_ind(res.pca,
                   col.ind = "cos2", # Colorer par le cos2
                   gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                   repel = TRUE     
)

p2 <- fviz_pca_var(res.pca,
                   col.var = "contrib", 
                   gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                   repel = TRUE     
)

grid.arrange(p1,p2,ncol=2,nrow=1)

########################## Data science step

set.seed(123)
gap_stat <- clusGap(t(data[,1:71]), FUN = kmeans,
                    K.max = 10, B = 5)
fviz_gap_stat(gap_stat)

final <- kmeans(t(data[,1:71]), 5, nstart = 25)
fviz_cluster(final, data = t(data[,1:71]))


data_2 <- data[,1:71]
data_2$cluster <- final$cluster

aux <- as.list(data_2)
for(cluster in aux[1:length(data_2)]){
  plot(hist(final$cluster) , cluster ) #hist(cluster)
}

hist(final$cluster)

# Set up data matrices.
for(i in c(1:3))
  final$cluster[i]



train.predictors <- model.matrix(Behavior~., select(train, -c(class)))[,-1] %>% scale()
train.response   <- as.numeric(train$Behavior=="S/C")
test.predictors <- model.matrix(Behavior~., select(test, -c(class)))[,-1] %>% scale()
test.response <- as.numeric(test$Behavior=="S/C")

# Find optimal lambda via cross validation and fit model.
cv.lambda <- cv.glmnet(train.predictors, train.response, family="binomial")$lambda.min
cv.lambda
pen.log.mod <- glmnet(train.predictors, train.response, family="binomial", lambda=cv.lambda)
pen.log.mod


# Make predictions and test error rate.
resp <- predict(pen.log.mod, test.predictors)
resp
probs <- exp(resp)/(1+exp(resp))
preds <- ifelse(probs>=0.5, 1, 0)
print("Prediction Accuracy:")
print(mean(preds==test.response))




set.seed(123)
trainIndex <- createDataPartition(data$Behavior, p = .7, 
                                  list = FALSE, 
                                  times = 1)


hftrain<-data[trainIndex,]
hftest<-data[-trainIndex,]

full_model<-glm(Behavior~.,data=train,family=binomial("logit"), 
                maxit = 100)
summary(full_model)

step_model <- step(full_model, 
                   direction = "backward")
summary(step_model)

hftest$prob<-predict(step_model,hftest,type='response')
hftest$Behavior_pred <-ifelse(hftest$prob>=0.5,1,0)
table(hftest$Behavior_pred,hftest$Behavior)
plot(roc(hftest$Behavior,hftest$Behavior_pred),col='red')
auc(roc(hftest$Behavior,hftest$Behavior_pred))

hftest %>% 
  arrange(prob) %>% 
  mutate(rank=rank(prob),Event=ifelse(prob>=0.5,'SC','CS')) %>% 
  ggplot(aes(rank,prob))+
  geom_point(aes(color=Event),size=4.5)+
  geom_smooth(method = "glm", 
              method.args = list(family = "binomial")) +
  ggtitle('Stepwise Logistic Regression',subtitle='Predicting Behavior')





# Variables pr�dictives
X <- model.matrix(Behavior~., data= train)[,-1]
# Variable de r�sultat (Target)
Y<- train$Behavior

glmnet(X, Y, alpha = 1 , lambda = NULL)

set.seed(123)
model_cross_valid <- cv.glmnet(X, Y, alpha = 1)
best_lambda <-model_cross_valid$lambda.min
best_lambda
plot(model_cross_valid)
model <- glmnet(X, Y, alpha = 1, lambda = best_lambda)
# affichons les coefficients de la r�gression
coef(model)
X_test <- model.matrix(Behavior ~., test)[,-1]
pred_lasso <- predict(model, s=best_lambda, newx= X_test)

#Ridge regression
set.seed(123)
model_cross_valid <- cv.glmnet(X, Y, alpha = 0)
# on affiche la meilleure valeur de lamnbda
best_lambda <-model_cross_valid$lambda.min
best_lambda

model <- glmnet(X, Y, alpha = 0, lambda = best_lambda)
X_test <- model.matrix(Behavior ~., test)[,-1]
pred_Ridge <- predict(model,s=best_lambda,newx= X_test)


RMSE_lasso <- sqrt(mean(pred_lasso-test$Behavior)^2)
RMSE_Ridge <- sqrt(mean(pred_Ridge-test$Behavior)^2)
SSE_lasso <- sum((pred_lasso - test$Behavior)^2)
SSE_Ridge <- sum((pred_Ridge - test$Behavior)^2)
SST <- sum((test$Behavior - mean(test$Behavior))^2)
#RMSE et R squared
results_RMSE <- data.frame(RMSE_RIDGE=RMSE_Ridge,
                           RMSE_LASSO=RMSE_lasso,
                           R_SQUARED_RIDGE=1-(SSE_Ridge/SST),
                           R_SQUARED_lasso=1-(SSE_lasso/SST))
results_RMSE







data_Split <- sort(sample(nrow(data),nrow(data)*.70))
train <- data[data_Split,]
test <- data[-data_Split,]



# Variables pr�dictives
X <- model.matrix(Behavior~., data= train)[,-1]
# Variable de r�sultat (Target)
Y<- train$Behavior

ridge2 <- glmnet(X,Y,family="binomial",standardize=FALSE,alpha=0)
plot(ridge2$lambda.min)





plot(ridge , xvar = "lambda")



summary()


predict.kmeans <- function(newdata, object){
  centers <- object$centers
  n_centers <- nrow(centers)
  dist_mat <- as.matrix(dist(rbind(centers, newdata)))
  dist_mat <- dist_mat[-seq(n_centers), seq(n_centers)]
  max.col(-dist_mat)}

set.seed(123)
gap_stat <- clusGap(t(data[,1:71]), FUN = kmeans,
                    K.max = 10, B = 5)
fviz_gap_stat(gap_stat)

final <- kmeans(t(data[,1:71]), 5, nstart = 25)
fviz_cluster(final, data = t(data[,1:71]))

par(mfrow = c(3,3))

plot_1 <- ggplot(data, aes(eval(parse(text = proteins[1])))) +
  geom_histogram(aes(fill=..count..), color = "black", alpha = 0.5) +
  scale_fill_gradient("Count", low="green", high="red") +
  labs(title = proteins[i],
       x = "Expression Level",
       y = "Count") +
  theme_light()
cat("#### ", proteins[i], "\n")
print(plot)
cat('\n\n')


plot_2 <- ggplot(data, aes(eval(parse(text = proteins[15])))) +
  geom_histogram(aes(fill=..count..), color = "black", alpha = 0.5) +
  scale_fill_gradient("Count", low="green", high="red") +
  labs(title = proteins[i],
       x = "Expression Level",
       y = "Count") +
  theme_light()
cat("#### ", proteins[i], "\n")
print(plot)
cat('\n\n')

grid.arrange(plot_1,plot_2)



p1 <- ggplot(data_1, aes(pCASP9_N)) +
  geom_histogram(aes(fill=..count..), color = "black", alpha = 0.5) +
  scale_fill_gradient("Count", low="darkblue", high="yellow") +
  labs(title = "Protein: pCASP9_N",
       x = "Expression Level",
       y = "Count")


p2 <- ggplot(data_1, aes(x = class, y = pCASP9_N)) +
  geom_boxplot(aes(fill = class), alpha = 0.7) +
  scale_fill_manual(values=c("#E59866", "darkblue", "#BA4A00", "#BDC3C7","#717D7E",
                             "#616A6B","#DC7633", "purple"))+
  stat_summary(fun.y=mean, colour="blue", geom="point", size=2) 

grid.arrange(p1,p2,nrow=2,ncol=1)










