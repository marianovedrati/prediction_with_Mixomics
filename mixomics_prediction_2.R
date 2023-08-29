###############DIABLO PER FARE PREVISIONI
df <- read.csv2("count_norm_importance.csv",row.names=1)


DESC_table <- read.csv2("DESC_SEX_gemma.csv")


m<-match(rownames(df), DESC_table$Sample.ID)
DESC_table<-DESC_table[m,]

citochine <- read.csv2("Desc_tab_citochine.csv")




#####tengo solo i campioni per qui ho le anlisi delle citochine

k<-citochine$Sample.ID %in% DESC_table$Sample.ID
citochine_k<-citochine[k,]
# verifico l'ordine delle samples
citochine_k$Sample.ID==DESC_table$Sample.ID
# ordino le semples nello stesso modo
m<-match(citochine_k$Sample.ID, DESC_table$Sample.ID)
DESC_table<-DESC_table[m,]
DESC_table$Sample.ID==citochine_k$Sample.ID




# ordino le semples nello stesso modo

df<-df[rownames(df) %in% DESC_table$Sample.ID,]
n <- match(DESC_table$Sample.ID, rownames(df))
df <- df[n,]
rownames(df)==citochine_k$Sample.ID #verifico che l'ordine sia corretto


###########SPLITTO I DATI

install.packages("caTools")
library(caTools)

ind <- sample.split(Y =DESC_table$Sample.ID, SplitRatio = 0.7)



#subsetting into Train data
DESC_train = DESC_table[ind,]

df_train <- df[ind,]

cito_train <- citochine_k[ind,]


#subsetting into Test data
DESC_test = DESC_table[!ind,]


df_test <- df[!ind,]

cito_test <- citochine_k[!ind,]



#######

rownames(cito_train) <- cito_train$Sample.ID
cito_train <- cito_train[,12:20]
cito_train<- as.matrix(cito_train)



#######miXomics dopo avere applicato filtro
library(mixOmics)
#help("mixOmics")


X <- list(gene_expression = (
  df_train), 
  citochine = cito_train) 

Y <- as.factor(DESC_train$cc_21)




##DESIGN MATRIX
design <- matrix(0.1, ncol = length(X), nrow = length(X), 
                 dimnames = list(names(X), names(X)))
diag(design) <- 0
design 


###CONTROLLO CON UNA SPARSEPLS LA CORRELAZIONE TRA I VARI DATASETA PER SECGLIERE DESIGN PIù GIUST E ACCURATO
res1.pls.EXP_CIT <- pls(X$gene_expression,X$citochine, ncomp = 1)
cor(res1.pls.EXP_CIT$variates$X, res1.pls.EXP_CIT$variates$Y)

# riscrivo la matrice di correlazione del design in seguito ai riultati della PLS
design <- matrix(0.58, ncol = length(X), nrow = length(X), 
                 dimnames = list(names(X), names(X)))
diag(design) <- 0
design




###nscelta numero di componti
diablo.cc <- block.plsda(X, Y, ncomp = 6, design = design)
diablo.cc$design
set.seed(35)
perf.diablo.cc <- perf(diablo.cc,validation = "Mfold",  folds = 7, nrepeat = 20)
plot(perf.diablo.cc)
perf.diablo.cc$choice.ncomp$WeightedVote



ncomp <- perf.diablo.cc$choice.ncomp$WeightedVote["Overall.BER","centroids.dist"]


###abbiamo valutato il numero di componenti principali da prendere e dopo aver fatto perf() function,
#abbiamo deciso di prendere 2 componenti.




#Number of variables to be selected

set.seed(40)

test.keepX <- list(gene_expression = c(5:10,seq(11,69,5)), 
                   citochine = c(1:9))
tune.diablo.cc <- tune.block.splsda(X, Y ,ncomp = 2,
                                    test.keepX = test.keepX, design = design,
                                    validation = "Mfold", folds = 5, nrepeat = 5, 
                                    dist = "centroids.dist")

list.keepX <- tune.diablo.cc$choice.keepX
list.keepX


#$gene_expression
#[1] 61 5

#$citochine
#[1] 3 1 

####scelgo il modello finale
diablo.cc <- block.splsda(X, Y, ncomp = 2, 
                          keepX = list.keepX, design = design)
diablo.cc$design
selectVar(diablo.cc)

#grafico e tipi di plot

plotDiablo(diablo.cc,ncomp = 2)
plotIndiv(diablo.cc, legend = TRUE)

plotArrow(diablo.cc)
plotVar(diablo.cc)

circosPlot(diablo.cc)
network(diablo.cc)
plotLoadings(diablo.cc)



###heatmap con cim function
cimDiablo(diablo.cc,margins = c(2,15),trim= F)







##############dopo aver costruito il modello su cui fittare i miei dati provo a fare la predizione


# Prepare test set data: here one block (proteins) is missing
df_test <- as.matrix(df_test)
rownames(cito_test) <- cito_test$Sample.ID
cito_test <- cito_test[,12:20]
cito_test<- as.matrix(cito_test)





X.test <- list(gene_expression = 
                 (df_test), 
               citochine = cito_test) 


predict.diablo.tcga <- predict(diablo.cc, newdata = X.test)

confusion.mat.tcga <- get.confusion_matrix(truth = DESC_test$cc_21, 
                                           predicted = predict.diablo.tcga$WeightedVote$centroids.dist[,2])
confusion.mat.tcga






# The warning message will inform us that one block is missing

#predict.diablo # List the different outputs


get.BER(confusion.mat.tcga)



