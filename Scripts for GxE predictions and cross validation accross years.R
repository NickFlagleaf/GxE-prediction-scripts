rm(list=ls())

#Install the required packages
pckgs<-installed.packages()
if(!"BGGE"%in%pckgs){install.packages("BGGE")}
if(!"rrBLUP"%in%pckgs){install.packages("rrBLUP")}
if(!"foreach"%in%pckgs){install.packages("foreach")}
if(!"iterators"%in%pckgs){install.packages("iterators")}
if(!"parallel"%in%pckgs){install.packages("parallel")}
if(!"doParallel"%in%pckgs){install.packages("doParallel")}
if(!"reshape2"%in%pckgs){install.packages("reshape2")}
if(!"corpcor"%in%pckgs){install.packages("corpcor")}
library(BGGE)
library(rrBLUP)
library(foreach)
library(iterators)
library(parallel)
library(doParallel)
library(reshape2)
library(corpcor)

# Several data sources must be used for this analysis but cannot be made directly available here. However the objects required are detailed:
# 'data' is a data frame of the multi environment trial data with the column names 'Variety','Harv_Year','Site' and 'Env' as the first four columns. 
#       The 14th to the 19th column have the quality trait data.
# 'A' is the additive kinship matrix based on pedigree data. Row and column names must be in data$Variety 
# 'markers.set' is a matrix of SNP marker data where row names must be in data$Variety and column names are marker names.
#               Marker values are codes as allele dosage, ie. 0,1 or 2.
# 'all_soil_properties' is a dataframe of soil environmental covariables. Row names must be in data$Env for each trial environment and column names are the soil properties.
# 'all.site.weather' is a dataframe of weather environmental covariables. The first column is the environment names that must be in data$Env. Colum names are the monthly weather variable names.


#Two custom functions are also used:

#Make ERM kernel:
Kernel.W <-function(W,type_kernel,h=1){
  if(type_kernel=="linear"){
    K<-tcrossprod(W)/ncol(W)
    return(K)
  }
  if(type_kernel=="gauss"){
    D2<-as.matrix(dist(W)^2)
    K<-exp(-h*D2/median(D2))
    return(K)
  }
}


#Function for constructing GxE kernels with the trial design matrices:
GE.function<-function (pheno_geno,KG,KE=NULL,model) {
  K<-KG
  IDs<-as.character(unique(pheno_geno$GID))
  K<-K[IDs,IDs]
  pheno_geno$GID <-factor(x=pheno_geno$GID,levels=rownames(K), ordered=T)
  Zg <- model.matrix(~factor(pheno_geno$GID)-1)  # genotype design matrix
  Ze <- model.matrix(~factor(pheno_geno$Location)-1)  # environment design matrix 
  ne<-apply(Ze,2,sum)
  K1 <- Zg%*% K%*%t(Zg)
  ZEZE<-tcrossprod(Ze)
  
  if (model=="MM"){
    ETA<-list(list(Kernel=ZEZE,Type="D"),list(Kernel=K1,Type="D"))
    result<-list(ETA,ne)
    return(result)
  }
  if (model=="MDs"){
    K2<-K1*ZEZE
    ETA<-list(list(Kernel=ZEZE,Type="D"),list(Kernel=K1,Type="D"),list(Kernel=K2,Type="D"))
    result<-list(ETA,ne)
    return(result)
  }
  if (model=="EMM"){
    ETA<-list(list(Kernel=ZEZE,Type="D"),list(Kernel=K1,Type="D"))
    for(j in 1:length(KE)){
      K<-KE[[j]]
      IDs<-as.character(unique(pheno_geno$Location))
      K<-K[IDs,IDs]
      pheno_geno$Location <-factor(x=pheno_geno$Location,levels=rownames(K), ordered=T)
      Ze <- model.matrix(~factor(pheno_geno$Location)-1)  # environment design mat
      K_env<-Ze%*%K%*%t(Ze)
      temp<-list(Kernel=K_env,Type="D")
      ETA[[2+j]]<-temp
    }
    result<-list(ETA,ne)
    return(result)
  }
  if (model=="EMDs"){
    K2<-K1*ZEZE
    ETA<-list(list(Kernel=ZEZE,Type="D"),list(Kernel=K1,Type="D"),list(Kernel=K2,Type="D"))
    for(j in 1:length(KE)){
      K<-KE[[j]]
      IDs<-as.character(unique(pheno_geno$Location))
      K<-K[IDs,IDs]
      pheno_geno$Location <-factor(x=pheno_geno$Location,levels=rownames(K), ordered=T)
      Ze <- model.matrix(~factor(pheno_geno$Location)-1)  # environment design mat
      K_env<-Ze%*%K%*%t(Ze)
      temp<-list(Kernel=K_env,Type="D")
      ETA[[3+j]]<-temp
    }
    result<-list(ETA,ne)
    return(result)
  }
  
  if (model=="RNMM"){
    ETA<-list(list(Kernel=ZEZE,Type="D"),list(Kernel=K1,Type="D"))
    for(j in 1:length(KE)){
      K<-KE[[j]]
      IDs<-as.character(unique(pheno_geno$Location))
      K<-K[IDs,IDs]
      pheno_geno$Location <-factor(x=pheno_geno$Location,levels=rownames(K), ordered=T)
      Ze <- model.matrix(~factor(pheno_geno$Location)-1)  # environment design mat
      K_env<-Ze%*%K%*%t(Ze)
      temp<-list(Kernel=K_env,Type="D")
      ETA[[1+2*j]]<-temp
      K_wg<-K1*K_env
      temp<-list(Kernel=K_wg,Type="D")
      ETA[[2+2*j]]<-temp
    }
    result<-list(ETA,ne)
    return(result)
  }
  
  
  if (model=="RNMDs"){
    K2<-K1*ZEZE
    ETA<-list(list(Kernel=ZEZE,Type="D"),list(Kernel=K1,Type="D"),list(Kernel=K2,Type="D"))
    for(j in 1:length(KE)){
      K<-KE[[j]]
      IDs<-as.character(unique(pheno_geno$Location))
      K<-K[IDs,IDs]
      pheno_geno$Location <-factor(x=pheno_geno$Location,levels=rownames(K), ordered=T)
      Ze <- model.matrix(~factor(pheno_geno$Location)-1)  # environment design mat
      K_env<-Ze%*%K%*%t(Ze)
      temp<-list(Kernel=K_env,Type="D")
      ETA[[2*j+2]]<-temp
      K2<-K1*K_env
      temp<-list(Kernel=K2,Type="D")
      ETA[[2*j+3]]<-temp
      
    }
    ne<-apply(Ze,2,sum)
    result<-list(ETA,ne)
    return(result)
    
  }
}






traits<-data[,14:19]
colnames(traits)<-c("Protein","HFN","Specific weight","Zeleny","Chopin W","Chopin P/L") #Specify trait names

#Make GRM and combine with A to get K
G<-make.positive.definite(A.mat(markers.set)) #Linear G matrix
G<-G[colnames(G)%in%colnames(A),colnames(G)%in%colnames(A)] #Remove varieties not in pedigree A
Agg=matrix(NA,ncol(G),nrow(G))#Make subset of A that's in G
Agg=A[colnames(G),colnames(G)]
meanG=mean(G)
meandiagG=mean(diag(G))
meanAgg=mean(Agg)
meandiagAgg=mean(diag(Agg))
b=(meandiagAgg-meanAgg)/(meandiagG-meanG)
a=meandiagAgg-meandiagG*b
cat("a=",a,"b=",b) 
G<-a+b*G #Adjust G 

#Make weigted K matrix with A and G:
w=0.25 #only one value for weighting was used
K<-matrix(NA,nrow=nrow(A),ncol=ncol(A),dimnames = dimnames(A))
K[colnames(G),colnames(G)]<-(w*Agg)+((1-w)*G)
K[is.na(K)]<-A[is.na(K)]
K<-make.positive.definite(K)
K_G<-K


all.site.weather<-all.site.weather[,colnames(all.site.weather)[-grep(x=colnames(all.site.weather),pattern = "Pressure")]]

envs <- unique(data$Env)

#Check the environmental data matched the names in the dataset
all_soil_properties<- as.data.frame(all_soil_properties)
all_soil_properties<- all_soil_propoerties[all_soil_properties$Env %in% envs,]
all.site.weather<- all.site.weather[all.site.weather$Env %in% envs,]
all.site.weather <- all.site.weather[complete.cases(all.site.weather),] # resulting in 377 (139 for the last 5 years)

# building the matrix of environmental covariates (EC)
EC <- merge(all.site.weather,all_soil_properties,by='Env')
dim(EC)
rownames(EC) <- EC$Env


EC <- as.matrix(EC[,-1])
dim(EC) # 315 environments]


# updated selected environments (with EC)
envs <- rownames(EC)

all.traits.out.of.fold.preds<-list() #Make list to output predictions to
#For each trait ----
for (i in 1:ncol(traits)) { 
  print(paste("Starting",colnames(traits)[i]))
  envs <- rownames(EC)
  data.nona      <- droplevels(data[data$Env %in% envs,])
  #remove NAs for trait and get only varieties in markers dataset
  data.nona      <- data.nona[!is.na(data.nona[,13+i]),]
  Y.data         <- data.nona[data.nona$Variety%in%rownames(K),c(4,1,13+i)]
  Y.data$Env     <- factor(as.character(Y.data$Env)) #Clean up factor levels
  Y.data$Variety <- factor(as.character(Y.data$Variety)) #cleanup factor levels
  
  colnames(Y.data)[3]<-"trait"
  years<-as.character(substr(Y.data$Env,1,4))
  sites<-as.character(substr(Y.data$Env,5,1000))
  Y.data<-cbind(Y.data,years,sites)
  
  # step 1: creating the W matrix
  W_matrix <- scale(EC[rownames(EC) %in% envs,],center = T,scale = T)
  K_E<-Kernel.W(W=W_matrix,type_kernel="gauss")	
  
  print("building kernels")
  pheno_geno<-Y.data
  colnames(pheno_geno)<-c("Location","GID","y","years","sites")
  common<-intersect(unique(pheno_geno$Location),unique(rownames(K_E)))
  pheno_geno<-pheno_geno[pheno_geno$Location%in%common,]
  K_E<-K_E[rownames(K_E)%in%common,colnames(K_E)%in%common]
  K_E<-list(K_E)
  
  cat("|")
  M1<-GE.function(pheno_geno=pheno_geno,KG=K_G,model="MM")
  cat("|")
  M2<-GE.function(pheno_geno=pheno_geno,KG=K_G,model="MDs")
  #cat("|")
  M3<-GE.function(pheno_geno=pheno_geno,K=K_G,KE=K_E,model="RNMM")
  cat("|")
  M4<-GE.function(pheno_geno=pheno_geno,K=K_G,KE=K_E,model="RNMDs")
  
  #Cross validate across years in parallel----
  cl <- makeCluster(length(levels(pheno_geno$years)),outfile=paste(Sys.time(),"out.txt"))
  registerDoParallel(cl)
  all.year.out.of.fold.preds<-foreach(j = 1:length(levels(pheno_geno$years)),.combine="rbind",.packages = c("BGGE"),.verbose=T) %dopar% {
    
    print(paste("Cross validation for",colnames(traits)[i],levels(years)[j]))
    test<-which(pheno_geno$years==levels(pheno_geno$years)[j])#make test set index
    yna<-pheno_geno$y
    yna[test]<-NA #mask data for the test year
    yna[pheno_geno$GID%in%pheno_geno$GID[test]]<-NA #mask data for genotypes in common with the those in the test year 
    
    #Fit each model
    print(paste(colnames(traits)[i],levels(pheno_geno$years)[j],"fitting M1"))
    M1.fit<-BGGE(y=yna,K=M1[[1]],ne=M1[[2]],ite=10000,burn=2000,thin=2)
    print(paste(colnames(traits)[i],levels(pheno_geno$years)[j],"fitting M2"))
    M2.fit<-BGGE(y=yna,K=M2[[1]],ne=M2[[2]],ite=10000,burn=2000,thin=2)
    print(paste(colnames(traits)[i],levels(pheno_geno$years)[j],"fitting M3"))
    M3.fit<-BGGE(y=yna,K=M3[[1]],ne=M3[[2]],ite=10000,burn=2000,thin=2)
    print(paste(colnames(traits)[i],levels(pheno_geno$years)[j],"fitting M4"))
    M4.fit<-BGGE(y=yna,K=M4[[1]],ne=M4[[2]],ite=10000,burn=2000,thin=2)
    
    #Get out of fold predictions
    out.of.fold.preds<-pheno_geno[test,]
    out.of.fold.preds<-cbind(out.of.fold.preds,M1.fit$yHat[test],
                             M2.fit$yHat[test],
                             M3.fit$yHat[test],
                             M4.fit$yHat[test])
    return(out.of.fold.preds)
  }
  stopCluster(cl)
  stopImplicitCluster()
  save(all.year.out.of.fold.preds,file=paste(colnames(traits)[i],"out.of.fold.preds.RData")) #save as you go
  all.traits.out.of.fold.preds[[i]]<-all.year.out.of.fold.preds
  print(paste("Finished",colnames(traits)[i],"!!!"))
}
save(all.traits.out.of.fold.preds,file="all.traits.out.of.fold.preds.RData") #save everything
write.csv(all.traits.out.of.fold.preds,file="All RKHS out of fold predictions.csv")



all.traits.pred.r<-list()
all.traits.year.means<-list()
#Get all within site prediction acuracies----
  for (i in 1:length(all.traits.out.of.fold.preds)) {
    sites<-unique(all.traits.out.of.fold.preds[[i]]$Location)
    all.sites.r<-matrix(NA,nrow=length(sites),ncol=6)
    colnames(all.sites.r)<-c("Trait","Env","M1","M2","M3","M4")
    all.sites.r[,1]<-colnames(traits)[i]
    all.sites.r[,2]<-as.character(sites)
    
    for(k in 1:length(sites)){
      is.site<-all.traits.out.of.fold.preds[[i]]$Location==sites[k]
      if(sum(all.traits.out.of.fold.preds[[i]]$Location[is.site]==sites[k])>6){
        all.sites.r[k,3]<-cor(all.traits.out.of.fold.preds[[i]]$y[is.site],all.traits.out.of.fold.preds[[i]]$`M1.fit$yHat[test]`[is.site],
                              use="pairwise.complete.obs")
        all.sites.r[k,4]<-cor(all.traits.out.of.fold.preds[[i]]$y[is.site],all.traits.out.of.fold.preds[[i]]$`M2.fit$yHat[test]`[is.site],
                              use="pairwise.complete.obs")
        all.sites.r[k,5]<-cor(all.traits.out.of.fold.preds[[i]]$y[is.site],all.traits.out.of.fold.preds[[i]]$`M3.fit$yHat[test]`[is.site],
                              use="pairwise.complete.obs")
        all.sites.r[k,6]<-cor(all.traits.out.of.fold.preds[[i]]$y[is.site],all.traits.out.of.fold.preds[[i]]$`M4.fit$yHat[test]`[is.site],
                              use="pairwise.complete.obs")
      }
    }
    all.traits.pred.r[[i]]<-all.sites.r
    all.traits.year.means[[i]]<-aggregate(all.traits.out.of.fold.preds[[i]],by = list(all.traits.out.of.fold.preds[[i]]$years),FUN = function(x) mean(na.omit(x)))
  }
#print mean within site prediction accuracies:
lapply(all.traits.pred.r,function(l) apply(l,2,function(x) mean(na.omit(as.numeric(x[3:6])))))






##########Run for all years---- ######################
cl <- makeCluster(6,outfile=paste(Sys.time(),"full.mods.out.txt"))
registerDoParallel(cl)
all.trait.fitted.models<-foreach(i = 1:ncol(traits),.packages = c("BGGE"),.verbose=T) %dopar% {
  print(paste("Starting",colnames(traits)[i]))
  envs <- rownames(EC)
  
  data.nona      <- droplevels(data[data$Env %in% envs,])
  #remove NAs for trait and get only varieties in markers dataset
  data.nona      <- data.nona[!is.na(data.nona[,13+i]),]
  Y.data         <- data.nona[data.nona$Variety%in%rownames(K),c(4,1,13+i)]
  Y.data$Env     <- factor(as.character(Y.data$Env)) #Clean up factor levels
  Y.data$Variety <- factor(as.character(Y.data$Variety)) #cleanup factor levels
  
  colnames(Y.data)[3]<-"trait"
  years<-as.character(substr(Y.data$Env,1,4))
  sites<-as.character(substr(Y.data$Env,5,1000))
  Y.data<-cbind(Y.data,years,sites)
  
  # step 1: creating the W matrix
  W_matrix <- scale(EC[rownames(EC) %in% envs,],center = T,scale = T)
  K_E<-Kernel.W(W=W_matrix,type_kernel="gauss")	
  
  
  print("building kernels")
  pheno_geno<-Y.data
  colnames(pheno_geno)<-c("Location","GID","y","years","sites")
  common<-intersect(unique(pheno_geno$Location),unique(rownames(K_E)))
  pheno_geno<-pheno_geno[pheno_geno$Location%in%common,]
  K_E<-K_E[rownames(K_E)%in%common,colnames(K_E)%in%common]
  K_E<-list(K_E)
  cat("|")
  M1<-GE.function(pheno_geno=pheno_geno,KG=K_G,model="MM")
  cat("|")
  M2<-GE.function(pheno_geno=pheno_geno,KG=K_G,model="MDs")
  cat("|")
  M3<-GE.function(pheno_geno=pheno_geno,K=K_G,KE=K_E,model="RNMM")
  cat("|")
  M4<-GE.function(pheno_geno=pheno_geno,K=K_G,KE=K_E,model="RNMDs")
  
  print(paste("Fitting models",colnames(traits)[i]))
  cat("M1")
  M1.fit<-BGGE(y=pheno_geno$y,K=M1[[1]],ne=M1[[2]],ite=10000,burn=2000,thin=2)
  cat("M2")
  M2.fit<-BGGE(y=pheno_geno$y,K=M2[[1]],ne=M2[[2]],ite=10000,burn=2000,thin=2)
  cat("M3")
  M3.fit<-BGGE(y=pheno_geno$y,K=M3[[1]],ne=M3[[2]],ite=10000,burn=2000,thin=2)
  cat("M4")
  M4.fit<-BGGE(y=pheno_geno$y,K=M4[[1]],ne=M4[[2]],ite=10000,burn=2000,thin=2)
  cat("Finished model fitting!")
  all.model.fits<-list(MM=M1.fit,MDs=M2.fit,RNMM=M3.fit,RNMDs=M4.fit)
  return(all.model.fits)
}
save(all.trait.fitted.models,file="all.trait.fitted.models.RData")



#get variance components for each trait and model:
all.trait.var.comps<-list() #make list for all traits
for (i in 1:length(all.trait.fitted.models)) {
  #Get variance components for each model:
  var.comps<-matrix(NA,ncol=4,nrow=6)
  colnames(var.comps)<-c("MM","MDs","RNMM","RNMDs")  
  rownames(var.comps)<-c("E","gE","w","g","gw","Residuals")
  var.comps["E","MM"]<-all.trait.fitted.models[[i]]$MM$K$K1$varu
  var.comps["g","MM"]<-all.trait.fitted.models[[i]]$MM$K$K2$varu
  var.comps["Residuals","MM"]<-all.trait.fitted.models[[i]]$MM$varE
  
  var.comps["E","MDs"]<-all.trait.fitted.models[[i]]$MDs$K$K1$varu
  var.comps["g","MDs"]<-all.trait.fitted.models[[i]]$MDs$K$K2$varu
  var.comps["gE","MDs"]<-all.trait.fitted.models[[i]]$MDs$K$K3$varu
  var.comps["Residuals","MDs"]<-all.trait.fitted.models[[i]]$MM$varE
  
  var.comps["E","RNMM"]<-all.trait.fitted.models[[i]]$RNMM$K$K1$varu
  var.comps["g","RNMM"]<-all.trait.fitted.models[[i]]$RNMM$K$K2$varu
  var.comps["w","RNMM"]<-all.trait.fitted.models[[i]]$RNMM$K$K3$varu
  var.comps["gw","RNMM"]<-all.trait.fitted.models[[i]]$RNMM$K$K4$varu
  var.comps["Residuals","RNMM"]<-all.trait.fitted.models[[i]]$RNMM$varE
  
  var.comps["E","RNMDs"]<-all.trait.fitted.models[[i]]$RNMDs$K$K1$varu
  var.comps["g","RNMDs"]<-all.trait.fitted.models[[i]]$RNMDs$K$K2$varu
  var.comps["gE","RNMDs"]<-all.trait.fitted.models[[i]]$RNMDs$K$K3$varu
  var.comps["w","RNMDs"]<-all.trait.fitted.models[[i]]$RNMDs$K$K4$varu
  var.comps["gw","RNMDs"]<-all.trait.fitted.models[[i]]$RNMDs$K$K5$varu
  var.comps["Residuals","RNMDs"]<-all.trait.fitted.models[[i]]$RNMDs$varE
  
  all.trait.var.comps[[i]]<-var.comps #Add to list
}
names(all.trait.var.comps)<-colnames(traits)

for (i in 1:length(all.trait.var.comps)) { #fill in 0s
  all.trait.var.comps[[i]][is.na(all.trait.var.comps[[i]])]<-0
}

for (i in 1:length(all.trait.var.comps)){ #convert to %
  all.trait.var.comps[[i]]<-apply(all.trait.var.comps[[i]],2,function(x) (x/sum(x))*100 )
}





