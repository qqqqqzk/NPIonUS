 library(dlnm); library(mgcv) ; library(MASS)
source('crosspred_re_refind.R', encoding = 'UTF-8')
df <- read.csv("Rt0427.csv")

# delete NPI i mplemented frequency lower than 5 times
aa <- c(0:30); bb <- c(4, 8, 10, 14, 16, 23, 25, 26, 27, 29); cc <- aa[!aa %in% bb]

# subset
for (i in cc) {assign(paste0("p", i), df[(i*15+4):(i*15+18)])}

# construct crossbasis
for (i in cc) {assign(paste0("cp", i), crossbasis(get(paste0("p", i)),lag=c(0,14),argvar=list('lin'),arglag=list('ps', df=5)))}

# cbPen
for (i in cc) {assign(paste0("cpp", i), cbPen(get(paste0("cp", i))))}

# GAM model
system.time({
  model <- gam(df$dif_Rt~cp0+cp1+cp2+cp3+cp5+cp6+cp7+cp9+cp11+cp12+cp13+cp15+cp17+cp18+cp19+cp20+cp21+cp22+cp24+cp28+cp30,
               family=gaussian,
               paraPen=list(cpp0=cp0, cpp1=cp1, cpp2=cp2, cpp3=cp3, cpp5=cp5, cpp6=cp6, cpp7=cp7,cpp9=cp9, cpp11=cp11, cpp12=cp12, cpp13=cp13, cpp15=cp15, cpp17=cp17, cpp18=cp18, cpp19=cp19, cpp20=cp20, cpp21=cp21, cpp22=cp22, cpp24=cp24, cpp28=cp28, cpp30=cp30), 
               method='REML')
})

# model performance
AIC(model)
summary(model)
cor(df$dif_Rt, fitted(model))^2

# crosspred
pcp0  <- crosspred(cp0, model,by=1,bylag=1,cen=0,cumul=T);pcp1  <- crosspred(cp1, model,by=1,bylag=1,cen=0,cumul=T);pcp2  <- crosspred(cp2, model,by=1,bylag=1,cen=0,cumul=T);pcp3  <- crosspred(cp3, model,by=1,bylag=1,cen=0,cumul=T);pcp5  <- crosspred(cp5, model,by=1,bylag=1,cen=0,cumul=T);pcp6  <- crosspred(cp6, model,by=1,bylag=1,cen=0,cumul=T);pcp7  <- crosspred(cp7, model,by=1,bylag=1,cen=0,cumul=T);pcp9  <- crosspred(cp9, model,by=1,bylag=1,cen=0,cumul=T);pcp11 <- crosspred(cp11,model,by=1,bylag=1,cen=0,cumul=T);pcp12 <- crosspred(cp12,model,by=1,bylag=1,cen=0,cumul=T);pcp13 <- crosspred(cp13,model,by=1,bylag=1,cen=0,cumul=T);pcp15 <- crosspred(cp15,model,by=1,bylag=1,cen=0,cumul=T);pcp17 <- crosspred(cp17,model,by=1,bylag=1,cen=0,cumul=T);pcp18 <- crosspred(cp18,model,by=1,bylag=1,cen=0,cumul=T);pcp19 <- crosspred(cp19,model,by=1,bylag=1,cen=0,cumul=T);pcp20 <- crosspred(cp20,model,by=1,bylag=1,cen=0,cumul=T);pcp21 <- crosspred(cp21,model,by=1,bylag=1,cen=0,cumul=T);pcp22 <- crosspred(cp22,model,by=1,bylag=1,cen=0,cumul=T);pcp24 <- crosspred(cp24,model,by=1,bylag=1,cen=0,cumul=T);pcp28 <- crosspred(cp28,model,by=1,bylag=1,cen=0,cumul=T);pcp30 <- crosspred(cp30,model,by=1,bylag=1,cen=0,cumul=T);

# export matfit
a <- as.data.frame(matrix(nrow=0,ncol=15))
for (i in cc) {a <- rbind(a, get(paste0('pcp',i))$matfit[2,])}
colnames(a)=c('lag0', 'lag1', 'lag2', 'lag3', 'lag4', 
              'lag5', 'lag6', 'lag7', 'lag8', 'lag9', 
              'lag10', 'lag11', 'lag12', 'lag13', 'lag14')

# export matse
b <- as.data.frame(matrix(nrow=0,ncol=15))
for (i in cc) {b <- rbind(b, get(paste0('pcp',i))$matse[2,])}
colnames(b)=c('lag0', 'lag1', 'lag2', 'lag3', 'lag4', 
              'lag5', 'lag6', 'lag7', 'lag8', 'lag9', 
              'lag10', 'lag11', 'lag12', 'lag13', 'lag14')

# export cumfit
c <- as.data.frame(matrix(nrow=0,ncol=15))
for (i in cc) {c <- rbind(c, get(paste0('pcp',i))$cumfit[2,])}
colnames(c)=c('lag0', 'lag1', 'lag2', 'lag3', 'lag4', 
              'lag5', 'lag6', 'lag7', 'lag8', 'lag9', 
              'lag10', 'lag11', 'lag12', 'lag13', 'lag14')

# export cumse
d <- as.data.frame(matrix(nrow=0,ncol=15))
for (i in cc) {d <- rbind(d, get(paste0('pcp',i))$cumse[2,])}
colnames(d)=c('lag0', 'lag1', 'lag2', 'lag3', 'lag4', 
              'lag5', 'lag6', 'lag7', 'lag8', 'lag9', 
              'lag10', 'lag11', 'lag12', 'lag13', 'lag14')

# export all
e <- as.data.frame(matrix(nrow=0,ncol=4))
for (i in cc) {e <- rbind(e, data.frame(get(paste0('pcp',i))$allfit[2], get(paste0('pcp',i))$allse[2], 
                                          get(paste0('pcp',i))$alllow[2], get(paste0('pcp',i))$allhigh[2]))} 
colnames(e) <- c('allfit', 'allse', 'alllow', 'allhigh') # allfit + 1.96*allse = allhigh # allfit=sum(pcp9$matfit[2,])
rownames(e) <- c(cc);e

# Save
write.csv(a, 'matfit0513.csv')
write.csv(b, 'matse0513.csv')
write.csv(c, 'cumfit0513.csv')
write.csv(d, 'cumse0513.csv')
write.csv(e, 'all0513.csv')