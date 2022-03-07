library(EpiEstim)
library(plyr)
us0427 <- read.csv("us0427.csv")
# generate subset for each state
us0427$index <- as.Date(us0427$index)
#============================= Rt 2019-2020 =============================
us0427 <- rename(us0427, c(inc_splines = 'I')) 
for (i in unique(us0427$REGION)) {
  assign((i), subset(us0427, REGION==i, select=c('index', 'I')))
}
for (i in unique(us0427$REGION)) {
  assign(paste0('Rt', i), estimate_R(get(i)[complete.cases(get(i)), ], method = "parametric_si",
                                     config = make_config(list(mean_si = 2.85, std_si = 0.93))))
}
a <- as.data.frame(matrix(nrow=0,ncol=0))  # 空白data.frame
for (i in unique(us0427$REGION)) {
a <- rbind(a, data.frame(rep(i, length(get(paste0('Rt', i))$R$`Mean(R)`)), 
                    get(paste0('Rt', i))$R$`Mean(R)`))
} 

#============================= Rt previous =============================
us0427 <- rename(us0427, c(avg = 'I'))
for (i in unique(us0427$REGION)) {
  assign((i), subset(us0427, REGION==i, select=c('index', 'I')))
}

for (i in unique(us0427$REGION)) {
  assign(paste0('Rt', i), estimate_R(get(i)[complete.cases(get(i)), ], method = "parametric_si",
                                     config = make_config(list(mean_si = 2.85, std_si = 0.93))))
} 
b <- as.data.frame(matrix(nrow=0,ncol=0))  # 空白data.frame
for (i in unique(us0427$REGION)) {
  b <- rbind(b, data.frame(rep(i, length(get(paste0('Rt', i))$R$`Mean(R)`)), 
                           get(paste0('Rt', i))$R$`Mean(R)`))
}
# export
write.csv(b, file='rt_p.csv')
write.csv(a, file='rt_2020.csv')