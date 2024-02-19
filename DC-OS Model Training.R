#BiocManager::install("GenomicFeatures")
#devtools::install_github('yikeshu0611/ggDCA')
#install.packages("survminer")
library(readr)  
library(survival)
library(ggplot2)
library(ggpubr)
library(magrittr)
library(survminer)
library(dplyr)
#train.data <- read.csv(file = "R6.CSV",header = T)
train.data <- read.csv(file = "imputed_data.CSV",header = T)
#data1 <- na.omit(data0)
groupedclin <-as.data.frame(train.data)
colnames(groupedclin)
#####首先告诉软件，哪些是你的连续变量（必须），但可能引入NAs#####
for (i in names(groupedclin)[c(67:94)]) {groupedclin[,i] <- as.numeric(groupedclin[,i])} # data[,"sex"]<- as.factor(data[,"sex"])
str(groupedclin)
#########查看缺失数据############
##检查是否存在NA值  
na_values <- is.na(groupedclin)  
print(na_values)  
##检查是否存在Inf值(异常字符)  
inf_values <- is.infinite(unlist(groupedclin)) #将列表转换为向量，然后检查向量中是否存在无穷大的值 
print(inf_values) 
data1 <-as.data.frame(groupedclin)
complete.cases(data1)#查看完整数据有多少
table(complete.cases(data1))#缺失数据比例
unique(unlist(lapply(data0,function(x) which(is.na(x)))))#显示缺失数据所在行数
missing.percent <- unlist(lapply(groupedclin,function(x) sum(is.na(x))))/nrow(groupedclin)#查看每个变量缺失比例（将超过10%缺失的删除）
missing.percent
library(VIM)#缺失值的可视化
matrixplot(groupedclin)#各变量分别看
aggr(groupedclin)#整体可视化
aggr(groupedclin, numbers = TRUE)
##删除缺失值##
#groupedclin <- na.omit(groupedclin)
table(groupedclin$V4)
############COX法单因素分析，优势是无需分组转换###########
cox4 <- coxph(Surv(V68, V66) ~ V48, data=groupedclin)
summary(cox4)
#groupedclin$V32 <- as.numeric(groupedclin$V32) 
sfit <- survfit(Surv(V68, V66) ~ V48, data=groupedclin)
#log rank检验
survdiff(Surv(V68, V66) ~ V48, data=groupedclin) 
#groupedclin$V31 <- as.factor(groupedclin$V31)
#ggsurvplot(sfit)
plot(sfit)
#########-----------###################
#批量进行单因素Cox生存分析
covariates <- c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"
                ,"V11","V12","V13","V14","V15","V16","V17","V18","V19","V20",
                "V21","V22","V23","V24","V25","V26","V27","V28","V29","V30",
                "V31","V32","V33","V34","V35","V36","V37","V38","V39","V40",
                "V41","V42","V43","V44","V45","V46","V47","V48","V49","V50",
                "V51","V52","V53","V54","V55","V56","V57","V58","V59","V60",
                "V61","V62","V63","V64","V69","V70",
                "V71","V72","V73","V74","V75","V76","V77","V78","V79","V80",
                "V81","V82","V83","V84","V85","V86","V87","V88","V89","V90",
                "V91","V92","V93","V94")
univ_formulas <- sapply(covariates,function(x) as.formula(paste('Surv(V68, V66)~', x)))
#univ_formulas
attach(groupedclin) 
#univ_models <- lapply( univ_formulas, function(x){coxph(x, data = info)})
univ_models <- lapply(univ_formulas, function(x){coxph(x, data = groupedclin)})
#univ_results
#要整理输出，注意要一个因素一行才行，像treatmentArm这样子chr类型的，会被当成不同种类而有不同的p-value（见下文多因素的输出）
#chr类型会被堪称不同种类，所以要变成numeric（看age）
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"],2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)","wald.test", "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))
as.data.frame(res)
#########-----------#######
####----------------------多因素Cox回归----------------------------####
#####---单变量初筛建模---#####
N_cox1 <- coxph(Surv(V68, V66) ~ V2+V6+V11+V13+V19+V21+V23+V25+V31+V35+V37+V38
                +V43+V44+V45+V46+V48+V49+V50+V51+V52+V53+V54
                +V55+V57+V58+V59+V60+V63+V64+V69+V70+V72+V73
                +V75+V76+V79+V80+V81+V82+V83+V85+V87+V89+V90
                +V91+V92+V94, data= groupedclin)
summary(N_cox1)#展示结果cox
#--COX筛选变量再确认--#
Univariate_cox <- coxph(Surv(V68, V66)~ V11+V25+V31+V43+V44+V64+V82+V91+V92+V94, data= groupedclin)
summary(Univariate_cox)#展示结果cox
library(gtsummary)
tbl_regression(Univariate_cox,exponentiate = TRUE) #一键列表，同下
#直接输出表格#
library(tableone)
ShowRegTable(Univariate_cox)
#write.csv(ShowRegTable(N_cox),file = "cox.csv")
#计算C-Index#
sum.surv <- summary(Univariate_cox)
c.index.univariate <- sum.surv$concordance
c.index.univariate
####绘制森林图####
# 加载所需的R包 # 
library(survival)  
library(survminer)  
# 拟合Cox比例风险模型 #  
#Univariate_cox <- coxph(Surv(V68, V66)~ V11+V31+V43+V44+V64+V82+V91+V92
#                        +V94, data= groupedclin)  
# 绘制森林图 # 
ggforest(Univariate_cox, # 使用coxph函数的结果  
         main = "Hazard ratio",  
         cpositions = c(0.02,-0.15, 0.25),  
         fontsize = 0.8,  
         refLabel = "reference",   
         noDigits = 2)
#####---AJCC 8th edition建模---#####
AJCC_cox <- coxph(Surv(V68, V66)~ V45+V46+V47+V48, data= groupedclin)
summary(AJCC_cox)#展示结果cox
#列表显示#
library(gtsummary)
tbl_regression(AJCC_cox,exponentiate = TRUE) #一键列表，同下
#计算C-Index#
sum.surv <- summary(AJCC_cox)
c.index.AJCC <- sum.surv$concordance
c.index.AJCC
####绘制森林图####
library(survival)  
library(survminer)  
# 绘制森林图 # 
ggforest(AJCC_cox, # 使用coxph函数的结果  
         main = "Hazard ratio",  
         cpositions = c(0.02,-0.15, 0.25),  
         fontsize = 0.8,  
         refLabel = "reference",   
         noDigits = 2)
#####---全变量纳入建模---#####
data2 <- train.data[, -c(65, 67)]
head(data2)
N_cox2 <- coxph(Surv(V68, V66)~ ., data= data2)
summary(N_cox2)#展示结果cox
#Full-cox model建模#
#Full_cox <- coxph(Surv(V68, V66)~ V11+V24+V30+V31+V35+V36+V37+V43+V56+V64+V77+V78+V82+V91+V94, data= groupedclin)
Full_cox <- coxph(Surv(V68, V66)~ V11+V30+V31+V37+V43+V56+V64+V77+V78+V82+V91+V94, data= groupedclin)
summary(Full_cox)#展示结果cox
library(gtsummary)
tbl_regression(Full_cox,exponentiate = TRUE) #一键列表，同下
#计算C-Index#
sum.surv <- summary(Full_cox)
c.index.full <- sum.surv$concordance
c.index.full
####绘制森林图####
# 加载所需的R包 # 
library(survival)  
library(survminer)  
# 拟合Cox比例风险模型 #  
#Full_cox <- coxph(Surv(V68, V66)~ V11+V30+V31+V37+V43+V56+V64
#+V77+V78+V82+V91+V94, data= groupedclin)  
# 绘制森林图 # 
ggforest(Full_cox, # 使用coxph函数的结果  
         main = "Hazard ratio",  
         cpositions = c(0.02,-0.15, 0.25),  
         fontsize = 0.8,  
         refLabel = "reference",   
         noDigits = 2)
#####################--Lasso回归初筛变量--##################
library(glmnet)
library(survival)
library(foreign)
#1、首先告诉软件，哪些是你的分类变量,分类变量全部转换为数值型，注意：异常值类型转换时可能会出现NAs#
for (i in names(groupedclin)[c(1:64)]) {groupedclin[,i] <- as.numeric(groupedclin[,i])} # data[,"sex"] <- as.factor(data[,"sex"])
str(groupedclin)
#groupedclin <- na.omit(groupedclin) ##必须不含NAs值
table(groupedclin$V48) #logist回归模型每个变量的阳性结局样本数至少是其10倍（EPV，event per variable），估算纳入模型的最大变量数
##########——Lasso回归建模——###########
#2、新建data frame
#定义x 
x = as.matrix(groupedclin[,c(1:64,69:94)])
#定义y
y= Surv(groupedclin$V68, groupedclin$V66==1)
######-----Lasso纳入全部变量-----######
fit <- glmnet(x,y,family="cox",alpha = 1)  
#解释偏差百分比以及相应的λ值
print(fit)
#解释偏差不再随着λ值的增加而减小，因此而停止
#画出收缩曲线图
#plot(fit,label = FALSE)
#plot(fit,label = TRUE,xvar = "lambda")#系数值如何随着λ的变化而变化
plot(fit,label = FALSE,xvar = "lambda")#删除变量标记
plot(fit,label = FALSE,xvar = "lambda")#偏差与系数之间的关系图
plot(fit,label = FALSE,xvar = "dev")
#指定lamda给出相应的筛选变量
lasso.coef <- predict(fit, type = "coefficients",s = 0.205300) # s指定
lasso.coef
####---lasso cross-validation（内部交叉验证）获取最优λ值---####
#glmnet包在使用cv.glmnet()估计λ值时，默认使用10折交叉验证。在K折交叉验证中，
#数据被划分成k个相同的子集（折），#每次使用k-1个子集拟合模型，然后使用剩下的
#那个子集做测试集，最后将k次拟合的结果综合起来（一般取平均数），确定最后的参数。
#在这个方法中，每个子集只有一次用作测试集。在glmnet包中使用K折交叉验证非常容易，
#结果包括每次拟合的λ值和响应的MSE。默认设置为α=1。
#进行交叉验证，选择最优的惩罚系数lambada
##############################
set.seed(123)
#cv.fit <- cv.glmnet(x,y,family="cox",alpha = 1)
#画出收缩系数图（偏差）
cv.fit <- cv.glmnet(x,y,family="cox",
                    alpha = 1,
                    type.measure = "deviance",
                    nfolds = 5)#deviance和λ的关系
plot(cv.fit)
###########绘制C-index和log(λ)的关系########## 
library(ggplot2)
#cv.fit <- cv.glmnet(x, y, family = "cox", alpha = 1, type.measure = "deviance", nfolds = 5)
cv.fit1 <- cv.glmnet(x, y, family = "cox", type.measure = "C")
#-Extract plotting data-#
plot_data <- data.frame(lambda = log(cv.fit1$lambda), cvm = cv.fit1$cvm)
#-Calculate confidence intervals-#
conf_intervals <- cv.fit1$cvm + cv.fit1$cvsd
plot_data$upper <- conf_intervals
plot_data$lower <- cv.fit1$cvm - cv.fit1$cvsd
#-Plot custom graph with confidence intervals(用cv.fit线标记)-#
ggplot(plot_data, aes(x = lambda, y = cvm)) +
  geom_point(color = "red", fill = "red", shape = 21, size = 2) +
  geom_vline(xintercept = log(cv.fit$lambda.min), linetype = "dashed", color = "gray") +
  geom_vline(xintercept = log(cv.fit$lambda.1se), linetype = "dashed", color = "gray") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.05, color = "gray") +
  labs(x = "log(λ)", y = "C-index", title = "") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black")) +
  scale_x_continuous(breaks = seq(-7, max(plot_data$lambda), by = 1)) +
  scale_y_continuous(breaks = seq(0, max(plot_data$cvm), by = 0.05)) +
  xlab("log(λ)") +
  ylab("C-index")
####提取变量（最小lambda，最佳）####
cv.fit$lambda.min
coef(fit,s=cv.fit$lambda.min) #用cv.fit/fit模型筛选变量一致，s代表标准
#predict(cv.fit,type='coefficient',s=cv.fit$lambda.min) #功能同上
#计算OR值#
# 使用的lambda.min值来拟合模型#  
fit_optimal_lambda <- glmnet(x, y, family = "cox", alpha = 1, lambda = cv.fit$lambda.min)  
#获取系数#  
coefs <- coef(fit_optimal_lambda)  
#计算OR#
exp(coefs)
####提取变量（1倍标准误，最优）####
cv.fit$lambda.1se
coef(fit,s=cv.fit$lambda.1se)
#predict(cv.fit,type='coefficient',s=cv.fit$lambda.1se) #功能同上
#计算OR值，注意：glmnet函数不能计算95%置信区间#
exp(coef(cv.fit)) #OR值
print(cv.fit) #变量个数汇总
#指定lamda给出相应的筛选变量
#predict(cv.fit,type='coefficient',s=0.1)
####把最佳和最优线加到收缩曲线图中####
plot(fit,xvar="lambda",label=FALSE)
abline(v=log(c(cv.fit$lambda.min,cv.fit$lambda.1se)),lty=2)
####Lasso-cox变量筛查####
##Lasso-cox建模##
Lasso_cox <- coxph(Surv(V68, V66) ~ V43+V44+V45+V48+V57+V64+V73+V82+V91+V94, data= groupedclin)
summary(Lasso_cox)#展示结果cox
##Lasso-cox筛选变量再确认##
Lasso_cox <- coxph(Surv(V68, V66)~ V43+V44+V57+V64+V73+V82+V91+V94, data= groupedclin)
summary(Lasso_cox)#展示结果cox
library(gtsummary)
tbl_regression(Lasso_cox,exponentiate = TRUE) #一键成表，同下
#计算C-Index#
sum.surv <- summary(Lasso_cox)
c.index.lasso <- sum.surv$concordance
c.index.lasso
####绘制森林图####
# 加载所需的R包 # 
library(survival)  
library(survminer)  
# 拟合Cox比例风险模型 #  
#Lasso_cox <- coxph(Surv(V68, V66)~ V43+V44+V57+V64+V73+V82+V91+V94, data= groupedclin) 
# 绘制森林图 # 
ggforest(Lasso_cox, # 使用coxph函数的结果  
         main = "Hazard ratio",  
         cpositions = c(0.02,-0.15, 0.25),  
         fontsize = 0.8,  
         refLabel = "reference",   
         noDigits = 2)
############################---COX模型评价---############################
Lasso_cox <- coxph(Surv(V68, V66)~ V43+V44+V57+V64+V73+V82+V91+V94, data= groupedclin)
Univariate_cox <- coxph(Surv(V68, V66)~ V11+V25+V31+V43+V44+V64+V82+V91+V92+V94, data= groupedclin)
Full_cox <- coxph(Surv(V68, V66)~ V11+V30+V31+V37+V43+V56+V64+V77+V78+V82+V91+V94, data= groupedclin)
AJCC_cox <- coxph(Surv(V68, V66)~ V45+V46+V47+V48, data= groupedclin)
###----模型间的比较----###
library(lmtest)
anova(Lasso_cox,Full_cox,test="Chisq") #模型有差异
anova(Lasso_cox,Univariate_cox,test="Chisq")
anova(Lasso_cox,AJCC_cox,test="Chisq")
#AIC越小，模型越优
AIC(Lasso_cox,Univariate_cox,Full_cox,AJCC_cox) 
#Likelihood ratio test，模型相似性检验
lrtest(Lasso_cox,Full_cox) 
lrtest(Lasso_cox,Univariate_cox)
lrtest(Lasso_cox,AJCC_cox)
##--模型NRI比较（C-index补充）--##
library(survival)  
# 使用Cox模型拟合数据  
#Univariate_cox <- coxph(Surv(V68, V66)~ V11+V31+V43+V44+V64+V82+V91+V92+V94, data= groupedclin)
#Full_cox <- coxph(Surv(V68, V66)~ V11+V30+V31+V37+V43+V56+V64+V77+V78+V82+V91+V94, data= groupedclin)
#Lasso_cox <- coxph(Surv(V68, V66) ~ V43+V44+V45+V48+V57+V64+V73+V82+V91+V94, data= groupedclin) 
# 计算风险分数  
Univariate_risk <- predict(Univariate_cox, type = "risk") 
Full_risk <- predict(Full_cox, type = "risk")  
Lasso_risk <- predict(Lasso_cox, type = "risk")  
# 计算ROC曲线  
Univariate_roc <- roc(groupedclin$V66, Univariate_risk) 
Full_roc <- roc(groupedclin$V66, Full_risk)  
Lasso_roc <- roc(groupedclin$V66, Lasso_risk)  
# 计算NRI（需要使用自定义calc_nri函数）  
calc_nri <- function(roc1, roc2) {  
  auc1 <- auc(roc1)  
  auc2 <- auc(roc2)  
  nri <- auc2 - auc1  
  return(nri)  
}  
# 计算NRI  
NRI1 <- calc_nri(Full_roc, Lasso_roc)  
print(NRI1)
NRI2 <- calc_nri(Univariate_roc, Lasso_roc)  
print(NRI2)
####--计算IDI（NRI）并作图--####
# 加载survIDINRI包
library(survIDINRI)
##--结局变量转换--##
groupedclin$V66 <- as.numeric(groupedclin$V66==1)
indat <- groupedclin[,c("V68","V66")]
x_new <- groupedclin[,c("V43","V44","V57","V64","V73","V82","V91","V94")] #Lasso_cox
x_Univariate <- groupedclin[,c("V11","V25","V31","V43","V44","V64","V82","V91","V92","V94")] #Univariate_cox
x_Full <- groupedclin[,c("V11","V30","V31","V37","V43","V56","V64","V77","V78","V82","V91","V94")] #Full_cox
x_AJCC <- groupedclin[,c("V45","V46","V47","V48")] #AJCC_cox
##--Lasso_cox vs. Univariate_cox--##
x1 <- IDI.INF(indat,covs0 = as.matrix(x_Univariate),
              covs1 = as.matrix(x_new), 
              t0=12,  #注意选取Cutoff值，可以调整
              npert=100)
IDI.INF.OUT(x1)
IDI.INF.GRAPH(x1) ## M1为IDI值，M2为NRI值，M3表示中位数差异及相应的置信区间
##--Lasso_cox vs. Full_cox--##
x2 <- IDI.INF(indat,covs0 = as.matrix(x_Full),
              covs1 = as.matrix(x_new), 
              t0=12, 
              npert=100)
IDI.INF.OUT(x2)
IDI.INF.GRAPH(x2) ## M1为IDI值，M2为NRI值，M3表示中位数差异及相应的置信区间
##--Lasso_cox vs. AJCC_cox--##
x3 <- IDI.INF(indat,covs0 = as.matrix(x_AJCC),
              covs1 = as.matrix(x_new), 
              t0=12, 
              npert=100)
IDI.INF.OUT(x3)
IDI.INF.GRAPH(x3) ## M1为IDI值，M2为NRI值，M3表示中位数差异及相应的置信区间
##################--生存DCA曲线--##################
##--方法一：用ggDCA包，美观差一些--##
library(rmda)
library(ggDCA)
library(ggplot2)
library(rms)
library(caret)
# 使用Cox模型拟合数据  
Lasso_cox <- coxph(Surv(V68, V66)~ V43+V44+V57+V64+V73+V82+V91+V94, data= groupedclin)
Univariate_cox <- coxph(Surv(V68, V66)~ V11+V25+V31+V43+V44+V64+V82+V91+V92+V94, data= groupedclin)
Full_cox <- coxph(Surv(V68, V66)~ V11+V30+V31+V37+V43+V56+V64+V77+V78+V82+V91+V94, data= groupedclin)
AJCC_cox <- coxph(Surv(V68, V66)~ V45+V46+V47+V48, data= groupedclin)
# 组合多Cox模型数据，默认中位生存时间,线宽度lwd = 0.1,可调节，注意顺序
dca_cph <- dca(Lasso_cox,Univariate_cox,Full_cox,AJCC_cox, model.names = c("Lasso-cox","Univariate-cox", "Full-cox", "AJCC-cox"))
ggplot(dca_cph, lwd = 0.1)
# 单个模型多个时间点(生存时间的下四分位数、中位数、上四分位数)
times = round(quantile(groupedclin$V68, c(0.25, 0.5, 0.75)), 2)
dca_cph <- dca(Lasso_cox, model.names = "Lasso-cox", times = times)
ggplot(dca_cph)
# 多个模型单个时间点
dca_cph <- dca(Lasso_cox,Univariate_cox,Full_cox,AJCC_cox, model.names = c("Lasso-cox","Univariate-cox", "Full-cox", "AJCC-cox")
               , times = 12)
ggplot(dca_cph, lwd = 0.5)
# 多个模型多个时间点
dca_cph <- dca(Lasso_cox,Univariate_cox,Full_cox,AJCC_cox, model.names = c("Lasso-cox","Univariate-cox", "Full-cox", "AJCC-cox")
               , times = c(12, 36, 60))
ggplot(dca_cph)
##--方法二：用dcurves包--##
library(dcurves)
library(foreign)
library(survival)
library(dplyr)
library(tidyr)
##COX回归建模，建模无需因子转换##
Lasso_cox <- coxph(Surv(V68, V66)~ V43+V44+V57+V64+V73+V82+V91+V94, data= groupedclin)
Univariate_cox <- coxph(Surv(V68, V66)~ V11+V25+V31+V43+V44+V64+V82+V91+V92+V94, data= groupedclin)
Full_cox <- coxph(Surv(V68, V66)~ V11+V30+V31+V37+V43+V56+V64+V77+V78+V82+V91+V94, data= groupedclin)
AJCC_cox <- coxph(Surv(V68, V66)~ V45+V46+V47+V48, data= groupedclin)
##算出每个模型的相应time(12,36,60..)的生存率，注意模型对应##
groupedclin$pr_Lasso_cox = c(1- (summary(survfit(Full_cox, newdata=groupedclin), times= 28)$surv))
groupedclin$pr_Univariate_cox = c(1- (summary(survfit(Lasso_cox, newdata=groupedclin), times= 28)$surv))
groupedclin$pr_Full_cox = c(1- (summary(survfit(Univariate_cox, newdata=groupedclin), times= 28)$surv))
groupedclin$pr_AJCC_cox = c(1- (summary(survfit(AJCC_cox, newdata=groupedclin), times= 28)$surv))
##画出每个模型DCA曲线，注意time与前对应##
dca(Surv(V68, V66) ~ pr_Lasso_cox, 
    data = groupedclin,
    time = 30,
    thresholds = 1:50 / 100) %>%
  plot(smooth = T) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dca(Surv(V68, V66) ~ pr_Univariate_cox, 
    data = groupedclin,
    time = 30,
    thresholds = 1:50 / 100) %>%
  plot(smooth = T) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dca(Surv(V68, V66) ~ pr_Full_cox, 
    data = groupedclin,
    time = 30,
    thresholds = 1:50 / 100) %>%
  plot(smooth = T) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dca(Surv(V68, V66) ~ pr_AJCC_cox, 
    data = groupedclin,
    time = 30,
    thresholds = 1:50 / 100) %>%
  plot(smooth = T) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
##画出多个模型DCA曲线，注意time与前对应（12m-22,36m-25,60m-28）##
dca(Surv(V68, V66) ~ pr_Lasso_cox + pr_Univariate_cox + pr_Full_cox + pr_AJCC_cox,  
    data = groupedclin,
    time = 28,
    thresholds = 1:50 / 100) %>%
  plot(smooth = T) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
####--重要：画clinical impact曲线--####
cox_fit <- Lasso_cox
# Obtain predicted risk scores for each patient
risk_scores <- predict(cox_fit, newdata = groupedclin, type = "risk")
# Determine risk thresholds (e.g., using quantiles)
risk_thresholds <- quantile(risk_scores, probs = seq(0, 1, 0.1))
# Initialize vectors to store the results
total_high_risk <- numeric(length(risk_thresholds))
true_positives <- numeric(length(risk_thresholds))
# Loop through each risk threshold and calculate the total high risk and true positives
for (i in 1:length(risk_thresholds)) {
  threshold <- risk_thresholds[i]
  # Count the total number of patients deemed high risk
  total_high_risk[i] <- sum(risk_scores >= threshold)
  # Count the number of true positive cases among those deemed high risk
  true_positives[i] <- sum(risk_scores >= threshold & groupedclin$V66 == 1)
}
# Create a data frame with the results
clinical_impact <- data.frame(threshold = risk_thresholds,
                              total_high_risk = total_high_risk,
                              true_positives = true_positives)
# Calculate the false positives (total_high_risk - true_positives)
clinical_impact$false_positives <- clinical_impact$total_high_risk - clinical_impact$true_positives
# Calculate the cumulative true positive rate (cumulative sum of true positives divided by the total number of positive cases)
clinical_impact$cumulative_tpr <- cumsum(clinical_impact$true_positives) / sum(groupedclin$V66 == 1)
# Plot the clinical impact curve with confidence intervals and white background
ggplot(clinical_impact, aes(x = total_high_risk, y = true_positives)) +
  geom_line(color = "red") +
  geom_ribbon(aes(ymin = true_positives - 1.96 * sqrt(true_positives),
                  ymax = true_positives + 1.96 * sqrt(true_positives)),
              fill = "lightblue", alpha = 0.3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
  labs(x = "Number high risk", y = "Number high risk with ture event", title = "Clinical Impact Curve") +
  xlim(0, 1000) +
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"))

#######----检查多重共线性VIF----#######
library(car)
vif(Lasso_cox) #方差膨胀因子vif小于5，提示无共线性问题
#vif(Full_cox)
#vif(Univariate_cox)
#####----检查模型中变量间的相关性----#####
library(corrplot)
#需转换为连续变量
#for (i in names(train.data)[c(13,35,36,42)]) {train.data[,i] <- as.numeric(train.data[,i])} ##变量数值型转换
cr <- cor(groupedclin[,c(43,44,57,64,73,82,91,94)])
corrplot.mixed(cr)
corrplot(cr, method = "ellipse")
#####----重要：训练模型中timeROC曲线----#####
library(timeROC)
groupedclin$lp <- predict(Lasso_cox, newdata = groupedclin,type ="lp") # 以Lasso_cox为例
time_roc_train <- timeROC(
  T = groupedclin$V68,  #结局时间
  delta = groupedclin$V66,  #生存结局
  marker = groupedclin$lp,  #预测变量
  cause = 1,  #阳性结局指标值
  weighting="marginal",  #权重计算方法，默认marginal，采用KM计算删失分布
  times = c(12, 12*3, 12*5),
  ROC = TRUE,  #保存sensitivities 和 specificties值
  iid = TRUE  #iid = TRUE 才会保存置信区间，但是样本量大了后，耗时耗资源
)
#计算train-1年的AUC
train1y <- paste0("train 1 year AUC (95%CI) = ",
                  sprintf("%.3f",time_roc_train$AUC[1]) ," (",
                  sprintf("%.3f",confint(time_roc_train, level = 0.95)$CI_AUC[1,1]/100)," - ",
                  sprintf("%.3f",confint(time_roc_train, level = 0.95)$CI_AUC[1,2]/100),")")
#计算train-3年的AUC
train3y <- paste0("train 3 year AUC (95%CI) = ",
                  sprintf("%.3f",time_roc_train$AUC[2]) ," (",
                  sprintf("%.3f",confint(time_roc_train, level = 0.95)$CI_AUC[2,1]/100)," - ",
                  sprintf("%.3f",confint(time_roc_train, level = 0.95)$CI_AUC[2,2]/100),")")
#计算train-5年的AUC
train5y <- paste0("train 5 year AUC (95%CI) = ",
                  sprintf("%.3f",time_roc_train$AUC[3]) ," (",
                  sprintf("%.3f",confint(time_roc_train, level = 0.95)$CI_AUC[3,1]/100)," - ",
                  sprintf("%.3f",confint(time_roc_train, level = 0.95)$CI_AUC[3,2]/100),")")
#train AUC 多时间点合一
plot(title="ROC Curve",time_roc_train,col="DodgerBlue",time=12,lty=1,lwd = 2,family = "serif") #绘制ROC曲线
plot(time_roc_train,time=12*3,lty=1,lwd = 2,add=TRUE,col="LightSeaGreen",family = "serif")#add=TRUE指在前一条曲线基础上新增
plot(time_roc_train,time=12*5,lty=1,lwd = 2,add=TRUE,col="DarkOrange",family = "serif" )
legend("bottom", c("1-year AUC: ", "3-year AUC: ", "5-year AUC: "), 
       col=c("DodgerBlue", "LightSeaGreen", "DarkOrange"), lty=1, lwd=2, bty="n", bg="transparent")
##显示AUC值和95%CI##
legend("bottom",c(train1y,train3y,train5y),col=c("DodgerBlue","LightSeaGreen","DarkOrange"),lty=1,lwd=2)
##############################---建立Cox-Nomogram列线表模型---##############################################
library(foreign)
library(rms)
library(survival)
#rms包数据预处理--必须,数据必须无缺失##
#attach(groupedclin)
#首先告诉软件，哪些是你的分类变量。分类变量转为因子
for (i in names(groupedclin)[c(1:64)]) {groupedclin[,i] <- as.factor(groupedclin[,i])} # data[,"sex"] <- as.factor(data[,"sex"])
str(groupedclin)
dd <- datadist(groupedclin)
options(datadist = "dd")
##必须——检查响应变量的格式和列名##
str(groupedclin$V66) # 确认V66列包含事件情况
str(groupedclin$V68) # 确认V68列包含生存时间
names(groupedclin) # 查看数据集列名
# 更改响应变量的格式和列名
#groupedclin$V66 <- as.factor(groupedclin$V66) # 将事件情况转换为因子变量
#names(groupedclin)[c(64, 65)] <- c("Status", "Time") # 更改响应变量的列名
Lasso.nom.fit <- cph(Surv(V68, V66) ~ V43+V44+V57+V64+V73+V82+V91+V94, data = groupedclin, surv=T)
surv <- Survival(Lasso.nom.fit)
surv1 <- function(x)surv(1*12,lp=x)
surv2 <- function(x)surv(1*36,lp=x)
surv3 <- function(x)surv(1*60,lp=x)
nom <- nomogram(Lasso.nom.fit, fun = list(surv1,surv2,surv3),lp=T,
                funlabel = c("1-year OS","3-year OS","5-year OS"),
                maxscale = 100,
                fun.at = c('0.95','0.85','0.80','0.70','0.6','0.5','0.4','0.3','0.2','0.1'))
####################画出经典列线表######################

plot((nom),xfrac = 0.3)

####################动态Nomogram########################
#install.packages("DynNom")
library(DynNom)
DynNom(Lasso.nom.fit, groupedclin) ##注意，数据框的名字不要用data，会报错（遇到过这一问题）
#####################NomogramEx提取方程#################
#install.packages("C:\\Program Files\\R\\rms_6.7-0.tar.gz", repos = NULL, type = "source") #添加本地文件
#install.packages("nomogramEx")
library(nomogramEx)
nomogramEx(nomo=nom,np=3,digit=3)

##########--更换模型中的变量名画Nom图--############

## 更改响应变量的格式和列名
#names(groupedclin)[c(68, 66)] <- c("Time","Status") # 更改响应变量的列名
#names(groupedclin)[names(groupedclin) == "V43"] <- "Perineural invasion"
#names(groupedclin)[names(groupedclin) == "V44"] <- "Microvascular invasion"
#names(groupedclin)[names(groupedclin) == "V57"] <- "Postoperative hemorrhage"
#names(groupedclin)[names(groupedclin) == "V64"] <- "Tumor recurrence"
#names(groupedclin)[names(groupedclin) == "V73"] <- "AST (POD 1)"
#names(groupedclin)[names(groupedclin) == "V82"] <- "Intraoperative blood transfusion"
#names(groupedclin)[names(groupedclin) == "V91"] <- "Operative time"
#names(groupedclin)[names(groupedclin) == "V94"] <- "Age"
#write.csv(groupedclin, file = "Nom-name.csv") 
##变量全名读取，注意空格等特殊字符用.替代##
newgroupedclin <- read.csv(file = "Nom-name.CSV", header = TRUE)
library(foreign)
library(rms)
library(survival)
#首先告诉软件，哪些是你的分类变量。分类变量转为因子
for (i in names(newgroupedclin)[c(67:94)]) {newgroupedclin[,i] <- as.numeric(newgroupedclin[,i])} # data[,"sex"]<- as.factor(data[,"sex"])
for (i in names(newgroupedclin)[c(1:64)]) {newgroupedclin[,i] <- as.factor(newgroupedclin[,i])} 
newgroupedclin$Status <- as.factor(newgroupedclin$Status) # 将事件情况转换为因子变量
str(newgroupedclin)
dd <- datadist(newgroupedclin)
options(datadist = "dd")
##必须——检查响应变量的格式和列名##
newgroupedclin$Time <- as.numeric(newgroupedclin$Time)
newgroupedclin$Status <- as.numeric(newgroupedclin$Status) - 1
str(newgroupedclin$Status) # 确认V66列包含事件情况
str(newgroupedclin$Time) # 确认V68列包含生存时间
##使用Nom和DynNom画图##
library(DynNom)
Lasso.nom.fit <- cph(Surv(Time, Status) ~ Perineural.invasion +
                       Microvascular.invasion +
                       Postoperative.hemorrhage +
                       Tumor.recurrence +
                       AST..POD.1. +
                       Intraoperative.blood.transfusion +
                       Operative.time +
                       Age,
                     data = newgroupedclin, surv = TRUE)

surv <- Survival(Lasso.nom.fit)
surv1 <- function(x) surv(1*12, lp = x)
surv2 <- function(x) surv(1*36, lp = x)
surv3 <- function(x) surv(1*60, lp = x)
nom <- nomogram(Lasso.nom.fit, fun = list(surv1, surv2, surv3), lp = TRUE,
                funlabel = c("1-year OS", "3-year OS", "5-year OS"),
                maxscale = 100,
                fun.at = c('0.95', '0.85', '0.80', '0.70', '0.6', '0.5', '0.4', '0.3', '0.2', '0.1'))
## Plotting the classic nomogram ##
plot(nom, xfrac = 0.3)
### Dynamic Nomogram ###
DynNom(Lasso.nom.fit, newgroupedclin)
###NomogramEx提取方程###
library(nomogramEx)
nomogramEx(nomo=nom,np=3,digit=3)
#####################---Nomogram评价---#######################
###Lasso-Nomogram校准曲线（每个模型都可以做）#####
fit <- cph(Surv(V68, V66) ~ V43+V44+V57+V64+V73+V82+V91+V94, data = groupedclin, x = TRUE, y = TRUE, surv= TRUE)
##m要根据样本量来确定，由于校准曲线一般将所有的样本均分为5组(对应图中的5分节点)，而m代表每组样本量数，因此m*5等于或者近似等于样本量
##B代表最大再抽样的样本量，一般b=1000，足够使用
##u需要和之前模型中定义好的time.inc一致,以月为单位
##--1年生存曲线校准,u值可调--##
cal <- calibrate(fit, cmethod='KM',method="boot",u=26,m=300,B=100)
plot(cal,lwd=2,lty=1,
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)),
     xlim = c(0.4,1),ylim = c(0.4,1),
     xlab = "Predicted probability of 1-year",
     ylab = "Actual probability of 1-year",
     col = c(rgb(192,98,83,maxColorValue = 255)))
lines(cal[,c("mean.predicted","KM")],type="b",lwd=2,
      col = c(rgb(192,98,83,maxColorValue = 255)), pch=16) 
abline(0,1,lty=3,lwd=2,
       col= c(rgb(0,118,192,maxColorValue = 255))) 
##--3年生存曲线校准--##
cal <- calibrate(fit, cmethod='KM',method="boot",u=30,m=300,B=100)
plot(cal,lwd=2,lty=1,
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)),
     xlim = c(0.4,1),ylim = c(0.4,1),
     xlab = "Predicted probability of 3-year",
     ylab = "Actual probability of 3-year",
     col = c(rgb(192,98,83,maxColorValue = 255)))
lines(cal[,c("mean.predicted","KM")],type="b",lwd=2,
      col = c(rgb(192,98,83,maxColorValue = 255)), pch=16) 
abline(0,1,lty=3,lwd=2,
       col= c(rgb(0,118,192,maxColorValue = 255))) 
##--5年生存曲线校准--##
cal <- calibrate(fit, cmethod='KM',method="boot",u=34,m=300,B=100)
plot(cal,lwd=2,lty=1,
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)),
     xlim = c(0.4,1),ylim = c(0.4,1),
     xlab = "Predicted probability of 5-year",
     ylab = "Actual probability of 5-year",
     col = c(rgb(192,98,83,maxColorValue = 255)))
lines(cal[,c("mean.predicted","KM")],type="b",lwd=2,
      col = c(rgb(192,98,83,maxColorValue = 255)), pch=16) 
abline(0,1,lty=3,lwd=2,
       col= c(rgb(0,118,192,maxColorValue = 255))) 
##--多时间点拟合曲线--##
cal1 <- calibrate(fit, cmethod='KM',method="boot",u=26,m=300,B=100)
cal2 <- calibrate(fit, cmethod='KM',method="boot",u=30,m=300,B=100)
cal3 <- calibrate(fit, cmethod='KM',method="boot",u=34,m=300,B=100)
plot(cal1,lwd=2,lty=0,
     errbar.col=c("#2166AC"),
     xlim = c(0.4,1),ylim = c(0.4,1),
     xlab = "Predicted probability",
     ylab = "Actual probability",
     col = c("#2166AC"),
     cex.lab=1.2, cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal1[,c("mean.predicted","KM")],type="b",lwd=1,
      col = c("#2166AC"), pch=16) 
mtext("")
plot(cal2,lwd=2,lty=0,
     errbar.col=c("#B2182B"),
     xlim = c(0.4,1),ylim = c(0.4,1),
     xlab = "Predicted probability",
     ylab = "Actual probability",
     col = c("#B2182B"),
     add= T)
lines(cal2[,c("mean.predicted","KM")],type="b",lwd=1,
      col = c("#B2182B"), pch=16)
mtext("")
plot(cal3,lwd=2,lty=0,
     errbar.col=c("#E7B800"),
     xlim = c(0.4,1),ylim = c(0.4,1),
     xlab = "Predicted probability of 1-year",
     ylab = "Actual probability of 1-year",
     col = c("#E7B800"),
     add= T)
lines(cal2[,c("mean.predicted","KM")],type="b",lwd=1,
      col = c("#E7B800"), pch=16)
abline(0,1,lty=3,lwd=2,
       col= c("#224444"))  
legend("topleft",
       legend = c("1-year","3-year","5-year"),
       col = c("#2166AC","#B2182B","#E7B800"),
       lwd = 2,
       cex = 1.2,
       bty = "n")
###不同Nom模型生存曲线的比较#####
#install.packages("pec")
library(pec)
Lasso_nom <- cph(Surv(V68, V66) ~ V43+V44+V57+V64+V73+V82+V91+V94, data = groupedclin, surv=T)
Univariate_nom <- cph(Surv(V68, V66) ~ V11+V25+V31+V43+V44+V64+V82+V91+V92+V94, data = groupedclin, surv=T)
Full_nom <- cph(Surv(V68, V66) ~ V11+V30+V31+V37+V43+V56+V64+V77+V78+V82+V91+V94, data = groupedclin, surv=T)
AJCC_nom <- cph(Surv(V68, V66) ~ V45+V46+V47+V48, data = groupedclin, surv=T)
t <- c(12,36,60) ##定义t为1,3,5年生存率
Lasso.survprob <- predictSurvProb(Lasso_nom, newdata = groupedclin, times = t)
Univariate.survprob <- predictSurvProb(Univariate_nom, newdata = groupedclin, times = t)
Full.survprob <- predictSurvProb(Full_nom, newdata = groupedclin, times = t)
AJCC.survprob <- predictSurvProb(AJCC_nom, newdata = groupedclin, times = t)
head(Lasso.survprob)
head(Full.survprob)
head(Univariate.survprob)
##3个模型c-index随时间的变化曲线##
library(rms)
c_index <- cindex(list("Lasso_nom"= Lasso_nom,"Full_nom"= Full_nom,"Univariate_nom"= Univariate_nom,"AJCC_nom"= AJCC_nom),
                  formula = Surv(V68, V66) ~ V43+V44+V57+V64+V73+V82+V91+V94, data = groupedclin,
                  eval.times = seq(1,60,1))
plot(c_index,xlim = c(0,60))
##对3个模型c-index随时间的变化曲线进行交叉验证##
c_index <- cindex(list("Lasso_nom"= Lasso_nom,"Full_nom"= Full_nom,"Univariate_nom"= Univariate_nom,"AJCC_nom"= AJCC_nom),
                  formula = Surv(V68, V66) ~ V43+V44+V57+V64+V73+V82+V91+V94, data = groupedclin,
                  eval.times = seq(1,60,1),
                  splitMethod = "bootcv",
                  B=100)
plot(c_index,xlim = c(0,60))
####----预测校准度的检验(Hosmer-Lemeshow拟合优度检验,指定时间t，以月为单位)----####
library(rms)
##指定时间点t，以月为单位##
t <- c(60) 
Lasso.survprob <- predictSurvProb(Lasso_nom,newdata = groupedclin, times = t)
Univariate.survprob <- predictSurvProb(Univariate_nom, newdata = groupedclin, times = t)
Full.survprob <- predictSurvProb(Full_nom, newdata = groupedclin, times = t)
AJCC.survprob <- predictSurvProb(AJCC_nom, newdata = groupedclin, times = t)
#cox-nom模型生存状态预测比对#
library(ResourceSelection)
hoslem.test(groupedclin$V66,Lasso.survprob)
hoslem.test(groupedclin$V66,Univariate.survprob)
hoslem.test(groupedclin$V66,Full.survprob)
hoslem.test(groupedclin$V66,AJCC.survprob)
###########################--模型风险分层--################################
#devtools::install_github("selva86/InformationValue")
library(InformationValue)
library(pec)
library(rms)
library(survival) 
####-------真实值与预测值的比较------####
##--方法1：Lasso_cox模型--##
#Lasso_risk <- predict(Lasso_cox, type = "risk") 
##--预测矩阵--##
#groupedclin$V66 <- as.factor(groupedclin$V66) 
##--方法2：Lasso_nom模型，优势在于调整时间窗(可以画ROC曲线)--##
#指定时间t，以月为单位#
t <- c(60) 
Lasso.survprob <- predictSurvProb(Lasso_nom,newdata = groupedclin, times = t)
confusionMatrix(groupedclin$V66, Lasso.survprob)#混淆矩阵
optimalCutoff(groupedclin$V66, Lasso.survprob)#最优cutoff,模型在此时最优
misClassError(groupedclin$V66, Lasso.survprob)#误判率
library(pROC)
plot.roc(groupedclin$V66, Lasso.survprob,
         main = "ROC Curve", percent = TRUE,
         print.auc = TRUE,
         ci = TRUE, of = "thresholds",
         thresholds = "best",  #最优best_cutoff,模型在此时最优
         print.thres = "best",
         col = "#2E9FDF",
         ci.col = "#FF5733")
#################计算灵敏度、特异度和准确率###################
#--使用最佳阈值划分预测结果--# 
Lasso_risk <- predict(Lasso_cox, type = "risk")  
predicted_class <- ifelse(Lasso_risk > 0.6, 1, 0)  #以0.6为例
#计算混淆矩阵#  
confusionMatrix(groupedclin$V66, predicted_class)  
#计算灵敏度、特异度和准确率#  
sensitivity <- sum(predicted_class == 1 & groupedclin$V66 == 1) / sum(groupedclin$V66 == 1)  
specificity <- sum(predicted_class == 0 & groupedclin$V66 == 0) / sum(groupedclin$V66 == 0)  
accuracy <- sum(predicted_class == groupedclin$V66) / length(groupedclin$V66)  
#输出性能指标#  
cat("Sensitivity:", sensitivity, "\n")  
cat("Specificity:", specificity, "\n")  
cat("Accuracy:", accuracy, "\n")
##########--分层检验（用方法2确定截断值）--###########
#######先确定并标记optimal_cutoff=0.6#########
library(pROC)
roc <- roc(groupedclin$V66, Lasso.survprob)
# Define optimal cutoff (可以通过流行病学数据定义) #
optimal_cutoff <- 0.6
# Plot ROC curve and mark optimal cutoff  
plot(roc)
abline(v = optimal_cutoff, lty = 2, col = "red")
# Add annotation   
text(optimal_cutoff, 0.2, paste("Optimal cutoff =", optimal_cutoff),col = "red")
#####AUC计算并画图#####
groupedclin$V66 <- as.factor(groupedclin$V66)  
#-方法1：指定T-#
t <- 24
Lasso.survprob <- predictSurvProb(Lasso_nom,newdata = groupedclin, times = t)
groupedclin$prob <- Lasso.survprob
#-方法2：不用指定T-#
#Lasso.prob <- predict(Lasso_cox, type="risk")
#groupedclin$prob <- Lasso.prob
group1 <- groupedclin$prob <= 0.6 # 预后差组
group2 <- groupedclin$prob > 0.6 # 预后好组
roc1 <- roc(groupedclin$V66[group1], groupedclin$prob[group1])
roc2 <- roc(groupedclin$V66[group2], groupedclin$prob[group2])  
plot(roc1, col = "#FF5733")
lines(roc2, col = "#5DADE2")
auc1 <- auc(roc1)
auc2 <- auc(roc2)
##-做图-##
text(0.2, 0.375, paste("Risk stratification"), col = "black")
text(0.2, 0.3, paste("High risk to death, AUC =",round(auc1,2)), col = "#FF5733")
text(0.2, 0.2, paste("Low risk to death, AUC =",round(auc2,2)), col ="#5DADE2") 
##############作为分组变量做图#############
library(dplyr)  
# 向groupedclin数据集中添加predicted_class变量列  
new_groupedclin <- groupedclin %>%  
  mutate(predicted_class = ifelse(Lasso.survprob > 0.6, 1, 0))  
# 查看新的数据集  
head(new_groupedclin)
new_groupedclin <- new_groupedclin[,-95]
write.csv(new_groupedclin,file = "new_groupedclin.csv")
# 删除含有缺失值的行 
df <- read.csv(file = "new_groupedclin.csv",header = T)
# 使用surv_fit()函数拟合Cox生存模型，并根据风险等级进行分层  
fit <- survfit(Surv(V68, V66) ~ predicted_class, data = df)  
# 使用ggsurvplot()函数绘制分层Cox生存曲线  
ggsurvplot(fit,   
           pval = TRUE, conf.int = TRUE,  
           conf.int.style = "ribbon",  
           risk.table = TRUE,  
           ggtheme = theme_light() +  
             theme(panel.grid.major = element_blank(),  
                   panel.grid.minor = element_blank(),  
                   panel.background = element_blank()),  
           legend.labs = c("High", "Low"),  
           palette = c('#FF0000', '#2E9FDF'),
           xlim = c(0, 160))
###########结局占比#############
##--High risk group--## 
n1 <- length(groupedclin$V66[group1])
pos1 <- sum(groupedclin$V66[group1] == 1)
neg1 <- sum(groupedclin$V66[group1] == 0)
n1
pos1
neg1
##计算阳性比例##
prop_pos1 <- pos1 / n1   
prop_neg1 <- neg1 / n1
prop_pos1
prop_neg1
##--Low risk group--##
n2 <- length(groupedclin$V66[group2])
pos2 <- sum(groupedclin$V66[group2] == 1)
neg2 <- sum(groupedclin$V66[group2] == 0)
n2
pos2
neg2
##计算阳性比例##
prop_pos2 <- pos2 / n2   
prop_neg2 <- neg2 / n2
prop_pos2
prop_neg2
####也可依据风险分数（中位数/四分位）进行风险分层，详见New basic-02####
library(survival) 
library(survminer)  
Lasso_cox <- coxph(Surv(V68, V66) ~ V43+V44+V57+V64+V73+V82+V91+V94, data= groupedclin) 
# 获取模型的系数 #  
cox_coefficients <- Lasso_cox$coefficients  
# 计算风险分数risk_score #  
groupedclin$risk_score <- cox_coefficients[1] # 添加截距项  
for (i in 1:length(cox_coefficients[-1])) {  
  groupedclin$risk_score <- groupedclin$risk_score + cox_coefficients[i+1] * groupedclin[,paste0("V", i)]  
}  
# 输出风险分数 # 
print(groupedclin$risk_score)
####-计算风险分数的中位数-####
median_risk <- median(groupedclin$risk_score)  
# 创建风险等级变量进行2分层 #
groupedclin$risk_level <- ifelse(groupedclin$risk_score < median_risk, "低风险", "高风险")
####-使用surv_fit()函数拟合Cox生存模型，并根据风险等级进行分层-####
fit3 <- survfit(Surv(V68, V66) ~ risk_level, data = groupedclin)  
# 使用ggsurvplot()函数绘制分层Cox生存曲线  
ggsurvplot(fit3,    
           data = groupedclin, risk.table = TRUE,   
           title = "Cox生存曲线", xlab = "时间", ylab = "生存概率",   
           conf.int = FALSE, palette = c("#FF5733", "#5DADE2"),   
           risk.table.y.text = TRUE, risk.table.col = "risk_level")
#log rank检验
survdiff(Surv(V68, V66) ~ risk_level, data=groupedclin) 
#######################———生存曲线画图———##########################
library(survival)
library(survminer)
# 分组  
groupedclin$V42 <- as.factor(groupedclin$V42)  
group1 <- groupedclin$V42 == 1 # 壶腹癌组  
group2 <- groupedclin$V42 == 2 # 十二指肠癌组  
# 提取壶腹癌组的数据  
group1_data <- groupedclin[group1, c("V68", "V66", "V42")]  
group2_data <- groupedclin[group2, c("V68", "V66", "V42")]  
# 进行log rank检验（单因素）  
fit1 <- survfit(Surv(V68, V66) ~ 1, data = group1_data) # KM分析
fit2 <- survfit(Surv(V68, V66) ~ 1, data = group2_data) 
summary(fit1) #查看总体生存率#
summary(fit2)
surv_median(fit1) #计算中位生存期#
surv_median(fit2)
###---- 总体生存模型的比较（有无差异）----###
library(lmtest)
fit1_cox <- coxph(Surv(V68, V66) ~ 1, data= group1_data)
fit2_cox <- coxph(Surv(V68, V66) ~ 1, data= group2_data)
anova(fit1_cox,fit2_cox,test="Chisq") #模型有差异
#Likelihood ratio test，模型相似性检验
lrtest(fit1_cox,fit2_cox) 
#图形美化#
ggsurvplot(fit1, 
           pval = TRUE, conf.int =TRUE, #可信区间
           conf.int.style="ribbon", #风格
           risk.table = TRUE, #增加风险表
           ggtheme = theme_light() + #主题
             theme(panel.grid.major = element_blank(), #删除主网格线
                   panel.grid.minor = element_blank(), #删除次网格线
                   panel.background = element_blank()), #删除背景
           legend.labs = c("L"), #图例
           palette = c('#2E9FDF')) #配色

ggsurvplot(fit2, 
           pval = TRUE, conf.int =TRUE, #可信区间
           conf.int.style="ribbon", #风格
           risk.table = TRUE, #增加风险表
           ggtheme = theme_light() + #主题
             theme(panel.grid.major = element_blank(), #删除主网格线
                   panel.grid.minor = element_blank(), #删除次网格线
                   panel.background = element_blank()), #删除背景
           legend.labs = c("H"), #图例
           palette = c('#E7B800')) #配色
####——2张图组合起来展示（上下展示）--###
library(gridExtra)  
# 生成第一张生存曲线图  
ggsurvplot1 <- ggsurvplot(fit1, pval = TRUE, conf.int = TRUE, conf.int.style="ribbon", risk.table = TRUE, ggtheme = theme_light() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()), legend.labs = c("L"), palette = c('#2E9FDF'))  
# 生成第二张生存曲线图  
ggsurvplot2 <- ggsurvplot(fit2, pval = TRUE, conf.int = TRUE, conf.int.style="ribbon", risk.table = TRUE, ggtheme = theme_light() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()), legend.labs = c("H"), palette = c('#E7B800'))  
# 将两张图整合到一个坐标系中  
combined_plot <- grid.arrange(ggsurvplot1$plot, ggsurvplot2$plot, nrow = 2)  
# 显示整合后的图  
print(combined_plot)
###########——重要：2张图画到1张中--############
library(survminer)  
# 生成生存曲线数据  
#fit1 <- survfit(Surv(V68, V66) ~ 1, data = group1_data) # KM分析
#fit2 <- survfit(Surv(V68, V66) ~ 1, data = group2_data) 
# Create a named list of survival objects
surv_list <- list(L = fit1, H = fit2)
# Plot survival curves
ggsurvplot(surv_list, data = groupedclin, combine = TRUE,
           pval = TRUE, conf.int = TRUE, conf.int.style = "ribbon",
           risk.table = TRUE, ggtheme = theme_light() +
             theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank()),
           legend.labs = c("CPV", "DC"), palette = c('#2E9FDF', '#E7B800'),
           labels = c("CPV", "DC"))
###########——KM曲线（2分类）--############
#km检验
sfit <- survfit(Surv(V68, V66) ~ V39, data=groupedclin)   #KM分析
# 绘制生存曲线 # 
ggsurvplot(sfit, 
           pval = TRUE, conf.int =TRUE, #可信区间
           conf.int.style="ribbon", #风格
           risk.table = TRUE, #增加风险表
           ggtheme = theme_light() + #主题
             theme(panel.grid.major = element_blank(), #删除主网格线
                   panel.grid.minor = element_blank(), #删除次网格线
                   panel.background = element_blank()), #删除背景
           legend.labs = c("H","L"), #图例
           palette = c('#E7B800','#2E9FDF')) #配色
###########——KM曲线（3分类）--############
ggsurvplot(sfit, 
           pval = TRUE, conf.int = TRUE, #可信区间
           conf.int.style="ribbon", #风格
           risk.table = TRUE, #增加风险表
           ggtheme = theme_light() + #主题
             theme(panel.grid.major = element_blank(), #删除主网格线
                   panel.grid.minor = element_blank(), #删除次网格线
                   panel.background = element_blank()), #删除背景 + #主题
           legend.labs = c("High Risk","Medium Risk","Low Risk"), #图例
           palette = c('#2E9FDF','#E7B800','#FF0000')) #配色

################################--PERFECT--#######################################
