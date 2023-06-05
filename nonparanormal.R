library(tidyverse)
library(GNAR)
library(QuantPsyc)
library(energy)
library(MVN)
library(huge)
library(lubridate)
library(igraph)
library(pcalg) #make direction
library(Rgraphviz)
library(ggplot2)
library(Hmisc)
# Note that original Rgraphviz package is no longer available..
# Instead we should install package Rgraphvis in BiocManager package.
# install.packages("BiocManager")
# BiocManager::install("Rgraphviz")


setwd('/Users/qsoon/Desktop/PGM')

source('util.R')

## preprocessing
obs_airkorea <- read.csv('data/observatory_airkorea.csv') # table1 has 690 obs
obs_kma <- read.csv('data/observatory_kma.csv',fileEncoding = "euc-kr") 

# "2018-01-01 ~ 2020-05-08"
t1 <- readRDS("data/forecast_dailymean.rds")
# t2 <- readRDS("./hourly1.rds") (Not used yet.)
train_table <- t1[[1]]$train_airkorea
train_table <- t1[[1]]$train_airkorea %>% add_column(year=year(train_table$date), month=month(train_table$date), day=day(train_table$date))

rec_boundary_seoul <- obs_airkorea %>% filter(lat >= 37.413294) %>% filter(lat<=37.715133) %>% filter(log >= 126.734086) %>% filter(log <= 127.269311)

district <- c("중구", "종로구", "용산구", "광진구", "성동구",
              "중랑구", "동대문구", "성북구", "도봉구", "은평구",
              "서대문구", "마포구", "강서구", "구로구", "영등포구",
              "동작구", "관악구", "강남구", "서초구", "송파구",
              "강동구", "금천구", "강북구", "양천구", "노원구")
seoul_table <- rec_boundary_seoul[1:25,] %>% add_column(district)

train_table_seoul <- train_table %>% filter(code %in% seoul_table$station)

## MAKE DISTRICT ##
district_table <- train_table_seoul %>% select(code,date,PM25)
district_table <- district_table %>% group_by(code,date) %>% summarise(daily_PM25 = mean(PM25,na.rm=TRUE))
district_table_summary <- district_table %>% spread(code,daily_PM25)
colnames(district_table_summary) <- c("date", district)
save(district_table_summary , file='data/seoul_districtPM25.RData')

## MAKE TIMESERIES ##
t11 <- readRDS("data/forecast_dailymean.rds")
# t2 <- readRDS("./hourly1.rds") Not used yet.
train_table2 <- t11[[1]]$train_kma
train_table2 <- t11[[1]]$train_kma %>% add_column(year=year(train_table2$date), month=month(train_table2$date), day=day(train_table2$date))

#  서울기상관측소(108) : 종로구송월길 52서울기상관측소  
seoul_kma <- train_table2[train_table2$loc==108,] # "2018-01-01 ~ 2020-05-08", No missing values.
seoul_kma_removed <- seoul_kma %>% select(c('year','month','day','temp','windspeed','winddirection','humidity'))

seoul_kma_removed_loc <- seoul_kma_removed[-1]

train_mean_seoul <- train_table_seoul %>% group_by(date) %>% 
  summarize(daily_SO2 = mean(SO2,na.rm=TRUE), 
            daily_CO = mean(CO,na.rm=TRUE),
            daily_O3 = mean(O3,na.rm=TRUE),
            daily_NO2 = mean(NO2,na.rm=TRUE),
            daily_PM10 = mean(PM10,na.rm=TRUE),
            daily_PM25 = mean(PM25,na.rm=TRUE))

result <- train_mean_seoul %>% add_column(year=year(train_mean_seoul$date), month=month(train_mean_seoul$date), day=day(train_mean_seoul$date))

seoul_daily <- result[-1] %>% mutate(seoul_kma_removed_loc)%>% relocate(c('year','month','day'))

save(seoul_daily, file='data/seoul_daily.RData')


############################################################
## TASK1: forecasting PM2.5 in Seoul from climate variables
############################################################
load('data/seoul_daily.RData')
seoul_daily # daily mean meteorology + air quality data in Seoul
seoul_daily <- seoul_daily[1:669,]
h=21 # how many days do we want to forecast?
data <- seoul_daily[1:(nrow(seoul_daily)-h),]

joint_vectors <- data[-c(1:3)]
N <- ncol(joint_vectors) # number of variables (nodes)
M <- nrow(joint_vectors)

## perform Multivariate normality test
hist.data.frame(joint_vectors)
mult.norm(joint_vectors)$mult.test # reject H0
mvnorm.etest(joint_vectors, R=100) # reject H0
(mvn(data = joint_vectors, mvnTest = "mardia"))$multivariateNormality # reject H0

## Use GGM
out_joint <- huge(as.matrix(joint_vectors), method = "glasso", nlambda=5,
                  lambda.min.ratio=0.4) # larger lambda.min.ratio, sparser graph 

out_joint_select <- huge.select(out_joint)

path_mat <- matrix(out_joint_select$refit,nrow=N,ncol=N)
g <- graph.adjacency(path_mat >.1, mode="undirected", diag=FALSE)
plot(g, layout=layout.fruchterman.reingold)
E(g)

## Use nonparanormal
joint_vectors_npn <- huge.npn(as.matrix(joint_vectors), npn.func = "truncation")
out_joint_npn <- huge(as.matrix(joint_vectors_npn), method = "glasso", nlambda=5,
                      lambda.min.ratio = 0.4)
out_joint_npn_select <- huge.select(out_joint_npn)

path_mat_npn <- matrix(out_joint_npn_select$refit,nrow=N,ncol=N)
g_npn <- graph.adjacency(path_mat_npn >.1, mode="undirected", diag=FALSE)
plot(g_npn, layout=layout.fruchterman.reingold)
E(g_npn)

## GNAR prediction
# GGM
PM25data <- list()
PM25data$A <- path_mat

PM25.GNARnet <- matrixtoGNAR(PM25data$A)

pred.gnar.PM25 <- forecast_narima0(vts = ts(joint_vectors), h = h, N = N, 
                                   net = PM25.GNARnet, max.alpha = 5, 
                                   max.beta = 3, globalalpha = TRUE, centering=TRUE, strd = TRUE)

rownames(pred.gnar.PM25) <- colnames(joint_vectors)
colnames(pred.gnar.PM25) <- 1:h

rmse(t(pred.gnar.PM25)[,"daily_PM25"], tail(seoul_daily,h)$daily_PM25)

# nonparanormal
PM25data_npn <- list()
PM25data_npn$A <- path_mat_npn

PM25_npn.GNARnet <- matrixtoGNAR(PM25data_npn$A)

pred.gnar.PM25_npn <- forecast_narima0(vts = ts(joint_vectors), h = h, N = N, 
                                       net = PM25_npn.GNARnet, max.alpha = 5, 
                                       max.beta = 3, globalalpha = TRUE, centering=TRUE, strd = TRUE)

rownames(pred.gnar.PM25_npn) <- colnames(joint_vectors)
colnames(pred.gnar.PM25_npn) <- 1:h

rmse(t(pred.gnar.PM25_npn)[,"daily_PM25"], tail(seoul_daily,h)$daily_PM25)


#############################
## TASK2: Find PM2.5 pathway 
#############################
load('data/seoul_districtPM25.RData')
district_table_summary # daily mean PM2.5 data in 25 districts of Seoul

joint_vectors2 <- (district_table_summary %>% drop_na())[-1] # remove NA
N2 <- ncol(joint_vectors2) # number of districts (nodes)
M2 <- nrow(joint_vectors2)

## perform Multivariate normality test
mult.norm(joint_vectors2)$mult.test # reject H0
mvnorm.etest(joint_vectors2, R=100) # reject H0
(mvn(data = joint_vectors2, mvnTest = "mardia"))$multivariateNormality # reject H0


## Preprocessing
load("data/seoul_daily_month.RData")
# 1. Remove all NA's
data <- data %>% na.omit()
## 1) ALL
firsthalf_data <- data %>% filter(month <= 5 | month == 12)
secondhalf_data <- data %>% filter(month >= 6 & month <= 11)

### (1) first half / second half ##
## Label ##
# 13 districts were selected as clusters from the full graph.
label <- c(1,2,3,6,7,13,16,18,21,22,23,24,25)

selected_fdata <- firsthalf_data[-(1:4)] %>% select(all_of(label)) # 12 ~ 5
selected_sdata <- secondhalf_data[-(1:4)] %>% select(all_of(label)) # 6 ~ 11

construct_CPDAG(selected_fdata,'NPN : 12 ~5',lambda=0.96,isnpn=TRUE,test_alpha=0.1)
construct_CPDAG(selected_sdata,'NPN : 6 ~ 11',lambda=0.95,isnpn=TRUE,test_alpha=0.1)

# Compare with normal cases

construct_CPDAG(selected_fdata,'Normal : 12 ~5',lambda=0.96,isnpn=FALSE,test_alpha=0.1)
construct_CPDAG(selected_sdata,'Normal : 6 ~ 11',lambda=0.95,isnpn=FALSE,test_alpha=0.1)

## 2) before COVID
# prev_COVID <- data %>% filter(year==2018 | (year==2019 & month < 10))
# prev_COVID
# spring_data <- prev_COVID %>% filter(month <= 5 & month >=3)
# summer_data <- prev_COVID %>% filter(month <= 8 & month >=6)
# fall_data <- prev_COVID %>% filter(month <= 11 & month >=9)
# winter_data <-prev_COVID %>% filter(month <= 2 | month ==12)
# 
# firsthalf_data <-prev_COVID %>% filter(month <= 5 | month == 12)
# secondhalf_data <-prev_COVID%>% filter(month >= 6 & month <= 11)
# 
# selected_fdata <- firsthalf_data[-(1:4)] %>% select(all_of(label)) # 12 ~ 5
# selected_sdata <- secondhalf_data[-(1:4)] %>% select(all_of(label)) # 6 ~ 11
# 
# construct_CPDAG(selected_fdata,'NPN : 12 ~5 (before COVID)',lambda=0.96,isnpn=TRUE,test_alpha=0.1)
# 
# construct_CPDAG(selected_sdata,'NPN : 6 ~ 11 (before COVID)',lambda=0.96,isnpn=TRUE,test_alpha=0.1)




# 2. fill NA

#load the original data again
load("data/seoul_daily_month.RData")
sum(is.na(data))
na_idx <- which(is.na(data), arr.ind=TRUE)

# fill NA values with row's mean.
for (i in 1:nrow(na_idx)){
  a = as.numeric(na_idx[i,1]) # row
  b = as.numeric(na_idx[i,2]) # column
  data[a,b] <- mean(as.numeric(data[a,-(1:4)]),na.rm=TRUE)
}
sum(is.na(data)) # 0

## Preprocessing(before COVID)
prev_COVID <- data %>% filter(year==2018 | (year==2019 & month < 10))

firsthalf_data <-prev_COVID %>% filter(month <= 5 | month == 12)
secondhalf_data <-prev_COVID%>% filter(month >= 6 & month <= 11)
selected_fdata <- firsthalf_data[-(1:4)] %>% select(all_of(label)) # 12 ~ 5
selected_sdata <- secondhalf_data[-(1:4)] %>% select(all_of(label)) # 6 ~ 11

# selected_spring <- spring_data[-(1:4)] %>% select(all_of(label)) # 12 ~ 5
# selected_summer <- summer_data[-(1:4)] %>% select(all_of(label)) # 12 ~ 5
# selected_fall <- fall_data[-(1:4)] %>% select(all_of(label)) # 12 ~ 5
# selected_winter <- winter_data[-(1:4)] %>% select(all_of(label)) # 12 ~ 5


construct_CPDAG(selected_fdata,'NPN : 12 ~5 (before COVID)',lambda=0.96,isnpn=TRUE,test_alpha=0.1)
construct_CPDAG(selected_sdata,'NPN : 6 ~ 11 (before COVID)',lambda=0.95,isnpn=TRUE,test_alpha=0.1)


construct_CPDAG(selected_fdata,'Normal : 12 ~5 (before COVID)',lambda=0.96,isnpn=TRUE,test_alpha=0.1)
construct_CPDAG(selected_sdata,'Normal : 6 ~ 11 (before COVID)',lambda=0.95,isnpn=TRUE,test_alpha=0.1)