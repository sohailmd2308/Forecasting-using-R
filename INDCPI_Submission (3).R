library("fpp3")
library(tidyverse)
library(urca)
library(imputeTS)
library(forecast)
library(vars)

dataset = readr::read_csv("INDCPIALLMINMEI.csv")
#can values go negative or non zero and < 0?


finaldataset = dataset%>%mutate(DATE=yearmonth(DATE))%>%as_tsibble(index=DATE)
finaldataset %>% autoplot(INDCPIALLMINMEI)

finaldataset%>%filter(is.na(INDCPIALLMINMEI))

missingvalues=finaldataset%>%summarise(count = sum(is.na(INDCPIALLMINMEI)))
sum(missingvalues$count)

tail(finaldataset, n = 4)

traindataset = filter_index(finaldataset,"1960 Jan"~"2021 May")
traindataset %>% autoplot(INDCPIALLMINMEI)

traindataset%>%gg_season(INDCPIALLMINMEI)


#box-cox transform
lambdaIndCPI=traindataset%>%features(INDCPIALLMINMEI,features=guerrero)%>%
  pull(lambda_guerrero)
lambdaIndCPI #0.1389429

traindataset%>%autoplot(box_cox(INDCPIALLMINMEI,lambdaIndCPI))+ylab("Transformed data - India CPI")

# Box cox with 0
lambdaLOG=0
traindataset%>%autoplot(box_cox(INDCPIALLMINMEI,lambdaLOG))+
  ylab("Logarithmic Transformation")

# STL
traindataset%>%model(STL(INDCPIALLMINMEI))%>%
  components%>%autoplot()

#add new variable after transform      
traindataset=traindataset%>%
  mutate(transformCPI=box_cox(INDCPIALLMINMEI,lambdaIndCPI))

#checking for seasonal differencing
traindataset%>%features(transformCPI,unitroot_nsdiffs)
traindataset%>%features(transformCPI,unitroot_kpss)
traindataset%>%features(difference(transformCPI),unitroot_kpss)


# check for difference
transformTS = as.ts(dplyr::select(traindataset,transformCPI))
plot(transformTS)
summary(ur.df(transformTS,type="trend",lags=30,selectlags="AIC"))
summary(ur.df(transformTS,type="trend",lags=40,selectlags="AIC"))
summary(ur.df(transformTS,type="trend",lags=50,selectlags="AIC"))

summary(ur.df(transformTS,type="trend",lags=36,selectlags="AIC"))
#test statistic is greater than 1%,5% and 10% test sizes

#to check on seasonally differenced data
summary(ur.df(diff(transformTS,12),type="drift",lags=36,selectlags="AIC"))# lag is 36 at 36
summary(ur.df(diff(transformTS,12),type="drift",lags=50,selectlags="AIC"))# lag is 41 at 50
summary(ur.df(diff(transformTS,12),type="drift",lags=70,selectlags="AIC"))# lag is 65 at 70
summary(ur.df(diff(transformTS,12),type="drift",lags=90,selectlags="AIC"))# lag is 72 at 90
summary(ur.df(diff(transformTS,12),type="drift",lags=110,selectlags="AIC"))# lag is 72 at 110
summary(ur.df(diff(transformTS,12),type="drift",lags=200,selectlags="AIC"))# lag is 77 at 200

summary(ur.df(diff(transformTS,12),type="drift",lags=77,selectlags="AIC"))# lag is 77 at 200
#test statistic is greater than 1%,5% and 10% test sizes

# partial autocorrelation graph
traindataset%>%gg_tsdisplay(difference(transformCPI,12)%>%difference(),
                            plot_type="partial",lag_max=77)
#guess model
report(traindataset%>%model(ARIMA(transformCPI~0+pdq(1,1,2)+PDQ(3,1,1)))) #BIC = -4474.03


fit = traindataset%>%model(
  CPImodel1=ARIMA(transformCPI~0+pdq(3,1,3)+PDQ(1,1,1)), #best AIC
  CPImodel2=ARIMA(transformCPI~0+pdq(2,1,1)+PDQ(1,1,1)),
  CPImodel3=ARIMA(transformCPI~0+pdq(2,1,1)+PDQ(1,1,2)),
  CPImodel4=ARIMA(transformCPI~0+pdq(2,1,1)+PDQ(3,1,2)),
  CPImodel5=ARIMA(transformCPI~0+pdq(1,1,2)+PDQ(1,1,1)), #best BIC
  CPImodel6=ARIMA(transformCPI~0+pdq(2,1,2)+PDQ(1,1,1)),
  CPImodel7=ARIMA(transformCPI~0+pdq(2,1,2)+PDQ(1,1,2)),
  CPImodel8=ARIMA(transformCPI~0+pdq(1,1,2)+PDQ(3,1,1)),
  CPImodel9=ARIMA(transformCPI~0+pdq(3,1,1)+PDQ(2,1,1)),
  CPImodel10=ARIMA(transformCPI~0+pdq(2,1,1)+PDQ(3,1,2)))

glance(fit)$BIC
glance(fit)$AIC

#residuals
candidate = traindataset%>%model(ARIMA(transformCPI~0+pdq(1,1,2)+PDQ(1,1,1))) #final model = guess model on BIC
residuals(candidate)%>%gg_tsdisplay(.resid,plot_type="partial",
                                    lag_max=77) #most variations absorbed except 6,10,60 and few close ones

finalMODEL=traindataset%>%model(ARIMA(transformCPI~0+pdq(1,1,2)+PDQ(1,1,1)))

finalMODEL%>%gg_tsresiduals(lag=77)

#ljung box
finalMODEL=traindataset%>%model(ARIMA(box_cox(INDCPIALLMINMEI,lambdaIndCPI)~0+pdq(1,1,2)+PDQ(1,1,1)))
augment(finalMODEL)%>%features(.innov,ljung_box,lag=77,dof=5)

qchisq(0.95,72) #test statistic (84.8) < chisq (92.8), we are good to proceed

#forecast
finalMODEL%>%forecast(h=4)%>% autoplot(finaldataset%>%filter_index("2019 Jan"~"2021 Sep"), level = 95)

forecastSARIMA = finalMODEL%>%forecast(h=4)
forecastSARIMA%>%as.data.frame()%>%dplyr::select(DATE,.mean)

setsarima1 = filter_index(finaldataset,"2021 Jun"~"2021 Sep")%>%dplyr::select(DATE,INDCPIALLMINMEI)
setsarima2 = filter_index(finaldataset,"2021 Jun"~"2021 Sep")%>%dplyr::select(INDCPIALLMINMEI)
setsarima3 = forecastSARIMA%>%as.data.frame()%>%dplyr::select(.mean)
setsarima1$INDCPIALLMINMEI = paste(setsarima2$INDCPIALLMINMEI)
setsarima1$forecastedSARIMA = paste(setsarima3$.mean)
setsarima1


------------------------------------
#Neural networks
traindataset%>%gg_tsdisplay(difference(transformCPI,12)%>%difference(),lag_max=77,plot_type="partial")

ModelNN = traindataset%>%model(NNETAR(box_cox(INDCPIALLMINMEI,lambdaIndCPI)~AR(p=1,P=1)))
report(ModelNN)
ModelNN2 = traindataset%>%model(NNETAR(box_cox(INDCPIALLMINMEI,lambdaIndCPI)~AR(p=13,P=1)))
report(ModelNN2)
#ModelNN2 = traindataset%>%model(NNETAR(box_cox(INDCPIALLMINMEI,lambdaIndCPI)~AR(p=13,P=13)))
#report(ModelNN2)


ModelNN%>%forecast(h=4,bootstrap=TRUE)%>% 
  autoplot(finaldataset%>%filter_index("2019 Jan"~"2021 Sep"), level = 95)

ModelNN2%>%forecast(h=4,bootstrap=TRUE)%>% 
  autoplot(finaldataset%>%filter_index("2019 Jan"~"2021 Sep"), level = 95)

forecastNN = ModelNN%>%forecast(h=4,bootstrap=TRUE) #,times=500,bootstrap=TRUE
forecastNN%>%as.data.frame()%>%dplyr::select(DATE,.mean)

forecastNN2 = ModelNN2%>%forecast(h=4,bootstrap=TRUE) #,times=500,bootstrap=TRUE
forecastNN2%>%as.data.frame()%>%dplyr::select(DATE,.mean)

candidateNN = traindataset%>%model(NNETAR(box_cox(INDCPIALLMINMEI,lambdaIndCPI)~AR(p=13,P=1)))
residuals(candidateNN)%>%gg_tsdisplay(.resid,plot_type="partial",
                                    lag_max=77) 

#compare NN1 and NN2
finaldataset%>%filter_index("2019 Jan"~"2021 Sep")%>%autoplot(INDCPIALLMINMEI)+
  autolayer(ModelNN%>%forecast(h=4),colour="RED",level=NULL)+
  autolayer(ModelNN2%>%forecast(h=4),colour="BLUE",level=NULL)+
  ggtitle("NN p=1,P=1 in red; NN P=1,p=13 in Blue")+labs(y="CPI forecast comparisons among neural nets")

set1 = filter_index(finaldataset,"2021 Jun"~"2021 Sep")%>%dplyr::select(DATE,INDCPIALLMINMEI)
set2 = filter_index(finaldataset,"2021 Jun"~"2021 Sep")%>%dplyr::select(INDCPIALLMINMEI)
set3 = forecastNN%>%as.data.frame()%>%dplyr::select(.mean)
set4 = forecastNN2%>%as.data.frame()%>%dplyr::select(.mean)
set1$INDCPIALLMINMEI = paste(set2$INDCPIALLMINMEI)
set1$forecastNN = paste(set3$.mean)
set1$forecastNN2 = paste(set4$.mean)
set1


#compare NN2 and SARIMA
set1 = filter_index(finaldataset,"2021 Jun"~"2021 Sep")%>%dplyr::select(DATE,INDCPIALLMINMEI)
set2 = filter_index(finaldataset,"2021 Jun"~"2021 Sep")%>%dplyr::select(INDCPIALLMINMEI)
set3 = forecastSARIMA%>%as.data.frame()%>%dplyr::select(.mean)
set4 = forecastNN2%>%as.data.frame()%>%dplyr::select(.mean)
set1$INDCPIALLMINMEI = paste(set2$INDCPIALLMINMEI)
set1$forecastedSARIMA = paste(set3$.mean)
set1$forecastedNN2 = paste(set4$.mean)
set1


--------------------------- #final forecast
  
finallambdaIndCPI=finaldataset%>%features(INDCPIALLMINMEI,guerrero)%>%
  pull(lambda_guerrero)

finalSARIMA6step=finaldataset%>%
  model(ARIMA(box_cox(INDCPIALLMINMEI,finallambdaIndCPI)~0+pdq(1,1,2)+PDQ(1,1,1)))


finalSARIMA6step %>% forecast(h=6) %>%
  autoplot(filter_index(finaldataset,"2018 Jun"~"2021 Sep"),level=95)

forecasth6 = finalSARIMA6step%>%forecast(h=6)
forecasth6%>%as.data.frame()%>%dplyr::select(DATE,.mean)

setcombine = filter_index(finaldataset,"2021 Jun"~"2021 Sep")%>%dplyr::select(DATE,INDCPIALLMINMEI)
setcombine2 = filter_index(finaldataset,"2021 Jun"~"2021 Sep")%>%dplyr::select(INDCPIALLMINMEI)
setcombine$INDCPIALLMINMEI = paste(setcombine2$INDCPIALLMINMEI)
setcombine