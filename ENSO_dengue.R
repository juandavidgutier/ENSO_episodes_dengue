library(dplyr)
library(maptools)
library(sp)
library(spdep)
library(INLA)

#load data
url_path_data = "https://raw.githubusercontent.com/juandavidgutier/ENSO_episodes_dengue/master/dataset_dengue.csv"
colombia_dengue <- read.csv(url_path_data)
dim(colombia_dengue)
head(colombia_dengue)

colombia_dengue$consensus4 <- as.factor(colombia_dengue$consensus4)
colombia_dengue$consensus3 <- as.factor(colombia_dengue$consensus3)
colombia_dengue$consensus2 <- as.factor(colombia_dengue$consensus2)

#tertiles
colombia_dengue$q_SST12 <- dplyr::ntile(colombia_dengue$SST12, 3)  
colombia_dengue$q_SST3 <- dplyr::ntile(colombia_dengue$SST3, 3) 
colombia_dengue$q_SST34 <- dplyr::ntile(colombia_dengue$SST34, 3) 
colombia_dengue$q_SST4 <- dplyr::ntile(colombia_dengue$SST4, 3) 
colombia_dengue$q_ESOI <- dplyr::ntile(colombia_dengue$Equatorial_SOI, 3) 
colombia_dengue$q_SOI <- dplyr::ntile(colombia_dengue$SOI, 3) 
colombia_dengue$q_NATL <- dplyr::ntile(colombia_dengue$NATL, 3) 
colombia_dengue$q_SATL <- dplyr::ntile(colombia_dengue$SATL, 3) 
colombia_dengue$q_TROP <- dplyr::ntile(colombia_dengue$TROP, 3) 

library(tidyr)
colombia_dengue = colombia_dengue %>% drop_na()

#ids
colombia_dengue$ID.area = colombia_dengue$ID
colombia_dengue$ID.area1 = colombia_dengue$ID
colombia_dengue$ID.period = colombia_dengue$Period
colombia_dengue$ID.period1 = colombia_dengue$Period
colombia_dengue$ID.area.period = colombia_dengue$Code.DANE.period
colombia_dengue$ID.area.int <- colombia_dengue$ID.area
colombia_dengue$ID.period.int <- colombia_dengue$ID.period

#map
colombia <- readShapePoly("D:/msnm_dengue.shp")

#neighborhood matrix
nn_col <- poly2nb(colombia)
nb2INLA("colombia.graph", nn_col)
colombia.adj <- paste(getwd(),"/colombia.graph",sep="")

# INLA models
#consensus of four climate agencies
formula.all_c4 <-  Cases ~ 1 + consensus4 + q_SST12 + q_SST3 + q_SST34 + q_SST4 + q_ESOI + q_SOI + q_NATL + q_SATL + q_TROP +
                            f(ID.area, model="bym", graph=colombia.adj) +
                            f(ID.period, model="rw1") +
                            f(ID.area.period, model="iid") 

                            
model.all_c4 <- inla(formula.all_c4, family="nbinomial2", data=colombia_dengue, E=exp, verbose=TRUE,
                   control.predictor=list(compute=TRUE),
                   control.compute=list(dic=TRUE,cpo=TRUE))


result_model.all_c4 <-  as.matrix(model.all_c4$summary.fixed[2:4,1:5],3)
result_model.all_c4 [] <- exp(result_model.all_c4)
result_model.all_c4
model.all_c4$dic$dic


formula.lab_c4 <-  Cases_lab ~ 1 + consensus4 + q_SST12 + q_SST3 + q_SST34 + q_SST4 + q_ESOI + q_SOI + q_NATL + q_SATL + q_TROP +
                            f(ID.area, model="bym", graph=colombia.adj) +
                            f(ID.period, model="rw1") +
                            f(ID.area.period, model="iid") 


model.lab_c4 <- inla(formula.lab_c4, family="nbinomial2", data=colombia_dengue, E=exp_lab, verbose=TRUE,
                  control.predictor=list(compute=TRUE),
                  control.compute=list(dic=TRUE,cpo=TRUE))

result_model.lab_c4 <- as.matrix(model.lab_c4$summary.fixed[2:4,1:5],3)
result_model.lab_c4 [] <- exp(result_model.lab_c4)
result_model.lab_c4
model.lab_c4$dic$dic


#consensus of three climate agencies
formula.all_c3 <-  Cases ~ 1 + consensus3 + q_SST12 + q_SST3 + q_SST34 + q_SST4 + q_ESOI + q_SOI + q_NATL + q_SATL + q_TROP +
                              f(ID.area, model="bym", graph=colombia.adj) +
                              f(ID.period, model="rw1") +
                              f(ID.area.period, model="iid") 


model.all_c3 <- inla(formula.all_c3, family="nbinomial2", data=colombia_dengue, E=exp, verbose=TRUE,
                  control.predictor=list(compute=TRUE),
                  control.compute=list(dic=TRUE,cpo=TRUE))

result_model.all_c3 <-  as.matrix(model.all_c3$summary.fixed[2:4,1:5],3)
result_model.all_c3 [] <- exp(result_model.all_c3)
result_model.all_c3
model.all_c3$dic$dic


formula.lab_c3 <-  Cases_lab ~ 1 + consensus3 + q_SST12 + q_SST3 + q_SST34 + q_SST4 + q_ESOI + q_SOI + q_NATL + q_SATL + q_TROP +
                                f(ID.area, model="bym", graph=colombia.adj) +
                                f(ID.period, model="rw1") +
                                f(ID.area.period, model="iid") 


model.lab_c3 <- inla(formula.lab_c3, family="nbinomial2", data=colombia_dengue, E=exp_lab, verbose=TRUE,
                  control.predictor=list(compute=TRUE),
                  control.compute=list(dic=TRUE,cpo=TRUE))

result_model.lab_c3 <- as.matrix(model.lab_c3$summary.fixed[2:4,1:5],3)
result_model.lab_c3[] <- exp(result_model.lab_c3)
result_model.lab_c3
model.lab_c3$dic$dic


#consensus of two climate agencies
formula.all_c2 <-  Cases ~ 1 + consensus2 + q_SST12 + q_SST3 + q_SST34 + q_SST4 + q_ESOI + q_SOI + q_NATL + q_SATL + q_TROP +
                                f(ID.area, model="bym", graph=colombia.adj) +
                                f(ID.period, model="rw1") +
                                f(ID.area.period, model="iid") 


model.all_c2 <- inla(formula.all_c2, family="nbinomial2", data=colombia_dengue, E=exp, verbose=TRUE,
                     control.predictor=list(compute=TRUE),
                     control.compute=list(dic=TRUE,cpo=TRUE))

result_model.all_c2 <- as.matrix(model.all_c2$summary.fixed[2:3,1:5],3)
result_model.all_c2 [] <- exp(result_model.all_c2)
result_model.all_c2
model.all_c2$dic$dic


formula.lab_c2 <-  Cases_lab ~ 1 + consensus2 + q_SST12 + q_SST3 + q_SST34 + q_SST4 + q_ESOI + q_SOI + q_NATL + q_SATL + q_TROP +
                                  f(ID.area, model="bym", graph=colombia.adj) +
                                  f(ID.period, model="rw1") +
                                  f(ID.area.period, model="iid") 


model.lab_c2 <- inla(formula.lab_c2, family="nbinomial2", data=colombia_dengue, E=exp_lab, verbose=TRUE,
                     control.predictor=list(compute=TRUE),
                     control.compute=list(dic=TRUE,cpo=TRUE))

result_model.lab_c2 <- as.matrix(model.lab_c2$summary.fixed[2:3,1:5],3)
result_model.lab_c2 [] <- exp(result_model.lab_c2)
result_model.lab_c2
model.lab_c2$dic$dic







