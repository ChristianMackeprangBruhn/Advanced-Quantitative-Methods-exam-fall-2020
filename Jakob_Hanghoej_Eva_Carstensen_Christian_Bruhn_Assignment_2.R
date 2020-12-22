##############################################################################################################
###############################                ASSIGNMENT 2                  #################################
###############################          ADVANCED QUANTATIVE METHODS         #################################
###############################                                              #################################
###############################          BY EVA BJERRE CARSTENSEN,           #################################
###############################            CHRISTIAN BRUHN AND               #################################
###############################                JAKOB HANGHØJ                 #################################
##############################################################################################################

#rm(list=ls())

library(haven)
library(readxl)
library(tidyverse)
library(stargazer)
library(stringdist)
library(sf)
library(rgdal)
library(maps)
library(maptools)
library(ggsn)
library(scales)
library(gpclib)
library(mapproj)
library(stargazer)
library(lmtest)
library(lfe)
library(robustbase)
library(sandwich)
library(rdrobust)
library(geosphere)
library(fields)
library(data.table)
library(rdd)

setwd("/Users/jakobhanghoej/Desktop/Statskundskab/9. semester/Advanced Quantative Methods/01_Exam_AQM/R")

############################################################################################
######################              TABLE OF CONTENT                ########################
############################################################################################

# Line 41  : 2000-DATA, REPLICATION OF KRASNO & GREEN
# Line 480 : 2000 MAP
# Line 540 : 2008-DATA, DATA TIDYING AND MERGING
# Line 1223: 2008-DATA, ANALYSIS AND REGRESSION MODELS
# Line 1555: 2008-DATA, KEELE & TITIUNIK
# Line 1738: 2008 MAP

##############################################################################################################
###############################                  2000 DATA                   #################################
###############################        REPLICATION OF KRASNO & GREEN         #################################
###############################                                              #################################
##############################################################################################################

KG <- read_dta("KrasnoGreen_JOP_2008.dta")

############################################################################################
######################           CONSTRUCTING VARIABLES             ########################
############################################################################################
summary(KG)

#generates a variable concerning the contact rate which they define as (#polled*proportion reportedly contacted)/#polled
main$contrate <- (main$n0901*main$cont0901+main$n1001*main$cont1001+main$n1016*main$cont1016)/
  (main$n0901+main$n1001+main$n1016+1)

#generates a variable that is the log of total candidate visits
main$candvisits <- log(main$bush+main$cheney+main$gore+main$lieb+1)


KG <- KG %>% 
  mutate(contact_rate = (n0901*cont0901+n1001*cont1001+n1016*cont1016)/(n0901+n1001+n1016+1), #Contact rate which is (polled*proportion reportedly contacted)/polled
         candidate_visits = log(bush+cheney+gore+lieb+1), #log of total candidate visits
         allgrp_final = end_grp/1000, #GRP of all ads final three weeks, divided by 1000
         attackgrp_final = end_att/1000, #GRP of attack ads final three weeks, divided by 1000
         promogrp_final = end_pro/1000, #GRP of promotional ads final three weeks, divided by 1000 
         contrastgrp_final = end_cont/1000, #GRP of contrast ads final three weeks, divided by 1000 
         diffgrp_final = promogrp_final-attackgrp_final, #Difference between GRP for promotional and attack ads
  ) 

KG$state <- as.factor(as.character(KG$state)) #recode state variable to a factor variable

KG$prezg <- KG$prezg/100

############################################################################################
######################          F-TEST WITH FIXED EFFECTS           ########################
############################################################################################
fmodel1 <- lm(prezg~state - 1 , 
              data = KG)
fmodel2 <- lm(prezg ~ turn96 + turn92 + turn88 + senate + state - 1,
              data = KG)
fmodel3 <- lm(prezg ~ turn96 + turn92 + turn88 + senate + lncpop00 + seng + hseg + issg + othg + govg + state - 1,
              data = KG)
fmodel4 <- lm(prezg ~ turn96 + turn92 + turn88 + senate + lncpop00 + seng + hseg + issg + othg + govg + contact_rate + candidate_visits + state - 1,
              data = KG)

anova(fmodel1, fmodel2)$"Pr(>F)"[2]
anova(fmodel1, fmodel3)$"Pr(>F)"[2]
anova(fmodel1, fmodel4)$"Pr(>F)"[2]
anova(fmodel1, fmodel2, fmodel3, fmodel4)$"Pr(>F)"[2]

############################################################################################
######################            RESULTS FROM TABLE 2              ########################
############################################################################################
?bptest
?felm
model1 <- lm(turn00 ~ prezg,
             data = KG)
# Model 1
model1 <- lm(turn00 ~ prezg + senate,
             data = KG)
bptest(model1, turn00 ~ lncpop00+pq, 
       data = KG, studentize=FALSE) #Test for heteroskedasticity with Koenkers original test. But = TRUE is more robust?
# And why include population and non-voters?
bptest(model1) #This way its insignificant
coeftest(model1, vcov = vcovHC(model1, type="HC1")) # robust standard errors

summary(model1)

# Model 2
model2test <- lm(turn00 ~ prezg + senate + state - 1,
                 data = KG)
model2 <- felm(turn00 ~ prezg + senate | state,
               data = KG)
bptest(model2, turn00 ~ lncpop00 + pq,
       data = KG, studentize = FALSE)
bptest(model2) #Test for heteroskedasticity
coeftest(model2, vcov = vcovHC(model2, type="HC1")) # robust standard errors
?vcovHC

summary(model2)

# Model 3
model3test <- felm(turn00 ~ prezg + senate + turn96 + turn92 + turn88 | factor(state),
                   data = KG)

model3 <- felm(turn00 ~ prezg + senate + turn96 + turn92 + turn88 | state,
               data = KG)
bptest(model3, turn00 ~ lncpop00 + pq,
       data = KG, studentize = FALSE)
bptest(model3)
coeftest(model3, vcov = vcovHC(model3, type="HC1")) # robust standard errors
summary(model3)

# Model 4
model4 <- lm(turn00 ~ prezg + senate + turn96 + turn92 + turn88,
             data = KG)
bptest(model4, turn00 ~ lncpop00 + pq,
       data = KG, studentize = FALSE)
bptest(model4)
coeftest(model4, vcov = vcovHC(model4, type="HC1")) # robust standard errors
summary(model4)

# Model 5
model5 <- felm(turn00 ~ prezg + turn96 + turn92 + turn88 + seng + hseg + issg + othg + govg | state,
               data = KG)
bptest(model5, turn00 ~ lncpop00 + pq,
       data = KG, studentize = FALSE)
bptest(model5)
coeftest(model5, vcov = vcovHC(model5, type="HC1")) # robust standard errors
summary(model5)

# Model 5a
model5a <- lm(turn00 ~ prezg + turn96 + turn92 + turn88 + seng + hseg + issg + othg + govg,
              data = KG)
bptest(model5a, turn00 ~ lncpop00 + pq,
       data = KG, studentize = FALSE)
bptest(model5a)
coeftest(model5a, vcov = vcovHC(model5, type="HC1")) # robust standard errors
summary(model5a)

# Model 6                 
model6 <- felm(turn00 ~ prezg + turn96 + turn92 + turn88 + seng + hseg + issg + othg + govg + contact_rate + candidate_visits | state,
               data = KG)

bptest(model6, turn00 ~ lncpop00 + pq,
       data = KG, studentize = FALSE)
bptest(model6)
coeftest(model6, vcov = vcovHC(model6, type="HC1")) # robust standard errors
summary(model6)

# Model 7
model7 <- lm(turn00 ~ prezg + turn96 + turn92 + turn88 + seng + hseg + issg + othg + govg + contact_rate + candidate_visits,
             data = KG)
bptest(model7, turn00 ~ lncpop00 + pq,
       data = KG, studentize = FALSE)
bptest(model7)
coeftest(model7, vcov = vcovHC(model7, type="HC1")) # robust standard errors
summary(model7)

# Table 2
a <- c(coef(model1)[2], coef(model2)[1], coef(model3)[1], coef(model4)[2], coef(model5)[1], coef(model6)[1], coef(model7)[2])*100 # #100 because it is percentage points
table2 <- stargazer(model1, model2, model3, model4, model5, model6, model7, type = "text")
table2mini1 <- stargazer(model1, model2, model3, type = "text")
table2mini2 <- stargazer(model4, model5, model6, model7, type = "text")


stargazer(model1, model2, model3, model4, model5, model6, model7,
          type = "latex", # or "latex"
          model.names = FALSE, # removes OLS and felm
          dep.var.labels = c("Turnout '00"),
          covariate.labels = "Gross Rating Points",
          omit = c("senate", "turn96", "turn92", "turn88", "seng", "hseg", "issg", "othg", "govg", "contact_rate", "candidate_visits", "Constant"),
          add.lines = list(c("Midterm turnout",                               "Yes", "Yes", "Yes", "Yes", "Yes", "Yes", "Yes"),
                           c("State fixed effects",                           "", "Yes", "Yes", "No", "Yes", "Yes", "No"),
                           c("Past presidential turnout",                     "",  "", "Yes",  "Yes", "Yes", "Yes", "Yes"),
                           c("Control for other election ads",                "",  "",  "", "", "Yes", "Yes", "Yes"),
                           c("Control for candidate visits and contact rate", "",  "",  "", "", "", "Yes", "Yes")),
          star.cutoffs = c(0.05, 0.01, 0.001),
          digits = 2,
          keep.stat = c("adj.rsq", "n")) # Show just the adj. r-square and number of obs.)


############################################################################################
######################              INTERACTION TERMS               ########################
############################################################################################
## TEST FOR INTERACTION TERMS IN TABLE 2 FROM PAGE 256 ## <- OBS SKAL VI BRUGE FIXED EFFECTS?
# Interaction model 1
intact_model1 <- lm(turn00 ~ senate + bg*prezg,
                    data = KG)
summary(intact_model1)

# Interaction model 2
intact_model2 <- lm(turn00 ~ senate + bg*prezg + state - 1,
                    data = KG)
summary(intact_model2)

# Interaction model 3
intact_model3 <- lm(turn00 ~ senate + turn96 + turn92 + turn88 + bg * prezg + state - 1,
                    data = KG)
summary(intact_model3)

# Interaction model 4
intact_model4 <- lm(turn00 ~ senate + turn96 + turn92 + turn88 + bg * prezg,
                    data = KG)
summary(intact_model4)

# Interaction model 5
intact_model5 <- lm(turn00 ~ senate + turn96 + turn92 + turn88 + bg * prezg + seng + hseg + issg + othg + govg + state - 1,
                    data = KG)
summary(intact_model5)

# Interaction model 6
intact_model6 <- lm(turn00 ~ senate + turn96 + turn92 + turn88 + bg * prezg + seng + hseg + issg + othg + govg + contact_rate + candidate_visits + state - 1,
                    data = KG)
summary(intact_model6)

# Interaction model 7
intact_model7 <- lm(turn00 ~ senate + turn96 + turn92 + turn88 + bg * prezg + seng + hseg + issg + othg + govg + contact_rate + candidate_visits, 
                    data = KG)
summary(intact_model7)



############################################################################################
######################              ROBUSTNESS  TESTS               ########################
############################################################################################
##########
# TEST FOR MEDIA TONE FROM PAGE 257
##########
# Media model 1
media_model1 <- lm(turn00 ~ attackgrp_final + promogrp_final + contrastgrp_final + senate,
                   data = KG)
summary(media_model1)
simple_model1 <- lm(turn00 ~ senate,
                    data = KG)
summary(simple_model1)

anova(media_model1, model1)$"Pr(>F)"[2] # INS, test whether the effect of media tone is significantly different than the summarized GRP
anova(media_model1, simple_model1)$"Pr(>F)"[2] # SIG, test whether the effect of media tone is signicantly different than without GRPs

# Media model 2
media_model2 <- lm(turn00 ~ attackgrp_final + promogrp_final + contrastgrp_final + senate + state - 1,
                   data = KG)
summary(media_model2)
simple_model2 <- lm(turn00 ~ senate + state - 1,
                    data = KG)
summary(simple_model2)

###FIRST ANOVA IS SIGNIFICANT
anova(media_model2, model2)$"Pr(>F)"[2] #SIG, test whether the effect of media tone is significantly different than the summarized GRP
anova(media_model2, simple_model2)$"Pr(>F)"[2] #INS, test whether the effect of media tone is signicantly different from model without GRPs

# Media model 3
media_model3 <- lm(turn00 ~ attackgrp_final + promogrp_final + contrastgrp_final + senate + turn96 + turn92 + turn88 + state - 1,
                   data = KG)
summary(media_model3)
simple_model3 <- lm(turn00 ~ senate + turn96 + turn92 + turn88 + state - 1,
                    data = KG)
summary(simple_model3)

###FIRST ANOVA IS SIGNIFICANT
anova(media_model3, model3)$"Pr(>F)"[2] #SIG, test whether the effect of attacking media tone is significantly different than the summarized GRP
anova(media_model3, simple_model3)$"Pr(>F)"[2] #INS, test whether the effect of media tone is signicantly different than without GRPs

# Media model 4
media_model4 <- lm(turn00 ~ attackgrp_final + promogrp_final + contrastgrp_final + senate + turn96 + turn92 + turn88,
                   data = KG)
summary(media_model4)
simple_model4 <- lm(turn00 ~ senate + turn96 + turn92 + turn88,
                    data = KG)
summary(simple_model4)

anova(media_model4, model4)$"Pr(>F)"[2] #INS, test whether the effect of attacking media tone is significantly different than the summarized GRP
anova(media_model4, simple_model4)$"Pr(>F)"[2] #INS, test whether the effect of media tone is signicantly different than without GRPs

# Media model 5
media_model5 <- lm(turn00 ~ attackgrp_final + promogrp_final + contrastgrp_final + senate + turn96 + turn92 + turn88 + seng + hseg + issg + othg + govg + state - 1,
                   data = KG)
summary(media_model5)
simple_model5 <- lm(turn00 ~ senate + turn96 + turn92 + turn88 + seng + hseg + issg + othg + govg + state - 1,
                    data = KG)
summary(simple_model5)

anova(media_model5, model5)$"Pr(>F)"[2] #SIG, test whether the effect of attacking media tone is significantly different than the summarized GRP
anova(media_model5, simple_model5)$"Pr(>F)"[2] #INS, test whether the effect of media tone is signicantly different than without GRPs

# Media model 6
media_model6 <- lm(turn00 ~ attackgrp_final + promogrp_final + contrastgrp_final + senate + turn96 + turn92 + turn88 + seng + hseg + issg + othg + govg + contact_rate + candidate_visits + state - 1,
                   data = KG)
summary(media_model6)

simple_model6 <- lm(turn00 ~ senate + turn96 + turn92 + turn88 + seng + hseg + issg + othg + govg + contact_rate + candidate_visits + state - 1,
                    data = KG)
summary(simple_model6)

anova(media_model6, model6)$"Pr(>F)"[2] #SIG, test whether the effect of attacking media tone is significantly different than the summarized GRP
anova(media_model6, simple_model6)$"Pr(>F)"[2] #INS, test whether the effect of media tone is signicantly different than without GRPs

# Media model 7
media_model7 <- lm(turn00 ~ attackgrp_final + promogrp_final + contrastgrp_final + senate + turn96 + turn92 + turn88 + seng + hseg + issg + othg + govg + contact_rate + candidate_visits,
                   data = KG)
summary(media_model7)
simple_model7 <- lm(turn00 ~ senate + turn96 + turn92 + turn88 + seng + hseg + issg + othg + govg + contact_rate + candidate_visits,
                    data = KG)
summary(simple_model7)

anova(media_model7, model7)$"Pr(>F)"[2] #SIG, test whether the effect of attacking media tone is significantly different than the summarized GRP
anova(media_model7, simple_model7)$"Pr(>F)"[2] #INS, test whether the effect of media tone is signicantly different than without GRPs


############################################################
# TEST FOR DIFFERENCE BETWEEN PRO AND ATTACK FROM PAGE 257 #
############################################################
# Model 1
diff_model1 <- lm(turn00 ~ diffgrp_final + senate,
                  data = KG)
summary(diff_model1)

# Model 2
diff_model2 <- lm(turn00 ~ diffgrp_final + senate + state - 1,
                  data = KG)
summary(diff_model2)

# Model 3
diff_model3 <- lm(turn00 ~ diffgrp_final + senate + turn96 + turn92 + turn88 + state - 1,
                  data = KG)
summary(diff_model3)

# Model 4
diff_model4 <- lm(turn00 ~ diffgrp_final + senate + turn96 + turn92 + turn88,
                  data = KG)
summary(diff_model4)

# Model 5
diff_model5 <- lm(turn00 ~ diffgrp_final + turn96 + turn92 + turn88 + seng + hseg + issg + othg + govg + state - 1,
                  data = KG)
summary(diff_model5)

# Model 6                 
diff_model6 <- lm(turn00 ~ diffgrp_final + turn96 + turn92 + turn88 + seng + hseg + issg + othg + govg + contact_rate + candidate_visits + state - 1,
                  data = KG)
summary(diff_model6)

# Model 7
diff_model7 <- lm(turn00 ~ diffgrp_final + turn96 + turn92 + turn88 + seng + hseg + issg + othg + govg + contact_rate + candidate_visits,
                  data = KG)
summary(diff_model7)


############################################################################################
######################              REPLICATE TABLE 3               ########################
############################################################################################
## Greg have told us not to do Table 3.


############################################################################################
######################                  ADDING ADS                  ########################
############################################################################################

no_of_ads <- c(0.1, 0.01,1.82,2.33,1.42,0.5,1.44,0.97,1.23,0.49,4.03,0.07,0.37,0.02,1.59,2.78,1.82,2.87,2.87,2.08,
               0.1,0.0,1.59,2.34,0.37,0.01,2.04,0.85,2.15,0.85,0.68,0.94,0.0,0.47,2.4,0.0,0.0,1.45,0.68,1.72,
               0.19,0.47,1.39,0.0,1.61,0.06,0.0,1.98,0.0,2.86,2.1,
               2.86,3.02,0.73,2.34,2.4,1.42,0.37,2.15,1.42,1.82,0.41,0.1,0.0,0.01,0.0,0.01,0.0,0.07,0.37,1.61,
               1.67,0.02,2.25,0.07,3.28,0.01,0.0,0.01,0.02,0.0,0.02,1.45,0.68,1.19,0.95,0.94,0.73,0.01,0.0,3.56,2.04,
               0.01,2.0,0.02,2.25,1.98,0.0,1.92,0.0,0.0,1.72,
               1.42,1.39,0.0,0.0,0.0,0.0,0.01,0.01,0.0,0.01,0.01,0.0,0.0,1.61,3.56,2.38,2.04,3.02,2.7,2,1.45,
               1.98,0.01,0.0,0.07,0.01)

KG <- KG %>%
  arrange(state)
KG$no_ads_1000s <- no_of_ads
KG$no_ads <- no_of_ads * 1000


# Model 1
model1test <- lm(turn00 ~ no_ads + senate,
                 data = KG)
bptest(model1, turn00 ~ lncpop00+pq, 
       data = KG, studentize=FALSE) #Test for heteroskedasticity with Koenkers original test. But = TRUE is more robust?
# And why include population and non-voters?
bptest(model1) #This way its insignificant
coeftest(model1, vcov = vcovHC(model1, type="HC1")) # robust standard errors

summary(model1test)

# Model 2
model2test <- lm(turn00 ~ no_ads + senate + state - 1,
                 data = KG)
model2 <- felm(turn00 ~ prezg + senate | state,
               data = KG)
bptest(model2, turn00 ~ lncpop00 + pq,
       data = KG, studentize = FALSE)
bptest(model2) #Test for heteroskedasticity
coeftest(model2, vcov = vcovHC(model2, type="HC1")) # robust standard errors
?vcovHC

summary(model2test)

# Model 3
model3test <- felm(turn00 ~ no_ads + senate + turn96 + turn92 + turn88 | factor(state),
                   data = KG)

model3 <- felm(turn00 ~ prezg + senate + turn96 + turn92 + turn88 | state,
               data = KG)
bptest(model3, turn00 ~ lncpop00 + pq,
       data = KG, studentize = FALSE)
bptest(model3)
coeftest(model3, vcov = vcovHC(model3, type="HC1")) # robust standard errors
summary(model3)

# Model 4
model4 <- lm(turn00 ~ prezg + senate + turn96 + turn92 + turn88,
             data = KG)
bptest(model4, turn00 ~ lncpop00 + pq,
       data = KG, studentize = FALSE)
bptest(model4)
coeftest(model4, vcov = vcovHC(model4, type="HC1")) # robust standard errors
summary(model4)

# Model 5
model5 <- felm(turn00 ~ prezg + turn96 + turn92 + turn88 + seng + hseg + issg + othg + govg | state,
               data = KG)
bptest(model5, turn00 ~ lncpop00 + pq,
       data = KG, studentize = FALSE)
bptest(model5)
coeftest(model5, vcov = vcovHC(model5, type="HC1")) # robust standard errors
summary(model5)

# Model 5a
model5a <- lm(turn00 ~ prezg + turn96 + turn92 + turn88 + seng + hseg + issg + othg + govg,
              data = KG)
bptest(model5a, turn00 ~ lncpop00 + pq,
       data = KG, studentize = FALSE)
bptest(model5a)
coeftest(model5a, vcov = vcovHC(model5, type="HC1")) # robust standard errors
summary(model5a)

# Model 6                 
model6 <- felm(turn00 ~ prezg + turn96 + turn92 + turn88 + seng + hseg + issg + othg + govg + contact_rate + candidate_visits | state,
               data = KG)

bptest(model6, turn00 ~ lncpop00 + pq,
       data = KG, studentize = FALSE)
bptest(model6)
coeftest(model6, vcov = vcovHC(model6, type="HC1")) # robust standard errors
summary(model6)

# Model 7
model7 <- lm(turn00 ~ prezg + turn96 + turn92 + turn88 + seng + hseg + issg + othg + govg + contact_rate + candidate_visits,
             data = KG)
bptest(model7, turn00 ~ lncpop00 + pq,
       data = KG, studentize = FALSE)
bptest(model7)
coeftest(model7, vcov = vcovHC(model7, type="HC1")) # robust standard errors
summary(model7)


##############################################################################################################
###############################                  2000 DATA                   #################################
###############################                     MAP                      #################################
###############################                                              #################################
##############################################################################################################

######################################################
#          REPLICATE KRASNO AND GREEN'S MAP          #
######################################################
#Load Krasno & Green (2008) data
KG <- read_dta("KrasnoGreen_JOP_2008.dta")

##############################
#          RECODING          #
##############################
#Excel file for DMA codes for Krasno & Green's DMAs
dmacodes <- read_xlsx("/Users/christianmackeprangbruhn/Documents/DMACODES.xlsx")

dma_states <- KG %>%
  select(state, dma, end_grp)

#Merge Krasno & Green data with DMA codes. We do this in order to merge the data with the shape-file data
dma_ads <- merge(dmacode, dma_states, by="dma")

#####################################
#          READ SHAPE FLES          #
#####################################
gpclibPermit()

dma <- readOGR("DMAs.shp")

dma.df <- fortify(dma, region = "DMA")
dma.df <- rename(dma.df, DMA = id)

####################################################
#          MERGE SHAPE FILE WITH K&G DATA          #
####################################################
#We need to make a variable numeric
class(dma_ads$code)
class(dma.df$DMA)
dma.df <- transform(dma.df, DMA = as.numeric(DMA))

#Add values together
dma.df <- dma.df %>% 
  left_join(dma_ads, by = c("DMA" = "code")) #der går noget galt; vi mister observationer

###############################
#          PLOT DMAS          #
###############################
states_map <- map_data("state")

ggplot(dma.df, aes(x=long, y=lat, group=group, fill=end_grp)) + 
  geom_polygon(color="#666666", size=.5) +
  scale_fill_continuous(legend_title, low = "lightskyblue", high="deepskyblue4", na.value="white", limits=c(0,14000), breaks=c(2500,5000,7500,10000,12500)) +
  geom_polygon(data = states_map, aes(x = long, y = lat, group = group), fill=NA, color="black") + # states' borders
  coord_map() +
  theme_void() +
  ggsave("map_states_2000.png")


##############################################################################################################
###############################                  2008 DATA                   #################################
###############################           DATA TIDYING AND MERGING           #################################
###############################                                              #################################
##############################################################################################################

############################################################################################
######################                LOADING DATA                  ########################
############################################################################################
#voting data per county
load("countypres_2000-2016.RData")
countyvote_raw <- x
rm(x)

#censusdata
population <- read_csv("co-est00int-agesex-5yr.csv")

#ads data
wiscads_raw <- read_dta("WiscAds2008_Presidential.dta")

wiscads_GSH_raw <- read_dta("WiscAds2008_GSHData.dta")

#candidate appearance data
appearance_raw <- read_csv("appearances.csv")

#dmaindex
load("county_dma.RData")
dmaindex_raw <- x
rm(x)

DMA <- read_csv("DMA-zip.csv")

DMA <- DMA %>%
  distinct(FIPS, .keep_all = TRUE) %>% 
  rename(dmacode = "DMA CODE") %>% 
  mutate(FIPS = as.numeric(FIPS))

#senateelection
electiondata_raw <- read_dta("nanda_voting_county_2004-2018_01P.dta")
names(electiondata_raw)

############################################################################################
######################                DATA TIDYING                  ########################
############################################################################################
#####################################
#          SENATE ELECTION          #
#####################################
senateelection06 <- electiondata_raw %>% 
  filter(year == 2006) %>% 
  rename(FIPS = stcofips,
         votes06 = ballots_cast,
         pop06 = cvap) %>% 
  select(FIPS, votes06, pop06) %>% 
  mutate(FIPS = as.numeric(FIPS))

pop02 <- electiondata_raw %>% 
  filter(year == 2004) %>% 
  rename(FIPS = stcofips,
         pop02 = cvap) %>% 
  select(FIPS, pop02) %>% 
  mutate(FIPS = as.numeric(FIPS))

senateelection02 <- electiondata_raw %>% 
  filter(year == 2002) %>% 
  rename(FIPS = stcofips) %>% 
  mutate(votes02 = (sen_dem_votes + sen_rep_votes)*1.008, #ratio to independent candidates
         votes02 = as.integer(votes02)) %>%  
  select(FIPS, votes02) %>% 
  mutate(FIPS = as.numeric(FIPS))

senateelection02 <- left_join(senateelection02, pop02, by = "FIPS")

preselection08 <- electiondata_raw %>% 
  filter(year == 2008) %>% 
  rename(FIPS = stcofips) %>% 
  select(FIPS, ballots_cast, cvap) %>% 
  mutate(FIPS = as.numeric(FIPS))

#### Unfortunately, it wasn't possible to find senate data from 20 media markets on county level

#####################################
#      WISCONSIN PRES AD DATA       #
#####################################
class(wiscads_raw["Weeks_To"])
wiscads_raw$Weeks_To <- as.numeric(wiscads_raw$Weeks_To)

# only 3 weeks before election date
wiscads <- wiscads_raw %>% 
  filter(Weeks_To <= 3) %>% 
  select(market, STATE_1, date, Days_To, Weeks_To, EST_SPENDING, statdist, CAND_ID, party, sponsor, program, AD_TONE)

wiscads <- wiscads %>% filter(STATE_1 != "NATIONAL")

wiscads <- subset(wiscads, STATE_1 != "")

wisc_totalads <- wiscads %>% 
  mutate(contrastads = ifelse(AD_TONE == 1, 1, 0),
         promoads = ifelse(AD_TONE == 2, 1, 0),
         attackads = ifelse(AD_TONE == 3, 1, 0),
         totalads = contrastads + promoads + attackads)

wisc_totalads <- wisc_totalads %>% 
  group_by(market) %>% 
  summarize(contrastads = sum(contrastads),
            promoads = sum(promoads),
            attackads = sum(attackads),
            totalads = sum(totalads),
            state_po = unique(STATE_1))


#In order to merge between different data sources, we have to insure that the names are the same 
#first, we reduce difference in strings by changing all to capital letters
wisc_totalads$market <- str_to_upper(wisc_totalads$market)
#second, remove comma and whitespace
wisc_totalads$market <- str_trim(wisc_totalads$market)
wisc_totalads$state_po <- str_trim(wisc_totalads$state_po)
#third, remove spaces and everything after the first word
wisc_totalads$market <- gsub("[,]", "", wisc_totalads$market)
wisc_totalads$market <- gsub("[.]", "", wisc_totalads$market)
wisc_totalads$market <- gsub("[-]", " ", wisc_totalads$market)
wisc_totalads$market <- gsub("[/]", " ", wisc_totalads$market)
wisc_totalads$market <- gsub("? .*", "", wisc_totalads$market) #This deletes everything after the first space

#Some names are the same after removing everything after the first word, so we change these manually
wisc_totalads <- wisc_totalads %>% 
  mutate(market = case_when(totalads == 4 & contrastads == 4 & state_po == "CA" ~ "SAN FRANCISCO",
                            totalads == 1 & state_po == "CA" ~ "SAN DIEGO",
                            totalads == 2763 & state_po == "NC" ~ "GREENVILLE SC",
                            totalads == 3690 & state_po == "NC" ~ "GREENVILLE NC",
                            totalads == 2180 & state_po == "IA" ~ "SIOUX CITY",
                            totalads == 4 & market == "SIOUX" & state_po == "SD" ~ "SIOUX FALLS",
                            totalads == 3278 & state_po == "MO" ~ "ST LOUIS",
                            totalads == 1 & contrastads == 0 & state_po == "NY" ~ "ALBANY, NY",
                            totalads == 1 & contrastads == 0 & state_po == "SC" ~ "CHARLESTON, SC",
                            totalads == 2143 & contrastads == 135 & state_po == "MO" ~ "COLUMBIA, MO",
                            totalads == 4 & contrastads == 0 & state_po == "MS" ~ "COLUMBUS, AL",
                            totalads == 4529 & contrastads == 734 & state_po == "OH" ~ "COLUMBUS, OH",
                            totalads == 2578 & contrastads == 386 & state_po == "FL" ~ "FT MYERS",
                            totalads == 3050 & contrastads == 631 & state_po == "IN" ~ "FT WAYNE",
                            totalads == 436 & contrastads == 1 & state_po == "MI" ~ "GRAND RAPIDS",
                            totalads == 1 & contrastads == 1 & state_po == "TN"  ~ "JACKSON, TN",
                            totalads == 10 & contrastads == 5 & state_po == "LA" ~ "LAFAYETTE, LA",
                            totalads == 1368 & state_po == "ME" & contrastads == 1 ~ "PORTLAND, ME",
                            totalads == 4 & state_po == "TX" & market == "SAN" & attackads == 4 ~ "SAN ANTONIO",
                            totalads == 2 & attackads == 2 & market == "WHICHITA" & state_po == "TX" ~ "WICHITA FALLS",
                            TRUE ~ as.character((market))))

wisc_totalads <- wisc_totalads %>% 
  group_by(market, state_po) %>% 
  summarize(totalads = sum(totalads),
            contrastads = sum(contrastads),
            attackads = sum(attackads),
            promoads = sum(promoads)) %>% 
  ungroup()


#####################################
#            DMAINDEX               #
#####################################
dmaindex <- dmaindex_raw %>% 
  select(-CNTYTVHH) %>%
  rename(market = DMA)

#Create FIPS code
dmaindex <- dmaindex %>%
  mutate(STATEFIPS = STATEFP,
         STATEFIPS = ifelse(STATEFIPS < 10 & CNTYFP < 10, STATEFIPS*100, STATEFIPS),
         STATEFIPS = ifelse(STATEFIPS >= 10 & STATEFIPS < 100 & CNTYFP < 100 & CNTYFP >= 10, STATEFIPS*10, STATEFIPS),
         STATEFIPS = ifelse(STATEFIPS < 10 & CNTYFP >= 10 & CNTYFP < 100, STATEFIPS*10, STATEFIPS),
         STATEFIPS = ifelse(STATEFIPS >= 10 & STATEFIPS < 100 & CNTYFP < 10, STATEFIPS*10, STATEFIPS))

dmaindex <- dmaindex %>% 
  unite(FIPS, STATEFIPS, CNTYFP, sep = "") %>% 
  select(-STATEFP)

dmaindex$FIPS <- as.numeric(dmaindex$FIPS)

#Reducing differences between dmas
dmaindex$market <- str_trim(dmaindex$market)
dmaindex$state_po <- str_trim(dmaindex$STATE)
dmaindex$market <- gsub("[()]", " ", dmaindex$market)
dmaindex$market <- gsub("[(,)]", " ", dmaindex$market)
dmaindex$market <- gsub("[(.)]", " ", dmaindex$market)
dmaindex$market <- gsub("[/]", " ", dmaindex$market)
dmaindex$market <- gsub("[-]", " ", dmaindex$market)

dmaindex$market <- gsub("? .*", "", dmaindex$market) #This deletes everything after the first space


# finding states with only one DMA
one_dma_states <- dmaindex %>%
  select(market, state_po)

one_dma_states <- one_dma_states[!duplicated(one_dma_states[c("market","state_po")]),]

one_dma_states <- one_dma_states %>% 
  group_by(state_po) %>% 
  summarize(count = n())


#####################################
#       WISCONSIN GSH AD DATA       #
#####################################
wiscads_GSH_raw$Weeks_To <- as.numeric(wiscads_GSH_raw$Weeks_To)

# only 3 weeks before election date
wiscads_GSH <- wiscads_GSH_raw %>% 
  filter(Weeks_To <= 3) %>% 
  select(market, CAND_ID, AD_TONE, office)

wiscads_GSH$CAND_ID <- gsub("[/]", " ", wiscads_GSH$CAND_ID)
wiscads_GSH$state_po <- gsub("? .*", "", wiscads_GSH$CAND_ID) #This deletes everything after the first space

wiscads_GSH <- wiscads_GSH %>% 
  group_by(state_po) %>% 
  filter(n_distinct(market) > 1) %>%  #removes states with only one DMA
  ungroup() %>% 
  select(-CAND_ID)

#Create columns for different ad tones
wiscads_totalGSH <- wiscads_GSH %>% 
  mutate(contrastads = ifelse(AD_TONE == 1, 1, 0),
         promoads = ifelse(AD_TONE == 2, 1, 0),
         attackads = ifelse(AD_TONE == 3, 1, 0),
         totalads = contrastads + promoads + attackads)

#Sum by market
wiscads_totalGSH <- wiscads_totalGSH %>% 
  group_by(market, office) %>% 
  summarize(contrastads = sum(contrastads),
            promoads = sum(promoads),
            attackads = sum(attackads),
            totalads = sum(totalads),
            state_po = unique(state_po)) %>% 
  ungroup()


#first we try to reduce difference in strings by changing all to capital letters
wiscads_totalGSH$market <- str_to_upper(wiscads_totalGSH$market)
#then remove comma and whitespace
wiscads_totalGSH$market <- str_trim(wiscads_totalGSH$market)
wiscads_totalGSH$state_po <- str_trim(wiscads_totalGSH$state_po)
#then remove spaces and everything after the first word
wiscads_totalGSH$market <- gsub("[,]", "", wiscads_totalGSH$market)
wiscads_totalGSH$market <- gsub("[.]", "", wiscads_totalGSH$market)
wiscads_totalGSH$market <- gsub("[-]", " ", wiscads_totalGSH$market)
wiscads_totalGSH$market <- gsub("[/]", " ", wiscads_totalGSH$market)
wiscads_totalGSH$market <- gsub("? .*", "", wiscads_totalGSH$market) #This deletes everything after the first space

#Some names are the same after removing everything after the first word, so we change these manually
wiscads_totalGSH <- wiscads_totalGSH %>% 
  mutate(market = case_when(totalads == 412 & contrastads == 4 & attackads == 404 & state_po == "CA" ~ "SAN FRANCISCO",
                            totalads == 724 & contrastads == 261 & state_po == "CA" ~ "SAN DIEGO",
                            totalads == 7207 & contrastads == 646 & state_po == "SC" ~ "GREENVILLE SC",
                            totalads == 7207 & contrastads == 646 & state_po == "NC" ~ "GREENVILLE NC",
                            totalads == 111 & contrastads == 20 & state_po == "SC" ~ "GREENVILLE SC",
                            totalads == 111 & contrastads == 20 & state_po == "NC" ~ "GREENVILLE NC",
                            totalads == 3294 & contrastads == 245 & state_po == "NC" ~ "GREENVILLE NC",
                            market == "SIOUX" & state_po == "IA" ~ "SIOUX CITY",
                            market == "SIOUX" & state_po == "NE" ~ "SIOUX CITY",
                            market == "SIOUX" & state_po == "MN" ~ "SIOUX FALLS",
                            market == "SIOUX" & state_po == "NJ" ~ "SIOUX FALLS",
                            market == "SIOUX" & state_po == "SD" ~ "SIOUX FALLS",
                            market == "ST" & state_po == "MO" ~ "ST LOUIS",
                            TRUE ~ as.character((market))))

#Create columns for senate and house ads
wiscads_totalGSH <- wiscads_totalGSH %>% 
  filter(office != 1) %>% 
  mutate(senateads = ifelse(office == 2, totalads, 0),
         houseads = ifelse(office == 3, totalads, 0)) %>% 
  select(market, state_po, senateads, houseads)

#Some names are the same after removing everything after the first word, so we change these manually
wiscads_totalGSH <- wiscads_totalGSH %>% 
  mutate(market = case_when(senateads == 3702 & state_po == "NY" ~ "ALBANY, NY",
                            senateads == 4056 & state_po == "SC" ~ "CHARLESTON, SC",
                            senateads == 1645 & state_po == "MO" ~ "COLUMBIA, MO",
                            senateads == 397 & state_po == "MS" ~ "COLUMBUS, AL",
                            senateads == 4527 & state_po == "OH" ~ "COLUMBUS, OH",
                            senateads == 1641 & state_po == "FL" ~ "FT MYERS",
                            houseads == 1482 & state_po == "IN" ~ "FT WAYNE",
                            senateads == 1201 & state_po == "LA" ~ "LAFAYETTE, LA",
                            senateads == 2985 & state_po == "ME" ~ "PORTLAND, ME",
                            senateads == 1818 & state_po == "TX" ~ "SAN ANTONIO",
                            TRUE ~ as.character((market))))

#Sum by market
wiscads_totalGSH <- wiscads_totalGSH %>% 
  group_by(market, state_po) %>% 
  summarize(senateads = sum(senateads),
            houseads = sum(houseads))%>% 
  ungroup()

#####################################
#      CENSUS POPULATION 2008       #
#####################################
pop <- population %>% 
  select(`SUMLEV`:`AGEGRP`, POPESTIMATE2008) %>%
  filter(SEX == 0) %>%  #total
  rename(statecode = STATE,
         countycode = COUNTY,
         state = STNAME,
         county = CTYNAME,
         agegroup = AGEGRP)

pop <- pop %>% 
  group_by(statecode, countycode, state, county, agegroup) %>% 
  summarize(population = sum(POPESTIMATE2008)) %>% 
  filter(statecode != "02") # remove Alaska 

# find voter aged population
pop <- pop %>% 
  mutate(population = ifelse(agegroup == 4, round(population*0.4), population))

pop <- pop %>% 
  group_by(statecode, countycode, state, county) %>% 
  filter(agegroup >= 4) %>% 
  summarize(population = sum(population))

pop <- pop %>% 
  unite(FIPS, statecode:countycode, remove = TRUE, sep = "")

pop$FIPS <- as.numeric(pop$FIPS)


#####################################
#         COUNTY VOTING DATA        #
#####################################
############
#TURNOUT '08
############
countyvote08 <- countyvote_raw %>% 
  filter(year == "2008" & party != "NA") %>% 
  select(state_po:FIPS, totalvotes) %>% 
  distinct()

countyvote08 <- countyvote08 %>% 
  filter(totalvotes != 0)

# Alaska data is not by county so we delete it
countyvote08 <- countyvote08 %>% 
  filter(FIPS <= 2000 | FIPS > 3000)

#merge to turnout08
turnout08 <- left_join(countyvote08, pop, by = c("FIPS"))

turnout08 <- turnout08 %>% 
  select(-county.y) %>% 
  rename(county = county.x)

turnout08$FIPS <- as.numeric(as.character(turnout08$FIPS))

turnout08 <- turnout08 %>% 
  rename(pop08 = population,
         votes08 = totalvotes) %>% 
  arrange(FIPS)

#remove states with only one DMA
turnout08 <- left_join(turnout08, one_dma_states, by = "state_po", all = TRUE)

#remove states with only one DMA
turnout08 <- turnout08 %>% 
  filter(count != 1) %>% 
  select(FIPS, votes08, pop08)

############
#TURNOUT '04
############
countyvote04 <- countyvote_raw %>% 
  filter(year == "2004" & party != "NA") %>% 
  select(state_po:FIPS, totalvotes) %>% 
  distinct()

countyvote04 <- countyvote04 %>% 
  filter(totalvotes != 0)

# Alaska data is not by county so we delete it
countyvote04 <- countyvote04 %>% 
  filter(FIPS <= 2000 | FIPS > 3000) %>% 
  rename(votes04 = totalvotes)

pop04 <- population %>% 
  select(`SUMLEV`:`AGEGRP`, POPESTIMATE2004) %>%
  filter(SEX == 0) %>%  #total
  rename(statecode = STATE,
         countycode = COUNTY,
         state = STNAME,
         county = CTYNAME,
         agegroup = AGEGRP)

pop04 <- pop04 %>% 
  group_by(statecode, countycode, state, county, agegroup) %>% 
  summarize(population = sum(POPESTIMATE2004))

# find voter aged population
pop04 <- pop04 %>% 
  mutate(population = ifelse(agegroup == 4, round(population*0.4), population))

pop04 <- pop04 %>% 
  group_by(statecode, countycode, state, county) %>% 
  filter(agegroup >= 4) %>% 
  summarize(population = sum(population))

pop04 <- pop04 %>% 
  unite(FIPS, statecode:countycode, remove = TRUE, sep = "") %>% 
  select(-county) %>% 
  rename(pop04 = population)

pop04$FIPS <- as.numeric(pop04$FIPS)

# merge to turnout04
turnout04 <- left_join(countyvote04, pop04, by = "FIPS")

turnout04 <- turnout04 %>% 
  select(-state)

turnout04$FIPS <- as.numeric(as.character(turnout04$FIPS))

turnout04 <- left_join(turnout04, one_dma_states, by = "state_po", all = TRUE)

turnout04 <- turnout04 %>% 
  filter(count != 1) %>% 
  select(-count)

names(turnout04)

turnout04 <- turnout04 %>% 
  arrange(FIPS) %>% 
  select(FIPS, votes04, pop04)

############
#TURNOUT '00
############
countyvote00 <- countyvote_raw %>% 
  filter(year == "2000" & party != "NA") %>% 
  select(state_po:FIPS, totalvotes) %>% 
  distinct()

countyvote00 <- countyvote00 %>% 
  filter(totalvotes != 0)

# Alaska data is not by county so we delete it
countyvote00 <- countyvote00 %>% 
  filter(FIPS <= 2000 | FIPS > 3000) %>% 
  rename(votes00 = totalvotes)

pop00 <- population %>% 
  select(`SUMLEV`:`AGEGRP`, POPESTIMATE2000) %>%
  filter(SEX == 0) %>%  #total
  rename(statecode = STATE,
         countycode = COUNTY,
         state = STNAME,
         county = CTYNAME,
         agegroup = AGEGRP)

pop00 <- pop00 %>% 
  group_by(statecode, countycode, state, county, agegroup) %>% 
  summarize(population = sum(POPESTIMATE2000))

# find voter aged population
pop00 <- pop00 %>% 
  mutate(population = ifelse(agegroup == 4, round(population*0.4), population))

pop00 <- pop00 %>% 
  group_by(statecode, countycode, state, county) %>% 
  filter(agegroup >= 4) %>% 
  summarize(population = sum(population))

pop00 <- pop00 %>% 
  unite(FIPS, statecode:countycode, remove = TRUE, sep = "") %>% 
  select(-county) %>% 
  rename(pop00 = population)

pop00$FIPS <- as.numeric(pop00$FIPS)

#merge to turnout00
turnout00 <- left_join(countyvote00, pop00, by = "FIPS")

turnout00 <- turnout00 %>% 
  select(-state)

turnout00$FIPS <- as.numeric(as.character(turnout00$FIPS))

turnout00 <- left_join(turnout00, one_dma_states, by = "state_po", all = TRUE)

turnout00 <- turnout00 %>% 
  filter(count != 1) %>% 
  select(-count)

names(turnout00)

turnout00 <- turnout00 %>% 
  arrange(FIPS) %>% 
  select(FIPS, votes00, pop00)


#####################################
#       CANDIDATE APPEARANCE        #
#####################################
appearance <- appearance_raw %>%
  select(state, fips, republican, democratic, visits) %>% 
  rename(FIPS = fips)

# merge appearance data with turnout data
turnout_appearance <- left_join(turnout08, appearance, by = "FIPS")
turnout_appearance <- left_join(turnout_appearance, turnout04, by = "FIPS")
turnout_appearance <- left_join(turnout_appearance, turnout00, by = "FIPS")

turnout_appearance <- turnout_appearance %>% 
  mutate(republican = ifelse(is.na(republican), 0, republican),
         democratic = ifelse(is.na(democratic), 0, democratic),
         visits = ifelse(is.na(visits), 0, visits))

############################################################################################
######################                MERGING KEY                   ########################
############################################################################################
#####################################
#   CREATE A MERGING KEY FOR DMAs   #
#####################################
#we only focus on states and markets to begin with and clean these categories
wiscdma <- wisc_totalads %>% 
  select(state_po, market)

#first we try to reduce difference in strings by changing all to capital letters
wiscdma$market <- str_to_upper(wiscdma$market)
#then remove comma and whitespace
wiscdma$market <- str_trim(wiscdma$market)
wiscdma$state_po <- str_trim(wiscdma$state_po)
#then remove spaces and everything after the first word
wiscdma$market <- gsub("[,]", "", wiscdma$market)
wiscdma$market <- gsub("[.]", "", wiscdma$market)
wiscdma$market <- gsub("[-]", " ", wiscdma$market)
wiscdma$market <- gsub("[/]", " ", wiscdma$market)
wiscdma$market <- gsub("? .*", "", wiscdma$market) #This deletes everything after the first space

#then we gather by using leftjoin 
dmacounty_pre <- inner_join(wiscdma, dmaindex, by = c("market"))

dmacounty_pre <- dmacounty_pre %>% 
  rename(state_po = STATE) %>% 
  select(-state_po.y, -state_po.x)


#Some names are the same after removing everything after the first word, so we change these manually
dmacounty <- dmacounty_pre %>% 
  mutate(DMAINDEX = case_when(market == "FORT" ~ 108,
                              market == "IDAHO" ~ 165,
                              market == "JONESBORO" ~ 181,
                              market == "PALM" ~ 2,
                              market == "MANCHESTER" ~ 6,
                              market == "CHARLOTTESVILL" ~ 186,
                              market == "BISMARCK" ~ 155, #The name is Minot-Bismarck
                              market == "MYRTLE" ~ 110, #The name is Florence-Myrtle
                              TRUE ~ as.numeric((DMAINDEX))))

dmacounty$FIPS <- as.numeric(as.character(dmacounty$FIPS))

############################################################################################
######################               MERGE DATASETS                 ########################
############################################################################################
full_dataset_raw <- inner_join(dmacounty, turnout_appearance, by = "FIPS") #turnout '04 and '00 + candidate visits
full_dataset_raw <- left_join(full_dataset_raw, senateelection06, by = "FIPS") # turnout senate election '06
full_dataset_raw <- left_join(full_dataset_raw, preselection08, by = "FIPS") # turnout presidential '08
full_dataset_raw <- left_join(full_dataset_raw, senateelection02, by = "FIPS") # turnout senate election '02
full_dataset_raw <- left_join(full_dataset_raw, DMA, by = "FIPS") # dmacode for map

full_dataset_raw <- full_dataset_raw %>% 
  filter(!is.na(pop08),
         !is.na(pop00))

#Use senate data from 2002 if there isn't any from 2006
full_dataset_raw <- full_dataset_raw %>% 
  mutate(pop06 = ifelse(votes06 == 0, pop02, pop06),
         votes06 = ifelse(votes06 == 0, votes02, votes06)) %>% 
  select(-votes02, -pop02)

names(full_dataset_raw)

full_dataset_pre <- full_dataset_raw %>% 
  group_by(market, DMAINDEX, state_po) %>% 
  summarize(votes08 = sum(votes08),
            pop08 = sum(pop08),
            turnout08 = votes08/pop08,
            votes04 = sum(votes04),
            pop04 = sum(pop04),
            turnout04 = votes04/pop04,
            votes00 = sum(votes00),
            pop00 = sum(pop00),
            turnout00 = votes00/pop00,
            votes06 = sum(votes06),
            pop06 = sum(pop06),
            senate = votes06/pop06,
            votes08_nanda = sum(ballots_cast),
            pop = sum(cvap),
            nanda08 = votes08_nanda/pop,
            nanda_turnout08 = votes08/pop,
            visits_rep = sum(republican),
            visits_dem = sum(democratic),
            totalvisits = sum(visits),
            dmacode = unique(dmacode)) %>% 
  select(-votes04, -pop04, -votes00, -pop00, -votes08, -votes08_nanda, -pop08, -votes06, -pop06) %>% 
  ungroup()

full_dataset_pre <- full_dataset_pre %>% 
  mutate(nanda08 = ifelse(is.na(nanda08), nanda_turnout08, nanda08))

#Some names are the same after removing everything after the first word, so we change these manually
full_dataset_pre <- full_dataset_pre %>% 
  mutate(market = case_when(DMAINDEX == 5 ~ "SAN FRANCISCO",
                            DMAINDEX == 26 ~ "SAN DIEGO",
                            DMAINDEX == 35 ~ "GREENVILLE SC",
                            DMAINDEX == 103 ~ "GREENVILLE NC",
                            DMAINDEX == 140 ~ "SIOUX CITY",
                            DMAINDEX == 112 ~ "SIOUX FALLS",
                            DMAINDEX == 193 ~ "ST JOSEPH",
                            DMAINDEX == 22 ~ "ST LOUIS",
                            DMAINDEX == 148 ~ "ALBANY, NY",
                            DMAINDEX == 105 ~ "CHARLESTON, SC",
                            DMAINDEX == 139 ~ "COLUMBIA, MO",
                            DMAINDEX == 131 ~ "COLUMBUS, AL",
                            DMAINDEX == 34 ~ "COLUMBUS, OH",
                            DMAINDEX == 70 ~ "FT MYERS",
                            DMAINDEX == 104 ~ "FT WAYNE",
                            DMAINDEX == 38 ~ "GRAND RAPIDS",
                            DMAINDEX == 183 ~ "JACKSON, TN",
                            DMAINDEX == 125 ~ "LAFAYETTE, LA",
                            DMAINDEX == 76 ~ "PORTLAND, ME",
                            DMAINDEX == 37 ~ "SAN ANTONIO",
                            DMAINDEX == 142 ~ "WICHITA FALLS",
                            TRUE ~ as.character((market))))

full_dataset <- left_join(full_dataset_pre, wisc_totalads, by = c("market"))
full_dataset <- left_join(full_dataset, wiscads_totalGSH, by = c("market")) 
names(full_dataset)

full_dataset <- full_dataset %>% 
  filter(!is.na(totalads)) %>% 
  select(-state_po) %>% 
  rename(state_po = state_po.x)

#####################################
#     GENERATING DATASET FOR MAP    #
#####################################
full_dataset_map <- full_dataset %>% # generating data for map
  distinct(market, state_po, totalads, .keep_all = TRUE) %>% #Remove all duplicates
  subset(!state_po %in% c("CT ", "DC ", "DE ", "ND ", "SD ")) %>% #removes states with only one DMA
  select(-state_po.y) %>%         
  arrange(state_po) 

#####################################
#GENERATING DATASET FOR REGRESSIONS #
#####################################
full_dataset <- full_dataset_map %>% #generating data for models
  filter(!is.na(senate)) #removes markets with no senate data

n_distinct(full_dataset$state_po)
n_distinct(full_dataset$market)
unique(full_dataset$state_po)

# States that have been excluded
#KY, WI, MA, MD, NY, AK, HI, UT, RI #no ads or senate data
#CT, DC, DE, ND, SD #only one DMA

#####################################
#      CONSTRUCTING VARIABLES       #
#####################################
full_dataset <- full_dataset %>% 
  mutate(lncpop = log(pop),
         pq = nanda08*(1-nanda08),
         senateads = ifelse(is.na(senateads), 0, senateads), 
         houseads = ifelse(is.na(houseads), 0, houseads),
         contrastads = contrastads/10000,
         promoads = promoads/10000,
         attackads = attackads/10000,
         totalads = totalads/10000,
         senateads = senateads/10000,
         houseads = houseads/10000,
         diffpromo = (promoads - attackads)*10000)


##############################################################################################################
###############################                 2008-DATA                    #################################
###############################       ANALYSIS AND REGRESSION MODELS         #################################
###############################                                              #################################
##############################################################################################################

############################################################################################
######################              REGRESSION MODELS               ########################
############################################################################################
###########
# Model 1 #
###########
model1 <- lm(nanda08 ~ totalads + senate,
             data = full_dataset)

# Test for heteroskedasticity
bptest(model1, nanda08 ~ lncpop+pq, 
       data = full_dataset, studentize=FALSE) #HETEROSKEDASTIC
bptest(model1)

# Robust standard errors
coeftest(model1, vcov = vcovHC(model1, type="HC1")) 
model1r <- lmrob(nanda08 ~ totalads + senate,
                 data = full_dataset)

summary(model1)

###########
# Model 2 #
###########
model2 <- felm(nanda08 ~ totalads + senate 
               | factor(state_po),
               data = full_dataset) 

# Test for heteroskedasticity
bptest(model2, nanda08 ~ lncpop+pq, 
       data = full_dataset, studentize=FALSE) #HETEROSKEDASTIC

# Robust standard errors
coeftest(model2, vcov = vcovHC(model2, type="HC1")) # robust standard errors

summary(model2)


###########
# Model 3 #
###########
model3 <- felm(nanda08 ~ totalads + senate + turnout04 + turnout00 
               | factor(state_po),
               data = full_dataset)

# Test for heteroskedasticity
bptest(model3, nanda08 ~ lncpop+pq, 
       data = full_dataset, studentize=FALSE) #HETEROSKEDASTIC

# Robust standard errors
coeftest(model3, vcov = vcovHC(model3, type="HC1")) # robust standard errors

summary(model3)


###########
# Model 4 #
###########
model4 <- lm(nanda08 ~ totalads + senate + turnout04 + turnout00,
             data = full_dataset)

# Test for heteroskedasticity
bptest(model4, nanda08 ~ lncpop+pq, 
       data = full_dataset, studentize=FALSE) #HETEROSKEDASTIC

# Robust standard errors
coeftest(model4, vcov = vcovHC(model4, type="HC1")) # robust standard errors

summary(model4)


###########
# Model 5 #
###########
model5 <- felm(nanda08 ~ totalads + senate + turnout04 + turnout00 + senateads + houseads 
               | factor(state_po),
               data = full_dataset)

# Test for heteroskedasticity
bptest(model5, nanda08 ~ lncpop+pq, 
       data = full_dataset, studentize=FALSE) #HETEROSKEDASTIC

# Robust standard errors
coeftest(model5, vcov = vcovHC(model5, type="HC1")) # robust standard errors

summary(model5)


###########
# Model 6 #
###########
model6 <- felm(nanda08 ~ totalads + senate + turnout04 + turnout00 + senateads + houseads + totalvisits 
               | factor(state_po), 
               data = full_dataset)

# Test for heteroskedasticity
bptest(model6, nanda08 ~ lncpop+pq, 
       data = full_dataset, studentize=FALSE) #HETEROSKEDASTIC

# Robust standard errors
coeftest(model6, vcov = vcovHC(model6, type="HC1")) # robust standard errors

summary(model6)


###########
# Model 7 #
###########
model7 <- lm(nanda08 ~ totalads + senate + turnout04 + turnout00 + senateads + houseads + totalvisits,
             data = full_dataset)

# Test for heteroskedasticity
bptest(model7, nanda08 ~ lncpop+pq, 
       data = full_dataset, studentize=FALSE) #HETEROSKEDASTIC

# Robust standard errors
coeftest(model7, vcov = vcovHC(model7, type="HC1")) # robust standard errors

summary(model7)

##############
# All models #
##############

summary(model1)
summary(model2)
summary(model3)
summary(model4)
summary(model5)
summary(model6)
summary(model7)

stargazer(model1, model2, model3, model4, model5, model6, model7,
          type = "latex",
          model.names = FALSE, # removes OLS and felm
          dep.var.labels = c("Turnout '08"),
          covariate.labels = "Total ads (divided by 1000)",
          omit = c("senate", "turnout04", "turnout00", "senateads", "houseads", "totalvisits", "Constant"),
          add.lines = list(c("Midterm turnout",                   "Yes", "Yes", "Yes", "Yes", "Yes", "Yes", "Yes"),
                           c("State fixed effects",               "", "Yes", "Yes", "No", "Yes", "Yes", "No"),
                           c("Past presidential turnout",         "",  "", "Yes",  "Yes", "Yes", "Yes", "Yes"),
                           c("Control for other election ads",    "",  "",  "", "", "Yes", "Yes", "Yes"),
                           c("Control for candidate visits",      "",  "",  "", "", "", "Yes", "Yes")),
          star.cutoffs = c(0.05, 0.01, 0.001),
          digits = 2,
          keep.stat = c("adj.rsq", "n"))


############################################################################################
######################                ROBUSTNESS TESTS              ########################
############################################################################################
#####################################
#              AD TONE              #
#####################################
###########
# Model 1 #
###########
tone_model1 <- lm(nanda08 ~ attackads + promoads + contrastads + senate,
                  data = full_dataset)
summary(tone_model1)

base_model1 <- lm(nanda08 ~ senate,
                  data = full_dataset)
summary(base_model1)

anova(tone_model1, model1)$"Pr(>F)"[2] #INS, test whether the effect of media tone is significantly different than total ads
anova(tone_model1, base_model1)$"Pr(>F)"[2] #SIG, test whether the effect of media tone is signicantly different than without ads


###########
# Model 2 #
###########
tone_model2 <- lm(nanda08 ~ attackads + promoads + contrastads + senate + state_po - 1,
                  data = full_dataset)
summary(tone_model2)

base_model2 <- lm(nanda08 ~ senate + state_po - 1,
                  data = full_dataset)
summary(base_model2)

anova(tone_model2, model2)$"Pr(>F)"[2] #INS, test whether the effect of media tone is significantly different than total ads
anova(tone_model2, base_model2)$"Pr(>F)"[2] #INS, test whether the effect of media tone is signicantly different than without ads

###########
# Model 3 #
###########
tone_model3 <- lm(nanda08 ~ attackads + promoads + contrastads + senate + turnout04 + turnout00 + state_po - 1,
                  data = full_dataset)
summary(tone_model3)

base_model3 <- lm(nanda08 ~ senate + turnout04 + turnout00 + state_po - 1,
                  data = full_dataset)
summary(base_model3)

anova(tone_model3, model3)$"Pr(>F)"[2] #INS, test whether the effect of media tone is significantly different than total ads
anova(tone_model3, base_model3)$"Pr(>F)"[2] #INS, test whether the effect of media tone is signicantly different than without ads


###########
# Model 4 #
###########
tone_model4 <- lm(nanda08 ~ attackads + promoads + contrastads + senate + turnout04 + turnout00,
                  data = full_dataset)
summary(tone_model4)

base_model4 <- lm(nanda08 ~ senate + turnout04 + turnout00,
                  data = full_dataset)
summary(base_model4)

anova(tone_model4, model4)$"Pr(>F)"[2] #INS, test whether the effect of media tone is significantly different than total ads
anova(tone_model4, base_model4)$"Pr(>F)"[2] #SIG, test whether the effect of media tone is signicantly different than without ads


###########
# Model 5 #
###########
tone_model5 <- lm(nanda08 ~ attackads + promoads + contrastads + senate + turnout04 + turnout00 + senateads + houseads + state_po - 1,
                  data = full_dataset)
summary(tone_model5)

base_model5 <- lm(nanda08 ~ senate + turnout04 + turnout00 + senateads + houseads + state_po - 1,
                  data = full_dataset)
summary(base_model5)

anova(tone_model5, model5)$"Pr(>F)"[2] #INS, test whether the effect of media tone is significantly different than total ads
anova(tone_model5, base_model5)$"Pr(>F)"[2] #INS, test whether the effect of media tone is signicantly different than without ads


###########
# Model 6 #
###########
tone_model6 <- lm(nanda08 ~ attackads + promoads + contrastads + senate + turnout04 + turnout00 + senateads + houseads + totalvisits + state_po - 1,
                  data = full_dataset)
summary(tone_model6)

base_model6 <- lm(nanda08 ~ senate + turnout04 + turnout00 + senateads + houseads + totalvisits + state_po - 1,
                  data = full_dataset)
summary(base_model6)

anova(tone_model6, model6)$"Pr(>F)"[2] #INS, test whether the effect of media tone is significantly different than total ads
anova(tone_model6, base_model6)$"Pr(>F)"[2] #INS, test whether the effect of media tone is signicantly different than without ads


###########
# Model 7 #
###########
tone_model7 <- lm(nanda08 ~ attackads + promoads + contrastads + senate + turnout04 + turnout00 + senateads + houseads + totalvisits,
                  data = full_dataset)
summary(tone_model7)

base_model7 <- lm(nanda08 ~ senate + turnout04 + turnout00 + senateads + houseads + totalvisits,
                  data = full_dataset)
summary(base_model7)

anova(tone_model7, model7)$"Pr(>F)"[2] #INS, test whether the effect of media tone is significantly different than total ads
anova(tone_model7, base_model7)$"Pr(>F)"[2] #INS, test whether the effect of media tone is signicantly different than without ads

#####################################
#         DIFF PROMO-ATTACK         #
#####################################
###########
# Model 1 #
###########
model1_diff <- lm(nanda08 ~ diffpromo + senate,
                  data = full_dataset)

summary(model1_diff)

###########
# Model 2 #
###########
model2_diff <- felm(nanda08 ~ diffpromo + senate 
                    | factor(state_po),
                    data = full_dataset) 

summary(model2_diff)


###########
# Model 3 #
###########
model3_diff <- felm(nanda08 ~ diffpromo + senate + turnout04 + turnout00 
                    | factor(state_po),
                    data = full_dataset)

summary(model3_diff)


###########
# Model 4 #
###########
model4_diff <- lm(nanda08 ~ diffpromo + senate + turnout04 + turnout00,
                  data = full_dataset)

summary(model4_diff)


###########
# Model 5 #
###########
model5_diff <- felm(nanda08 ~ diffpromo + senate + turnout04 + turnout00 + senateads + houseads 
                    | factor(state_po),
                    data = full_dataset)

summary(model5_diff)


###########
# Model 6 #
###########
model6_diff <- felm(nanda08 ~ diffpromo + senate + turnout04 + turnout00 + senateads + houseads + totalvisits 
                    | factor(state_po), 
                    data = full_dataset)

summary(model6_diff)


###########
# Model 7 #
###########
model7_diff <- lm(nanda08 ~ diffpromo + senate + turnout04 + turnout00 + senateads + houseads + totalvisits,
                  data = full_dataset)

summary(model7_diff)




##############################################################################################################
###############################                  2008 DATA                   #################################
###############################     GRDD WITH DATA FROM KEELE & TITIUNIK     #################################
###############################                                              #################################
##############################################################################################################

############################################################################################
######################                LOADING DATA                  ########################
############################################################################################
KT <- read_dta("Voters_Final.dta")
BDR <- read.dbf("BorderSegmentPoints_Project.dbf")

KT <- KT %>% 
  filter(schoolid == 1) # Only include individuals who are in the same school district

############################################################################################
######################                DATA TIDYING                  ########################
############################################################################################
BDR <- BDR %>% 
  rename(latitude = POINT_Y, 
         longitude = POINT_X,
         BDR_uniqueid = StepID)

KT_geo <- KT %>% 
  select(uniqueid, latitude, longitude, treat_county) %>% 
  rename(KT_uniqueid = uniqueid)

BDR_geo <- BDR %>% 
  select(BDR_uniqueid, latitude, longitude)


combinations <- CJ(KT$uniqueid, BDR$BDR_uniqueid)
names(combinations) <- c("KT_uniqueid", "BDR_uniqueid")

geo <- combinations %>% 
  left_join(KT_geo, by = "KT_uniqueid")

geo <- geo %>% 
  left_join(BDR_geo, by = "BDR_uniqueid") %>% 
  rename(KT_latitude = latitude.x,
         KT_longitude = longitude.x,
         BDR_latitude = latitude.y,
         BDR_longitude = longitude.y)

#####################################
#       CALCULATE DISTANCES         #
#####################################
geo[, distance := distHaversine(matrix(c(KT_longitude, KT_latitude), ncol = 2),
                                matrix(c(BDR_longitude, BDR_latitude), ncol = 2))]

geodata <- geo %>% 
  group_by(KT_uniqueid) %>% 
  summarize(distance = min(distance),
            treat_county = unique(treat_county))

geodata <- geodata %>% 
  mutate(distance = ifelse(treat_county == 0, distance * -1, distance)) %>% 
  rename(uniqueid = KT_uniqueid)


KT_merged <- KT %>% 
  left_join(geodata, by = "uniqueid")


############################################################################################
######################           Pre-regression analysis            ########################
############################################################################################
DCdensity(KT_merged$distance, cutpoint = 0)


############################################################################################
######################                 REGRESSIONS                  ########################
############################################################################################
#####################################
#         Model 1 - linear          #
#####################################
#Regression with optimal banwidth
rd_model1 <- rdrobust(KT_merged$e2008g, KT_merged$distance, c = 0, p = 1, kernel = "triangular")
summary(rd_model1)
# OBS = 2044

rd_model1_c <- rdrobust(KT_merged$e2008g, KT_merged$distance, covs = KT_merged$income, c = 0, p = 1)
summary(rd_model1_c)

?rdplot

# Plot
pdf("/Users/jakobhanghoej/Desktop/Statskundskab/9. semester/Advanced Quantative Methods/01_Exam_AQM/New_data/Plots/Rdplot_1.pdf") #start print device
rdplot(KT_merged$e2008g, KT_merged$distance, 
       c = 0, 
       p = 1,
       h = 797.398, #optimal bandwidth from regression above
       x.lim = c(-800, 800), 
       y.lim = c(0, 1),
       x.label = "Distance to the border between DMAs (meters)",
       y.label = "Proportion of voters",
       title = "",
       col.dots = "grey",
       col.lines = "darkred",
       kernel = "triangular")
dev.off() #stopping print device

#####################################
#       Model 2 - 2 dg. poly       #
#####################################
rd_model2 <- rdrobust(KT_merged$e2008g, KT_merged$distance, c = 0, p = 2, kernel = "triangular")
summary(rd_model2)
# OBS = 5895

rd_model2_c <- rdrobust(KT_merged$e2008g, KT_merged$distance, covs = KT_merged$income, c = 0, p = 2)
summary(rd_model2_c)

# Plot
pdf("/Users/jakobhanghoej/Desktop/Statskundskab/9. semester/Advanced Quantative Methods/01_Exam_AQM/New_data/Plots/Rdplot_2.pdf") #start print device
rdplot(KT_merged$e2008g, KT_merged$distance, 
       c = 0, 
       p = 2, 
       h = 1672.799, #optimal bandwidth from regression above
       x.lim = c(-1800, 1800), 
       y.lim = c(0, 1),
       x.label = "Distance to the border between DMAs (meters)",
       y.label = "",
       title = "",
       col.dots = "grey",
       col.lines = "darkred",
       kernel = "triangular")
dev.off() #stopping print device

#####################################
#       Model 3 - 3 dg. poly       #
#####################################
rd_model3 <- rdrobust(KT_merged$e2008g, KT_merged$distance, c = 0, p = 3, kernel = "triangular")
summary(rd_model3)
# OBS = 10,085

rd_model3_c <- rdrobust(KT_merged$e2008g, KT_merged$distance, covs = KT_merged$income, c = 0, p = 3)
summary(rd_model3_c)

# Plot
pdf("/Users/jakobhanghoej/Desktop/Statskundskab/9. semester/Advanced Quantative Methods/01_Exam_AQM/New_data/Plots/Rdplot_3.pdf") #start print device
rdplot(KT_merged$e2008g, KT_merged$distance, 
       c = 0, 
       p = 3, 
       h = 2321.962, #optimal bandwidth from regression above
       x.lim = c(-2500, 2500), 
       y.lim = c(0, 1),
       x.label = "Distance to the border between DMAs (meters)",
       y.label = "Proportion of voters",
       title = "",
       col.dots = "grey",
       col.lines = "darkred", 
       kernel = "triangular")
dev.off()

#####################################
#       Model 4 - 4 dg. poly       #
#####################################
rd_model4 <- rdrobust(KT_merged$e2008g, KT_merged$distance, c = 0, p = 4, kernel = "triangular")
summary(rd_model4)
# OBS = 6,710
?rdplot
# Plot
pdf("/Users/jakobhanghoej/Desktop/Statskundskab/9. semester/Advanced Quantative Methods/01_Exam_AQM/New_data/Plots/Rdplot_4.pdf") #start print device
rdplot(KT_merged$e2008g, KT_merged$distance, 
       c = 0, 
       p = 4, 
       h = 1848.973, #optimal bandwidth from regression above
       x.lim = c(-2000, 2000), 
       y.lim = c(0, 1),
       x.label = "Distance to the border between DMAs (meters)",
       y.label = "",
       title = "",
       col.dots = "grey",
       col.lines = "darkred",
       kernel = "triangular")
dev.off()

summary(rd_model1)
summary(rd_model2)
summary(rd_model3)
summary(rd_model4)


##############################################################################################################
###############################                  2008 DATA                   #################################
###############################                   MAP DATA                   #################################
###############################                                              #################################
##############################################################################################################

##############################################
#          LOADING AND READING DATA          #
##############################################
mapdf <- full_dataset_map %>% 
  select(market, dmacode, state_po, totalads)

gpclibPermit()

dma <- readOGR("dma_2008/DMAs.shp")

dma.df <- fortify(dma, region = "DMA")
dma.df <- rename(dma.df, DMA = id)

#Merge with shape data
class(dma.df$DMA)
dma.df <- transform(dma.df, DMA = as.numeric(DMA)) #turn to numeric
class(mapdf$dmacode)
mapdf <- transform(mapdf, dmacode = as.numeric(dmacode)) #turn to numeric

dma_mapdf <- dma.df %>% 
  left_join(mapdf, by = c("DMA" = "dmacode"))

dma_mapdf <- transform(dma_mapdf, totalads = as.numeric(totalads)) #turn to numeric

###############################
#          PLOT DMAS          #
###############################
states_map <- map_data("state")

legend_title <- "Total ads shown"

ggplot(dma_mapdf, aes(x=long, y=lat, group=group, fill=totalads)) + 
  geom_polygon(color="#666666", size=.5) +
  scale_fill_continuous(legend_title, low = "lightskyblue", high="deepskyblue4", na.value="white", limits=c(0,7200), breaks=c(1000,3000,5000,7000)) +
  geom_polygon(data = states_map, aes(x = long, y = lat, group = group), fill=NA, color="black") +
  coord_map() +
  theme_void() +
  ggsave("map_states_2008c.png")

#########################################
#          FOCUS ON NEW JERSEY          #
#########################################
nj_stat <- states_map %>% 
  filter(region %in% c("new jersey"))

pa_nj_stat <- states_map %>% 
  filter(region %in% c("pennsylvania", "new jersey"))

pa_nj_dma <- dma_mapdf %>% 
  filter(state_po %in% c("PA ", "NJ "))

nj_dma <- pa_nj_dma %>% 
  filter(market == "PHILADELPHIA")
nj_dma$state_po <- "NJ"

#Sites
sites <- data.frame(longitude = c(-74.6), latitude = c(40.31))

ggplot(nj_dma, aes(x=long, y=lat, group=group, fill=totalads)) +
  geom_polygon(color="#666666", size=.5) +
  geom_polygon(data = nj_stat, aes(x = long, y = lat, group = group), fill=NA, color="black") + # states' borders
  scale_fill_continuous(legend_title, low = "white", high=muted("deepskyblue2"), na.value="white", limits=c(0,5035), breaks=c(0,5000)) +
  coord_map() +
  theme_void() +
  annotate("text", x = -74.57, y = 41.05, label = "Area of interest", size=5) +
  annotate("text", x=-76.1, y=41.2, label = "Philadelphia DMA", size=5) +
  annotate("text", x=-76.1, y=41.0, label = "New York DMA", size=5) +
  annotate("label", x = -75.45, y = 41.2, label = "     ", fill=muted("deepskyblue2")) +
  annotate("label", x = -75.45, y = 41.0, label = "     ") +
  annotate("label", x = -74.57, y = 40.32, label = "       \n       ",alpha=1,fill=NA) +
  annotate("segment", x = -74.57, xend = -74.57, y = 41, yend = 40.45, colour = "black", size=0.5, alpha=1, arrow=arrow()) +
  ggsave("New Jersey map.png")



