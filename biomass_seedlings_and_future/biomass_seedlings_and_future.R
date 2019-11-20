#Seedling Biomass, Adult Biomass, Demographic Model Future Biomass####
#1. Seedlings####
#WD, packages####
setwd("/Users/Ana/Dropbox/2017\ Palmyra\ Seed\ Predation\ Paper/Analyses/Draft_1_Mar_2019/biomass_seedlings_and_future")

library(tidyverse)
library(car)
library(nlstools)
library(ggplot2)
library(knitr)
library(directlabels)
library(lme4)
library(DHARMa)
library(RVAideMemoire)
library(MuMIn)
library(cowplot)
library(pbkrtest)
library(popbio)
#CN Seedling Biomass Eqn####
## Coconut Seedling Biomass Determination

#The seedlings in the permanent plots range in size and age, and thus, biomass. 
#To determine how much total seedling biomass is in each plot, we determined an 
#age-specific biomass relationship based on relationships determined between 
#leaf length, leaf area, leaf biomass, and age. We used a combination of data 
#from the literature, leaf images (leaf length to leaf area), and seedlings grown
#in a nursery experiment (leaf area to biomass) to determine the relationship 
#between biomass and age for a group of caged and uncaged seedlings tracked over
#a period of eight years (total leaf number and leaf lengths measured in each 
#year). We used this age-specific biomass relationship in combination with the 
#age distributions of seedlings in our plots to determine the total seedling 
#biomass in each plot. 

# Step 1: Validate relationship between lamina length and leaf area #
#(literature and leaf image data)

#From the literature, we determined that the fronds of Cocos ages 
#1 - 4 follow a relationship of the form:
  
#  1. Leaf Area = 0.8282*(rachis Length)^1.5662 (Figure 1 in R Markdown)

#These data were determined from 57 leaves ranging from 48 to 330 cm, 
#from plants ranging from age one to four years in a plantation in Campos 
#dos Goytacazes, RJ, Brazil (de Sousa et al. 2005). The rachis length was measured 
#in the field and the area of each leaf was determined using a LI-3100 area meter 
#(LICOR Inc., Lincoln, Nebraska, USA). This relationship was validated against 31
#additional leaves from four other plantations. For these additional leaves, the 
#expected leaf area from the model was compared to the observed leaf area to determine 
#the predictability of this model (Figure 2). 

#![Figure 1](Screen\ Shot\ 2018-06-14\ at\ 9.43.10\ AM.png) 
#![Figure 2](Screen\ Shot\ 2018-06-14\ at\ 9.43.20\ AM.png)

#We independently fit a power model to leaves collected and imaged from Palmyra 
#(Young et al. 2011) to determine if this relationship was similar to one 
#determined on leaves specific to Palmyra. We used images of 42 leaves that ranged 
#in length from 7.10 to 84.58 cm. We fit a power model to the relationship between 
#leaf length and leaf area. 

length <- read.csv("Cocos scan length to LA.csv")

length$area.cm.2 <- as.numeric(as.character(length$area.cm.2))
length <- length[which(length$diagonal.length.cm != "NA"),]
length <- length[which(length$area.cm.2 != "NA"), ]

#power relationship between area and length
m <- nls(area.cm.2 ~ diagonal.length.cm^b, data = length, start = list(b=1))
summary(m)

#goodness of fit via an R^2
RSS.p <- sum(residuals(m)^2)
TSS <- sum((length$area.cm.2 - mean(length$area.cm.2))^2)
1 - (RSS.p/TSS) #this is R^2 value for this model fit

#visualize
ggplot(length, aes(x = diagonal.length.cm, y = area.cm.2)) +
  geom_point() + 
  stat_smooth(method = 'nls', formula = 'y~x^b', start = list(b = 1), se=FALSE, color = "black") +
  geom_smooth(se = FALSE, method="lm", formula = y ~ I(x^1.5662), color = "black", lty = 2) +
  xlim(0, 100) +
  annotate("text", x = 95, y = 875, label = "Literature") +
  annotate("text", x = 95, y = 815, label = "Palmyra") +
  labs(title = "Power model fit of Palmyra seedling leaves", x = "Leaf length (cm)", y = "Leaf area (cm^2)") +
  theme_bw()

summary(m)$coefficients

#Our results indicate that there is a significant power relationship between leaf 
#diagonal length and leaf area (p-value < 0.001). The R^2 of this model fit is 0.67, 
#and the relationship is given as:
  
#2. Lamina Area = Lamina Length ^ 1.51 **is this lamina and not rachis length?**
  
#Furthermore, this relationship is very close to the relationship in the literature.

# Step 2: Determine a relationship between biomass and leaf area (nursery data - 
#Young et al. 2011)

#From the nursery experiments, we determined a relationship between biomass and 
#leaf area for a total of 52 seedlings. Dry biomass was taken on these seedlings, 
#including leaf and stem biomass, as well as specific leaf area (cm^2/dry g). We 
#determined total leaf area by multiplying specific leaf area by leaf dry biomass. 
#We fit a linear model to the relationship between total leaf area and total dry 
#biomass (leaves and stem).

nursery <- read.csv("Nursury_CN_LA_Biomass.csv")


row.has.na <- apply(nursery, 1, function(x){any(is.na(x))})
sum(row.has.na)
nursery <- nursery[!row.has.na,]

nursery$leaf.bio <- nursery$Lamina.weight + nursery$Stem.weights


model <- lm(leaf.bio ~ Lamina.Area, data=nursery)
summary(model)

par(mfrow=c(1,1))
plot(model, which=c(1,2))

par(mfrow=c(1,1))
plot(leaf.bio ~ Lamina.Area, data=nursery, main = "Total leaf biomass (stem and lamina) by leaf area", 
     xlab = "Leaf area (cm^2)", ylab = "Leaf Biomass (dry grams)")
abline(lm(leaf.bio ~ Lamina.Area, data=nursery))

coef(model)
#gives coefficients of a model of the form leaf.bio = 49.89 + 0.22*Lamina.Area

#Our data suggest a significant relationship between leaf area and leaf+stem biomass 
#(p-value < 0.001). The adjusted R^2 is 0.78 and the relationship is given as:
  
 # 3. Leaf Biomass (lamina + stem) = 49.89 + 0.22*Lamina Area

#Step 3: Determine age-specific biomass (cage seedling data - Young et al. 2013)

#The relationship determined from the literature for leaf length to area matched 
#the one determined by imaged leaves fairly closely, so we used the relationship from 
#the literature, since this was drawn from leaves with a broader size range. 

#We determined lamina area per frond on each seedling using equation 1, and added 
#the totals to get a total leaf area for each seedling. We used this total leaf area 
#to determine total leaf dry biomass using equation 3. 

#We then fit a power function predicting leaf dry weight by measurement time 
#(a yearly number between 1-8), #with 1 being the first measurement in 2007, 
#and 8 being the last measurement in 2014. 
#We first fit a model including a variable indicating whether the seedling 
#was initially caged or not, and then fit a model without this caged variable. 
#We chose best model fit based on model AICs. The complete dataset includes 570 observations
#of 204 individuals taken 1-5 times per individual depending on individual seedling survival. 

ages <- read.csv("cage_LA_Bio_age.csv")

ages$cage.id[ages$Cage == "CA"] <- 2
ages$cage.id[ages$Cage == "CO"] <-1

m.age.cage <- nls(leaf.dw.g ~ a*Time^b + cage.id, data=ages, start = list (a=2, b=2),
                  trace = T)

m.age <- nls(leaf.dw.g ~ a*Time^b, data= ages, start = list (a=2, b=2),
             trace = T)

AIC(m.age.cage, m.age)

summary(m.age)
RSS.p <- sum(residuals(m.age)^2)
TSS <- sum((ages$leaf.dw.g - mean(ages$leaf.dw.g))^2)
1 - (RSS.p/TSS) #R^2 of the model

ggplot(data = ages, 
       aes(x = Time, y = leaf.dw.g)) +
  geom_point() +
  geom_smooth(se = FALSE, method = "lm", formula = y ~I(x^2.39), color="black") +
  labs(title = "Seedling total dry biomass by age", x = "Sampling year (starting at 1)", y = "Total seedling dry biomass (g)") +
  theme_bw()

#The model fit suggests a significant power function (p-value < 0.001),
#with an R^2 of 0.49, and the relationship is given as:
  
#  4. Seedling Total Dry Biomass = 54*Time ^ 2.39

#PG Seedling biomass Eqn####
#for pisonia, we used ideas from Chave et al. 2005 to determine AGB for seedlings, with
#a relationship given as: 
#AGB = 0.0776 * (p*D^2*H)^0.940
#consider something more seedling-specific here, maybe drawn from data from nursery
#experiment if i can find it.

#because we did not measure heights for seedlings in our plots, we had to first tie
#the diameters we did measure to a possible diameter-height relationship for Pisonia
height <- read.csv("2018_DBH_ht.csv")

PG_height <- height[which(height$Species == "Pisonia "),]

plot(tree_height_m ~ Diam_cm, data = PG_height) #when we assess this, it appears that
#there is some sort of logistic growth for this tree. so we will use the coefficients
#from this model to verify our starting points in the nls function

coef(lm(logit(tree_height_m) ~ Diam_cm,
        data = PG_height))

PG_height
#this is the relationship with three variables phi1, phi2, phi3 that fits a logistic model to 
#the relationship between DBH and height. phi2 and phi3 correspond to the lm coefficents
#from above, and phi1 I played around with specifying until the model fits.
ht.dbh<-nls(tree_height_m~phi1/(1+exp(-(phi2+phi3*Diam_cm))),
            start=list(phi1=100,phi2=-3.27,phi3=.015),data=PG_height,trace=TRUE, 
            control = list(maxiter = 500))

#this is a linear model to compare to the logistic model
ht.dbh.lin <- lm(tree_height_m ~ Diam_cm, data = PG_height)
AIC(ht.dbh.lin, ht.dbh) #comparing shows logistic is better predictor

ht.dbh.resid <- nlsResiduals(ht.dbh) #assess residuals
plot(ht.dbh.resid) #look pretty good 
summary(ht.dbh)

coefficients(ht.dbh) #this will give us what the equation is

#DBH to height relationship:
#TreeHt(m) = 25.97/(1+exp(-(-2.25+0.06*Tree.DBH.cm)))

#Now we have to get some idea of what biomass will be for pisonia
#average wood density (p) for Pisonia genus given as: 0.3014 g/cm^3 (Climate Action Reserve)
 #from Chave et al. 2005, AGB in biomass 
#is given as: AGB (kg) = 0.0776 x (p*D^2*H)^0.940
#where D is in cm and H is in m

#so:
#AGB(kg) = 0.0776 * (0.3014*(D^2)*H)^0.94

#PF Seedling Biomass####
#from ashish et al 2015, all values in grams
#PF Biomass for seedlings (diam 0-4cm) is 91.05 +- 1.86g, 137.9 +- 0.76g, 95.76 +- 4.64g
#PF Biomass for saplings (diam 4-8 cm) is 1333.76 +- 16.73, 841.8 +- 5.27, 864.42+-26.67
#PF Biomass for larger saplings (diam 8-12cm) is 1610.23+-2.04, 1182.43+-14.17, 1033.47 +- 16.87

#Seedling count to Biomass####
counts_CN <- read.csv("seedling_size_dists.csv")

counts_CN <- counts_CN %>%
  gather(class, count, CN_Yr1:CN_Yr5)

ggplot(counts_CN, aes(x = class, y = count, fill= Rats)) + 
  geom_bar(stat="identity", position="dodge")

#by year biomass
counts_CN$biomass_g <- ifelse(counts_CN$class == "CN_Yr1", counts_CN$count*54,
                            ifelse(counts_CN$class == "CN_Yr3", counts_CN$count*54*(3^2.39),counts_CN$count*54*(5^2.39)
                                   ))
#species
counts_CN$Species <- counts_CN$class

#PF and PG counts
PG_PF <- read.csv("PG_PF_seedling_Diams_by_Year.csv")

ggplot(PG_PF, aes(x = Year, y = Diam_cm, color = Species))+
  geom_jitter()
#subset for PG
PG_counts <- PG_PF[which(PG_PF$Species == "PG"), ]

#make a height column based on eqn determined between height and DBH
PG_counts$height_m <- 25.97/(1+exp(-(-2.25+0.06*PG_counts$Diam_cm)))

#convert this to biomass given equation from Chave et al. 2005
PG_counts$biomass_g <- 1000*(0.076 * (0.3014*(PG_counts$Diam_cm^2)*PG_counts$height_m)^0.94)

PF_counts <- PG_PF[which(PG_PF$Species == "PF"), ]

#smallest size from Ashish et al 2015
sm <- mean(c(91.05, 137.9, 95.76))
sm
#med size
md <- mean(c(1333.76, 841.8, 864.42))
md
#large size
lg <- mean(c(1610.23, 1182.43, 1033.47))
lg
#convert counts to biomass
PF_counts$biomass_g <- ifelse(PF_counts$Diam_cm > 0 & PF_counts$Diam_cm < 4, sm,
                              ifelse(PF_counts$Diam_cm >= 4 & PF_counts$Diam_cm < 8, md,lg
                                     ))
#select variables to subset
myvars <- c("Species", "Plot", "Forest_Type", "Rats", "Year", "biomass_g")

#create new dataframes that only include important variables
CN <- counts_CN[myvars]
PG <- PG_counts[myvars]
PF <- PF_counts[myvars]

#create one data frame with the right columns for biomass data
all_biomass <- rbind(CN, PG, PF)
all_biomass$Year <- as.factor(all_biomass$Year)
all_biomass$Species <- as.factor(all_biomass$Species)
all_biomass$Rats <- factor(all_biomass$Rats, levels = c("pre", "post"))

#make a new column by species that is consistent across age groups
all_biomass$Species1 <- ifelse(all_biomass$Species == "CN_Yr1", "CN",
                              ifelse(all_biomass$Species == "CN_Yr3", "CN",
                                     ifelse(all_biomass$Species == "CN_Yr5", "CN",
                                            ifelse(all_biomass$Species == "PF", "PF", "PG"))))


#Seedling Biomass LMMs####
#make a sum of biomass by plot by year to do statistics on
B_total <- all_biomass %>%
  group_by(Plot, Rats, Year, Forest_Type) %>%
  summarize(sum_biomass = sum(biomass_g, na.rm=T))

B_total$Eradication <- ifelse(B_total$Rats == "pre", 1,2)

#Figure in paper#
(c <- ggplot(B_total, aes(x = Rats, y = (sum_biomass/300) *(.000001/1) * (1/0.0001), fill = Rats)) + 
  geom_boxplot(fill = "#00A08A") +   
  labs(x = "Rats", y = "Seedling biomass (Mg/ha)") +
  theme_bw() + ylim(0,20) +
  theme(axis.text = element_text(size = 35), axis.title=element_text(size=40), legend.position="none")+
  annotate("text", x = 1.5, y = 17, label = "*", size = 22))

B_cn <- all_biomass[which(all_biomass$Species1 == "CN"),]
(c_cnonly <- ggplot(B_cn, aes(x = Rats, y = (biomass_g/300) *(.000001/1) * (1/0.0001), fill = Rats)) + 
    geom_boxplot(fill = "#00A08A") +   
    labs(x = "Rats", y = "Juvenile biomass (Mg/ha)") +
    theme_bw() + ylim(0,15) +
    theme(axis.text = element_text(size = 35), axis.title=element_text(size=40), legend.position="none")+
    annotate("text", x = 1.5, y = 12, label = "*", size = 22))
#summary values:
B_total %>%
  group_by(Rats) %>%
  summarize(mean = mean(sum_biomass)/30000)
(9.48-3.23)/(0.689 - .134)
B_cn %>%
  group_by(Rats) %>%
  summarize(mean = mean(biomass_g)/30000)
(3.23/.134)
((.134 - 3.23)/.134)*100
B_total %>%
  group_by(Rats) %>%
  summarize(sd = sd(sum_biomass)/30000)
#base model of biomass
B_M1 <- lmer(sum_biomass ~ Eradication + (1|Year/Plot),
             data=B_total)

#assess biomass model fit
simulationOutput <- simulateResiduals(fittedModel = B_M1)
plot(simulationOutput, asFactor=T)
testDispersion(simulationOutput)  

#models with different random strucutres
B_M2 <- lmer(sum_biomass ~ Eradication + (1|Plot),
             data=B_total)

B_M3 <- lmer(sum_biomass ~ Eradication + (1|Year),
             data=B_total)
#assessing AIC, it seems that 1|Plot is best
AICc(B_M1, B_M2, B_M3)

#null model to compare to with no fixed effects
B_Mnull <- lmer(sum_biomass ~ 1 + (1|Plot),
                data = B_total)

#assessing AIC values for best compared to null
AICc(B_M2, B_Mnull)
AICc(B_M1, B_M2, B_M3, B_Mnull)

#summary of best model, with residuals and marginal R^2 for fixed effects
summary(B_M2)
plot(residuals(B_M2))
r.squaredGLMM(B_M2)

#bootstrapped p-value from many distributions of the models
PB_B <- PBmodcomp(B_M2, C_Mnull, nsim = 10000)

#this is the average biomass per plot with and without rats
all_biomass %>%
  group_by(Plot, Rats, Year) %>%
  summarize(sum_biomass = sum(biomass_g, na.rm=T)) %>%
  group_by(Rats) %>%
  summarize(biomass = mean(sum_biomass, na.rm=T))

#mean biomass of individual seedling with and without rats
all_biomass %>%
  group_by(Rats) %>%
  summarize(biomass = mean(biomass_g, na.rm=T))

#This is looking at average total biomass of each plot with and without rats
seed_bio <- B_total %>%
  group_by(Plot, Rats) %>%
  summarize(bio_g = mean(sum_biomass, na.rm=T))

#change column name for later merging with adult and projected DFs
colnames(seed_bio)[which(names(seed_bio) == "Rats")] <- "Erad"

#Seedling nums for responses####
all_biomass$Species1 <- factor(all_biomass$Species1, levels = c("CN", "PG", "PF"))

#juvenile biomass by species for paper
(seed_sp <- ggplot(all_biomass, aes(x = Rats, y = (biomass_g/300)*0.01, fill = Species1)) +
    geom_boxplot() + theme_bw() +
    scale_fill_manual(values = c("#00A08A", "#F2AD00", "#5BBCD6")) +
    labs(x = "Rats", y = "Juvenile biomass (Mg/ha)") +
    theme(axis.text = element_text(size = 35), axis.title=element_text(size=40)))

counts_CN$Species <- "CN"

all_biomass %>%
  group_by(Plot, Year, Species1) %>%
  summarize(sum = sum(biomass_g, na.rm=T)) %>%
  group_by(Year, Species1) %>%
  summarize(mean = mean(sum, na.rm=T))

mean(c(10616, 12791, 12405)) #11937.33
mean(c(8880, 9570, 9644)) #9364.667
mean(c(6635, 7017, 3248)) #5633.333

11937.33/sum(11937.33, 9364.667, 5633.333)
(251040 - 11937.33)/(sum(251040, 43006, 59264) - sum(11937.33, 9364.667, 5633.333))
21/4
21/10.5
251040/(251040 + 43006+59264)
251040/11937.33
43006/9364.667
59264/5633.333
11937.33/251040
9364.667/43006
5633.333/59264

251040/11937.33
(251040/11937.33)/(43006/9364.667)
59264/5633.333
?mean
#summarize total CN in plots
CN_counts <- counts_CN %>%
  group_by(Plot, Year, Rats, Species) %>%
  summarize(n = sum(count, na.rm=T))
#summarize total number PG and PF in plots
PG_PF_counts <- PG_PF %>%
  group_by(Plot, Year, Rats, Species) %>%
  tally()
#bind into a dataframe
all_counts <- CN_counts %>% full_join(PG_PF_counts)
all_counts$Rats <- factor(all_counts$Rats, levels = c("pre", "post"))

#add in zero counts for plotting and averages for plots/years w/o specific species
Plot <- c("Holei 1", "Holei 1", "Holei 1", "Holei 1", "Holei 2", "Kaula", "Kaula", "Kaula",
          "Kaula", "Kaula", "Kaula", "Papala", "Papala", "Papala", "Papala", "Papala", 
          "Papala", "Papala", "Papala", "Paradise", "Sand", "Sand", "Sand", "Sand")
Year <- c(2008, 2010, 2016, 2016, 2016, 2008, 2008, 2009, 2009, 2016, 2016, 2008, 2008,
          2009, 2009, 2010, 2010, 2016, 2016, 2016, 2008, 2009, 2010, 2016)
Rats <- c("pre", "pre", "post", "post", "post", "pre", "pre", "pre", "pre", "post", "post",
          "pre", "pre", "pre", "pre", "pre", "pre", "post", "post", "post", "pre", "pre",
          "pre", "post")
Species <- c("PG","PG", "PF", "PG", "PG", "PF", "PG", "PF", "PG", "PF", "PG", 
             "PF", "PG", "PF", "PG", "PF", "PG", "PF", "PG", "PF", "PF", "PF", "PF", "PF")
n <- c(0)
zeros <- as.data.frame(cbind(Plot, Year, Rats, Species, n))

zeros$Plot <- as.character(zeros$Plot)
zeros$Year <- as.integer(as.character(zeros$Year))
zeros$Species <- as.character(zeros$Species)
zeros$n <- as.numeric(as.character(zeros$n))
all_counts <- all_counts %>% full_join(zeros)
all_counts$Rats <- factor(all_counts$Rats, levels = c("pre", "post"))

#abundance by species graph not in paper
(seed_sp_abund <- ggplot(all_counts, aes(x = Rats, y=n, fill = Species))+
  geom_boxplot() + theme_bw() +
  scale_fill_manual(values = c("#00A08A", "#F2AD00", "#5BBCD6")) +
  labs(x = "Rats", y = "Juvenile abundance") +
  theme(axis.text = element_text(size = 35), axis.title=element_text(size=40)))

#this is ttotal seedling abundance in each plot with and without rats
abund <- all_counts %>%
  group_by(Plot, Rats) %>%
  tally()

abund_sp <- all_counts %>%
  group_by(Plot, Rats, Species) %>%
  tally()

abund%>%
  group_by(Rats) %>%
  summarize(mean = mean(nn))

abund_cn <- abund_sp[which(abund_sp$Species == "CN"),]

#figure in paper for all juvenile abundance following eradication
(abund_1 <- ggplot(abund, aes(x = Rats, y =n)) +
  geom_boxplot(fill = "#00A08A")+ theme_bw() +
  labs(x = "Rats", y = "Juvenile abundance") +
    ylim(0,940) +
    annotate("text", x = 1.5, y = 850, label = "*", size = 22) +
  theme(axis.text = element_text(size = 35), axis.title=element_text(size=40)))

(abund_cn <- ggplot(abund_cn, aes(x = Rats, y =n)) +
    geom_boxplot(fill = "#00A08A")+ theme_bw() +
    labs(x = "Rats", y = "Juvenile abundance") +
    ylim(0,940) +
    annotate("text", x = 1.5, y = 850, label = "*", size = 22) +
    theme(axis.text = element_text(size = 35), axis.title=element_text(size=40)))
#2. Adult Biomass####
setwd("/Users/Ana/Dropbox/2017\ Palmyra\ Seed\ Predation\ Paper/Analyses/Draft_1_Mar_2019")

#from Chave et al. 2005
#AGB_est = p x exp(-1.239 + (1.980*ln(D)) + (0.207*(ln(D))^2) - 0.0281*(ln(D))^3) 
#PG: average wood density (p) for Pisonia genus given as: 0.32 g/cm^3 (n=6 studies, 4 species) FIND THESE AGAIN!! AHHH
#TA: there is one Tournefortia in this dataset, a value of 0.47 (CHave et al. 2006 databse)
#PF: wood density we measured was .207994673, from the literature there is 0.33 g/cm^3. average is .27 FIND LIT VALUE!!
##CN: average wood density for CN is 0.50 g/cm^3 (n = 3 studies) FIND STUDIES
setwd("/Users/Ana/Dropbox/2017\ Palmyra\ Seed\ Predation\ Paper/Analyses/Draft_1_Mar_2019/")
trees <- read.csv("Adult_DBH_for_Biomass.csv")

PG <- trees[which(trees$Species == "PG"),]
PG$Bio <- 0.32 * exp(-1.239 + (1.980*log(PG$DBH)) + (0.207*(log(PG$DBH))^2) - 0.0281*(log(PG$DBH))^3)

CN <- trees[which(trees$Species == "CN"),]
CN$Bio <- 0.5 * exp(-1.239 + (1.980*log(CN$DBH)) + (0.207*(log(CN$DBH))^2) - 0.0281*(log(CN$DBH))^3)

PF <- trees[which(trees$Species == "PF"),]
PF$Bio <- 0.27 * exp(-1.239 + (1.980*log(PF$DBH)) + (0.207*(log(PF$DBH))^2) - 0.0281*(log(PF$DBH))^3)

TA <- trees[which(trees$Species == "TA"),]
TA$Bio <- 0.47 * exp(-1.239 + (1.980*log(TA$DBH)) + (0.207*(log(TA$DBH))^2) - 0.0281*(log(TA$DBH))^3)

trees1 <- rbind(PG, CN, PF, TA)
#total bio in adults by plot
tree_bio <- trees1 %>%
  group_by(Plot) %>%
  summarize(bio_g = sum(Bio*1000, na.rm=T))
  
#total bio in adults by species in each plot
tree_bio_sp <- trees1 %>%
  group_by(Plot, Species) %>%
  summarize(bio_g = sum(Bio*1000, na.rm=T))
tree_bio_sp$Species <- factor(tree_bio_sp$Species, levels = c("CN", "PG", "PF", "TA"))

#plot of by species 
(adult_sp <- ggplot(tree_bio_sp, aes(x = Species, y = (bio_g/300)* 0.01)) +
  geom_boxplot(fill = c("#00A08A", "#F2AD00", "#5BBCD6", "#F98400")) +theme_bw() +
  labs(x = "Species", y = "Adult biomass (Mg/ha)") +
theme(axis.text = element_text(size = 35), axis.title=element_text(size=40)))

#3. Projected Added Biomass from Dem Models####
#per seedling and per new adult biomass values#
adult_G <- (54*(12^2.39)) - (54*(11^2.39))
sdlg_G <- 54
#stages for model
stages <- c("seed", "seedling", "adult")
####pre-eradication demographic model Matrix#
A_pre <- matrix(c(0, 0, 4.14, 0.34, 0.6, 0, 0, 0.03, 0.97), nrow = 3, byrow = TRUE, dimnames = list(stages, stages))  
#elasticities and lambda for pre model
eigen.analysis(A_pre)
#post-eradication demographic model matrix
A_post <- matrix(c(0, 0, 11.05, 0.34, 0.6, 0, 0, 0.03, 0.97), nrow = 3, byrow = TRUE, dimnames = list(stages, stages))
#elasticities and lambda for post model
eigen.analysis(A_post)

#population vectors per plot to populate models####
#n_hol1 <- c(0, 14.33, 30) #Holei 1
#n_pap <- c(13.48, 68.33, 24) #Papala
#n_par <- c(0, 15.25, 4) #Paradise
#n_eas <- c(67.41, 46.5, 3) #Eastern
#n_san <- c(0, 1.67, 17) #Sand
#n_kau <- c(0, 1.67, 23) #Kaula
#n_hol2 <- c(0, 16, 0) #Holei 2

#a list of vectors for each population to input into the for loop below
pops <- list(n_hol1=c(0,14.33,30), n_pap=c(13.48, 68.33, 24), n_par=c(0, 15.25, 4), n_eas=c(67.41, 46.5, 3),
             nsan=c(0, 1.67, 17), n_kau=c(0, 1.67, 23), n_hol2=c(0, 16, 0))

#Function for the future total biomass contribution over 12 years#
future_biomass <- function(left_matrix, pop_vector){ #a function of a Leftkovich matrix and a population vector
  pop_proj <- pop.projection(left_matrix, pop_vector, 12) #projects population out 12 years by stagge
  
  stage_vect <- pop_proj$stage.vectors #these are the population sizes per stage for the 12 year model
  
  stage_vect_t <- as.data.frame(t(stage_vect)) #transpose this dataframe
  
  stage_vect_s <- stage_vect_t %>%
    mutate(seed_diff = seedling - lag(seedling, default = first(seedling))) #find the number of new seedlings per time period
  
  stage_vect_s$gB_s <- stage_vect_s$seed_diff * sdlg_G #convert seedling number to seedling biomass in grams
  
  stage_vect_a <- stage_vect_t %>%
    mutate(adult_diff = adult-lag(adult, default = first(adult))) #find the number of new adults per time period
  
  stage_vect_a$gB_a <- stage_vect_a$adult_diff * adult_G #convert new adults to biomass
  
  stage_vect <- stage_vect_a %>% full_join(stage_vect_s) #create one dataframe of these values
  
  stage_vect$gB <- stage_vect$gB_a + stage_vect$gB_s #create a column in the data frame that is the per time period biomass addition (seedlings + adults)
  
  return(sum(stage_vect$gB)) #the total for the 12 years of new biomass in seedlings and adults
}

summary_matrix <- matrix(-999, ncol = 2, nrow = length(pops)) #creates an empty matrix to populate per-plot biomass additions in both pre and post time periods
for(i in 1:length(pops)){
  summary_matrix[i,1] = future_biomass(A_pre, pops[[i]]) #this creates a matrix column of biomass per plot with pre matrix
  summary_matrix[i,2] = future_biomass(A_post, pops[[i]]) #this creates a matrix column of biomass per plot with post matrix 
}
colnames(summary_matrix) <- c("pre", "post") #names those columns pre and post

#convert to DF for wilcoxon signed rank test
summary_df <- as.data.frame(summary_matrix)

#comparing added biomass with and without rats for each plot population
wilcox.test(summary_df$pre, summary_df$post, paired = TRUE)

#add a plot column to the summary_df
summary_df$Plot <- c("Holei 1", "Papala", "Paradise", "Eastern", "Sand", "Kaula", "Holei2")

proj_bio <- gather(summary_df,key = "Erad", value = bio_g, pre:post)
proj_bio$Type <- "Projected"

#4. Add all Biomass Together####
#seedling
#adult
#future
#adding new columns to tree dataframe for merging
tree_bio$Type <- "Adult"
tree_bio2 <- tree_bio
tree_bio <- as.data.frame(rbind(tree_bio, tree_bio2))
tree_bio$Plot <- as.character(tree_bio$Plot)
tree_bio$Erad <- c("pre", "pre", "pre", "pre", "pre", "pre", "pre", "post", "post", "post", "post", "post", "post", "post")

#adding new columns to seedling dataframe for merging
seed_bio$Type <- "Seedling"
seed_bio$Plot <- as.character(seed_bio$Plot)
seed_bio$Erad <- as.character(seed_bio$Erad)

#adding new columns to projected dataframe for merging
proj_bio$Plot <- as.character(proj_bio$Plot)
proj_bio$Erad <- as.character(proj_bio$Erad)
proj_bio$bio_g <- as.numeric(as.character(proj_bio$bio_g))

#merging all biomass componenets
all_bio <- tree_bio %>% full_join(seed_bio)
all_bio <- all_bio %>% full_join(proj_bio)

#making a sum of al of this biomass
sum_all_bio <- all_bio %>%
  group_by(Plot, Erad) %>%
  summarize(sum = sum(bio_g))

#make this a wide dataframe
sum_all_wide <- spread(sum_all_bio, Erad, sum)
#take out the vectors of pre and post
pre <- sum_all_wide$pre
post <- sum_all_wide$post

#comparing all biomass TOTAL pre-post
wilcox.test(pre, post, paired = TRUE)

#plot that
ggplot(all_bio, aes(x = Erad, y = bio_g, color = Type)) +
  geom_point()

#looking at the sums from different components of biomass
sum_bio <- all_bio %>%
  group_by(Erad, Type) %>%
  summarize(mean = mean(bio_g, na.rm=T))

#12227915.41 #adult
#(12227915.41/300)*0.01
#172028.15 #post proj
#284269.01 #post seed
#172028.15/12227915.41
#42671.75 #pre proj
#19180.50 #pre seed

sum_bio$Type <- factor(sum_bio$Type, levels = c("Projected", "Seedling", "Adult"))
ggplot(sum_bio, aes(x = Erad, y = mean, fill = Type)) +
  geom_bar(stat="identity")

#subsetting out adult biomass
sum_bio_new <- sum_bio[which(sum_bio$Type != "Adult"),]
#centering new biomass around zero
sum_bio_new$centered <- sum_bio_new$mean + 12227915.41

#looking at new bio added
ggplot(sum_bio_new, aes(x = Erad, y = (mean/300)*.01, fill = Type))+
  geom_bar(stat="identity") 

#spread all biomass
all_bio_wide <- spread(all_bio, Type, bio_g)
#combining current and future biomass
all_bio_wide$current <- all_bio_wide$Adult + all_bio_wide$Seedling
#percent more biomass based on current total
all_bio_wide$change <- (all_bio_wide$Projected / all_bio_wide$current)*100
#percent increase in biomass
all_bio_wide$change2 <- ((all_bio_wide$Projected + all_bio_wide$Seedling)/all_bio_wide$Adult)*100
all_bio_wide$Erad <- factor(all_bio_wide$Erad, levels = c("pre", "post"))
#percent of total biomass added now in projected
ggplot(all_bio_wide, aes(x = Erad, y = change)) +
  geom_boxplot() +
  labs(x = "Rats", y = "Percent Increase")
#percent increase in biomass 
ggplot(all_bio_wide, aes(x = Erad, y = (Projected/300)*.01)) +
  geom_boxplot() +
  labs(x = "Rats", y = "Total Increase (Mg/ha)")

#subset out only projected biomass
all_bio_proj <- all_bio[which(all_bio$Type == "Projected"),]

#just looking at mean projected additional biomass
mean_proj <- all_bio_proj %>%
  group_by(Erad) %>%
  summarize(mean = mean(bio_g))
#(172028/300)*0.01
#(42672/300)*0.01
#standard error for this
se_proj <- all_bio_proj %>%
  group_by(Erad) %>%
  summarize(se = (sd(bio_g)/sqrt(6)))

all_bio_proj %>%
  group_by(Erad) %>%
  summarize(sd = sd(bio_g)/30000)
            
all_bio_proj %>%
  group_by(Erad) %>%
  summarize(se = (sd(bio_g)/sqrt(6))/300*0.01)
#making a dataframe for bar graph
bar_proj <- mean_proj %>% 
  inner_join(se_proj)
bar_proj$Erad <- factor(bar_proj$Erad, levels = c("pre", "post"))

#figure for paper for projected future biomass
(d <- ggplot(bar_proj, aes(x = Erad, y = ((mean/300)*.01))) +
  geom_bar(stat="identity", fill = "#007867") +
  geom_errorbar(aes(ymin = ((mean/300)*.01)-((se/300)*.01), ymax = ((mean/300)*.01)+((se/300)*.01)), width = 0.1)+
  labs(x = "Rats", y = "Adult Biomass (Mg/ha)") +
  theme_bw() +
  theme(axis.text = element_text(size = 35), axis.title=element_text(size=40), legend.position="none")+
  annotate("text", x = 1.5, y = 6.75, label = "*", size = 22) +
    geom_hline(yintercept = 0, size = 1) +
    geom_hline(yintercept =-0.08, size = 1)+
  ylim(-0.08,7.6))

#Paper Figures####
#All juvenile biomass not by species#
(c <- ggplot(B_total, aes(x = Rats, y = (sum_biomass/300) *(.000001/1) * (1/0.0001), fill = Rats)) + 
   geom_boxplot(fill = "#00A08A") +   
   labs(x = "Rats", y = "Juvenile biomass (Mg/ha)") +
   theme_bw() + ylim(0,20) +
   theme(axis.text = element_text(size = 35), axis.title=element_text(size=40), legend.position="none")+
   annotate("text", x = 1.5, y = 17, label = "*", size = 22))


#juvenile biomass by species
(seed_sp <- ggplot(all_biomass, aes(x = Rats, y = (biomass_g/300)*0.01, fill = Species1)) +
    geom_boxplot() + theme_bw() +
    scale_fill_manual(values = c("#00A08A", "#F2AD00", "#5BBCD6")) +
    labs(x = "Rats", y = "Juvenile biomass (Mg/ha)") +
    theme(axis.text = element_text(size = 35), axis.title=element_text(size=40)))


#all juvenile abundance not by species
(abund_1 <- ggplot(abund, aes(x = Rats, y =n)) +
    geom_boxplot(fill = "#00A08A")+ theme_bw() +
    labs(x = "Rats", y = "Juvenile abundance") +
    ylim(0,940) +
    annotate("text", x = 1.5, y = 850, label = "*", size = 22) +
    theme(axis.text = element_text(size = 35), axis.title=element_text(size=40)))

#juv abund CN only
(abund_cn_only <- ggplot(abund_cn, aes(x = Rats, y =n/300)) +
    geom_boxplot(fill = "#00A08A")+ theme_bw() +
    labs(x = "Rats", y = "Juvenile abundance") +
    ylim(0,3) +
    theme(axis.text = element_text(size = 35), axis.title=element_text(size=40)))

abund_cn %>%
  group_by(Rats) %>%
  summarize(mean = mean(n/300))

#adult biomass by species
(adult_sp <- ggplot(tree_bio_sp, aes(x = Species, y = (bio_g/300)* 0.01)) +
    geom_boxplot(fill = c("#00A08A", "#F2AD00", "#5BBCD6", "#F98400")) +theme_bw() +
    labs(x = "Species", y = "Adult biomass (Mg/ha)") +
    theme(axis.text = element_text(size = 35), axis.title=element_text(size=40)))

#projected future biomass
(d <- ggplot(bar_proj, aes(x = Erad, y = ((mean/300)*.01))) +
    geom_bar(stat="identity", fill = "#007867") +
    geom_errorbar(aes(ymin = ((mean/300)*.01)-((se/300)*.01), ymax = ((mean/300)*.01)+((se/300)*.01)), width = 0.1)+
    labs(x = "Rats", y = "Adult Biomass (Mg/ha)") +
    theme_bw() +
    theme(axis.text = element_text(size = 35), axis.title=element_text(size=40), legend.position="none")+
    annotate("text", x = 1.5, y = 6.75, label = "*", size = 22) +
    geom_hline(yintercept = 0, size = 1) +
    geom_hline(yintercept =-0.08, size = 1)+
    ylim(-0.08,7.6))

#this is the plot grid for the paper, which includes the two seed values from the
#seedling_nut_lmes.R mixed effects model r script. considering blending for final
plot_grid(b,a,abund_1, d,nrow=1,align="h")
#(21536890.54/300)*.01

#this is the by species biomass figure in paper
plot_grid(adult_sp,seed_sp,nrow=1,align="h")

