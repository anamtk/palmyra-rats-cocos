library(ggplot2) #plotting capability
library(gridExtra)
library(rcompanion) #for the cldList function to do posthoc dunn test for kruskal wallis
library(nlme)
library(car)
library(multcompView)
library(emmeans)
library(tidyverse)
library(lme4)
setwd('/Users/Ana/Dropbox/2017\ Palmyra\ Seed\ Predation\ Paper/Analyses/Draft-18-Nov-18/')
#Nut to seed transition####
nut_2_seed <- read.csv('nut_to_seedling_transition.csv')

nut_2_seed$Year <- as.factor(nut_2_seed$Year)
nut_2_seed <- nut_2_seed[which(nut_2_seed$seed_seedling_survival != 'NA'),]

lme <- lme(fixed = seed_seedling_survival ~ Year, random =~1|Plot, 
               data=nut_2_seed)

lme.null <- lme(fixed = seed_seedling_survival ~ 1, random =~1|Plot, 
                    data=nut_2_seed)

nagelkerke(lme,
           lme.null)
#from this, it appears that there is no significant by-year effect of seed to seedling transition. 

x <- residuals(lme)
plotNormalHistogram(x)


plot(fitted(lme),
     residuals(lme))
abline(0,0)

#this looks ok - a little odd, but not surprising given the small number of samples. 

#what is value?
nut_2_seed %>%
  summarize(mean = mean(seed_seedling_survival, na.rm=T))

#Seedling Survival####
seedling_sur <- read.csv("seedling_survival.csv")  
str(seedling_sur)
seedling_sur$Year <- as.factor(seedling_sur$Year)
seedling_sur <- seedling_sur[which(seedling_sur$Survival_Rate != 'NA'),]
str(seedling_sur)
lme <- lme(fixed =Survival_Rate ~ Year, random =~1|Plot, 
           data=seedling_sur)

lme.null <- lme(fixed = Survival_Rate ~ 1, random =~1|Plot, 
                data=seedling_sur)

nagelkerke(lme,
           lme.null)
#from this, it appears that there is no significant by-year effect of seed to seedling transition. 

x <- residuals(lme)
plotNormalHistogram(x)

plot(fitted(lme),
     residuals(lme))
abline(0,0)

seedling_sur %>%
  group_by(Plot, Year) %>%
  summarize(mean = mean(Survival_Rate, na.rm=T))

#Seedling Mortality and size Distributions####
#we needed to determine annual litter inputs, which required knowing whether
#1. seedling abundance changed (lmes) 2. seedling mortality changed
#(based on survival, the answer here is no) and 3. whether size distributions
#of seedlings changed, thus driving an increase in the per seedling 
#contribution to biomass and thus litter.

#we impoted data on seedling size distributions to see if they changed pre-post
counts_CN <- read.csv("seedling_size_dists.csv")
str(counts_CN)

counts_CN$total <- counts_CN$CN_Yr1 + counts_CN$CN_Yr3 + counts_CN$CN_Yr5
counts_CN$bio_1 <- counts_CN$CN_Yr1 * (54)
counts_CN$bio_3 <- counts_CN$CN_Yr3 * (54*3^2.39)
counts_CN$bio_5 <- counts_CN$CN_Yr5 * (54*5^2.39)
counts_CN$total_bio <- counts_CN$bio_1 + counts_CN$bio_3 + counts_CN$bio_5

#per seedling biomass contribution
counts_CN$per_seed_bio <- counts_CN$total_bio/counts_CN$total

counts_CN <- counts_CN[which(counts_CN$per_seed_bio != "NaN"),]
lme <- lme(fixed =per_seed_bio ~ Year, random =~1|Plot, 
           data=counts_CN)

lme.null <- lme(fixed = per_seed_bio ~ 1, random =~1|Plot, 
                data=counts_CN)

nagelkerke(lme,
           lme.null)

counts_CN %>%
  summarize(per_seed = mean(per_seed_bio, na.rm=T))
#seems that per seedling contribution did not significantly change,
#therefore, litter inputs are just the average per seedling
#biomass contribution * total number of seedlings dead

seedling_sur$seed_litter <- seedling_sur$Prev_Year_Dead * 563.8
seedling_sur%>% 
  group_by(Erad) %>%
  summarize(mean = mean(seed_litter, na.rm=T))
str(seedling_sur)

M1 <- lmer(seed_litter ~ Erad + (1|Plot),
                  data=seedling_sur)

Mnull <- lmer(seed_litter ~ 1 + (1|Plot),
           data=seedling_sur)

AIC(M1,Mnull)

summary(M1)
plot(residuals(M1))
r.squaredGLMM(M1)

PBmodcomp(M1, Mnull, nsim=10000)

seedling_sur$Erad <- factor(seedling_sur$Erad, levels = c("pre", "post"))
(l <- ggplot(seedling_sur, aes(x = Erad, y = ((seed_litter/300)*(.000001/1) * (1/0.0001)))) +
  geom_boxplot(fill = c("#d6604d", "#878787"), size = 0.6) + 
  labs(y = "Seedling litter deposition (t/ha/yr)")  + theme_bw() +
  theme(axis.text = element_text(size = 35), 
        axis.title = element_text(size = 40))+
    ylim(0,7) +
    annotate("text", x = 1.5, y = 6, label = "**", size = 22))

plot_grid(l,c,nrow=2,align="v")
seedling_sur %>%
  group_by(Erad) %>%
  summarize(mean = mean((seed_litter/300)*(.000001/1) * (1/0.0001), na.rm=T))
2.65/0.214
((79393-6407)/6407)*100
79393/6407
120/90
((120-90)/90)*100
grid.arrange(c, l, nrow=1)
str(seedling_sur)
#Seedling to Adult####
sdl_adult <- read.csv("sdlg_adult.csv")
str(sdl_adult)
sdl_adult$Year <- as.factor(sdl_adult$Year)

myvars <- c("Plot", "Eradication", "survivorship")
sdl_adult <- sdl_adult[myvars]
  
wide_sdl_adult <- spread(sdl_adult, Eradication, survivorship)

wilcox.test(wide_sdl_adult$pre, wide_sdl_adult$post, paired = T)

#no significant difference

sdl_adult %>%
  summarize(mean = mean(survivorship, na.rm=T))
