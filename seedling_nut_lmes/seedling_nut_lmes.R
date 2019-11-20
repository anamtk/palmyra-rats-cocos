library(ggplot2)
library(gridExtra)
#library(glmmTMB)
library(DHARMa)
#library(RVAideMemoire)
library(MuMIn)
library(lme4)
library(pbkrtest)
library(parallel)
library(MASS)
library(tidyverse)
library(cowplot)
#Mac
#setwd("/Users/Ana/Dropbox/2017\ Palmyra\ Seed\ Predation\ Paper/Analyses/Draft-18-Nov-18/")

setwd("/Users/Ana/Dropbox/2017\ Palmyra\ Seed\ Predation\ Paper/Analyses/Draft_1_Mar_2019/seedling_nut_lmes")
#PC 
#setwd("C:/Users/kuile/Dropbox/2017 Palmyra Seed Predation Paper/Analyses/Draft-18-Nov-18")

####NUT MODEL####
quad.nuts <- read.csv("Quadrat_Nuts_no_H2.csv")#import data minus Holei 2 since it has no adults

#manipulate factors
quad.nuts$predated <- quad.nuts$Immature.Predated + quad.nuts$Mature.Predated
quad.nuts$Rats <- ifelse(quad.nuts$Eradication=='pre', 1, 2)
quad.nuts$Year <- as.factor(quad.nuts$Year)

ggplot(quad.nuts, aes(x = Year, y = predated)) +
  geom_jitter()
ggplot(quad.nuts, aes(x = Year, y = Mature.Intact)) +
  geom_jitter()
ggplot(quad.nuts, aes(x = Year, y = Immature.Predated)) +
  geom_jitter()
ggplot(quad.nuts, aes(x = Year, y = Mature.Predated)) +
  geom_jitter()

#Proportion Predated in Trees####
quad.nuts$tot.nuts <- quad.nuts$Mature.Intact + quad.nuts$Immature.Intact + quad.nuts$Immature.Rotten +
                          quad.nuts$Mature.Rotten + quad.nuts$Immature.Predated + quad.nuts$Mature.Predated
quad.nuts$rat.prop <- quad.nuts$Immature.Predated/quad.nuts$tot.nuts
quad.nuts$rat.prop[is.na(quad.nuts$rat.prop)] <- 0

M.prop <- glmer(rat.prop ~ Rats + (1|Site),
                data = quad.nuts,
                family = "binomial")

summary(M.prop)
plot(residuals(M.prop))
library(effects)
plot(allEffects(M.prop))
library(lattice)
dotplot(ranef(M.prop, condVar = T))

simulationOutput <- simulateResiduals(fittedModel = M.prop)
plot(simulationOutput, asFactor = T)
testZeroInflation(simulationOutput)
testDispersion(simulationOutput)

M.prop.nrand <- glm(rat.prop ~ Rats,
                    data= quad.nuts,
                    family = "binomial")

AICc(M.prop, M.prop.nrand)

M.prop.null <- glmer(rat.prop ~ (1|Site),
                     data= quad.nuts,
                     family = "binomial")

AICc(M.prop, M.prop.null)
PBmodcomp(M.prop, M.prop.null, nsim = 100)
anova(M.prop, M.prop.null)
summary(M.prop)
plot(residuals(M.prop))
r.squaredGLMM(M.prop)

sim.prop = simulate(M.prop.null, nsim = 100) #simulations that simulate dists. of the null model. check nsim, but I think 10,000 is appropriate
LR = vector() #empty vector where the likelihood ratios simulated will go
for (i in 1:100) { #this simulates likelihood given both models
  try(
    {
      M.full = glmer(sim.prop[[i]] ~ Rats + (1|Site),
                        data=quad.nuts,
                        family = "binomial")
      M.res = glmer(sim.prop[[i]] ~ 1 + (1|Site),
                       data=quad.nuts,
                    family = "binomial")
    })
  if((!is.character(M.full))&(!is.character(M.res))){
    LR[i] = 2*(logLik(M.full) - logLik(M.res))
  }else{
    LR[i] = NA #this section was because of my weird distribution spitting errors in the first bit. 
  }
  if((i%%100) == 0){print(paste("Iteration ",i))} #this gives status
}
LR <- LR[!is.na(LR)] #removes NA values from likelihood ratios created by the fix in the "else" statement above
hist(LR) #gives a histogram of likelihood ratios
p <- (sum(LR > 2*(logLik(M.prop) - logLik(M.prop.null))))/length(LR) #this gives a p-value of the full vs. null model likelihood
p #p-value
#imprecision in bootstrap p-value
#2 x sqrt(pv_est * (1 - pv_est) / N)
c <- 2*sqrt(p*(1-p)/100) #this is the confidence interval around that. note that the /1000 should be fixed to be number of simulations
c

means <- quad.nuts %>%
  group_by(Rats) %>%
  summarize(mean = mean(rat.prop, na.rm=T))

se <- quad.nuts %>%
  group_by(Rats) %>%
  summarize(se = sd(rat.prop, na.rm=T)/sqrt(nrow(quad.nuts)))
se <- se$se

seed_graph <- as.data.frame(cbind(means,se))
str(seed_graph)
seed_graph$Erad <- ifelse(seed_graph$Rats==1, "pre", "post")
seed_graph$Erad <- factor(seed_graph$Erad, levels = c("pre", "post"))
(a <- ggplot(seed_graph, aes(x = Erad, y = mean)) + 
    geom_crossbar(aes(x=Erad, ymin=mean-se, ymax=mean+se, fill = Erad)) +
    scale_fill_manual(values = c("#00A08A", "#00A08A"))+ 
    labs(x = "Rats", y = "Proportion predated") +
    theme_bw() +
    theme(axis.text = element_text(size = 35), axis.title=element_text(size=40), legend.position="none") +
    annotate("text", x = 1.5, y = 0.2, label = "***", size = 22) +
    ylim(0,0.3))

#Mature Intact ####
#a model to test overdispersion and zero inflation
M1 <- glmer(Mature.Intact ~ Rats + (1|Site),
              data=quad.nuts, 
              family=poisson)

simulationOutput <- simulateResiduals(fittedModel = M1) #simulate residuals, look weird, switched to nb
plot(simulationOutput, asFactor=T) #plot these
testZeroInflation(simulationOutput) #test zero inflation
testDispersion(simulationOutput) #check that it is not overdispersed

M1.nb <- glmer.nb(Mature.Intact ~ Rats + (1|Site),
                  data=quad.nuts)

simulationOutput <- simulateResiduals(fittedModel = M1.nb) #simulate residuals
plot(simulationOutput, asFactor=T) #plot these
testZeroInflation(simulationOutput) #test zero inflation
testDispersion(simulationOutput) #check that it is not overdispersed

AICc(M1, M1.nb)

M2 <- glm.nb(Mature.Intact ~ Rats,
            data=quad.nuts)

AICc(M1, M1.nb, M2) #testing random structure, where M1 clearly performs better

M3 <- glmer.nb(Mature.Intact ~ (1|Site),
            data=quad.nuts)

AICc(M1, M1.nb, M2, M3) #ccomparing fixed effect set up

#summarizing best model
summary(M1.nb)
plot(residuals(M1.nb))
r.squaredGLMM(M1.nb)
(nc <- detectCores())
cl <- makeCluster(rep("localhost", nc))

refdist <- PBrefdist(M1.nb, M3, nsim = 10000, cl=cl)
pb_nuts <- PBmodcomp(M1.nb, M3, ref=refdist)

#stop parallel processing:
stopCluster(cl)

means <- quad.nuts %>%
  group_by(Rats) %>%
  summarize(mean = mean(Mature.Intact, na.rm=T))

se <- quad.nuts %>%
  group_by(Rats) %>%
  summarize(se = sd(Mature.Intact, na.rm=T)/sqrt(nrow(quad.nuts)))
se <- se$se

seed_graph <- as.data.frame(cbind(means,se))
str(seed_graph)
seed_graph$Erad <- ifelse(seed_graph$Rats==1, "pre", "post")
seed_graph$Erad <- factor(seed_graph$Erad, levels = c("pre", "post"))
(a <- ggplot(seed_graph, aes(x = Erad, y = mean)) + 
  geom_crossbar(aes(x=Erad, ymin=mean-se, ymax=mean+se, fill = Erad)) +
  scale_fill_manual(values = c("#00A08A", "#00A08A"))+ 
  labs(x = "Rats", y = "Seed abundance (seeds/square meter)") +
  theme_bw() +
  theme(axis.text = element_text(size = 35), axis.title=element_text(size=40), legend.position="none") +
  annotate("text", x = 1.5, y = 0.4, label = "***", size = 22) +
  ylim(0,0.45))

#Total Predated ####
#a model to test overdispersion and zero inflation
M1 <- glmer(predated ~ Rats + (1|Site),
            data=quad.nuts, 
            family=poisson)

simulationOutput <- simulateResiduals(fittedModel = M1) #simulate residuals, look good
plot(simulationOutput, asFactor=T) #plot these
testZeroInflation(simulationOutput) #test zero inflation
testDispersion(simulationOutput) #check that it is not overdispersed

M2 <- glm(predated ~ Rats,
             data=quad.nuts)

AICc(M1, M2) #testing random structure, where M1 clearly performs better

M3 <- glmer(predated ~ (1|Site),
               data=quad.nuts,
            family = poisson)

AICc(M1, M2, M3) #ccomparing fixed effect set up

#summarizing best model
summary(M1)
plot(residuals(M1))
r.squaredGLMM(M1)
(nc <- detectCores())
cl <- makeCluster(rep("localhost", nc))

refdist <- PBrefdist(M1, M3, nsim = 10000, cl=cl)
pb_nuts <- PBmodcomp(M1, M3, ref=refdist)
pb_nuts
#stop parallel processing:
stopCluster(cl)

means <- quad.nuts %>%
  group_by(Rats) %>%
  summarize(mean = mean(predated, na.rm=T))

se <- quad.nuts %>%
  group_by(Rats) %>%
  summarize(se = sd(predated, na.rm=T)/sqrt(nrow(quad.nuts)))
se <- se$se

seed_graph <- as.data.frame(cbind(means,se))
str(seed_graph)
seed_graph$Erad <- ifelse(seed_graph$Rats==1, "pre", "post")
seed_graph$Erad <- factor(seed_graph$Erad, levels = c("pre", "post"))
(b <- ggplot(seed_graph, aes(x = Erad, y = mean)) + 
    geom_crossbar(aes(x=Erad, ymin=mean-se, ymax=mean+se, fill = Erad)) +
    scale_fill_manual(values = c("#F2AD00", "#F2AD00"))+ 
    labs(x = "Rats", y = "Seed abundance (seeds/square meter)") +
    theme_bw() +
    theme(axis.text = element_text(size = 35), axis.title=element_text(size=40), legend.position="none") +
    annotate("text", x = 1.5, y = 0.4, label = "***", size = 22) +
    ylim(0,0.45))

plot_grid(b,a,nrow=1,align="h")

#Pre-dispersal predation####
M1 <- glmer(Immature.Predated ~ Rats + (1|Site),
            data=quad.nuts, 
            family=poisson)
?glmer
simulationOutput <- simulateResiduals(fittedModel = M1) #simulate residuals, look good
plot(simulationOutput, asFactor=T) #plot these
testZeroInflation(simulationOutput) #test zero inflation
testDispersion(simulationOutput) #check that it is not overdispersed

M2 <- glm(Immature.Predated ~ Rats,
          data=quad.nuts)

AICc(M1, M2) #testing random structure, where M1 clearly performs better

M3 <- glmer(Immature.Predated ~ (1|Site),
            data=quad.nuts,
            family = poisson)

AICc(M1, M2, M3) #ccomparing fixed effect set up

#summarizing best model
summary(M1)
plot(residuals(M1))
r.squaredGLMM(M1)
(nc <- detectCores())
cl <- makeCluster(rep("localhost", nc))

refdist <- PBrefdist(M1, M3, nsim = 10000, cl=cl)
pb_nuts <- PBmodcomp(M1, M3, ref=refdist)
pb_nuts
#stop parallel processing:
stopCluster(cl)

means <- quad.nuts %>%
  group_by(Rats) %>%
  summarize(mean = mean(Immature.Predated, na.rm=T))

quad.nuts %>%
  group_by(Eradication) %>%
  summarize(mean = mean(Immature.Predated, na.rm=T))

quad.nuts %>%
  group_by(Eradication) %>%
  summarize(sd = sd(Immature.Predated, na.rm=T))

se <- quad.nuts %>%
  group_by(Rats) %>%
  summarize(se = sd(Immature.Predated, na.rm=T)/sqrt(nrow(quad.nuts)))
se <- se$se

((.306-.0139)/.306)*100
seed_graph <- as.data.frame(cbind(means,se))
str(seed_graph)
seed_graph$Erad <- ifelse(seed_graph$Rats==1, "pre", "post")
seed_graph$Erad <- factor(seed_graph$Erad, levels = c("pre", "post"))
(b <- ggplot(seed_graph, aes(x = Erad, y = mean)) + 
    geom_crossbar(aes(x=Erad, ymin=mean-se, ymax=mean+se, fill = Erad)) +
    scale_fill_manual(values = c("#F2AD00", "#F2AD00"))+ 
    labs(x = "Rats", y = "Seed abundance (seeds/square meter)") +
    theme_bw() +
    theme(axis.text = element_text(size = 35), axis.title=element_text(size=40), legend.position="none") +
    annotate("text", x = 1.5, y = 0.4, label = "***", size = 22) +
    ylim(0,0.45))

plot_grid(b,a,nrow=1,align="h")

quad.nuts%>%
  group_by(Year) %>%
  summarize(mean = mean)
#Mature Predated####
M1 <- glmer(Mature.Predated ~ Rats + (1|Site),
            data=quad.nuts, 
            family="poisson", 
            control = lmerControl(optimizer ="Nelder_Mead"))

simulationOutput <- simulateResiduals(fittedModel = M1) #simulate residuals, look good
plot(simulationOutput, asFactor=T) #plot these
testZeroInflation(simulationOutput) #test zero inflation
testDispersion(simulationOutput) #check that it is not overdispersed

quad.nuts %>%
  group_by(Year) %>%
  summarize(sum = sum(Mature.Predated))
quad.nuts%>%
  group_by(Year) %>%
  tally()

quad.nuts %>%
  group_by(Year) %>%
  summarize(sum = sum(Immature.Predated))
########SUPPLEMENTARY SEEDLING MODEL#####

seedlings <- read.csv("2007-2017\ Seedlings_all_species.csv") #seedlings of all species counts

#manipulate data
seedlings$CN_Count <- as.numeric(as.character(seedlings$CN_Count))
seedlings$PF_Count <- as.numeric(as.character(seedlings$PF_Count))
seedlings$PG_Count <- as.numeric(as.character(seedlings$PG_Count))
seedlings$Year <- as.factor(seedlings$Year)
seedlings$Eradication <- ifelse(seedlings$Rats=='pre', 1, 2)
seedlings$total_seedlings <- seedlings$CN_Count + seedlings$PF_Count + seedlings$PG_Count

ggplot(seedlings, aes(x = Year, y = total_seedlings)) +
  geom_point()

seedlings %>%
  group_by(Rats) %>%
  summarize(mean= mean(CN_Count, na.rm=T))

seedlings_full <- seedlings[which(seedlings$total_seedlings != "NA") ,]
str(seedlings_full)
#full model for DHARMa diagnostics, which are the same as above
M1 <- glmer(total_seedlings ~ Eradication + (1|Year/Plot),
               data=seedlings_full, 
               family=poisson)

simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput, asFactor=T)
testDispersion(simulationOutput)

M1 <- glmer(total_seedlings ~ Eradication + (1|Year/Plot),
               data=seedlings_full, 
               family=poisson)

M2 <- glmer(total_seedlings ~ Eradication + (1|Year),
               data=seedlings_full, 
               family=poisson)

M3 <- glmer(total_seedlings ~ Eradication + (1|Plot),
               data=seedlings_full, 
               family=poisson)

AICc(M1,M2,M3) #compare randoms, with M1 outperforming

#null to compare these to for fixed effects
Mnull <- glmer(total_seedlings ~ 1 + (1|Year/Plot),
               data=seedlings_full, 
               family=poisson)

AIC(M1,M2, M3, Mnull) #compare

summary(M1)
plot(residuals(M1))
r.squaredGLMM(M1)

(nc <- detectCores())
cl <- makeCluster(rep("localhost", nc))

refdist <- PBrefdist(M1, Mnull, nsim = 10000, cl=cl)
refdist
pb_sdlg <- PBmodcomp(M1, Mnull, ref=refdist)
pb_sdlg
#stop parallel processing:
stopCluster(cl)

seedlings_full %>%
  group_by(Eradication) %>%
  summarize(mean = mean(total_seedlings, na.rm=T))

#Seedling repeated measures ANOVA####

#Parametric##
##resource: https://www.r-bloggers.com/beware-the-friedman-test/ and
#the book "Serious Stats"
#this is just another way to call a repeated measures ANOVA. we are asking, do
#total seedlings vary by year, when grouped by the repeated measure of Plot
library(nlme)
library(emmeans)
library(rcompanion)
lme.raw <- lme(fixed = total_seedlings ~ Year, random =~1|Plot, 
               data=seedlings_full)

lme.raw.null <- lme(fixed = total_seedlings ~ 1, random =~1|Plot, 
                    data=seedlings_full)
??nagelkerke
nagelkerke(lme.raw,
           lme.raw.null)

model.emm <- emmeans(lme.raw, "Year")

#pair-wise comparisons of all years
pairs(model.emm)
#from this, it appears that the p-values follow the same trend as above, but
#are higher, suggesting perhaps that this approach is more conservative?

x <- residuals(lme.raw)
plotNormalHistogram(x)
#WOW, NOT normal!

plot(fitted(lme.raw),
     residuals(lme.raw))
abline(0,0)
#hahahahah this is terrible!!! Let us try rank transformation of this

#Rank-Transformed####
seedlings_full$r.Seeds <- rank(seedlings_full$total_seedlings)

lme.rank <- lme(fixed = r.Seeds ~ Year, random =~1|Plot, 
                data=seedlings_full)

lme.rank.null <- lme(fixed = r.Seeds ~ 1, random =~1|Plot, 
                     data=seedlings_full)

nagelkerke(lme.rank, 
           lme.rank.null)

emodel.emm <- emmeans(lme.rank, "Year")

#pair-wise comparisons of all years
pairs(model.emm)

model.pairs <- pairs(model.emm)
model_table <- write.table("model.pairs")
marginal = emmeans(lme.rank, 
                   ~ Year)

CLD(marginal,
    alpha   = 0.05, 
    Letters = letters,     ### Use lower-case letters for .group
    adjust  = "tukey")

x <- residuals(lme.rank)
plotNormalHistogram(x)
#these residuals look normal!

#check for homoskedastictiy
plot(fitted(lme.rank),
     residuals(lme.rank))
abline(0,0)
#these also look good!

#Precip####
precip <- read.csv("Palmyra\ Annual\ Precipitation\ 2002-2016.csv") 
#precip <- precip[1:15,]
precip <- precip[6:15,]
str(precip)
precip$Year <- as.numeric(as.character(precip$Year))
precip$Eradication <- ifelse(precip$Year <2011, "pre", "post")

precip %>%
  group_by(Eradication) %>%
  summarize(mean = mean(total.cm))

precip.pre <- precip[which(precip$Eradication == "pre"),]
precip.post <- precip[which(precip$Eradication == "post"),]
wilcox.test(precip.pre$total.cm, precip.post$total.cm, paired = F)
