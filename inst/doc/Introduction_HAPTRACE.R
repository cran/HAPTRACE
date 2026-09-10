## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(HAPTRACE)
set.seed(1)

## -----------------------------------------------------------------------------
MAPfile <- generateMAP(nChr = 10, n_markers = 1000, len_Chr = 100)
QTLfile <- QTLeffects(n_QTL_Chr = 10, nChr = 10, len_Chr = 100, effect_distribution = "gamma")

MAP_QTL_Pop <- mapQTLfile(map = MAPfile, QTL = QTLfile)
Effect <- as.vector(as.numeric(MAP_QTL_Pop$Effect))


## -----------------------------------------------------------------------------
Historical_Population <- generateHP(n_generations = 5, n_animals = 2000, 
                                   map = MAP_QTL_Pop, nSire = 200, nDam = 500, 
                                   nProgeny = 1000, mutationRate = 2.5 * 10^-5,
                                   SelType = "random", Effect = Effect, nChr = 10, 
                                   h2 = 0.3, trait_mean = 40, VarE = 0.6,recL = 10, 
                                   minLength = 200)

Sim_HP <- Historical_Population$Haplotype
str(Sim_HP)
class(Sim_HP)
Sim_HP[1:10, 1:10]

## -----------------------------------------------------------------------------
Population_A <- Sim_HP[1:(nrow(Sim_HP)/2),]
Population_B <- Sim_HP[(nrow(Sim_HP)/2 + 1):nrow(Sim_HP), ]

EffectA <- Rescale(Haplotype = Population_A, Effect = Effect, Phenosd = 3, h2 = 0.3)
Effect_Sim <- EffectA$QTLeffect

## -----------------------------------------------------------------------------
# Assignment of population code
Population_A <- PopulationID(Population = Population_A, PopCode = "P1")
Population_B <- PopulationID(Population = Population_B, PopCode = "P2")

# Recode haplotype of populations A and B with 1 and 2, respectively
Population_A <- Recode_Haplotype(df = Population_A, popBase = 1)
Population_B <- Recode_Haplotype(df = Population_B, popBase = 2)

## -----------------------------------------------------------------------------
# Population A
popA <- GenOnePopulation(ngenerations = 1, map = MAP_QTL_Pop,
                         populationSim = Population_A, nSire = 20, nDam = 100,
                         mutationRate = 2.5 * 10^-5, SelType = "random", 
                         h2 = 0.3, trait_mean = 10, VarE = 0.6,Effect = Effect_Sim, 
                         recL = 10, nChr = 10, minLength = 200, IndStartVal=1,
                         prefixID = "A", nProgeny = 300)

# Population B
popB <- GenOnePopulation(ngenerations = 1, map = MAP_QTL_Pop,
                         populationSim = Population_B, nSire = 10, nDam = 100,
                         mutationRate = 2.5 * 10^-5, SelType = "random", 
                         h2 = 0.3, trait_mean = 30, VarE = 0.6 ,Effect = Effect_Sim,
                         recL = 10,nChr = 10, minLength = 200, IndStartVal=1,
                         prefixID = "B", nProgeny = 300)

# Simulation of a crossbred population derived from populations A and B

## -----------------------------------------------------------------------------
# First cross between populations A and B
F1Cross <- TwoPopCross(map = MAP_QTL_Pop, populationSim_A = popA$Haplotype,
                       populationSim_B = popB$Haplotype, nSireA = 20, nDamA = 100,
                       nSireB = 10, nDamB = 100, Sire_From = "A",
                       mutationRate = 2.5 * 10^-5, SelType = "random", h2 = 0.3,
                       trait_mean = 20,VarE_PopA = 0.6, VarE_PopB = 0.6, 
                       Effect_PopA = Effect_Sim, Effect_PopB = Effect_Sim,
                       IndStartVal = 1, prefixID = "F1C", nProgeny = 300, recL = 10, 
                       nChr = 10, minLength = 200)

## -----------------------------------------------------------------------------
# Simulation of an admixed population
F3 <- GenOnePopulation(ngenerations = 2, map = MAP_QTL_Pop,
                       populationSim = F1Cross$Haplotype,
                       nSire = 20,nDam = 100,mutationRate = 2.5 * 10^-5,
                       SelType = "random",h2 = 0.3, trait_mean = 20,
                       VarE = 0.6, Effect = Effect_Sim, recL = 10, 
                       nChr = 10, minLength = 200, IndStartVal = 1,
                       prefixID = "C2", nProgeny = 300)

## -----------------------------------------------------------------------------
# Extracting pedigree and phenotype
F3Pheno = Ext_Pheno(df = F3$population)
head(F3Pheno)
F3Ped_Pheno = Ped_Pheno(Pedigree = F3$parents, Phenotype = F3Pheno)
head(F3Ped_Pheno)

## -----------------------------------------------------------------------------
# Extracting genotype information for F1Cross and F3 populations
CBgenotype <- rbind(MakeGeno(Coded_to_Haplo(F1Cross$Haplotype)),
                    MakeGeno(Coded_to_Haplo(F3$Haplotype)))
CBgenotype[1:10, 1:10]
# Coded haplotype for F1Cross and F3 populations
CBHaplotype <- rbind(F1Cross$Haplotype,F3$Haplotype)
CBHaplotype[1:10, 1:10]

## ----fig.cap=paste("Figure 1: Haplotype plot for crossbred animals obtained from `plotHaplo` function"), fig.width = 6, fig.height = 6----
# Haplotype plot
smallF3 <- F3$Haplotype[1:100, ]
F3pAplot <- plotHaplo(Haplotype = smallF3, 
                     colors = c("#CC79A7", "#0072B2"), 
                     map = MAP_QTL_Pop, nChr = 10, 
                     legendlabels = c("A", "B"))


## ----fig.cap=paste("Figure 2: Admixture plot for crossbred animals obtained from `Admixplot` function"), fig.width = 6, fig.height = 6----
# Admixture plot
TBC <- trueBreedComp(Haplotype = smallF3)

# Arrange the population
TBC_data <- TBC[c(2,1), ]

Admixplot(trueBC = TBC_data, 
          color = c("#CC79A7", "#0072B2"),
          legendlabels = c("A", "B"))

