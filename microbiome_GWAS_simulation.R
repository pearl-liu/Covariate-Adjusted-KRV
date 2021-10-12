#!/usr/local/bin/Rscript

## Investigate type I error and power in covariate-adjusted KRV tests 
## Simulate realistic genetic and microbiome data, and introduce population structure as a confounder 

#ptm <- proc.time()  # Time the simulation

#### Load libraries and data sets ####
library(argparse)
library(MiRKAT)
library(GUniFrac) 
library(SKAT)
library(dirmult)
library(doBy)
library(svd)


source('functions.R')

data("throat.meta")
data("throat.otu.tab")
data("throat.tree")

# Estimated parameters for Dirichlet-multinomial distribution based on the Charleson data set
DM_parameter_est = readRDS("throat.otu.tab.DM.parameter.estimates.rds")

# simulated SNP data over 1Mb using cosi2 program
eur_afr_haplo_mt = readRDS("10000eur_10000afr_1Mb_haplotype_matrix.rds")
eur_afr_haplo_snp_pos = readRDS("10000eur_10000afr_1Mb_haplotype_SNP_position.rds")
eur_afr_haplo_maf = readRDS("10000eur_10000afr_1Mb_haplotype_SNP_MAF.rds")

# labels of OTUs belonging to a common phylogenetic cluster
com_cluster = readRDS('throat_OTU_labels_node1148.rds')

eur_haplo_mt <- eur_afr_haplo_mt[1:10000,]
afr_haplo_mt <- eur_afr_haplo_mt[10001:20000,]

rm(eur_afr_haplo_mt)
gc()


#### Simulation overview ####
## Simulate abundance data with 856 OTUs using estimated parameters,
## assuming 1000 total counts per sample,
## for n individuals.

## Randomly pair haplotypes from 10,000 haplotypes to form the genotypes of the n individuals. 
## Test whether the common SNPs (MAF >= 0.05) within a 8kb subregion 
## of the 1 Mb region is associated with microbiome composition. 

## Use Bray-Kurtis and different UniFrac kernels for microbiome data
## Use linear kernel for genotype data



#### Load command line arguments  ####
parser <- ArgumentParser()
parser$add_argument("--n", type = "double",
                    help = "sample size")
parser$add_argument("--n_sim", type = "double",
                    help = "number of simulations")
parser$add_argument("--n_pop", type = "double",
                    help = "number of populations for genetic data")
parser$add_argument("--seed_num", type = "double",
                    help = "seed number for simulation")
parser$add_argument("--output_file", type = "character",
                    help = "output file name")
parser$add_argument("--power", type = "double",
                    help = "indicator for power(1) vs. type I error(0) analysis")
parser$add_argument("--power_situation", type = "double",
                    help = "specify power situation(0,1,2,3)")

args <- parser$parse_args()


n <- args$n # sample size
n_sim <- args$n_sim # number of simulations
n_pop = args$n_pop # specify the number of populations for genetic data
seed_num = args$seed_num # seed number for simulation
power_ind <- args$power # 1: power analysis, 0: type I error analysis
power_situ <-args$power_situation # 0 (type I error analysis), 1, 2 or 3

#### Baseline parameters ####
pop_size <- nrow(eur_haplo_mt)

## Create population stratification for microbiome data:
# Increase the counts of the 10 most common OTUs by 'fac' in AFR population
fac <- 1.4
afr_pi_fac <- rep(1, length(DM_parameter_est$pi))
afr_pi_fac[which.maxn(DM_parameter_est$pi, 10)] <- fac

#### Main simulation function ####
sim_func_KRV <- function(i) {
  if (n_pop == 2) {
    ## Form the genotypes of the n individuals (n/2 eur, n/2 afr)
    geno_matrix_list <- lapply(list(eur_haplo_mt,
                                    afr_haplo_mt), 
                               function(x) {pair_geno(n/2, pop_size, x)})
    
    # First n/2 individuals are EUR, second n/2 individuals are AFR.
    geno_matrix <- do.call(rbind, geno_matrix_list)
  }
  else if (n_pop == 1) {
    ## Form n euro individuals
    geno_matrix <- pair_geno(n, pop_size, eur_haplo_mt)
  }
  
  # construct a linear kernel based on the genotype matrix
  geno_matrix_af <- apply(geno_matrix, 2, mean)/2
  geno_matrix_maf <- sapply(geno_matrix_af, function(x){min(x, 1-x)})
  snp_id <- which(eur_afr_haplo_snp_pos <= 8000 &    # select common SNPs within a specific region
                    geno_matrix_maf >= 0.05) 
  geno.K <- geno_matrix[,snp_id] %*% t(geno_matrix[,snp_id])  
  
  
  # simulate microbiome count data and introduce population structure
  sim_otu_data <- simPop(J=n, K=856, n = 1000,   
                         pi = DM_parameter_est$pi, 
                         theta = DM_parameter_est$theta)
  sim_otu_tab <- sim_otu_data[[3]]
  
  afr_fac_mt <- matrix(rep(afr_pi_fac,n/2),n/2,856,byrow = T)
  
  sim_otu_tab[(n/2+1):n,] <- sim_otu_tab[(n/2+1):n,] * afr_fac_mt
  
  if (power_ind ==  1) {  # if performing power analysis, introduce genetic effect on otu counts
  
  # effect of one common snp on otu counts
  snp_idx <- sample(snp_id, 1)
  snp_fac1 <- rep(1, n) + 1.9*geno_matrix[,snp_idx]    # effect size (1.9) can be modified
  snp_fac2 <- rep(1, n) + 1.9*geno_matrix[,snp_idx]
  snp_shift <- rep(0, n) + 2*geno_matrix[,snp_idx]
  
  if (power_situ == 1) {  
  ## change the abundance of the 11th-20th most common OTUs
   max11.20 <- which.maxn(DM_parameter_est$pi, 20)[11:20]
   sim_otu_tab[,max11.20] <-sim_otu_tab[,max11.20] * snp_fac1 
  } else if (power_situ == 2) {
  ## change the abundance of OTUs from a common cluster 
   cluster_idx <- which(colnames(throat.otu.tab) %in% com_cluster)
   sim_otu_tab[,cluster_idx] <-sim_otu_tab[,cluster_idx] * snp_fac2 
  } else if (power_situ == 3) {
  ## change the abundance of 5 rare OTUs
    min1.40 <- which.minn(DM_parameter_est$pi, 40)
    min5 <- sample(min1.40, 5)
    sim_otu_tab[,min5] <-sim_otu_tab[,min5] + snp_shift
  }
  }
  
  # round the counts to integer
  sim_otu_tab <- round(sim_otu_tab)
  
  # rarefy otu matrix
  sim_otu_tab <- Rarefy(sim_otu_tab)$otu.tab.rff
  
  colnames(sim_otu_tab)= colnames(throat.otu.tab)
  
  # Construct distance matrices
  otu.D_bc <- as.matrix(vegdist(sim_otu_tab, method="bray")) 
  
  unifracs <- GUniFrac(sim_otu_tab, throat.tree, alpha=c(0, 0.5, 1))$unifracs
  otu.D_uu <- unifracs[,,"d_UW"]
  otu.D_uw <- unifracs[,,"d_1"]
  otu.D_ug <- unifracs[,,"d_0.5"]
  
  # convert distance matrices to kernel matrices
  otu.K_bc <- D2K(otu.D_bc)
  otu.K_uu <- D2K(otu.D_uu)
  otu.K_uw <- D2K(otu.D_uw)
  otu.K_ug <- D2K(otu.D_ug)
  
  otu.K.list <- list(otu.K_bc, otu.K_uu, otu.K_uw, otu.K_ug)
  
  ## Perform KRV with no adjustment
  p_value_noadj <- mapply(function(x){
    KRV(kernels.otu = x, kernel.y = geno.K)
  }, otu.K.list)
  
  ## Perform adjusted KRV
  # Perform PCA on the genotype matrix over the entire 1 Mb
  geno_matrix_centered <- scale(geno_matrix[,-snp_id], center = T,
                                scale = F)
  geno_pca <- propack.svd(geno_matrix_centered, neig = 1)
  
  # covariate matrix containing the first PC
  X <- cbind(1, geno_pca$u)
  
  p_value_adj <- mapply(function(x) {
                   KRV(kernels.otu = x, kernel.y = geno.K, X = X,
                       adjust.type = 'both') },
                             otu.K.list)
  
  
  if (power_ind ==  1) {  # perform analysis with competing methods
  ## Obtain 1st PC from microbiome kernel and 1st PC from genotype kernel
   # and then perform linear regression
    PC.m.list <- lapply(otu.K.list, function(x) {K_pc(x, 1)})
    PC.g <- K_pc(geno.K, 1)
    p_value_reg <- mapply(function(x){
      model1 <- lm(x ~ PC.g + geno_pca$u)
      return(summary(model1)$coefficients[2,4])
    },  PC.m.list)
  
  ## Obtain 1st PC from microbiome kernel and perform SKAT
    p_value_skat <- mapply(function(x){
      pheno.res <- SKAT_Null_Model(x ~ geno_pca$u)
      return(SKAT(geno_matrix[, snp_id], pheno.res,
                kernel = "linear")$p.value)
    }, PC.m.list)
  
  }
  
  if (power_ind == 0) {
    return(c(p_value_noadj, p_value_adj)) 
  } else {
    return(c(p_value_adj, p_value_reg, p_value_skat))
  }
}


#### Conduct simulation ####

if (power_ind == 0) {
   p_value_df <- matrix(NA, n_sim, 8) 
} else {
   p_value_df <- matrix(NA, n_sim, 12)
}

set.seed(seed_num)
for (i in 1:n_sim) {
  p_value_df[i,] <- sim_func_KRV(i)
}

if (power_ind == 0) {
  colnames(p_value_df) = c("BrayCurtis_noadj", "UnweightedUniFrac_noadj", 
                           "WeightedUniFrac_noadj", "GeneralizedUniFrac_noadj",
                           "BrayCurtis_adj", "UnweightedUniFrac_adj", 
                           "WeightedUniFrac_adj", "GeneralizedUniFrac_adj")
} else {
  colnames(p_value_df) = c("BrayCurtis_adjKRV", "UnweightedUniFrac_adjKRV", 
                             "WeightedUniFrac_adjKRV", "GeneralizedUniFrac_adjKRV",
                             "BrayCurtis_regression", "UnweightedUniFrac_regression", 
                             "WeightedUniFrac_regression", "GeneralizedUniFrac_regression",
                             "BrayCurtis_SKAT", "UnweightedUniFrac_SKAT", 
                             "WeightedUniFrac_SKAT", "GeneralizedUniFrac_SKAT"
                             )
}

write.csv(p_value_df, file=args$output_file, row.names = F)


#elapsed <- proc.time() - ptm
#print(elapsed)






