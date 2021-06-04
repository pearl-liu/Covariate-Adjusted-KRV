#### Functions used in simulations and demo ####
library(ggplot2)
library(reshape2)

# Randomly pair haplotypes from the 10,000 haplotypes to form the 
# genotypes of the n individuals. 
pair_geno <- function(n, pop_size, haplo_mt) {
  haplo_idx <- sample(pop_size, 2*n)
  haplo_matrix <- haplo_mt[haplo_idx, ]
  haplo_mt_group <- rep(1:n, each=2)
  geno_matrix <- rowsum(haplo_matrix, 
                        haplo_mt_group)
  return(geno_matrix)
}


## obtain i-th PC from kernel PCA
K_pc <- function(K, i) {
  eK <- eigen(K, symmetric=TRUE)
  pc <- eK$vector[, which.maxn(eK$values, i)[i]]
  return(pc)
}

## Output empirical type I error rates or powers based on p-value results
result_eval <- function(power_ind, file_name, sig_level) {
  p_val_list <-  read.csv(file_name)
  if (power_ind == 0) { # evaluate type I error
    emp_type_I_error <- data.frame("Kernel"=c('Bray-Curtis', 'Unweighted UniFrac', 
                                              'Weighted UniFrac', 'Generalized UniFrac'),
                                   "alpha"=sig_level,
                                   "UnadjustedKRV"=NA, "AdjustedKRV"=NA)
    for (i in 1:4) {
      emp_type_I_error$UnadjustedKRV[i] <-  sum(p_val_list[,i] < sig_level)/nrow(p_val_list)
      emp_type_I_error$AdjustedKRV[i] <-  sum(p_val_list[,i+4] < sig_level)/nrow(p_val_list)
    }
    return(emp_type_I_error)
  } else {  # evaluate power
    emp_power <- data.frame('Method'=c('AdjustedKRV','LinearRegression','SKAT'),
                                      'BrayCurtis'=NA, 'UnweightedUniFrac'=NA, 
                                       'WeightedUniFrac'=NA, 'GeneralizedUniFrac'=NA)
    for (i in 1:4) {
      emp_power[1,i+1] <-  sum(p_val_list[,i] < sig_level)/nrow(p_val_list)
      emp_power[2,i+1] <-  sum(p_val_list[,i+4] < sig_level)/nrow(p_val_list)
      emp_power[3,i+1] <-  sum(p_val_list[,i+8] < sig_level)/nrow(p_val_list)
    }
    return(emp_power)
  }
}

## Make a plot comparing power between different methods
addline_format <- function(x,...){
  gsub('\\s','\n',x)
}

power_plot <- function(power_result, power_scenario,
                       kernel_name_output) {
  power_result_m <- melt(power_result, id.vars = 'Method')
  colnames(power_result_m)[2:3] <- c('Kernel', 'Power') 
  power_result_m$Kernel <- factor(power_result_m$Kernel, 
                                  levels = c("BrayCurtis", "UnweightedUniFrac", 
                                             "WeightedUniFrac", "GeneralizedUniFrac"),
                                  labels =  kernel_name_output)
  output_plot <- ggplot(power_result_m, aes(Kernel, Power, fill = Method)) +
    geom_bar(stat = 'identity', position = 'dodge') +
    theme_bw()+
    theme(
      plot.title = element_text(hjust = 0.5),
      #  legend.title=element_text(size=20),
      #  legend.text=element_text(size=18),
      legend.position = "bottom") +
    #  legend.text.align = 0,
    #  axis.title=element_text(size=20),
    #  axis.text=element_text(size=17)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0,1)) +
    scale_fill_manual(values = c('blue','orange','green'),
                      labels = c('Adjusted KRV', 'Linear Regression', 'SKAT'))+
    labs(title = paste0('Power Scenario ', power_scenario)) +
    scale_x_discrete(labels = addline_format(kernel_name_output))
  
  return(output_plot)
}


