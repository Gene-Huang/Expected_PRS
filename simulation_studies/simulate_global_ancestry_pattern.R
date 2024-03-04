
# generate an n by 3 matrix representing a 3-way admixed population
# each row provides the proportion global ancestries of an individual. 
# multiple admixture patterns are implemented. 
simulate_global_ancestry_pattern <- function(n, 
                                             pattern = "realistic_admixture"){
  stopifnot(is.element(pattern, c("realistic_admixture", 
                                  "fixed_admixture", 
                                  "a1", 
                                  "a2", 
                                  "mixed_a1_a2")))
  
  if (pattern == "realistic_admixture"){
    a1_prop <- runif(n, 0,1)
    a2_prop <- runif(n,  rep(0, n),1 - a1_prop)
    a3_prop <- rep(1, n) - a1_prop - a2_prop
    admixture_prop <- data.frame(a1 = a1_prop, 
                                 a2 = a2_prop, 
                                 a3 = a3_prop)
    
  }
  
  if (pattern == "fixed_admixture"){
    admixture_prop <- data.frame(a1 = rep(0.4, n), 
                                 a2 = rep(0.3, n),
                                 a3 = rep(0.3, n))
  }
  
  if (pattern == "a1"){
    admixture_prop <- data.frame(a1= rep(1, n), 
                                 a2 =0, 
                                 a3 = 0)
  }
  
  if (pattern == "a2"){
    admixture_prop <- data.frame(a1= rep(0, n), 
                                 a2 =1, 
                                 a3 = 0)
  }
  
  if (pattern == "mixed_a1_a2"){
    inds_a1 <- sample(1:n, size = 0.6*n)
    a1 <- rep(0, n)
    a1[inds_a1] <- 1
    a2 <- rep(1,n)
    a2[inds_a1] <- 0
    admixture_prop <- data.frame(a1= a1, 
                                 a2 =a2, 
                                 a3 = 0)   
  }  
  
  rownames(admixture_prop) <- paste0("person_", 1:n)
  return(admixture_prop)
}





