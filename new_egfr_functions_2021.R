## New Race Equations

## Adapted from equations on {nephro} package, updated for release of paper by Inker et al, NEJM 2021 (https://www.nejm.org/doi/full/10.1056/NEJMoa2102953)

NoRace_eGFRcr <- function(creatinine, age, sex) {
# Creatinine in micromol/l, age in years, sex == 0 is Female | sex == 1 is Male  
  if (!is.null(creatinine) & !is.null(sex) & !is.null(age))
  {
    creatinine <- as.numeric(creatinine) * 0.01131222
    sex <- as.numeric(sex)
    age <- as.numeric(age)
    
    n <- length(creatinine)
    
    if (length(sex) == n & length(age) == n)
    {
      # Identify missing data and store the index
      idx <- c(1:n)[is.na(creatinine) | is.na(sex) | is.na(age)]
      
      # Replace missing data with fake data to avoid problems with formulas
      creatinine[is.na(creatinine)] <- 10
      sex[is.na(sex)] <- 10
      age[is.na(age)] <- 10     
      
      # New eGFR equation for creatinine   
      k <- a <- numeric(n)
      k[sex==0] <- 0.7
      k[sex==1] <- 0.9
      a[sex==0] <- -0.241
      a[sex==1] <- -0.302
      one <- rep(1,n)
      eGFR <- apply(cbind(creatinine/k,one),1,min,na.rm=T)^a * apply(cbind(creatinine/k,one),1,max,na.rm=T)^-1.20 * 0.9938^age
      eGFR[sex==0] <- eGFR[sex==0] * 1.012
      
      # Restore missing data at the indexed positions
      eGFR[idx] <- NA
      
      # Output
      142 * eGFR
    } else
      stop ("Different number of observations between variables")
  } else
    stop ("Some variables are not defined")
}


NoRace_eGFRcrcys <- function(creatinine, cystatin, age, sex) {
  # Creatinine in micromol/l, age in years, sex == 0 is Female | sex == 1 is Male  
  if (!is.null(creatinine) & !is.null(cystatin) & !is.null(sex) & !is.null(age))
  {
    creatinine <- as.numeric(creatinine) * 0.01131222
    cystatin <- as.numeric(cystatin)
    sex <- as.numeric(sex)
    age <- as.numeric(age)
    n <- length(creatinine)
    
    if (length(sex) == n & length(age) == n & length(cystatin) == n)
    {
      # Identify missing data and store the index
      idx <- c(1:n)[is.na(creatinine) | is.na(cystatin) | is.na(sex) | is.na(age)]
      
      # Replace missing data with fake data to avoid problems with formulas
      creatinine[is.na(creatinine)] <- 10
      cystatin[is.na(cystatin)] <- 10
      sex[is.na(sex)] <- 10
      age[is.na(age)] <- 10     
      
      # New eGFR equation for creatinine and cystatin combination  
      k_sex <- rep(1,n)
      k_sex[sex==0] <- 0.969
      k <- rep(0.9,n)
      k[sex==0] <- 0.7
      a <- rep(-0.144,n)
      a[sex==0] <- -0.219
      one <- rep(1,n)
      CR <- cbind(creatinine/k,one)
      CY <- cbind(cystatin/0.8,one)                 
      eGFR <- 135 * apply(CR,1,min,na.rm=T)^a * apply(CR,1,max,na.rm=T)^-0.544 * apply(CY,1,min,na.rm=T)^-0.323 * apply(CY,1,max,na.rm=T)^(-0.778) * 0.9961^age * k_sex
      # Restore missing data at the indexed positions
      eGFR[idx] <- NA
      
      # Output
      eGFR
    } else
      stop ("Different number of observations between variables")
  } else
    stop ("Some variables are not defined")
}