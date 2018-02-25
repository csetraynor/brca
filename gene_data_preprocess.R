library(tidyverse)
#obtain gene expression data
de <- read_tsv("brca_metabric/data_expression.txt", col_names = TRUE)[-2]

#do you need to sample genes? 
d <- de[sample(nrow(de), 100), ]
  
# first remember the names
n <- d$Hugo_Symbol
# transpose all but the first column (name)
df.aree <- as.data.frame(t(d[,-1]))
colnames(df.aree) <- n
df.aree$myfactor <- factor(row.names(df.aree))
str(df.aree) # Check the column types
#Clinical data
md <- read_tsv("brca_metabric/data_clinical_patient.txt", skip = 4)
colnames(md) <- tolower(colnames(md))
Y <- md[md$patient_id %in% df.aree$myfactor,] 

Y %>%
  VIM::aggr(prop = FALSE, combined = TRUE, numbers = TRUE, sortVars = TRUE, sortCombs = TRUE)
# 
# y <- Y %>% select(patient_id, os_status, os_months)
require(impute)
X <- tbl_df(impute.knn(as.matrix(df.aree[-101]), k = 10)$data)
X <- X %>% cbind(patient_id = df.aree$myfactor) 
pre_data <- left_join(Y, X)[-1]


# Get metastasis data for the case list
md <- read_tsv("msk_impact_2017_clinical_data.tsv")
#select only breast cancer
mdg <- md %>% filter(`primary tumor site` == "Breast") %>% 
  rename(sample_id = `sample id`,
         patient_id = `patient id`) 
mdg <- mdg[mdg$patient_id %in% Y$patient_id, ]

#--- This function will take a formula object as input --- #
gen_stan_data <- function(data, formula = as.formula(~1)) {
  if(!inherits(formula, 'formula'))
    formula <- as.formula(formula)
  
  observed_data <- data %>%
    dplyr::filter(os_status == "DECEASED")
  
  censored_data <- data %>%
    dplyr::filter(os_status != "DECEASED")
  
  Xobs_biom <- observed_data %>%
    model.matrix(formula, data = .)
  
  Xcen_biom <- censored_data %>%
    model.matrix(formula, data = .)
  
  assertthat::assert_that(ncol(Xcen_biom) == ncol(Xobs_biom))
  M_biom <- ncol(Xcen_biom)
  
  if (M_biom > 1){
    if("(Intercept)" %in% colnames(Xobs_biom))
      Xobs_biom <- array(Xobs_biom[,-1], dim = c(nrow(observed_data), M_biom - 1))
    if("(Intercept)" %in% colnames(Xcen_biom))
      Xcen_biom <- array(Xcen_biom[,-1], dim = c(nrow(censored_data), M_biom - 1))
    assertthat::assert_that(ncol(Xcen_biom) == ncol(Xobs_biom))
    M_biom <- ncol(Xcen_biom)
  }
  
  stan_data <- list(
    Nobs = nrow(observed_data),
    Ncen = nrow(censored_data),
    yobs = observed_data$os_months,
    ycen = censored_data$os_months,
    M_biom = M_biom,
    Xcen_biom = array(Xcen_biom, dim = c(nrow(censored_data), M_biom)),
    Xobs_biom = array(Xobs_biom, dim = c(nrow(observed_data), M_biom))
  )
}

gen_inits <- function(M_biom){
  function()
    list(
      alpha_raw = 0.01*rnorm(1),
      mu = rnorm(1),
      tau_s1_biom_raw = 0.1*abs(rnorm(1)),
      tau_s2_biom_raw = 0.1*abs(rnorm(1)),
      tau2_biom_raw = array(abs(rnorm(M_biom)), dim = c(M_biom)),
      tau1_biom_raw = array(abs(rnorm(M_biom)), dim = c(M_biom)),
      beta_biom_raw = array(rnorm(M_biom), dim = c(M_biom))
    )
}
###------ Run Stan --------##
x <- colnames(pre_data)[3:102]
nChains <- 1
stanfile <- "pem/wei_hs.stan"
testfit <- rstan::stan(stanfile,
                       data = gen_stan_data(pre_data), 
                as.formula(paste("~", paste(x, collapse = "+"))),
                       init = gen_inits,
                       iter = 4,
                       chains = 1)
