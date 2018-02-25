library(tidyverse)
library(caret)
library(cgdsr)

#download data
# Create CGDS object
mycgds = CGDS("http://www.cbioportal.org/public-portal/")
test(mycgds)
# Get list of cancer studies at server
list_cancer = getCancerStudies(mycgds)
# Get available case lists (collection of samples) for a given cancer study
mycancerstudy = getCancerStudies(mycgds)[23,1]
mycaselist = getCaseLists(mycgds,mycancerstudy)[2,1]
# Get medical data for the case list
# md <- read_tsv("brca_tcga_pub2015/data_clinical_patient.txt", col_names = TRUE, skip = 4)
#Definitive
md <- tbl_df(getClinicalData(mycgds, mycaselist))

# # Get genetic profile
# mygeneticprofile = getGeneticProfiles(mycgds,mycancerstudy)[2,1]

# Only 423 
ge <- read_tsv("brca_tcga_pub2015/data_expression_median.txt", col_names = TRUE)

# Get sample names in rownames
md <- md %>% tibble::rownames_to_column("num_id") 
ge <- tibble::rownames_to_column(ge, var = "sample_id")

# Get an easier code to read
names(md) <- tolower(names(md))
names(ge) <- tolower(names(ge))

#---- Data Cleaning ----#
#convert missig values into NA
convert_blank_to_na <- function(x) {
  if(!purrr::is_character(x)){
    warning('input variate not character - return original')
    return(x)
  } else {
    ifelse(x == '', NA, x)
  }
}
md <- md %>% dplyr::mutate_all(funs(convert_blank_to_na))

VIM::marginplot(md[c("metastatic_site","os_months")])

md %>%
  VIM::aggr(prop = FALSE, combined = TRUE, numbers = TRUE, sortVars = TRUE, sortCombs = TRUE)



# Only 423 
dme <- read_tsv("brca_igr_2015/data_mutations_extended.txt", col_names = TRUE, skip = 1)















Xdummies <- caret::dummyVars( ~ .,data =  ge[,2:4]) #for computational reason we pick two genes.
X <- predict(Xdummies, newdata =  ge[,2:4])
ge <- X %>%
  select(-contains("NaN")) %>%
  cbind(ge$`sample id`)


ge[]
msk <- left_join(md, ge, by = "Sample ID")


glio <- msk %>% filter(`Cancer Type` == "Glioma")
#obtain id
ge <- read_tsv("msk_impact_2017/data_mutations_uniprot.txt")
id_ge <- ge$Tumor_Sample_Barcode
#obtain sample id
Y <- md[md$`Sample ID` %in% id_ge,c(1:6)] 







