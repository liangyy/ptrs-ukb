# ARGS1: input query header including sample filters and covariates
# ARGS2: output query yaml

args = commandArgs(trailingOnly = TRUE)

library(dplyr)
source('../../code/rlib_doc.R')

# read in data dictionary for uk biobank
dict = read.csv("~/Downloads/Data_Dictionary_Showcase.csv")

# read in Martin et al trait list
qtrait_sub = read.csv('t2d_traits.csv')
# qtrait_sub$UKBB.code = as.character(qtrait_sub$UKBB.code)
ethnicityid = 21000

# extract the target trait from dictory and construct query key:value pairs
myextract = dict %>% filter(FieldID %in% c(qtrait_sub$UKBB.code, ethnicityid))
myextract = myextract %>% left_join(
  rbind(
    qtrait_sub %>% select(UKBB.code, Trait),
    data.frame(UKBB.code = ethnicityid, Trait = 'ethnicity')
  ), by = c('FieldID' = 'UKBB.code'))
mydata = meta_info_to_query(myextract$Trait, myextract$FieldID, myextract$Instances, myextract$Array)


# check whether the values in query exists in current ukbREST instance
## load fields existing in ukbREST instance 
good_fields = readRDS('../../output/fields_on_nucleus.rds')
## and check/filter non-existing ones
for(n in names(mydata)) {
  if(! mydata[[n]] %in% good_fields) {
    mydata[[n]] = NULL
  }
}


# load in query with sample filters and covars
# and add these new fields in
current_yaml = yaml::read_yaml(args[1])  
for(n in names(mydata)) {
  current_yaml$data[[n]] = mydata[[n]]
}
yaml::write_yaml(current_yaml, args[2])
