library(data.table)
library(tidyverse)
flow<- fread("C:/Users/vvijayak/Documents/GitHub/supply-chain-inputs/cases/BASECASE_Centralonly-lowercap_stationonly_e5/outputs/flow_SERA_2.tsv")
listofstat<-fread("C:/Users/vvijayak/Downloads/listofstations.csv")
listofstat$id<-1
comb<-left_join(flow,listofstat,by="Infrastucture_ID")
comb2<-comb%>% mutate_at(11, ~replace_na(.,0))
comb_filter<-comb2%>%filter(id!=1)
# write flow output file
write.table(comb_filter,file="C:/Users/vvijayak/Documents/GitHub/supply-chain-inputs/cases/BASECASE_Centralonly-lowercap_stationonly_e5/outputs/flow_SERA_2_.tsv", sep = "\t",row.names = F, quote = F)
