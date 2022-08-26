rm(list=ls())
library(data.table)
setwd("C:/Users/chunter/Documents/Projects/SERA/Infrastructure-Optimization/Supply-Dynamics/supply-chain/inputs/prices")

usa <- fread(paste0("energy-prices-reference-all-eia-real-dollars.tsv"), header = TRUE, sep = "\t", stringsAsFactors = TRUE)

canada <- usa[, .(`Price [$/unit]` = mean(`Price [$/unit]`), "Zone" = "Canada", "Billable?" = "True"), by = c("Material","Year")]
canada[Material == "Electricity (Industrial) [kWh]" | Material == "Electricity (Commercial) [kWh]", `Price [$/unit]` := `Price [$/unit]` * (7.31/8.71)]

write.table(canada, file = "canada-energy-prices-approximation.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
