library(dplyr)
library(tidyverse)
library(data.table)
library(readxl)
library("writexl")

setwd("./DATASET")

f <- list.files(pattern="tsv$")
TOTALS <- do.call(bind_rows, lapply(f, function(file) read.csv(file,header = TRUE, sep = "\t")))

TOTALS <- as.data.frame(TOTALS)


write_xlsx(TOTALS,"./TOTAL_String.xlsx")
write.table(TOTALS, file="./TOTAL_String_tab.csv",sep = "\t")








#TOTALS <- do.call(bind_rows, lapply(f, function(file) read.table(file,header = TRUE, sep = "\t",col.names = )))



#TOTALSread_table <- do.call(bind_rows, lapply(f, function(file) read_table(file,col_names = TRUE)))


