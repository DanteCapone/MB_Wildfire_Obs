librarian::shelf(dataRetrieval, ggplot2, tidyverse)
# Choptank River near Greensboro, MD
siteNumber <- "11161000"
ChoptankInfo <- readNWISsite(siteNumber)
parameterCd <- "00060"

# Raw daily data:
SLR_outflow <- readNWISdv(
  siteNumber, parameterCd,
  "1980-01-01", "2024-04-01"
)

pCode <- readNWISpCode(parameterCd)
names(pCode)
unique(pCode$parameter_units)


#rename columns
SLR_outflow = renameNWISColumns(SLR_outflow)


## Visualize

SLR_outflow %>% 
  ggplot(aes(Date,Flow))+
  geom_point()
here()
write.csv(SLR_outflow,here("san_lorenzo_daily_outflow.csv"))

# -------------------------------------------------------------------------

specificCond <- readWQPqw(
  siteNumbers = "WIDNR_WQX-10032762",
  parameterCd = "Specific conductance",
  startDate = "2011-05-01",
  endDate = "2011-09-30"
)
