# Run this part:
library(devtools)

set.seed(100)

nProd <- 10
nFirms <- 100
# Max: 26

prodNames <- letters[1:nProd]

outputPriceMat <- matrix(runif(nProd * nFirms), ncol = nProd)
colnames(outputPriceMat) <- paste0("P", prodNames)
outputQuantMat <- matrix(runif(nProd * nFirms), ncol = nProd)
colnames(outputQuantMat) <- paste0("Q", prodNames)

outputPriceMat <- round(outputPriceMat + 1, 5)
outputQuantMat <- round(outputQuantMat + 1, 5)

# write.csv(as.data.frame(cbind(outputPriceMat, outputQuantMat) ),
#  #file = "~/tests/testthat/priceQuantMat.txt",
#  file = "~/svn/micEcon/pkg/micEconIndex/tests/testthat/priceQuantMat.txt",
#  row.names = FALSE)

priceQuantMat <- as.data.frame(cbind(outputPriceMat, outputQuantMat) )
# devtools::use_data(, internal = TRUE, overwrite = TRUE)

# Then pause and run TFPIP/calcTFPIP.R according to instructions there
# And then run the below:

load("data-raw/TFPIPcheck.Rdata", verbose = TRUE)
devtools::use_data(TFPIPresult, priceQuantMat, internal = TRUE, overwrite = TRUE)


#quantityIndex( paste0("P", prodNames),
#  paste0("Q", prodNames), 1,
# as.data.frame(cbind(outputPriceMat, outputQuantMat) ))






#eg1-dta.txt        DATA FILE NAME
#eg1-out.txt        OUTPUT FILE NAME
#5              NUMBER OF OBSERVATIONS
#2              NUMBER OF OUTPUTS
#3              NUMBER OF INPUTS
#0              0=TORNQVIST AND 1=FISHER
#0	         0=NON-TRANSITIVE AND 1=TRANSITIVE
