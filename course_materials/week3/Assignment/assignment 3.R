# week 3 assignment

# 1. Try creating risk ratios with these new case-control data (simulated)
nairobi_cases <- read.csv("cases_nairobi.csv")

# 2. Identify 'best' IDW model using CV-MSE values for different 
# powers of IDW function using the BF malaria data
# HINT: you can do this manually, or you could wrap this into a function..


# 3. Open this dataset of hookworm in Uganda, 
# and compare the best IDW surface to a kriged prevalence across 
# the window. HINT: when using the logit transform, 
# you may need to add a small amount
# (i.e 0.001) to any values of 0 prevalence
HK<-read.csv("tanzania_uganda_hkprev.csv")
