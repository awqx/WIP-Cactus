# Models evaluated by standards from Tropsha 2010, DOI: 10.1002/minf201000061
# "Best Practices for QSAR Model Development, Validation, and Exploitation"


# Set working directory as needed
# setwd("~/SREP LAB/qsar")

source("./model.code/eval.functions.R")

# Importing Models' Results -----------------------------------------------

cube <- readRDS("./models/cubist/compiled.results.RDS")
svm  <- readRDS("./models/svm/polysvm.tst.results.RDS")
rf   <- readRDS("./models/rforest/rf.results.RDS") 
glm  <- readRDS("./models/glmnet/glm.results.RDS")
pls  <- readRDS("./models/pls/pls.results.RDS")

# Evaluating Standards ----------------------------------------------------

# 0.85 <= k <= 1.15
find.k(cube) # 0.985
find.k(svm) # 0.954
find.k(rf) # 0.965
find.k(glm) # 1.002
find.k(pls) # 0.951

# 0.85 <= k-prime <- 1.15
find.kprime(cube) # 0.972
find.kprime(svm) # 0.996
find.kprime(rf) # 0.982
find.kprime(glm) # 0.956
find.kprime(pls) # 1.009

# Q2-F1 > 0.5
find.q2.f1(cube) # 0.659
find.q2.f1(svm) # 0.601
find.q2.f1(rf) #0.678
find.q2.f1(glm) # 0.694
find.q2.f1(pls) # 0.683

# (R2-R20-prime)/R2 < 0.1
proc.time()
(find.r2(cube) - find.r20prime(cube))/find.r2(cube) # 0.040
(find.r2(svm) - find.r20prime(svm))/find.r2(svm)  # 0.020
(find.r2(rf) - find.r20prime(rf))/find.r2(rf)  # 0.032
(find.r2(glm) - find.r20prime(glm))/find.r2(glm) # 0.009
(find.r2(rf) - find.r20prime(pls))/find.r2(pls)  # 0.024

# abs(R20- - R20prime) < 0.3
abs(find.r20(cube) - find.r20prime(cube)) #0.007
abs(find.r20(svm) - find.r20prime(svm)) # 0.080
abs(find.r20(rf) - find.r20prime(rf)) #0.312
abs(find.r20(glm) - find.r20prime(glm)) # 0.187
abs(find.r20(pls) - find.r20prime(pls)) #0.113


# Q2 ----------------------------------------------------------------------

# On one hand, Xu 2015 provides a single equation for evaluating Q^2
# On the other hand, it seems to be more standard to actually create and 
# evaluate models. 
# 10 fold CV used on the training set only
# Data changes betweem different CV methods: sparse, df, or mat

View(cube)
