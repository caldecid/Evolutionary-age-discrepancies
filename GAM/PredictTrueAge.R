# install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(cmdstanr)
# install_cmdstan(overwrite = TRUE, cores = 1) # Only once to configure cmdstanr
library(brms)
library(mgcv)


# Read data
A <- read.table('ages_training_predictors.csv', header = TRUE, sep = ',')


# Log transform ages 
A$true_age_log <- log10(A$true.age)
A$estimated_age_log <- log10(A$Estimated.age)
A$mode <- as.factor(A$mode)


# Divide data set in test and training part
set.seed(42)
Train <- sample(1:nrow(A), nrow(A) * 0.9)
Test <- c(1:nrow(A))[-Train]


# Maximum likelihood GAM
########################
Gam <- bam(true_age_log ~ mode + 
             sister +
             s(estimated_age_log, by = mode) +
             s(div, by = mode) +
             s(turnover, by = mode) +
             s(root.age, by = mode) +
             s(number.sp, by = mode),
           data = A[Train, ],
           samfrac = 0.1)
summary(Gam)

# True vs predicted age
PredGam <- predict(Gam, newdata = A[Test, ])
plot(A$true.age[Test], 10**PredGam, col = A$mode[Test])



# Bayesian GAM (lasts ~2 hours)
###############################
GamBrm <- brm(true_age_log ~ mode + 
                s(estimated_age_log, by = mode) +
                s(div, by = mode) +
                s(turnover, by = mode) +
                s(root.age, by = mode) +
                s(number.sp, by = mode),
              data = A[Train, ], 
              family = 'gaussian',
              chains = 1, cores = 1,
              backend = 'cmdstanr',   # Not sure whether this works on Windows
              threads = threading(2), # Multithreading, will this works on Win?
              control = list(adapt_delta = 0.8, max_treedepth = 10),
              iter = 2000, warmup = 1000, thin = 1)

summary(GamBrm)
bayes_R2(GamBrm)

# Plot effect of predictors
# (Maybe you get some nicer axis labels?)
Eff <- conditional_effects(GamBrm)
plot(Eff)


# Save and load Bayesian GAM to safe time
saveRDS(GamBrm, 'GamBrm_Log10transf.rds')
GamBrm <- readRDS('GamBrm_Log10transf.rds')


# Corrected age based on BGAM
PredBGAM <- predict(GamBrm, newdata = A)


# Means square error before correction
raw_mse <- mean((A[Test, 'true_age_log'] - A[Test, 'estimated_age_log'])**2)
# Means square error after correction
test_mse <- mean((A[Test, 'true_age_log'] - PredBGAM[, 1])**2)

TestBud <- A$mode[Test] == 0
# Means square error before correction for budding speciation
mean((A[Test, 'true_age_log'][TestBud] - A[Test, 'estimated_age_log'][TestBud])**2)
# Means square error after correction for budding speciation
mean((A[Test, 'true_age_log'][TestBud] - PredBGAM[, 1][TestBud])**2)
# Means square error before correction for cladogenetic speciation
mean((A[Test, 'true_age_log'][!TestBud] - A[Test, 'estimated_age_log'][!TestBud])**2)
# Means square error after correction for cladogenetic speciation
mean((A[Test, 'true_age_log'][!TestBud] - PredBGAM[, 1][!TestBud])**2)


# Coverage
InsideCredInt <- sapply(1:nrow(A), function(x) 
  PredBGAM[x, 3] < A$true_age_log[x] && PredBGAM[x, 4] > A$true_age_log[x] )
sum(InsideCredInt) / length(InsideCredInt)
sum(InsideCredInt[Test]) / length(Test) # Coverage test data

# Plot true vs predicted age
Col <- rep('deepskyblue', nrow(A))
Col[Test] <- 'orange2'
par(las = 1, mar = c(4, 4, 0.5, 0.5))
plot(A$true_age_log, PredBGAM[, 1], 
     pch = 19, cex = 0.5, col = Col,
     xlab = 'True age', ylab = 'Predicted age')
for (i in 1:nrow(A)) {
  lines(rep(A$true_age_log[i], 2), PredBGAM[i, 3:4],
        col = Col[i], lwd = 0.3)
}
abline(a = 0, b = 1, lty = 2)

