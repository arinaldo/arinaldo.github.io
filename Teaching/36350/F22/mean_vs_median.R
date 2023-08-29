# For consistency, use the same seed as in lab 1
set.seed(6433635)
huge.bin.draws <- rbinom(200e6, 10e6, 50e-8)

# crude benchmark of mean vs median
system.time(mean(huge.bin.draws))
system.time(median(huge.bin.draws))

# this is much faster:a
system.time(sum(huge.bin.draws)/length(huge.bin.draws))




# better benchmarking - same results
install.packages(c("bench"))
library("bench")
set.seed(6433635)
huge.bin.draws <- rbinom(200e6, 10e6, 50e-8)

bnch_mean <- bench::mark(
    mean(huge.bin.draws),
    sum(huge.bin.draws)/length(huge.bin.draws)
    )
# Display results
bnch_mean

bnch_median <- bench::mark(
+ median(huge.bin.draws) )
bnch_median
