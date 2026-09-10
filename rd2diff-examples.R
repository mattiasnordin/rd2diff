# load the file containing the functions 
# (have to be in the correct directory)
source("rd2diff.R")

G <- 20
ng <- 50 + rbinom(G, 400, .5)
n <- sum(ng)

# Create basic dataset with a running variable with equal spacing
# (that is, ordinal running variable)
df1 <- data.frame(x=rep(1:G, ng))
cut1 <- 12
df1$d <- ifelse(df1$x >= cut1, 1, 0)
df1$y <- df1$x + .5 * df1$d + rnorm(n)

# Two different ways of running rd2diff
rd2diff(formula=y ~ x, cutoff=cut1, data=df1)
rd2diff(df1$y, df1$x, cutoff=cut1)

# Restrict the bandwidth such that xc-3 <= xg <= xc+2
rd2diff(formula=y ~ x, cutoff=cut1, data=df1, h=c(3, 2))

# Create basic dataset with a running variable with unequal spacing
df2 <- data.frame(x=rep(rnorm(G), ng))
cut2 <- 0
df2$d <- ifelse(df2$x >= cut2, 1, 0)
df2$y <- df2$x + .5 * df2$d + .1 * rnorm(n)

rd2diff(formula=y ~ x, cutoff=cut2, data=df2)
rd2diff(df2$y, df2$x, cutoff=cut2)

# Create dataset with ordinal running variable and covariate
df3 <- data.frame(x=rep(1:G, ng))
cut3 <- 12
df3$d <- ifelse(df3$x >= cut3, 1, 0)
df3$z <- rnorm(n) + df3$d
df3$y <- df3$x+ df3$z + .5 * rnorm(n)
rd2diff(formula=y ~ x, cutoff=cut3, data=df3)
rd2diff(formula=y ~ x | z, cutoff=cut3, data=df3)

# Aggregate data for first example
dfagg <-  do.call(data.frame, aggregate(
  . ~ x, data=df1[, c("x", "y")],
  FUN = function(x) c(mean = mean(x), var = var(x), n=length(x))))

rd2diff.aggregate(dfagg$y.mean, dfagg$x, cutoff=cut1,
                  std_dev_depvar_means=sqrt(dfagg$y.var),
                  ng=dfagg$y.n)
