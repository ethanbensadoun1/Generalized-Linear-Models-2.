Quasi-likelihood Methods and Overdispersion for Ray Allen 3-Pointer Data
Set
================
Ethan Bensadoun
2024-08-15

``` r
Ray3s <- read.csv("https://users.stat.ufl.edu/~aa/glm/data/Basketball.dat", 
                  header= FALSE, sep = "",stringsAsFactors = FALSE)
Ray3s <- Ray3s[-1,]
names(Ray3s) <- c("Game", "Made", "Attempts")
rownames(Ray3s) <- NULL
Ray3s
```

    ##    Game Made Attempts
    ## 1     1    0        4
    ## 2     2    7        9
    ## 3     3    4       11
    ## 4     4    3        6
    ## 5     5    5        6
    ## 6     6    2        7
    ## 7     7    3        7
    ## 8     8    0        1
    ## 9     9    1        8
    ## 10   10    6        9
    ## 11   11    0        5
    ## 12   12    2        5
    ## 13   13    0        5
    ## 14   14    2        4
    ## 15   15    5        7
    ## 16   16    1        3
    ## 17   17    3        7
    ## 18   18    0        2
    ## 19   19    8       11
    ## 20   20    0        8
    ## 21   21    0        4
    ## 22   22    0        4
    ## 23   23    2        5
    ## 24   24    2        7

``` r
# Here I needed to upload the three-point shooting data of Ray Allen, but making
# sure that each column is seperated into the 3 variables
```

Now that the data is successfully uploaded, we can take a look at the
problem:

- Let pi_i = Pr(Y_ij = 1) = 1 - Pr(y_ij = 0) and y_ij = sum_j (y_ij /
  n_i) be the sample proportion of successes.
- In this case, pi_i represents the probability of successful
  three-point shots by Ray in game i.
- For the purpose of this problem, we can make the necessary assumption
  of independence as indicated in the problem, {yi} is independent.
- From the problem we also know the following:
  - (n_i)(y_i) = number three points shots made out of ni
  - Attempts is a bin(n_i,p_i) variate
  - yi is the observed proportion of a successful three-point shot

## (a)

Here we can fit the model where pi = Beta0,

``` r
# I tried running the GLM model, but it didn't work since the variables in the 
# data frame were not numeric, thus, we need to convert them:
Ray3s$Made <- as.numeric(as.character(Ray3s$Made))
Ray3s$Attempts <- as.numeric(as.character(Ray3s$Attempts))
# We can follow the form success/attempts, where the predictor is an intercept
Ray3s_fit <- glm(cbind(Made, Attempts - Made) ~ 1, 
                 family = binomial, data = Ray3s)

Ray3s_fit_summary <- summary(Ray3s_fit)
Ray3s_fit_summary
```

    ## 
    ## Call:
    ## glm(formula = cbind(Made, Attempts - Made) ~ 1, family = binomial, 
    ##     data = Ray3s)
    ## 
    ## Coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)   
    ## (Intercept)  -0.4633     0.1706  -2.716   0.0066 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 58.307  on 23  degrees of freedom
    ## Residual deviance: 58.307  on 23  degrees of freedom
    ## AIC: 96.384
    ## 
    ## Number of Fisher Scoring iterations: 4

``` r
Beta0_hat <- coef(Ray3s_fit)["(Intercept)"]
Beta0_hat
```

    ## (Intercept) 
    ##  -0.4632847

``` r
standard_error <- Ray3s_fit_summary$coefficients["(Intercept)", "Std. Error"]
standard_error
```

    ## [1] 0.170567

From the summary table above, we find that the estimated intercept,
Beta0_hat = -0.4633 and the standard error, SE = 0.1706.

- The fact that Beta0_hat is negative implies that the log-odds of a
  successful three-point shot are less than zero.

## (b)

A factor that can cause overdisperion in the case of Ray Allen’s
three-point stats, is that some of the shots made may not actually be
“independent” in nature. Take for example a high-scoring game, where a
player may be having a hot streak or the inveserse if he is not
performing very well due to health conditions. This may cause
overdisperion, as the three-point shots made given the number of
attempts are actually dependent on each other.

``` r
# Here we can calculate the Pearson chi-square statistic which assesses the 
# goodness-of-fit of the Ray3s_fit model above
X2 <- sum(residuals(Ray3s_fit,type="pearson")^2)
X2
```

    ## [1] 46.82661

``` r
n <- 24 # this is the total number of games played by Ray in the dataset 
p <- 1 # this is the number of parameters (in this case just Beta0 hat)
df_residual <- n - p
df_residual
```

    ## [1] 23

``` r
# Now we can calculate phi (i.e., the overdispersion constant)
phi <- X2/df_residual
phi
```

    ## [1] 2.035939

Notice that phi \> 1 represents overdispersion of the model in this
case, since 2.035939 \> 1.

Now for the Adjusted standard error:

``` r
Adjusted_SE <- sqrt(phi) * standard_error
Adjusted_SE
```

    ## [1] 0.2433758

As a result, we can say that 0.2433758, is a much more plausible
standard error for Beta0 hat in this prediction equation.

Now we can find and compare the 95% C.Is for Beta0:

``` r
# confidence interval for Beta0 hat , before the model's SE is adjusted:
C.I_original <- confint(Ray3s_fit)
```

    ## Waiting for profiling to be done...

``` r
C.I_original
```

    ##      2.5 %     97.5 % 
    ## -0.8026429 -0.1324412

``` r
# now we need to set the lower and upper bounds of the C.I given the adjusted SE:
lwr <- Beta0_hat - qnorm(0.975) * Adjusted_SE
upr <- Beta0_hat + qnorm(0.975) * Adjusted_SE
C.I_adjusted <- c(lwr, upr)
C.I_adjusted
```

    ## (Intercept) (Intercept) 
    ## -0.94029248  0.01372313

C.I_original represents the 95% confidence interval of the estimated
intercept in the model, under the assumption of no overdispersion. Note
that, as calculated in part (a) of the problem, Beta0_hat = -0.4633,
which falls in the C.I. (-0.8026429 -0.1324412). Thus, we can say that
the true average log-odds of making a Ray Allen shooting a three-point
shot, seems to be true 95% of the time within the range that was
calculated.

As expected, the adjusted confidence interval is wider than the
original. This is due to the extra-binomial variation from phi. Since
the interval is larger, (-0.94029248 0.01372313 ), we can say that there
is increased uncertainty about what Beta0_hat, which is more realisitic
when prediciting the accuracy of the intercept. In addition, the fact
that the adjusted C.I is wider can mean that the data exhibit more
variability than what the binomial Ray3s_fit model reports.
