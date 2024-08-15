Heterogeneity of Poissoin and Negative Binomial GLMs
================
Ethan Bensadoun
2024-08-15

## (a)

In terms of the problem, we note the following:

- y_i denotes the response for the subject i

- x_i = 1 for blacks, and x_i =0 for whites.

``` r
Homicides$race <- as.numeric(as.character(Homicides$race))
Homicides$count <- as.numeric(as.character(Homicides$count))

Homicides_fit <- glm(count ~ race, data = Homicides, family = poisson)
summary(Homicides_fit)
```

    ## 
    ## Call:
    ## glm(formula = count ~ race, family = poisson, data = Homicides)
    ## 
    ## Coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept) -2.38321    0.09713  -24.54   <2e-16 ***
    ## race         1.73314    0.14657   11.82   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for poisson family taken to be 1)
    ## 
    ##     Null deviance: 962.80  on 1307  degrees of freedom
    ## Residual deviance: 844.71  on 1306  degrees of freedom
    ## AIC: 1122
    ## 
    ## Number of Fisher Scoring iterations: 6

``` r
Beta1_hat <- coef(Homicides_fit)["race"]
```

The coefficient of interest in this case is the one for race, i.e.,
Beta1_hat = 1.73314. Note that for the Poisson distribution, the
canonical link function is the log. The log link function, in this
context, represents the difference of logs between blacks and whites in
the expected counts of homicide victims in the 12 month period.

``` r
# Now we exponentiate the coefficient in order to transform the log-difference 
# back to the original scale
exp_Beta1_hat <- exp(Beta1_hat)
```

Thus, we can interpret exp_Beta1_hat = 5.66, as saying the following:

- Within the past 12 months, blacks are 5.66 more likely to have known
  someone who was the victim of a homicide.

## (b)

- Since we have fit a Poisson GLM, an inadequacy that may come up stems
  from the fact that the mean and variance are equal in this type of
  distribution. Thus, overdispersion (phi \>1) will occur when the
  variance is greater than the mean. This will be due to heterogeneity
  in the data.
- There may be some auxiliary variables that are unaccounted for in the
  data causing an increased amount of variability. Thus, due to
  variables that are not accounteed for, the variance in the number of
  known homicide victims may be different across individuals, which
  cannot be explained by Beta1, i.e., race.
- Some of these hetergenous factors may include: age, gender, geographic
  location, education level, etc…

``` r
Homicides_fit_nb <- glm.nb(count ~ race, data = Homicides)
summary(Homicides_fit_nb)
```

    ## 
    ## Call:
    ## glm.nb(formula = count ~ race, data = Homicides, init.theta = 0.2023119205, 
    ##     link = log)
    ## 
    ## Coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)  -2.3832     0.1172 -20.335  < 2e-16 ***
    ## race          1.7331     0.2385   7.268 3.66e-13 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(0.2023) family taken to be 1)
    ## 
    ##     Null deviance: 471.57  on 1307  degrees of freedom
    ## Residual deviance: 412.60  on 1306  degrees of freedom
    ## AIC: 1001.8
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  0.2023 
    ##           Std. Err.:  0.0409 
    ## 
    ##  2 x log-likelihood:  -995.7980

``` r
theta <- 0.2023
# the estimated dispersion parameter for a neg.bin model can be calculated by 
# the following formula: 1/theta = alpha
alpha <- 1/theta
race_std_error <- 0.2385
```

Since, under this parametrization overdispersion relative to the Poisson
is indicated by theta \> 0. Note that the greater theta is, the greater
the overdispersion is relative to the Poisson. In this case, we have
that alpha is 4.943154 which is also much greater than 0. The standard
errors for the negative binomial are also much larger. That is, the std
error for Beta1 hat coeff, is 0.2385, whereas for the poisson glm it is
0.14657. The fact that the std errors are larger for the negative
binomial is more realistic. In addition the AIC score is smaller
(1001.8) when compared to the Poisson glm AIC score (1122).

Overall, assuming a Poisson distribution to analyze this data is not
consistent with the results we got after fitting the data of in both
models. In conclusion, there is indeed overdispersion in the Poisson GLM
model, which is proven by fitting and analyzing the negative binomial
glm model whch is much better at capturing the extra variability in the
data that the Poisson glm cannot.

\*\* maybe need to add part about: estimate how the variance depends on
the mean. \*\*

## (c)

``` r
# Here we can find the 95 confidence intervals for both models:
exp(confint(Homicides_fit)["race", ])
```

    ## Waiting for profiling to be done...

    ##    2.5 %   97.5 % 
    ## 4.236333 7.532534

``` r
exp(confint(Homicides_fit_nb)["race", ])
```

    ## Waiting for profiling to be done...

    ##    2.5 %   97.5 % 
    ## 3.577823 9.131638

``` r
# Note that we need to exponentiate the confidence intervals because we want to 
# have the log rate ratio. 
```

As a result, both confidence intervals give you estimated range of
values that is believed to contain the true rate ratio with 95%
confidence. The exponeniated CIs represent the ratio of of the mean
count of known murder victims for blacks to whites (beta1 hat coeffs for
both models.). In both models the beta1_hat coeff is equal to 1.7331.
Thus, this falls in the 95% CI in both cases. We can conclude that,
since there is overdispersion in the data, as proven in part (b), the
negative binomial GLM provides a more accurate interval because it
accounts for this extra variability, despite being wider. Even though,
the Poisson glm C.I is narrower, we can say that it underestimates the
true variability, thus, it is not the good choice in this case.

## (d)

Now we need to fit the models using the QL method:

``` r
quasi_Homicides_fit <- glm(count ~ race, data = Homicides, family = quasipoisson)
summary(quasi_Homicides_fit)
```

    ## 
    ## Call:
    ## glm(formula = count ~ race, family = quasipoisson, data = Homicides)
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  -2.3832     0.1283  -18.57   <2e-16 ***
    ## race          1.7331     0.1937    8.95   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for quasipoisson family taken to be 1.745694)
    ## 
    ##     Null deviance: 962.80  on 1307  degrees of freedom
    ## Residual deviance: 844.71  on 1306  degrees of freedom
    ## AIC: NA
    ## 
    ## Number of Fisher Scoring iterations: 6

``` r
quasi_Homicides_fit_nb = glm(count ~ race, family=quasi(link="log",variance="mu^2"),data=Homicides)
summary(quasi_Homicides_fit_nb)
```

    ## 
    ## Call:
    ## glm(formula = count ~ race, family = quasi(link = "log", variance = "mu^2"), 
    ##     data = Homicides)
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  -2.3832     0.1200 -19.861  < 2e-16 ***
    ## race          1.7331     0.3442   5.036 5.43e-07 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for quasi family taken to be 16.54457)
    ## 
    ##     Null deviance: 1839.0  on 1307  degrees of freedom
    ## Residual deviance: 1870.9  on 1306  degrees of freedom
    ## AIC: NA
    ## 
    ## Number of Fisher Scoring iterations: 10

We have the same coefficient for race in both cases, as previously seen.
However, the standard errors are much larger in the quasi-likelihood
models which is indicative of a higher variability due to overdispersion
in the data. Since the model adjusts itself due to overdispersion the
result is indeed larger standard errors. In addition, we can see that
the AIC scores in both cases is NA. This is expected in this case
because there isn’t a full likelihood function (unlike in the previous
cases). We can also see that the dispersion paramter for the quasi
family is taken to be 16.54457. These are much larger, which araises
from the fact that it is harder to model.
