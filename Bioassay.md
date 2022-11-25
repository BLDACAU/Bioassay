---
author:
- Jooyoung Lee
editor: visual
title: "Statistical Analysis of Bioassay Data: Data Example"
toc-title: Table of contents
---

## R packages

::: cell
``` {.r .cell-code}
## install packages and loading
#install.packages("drc")
library(drc)
```


## Data set

::: cell
``` {.r .cell-code}
# load data
dat <- read.table("~/Dropbox/CAU/Workshop/2022/Bioassay_112522/dat.txt", header=T)
head(dat)
```

::: {.cell-output .cell-output-stdout}
      response dose preparations
    1      161 1.00            S
    2      151 1.00            T
    3      162 1.50            T
    4      194 2.25            S
    5      176 1.50            S
    6      193 2.25            T
:::
:::

## Parallel Line Model

::: cell
``` {.r .cell-code}
dat$preparations <- factor(dat$preparations, levels=c("S", "T"))
dat$dose <- as.numeric(log(dat$dose))
dat$response <- as.numeric(dat$response)
dat$doselevel <- factor(dat$dose)

model1 <- lm(response ~  preparations + dose + dose*preparations + doselevel*preparations, data=dat)

# preparations: preparations
# dose: linear regression 
# dose*preparations: Non-parallelism
# doselevel*preparations: Non-linearity

anova(model1)
```

::: {.cell-output .cell-output-stdout}
    Analysis of Variance Table

    Response: response
                           Df Sum Sq Mean Sq  F value    Pr(>F)    
    preparations            1   11.1    11.1   0.2258    0.6381    
    dose                    1 7958.5  7958.5 161.7175 1.292e-13 ***
    doselevel               1   31.2    31.2   0.6339    0.4322    
    preparations:dose       1   73.5    73.5   1.4931    0.2313    
    preparations:doselevel  1    5.4     5.4   0.1092    0.7434    
    Residuals              30 1476.4    49.2                       
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
:::

``` {.r .cell-code}
model2 <- lm(response ~  preparations + dose, data=dat)
summary(model2)
```

::: {.cell-output .cell-output-stdout}

    Call:
    lm(formula = response ~ preparations + dose, data = dat)

    Residuals:
         Min       1Q   Median       3Q      Max 
    -13.4444  -4.6944  -0.3064   3.3399  22.3140 

    Coefficients:
                  Estimate Std. Error t value Pr(>|t|)    
    (Intercept)    159.686      2.095   76.22  < 2e-16 ***
    preparationsT   -2.103      2.312   -0.91     0.37    
    dose            44.053      3.424   12.87 2.08e-14 ***
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    Residual standard error: 6.933 on 33 degrees of freedom
    Multiple R-squared:  0.834, Adjusted R-squared:  0.8239 
    F-statistic: 82.89 on 2 and 33 DF,  p-value: 1.356e-13
:::

``` {.r .cell-code}
# the log of relative potency
logrho <- model2$coefficients[2]/model2$coefficients[3]
logrho
```

::: {.cell-output .cell-output-stdout}
    preparationsT 
      -0.04774804 
:::

``` {.r .cell-code}
# the relative potency
exp(logrho)
```

::: {.cell-output .cell-output-stdout}
    preparationsT 
         0.953374 
:::
:::

## Slope ratio Model

::: cell
``` {.r .cell-code}
model4 <- lm(response ~   dose + dose:preparations , data=dat)

# the log of relative potency
logrho <- (model4$coefficients[2]+model4$coefficients[3])/model4$coefficients[2]
logrho
```

::: {.cell-output .cell-output-stdout}
        dose 
    1.012457 
:::

``` {.r .cell-code}
# the relative potency
exp(logrho)
```

::: {.cell-output .cell-output-stdout}
        dose 
    2.752356 
:::
:::

## The four-parameter logistic model

::: cell
``` {.r .cell-code}
data(acidiq)
head(acidiq)
```

::: {.cell-output .cell-output-stdout}
      dose pct       rgr
    1    0 999 0.2904680
    2    0 999 0.2834756
    3    0 999 0.2925334
    4    0 999 0.3141099
    5    0 999 0.3095383
    6    0 999 0.3199928
:::

``` {.r .cell-code}
## same upper, lower asymptote, and slope factor model
acidiq.model1 <- drm(rgr ~ dose, pct, data = acidiq, fct = LL.4(), pmodels = list(~1, ~1, ~1, ~factor(pct) - 1))
```

::: {.cell-output .cell-output-stdout}
    Control measurements detected for level: 999
:::

``` {.r .cell-code}
modelFit(acidiq.model1)
```

::: {.cell-output .cell-output-stdout}
    Lack-of-fit test

              ModelDf      RSS Df F value p value
    ANOVA         123 0.023854                   
    DRC model     170 0.056393 47  3.5698  0.0000
:::

``` {.r .cell-code}
summary(acidiq.model1)
```

::: {.cell-output .cell-output-stdout}

    Model fitted: Log-logistic (ED50 as parameter) (4 parms)

    Parameter estimates:

                    Estimate Std. Error  t-value   p-value    
    b:(Intercept) 2.0513e+00 1.1455e-01  17.9079 < 2.2e-16 ***
    c:(Intercept) 2.9822e-02 3.3819e-03   8.8182 1.337e-15 ***
    d:(Intercept) 3.0083e-01 2.8216e-03 106.6169 < 2.2e-16 ***
    e:100         3.2329e+02 2.1757e+01  14.8594 < 2.2e-16 ***
    e:83          3.7806e+02 2.3425e+01  16.1394 < 2.2e-16 ***
    e:67          4.9157e+02 2.8812e+01  17.0613 < 2.2e-16 ***
    e:50          5.1889e+02 2.9700e+01  17.4709 < 2.2e-16 ***
    e:33          5.3173e+02 3.1125e+01  17.0836 < 2.2e-16 ***
    e:17          3.8530e+02 2.1690e+01  17.7640 < 2.2e-16 ***
    e:0           3.5125e+02 2.0075e+01  17.4969 < 2.2e-16 ***
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    Residual standard error:

     0.01821322 (170 degrees of freedom)
:::

``` {.r .cell-code}
acidiq.model2 <- drm(rgr ~ dose, pct, data = acidiq, fct = LL.4(),
pmodels = list(~factor(pct), ~factor(pct), ~factor(pct), ~factor(pct) - 1))
```

::: {.cell-output .cell-output-stdout}
    Control measurements detected for level: 999
:::

``` {.r .cell-code}
summary(acidiq.model2)
```

::: {.cell-output .cell-output-stdout}

    Model fitted: Log-logistic (ED50 as parameter) (4 parms)

    Parameter estimates:

            Estimate Std. Error t-value   p-value    
    b:100 1.2787e+00 1.3309e-01  9.6075 < 2.2e-16 ***
    b:83  1.5339e+00 2.0639e-01  7.4321 7.222e-12 ***
    b:67  1.8407e+00 2.2199e-01  8.2916 5.508e-14 ***
    b:50  2.3613e+00 2.7409e-01  8.6151 8.402e-15 ***
    b:33  2.8670e+00 3.5189e-01  8.1475 1.266e-13 ***
    b:17  2.0412e+00 2.5656e-01  7.9560 3.792e-13 ***
    b:0   2.5446e+00 3.4011e-01  7.4817 5.485e-12 ***
    c:100 9.6130e-03 1.0263e-02  0.9367 0.3504093    
    c:83  1.2524e-02 1.0168e-02  1.2317 0.2199561    
    c:67  1.9499e-02 8.8396e-03  2.2059 0.0288907 *  
    c:50  2.6577e-02 7.0401e-03  3.7751 0.0002288 ***
    c:33  5.7116e-02 6.0049e-03  9.5117 < 2.2e-16 ***
    c:17  2.8553e-02 6.8025e-03  4.1974 4.577e-05 ***
    c:0   3.4560e-02 6.1532e-03  5.6165 9.028e-08 ***
    d:100 2.9184e-01 4.1646e-03 70.0756 < 2.2e-16 ***
    d:83  3.0365e-01 9.6544e-03 31.4519 < 2.2e-16 ***
    d:67  3.1383e-01 7.3778e-03 42.5379 < 2.2e-16 ***
    d:50  2.9464e-01 5.9817e-03 49.2564 < 2.2e-16 ***
    d:33  3.0838e-01 5.7283e-03 53.8353 < 2.2e-16 ***
    d:17  3.2815e-01 8.6033e-03 38.1420 < 2.2e-16 ***
    d:0   2.9641e-01 6.6422e-03 44.6253 < 2.2e-16 ***
    e:100 3.8983e+02 3.7866e+01 10.2951 < 2.2e-16 ***
    e:83  4.1741e+02 3.5719e+01 11.6859 < 2.2e-16 ***
    e:67  4.8677e+02 3.3648e+01 14.4663 < 2.2e-16 ***
    e:50  5.4376e+02 3.1871e+01 17.0609 < 2.2e-16 ***
    e:33  4.4051e+02 2.3735e+01 18.5596 < 2.2e-16 ***
    e:17  3.3975e+02 2.0353e+01 16.6931 < 2.2e-16 ***
    e:0   3.4727e+02 1.9458e+01 17.8470 < 2.2e-16 ***
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    Residual standard error:

     0.0149249 (152 degrees of freedom)
:::

``` {.r .cell-code}
anova(acidiq.model2, acidiq.model1)
```

::: {.cell-output .cell-output-stdout}

    1st model
     fct:     LL.4()
     pmodels: ~1, ~1, ~1, ~factor(pct) - 1
    2nd model
     fct:     LL.4()
     pmodels: ~factor(pct), ~factor(pct), ~factor(pct), ~factor(pct) - 1
:::

::: {.cell-output .cell-output-stdout}
    ANOVA table

              ModelDf      RSS Df F value p value
    2nd model     170 0.056393                   
    1st model     152 0.033858 18  5.6201  0.0000
:::
:::
