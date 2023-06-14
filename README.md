# MST_uncured
R-code for "A two-sample comparison of mean survival times of uncured sub-populations"

Comparing the survival times among two groups is a common problem in time-to-event analysis, for example if one would like to understand
whether one medical treatment is superior to another. In the standard survival analysis setting, there has been a lot of discussion on how to
quantify such difference and what can be an intuitive, easily interpretable, summary measure. In the presence of subjects that are immune to the event
of interest (‘cured’), we illustrate that it is not appropriate to just compare the overall survival functions. Instead, it is more informative to compare the
cure fractions and the survival of the uncured sub-populations separately from each other. Our research is mainly driven by the question: if the cure
fraction is similar for two available treatments, how else can we determine which is preferable? To this end, we propose the use of the mean survival times  of
the uncured subpopulations of both  groups (MSTu). 

We first consider a nonparametric approach and estimate the MSTu for both groups.  Secondly, we consider a semiparametric approach to adjust the mean survival time of the uncured for potential confounders. For each group we employ the widely used logistic-Cox mixture cure model and estimate the MSTu conditional on a given covariate value. An asymptotic and
a permutation-based method have been developed for construction of confidence intervals and hypothesis testing on the difference in the MSTu and conditional MSTu respectively.
The methods are illustrated with medical data of leukemia patients.

For more details about the method, the asymptotic properties and the finite sample behaviour see Dobler and Musta (2023a,2023b).

Here we demonstrate the method for the melanoma data (ECOG phase III clinical trial e1684) from the smcure package (Cai et al. (2012)). The purpose of this study was to evaluate the effect of treatment (high dose interferon alpha-2b regimen) as the postoperative adjuvant therapy. The event time is the time from initial treatment to recurrence of melanoma and three covariates have been considered: age (continuous variable centered to the mean), gender (0=male and 1=female) and treatment (0=control and 1=treatment). The data consists of 284 observations (after deleting missing data) out of which 196 had recurrence of the melanoma cancer (around 30% censoring). The Kaplan Meier estimator of the survival function can be found in 'KM_estimator.png'.

One can run the example via '2step_logCox-MCM.R'. This produces estimators for the parameters of the incidence, the regression coefficients and the baseline survival of the Cox component (latency). The bootstrap procedure is used to estimate the variance of the parameter estimates and the resulting p-values. The R file uses the packages 'smcure' for the data and 'np' for bandwidth selection. Additional functions are defined in the file 'functions_2step_logCox-MCM.R'.

References

- D. Dobler and E. Musta. A two-sample comparison of mean survival times of uncured sub-populations – Part I: Nonparametric analyses
- D. Dobler and E. Musta. A two-sample comparison of mean survival times of uncured sub-populations – Part II: Semiparametric analyses 
