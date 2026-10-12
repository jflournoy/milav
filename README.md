[![DOI](https://zenodo.org/badge/129271129.svg)](https://zenodo.org/badge/latestdoi/129271129)

# milav: Measurement Invariance in lavaan

This package exists mostly to collect functions I wrote for performing measurement invariance model comparisons for the Life and Time study.

This package also showcases the [results of simulations of changes in fit indices](https://jflournoy.github.io/milav/articles/simulate-mi-tests.html) in order to evaluate reasonable cutoffs based on the 4-group, 4-wave structure of the data in that study.

## What the simulation shows

Invariance tests usually compare a constrained model to a freer one and reject invariance when a fit index drops by more than a fixed amount, most often the cutoffs from Cheung and Rensvold (2002): ΔCFI ≤ −.01, ΔMcDonald's noncentrality index (MFI) ≤ −.02, and Δγ̂ ≤ −.001. Those cutoffs were not derived for this design, and later work showed that fit-index cutoffs shift with measurement quality (Kang et al., 2016; McNeish et al., 2018). So I simulated the null directly.

**Design.** Data were generated from a single-factor model that is strictly invariant across 4 groups (n = 250 each). Each dataset was fit twice, with loadings, intercepts, and residual variances either free across groups or constrained equal, and the change in fit was recorded. Because invariance holds in the generating model, the lower tail of each distribution is the cutoff that would give a known false-rejection rate. Fifteen conditions crossed 4, 6, or 8 items with loadings drawn uniformly within bands from .50–.59 to .90–.99, with 1,000 replications each.

**Results.** Each range below is the 1st percentile under the null, across the 15 conditions.

- **ΔCFI depends heavily on measurement quality.** Its 1st percentile runs from −.0127 (4 items, loadings .50–.59) to −.0011 (strong loadings, 6–8 items). With few weak items, the −.01 cutoff rejects a true invariance hypothesis more than 1% of the time. With strong items the null distribution sits an order of magnitude inside the cutoff, so a fixed −.01 is not calibrated to anything.
- **ΔMFI and Δγ̂ are far more stable.** ΔMFI runs from −.0153 to −.0074, and Δγ̂ from −.0038 to −.0018. The −.02 MFI cutoff is conservative in every condition. The −.001 γ̂ cutoff is too strict: the 5th percentile of Δγ̂ fell below −.001 in every condition, so it rejects true invariance more than 5% of the time throughout.
- **Adjusted γ̂ is not usable as a difference test here.** Its median drifts positive at loadings of .80 and above.

Pooled over all conditions, the 1st percentiles are ΔMFI −.012, ΔCFI −.006, and Δγ̂ −.0026.

**Scope.** These numbers apply to this design: one factor, 4 groups of 250, and the strict constraint set described above. Other sample sizes, group counts, or model structures will move the distributions, so rerun the vignette with your own structure rather than borrowing these cutoffs.

## Status

Written in 2018 for the Life and Time Study (see Thalmayer et al., 2019) and not updated since. It has not been tested against current lavaan releases.

```r
remotes::install_github("jflournoy/milav")
```

## Citation

Flournoy, J. C. (2018). *milav: Measurement invariance tests in lavaan* (Version 0.1.0-alpha) [R package]. Zenodo. https://doi.org/10.5281/zenodo.1294089

## References

Cheung, G. W., & Rensvold, R. B. (2002). Evaluating goodness-of-fit indexes for testing measurement invariance. *Structural Equation Modeling, 9*(2), 233–255. https://doi.org/10.1207/S15328007SEM0902_5

Kang, Y., McNeish, D. M., & Hancock, G. R. (2016). The role of measurement quality on practical guidelines for assessing measurement and structural invariance. *Educational and Psychological Measurement, 76*(4), 533–561. https://doi.org/10.1177/0013164415603764

McNeish, D., An, J., & Hancock, G. R. (2018). The thorny relation between measurement quality and fit index cutoffs in latent variable models. *Journal of Personality Assessment, 100*(1), 43–52. https://doi.org/10.1080/00223891.2017.1281286

Thalmayer, A. G., Saucier, G., Srivastava, S., Flournoy, J. C., & Costello, C. K. (2019). Ethics-relevant values in adulthood: Longitudinal findings from the Life and Time Study. *Journal of Personality, 87*(6), 1119–1135. https://doi.org/10.1111/jopy.12462
