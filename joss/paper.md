---
title: "OrbDot: A Python package for studying the secular evolution of exoplanet orbits"
tags:
  - Python
  - astronomy
  - exoplanets
  - exoplanet transits
  - exoplanet radial velocities
  - orbital evolution
  - nested sampling
  - model fitting
authors:
  - name: Simone R. Hagey
    orcid: 0000-0001-8072-0590
    affiliation: 1
  - name: Aaron Boley
    orcid: 0000-0002-0574-4418
    affiliation: 1
affiliations:
  - name: Department of Physics and Astronomy, The University of British Columbia, 6224 Agricultural Road Vancouver, BC, V6T 1Z1, Canada
    index: 1
date: 1 September 2025
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
aas-doi: 10.3847/1538-3881/aded15
aas-journal: Astronomical Journal
---

# Summary

Gradual changes in exoplanet orbits, known as secular variations, can be detected through observations of transits, eclipses, and radial velocities that span multiple decades in time. Their detection and characterization enable the study of a wide range of dynamical phenomena, such as orbital decay and precession, which operate on timescales of millions of years. Under certain conditions, measurements of secular variations can even probe the interior structure of exoplanets, providing a unique tool for understanding exoplanet formation and evolution.

The necessity to search over many orbital epochs coupled with an ever-growing archive of exoplanet observations creates a need for fast and flexible open-source software that can reliably detect gradual changes in exoplanet orbits. `OrbDot` addresses this need by offering robust tools for fitting secular evolution models to exoplanet transit and eclipse mid-times, transit durations, and radial velocity data.

A key advantage of `OrbDot` is its ability to fit multiple types of data simultaneously, which can help to break parameter degeneracies. It also excels at assisting in result interpretation by generating reports on model comparisons and assessments of various physical effects in the context of the models and their corresponding theory. For example, analysis reports could determine key parameters for assessing tidal energy dissipation, apsidal precession mechanisms, variations due to systemic proper motion, and the dynamical effects of non-resonant companion objects, depending on the applied models.

`OrbDot` remains highly efficient with multiple data types and a high number of free parameters, as it utilizes powerful nested sampling algorithms of the `nestle` [@nestle; @Skilling:2006] and PyMultiNest [@Buchner:2014; @Feroz:2009] packages. The intricacies of the implementation are abstracted such that the `OrbDot` input files are simple and the method calls require only a list of free parameters, along with the desired model for fitting.

Extensive documentation, including examples, is hosted on [ReadTheDocs](https://orbdot.readthedocs.io).

The examples demonstrate that `OrbDot` can quickly reproduce literature results using only a few lines of code. Readers may be especially interested in the `OrbDot` example analysis of the transit and eclipse mid-times of Hot Jupiter WASP-12 b, which is well-known for showing strong evidence for orbital decay.

A complementary case study of TrES-1 b [@Hagey:2025] illustrates the full capabilities of `OrbDot`, placing it in a broader scientific context. Moreover, an early version of this code was used for the orbital analysis of the Hot Neptune LTT-9779 b, published in @Edwards:2023.

# Statement of need

Many exoplanet systems now have transit and radial velocity data spanning over a decade, enabling studies of secular variations. While tools for analyzing short-term transit variations exist, there is a lack of open-source software dedicated to long-term orbital evolution. This does not, however, reflect a lack of interest, as the number of such studies is growing rapidly.

`OrbDot` lowers the barrier to entry for researchers at all levels, including undergraduates, by making advanced statistical methods accessible without requiring extensive computational experience. Despite its ease of use, `OrbDot` is not intended to be a black box. Rather, with extensive documentation, examples, and accessible source code, it is presented to the community with transparency that lends itself to community contributions and independent verification of results. It is designed to be easily extended, as the nested sampling framework supports custom log-likelihood models with free parameters that are part of the `OrbDot` ecosystem. This ensures that the software may evolve to meet the needs of the research community.

# Similar software

Some existing tools have features that overlap with `OrbDot`’s functionality, but none provides its full suite of capabilities. The most similar codes, `Susie` [@susie; @Barker:2024] and `PdotQuest` [@Wang:2024], which are not fully packaged software, are designed to fit secular evolution models to transit and eclipse timing data, but not to radial velocities. `OrbDot` also has greater flexibility in model fitting than `Susie` and `PdotQuest`. For example, `Susie` employs simple least-squares fitting, and while `PdotQuest` [@Wang:2024] uses MCMC, it currently supports only the orbital decay model, making both codes narrower in scope than `OrbDot`. Moreover, `OrbDot` includes tools for theoretical interpretation.

The well-known package `TTVFast` [@Deck:2014] is highly robust and capable of modeling both transit and radial velocity data, but it is focused on short-term timing variations driven by multi-planet dynamics near mean-motion resonances. `OrbDot`, in contrast, explores long-term secular models. Similarly, general-purpose frameworks such as `juliet` [@Espinoza:2019], `exoplanet` [@ForemanMackey:2021], `EXOFAST` [@Eastman:2019], `ExoStriker` [@Trifonov:2019], and `allesfitter` [@Gunther:2021] can jointly model transit and RV data, but their TTV models are also restricted to short-term variations, with an emphasis on transit light curve fitting rather than directly using transit and eclipse mid-times to constrain models.

The RV-focused packages `RadVel` [@Fulton:2018] and `Kima` [@Faria:2018] provide flexible modeling of radial velocity datasets, but do not have a framework for incorporating transit and eclipse data into the modeling. The codes also have limited options for studying long-term trends. For example, `RadVel` [@Fulton:2018] includes linear and quadratic RV terms, which `OrbDot` also supports, but it does not model evolving orbital elements. `Kima` [@Faria:2018] incorporates an apsidal precession model, but only for circumbinary planets -- a niche application.

In summary, `OrbDot` is a fully packaged, documented, and maintained software suite designed specifically for studies of secular orbital evolution. It unifies transit, eclipse, and RV data, bypasses light-curve fitting to focus directly on timing measurements, implements robust Bayesian inference with nested sampling, and provides a flexible framework for selecting models, priors, and parameterizations. In addition, it integrates tools for theoretical analysis, making `OrbDot` the first software to systematically support data-driven, comprehensive studies of long-term exoplanet orbital dynamics.

# Acknowledgements

This work was supported, in part, by an NSERC Discovery Grant (DG-2020-04635) and the University of British Columbia. SH's contribution was further supported, in part, by an NSERC PGS-D and a Li Tze Fong Fellowship.

# References
