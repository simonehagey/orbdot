---
title: "``OrbDot``: A Python package for studying the secular evolution of exoplanet orbits"
tags:
  - python
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
    orcid: 0000-0000-0000-0000
    affiliation: 1
affiliations:
  - name: Department of Physics and Astronomy, The University of British Columbia, 6224 Agricultural Road Vancouver, BC V6T 1Z1, Canada
    index: 1
date: 4 April 2025
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

Gradual changes in exoplanet orbits can be detected through observations that span multiple decades in time. Called secular variations, these changes manifest as long-term trends in observed transits, eclipses, and radial velocities. Their detection and characterization enable the study of a wide range of dynamical phenomena, such as orbital decay and precession, which operate on timescales of millions of years. Under certain conditions, measurements of secular variations can even probe the interior structure of exoplanets, providing a unique tool for understanding exoplanet formation, evolution, and dynamics.

The necessity to search over long periods of time coupled with an ever-growing archive of exoplanet observations creates a need for fast and flexible open-source software that can reliably detect gradual changes in exoplanet orbits. `OrbDot` addresses this need by offering a robust set of tools for fitting secular evolution models to exoplanet transit and eclipse mid-times, transit durations, and radial velocity data.

A key advantage of `OrbDot` is its ability to fit multiple data types simultaneously, which can help to break parameter degeneracies. `OrbDot` also excels at assisting in result interpretation by generating automatic reports of model comparisons, as well as assessing various physical effects in the context of models and their corresponding theory. For example, analysis reports could determine key parameters for assessing tidal energy dissipation, apsidal precession mechanisms, variations due to systemic proper motion, and non-resonant companion objects, depending on the applied models.

`OrbDot` remains highly efficient with multiple data types and a high number of free parameters, as it utilizes the powerful nested sampling algorithms of the `Nestle` [@nestle; @Skilling:2006] and `PyMultiNest` [@pymultinest; @Buchner:2014; @Feroz:2009] packages. The intricacies of the implementation are abstracted such that the `OrbDot` input files are simple and the method calls require only a list of free parameters, along with the desired model for fitting.

Extensive documentation, including examples, is hosted on [ReadTheDocs](https://`orbdot.readthedocs.io).

The examples demonstrate that `OrbDot` can reproduce literature results quickly using only a few lines of code. Readers may be especially interested in the `OrbDot` example analysis of the transit and eclipse mid-times of Hot Jupiter WASP-12 b, which is well-known for showing strong evidence for orbital decay.
  
In addition to the documentation, there is a complementary case study of Hot Jupiter TrES-1 b (CITE) that demonstrates the full capabilities of the `OrbDot` model fitting functions and interpretive tools, placing the need for this software in deeper scientific context. The TrES-1 b example may further leave the reader with ideas for how to utilize `OrbDot` in their own work. Moreover, an early version of this code was used for the orbital analysis of the Hot Neptune LTT-9779 b, published in @Edwards:2023.

# Statement of need

Many star-planet systems now have transit and radial velocity observations that span ten years or more. The exoplanet field is also still in a phase of significant growth, especially for population-level studies. Streamlined tools are greatly needed to enable access to the growing data and to foster individual as well as group research efforts. Although several Python packages exist for analyzing short-term transit variations, there is a lack of open-source tools for analyzing secular variations. This does not, however, reflect a lack of interest, as the number of such studies is growing rapidly.

The primary purpose of releasing this package is to assist researchers with their data analysis so that they may focus on the science. It makes advanced statistical methods accessible even to those with limited computational experience, lowering the barrier to entry for engaging in this research. It is suitable for researchers at all stages, including undergraduates.

Despite its ease of use, `OrbDot` is not intended to be a black box. Rather, with extensive documentation, examples, and accessible source code, it is presented to the community with transparency that lends itself to community contributions and independent verification of results. It is designed to be easily extendable, as the nested sampling algorithms are implemented in a way that the model fitting routines can apply to any log-likelihood with free parameters that are part of the `OrbDot` ecosystem. This ensures that the tools may evolve to meet the needs of the research community.

# Acknowledgements

This work was supported, in part, by an NSERC Discovery Grant (DG-2020-04635) and the University of British Columbia. SH's contribution was further supported, in part, by an NSERC PGS-D and a Li Tze Fong Fellowship.

# References