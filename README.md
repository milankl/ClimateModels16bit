[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4015133.svg)](https://doi.org/10.5281/zenodo.4015133)
# Number formats, error mitigation and scope for 16-bit arithmetics in weather and climate modelling analysed with a shallow water model

M Klöwer¹, PD Düben² and TN Palmer¹, 2020.

¹Atmospheric, Oceanic and Planeteray Physics, University of Oxford, UK.   
²European Centre for Medium-Range Weather-Forecasts, Reading, UK.

accepted in [Journal of Advances in Modeling the Earth System (JAMES)](https://agupubs.onlinelibrary.wiley.com/journal/19422466). This repository contains scripts to reproduce the results in this publication.

Additional software that was used
- [ShallowWaters.jl](https://github.com/milankl/ShallowWaters.jl)
- [SoftPosit.jl](https://github.com/milankl/SoftPosit.jl)
- [BFloat16s.jl](https://github.com/JuliaComputing/BFloat16s.jl)
- [Float8s.jl](https://github.com/milankl/Float8s.jl)

## Abstract

The need for high precision calculations with 64-bit or 32-bit floating-point
arithmetic for weather and climate models is questioned.
Lower precision numbers can accelerate simulations and are increasingly supported
by modern computing hardware. This paper investigates the potential of 16-bit
arithmetic when applied within a shallow water model that serves as a medium
complexity weather or climate application. There are several 16-bit number
formats that can potentially be used (IEEE half precision, BFloat16, posits,
integer and fixed-point). It is evident that a simple change to 16-bit arithmetic
will not be possible for complex weather and climate applications as it will
degrade model results by intolerable rounding errors that cause a stalling of
model dynamics or model instabilities. However, if the posit number format is
used as an alternative to the standard floating-point numbers the model degradation
can be significantly reduced. Furthermore, a number of mitigation methods,
such as rescaling, reordering and mixed-precision, are available to make model
simulations resilient against a precision reduction. If mitigation methods are applied,
16-bit floating-point arithmetic can be used successfully within the shallow water model.
The results show the potential of 16-bit formats for at least parts of complex weather
and climate models where rounding errors would be entirely masked by initial
condition, model or discretization error.
