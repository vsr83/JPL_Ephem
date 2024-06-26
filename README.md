# JPL_Ephem
Attempt to reimplement the numerical integration used to construct JPL Ephemerides.

The implementation work is still in early stages. The implementation is heavily inspired by DE118i written by Steve Moshier [3]. Most of the routines have been tested by comparisons to DE118i.

The numerical integration is performed for the acceleration equations from [1], [2] for point masses in the isotropic, parametrized post-Newtonian (PPN) n-body metric taking into account relativistic effects. In addition, figure effects from the Moon and Earth oblateness are taken into account as well as the effect of Earth tides.

The integration scheme used is the 8th order Adams-Bashforth-Moulton method with the initial 8 steps computed with 4th order Runge-Kutta. From testing, it seems that using timesteps shorter than 0.05 days (1.2 hours) does not improve the accuracy.

The current implementation has not been written with performance in mind. The Moon libration and oblateness computations use coefficients taken from several different sources and do not take into account the tidal and spin distortions of the Moon. Also, unlike in the references, the OSV for the Sun is included in the integration. The integration also relies on 64-bit floating point arithmetic.

The position error w.r.t. JPL Horizons over 100-year integration with 864-second timesteps is shown below for a integration with and without 67 major asteroids, respectively.
[![Integration error with asteroids](error_with_asteroids.png)](error_with_asteroids.png)
[![Integration error without asteroids](error_without_asteroids.png)](error_without_asteroids.png)


## References
1. Newhall, Standish, Williams - DE 102: a numerically integrated ephemeris of the Moon and planets spanning forty-four centuries, Astronomy and Astrophysics, 125, 150-167, 1983 [link](https://adsabs.harvard.edu/full/1983A%26A...125..150N).
2. Urban, Seidelmann - Explanatory Supplement to the Astronomical Almanac, 3rd edition, University Science Books, 2013. [link](https://www.amazon.com/Explanatory-Supplement-Astronomical-Almanac-Urban/dp/1891389858)
3. Steve Moshier, DE118i available at [link](http://www.moshier.net/de118i-2.zip).
