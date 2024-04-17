# JPL_Ephem
Attempt to reimplement the numerical integration used to construct JPL Ephemerides.

The implementation work is still in early stages. The implementation is heavily inspired by DE118i written by Steve Moshier [3].

The implementation is missing a proper integration scheme. In addition, the implementation has not been written for performance. The Moon libration and oblateness computations use coefficients taken from several different sources. Also, unlike in the references, the OSV for the Sun is included in the integration.

The position error w.r.t. JPL Horizons over 100-year integration with 864-second timesteps is shown below for a integration with and without 67 major asteroids, respectively.
[![Integration error with asteroids](error_with_asteroids.png)](error_with_asteroids.png)
[![Integration error without asteroids](error_without_asteroids.png)](error_without_asteroids.png)


## References
1. Newhall, Standish, Williams - DE 102: a numerically integrated ephemeris of the Moon and planets spanning forty-four centuries, Astronomy and Astrophysics, 125, 150-167, 1983 [link](https://adsabs.harvard.edu/full/1983A%26A...125..150N).
2. Urban, Seidelmann - Explanatory Supplement to the Astronomical Almanac, 3rd edition, University Science Books, 2013. [link](https://www.amazon.com/Explanatory-Supplement-Astronomical-Almanac-Urban/dp/1891389858)
3. Steve Moshier, DE118i available at [link](http://www.moshier.net/de118i-2.zip).
