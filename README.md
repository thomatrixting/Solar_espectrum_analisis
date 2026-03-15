# Solar_espectrum_analisis
This project is an analysis of the light spectrum emitted by the Sun based on data from SOURCE mission [1], with the intent to estimate de sun surface temperature ussing plank ecuations.

The measured spectral irradiance is modeled as the Planck spectral radiance scaled by the Sun's solid angle as seen from Earth:

$$I_\lambda(\lambda) = \Omega_{\odot}\, B_\lambda(\lambda, T)$$

where $\Omega_{\odot}$ is the Sun's solid angle and $B_\lambda$ is the Planck function:

$$B_\lambda(\lambda, T) = \frac{2 h c^2}{\lambda^5} \frac{1}{e^{\frac{h c}{\lambda k_B T}} - 1}$$

Obtainig a value of $T=5536.82\pm 0.02$

### References

[1] GES DISC, “GES DISC Dataset: SORCE SIM Level 3 Solar Spectral Irradiance Daily Means V027 (SOR3SIMD 027),” NASA Goddard Earth Sciences Data and Information Services Center (GES DISC). Disponible en:  (accedido el 15 de marzo de 2026)