# TMTL
This code was developed to calculate Thermal Lens and Thermal mirror intensities, taking into account sample's and the 
air signal. TLLAM contains the model as proposed by Malacarne et al[1].

# TLLAM
The code uses GSL special function, integrations routines and OpenMP in order to evaluate several integrand points in parallel. The code should be compiled like this:
gcc -O3 -static tllam.c -fopenmp -lgsl -lgslcblas -lm -o tllam

#References
1. Malacarne, L. C., Astrath, N. G. C., Lukasievicz, G. V. B., Lenzi, E. K., Baesso, M. L., & Bialkowski, S. E. (2011). Time-Resolved Thermal Lens and Thermal Mirror Spectroscopy with Sample—Fluid Heat Coupling: A Complete Model for Material Characterization. Applied Spectroscopy, 65(1), 99–104. https://doi.org/10.1366/10-06096
