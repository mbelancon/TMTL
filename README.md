# TMTL
This code was developed to calculate Thermal Lens and Thermal mirror intensities, taking into account sample's and the 
air signal. TLLAM contains the model as proposed by Malacarne et al[1].

# TLLAM
The code uses GSL special function, integrations routines and OpenMP in order to evaluate several integrand points in parallel. 

The code should be compiled like this:
gcc -O3 -static tllam.c -fopenmp -lgsl -lgslcblas -lm -o tllam

and used like this:
./tllam datatofit ds ths thf m w0 v l df cs cf dens denf

where datatofit is the file containing the data to be fitter, which should be in the form:
Time  Intensity

0.	  1.

0.002	1.0121889283900456

0.004	1.0218384966988319

0.006	1.0289487049263586

....  .....

ds is the thermal diffusivity of the sample
ths is the amplitude of the sample signal
thf is the amplitude of the fluid signal
m, w0 and v are geometrical parameters of the experimental setup, as described in the reference [1].
l is the sample's thickness
df is the thermal diffusivity of the fluid
cs is the specific heat of the sample
cf is the specific hear of the fluid
dens is the mass density of the sample
denf is the mass density of the fluid


#References
1. Malacarne, L. C., Astrath, N. G. C., Lukasievicz, G. V. B., Lenzi, E. K., Baesso, M. L., & Bialkowski, S. E. (2011). Time-Resolved Thermal Lens and Thermal Mirror Spectroscopy with Sample—Fluid Heat Coupling: A Complete Model for Material Characterization. Applied Spectroscopy, 65(1), 99–104. https://doi.org/10.1366/10-06096
