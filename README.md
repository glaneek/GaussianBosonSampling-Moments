# GaussianBosonSampling-Moments

GBS-M is a collection of scripts that demonstrates how a GBS system can be characterised through moments, procided GBS experimental data.
It makes use of the strawberryfield python package: https://strawberryfields.readthedocs.io/en/stable/

## Theory

### Modeling the interferometer

A linear optical tranformation can be mathematically represented by a Symplectic matrix. For example a unitary tranformation 

![equation](https://latex.codecogs.com/gif.latex?U%3D%20%5Cbegin%7Bpmatrix%7D%20A%20%26%20B%20%5C%5C%20B%5E*%20%26%20A%5E*%20%5Cend%7Bpmatrix%7D)

can be represented by the corresponding symplectic matrix

![equation](https://latex.codecogs.com/gif.latex?S%3D%5Cfrac%7B1%7D%7B2%7D%20%5Cbegin%7Bpmatrix%7D%20A&plus;B&plus;A%5E*&plus;B%5E*%20%26%20-i%28A-A%5E*&plus;B%5E*-B%29%5C%5C%20i%28A-B%5E*&plus;B-A%5E*%29%20%26%20A&plus;A%5E*-B-B%5E*%20%5Cend%7Bpmatrix%7D)

A tranformation from a combination of linear optical components is given by a direct sum of the corresponding sympecitc matrices. For example the following tranformation 

![alt text](https://github.com/glaneek/GaussianBosonSampling-Moments/blob/master/tranformationExample.PNG?raw=true)

can be described as a single 4x4 symplectic matrix

![equation](https://latex.codecogs.com/gif.latex?S_%7BBS%7D%5Ccdot%28S_%7B%5Czeta%7D%5Coplus%20%5Cmathbb%7BI%7D%29%3D%20S_%7BBS%7D%5Ccdot%20%5Cbegin%7Bpmatrix%7D%20%5Ccosh%20r%20&plus;%20%5Ccos%5Cphi%5Csinh%20r%20%26%20%5Csin%5Cphi%5Csinh%20r%20%26%200%20%26%200%20%5C%5C%20-%5Csin%5Cphi%5Csinh%20r%20%26%20%5Ccosh%20r%20-%20%5Ccos%5Cphi%5Csinh%20r%20%26%200%20%26%200%5C%5C%200%20%26%200%20%26%201%20%26%200%5C%5C%200%20%26%200%20%26%200%20%26%201%20%5Cend%7Bpmatrix%7D)

In this manner the interferometer of the GBS can be represented as a single matrix S.

### Moments

The first and second order moments along the quatratures are then tranformed by the interferometer as follows

* Frist order: ![equation](https://latex.codecogs.com/gif.latex?%5CVec%7Bd%7D%27%3DS%5CVec%7Bd%7D)
* Second order: ![equation](https://latex.codecogs.com/gif.latex?V%27%3DSVS%5ET)

This leads to first and second order moments of the photon numbers in terms of the second order, tranformed, quadrature covariance matrix (V') such that

* First order:  ![equation](https://latex.codecogs.com/gif.latex?%5Clangle%20%5Chat%7BN%7D_n%20%5Crangle%20%3D%5Cfrac%7B1%7D%7B2%7D%28V%27_%7Bnn%7D&plus;V%27_%7Bn&plus;1n&plus;1%7D-1%29)

* Second order:  ![equation](https://latex.codecogs.com/gif.latex?%5Clangle%20N_n%20N_m%20%5Crangle%3D%20%5Cleft%5B%20%5Csum%5EW_%7Bnnmm%7D&plus;%20%5Csum%5EW_%7Bnn%28m&plus;1%29%28m&plus;1%29%7D&plus;%20%5Csum%5EW_%7B%28n&plus;1%29%28n&plus;1%29mm%7D&plus;%20%5Csum%5EW_%7B%28n&plus;1%29%28n&plus;1%29%28m&plus;1%29%28m&plus;1%29%7D-%20%5Csum_%7Bnm%7D%5ED&plus;1%20%5Cright%5D%20%5Cleft%28%5Cfrac%7B%282%5Cpi%29%5E%7B2n%7D%7D%7B%5Ctext%7Bdet%7DV%5ET%7D%5Cright%29%5E%7B%5Cfrac%7B1%7D%7B2%7D%7D%5C%5C)

where

![equation](https://latex.codecogs.com/gif.latex?%5Csum_%7Bnm%7D%5ED%3DV%27_%7Bnn%7D&plus;V%27_%7Bn&plus;1n&plus;1%7D&plus;V%27_%7Bmm%7D&plus;V%27_%7Bm&plus;1m&plus;1%7D) 

This allows the analytical construction of photon number covariance matrix for a well defined interferometer

![equation](https://latex.codecogs.com/gif.latex?%5Clangle%20%5Chat%7BN%7D_i%20%5Chat%7BN%7D_j%5Crangle-%5Clangle%20%5Chat%7BN%7D_i%20%5Crangle%20%5Clangle%20%5Chat%7BN%7D_j%20%5Crangle%3D%5Cfrac%7B1%7D%7B2%7D%28V%5E2_%7Bi%2Cj%7D&plus;V%5E2_%7Bi&plus;1%2Cj&plus;1%7D&plus;V%5E2_%7Bi%2Cj&plus;1%7D&plus;V%5E2_%7Bi&plus;1%2Cj%7D%29)

which can be compared to experimental moments calculated from the GBS data using the following statistics formula

![equation](https://latex.codecogs.com/gif.latex?%5Cmathrm%7BE%7D%28g%28x_1%2C...%2Cx_n%29%29%3D%5Csum_%7Bx_1%7D%5Ccdot%5Ccdot%5Ccdot%5Csum_%7Bx_n%7Dg%28x_1%2C...%2Cx_n%29%5Cmathrm%7BP%7D%28x_1%3DX_1%2C...%2Cx_n%3DX_n%29)

### Comparison of experimental and analytical moments

To compare the two covariance martices, the Residual Sum of Squares is considered such that 

![equation](https://latex.codecogs.com/gif.latex?RSS%3D%5Csum_%7B%5Csigma%5Cin%20cov%5BN_i%2CN_j%5D%7D%5Cleft%5B%28%5CDelta%5CTilde%7BN%7D_%7B%5Csigma%7D%29%5E2-%28%5CDelta%20N_%7B%5Csigma%7D%29%5E2%29%5Cright%5D%5E2)

where

* ![equation](https://latex.codecogs.com/gif.latex?%5CDelta%20N_%7B%5Csigma%7D) is the expimental covariance matrix
* ![equation](https://latex.codecogs.com/gif.latex?%5CDelta%5CTilde%7BN%7D_%7B%5Csigma%7D) is the analytical covariance matrix.

### System Characterisation

The closer the analytical system is characterised to the atual experimental set up, the lower the RSS value. i.e the RSS value can be thought of a function to be minimized with GBS data, interferometers and mode squeezing parameters as arguments.

Therefore, starting from ideal optical component parameters, we try to minimize the RSS value by twitching their values. The optimal parameter values are then the true values of the experimental set up. 

## Scripts

* GBS_functions.py: This script holds every required funcion for the rest of the scripts. Make sure you 
		  import everything. Every function has its own discription in line. In general the functions are
		  it separated into four labelled categories: 
	* Strawberryfield: functions to model a system and get probabilities
	* Symplectic Matricesl: functions to model a system and get photon number covariances through the symplectic matrices formalism
	* Sampling Simulation : it samples from the strawberryfield exact PDF to simulate experimental data, of finite sample size
	* General: functions that extract data from files.


## Run an optimisation

* GBS_model_generator.py: Define your model here! beamsplitters, phaseshifters, squeezers. It models the system
 			with the strawberry field package and saves the probabilities for corresponding states
			in a named file. Errori identification: Manually introduce an error value to one of 
			the optical components. 

* GBS_sampling.py: Gets the exact probabilities-file that was generated from GBS_model_generator.py, and samples
		 to simulate experimental data.

* GBS_optimizer.py: The optimizer will sample, calculate moments from experimental data, and compare values
		  to the symplectic theoretical ones, whith a guessed error. It will repeat the process
	          with a different error value. It slowly converges to the true error value.

## Example Run

The default system that is simulated by the scripts is a 4 mode interferoemter and 8 beam splitters:

![alt text](https://github.com/glaneek/GaussianBosonSampling-Moments/blob/master/system.png?raw=true)

By considering statistics of outputs states between |0,0,0,0> and |12,12,12,12> the BS values where identified with +/-0.5% which is almost an order of o magnitude compared to the typical provided manufacturing error bound.sy

![alt text](https://github.com/glaneek/GaussianBosonSampling-Moments/blob/master/characterising8BS.png?raw=true)



