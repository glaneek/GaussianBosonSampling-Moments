# GaussianBosonSampling-Moments

GBS-M is a collection of scripts that demonstrates how a GBS system can be characterised through moments.
It makes use of the strawberryfield python package: https://strawberryfields.readthedocs.io/en/stable/  

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
