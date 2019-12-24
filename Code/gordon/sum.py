from brian2 import *
codegen.target = 'cython'

### General parameters
duration = 50*ms       # Total simulation time
sim_dt = 25*ms            # Integrator/sampling step

# Model definition
defaultclock.dt = sim_dt     # Set the integration time

astro_eqs = '''
C : 1
g : 1
ddd : 1
eee : 1
coupling : 1
nb_connections : 1
'''

N_astro = 3 # Total number of astrocytes in the network
astrocytes = NeuronGroup(N_astro, astro_eqs)
astrocytes.g = [10,20,40]
astro_mon = StateMonitor(astrocytes, variables=['C','ddd','eee', 'coupling', 'nb_connections'], record=True)

# Diffusion between astrocytes
astro_to_astro_eqs = '''
# Bad idea to have pre on both sides of an equation
# Summation over all connections (including self-connections)
#eee_pre = g_pre : 1 (summed)  # wrong answer
eee_pre = g_post : 1 (summed)  # correct answer (self-connecton is important)
ddd_post = g_pre : 1 (summed)  # correct answer (self-connecton is important)

# How do I use the number of connections in the Neuron without using (summed)? 
# Does not seem possible. Perhaps it is never required?
nb_connections_pre = N_pre : 1 (summed)

#sum = g_pre : 1  (summed) # Not allowed. 

# Correct computation of diffusion terms. Same result with and without self-connections
# Correct computation of diffusion terms. Same result with and without self-connections
coupling_pre = g_post - g_pre : 1 (summed)
'''

astro_to_astro = Synapses(astrocytes, astrocytes, model=astro_to_astro_eqs)
#syn_mon = StateMonitor(astro_to_astro, variables=['fff'], record=True)
astro_to_astro.connect()

# Simulation run
run(duration, report='text')

print("C= ", astro_mon.C)
print("ddd= ", astro_mon.ddd)
print("eee= ", astro_mon.eee)
print("coupling= ", astro_mon.coupling)
print("nb_connections= ", astro_mon.nb_connections)
#print("fff= ", syn_mon.fff)
