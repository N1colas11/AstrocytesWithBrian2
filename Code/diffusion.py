# My code still does not work

# 2019-10-19,9.32am, working code. 

from brian2 import *

set_device('cpp_standalone', directory="cpp")  # Use fast "C++ standalone mode"


# link two compartments through diffusion. Single synapse. 
# Use a Poisson input on one of the branches. 

# f(t) should be Poisson

# du/dt = D (u-v) + f(t)
# dv/dt = D (v-u)

# To do this with Brian2, create a Neuron with equations: 

#Parameters: 
D = 1*Hz;
#freq = 50*Hz
sim_dt = 1.e-2*second
defaultclock.dt = sim_dt     # Set the integration time
duration = 5*second
source = 1*Hz
w = 0.2 # Amount added by each incoming spike. Simulate Poisson Input

comp_eqs = """
du/dt = diff : 1
diff : Hz
freq : Hz
"""

N_comp = 2
# When exceeded, the threshold creates an input spike
compartments = NeuronGroup(N_comp, comp_eqs, threshold='rand()< freq*dt', method='heun')
compartments.u = [1.,0.]  # initial condition for both Neurons
# could use a parameter instead of a state variable for freq if constant
compartments.freq = np.asarray([0.2,0.2],dtype=float)*Hz  # asarray not required

synapse_eqs = """
    delta_u = u_post-u_pre : 1
    diff_post = -D * delta_u : Hz (summed) 
"""

# My code is not producing the correct results. MDP's is. 
comp_to_comp = Synapses(compartments, compartments, 
                        model=synapse_eqs,
                        on_pre='u_pre += w',
                        on_post='u_post += w')

# Connect i (pre) to j (post). Not correct. We must connect in both direction. 
#comp_to_comp.connect(i=0, j=1)

# Connects i (pre) to j (post) AND i (post) to j (pre)
comp_to_comp.connect()

# You want to record. It is possible to record at lower interval frequency than that imposed by dt
monitor = StateMonitor(compartments,'u',record=True)

run(duration, report='text')
print("*** End of simulation ***")
#----------------------------------------------------------------------

print(monitor.t_)
plt.plot(monitor.t_,monitor.u[0],label='C1')
plt.plot(monitor.t_,monitor.u[1],label='C2')
plt.legend()
plt.show()
#----------------------------------------------------------------------
