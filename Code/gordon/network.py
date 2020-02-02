from brian2 import *

duration = 1*ms          # Total simulation time
sim_dt = 0.1*ms            # Integrator/sampling step

# Model definition
defaultclock.dt = sim_dt     # Set the integration time

astro_eqs1 = '''
   dC/dt = 1*Hz : 1
   A : 1
'''

astro_eqs2 = '''
   dD/dt = 1*Hz : 1
   B : 1
'''

N_astro = 4 # Total number of astrocytes in the network
astrocytes1 = NeuronGroup(N_astro, astro_eqs1, method='euler', name="ng1", order=2)
astrocytes2 = NeuronGroup(N_astro, astro_eqs2, method='euler', name="ng2", order=4)
astrocytes1.C = 1
astrocytes2.D = 2

# Initital Conditions

# How to make A_pre update before the state updater?
synapse_eqs1 = '''
   A_pre = C_pre : 1 (summed)  
   dE/dt = 1.0*Hz * A_pre : 1
'''

synapse_eqs2 = '''
   B_pre = D_pre : 1 (summed)
   dF/dt = 1.0*Hz : 1
'''

synapse1 = Synapses(astrocytes1, astrocytes2, model=synapse_eqs1, method='euler', name="sg1", order=1)
synapse2 = Synapses(astrocytes2, astrocytes1, model=synapse_eqs2, method='euler', name="sg2", order=3)
synapse1.connect()
synapse2.connect()


# None of the expressions below is changing "when" in scheduling_summary() . Why? 
#astrocytes2.when = 'groups'
#astrocytes1.when = 'start'
#synapse2.when = 'start'
#synapse1.when = 'end'

print('==> ', dir(astrocytes1))
print('==> ', dir(synapse1))
print('==> astrocytes1 variables: ', astrocytes1.variables)
print('==> synapse1 variables: ', dir(synapse1.variables))
for k,v in synapse1.variables.items():
	print("k,v: ", k,v)
print('==> synapse1 variables: ', synapse1.variables.items())
#quit()

astro1_mon = StateMonitor(astrocytes1, variables=['C'],record=True, name='ngm1')
astro2_mon = StateMonitor(astrocytes2, variables=['D'],record=True, name='ngm2')
syn1_mon   = StateMonitor(synapse1, variables=['E'],record=True)
syn2_mon   = StateMonitor(synapse2, variables=['F'],record=True)

#print(scheduling_summary())
#print("==> ", Network.schedule)
#run(duration, report='text')
#quit()
#----------------------------------------------------------------------

# sg1, ng1, sg2, ng2
#print("==> astrocytes1: ", dir(astrocytes1))
#astrocytes1.when = 'synapses'
#astrocytes2.when = 'synapses'
#astrocytes1.state_updater.when = 'synapses'
#astrocytes1.state_updater.order = -2
#astrocytes2.state_updater.order = -1
#print('==> ', dir(synapse1))
#print('==> summedvariableupdater', dir(synapse1.summed_updaters))
#print('==> ', dir(synapse1.subexpression_updater))
#print()
#print('==> ', synapse1.equations)
#synapse1.when = 'groups'
#synapse1.order = 0
#synapse2.when = 'start'
#synapse2.order = 1



net = Network(collect())
#net.add(synapse1)
#net.add(synapse2)
#net.add(astrocytes2)
#net.add(astrocytes1)

print('\n==> net: ', dir(net), '\n')
print('\n==> net.get_states: ', net.get_states(), '\n')
print('\n==> net.objects: ', net.objects, '\n')
print('\n==> Synapse: ', dir(synapse1), '\n')
print()
for o in net.objects: 
	print("object: ", o)
print(net.objects)
for s in net.get_states():
    print("state: ", s)
print()

print("when ng1= ", astrocytes1.when)
print("when ng2= ", astrocytes2.when)
print("when sg1= ", synapse1.when)
print("when sg2= ", synapse2.when)

print(net.scheduling_summary())
print(net.schedule)
#print(net.get_states())

net.run(duration, report='text')

