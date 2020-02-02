import numpy as np

def printData(astro_mon, b_mon, syn1_mon, syn2_mon):
	print("C= ", astro_mon.C)
	print("V= ", b_mon.V)  
	#print("CC0= ", b_mon.CC0)
	#print("coupling_C= ", astro_mon.coupling_C)
	#print("coupling_electro= ", astro_mon.coupling_electro)
	#print("electro_diffusion= ", astro_mon.electro_diffusion)
	#print("VVT=-2*V/V_T ", b_mon.VVT)
	#print("V0= ", astro_mon.V0)
	#print("V_T= ", V_T)
	print("A000_nb_connections: ", astro_mon.A000_nb_connections)
	print("A1= ", astro_mon.A1)
	print("B1= ", astro_mon.B1)
	print("C1= ", astro_mon.C1)
	print("CC0= ", astro_mon.CC0)
	print("V0= ", astro_mon.V0)
	#print("s= ", b_mon.s)
	print("dR2= ", b_mon.dR2)
	print("L= ", b_mon.L)
	print("syn2::CC= ", syn2_mon.CC)
	print("syn2::Q= ", syn2_mon.Q)
	print("syn1::P= ", syn1_mon.P)


#----------------------------------------------------------------------
def plots(astro_mon):
	################################################################################
	# Analysis and plotting
	################################################################################
	fig, (ax1,ax2) = plt.subplots(nrows=1, ncols=2, figsize=(6.26894 * 2, 6.26894 * 0.66),
                       	gridspec_kw={'left': 0.1, 'bottom': 0.12})
	
	ax1.plot(astro_mon.t/second,astro_mon.C[0]/umolar,label='C1')
	ax1.plot(astro_mon.t/second,astro_mon.C[1]/umolar,label='C2')
	ax1.legend()
	ax1.set(xlabel='time (ms)',ylabel='Ca concentration (umolar)')
	
	plt.savefig("plot.pdf")
	#plt.show()
#----------------------------------------------------------------------
def connectionMatrices(astrocytes1, synapses1, synapses2):
	matrix = np.zeros([len(astrocytes1), len(astrocytes1)], dtype=bool)
	matrix[synapses1.i[:], synapses1.j[:]] = True
	print("Connection matrix 1:")
	print(matrix)

	# Matrix2 is a permutation matrix. Is this significant?
	matrix2 = np.zeros([len(astrocytes1), len(astrocytes1)], dtype=bool)
	matrix2[synapses2.i[:], synapses2.j[:]] = True
	print("Connection matrix 2:")
	print(matrix2)
#-------------------------------------------
def lowLevelInfo(synapses1):
	for k,v in synapses1.variables.items():
		#print("k: ", k)
		#print("v: ", v)
		#print("..dir(k)= ", dir(k))
		#print("..dir(v)= ", dir(k))
		pass

	#print("..dir(synapses1)= ", dir(synapses1))
	#print("..dir(synapses1.subexpression_updater)= ", dir(synapses1.subexpression_updater))
	#print("..dir(synapses1.summed_updaters)= ", dir(synapses1.summed_updaters))
	#print()
#----------------------------------------------------------------------
