
def printData(astro_mon, b_mon, syn2_mon):
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

