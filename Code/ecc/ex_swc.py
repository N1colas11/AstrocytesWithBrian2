#------------------------------------------------------------------------------------------------------#
# import allen institute package for managing .swc files
# documentation for download and use at --> https://alleninstitute.github.io/AllenSDK/install.html
# documentation for swc specific functions --> https://alleninstitute.github.io/AllenSDK/cell_types.html
import allensdk.core.swc as swc
#------------------------------------------------------------------------------------------------------#

import matplotlib.pylab as plt

################################################################################
# Morphology
# Morphological neuron reconstructions are available for download as SWC files. The SWC file format is a white-space delimited text file with a standard set of headers. The file lists a set of 3D neuronal compartments, each of which has:
# Column	Data Type	Description
# id	string	compartment ID
# type	integer	compartment type
# x	float	3D compartment position (x)
# y	float	3D compartment position (y)
# z	float	3D compartment position (z)
# radius	float	compartment radius
# parent	string	parent compartment ID
# Comment lines begin with a '#'. Reconstructions in the Allen Cell Types Database can contain the following compartment types:
################################################################################
file_name = 'ex_astro.swc'
morphology = swc.read_swc(file_name)


# the compartment list has all of the nodes in the file
# print soma compartment data
if True:
    import pprint
    #pprint.pprint(morphology.compartment_list)

    fig, axes = plt.subplots(1, 2, sharey=True, sharex=True)
    # Make a line drawing of x-y and y-z views
    for i,n in enumerate(morphology.compartment_list):
        pprint.pprint(morphology.compartment_list[i])
        morphology.compartment_list[i]['rho'] = .5
        #print("compartment ", i," radius ", morphology.compartment_list[i]['radius'])
        print("compartment ", i," rho ", morphology.compartment_list[i]['rho'])
        for c in morphology.children_of(n):
            axes[0].plot([n['x'], c['x']], [n['y'], c['y']], color='black')
            axes[1].plot([n['z'], c['z']], [n['y'], c['y']], color='black')

    axes[0].set_ylabel('y')
    axes[0].set_xlabel('x')
    axes[1].set_xlabel('z')
    plt.show()
    plt.savefig("morpho_ex.pdf")





