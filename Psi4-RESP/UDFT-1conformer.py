import psi4
import resp

# Initialize three different conformations of ethanol
geometry = """
0   2
C    0.00000000  0.00000000  0.00000000
C    1.48805540 -0.00728176  0.39653260
O    2.04971655  1.37648153  0.25604810
H    1.58679428 -0.33618761  1.43102358
H    2.03441010 -0.68906454 -0.25521028
F   -0.40814044 -1.00553466  0.10208540
F   -0.54635470  0.68178278  0.65174288
F   -0.09873888  0.32890585 -1.03449097
"""
mol1 = psi4.geometry(geometry)
mol1.update_geometry()
mol1.set_name('conformer1')

molecules = [mol1]

# Specify options

psi4.set_options({'reference': 'uks',
                  'scf_type' : 'direct'})
options = {'VDW_SCALE_FACTORS' : [1.4, 1.6, 1.8, 2.0],
           'VDW_POINT_DENSITY'  : 1.0,
           'RESP_A'             : 0.0005,
           'RESP_B'             : 0.1,
           'RESTRAINT'          : True,
           'IHFREE'             : False,
           'CONSTRAINT_GROUP'   : [[7,8,9]],  # constrain Fluorine atom
           "METHOD_ESP"         : "b3lyp",
           "BASIS_ESP"          : "CC-PVTZ"
           }

# Call for first stage fit
charges1 = resp.resp(molecules, options)

print("Restrained Electrostatic Potential Charges")
print(charges1[1])

options['RESP_A'] = 0.001
resp.set_stage2_constraint(molecules[0], charges1[1], options)

options['grid'] = []
options['esp'] = []
# Add constraint for atoms fixed in second stage fit
for mol in range(len(molecules)):
    options['grid'].append('%i_%s_grid.dat' %(mol+1, molecules[mol].name()))
    options['esp'].append('%i_%s_grid_esp.dat' %(mol+1, molecules[mol].name()))

# Call for second stage fit
charges2 = resp.resp(molecules, options)
print("\nStage Two\n")
print("RESP Charges")
print(charges2[1])
