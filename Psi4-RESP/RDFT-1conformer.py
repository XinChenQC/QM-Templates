import psi4
import resp

# Initialize three different conformations of ethanol
geometry = """
0   1
    C   -0.2753000000       -2.2228000000        0.0000000000     
    C    1.0407000000       -1.7594000000        0.0000000000     
    C   -1.0775000000        0.0564000000        0.0000000000     
    C   -1.3343000000       -1.3149000000        0.0000000000     
    H   -0.4776000000       -3.3036000000        0.0000000000     
    H    1.8753000000       -2.4754000000        0.0000000000     
    H   -1.9127000000        0.7719000000        0.0000000000     
    H   -2.3715000000       -1.6800000000        0.0000000000     
    C    0.5221000000        2.0333000000        0.0000000000     
    H    1.5322000000        2.3862000000        0.0000000000     
    C    1.2975000000       -0.3886000000        0.0000000000     
    H    2.3254000000       -0.0916000000        0.0000000000     
    C    0.2381000000        0.5197000000        0.0000000000     
    O   -0.4298000000        2.8564000000        0.0000000000     
"""

psi4.set_options({'reference': 'rhf',
                  'scf_type' : 'direct'})

mol1 = psi4.geometry(geometry)
mol1.update_geometry()
mol1.set_name('conformer1')
molecules = [mol1]

# Specify options
options = {'VDW_SCALE_FACTORS' : [1.4, 1.6, 1.8, 2.0],
           'VDW_POINT_DENSITY'  : 1.0,
           'RESP_A'             : 0.0005,
           'RESP_B'             : 0.1,
           'RESTRAINT'          : True,
           'IHFREE'             : False,
           'CONSTRAINT_GROUP'   : [[4,2],[8,6],[11,3],[12,7]],
           "METHOD_ESP"         : "b3lyp",
           "BASIS_ESP"          : "6-31+g*"
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
