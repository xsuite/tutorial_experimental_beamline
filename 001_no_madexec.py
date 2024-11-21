# Try to reproduce this one
# https://gitlab.cern.ch/acc-models/acc-models-ea/-/blob/en-ea-le/H6/scenarios/positive-120gev-fm-focus/focus.py?ref_type=heads

import xtrack as xt
import matplotlib.pyplot as plt
plt.close('all')

# Load lattice and strengths
env = xt.load_madx_lattice('h6fm04.seq')
env.vars.load_madx('h6-fm.str')
line = env['h6']

# Associate a reference particle
line.particle_ref = xt.Particles(p0c=120e9, mass0=xt.PROTON_MASS_EV)

# Survey
sv = line.survey()
sv.plot(projection='ZX')
sv.plot(projection='ZY', axis=None, element_width=2)

# Inspect the lattice
tt = line.get_table(attr=True)
tt_bend = tt.rows[tt.element_type=='Bend']
tt_quad = tt.rows[tt.element_type=='Quadrupole']

tt_bend.cols['s', 'length', 'angle_rad', 'k0l', 'rot_s_rad'].show()
tt_quad.cols['s', 'length', 'k1l', 'rot_s_rad'].show()


tw = line.twiss(betx=10, bety=10, alfx=0, alfy=0)