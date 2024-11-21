# Try to reproduce this one
# https://gitlab.cern.ch/acc-models/acc-models-ea/-/blob/en-ea-le/H6/scenarios/positive-120gev-fm-focus/focus.py?ref_type=heads

import numpy as np
import xtrack as xt
import matplotlib.pyplot as plt
plt.close('all')

# Load lattice and strengths
env = xt.load_madx_lattice('h6fm04.seq')
env.vars.load_madx('h6-fm.str')
line = env['h6']
line.replace_all_repeated_elements()

# Associate a reference particle
line.particle_ref = xt.Particles(p0c=120e9, mass0=xt.PROTON_MASS_EV)

# Survey
sv = line.survey()
sv.plot(projection='ZX')
sv.plot(projection='ZY', axis=None, element_width=2)

# Inspect the lattice using tables
tt = line.get_table(attr=True)
tt_bend = tt.rows[tt.element_type=='Bend']
tt_quad = tt.rows[tt.element_type=='Quadrupole']

tt_bend.cols['s', 'length', 'angle_rad', 'k0l', 'rot_s_rad'].show()
tt_quad.cols['s', 'length', 'k1l', 'rot_s_rad'].show()

# Inspect one element
line.info('mbxhc.x0410124')

# Inspect one variable
line.info('kb4', limit=None)

# Twiss
tw = line.twiss(betx=10, bety=10, alfx=0, alfy=0)
twplt = tw.plot()
twplt.ylim(left_lo=0, left_hi=2e3, right_lo=-10, right_hi=3)


# Add markers for focus points
s_focus = np.linspace(490., 540., 9)
line.discard_tracker()
for ii, ss in enumerate(s_focus):
    line.insert_element(element=xt.Marker(), name=f'ff_{ii+1}', at_s=ss)

# Match waist at focus points
opt = line.match(
    solve=False,
    betx=10, bety=10, alfx=0, alfy=0,
    vary=xt.VaryList(['kq10', 'kq11', 'kq12', 'kq13', 'kq14', 'kq15', 'kq16'], step=1e-3),
    targets=[
        xt.TargetSet(at='ff_5', betx=xt.LessThan(15), bety=xt.LessThan(12),
                                alfx=0, alfy=0, dx=0, dy=0),
        xt.TargetSet(at='xced.x0410440', betx=xt.LessThan(750), bety=xt.LessThan(750),
                                         alfx=0, alfy=0),
        xt.TargetSet(at='xcsv.x0410384', betx=xt.LessThan(40), bety=xt.LessThan(40),
                                         alfx=0),
    ]
)
opt.step(20)
opt.targets['ff_5_dx'].weight = 10000
opt.targets['ff_5_dy'].weight = 10000
opt.step(20)

tw = line.twiss(betx=10, bety=10, alfx=0, alfy=0)
twplt = tw.plot()
twplt.ylim(left_lo=0, left_hi=2e3, right_lo=-10, right_hi=3)
plt.axvline(x=tw['s', 'ff_5'], color='k', linestyle='--')

bsz = tw.get_beam_covariance(gemitt_x=2e-7, gemitt_y=5e-8)
bsz.rows['ff_5'].cols['s', 'sigma_x', 'sigma_y'].show()

plt.figure(10)
ax = plt.gca()
bsplt = tw.plot(lattice_only=True, ax=ax)
ax = bsplt.left
ax.plot(bsz.s, 1e3*bsz.sigma_x, label=r'$\sigma_x$')
ax.plot(bsz.s, 1e3*bsz.sigma_y, label=r'$\sigma_y$')
ax.set_ylabel(r'$\sigma$ [mm]')
ax.legend(loc='upper right')
bsplt.lattice.set_ylim(-5, 1.2)
ax.set_ylim(0, 15)
plt.axvline(x=tw['s', 'ff_5'], color='k', linestyle='--')

plt.show()
