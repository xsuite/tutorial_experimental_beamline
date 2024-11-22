from cpymad.madx import Madx
import numpy as np
import xtrack as xt
import xobjects as xo

mad = Madx()
mad.input('''

call, file = "../acc-models-ea/H6/h6fm04.seq";
call, file = "../acc-models-ea/H6/scenarios/positive-120gev-fm-focus/h6-fm.str";


beam, particle=proton, sequence=H6, PC=120.0,
      ex = 2e-07,
      ey = 5e-08;

use, sequence=H6;

set,  format="15.9f";

twiss, chrom=true, ripken=true, rmatrix=true, betx=10, alfx=0, bety=10, alfy=0, deltap=0.01;
''')
tw_mad = xt.Table(mad.table.twiss)
lmad = xt.Line.from_madx_sequence(mad.sequence.h6, deferred_expressions=True)
lmad.particle_ref = xt.Particles(p0c=120e9, mass0=xt.PROTON_MASS_EV)

tt = lmad.get_table(attr=True)
tt_bend = tt.rows[tt.element_type=='Bend']

lmad.discard_tracker()
for nn, ll in zip(tt_bend.name, tt_bend.length):
    lmad.element_dict[nn] = xt.Drift(length=ll)
    # lmad[nn].k0 = 0
    # lmad[nn].h = 0

tw1 = lmad.twiss(betx=10, bety=10)

ss = lmad.to_madx_sequence('h6')
with open('sss.seq', 'w') as fid:
    fid.write(ss)
mad2 = Madx()
mad2.call('sss.seq')
mad2.beam()
mad2.input('''
beam, particle=proton, sequence=H6, PC=120.0,
      ex = 2e-07,
      ey = 5e-08;

use, sequence=H6;

set,  format="15.9f";

twiss, chrom=true, ripken=true, rmatrix=true, betx=10, alfx=0, bety=10, alfy=0;
''')
tw_mad2 = xt.Table(mad2.table.twiss)
line2 = xt.Line.from_madx_sequence(mad2.sequence.h6, deferred_expressions=True)
line2.particle_ref = xt.Particles(p0c=120e9, mass0=xt.PROTON_MASS_EV)
tw2 = line2.twiss(betx=10, bety=10)

import matplotlib.pyplot as plt
plt.close('all')
plt.figure(1)
plt.plot(tw1.s, tw1.betx, label='l1', alpha=0.5)
plt.plot(tw2.s, tw2.betx, label='l2', alpha=0.5)
plt.plot(tw_mad2.s, tw_mad2.betx, label='mad2', alpha=0.5)

plt.show()



prrrr

env = xt.load_madx_lattice('h6fm04.seq')
env.vars.load_madx('h6-fm.str')
line = env['h6']
line.particle_ref = xt.Particles(p0c=120e9, mass0=xt.PROTON_MASS_EV)
tt = line.get_table(attr=True)
line.configure_bend_model(edge='full')

tw = line.twiss(betx=10, bety=10)

check_at = [
 't4..centre',
 'vxsv.x0410104',
 'begin.vac',
 'xwca.x0410404',
 'xsci.x0410475',
 'xemc.x0410476',
 'xdwc.x0410488',
 'h6a',
 'h6b',
 'h6c']

tw_check = tw.rows[check_at]
tw_mad_check = tw_mad.rows[[nn+':1' for nn in check_at]]

# xo.assert_allclose(tw_check.betx, tw_mad_check.betx, rtol=1e-5, atol=0)

import matplotlib.pyplot as plt
plt.close('all')

plt.plot(tw_mad.s, tw_mad.betx, label='mad')
plt.plot(tw.s, tw.betx, label='env')
plt.plot(t1.s, t1.betx, label='cpymad')
plt.show()