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

twiss, chrom=true, rmatrix=true, betx=10, alfx=0, bety=10, alfy=0;
''')
tw_mad = xt.Table(mad.table.twiss)

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

xo.assert_allclose(tw_check.betx, tw_mad_check.betx, rtol=2e-5, atol=0)
xo.assert_allclose(tw_check.bety, tw_mad_check.bety, rtol=2e-5, atol=0)
xo.assert_allclose(tw_check.dx, tw_mad_check.dx, rtol=2e-5, atol=1e-4)
xo.assert_allclose(tw_check.dy, tw_mad_check.dy, rtol=2e-5, atol=1e-4)

trm = tw.get_R_matrix_table()
trm_check = trm.rows[check_at]

for ii in range(6):
    for jj in range(6):
        rterm_mad = tw_mad_check[f're{ii+1}{jj+1}']
        rterm_xs = trm_check[f'r{ii+1}{jj+1}']
        atol=1e-4*np.max(np.abs(rterm_mad))
        if atol<1e-14:
            atol=1e-14
        xo.assert_allclose(rterm_xs, rterm_mad, rtol=2e-5, atol=atol)

import matplotlib.pyplot as plt
plt.close('all')
plt.figure(1)
plt.plot(tw_mad.s, tw_mad.betx, label='mad')
plt.plot(tw.s, tw.betx, label='env')

plt.figure(2)
plt.plot(tw_mad.s, tw_mad.re11, label='mad')
plt.plot(trm.s, trm.r11, label='env')

plt.figure(3)
plt.plot(tw_mad.s, tw_mad.re12, label='mad')
plt.plot(trm.s, trm.r12, label='env')

plt.figure(4)
plt.plot(tw_mad.s, tw_mad.re21, label='mad')
plt.plot(trm.s, trm.r21, label='env')

plt.figure(5)
plt.plot(tw_mad.s, tw_mad.re22, label='mad')
plt.plot(trm.s, trm.r22, label='env')



plt.show()