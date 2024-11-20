# Try to reproduce this one
# https://gitlab.cern.ch/acc-models/acc-models-ea/-/blob/en-ea-le/H6/scenarios/positive-120gev-fm-focus/focus.py?ref_type=heads
from cpymad.madx import Madx

mad = Madx()
mad.input('''

call, file = "../acc-models-ea/H6/scenarios/positive-120gev-fm-focus/h6-fm.str";
call, file = "../acc-models-ea/H6/h6fm04.seq";

beam, particle=proton, sequence=H6, PC=120.0,
      ex = 2e-07,
      ey = 5e-08;

use, sequence=H6;

set,  format="15.9f";

twiss, chrom=true, rmatrix=true, save, file="h6_nominal.tfs", betx=10, alfx=0, bety=10, alfy=0, deltap=0.01;
''')

import xtrack as xt
import matplotlib.pyplot as plt
plt.close('all')

line = xt.Line.from_madx_sequence(mad.sequence.h6, deferred_expressions=True)
line.particle_ref = xt.Particles(p0c=120e9, mass0=xt.PROTON_MASS_EV)

sv = line.survey()
sv.plot(projection='ZX')
sv.plot(projection='ZY', axis=None, element_width=2)

tw = line.twiss(betx=10, bety=10, alfx=0, alfy=0)
twplt = tw.plot()
twplt.ylim(left_lo=0, left_hi=2e3, right_lo=-10, right_hi=3)

tw_1mm = line.twiss(betx=10, bety=10, alfx=0, alfy=0, x=1e-3)
tw_1mrad = line.twiss(betx=10, bety=10, alfx=0, alfy=0, px=1e-3)
tw_1percent = line.twiss(betx=10, bety=10, alfx=0, alfy=0, delta=0.01)

# Add to general twiss table
tw['x_1mm'] = tw_1mm.x
tw['x_1mrad'] = tw_1mrad.x
tw['x_1percent'] = tw_1percent.x
tw.plot('x_1mm x_1mrad x_1percent')

# Add marker for focus point
s_focus = 500.

line.discard_tracker()
line.insert_element(element=xt.Marker(), name='focus', at_s=s_focus)

opt = line.match(
    solve=False,
    betx=10, bety=10, alfx=0, alfy=0,
    vary=xt.VaryList(['kq10', 'kq11', 'kq12', 'kq13', 'kq14', 'kq15', 'kq16']),
    targets=[
        xt.TargetSet(at='focus', betx=xt.LessThan(0.5), bety=xt.LessThan(0.5),
                                 alfx=0, alfy=0, dx=0, dy=0),
        xt.TargetSet(at='xced.x0410440', betx=xt.LessThan(750), bety=xt.LessThan(750),
                                         alfx=0, alfy=0),
        xt.TargetSet(at='xcsv.x0410384', betx=xt.LessThan(40), bety=xt.LessThan(40),
                                         alfx=0),
    ]
)


