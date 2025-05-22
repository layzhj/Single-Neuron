from neuron import h, load_mechanisms
import os
import matplotlib.pyplot as plt
import pandas as pd
from cell.PyramidalCell import PyramidalCell
from litus.bls import BilayerSonophore

h.load_file('stdrun.hoc')
load_mechanisms('./mechanisms/')

trans_cm = 1e2  # from f/m2 to uf/cm2
trans_m = 1e6   # from m to um
trans_qm = 1e-5  # from uf/cm2*mv to f/m2*v

cell = PyramidalCell(0, k_cal=5/3, k_kca=50)

Qm = cell.soma(0.5).v * cell.soma(0.5).cm * trans_qm
bls = BilayerSonophore(a=32e-9, cm0= 1/trans_cm, Qm0=Qm)
print(bls.Delta * trans_m)

stim = h.DcDt(cell.soma(0.5))
stim._ref_c = cell.soma(0.5)._ref_cm
stim.A = 500e3
stim.f = 100
stim.Delta = bls.Delta * trans_m
stim.LJ_C = bls.LJ_approx['C']
stim.LJ_alpha = bls.LJ_approx['x0'] * trans_m
stim.m = bls.LJ_approx['nrep']
stim.n = bls.LJ_approx['nattr']
stim.tbegin = 0
stim.tdur = 300
stim.PRF = 1
stim.DC = 0.05

# FOR FIRNG
iclamp = h.IClamp(cell.soma(0.5))
iclamp.delay = 100
iclamp.dur = 200
iclamp.amp = 0.3

t_vec = h.Vector().record(h._ref_t, 0.1)
v_vec = h.Vector().record(cell.soma(0.5)._ref_v, 0.1)
c_vec = h.Vector().record(cell.soma(0.5)._ref_c, 0.1)

h.cvode_active(1)
cv = h.CVode()
cv.atolscale("DcDt.ng", 1e-22)
cv.atolscale("DcDt.U", 1)
cv.atolscale("DcDt.Z", 1e-6)

cv.solve(400)
# h.tstop = 400
# h.run()

save_path = './sonic'
os.makedirs(save_path, exist_ok=True)

df = pd.DataFrame({'t': t_vec, 'v': v_vec})
df.to_csv(save_path+'/pyramidal.csv')

print(t_vec.as_numpy())
print(v_vec.as_numpy())

fig, ax = plt.subplots(2, 1)

ax[0].plot(t_vec, v_vec, color='orange')
ax[0].set_ylabel('Voltage (mV)')

ax[1].plot(t_vec, c_vec, color='orange')
ax[1].set_ylabel('Capacitance ($uF/cm^{2}$)')
ax[1].set_xlabel('Time (ms)')

plt.plot(t_vec, v_vec)
plt.show()
