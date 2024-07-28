import os
from obspy import read
import matplotlib.pyplot as plt
import numpy as np
import subprocess

raw_path = './IU.CCM.BHE'
vel_path = './vel'
none_path = './none'
os.system('rm -r vel none')

# remove response to vel
s_vel = "r " + raw_path + "\n"
s_vel += "rmean; rtr; taper\n"
# ######## Notice the Resp name and the filter band  ########
s_vel += "trans from polezero subtype " + '../pz/IU.CCM.BHE.PZ' +" to vel freq 0.004 0.005 4 4.5 \n"
# unit to cm/s
s_vel += "mul 1.0e2\n"
# ######## Notice the Resp name and the filter band  ########
s_vel += "w " + vel_path + "\n"
s_vel += "q \n"
subprocess.Popen(['sac'], stdin=subprocess.PIPE).communicate(s_vel.encode())

# remove response to none
s_none = "r " + raw_path + "\n"
s_none += "rmean; rtr; taper\n"
# ######## Notice the Resp name and the filter band  ########
s_none += "trans from polezero subtype " + '../pz/IU.CCM.BHE.PZ' +" to none freq 0.004 0.005 4 4.5 \n"
# unit to cm
s_none += "mul 1.0e2\n"
# ######## Notice the Resp name and the filter band  ########
s_none += "w " + none_path + "\n"
s_none += "q \n"
subprocess.Popen(['sac'], stdin=subprocess.PIPE).communicate(s_none.encode())

# read data
tr_raw = read(raw_path)[0]
tr_vel = read(vel_path)[0]
tr_none = read(none_path)[0]

# for plot
data_raw, data_vel, data_none = tr_raw.data, tr_vel.data, tr_none.data
sampling_rate = tr_raw.stats.sampling_rate
t = np.arange(0, len(data_raw)) / sampling_rate
p_loc = tr_raw.stats.sac.t1 - tr_raw.stats.sac.b
s_loc = tr_raw.stats.sac.t2 - tr_raw.stats.sac.b

fig, ax = plt.subplots(3, 1, figsize=(10, 6), sharex=True)
figure_path = './test_response.png'

ax[0].plot(t, data_raw, color='black')
ax[0].set_xlim(0, np.max(t))
ax[0].axvline(x=p_loc, color='blue', linestyle='--')
ax[0].text(p_loc, np.min(data_raw), 'P', ha='right', color='blue', fontsize=15)
ax[0].axvline(x=s_loc, color='red', linestyle='--')
ax[0].text(s_loc, np.min(data_raw), 'S', ha='right', color='red', fontsize=15)

ax[1].plot(t, data_vel, color='black')
ax[1].set_xlim(0, np.max(t))
ax[1].axvline(x=p_loc, color='blue', linestyle='--')
ax[1].text(p_loc, np.min(data_vel), 'P', ha='right', color='blue', fontsize=15)
ax[1].axvline(x=s_loc, color='red', linestyle='--')
ax[1].text(s_loc, np.min(data_vel), 'S', ha='right', color='red', fontsize=15) 

ax[2].plot(t, data_none, color='black')
ax[2].set_xlim(0, np.max(t))
ax[2].axvline(x=p_loc, color='blue', linestyle='--')
ax[2].text(p_loc, np.min(data_none), 'P', ha='right', color='blue', fontsize=15)
ax[2].axvline(x=s_loc, color='red', linestyle='--')
ax[2].text(s_loc, np.min(data_none), 'S', ha='right', color='red', fontsize=15) 


plt.figtext(0.25, 0.95, "raw", ha="center", va="center", fontsize=15, color='blue')
plt.figtext(0.25, 0.64, "vel", ha="center", va="center", fontsize=15, color='blue')
plt.figtext(0.25, 0.32, "none", ha="center", va="center", fontsize=15, color='blue')    

plt.tight_layout()
plt.savefig(figure_path)
plt.close()