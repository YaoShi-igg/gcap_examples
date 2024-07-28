'''Pick the arrival times.
If run it, use picked arrival for weight.dat, else, use theoritical.
'''
import os, glob
from obspy import read
import matplotlib.pyplot as plt
import numpy as np

###################################################params ####################################
# windows or linux
# platform = 'windows'
platform = 'windows'
# i/o paths
data_dir = ('./event_data_rtz')
event_index = 5
chns = ['r', 't', 'z']
# index of the event
event_index = 0
filter_range = [0.02, 0.2]
###################################################params ####################################

def Pick(event):
    # Left mouse button
    if event.button == 1:
        global phase
        phase.append(event.xdata)
        print("Phase arrivals:", event.xdata)
    # scroll wheel up
    elif event.button == 4: 
        phase=[]
        print('Clear picks, click left button again')
    # Right mouse button
    elif event.button == 3:
        plt.close()
        if len(phase) == 0:
            print('No picks')

def plot_waveforms(tr_z, tr_r, tr_t, sta_code, filter_range):
    """
    @func: Plot waveforms.
    @params: Dir of event, code of station and directory for figures and filter range.
    """
    # bandpass filter
    tr_z.filter('bandpass', freqmin=filter_range[0], freqmax=filter_range[1])
    tr_r.filter('bandpass', freqmin=filter_range[0], freqmax=filter_range[1])
    tr_t.filter('bandpass', freqmin=filter_range[0], freqmax=filter_range[1])

    # read the data
    data_z, data_r, data_t = tr_z.data, tr_r.data, tr_t.data
    sampling_rate = tr_z.stats.sampling_rate
    t = np.arange(0, len(data_z)) / sampling_rate
    # P and S arrival location on the waveform
    p_loc = tr_z.stats.sac.user9 - tr_z.stats.sac.b
    s_loc = tr_z.stats.sac.user8 - tr_z.stats.sac.b

    # Plot the figure
    fig, ax = plt.subplots(3, 1, figsize=(10, 6), sharex=True)
    # plot the seismic record
    ax[0].plot(t, data_z, color='black')
    ax[0].set_xlim(0, np.max(t))
    ax[0].set_ylabel('Velocity (cm/s)')
    ax[0].axvline(x=p_loc, color='blue', linestyle='--')
    ax[0].text(p_loc, np.min(data_z), 'P', ha='right', color='blue', fontsize=15)
    ax[0].axvline(x=s_loc, color='red', linestyle='--')
    ax[0].text(s_loc, np.min(data_z), 'S', ha='right', color='red', fontsize=15)

    ax[1].plot(t, data_r, color='black')
    ax[1].set_xlim(0, np.max(t))
    ax[1].set_ylabel('Velocity (cm/s)')
    ax[1].axvline(x=p_loc, color='blue', linestyle='--')
    ax[1].text(p_loc, np.min(data_r), 'P', ha='right', color='blue', fontsize=15)
    ax[1].axvline(x=s_loc, color='red', linestyle='--')
    ax[1].text(s_loc, np.min(data_r), 'S', ha='right', color='red', fontsize=15) 

    ax[2].plot(t, data_t, color='black')
    ax[2].set_xlim(0, np.max(t))
    ax[2].set_xlabel('Time (s)')
    ax[2].set_ylabel('Velocity (cm/s)')
    ax[2].axvline(x=p_loc, color='blue', linestyle='--')
    ax[2].text(p_loc, np.min(data_t), 'P', ha='right', color='blue', fontsize=15)
    ax[2].axvline(x=s_loc, color='red', linestyle='--')
    ax[2].text(s_loc, np.min(data_t), 'S', ha='right', color='red', fontsize=15)     

    plt.tight_layout()
    plt.suptitle(sta_code, ha='center', va='top', fontsize=15, color='blue')

    return fig

if __name__ == '__main__':

    event_list = os.listdir(data_dir)
    event = event_list[event_index]

    waveform_dir = os.path.join(data_dir, event)
    waveform_list = os.listdir(waveform_dir)
    # print(waveform_list, len(waveform_list))

    sta_code_list = ['.'.join(waveform.split('.')[:2]) for waveform in waveform_list]
    sta_code_list = list(set(sta_code_list))
    if 'weight.dat' in sta_code_list:
        sta_code_list.remove('weight.dat')
    print('%s stations for this event %s' % (len(sta_code_list), event))

    for sta_code in sta_code_list:
        print(sta_code)
        # read the data
        filepath_z = glob.glob(waveform_dir + '/' + sta_code + '*' + chns[0])[0]
        filepath_r = glob.glob(waveform_dir + '/' + sta_code + '*' + chns[1])[0]
        filepath_t = glob.glob(waveform_dir + '/' + sta_code + '*' + chns[2])[0]
        tr_z, tr_r, tr_t = read(filepath_z)[0], read(filepath_r)[0], read(filepath_t)[0]
        fig = plot_waveforms(tr_z=tr_z, tr_r=tr_r, tr_t=tr_t, sta_code=sta_code, filter_range=filter_range)

        # save the picks
        phase=[]
        fig.canvas.mpl_connect('button_press_event', Pick)
        plt.show()
        tr_z.stats.sac['t1'] = phase[0] + tr_z.stats.sac['b']
        tr_r.stats.sac['t1'] = phase[0] + tr_z.stats.sac['b']
        tr_t.stats.sac['t1'] = phase[0] + tr_z.stats.sac['b']
        tr_z.stats.sac['t2'] = phase[1] + tr_z.stats.sac['b']
        tr_r.stats.sac['t2'] = phase[1] + tr_z.stats.sac['b']
        tr_t.stats.sac['t2'] = phase[1] + tr_z.stats.sac['b']
        # print(phase[0], tr_z.stats.sac['b'], tr_z.stats.sac['user9'])
        # clear picks
        phase=[]
        # save the change
        tr_z.write(filepath_z, format='SAC')
        tr_r.write(filepath_r, format='SAC')
        tr_t.write(filepath_t, format='SAC')
