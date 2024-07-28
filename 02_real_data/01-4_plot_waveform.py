"""Plot nez and rtz waveforms.
"""
import os, glob, shutil
from obspy import read
import numpy as np
import matplotlib.pyplot as plt

################## params ##################
data_dir = './event_data_rtz'
out_dir   = './figures/rtz'
chns = ['z', 'r', 't']
# data_dir = './event_data_nez'
# out_dir   = './figures/nez'
# chns = ['Z', 'N', 'E']
if os.path.exists(out_dir): shutil.rmtree(out_dir)
os.makedirs(out_dir)
filter_range = [0.02, 0.2]
################## params ##################

def plot_waveforms(waveform_dir, sta_code, chns, figure_dir, filter_range):
    """
    @func: Plot waveforms.
    @params: Dir of event, code of station and directory for figures and filter range.
    """
    # read the data
    filepath_z = glob.glob(waveform_dir + '/' + sta_code + '*' + chns[0])[0]
    filepath_r = glob.glob(waveform_dir + '/' + sta_code + '*' + chns[1])[0]
    filepath_t = glob.glob(waveform_dir + '/' + sta_code + '*' + chns[2])[0]
    tr_z, tr_r, tr_t = read(filepath_z)[0], read(filepath_r)[0], read(filepath_t)[0]
    # bandpass filter
    tr_z.filter('bandpass', freqmin=filter_range[0], freqmax=filter_range[1])
    tr_r.filter('bandpass', freqmin=filter_range[0], freqmax=filter_range[1])
    tr_t.filter('bandpass', freqmin=filter_range[0], freqmax=filter_range[1])

    # read the data
    data_z, data_r, data_t = tr_z.data, tr_r.data, tr_t.data
    sampling_rate = tr_z.stats.sampling_rate
    t = np.arange(0, len(data_z)) / sampling_rate
    # P and S arrival location on the waveform
    p_loc = tr_z.stats.sac.t1 - tr_z.stats.sac.b
    s_loc = tr_z.stats.sac.t2 - tr_z.stats.sac.b

    # Plot the figure
    fig, ax = plt.subplots(3, 1, figsize=(10, 6), sharex=True)
    figure_name = '.'.join([sta_code, 'pdf'])
    figure_path = os.path.join(figure_dir, figure_name)

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
    plt.savefig(figure_path)
    print(waveform_dir, figure_path)
    plt.close()

    return         

if __name__ == '__main__':

    event_list = sorted(os.listdir(data_dir))
    for event in event_list:
        # directory for figures
        figure_dir = os.path.join(out_dir, event)
        if not os.path.exists(figure_dir): os.makedirs(figure_dir)
        # get the list of station id
        waveform_dir = os.path.join(data_dir, event)
        waveform_list = os.listdir(waveform_dir)
        sta_code_list = ['.'.join(waveform.split('.')[:2]) for waveform in waveform_list]
        sta_code_list = list(set(sta_code_list))

        for sta_code in sta_code_list:
            plot_waveforms(waveform_dir=waveform_dir, sta_code=sta_code, chns=chns,
                           figure_dir=figure_dir, filter_range=filter_range)
