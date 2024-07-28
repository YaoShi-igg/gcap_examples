"""Rotate nez seismic records to rtz.
"""
import os, glob, shutil
from obspy import read
from obspy.signal.rotate import rotate_ne_rt

##########################params##########################
nez_dir = './event_data_nez'
rtz_dir = './event_data_rtz'
if os.path.exists(rtz_dir): shutil.rmtree(rtz_dir)
os.makedirs(rtz_dir)
##########################params##########################

event_list = os.listdir(nez_dir)
for event in event_list:
    waveform_dir = os.path.join(nez_dir, event)
    waveform_list = os.listdir(waveform_dir)
    sta_code_list = ['.'.join(waveform.split('.')[:2]) for waveform in waveform_list]
    sta_code_list = list(set(sta_code_list))
    
    for sta_code in sta_code_list:
        # read waveform
        waveform_n = glob.glob(waveform_dir + '/' + sta_code + '.??N')[0]
        waveform_e = glob.glob(waveform_dir + '/' + sta_code + '.??E')[0]
        waveform_z = glob.glob(waveform_dir + '/' + sta_code + '.??Z')[0]
        tr_n, tr_e= read(waveform_n)[0], read(waveform_e)[0]
        # rotate the waveform
        data_r, data_t = rotate_ne_rt(n=tr_n.data, e=tr_e.data, ba=tr_n.stats.sac.baz)
        # Save to rtz dir
        waveform_out_dir = os.path.join(rtz_dir, event)
        if not os.path.exists(waveform_out_dir): os.makedirs(waveform_out_dir)
        waveform_out_r = os.path.join(waveform_out_dir, sta_code + '.r')
        waveform_out_t = os.path.join(waveform_out_dir, sta_code + '.t')
        waveform_out_z = os.path.join(waveform_out_dir, sta_code + '.z')

        # Save the R and T data to trace 
        tr_r, tr_t = tr_n.copy(), tr_e.copy()
        tr_r.data, tr_t.data = data_r, data_t

        tr_r.write(waveform_out_r, format='SAC')
        tr_t.write(waveform_out_t, format='SAC')
        os.system('cp ' + waveform_z + ' ' + waveform_out_z)

