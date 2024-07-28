"""Preprocess for event data.
"""
# 1. Modify the ref.time to earthquake orign time
# 2. Add the station, earthquake and arrival infomation to the header
# 3. Downsampling to 10hz, rmean, rtrend, taper
# 4. Remove response by pz file and modify unit to cm/s
import os, shutil
import subprocess
from obspy import read, UTCDateTime
from obspy.geodetics import locations2degrees
from obspy.io.sac.util import obspy_to_sac_header
from obspy.taup import TauPyModel
from obspy.taup.taup_create import build_taup_model
build_taup_model("./input/cus.nd") 
# here use ak135 for test, for real data, we can use local model
model = TauPyModel(model="cus") 
##########################params##########################
# path for raw_data
raw_data = './raw_data'
# path for out data
out_data = './event_data_nez'
if os.path.exists(out_data): shutil.rmtree(out_data)
os.makedirs(out_data)
# dir for pz file
pz_dir = './pz'
# read catalog
cat_path = './input/catalog.dat'
f = open(cat_path); cat_lines = f.read().splitlines(); f.close
# read stations
sta_path = './input/station.dat'
g = open(sta_path); sta_lines = g.read().splitlines(); g.close()
# sampling rate for dest data
dest_sampling_rate = 10 #Hz
##########################params##########################
def dtime2str(utc_time):
    """
    @func: Convert utc time to string.
    @param: Utc time.
    @return: String for utc time. 
    """
    date = ''.join(str(utc_time).split('T')[0].split('-'))
    time = ''.join(str(utc_time).split('T')[1].split(':'))[:10]
    time_all = date + time
    str_time = ''.join(time_all.split('.'))

    return str_time

def get_catalog(event_name):
    """
    @func: Get lat, lon, dep and mag information from catalog file
    @params: Name of the event (YYYYMMDDHHMMSS)
    @return: Origin, latitude, longitude, depth and magnetitude of the event
    """
    global cat_lines
    for cat_line in cat_lines: 
        cat = cat_line.split()
        orign, evlo, evla, evdp, mag = cat[0], float(cat[1]), float(cat[2]), float(cat[3]), float(cat[4])
        if event_name == dtime2str(orign):
            return UTCDateTime(orign), evlo, evla, evdp, mag 
        
def get_sta(sta_id):
    """
    @func: Get stla, stlo, stel info from station file
    @params: Network and station
    @return: stla, stlo, stel
    """
    global sta_lines
    for sta_line in sta_lines:
        sta = sta_line.split()
        sta_code, stlo, stla, stel = sta[0], float(sta[1]), float(sta[2]), float(sta[3])
        if sta_code == sta_id:          
            # print(stlo, stla, stel)
            return stlo, stla, stel * 1e3
        
def cal_arrival(dist_degree, dep):
    """
    @func: Calculate P and S arrival
    @params: Dstance and depth
    @return: P and S arrival
    """
    # arrivals for p wave
    p_arrivals = model.get_travel_times(source_depth_in_km=dep, 
                                        distance_in_degree=dist_degree, phase_list=['P', 'Pn', 'Pg', 'p', 'Pdiff'])             
    # get the first arrival time of p wave
    p_arr = min(p_arrival.time for p_arrival in p_arrivals)
    # arrivals for s wave
    s_arrivals = model.get_travel_times(source_depth_in_km=dep, 
                                        distance_in_degree=dist_degree, phase_list=['S', 'Sn', 'Sg', 's', 'Sdiff'])              
    # get the first arrival time of p wave
    s_arr = min(s_arrival.time for s_arrival in s_arrivals)

    return round(p_arr, 2), round(s_arr, 2)

def ch_head(tr, orign, evla, evlo, evdp, mag, stla, stlo, stel, p_arr, s_arr):
    """
    @func: Change head for event waveform.
    @params: Trace for waveform, orign, lon, lat, depth(km), mag of event, 
             lon, lat, evel(m) of station, P and S arrival
    @return:  
    """
    tr.stats.sac = obspy_to_sac_header(tr.stats)
    # open the auto cal
    tr.stats.sac.lcalda = True
    # ref.time
    tr.stats.sac['nzyear'] = UTCDateTime(orign).year
    tr.stats.sac['nzjday'] = UTCDateTime(orign).julday
    tr.stats.sac['nzhour'] = UTCDateTime(orign).hour
    tr.stats.sac['nzmin']  = UTCDateTime(orign).minute
    tr.stats.sac['nzsec']  = UTCDateTime(orign).second
    tr.stats.sac['nzmsec'] = int(UTCDateTime(orign).microsecond/1e3)
    tr.stats.sac['o'] = 0
    # start time
    tr.stats.sac['b'] = '0'
    # station info
    tr.stats.sac.stla = stla
    tr.stats.sac.stlo = stlo
    tr.stats.sac.stel = stel
    # event info
    tr.stats.sac.evla = evla
    tr.stats.sac.evlo = evlo
    tr.stats.sac.evdp = evdp
    tr.stats.sac.mag = mag 
    # arrival info
    tr.stats.sac.t1 = p_arr
    tr.stats.sac.t2 = s_arr

    return

def match_pz(pz_dir, waveform):
    """
    @func: Mathch the waveform file with pz file
    @params: Directory for pz file and name of waveform
    @return: Path for pz file
    """
    pz_list = os.listdir(pz_dir)
    net = waveform.split('.')[0]
    sta = waveform.split('.')[1]
    chn = waveform.split('.')[2]
    # print(net, sta, chn)
    for pz in pz_list:
        if pz.split('.')[0] == net and pz.split('.')[1] == sta and pz.split('.')[2] == chn:
            return pz_dir + '/' + pz

def rm_response(file_path, pz):
    """
    @func: Remove resoponse for waveform
    @params: Path for waveform and pz file
    @return: None
    """
    #remove response
    os.putenv("SAC_DISPLAY_COPYRIGHT", '0')
    s = "r " + file_path + "\n"
    s += "rmean; rtr; taper\n"
    # ######## Notice the Resp name and the filter band  ########
    s += "trans from polezero subtype " + pz +" to vel freq 0.004 0.005 4 4.5 \n"
    # unit to cm/s (* 1e9 / 1e7)
    s += "mul 1.0e2\n"
    # ######## Notice the Resp name and the filter band  ########
    s += "w " + file_path + "\n"
    s += "q \n"
    subprocess.Popen(['sac'], stdin=subprocess.PIPE).communicate(s.encode())

if __name__ == '__main__':

    event_list = os.listdir(raw_data)
    for event in event_list:
        # get information for event
        orign, evlo, evla, evdp, mag  = get_catalog(event_name=event)
        waveform_dir = os.path.join(raw_data, event)
        waveform_list = os.listdir(waveform_dir)

        for waveform in waveform_list:
            sta_id = '.'.join(waveform.split('.')[:2])
            waveform_path = os.path.join(waveform_dir, waveform)
            # get information for station
            stlo, stla, stel = get_sta(sta_id=sta_id)
            # get arrivals
            dist_degree = locations2degrees(stla, stlo, evla, evlo)
            p_arr, s_arr = cal_arrival(dist_degree=dist_degree, dep=evdp)

            tr = read(waveform_path)[0]
            # downsamlpling to 10hz
            if tr.stats.sampling_rate != dest_sampling_rate:
                decimate_factor = int(tr.stats.sampling_rate / dest_sampling_rate)
                tr.decimate(decimate_factor, strict_length=False, no_filter=True)
            # change head
            ch_head(tr=tr, orign=orign, evla=evla, evlo=evlo, evdp=evdp, mag=mag, 
                    stla=stla, stlo=stlo, stel=stel, p_arr=p_arr, s_arr=s_arr)

            # Save the waveform data to out path
            waveform_out_dir = os.path.join(out_data, event)
            if not os.path.exists(waveform_out_dir): os.makedirs(waveform_out_dir)
            waveform_out_path = os.path.join(waveform_out_dir, waveform)
            tr.write(waveform_out_path, format='SAC')

            # remove the response
            # match the resp file
            pz = match_pz(pz_dir=pz_dir, waveform=waveform)
            # remove the response
            rm_response(file_path=waveform_out_path, pz=pz)
