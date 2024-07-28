"""Download waveforms, response for event by catalogs, write station_info to station file.
"""
import re, os, shutil
from obspy.clients.fdsn import Client
from obspy import read_events

############################### i/o paths ########################
cat_path = './input/catalog.xml'
data_dir = './raw_data'
if os.path.exists(data_dir): shutil.rmtree(data_dir)
os.makedirs(data_dir)
pz_dir = './pz'
if os.path.exists(pz_dir): shutil.rmtree(pz_dir)
os.makedirs(pz_dir)
sta_file = open('./input/station.dat', 'w')
data_f = open('./data_fail.dat', 'w')
############################### i/o paths ########################

############################### params ##############################
client = Client("IRIS")
# download waveforms 10s before orign and 100s after orign
time_range = [0, 100]
# specify stations
# distance between 0-3 degrees
radius = [1, 3]
net = 'IU,NM'
sta = '*'
loc = '*'
chn = 'BH*'
############################### params ##############################
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

def get_stations(t0, t1, evla, evlo, radius):
    """
    @func: Get station list in specific distance range.
    @params: Start and end time, location of source and distance range.
    @return: List of stations. 
    """
    global net, sta, loc, chn
    # get the station 
    inventory = client.get_stations(network=net, station=sta, 
                    channel=chn, location=loc, starttime=t0, endtime=t1,
                    latitude=evla, longitude=evlo, minradius=radius[0], maxradius=radius[1])
    sta_list = inventory.get_contents()['stations']
    sta_list = [re.sub(r'\([^)]*\)', '', entry).strip() for entry in sta_list]
    print('Get stations within %d-%d degrees from source.' % (radius[0], radius[1]))

    return sta_list

def get_waveforms(orign, t0, t1, sta_list):
    """
    @func: Download event waveforms for specific stattions and time range.
    @params: Start and end time, list of stations.
    @return
    """
    global loc, chn
    # set the dir name
    event_id = dtime2str(orign)
    waveform_dir = os.path.join(data_dir, event_id)
    if not os.path.exists(waveform_dir): os.makedirs(waveform_dir)
    print('Download waveforms for event: %s' % (event_id))

    for sta_id in sta_list:
        net, sta = sta_id.split('.')[0], sta_id.split('.')[1]
        try:
            # get waveforms
            st = client.get_waveforms(network=net, station=sta, location=loc, 
                                      channel=chn, starttime=t0, endtime=t1)            
            for i in range(len(st)):
                # Save waveforms
                waveform_name = '.'.join([net, sta, st[i].stats.channel])
                out_path = os.path.join(waveform_dir, waveform_name)
                st[i].write(out_path, format="SAC")

                # get PZ file for every channel
                response = client.get_stations(network=st[i].stats.network, station=st[i].stats.station,
                                location=st[i].stats.location, channel=st[i].stats.channel,
                                starttime=t0, endtime=t1, level="response")
                # Save pz file
                pz_name = '.'.join([st[i].stats.network, st[i].stats.station, st[i].stats.channel, 'PZ'])
                pz_path = os.path.join(pz_dir, pz_name)
                response.write(pz_path, format="SACPZ")

            # extract station information and save to file
            stla = response.networks[0].stations[0].latitude
            stlo = response.networks[0].stations[0].longitude
            stel = response.networks[0].stations[0].elevation / 1e3
            sta_file.write('%s %.2f %.2f %.2f\n' % (sta_id, stlo, stla,stel))
            print("Get the data and response for %s.%s" % (net, sta))
        except:
            print('%s.%s has no data' % (net, sta))
            data_f.write('{},{},{}\n'.format(event_id, net, sta))
    sta_file.close()
    data_f.close()

    return 

if __name__ == '__main__':

    # load the catalog(in the format of catalog)
    cat = read_events(cat_path)
    for ev in cat:
        # read the information of event
        info_event = ev.origins
        orign_time = info_event[0].time
        evla = info_event[0].latitude
        evlo = info_event[0].longitude 
        # time before and after the event
        t0 = orign_time + time_range[0]
        t1 = orign_time + time_range[1]

        # get list of stations
        sta_list = get_stations(t0=t0, t1=t1, evla=evla, evlo=evlo, radius=radius)
        # Download waveforms
        get_waveforms(orign=orign_time, t0=t0, t1=t1, sta_list=sta_list)

