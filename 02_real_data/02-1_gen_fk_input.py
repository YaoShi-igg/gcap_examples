# Generate input file for fk
# file format for fk
# source_id receiver_id dist azimuth evdp mag 
from obspy.geodetics import gps2dist_azimuth
import json

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

def read_source(source_file):
    """
    @func: Read information for source.
    @params: Path for source file.
    @return: Dictionary for source.
    """
    source_dict = {}
    f = open(source_file, 'r')
    lines = f.read().splitlines()
    for line in lines:
        line = line.split()
        evid, evlo, evla, evdp, mag = line[0], float(line[1]), float(line[2]), float(line[3]), float(line[4])
        evid = dtime2str(evid)
        source_dict[evid] = {'evlo': evlo, 'evla': evla, 'evdp': evdp, 'mag':  mag}
    return source_dict

def read_receiver(receiver_file):
    """
    @func: Read information for receiver.
    @params: Path for receiver file.
    @return: Dictionary for receiver.
    """
    receiver_dict = {}
    f = open(receiver_file, 'r')
    lines = f.read().splitlines()
    for line in lines:
        line = line.split()
        sta_id, stlo, stla = line[0], float(line[1]), float(line[2])
        receiver_dict[sta_id] = {'stlo': stlo, 'stla': stla}
    return receiver_dict

def gen_fk_input(source, receiver):
    """
    @func: Convert source and reivers to the format for fk.
    @params: Dictionary for source and receiver.
    @return: Save input file for fk and generate input file for syn.
    """
    source_id_list   = list(source.keys())
    receiver_id_list = list(receiver.keys())
    # dcitionary for syn and save to json file
    syn_dict = {}

    f = open(input_fk, 'w')
    for source_id in source_id_list:
        evdp, mag = round(source[source_id]['evdp']), source[source_id]['mag']
        syn_dict[source_id] = {'evdp': evdp, 'mag': mag}
        syn_dict[source_id]['receiver'] = {}

        for receiver_id in receiver_id_list:
            # meters, degrees, degrees
            dist, az, baz = gps2dist_azimuth(lat1=source[source_id]['evla'], 
                                             lon1=source[source_id]['evlo'], 
                                             lat2=receiver[receiver_id]['stla'], 
                                             lon2=receiver[receiver_id]['stlo'], 
                                             a=r)
            dist, az = round(dist / 1e3), round(az, 1)
            syn_dict[source_id]['receiver'][receiver_id]={'dist': dist, 'az': az}                        
            # write to file 
            f.write('%s %s %d %.1f %d %.1f\n' % 
                    (source_id, receiver_id, dist, az, evdp, mag))

    return syn_dict

if __name__ == '__main__':

    # input file
    source_file   = './input/catalog.dat'
    receiver_file = './input/station.dat'
    # out_file
    input_fk  = './input/fk.in'
    input_syn = './input/syn.json'
    # earth_radius, meters, same to setting in fk
    r = 6371000

    # read source and receiver
    source_dict = read_source(source_file=source_file)
    receiver_dict = read_receiver(receiver_file=receiver_file)

    # Convert the format and save to fk.in
    syn_dict = gen_fk_input(source=source_dict, receiver=receiver_dict)
    with open(input_syn, 'w') as f:
        json.dump(syn_dict, f)
