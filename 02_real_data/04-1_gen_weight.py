# Generate weight file for gcap
# fromat
# StationName dist w1 w2 w3 w4 w5 tp ts
# w1 to w5 are the weights for PnlZ，PnlR，Z，R，T
# tp for P wave arrival, but ts for initial shift of surface wave
import os, glob
import json
from obspy import read

######################### i/o paths ###########################
data_dir = "./event_data_rtz"
input_syn = "./input/syn.json"
with open(input_syn, 'r') as f:
    syn_dict = json.load(f)
######################### i/o paths ###########################
# params for default
w1, w2, w3, w4, w5 = 1, 1, 1, 1, 1
ts = 0

# get epicentral distance and P arrival from input file and waveform
source_list = list(syn_dict.keys())
for source in source_list:
    source_dir = os.path.join(data_dir, source)
    weight_file = os.path.join(source_dir, 'weight.dat')
    # get the receivers
    receiver_list = list(syn_dict[source]['receiver'].keys())

    with open(weight_file, 'w') as f:
        for receiver in receiver_list:
            dist = syn_dict[source]['receiver'][receiver]['dist']
            # get P wave arrival from sachead of waveform file
            waveform_z = glob.glob(source_dir + '/' + receiver + '.z')[0]
            tr_z = read(waveform_z)[0]
            tp = tr_z.stats.sac.t1
            f.write('%s %d %d %d %d %d %d %.1f %d\n' %
                    (receiver, dist, w1, w2, w3, w4, w5, tp, ts))
