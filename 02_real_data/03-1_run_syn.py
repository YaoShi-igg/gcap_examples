# Green's function convolution with source time function
# Combine different green function to get RTZ seismograms based on focal mechanism

import os, shutil
import json

######################### i/o paths ###########################
glib = "./glib"
seislib = "./syn_seismograms"
if os.path.exists(seislib): shutil.rmtree(seislib)
os.makedirs(seislib)
input_syn = "./input/syn.json"
with open(input_syn, 'r') as f:
    syn_dict = json.load(f)
######################### i/o paths ###########################

########################## params #############################
# model for synthesizing green function
model='cus'
# params for differen source
# mag(dyne-cm) for iso; mag(dyne), strike and dip for single force
# mag(mw), strike, dip and rake for dc, mag/Mxx/Mxy/Mxz/Myy/Myz/Mzz for moment tensor
# here for dc and mag in syn_dict
strike, dip, rake = '125', '69', '-17'
# bandpass filter, generally no need to filter
# f1=0.1; f2=1; n=4
# source time function (D for trapezoid; S for sac file)
d = '1'
########################## params #############################
source_list = list(syn_dict.keys())
for source in source_list:
    # make directory for source and run syn in it
    source_dir = os.path.join(seislib, source)
    if not os.path.exists(source_dir): os.makedirs(source_dir)
    evdp, mag = str(syn_dict[source]['evdp']), str(syn_dict[source]['mag'])
    
    # get the receivers
    receiver_list = list(syn_dict[source]['receiver'].keys())
    for receiver in receiver_list:
        # out name for sac file
        out_name = receiver + '.z'
        dist = str(syn_dict[source]['receiver'][receiver]['dist'])
        az = str(syn_dict[source]['receiver'][receiver]['az'])

        # parameters for syn
        M = '-M' + '/'.join([mag, strike, dip, rake])
        D = '-D' + str(d)
        A = '-A' + az
        G = '-G' + '/'.join([glib, model, model + '_' + evdp, dist + '.grn.0'])
        O = '-O' + out_name
        cmd = ' '.join(['syn', M, D, A, G, O])
        print(cmd)
        os.system(cmd)
    os.system('mv ??.*.[rtz] ' + source_dir)
        # print('STATUS: {}'.format(syn_output))
