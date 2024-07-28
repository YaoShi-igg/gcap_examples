from obspy import read

data_n = read('./event_data_nez/20080418093658680/IU.CCM.BHN')[0].data
data_r = read('./event_data_rtz/20080418093658680/IU.CCM.r')[0].data
data_e = read('./event_data_nez/20080418093658680/IU.CCM.BHE')[0].data
data_t = read('./event_data_rtz/20080418093658680/IU.CCM.t')[0].data

print(max(data_n), max(data_e), max(data_r), max(data_t))
