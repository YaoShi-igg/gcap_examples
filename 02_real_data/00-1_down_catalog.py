"""Downdload catalogs.
"""
from obspy import UTCDateTime
from obspy.clients.fdsn import Client

########################## params #############################
client = Client("IRIS")
# the path to the file for catalogs
cat_path = './input/catalog.dat'
# parameters for events
t1 = UTCDateTime('2008-04-18T00:00:00')
t2 = UTCDateTime('2008-04-18T10:00:00')
# domain range
lat_range = [38, 39] 
lon_range = [-88, -87]
mag_range = [4, 6]
cat_params = {'starttime' : t1, 'endtime' : t2,
              'minlatitude' : lat_range[0],  'maxlatitude' : lat_range[1],
              'minlongitude' : lon_range[0], 'maxlongitude' : lon_range[1],
              'minmagnitude' : mag_range[0], 'maxmagnitude' : mag_range[1]}
########################## params #############################

def save_cat(cat_path, cat):
    """
    @func: Save the catalog from client.get_events to file
    @params: Path for out file and catalog(format: obspy.core.event.catalog.Catalog) 
    @return: None
    """
    cat_file = open(cat_path, 'w')
    for ev in cat:
        info_event = ev.origins
        orign_time = info_event[0].time
        lat = info_event[0].latitude
        lon = info_event[0].longitude
        # convert to kilometers
        dep = info_event[0].depth / 1e3
        info_mag = ev.magnitudes
        mag = info_mag[0].mag  
        cat_file.write('%s %.2f %.2f %d %.2f\n' % (orign_time, lon, lat, dep, mag))
    cat_file.close()

    return

if __name__ == '__main__': 

    # Download the catalog
    cat = client.get_events(**cat_params)
    # save the catalog in dat and xml
    cat.write('./input/catalog.xml', format='QUAKEML')
    save_cat(cat_path=cat_path, cat=cat)
