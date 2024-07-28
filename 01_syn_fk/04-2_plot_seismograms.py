# Plot seismograms(RTZ).
import os, glob
from obspy import read
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib import mlab
from obspy.imaging.cm import obspy_sequential
import math 
################## i/o paths ##################
seis_lib = './seismograms'
out_dir   = './figures/seismograms'
################## i/o paths ##################

def get_sta_id_list(seis_list):
    """
    @func: Get list of station id from source.
    @params: List of seismograms
    """
    sta_id_list = []
    for seis in seis_list:
        sta_id = '.'.join(seis.split('.')[:2])
        if sta_id not in sta_id_list and 'gmt' not in sta_id and 'weight' not in sta_id:
            sta_id_list.append(sta_id)

    return sta_id_list
def _nearest_pow_2(x):
    """
    Find power of two nearest to x
    :type x: float
    :param x: Number
    :rtype: int
    :return: Nearest power of 2 to x
    """
    a = math.pow(2, math.ceil(np.log2(x)))
    b = math.pow(2, math.floor(np.log2(x)))
    if abs(a - x) < abs(b - x):
        return a
    else:
        return b

def spectrogram(data, samp_rate, axes, per_lap=0.9, wlen=None, 
                dbscale=False, mult=8.0, cmap=obspy_sequential, 
                zorder=None, clip=[0.0, 1.0]):
    """
    Computes and plots spectrogram of the input data.

    :param data: Input data
    :type samp_rate: float
    :param samp_rate: Samplerate in Hz
    :type axes: :class:`matplotlib.axes.Axes`
    :param axes: Plot into given axes, this deactivates the fmt and
        outfile option.
    :type per_lap: float
    :param per_lap: Percentage of overlap of sliding window, ranging from 0
        to 1. High overlaps take a long time to compute.
    :type wlen: int or float
    :param wlen: Window length for fft in seconds. If this parameter is too
        small, the calculation will take forever. If None, it defaults to a
        window length matching 128 samples.
    :type dbscale: bool
    :param dbscale: If True 10 * log10 of color values is taken, if False the
        sqrt is taken.
    :type mult: float
    :param mult: Pad zeros to length mult * wlen. This will make the
        spectrogram smoother.
    :type cmap: :class:`matplotlib.colors.Colormap`
    :param cmap: Specify a custom colormap instance. If not specified, then the
        default ObsPy sequential colormap is used.
    :type zorder: float
    :param zorder: Specify the zorder of the plot. Only of importance if other
        plots in the same axes are executed.
    :type clip: [float, float]
    :param clip: adjust colormap to clip at lower and/or upper end. The given
        percentages of the amplitude range (linear or logarithmic depending
        on option `dbscale`) are clipped.
    """
    import matplotlib.pyplot as plt
    # enforce float for samp_rate
    samp_rate = float(samp_rate)

    # set wlen from samp_rate if not specified otherwise
    if not wlen:
        wlen = 128 / samp_rate
    npts = len(data)

    # nfft needs to be an integer, otherwise a deprecation will be raised
    # XXX add condition for too many windows => calculation takes for ever
    nfft = int(_nearest_pow_2(wlen * samp_rate))

    if npts < nfft:
        msg = (f'Input signal too short ({npts} samples, window length '
               f'{wlen} seconds, nfft {nfft} samples, sampling rate '
               f'{samp_rate} Hz)')
        raise ValueError(msg)

    if mult is not None:
        mult = int(_nearest_pow_2(mult))
        mult = mult * nfft
    nlap = int(nfft * float(per_lap))

    data = data - data.mean()

    # Here we call not plt.specgram as this already produces a plot
    # matplotlib.mlab.specgram should be faster as it computes only the
    # arrays
    # XXX mlab.specgram uses fft, would be better and faster use rfft
    specgram, freq, time = mlab.specgram(data, Fs=samp_rate, NFFT=nfft,
                                         pad_to=mult, noverlap=nlap)
    # print(time)

    if len(time) < 2:
        msg = (f'Input signal too short ({npts} samples, window length '
               f'{wlen} seconds, nfft {nfft} samples, {nlap} samples window '
               f'overlap, sampling rate {samp_rate} Hz)')
        raise ValueError(msg)

    # db scale and remove zero/offset for amplitude
    if dbscale:
        specgram = 10 * np.log10(specgram[1:, :])
    else:
        specgram = np.sqrt(specgram[1:, :])
    freq = freq[1:]

    vmin, vmax = clip
    if vmin < 0 or vmax > 1 or vmin >= vmax:
        msg = "Invalid parameters for clip option."
        raise ValueError(msg)
    _range = float(specgram.max() - specgram.min())
    vmin = specgram.min() + vmin * _range
    vmax = specgram.min() + vmax * _range
    norm = Normalize(vmin, vmax, clip=True)

    # calculate half bin width
    halfbin_time = (time[1] - time[0]) / 2.0
    halfbin_freq = (freq[1] - freq[0]) / 2.0

    # argument None is not allowed for kwargs on matplotlib python 3.3
    kwargs = {k: v for k, v in (('cmap', cmap), ('zorder', zorder))
              if v is not None}

    # pcolor expects one bin more at the right end
    freq = np.concatenate((freq, [freq[-1] + 2 * halfbin_freq]))
    time = np.concatenate((time, [time[-1] + 2 * halfbin_time]))
    # center bin
    time -= halfbin_time
    freq -= halfbin_freq

    # rasterized=True, delete the grid lines
    cax = axes.pcolormesh(time, freq, specgram, norm=norm, rasterized=True, **kwargs)
    axes.set_yscale('log')
    axes.set_ylabel('Frequency (Hz)')
    axes.grid(False)

    # set the colorbar(Normalized)
    # cbar = plt.colorbar(cax, ax=axes, label='Normalized Amplititude', location='bottom', shrink=0.25)
    # cbar.set_ticks(np.linspace(0, vmax, 6))
    # cbar.set_ticklabels(['0', '0.2', '0.4', '0.6', '0.8', '1.0'])

    return

def plot_seismograms(source, sta_id, figure_dir, filter_range):
    """
    @func: Plot seismograms.
    @params: Name of source, id of station and directory for figures.
    """
    # read the data
    seis_dir = os.path.join(seis_lib, source)
    print(seis_dir + '/' + sta_id + '.z')
    filepath_z = glob.glob(seis_dir + '/' + sta_id + '.z')[0]
    filepath_r = glob.glob(seis_dir + '/' + sta_id + '.r')[0]
    filepath_t = glob.glob(seis_dir + '/' + sta_id + '.t')[0]
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
    fig, ax = plt.subplots(4, 1, figsize=(10, 6), sharex=True)
    figure_name = '.'.join([sta_id, 'pdf'])
    figure_path = os.path.join(figure_dir, figure_name)

    # plot the seismic record
    ax[0].plot(t, data_z, color='black')
    ax[0].set_xlim(0, np.max(t))
    ax[0].set_ylabel('Velocity (cm/s)')
    ax[0].axvline(x=p_loc, color='blue', linestyle='--')
    ax[0].text(p_loc, np.min(data_z), 'P', ha='right', color='blue', fontsize=15)
    ax[0].axvline(x=s_loc, color='red', linestyle='--')
    ax[0].text(s_loc, np.min(data_z), 'S', ha='right', color='red', fontsize=15)

    # plot the spectrogram for R channel data
    spectrogram(data=data_z, samp_rate=sampling_rate, axes=ax[1], cmap='hot')

    ax[2].plot(t, data_r, color='black')
    ax[2].set_xlim(0, np.max(t))
    ax[2].set_ylabel('Velocity (cm/s)')
    ax[2].axvline(x=p_loc, color='blue', linestyle='--')
    ax[2].text(p_loc, np.min(data_r), 'P', ha='right', color='blue', fontsize=15)
    ax[2].axvline(x=s_loc, color='red', linestyle='--')
    ax[2].text(s_loc, np.min(data_r), 'S', ha='right', color='red', fontsize=15) 

    ax[3].plot(t, data_t, color='black')
    ax[3].set_xlim(0, np.max(t))
    ax[3].set_xlabel('Time (s)')
    ax[3].set_ylabel('Velocity (cm/s)')
    ax[3].axvline(x=p_loc, color='blue', linestyle='--')
    ax[3].text(p_loc, np.min(data_t), 'P', ha='right', color='blue', fontsize=15)
    ax[3].axvline(x=s_loc, color='red', linestyle='--')
    ax[3].text(s_loc, np.min(data_t), 'S', ha='right', color='red', fontsize=15)     

    plt.tight_layout()
    plt.savefig(figure_path)
    plt.close()

    return

if __name__ == '__main__':

    filter_range = [0.02, 0.2]
    source_list = sorted(os.listdir(seis_lib))
    for source in source_list:
        # directory for figures
        figure_dir = os.path.join(out_dir, source)
        if not os.path.exists(figure_dir): os.makedirs(figure_dir)
        # get the list of station id
        seis_dir = os.path.join(seis_lib, source)
        seis_list = os.listdir(seis_dir)
        sta_id_list = get_sta_id_list(seis_list=seis_list)

        for sta_id in sta_id_list:
            plot_seismograms(source=source, sta_id=sta_id, figure_dir=figure_dir, filter_range=filter_range)
