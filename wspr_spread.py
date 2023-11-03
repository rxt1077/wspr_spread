import struct
import sys
import numpy as np
import wspr
import wave
import matplotlib.pyplot as plt
import argparse
import csv

# constants used in analysis
# we're working with a 45000 sample FFT that covers 375 Hz
# and is sampled at 375 Hz
SAMPLES = 46080
SAMPLE_RATE = 375
BW = 375
# DF allows us to convert from freqs to an index
DF = SAMPLE_RATE / SAMPLES
# all our calcs on g() are around the middle 
MIDDLE_INDEX = round(SAMPLES / 2)

# the region of interest is from -4 Hz to 4 Hz
REGION_START = round(MIDDLE_INDEX - 4 / DF)
REGION_END = round(MIDDLE_INDEX + 4 / DF)
REGION_MIDDLE = round((REGION_END - REGION_START) / 2)

# the leftmost 2 Hz of the region of interest
LEFT_NOISE_START = round(REGION_MIDDLE - 4 / DF)
LEFT_NOISE_END = round(REGION_MIDDLE - 2 / DF)

# the rightmost 2 Hz of the region of interest
RIGHT_NOISE_START = round(REGION_MIDDLE + 2 / DF)
RIGHT_NOISE_END = round(REGION_MIDDLE + 4 / DF)

# the power region is the middle 2 Hz of the region of interest
POWER_START = round(REGION_MIDDLE - 1 / DF)
POWER_END = round(REGION_MIDDLE + 1 / DF)

class WSPRSpreadError(Exception):
    """Basic exception class for errors"""

def plot_fft(filename, title, freqs, fft):
    """Utility function to plot things"""

    fig = plt.figure()
    plt.title(title)
    plt.xlabel('Freq (Hz)')
    plt.ylabel('Magnitude')
    plt.plot(freqs, fft)
    plt.savefig(filename)
    plt.close(fig)

def read_wav_file(fn):
    """Reads a WAV file the same way readwavfile() does in wsprd.c
    Currently only 2 minute mode is supported"""

    npoints = 114 * 12000
    nfft2 = 46080
    nh2 = nfft2 / 2
    nfft1 = nfft2 * 32
    df = 12000.0 / nfft1
    i0 = int(1500.0 / df + 0.5)

    realin = np.zeros(nfft1, dtype=float)

    # read WAV file, check header, and scale values
    with wave.open(fn, 'rb') as fp:
        params = fp.getparams()
        if params.nchannels != 1:
            raise WSPRSpreadError(f"Invalid number of WAV channels: {params.nchannels}")
        if params.sampwidth != 2:
            raise WSPRSpreadError(f"Invalid WAV sample width: {params.sampwidth}")
        if params.framerate != 12000:
            raise WSPRSpreadError(f"Invalid WAV framerate: {params.framerate}")
        if params.nframes != 1440000:
            raise WSPRSpreadError(f"Invalid number of frames in WAV: {params.nframes}")
        if params.comptype != 'NONE':
            raise WSPRSpreadError("Only uncompressed WAVs are supported")
        for i in range(0, npoints):
            realin[i] = int.from_bytes(fp.readframes(1), "little", signed=True) / 32768.0

    # real to complex fft yielding 737281 complex numbers
    fftout = np.fft.rfft(realin)

    fftin = np.zeros(nfft2, dtype=complex)

    # i0 is the index of 1500 Hz in the DFT
    # nh2 is half of nfft2 (23040)
    for i in range(0, nfft2):
        j = i0 + i
        if i > nh2:
            j -= nfft2
        fftin[i] = fftout[j]
    # the first half of fftin is fftout[184320 - 207360] (1500 Hz to 1688 Hz)
    # the second half of fftin is fftout[161281 - 184319] (1313 Hz to 1500 hz)
    # fftin holds 375 Hz of bandwidth around 1500 Hz

    # inverse fft of that 375 Hz subsection
    # looks like the FFTW inverse transform is unscaled:
    # https://www.fftw.org/fftw3_doc/The-1d-Discrete-Fourier-Transform-_0028DFT_0029.html
    # https://numpy.org/doc/stable/reference/routines.fft.html#normalization
    fftout = np.fft.ifft(fftin, norm="forward")

    # wsprd does a bit of scaling
    iq = fftout / 1000

    return iq

def reference(f0, symbols, callsign):
    """Creates a reference WSPR signal vector centered on f0
    Based on subtract_signal2() from wsprd.c.
    wsprd.c uses single precision floats, if you want similar results
    you need to pay _very_ close attention to the precision of your
    variables."""

    dt = np.float32(1.0 / SAMPLE_RATE)
    df = np.float32(SAMPLE_RATE / 256.0)
    nspersym = 256
    twopidt = np.float32(np.float32(2.0) * np.pi * dt)
    phi = np.float32(0.0)
    tone_separation = np.float32(12000 / 8192)
    #tone_separation = np.float32(1.5)

    reference = np.zeros(SAMPLES, dtype=complex)
    with open(f"wsprd/{callsign}_ref_params.csv", newline='') as csvfile:
        csv_reader = csv.reader(csvfile)

        for i in range(0, 162):
            # NOTE: no freq drift
            dphi = np.float32(twopidt * (np.float32(f0) + (np.float32(symbols[i]) - tone_separation) * df))
            for j in range(0, nspersym):
                ii = nspersym * i + j

                row = next(csv_reader)
                wsprd_f0 = np.float32(row[0])
                wsprd_twopidt = np.float32(row[1])
                wsprd_df = np.float32(row[2])
                wsprd_i = int(row[3])
                wsprd_cs = np.float32(row[4])
                wsprd_dphi = np.float32(row[5])
                wsprd_ii = int(row[6])
                wsprd_phi = np.float32(row[7])
                wsprd_drift0 = np.float32(row[8])
                wsprd_component = np.float32(row[9])

#                if wsprd_i != i:
#                    print(f"{callsign} i mismatch {wsprd_i} {i}")
#                if wsprd_ii != ii:
#                    print(f"{callsign} ii mismatch {wsprd_ii} {ii}")
#                if wsprd_f0 != f0:
#                    print(f"{callsign} f0 mismatch {wsprd_f0} {f0}")
#                if wsprd_twopidt != twopidt:
#                    print(f"{callsign} twopidt mismatch {wsprd_twopidt} {twopidt}")
#                if wsprd_df != df:
#                    print(f"{callsign} df mismatch {wsprd_df} {df}")
#                if wsprd_cs != symbols[i]:
#                    print(f"{callsign} cs mimatch {wsprd_cs} {symbols[i]}")
                if wsprd_dphi != dphi:
                    print(f"{callsign} dphi mismatch {wsprd_dphi} {dphi}")
                    print(f"{wsprd_twopidt} * ({wsprd_f0} + ({wsprd_cs} - 1.5) * {wsprd_df})")
                    new_dphi = np.float32(wsprd_twopidt * (f0 + (wsprd_cs - np.float32(1.5)) * wsprd_df))
                    print(new_dphi)
                    print(f"{callsign} {wsprd_dphi} {new_dphi}")
                    print(f"{callsign} {wsprd_drift0} {wsprd_component}")
                #    dphi = wsprd_dphi
#                if wsprd_phi != phi:
#                    print(f"{callsign} phi mismatch {wsprd_phi} {phi}, fixing")
#                    print(f"{callsign} {type(phi)}")
                    #phi = wsprd_phi

                reference[ii] = np.cos(phi) + 1j * np.sin(phi)
                phi += dphi

    with open(f"wsprd/{callsign}_refi_refq.csv", newline='') as csvfile:
        average_diff = 0.0
        for i, row in enumerate(csv.reader(csvfile)):
            wsprd_refi = float(row[0])
            wsprd_refq = float(row[1])
            our_refi = reference[i].real
            our_refq = reference[i].imag

            diff_i = wsprd_refi - our_refi
            diff_q = wsprd_refq - our_refq
            diff = abs(diff_i) + abs(diff_q)
            average_diff += diff
            #if (i % 100) == 0:
            #    print(f"{callsign} reference {i} {diff_i} {diff_q}")

            #reference[i] = wsprd_refi + 1j * wsprd_refq
        average_diff /= 162 * nspersym
        print(f"{callsign} average_diff {average_diff}")

    return reference

def calc_g(signal, iq_dat, r):
    """"From wsprd.c:
    s(t) * conjugate(r(t))
    beginning of first symbol in reference signal is at i=0
    beginning of first symbol in received data is at shift0.
    filter transient lasts nfilt samples
    leave nfilt zeros as a pad at the beginning of the unfiltered reference signal
    """

    nsym = 162
    nspersym = 256
    nfilt = 360

    g = np.zeros(SAMPLES, dtype=complex)
    shift0 = round((signal['dt'] + 1.0) * SAMPLE_RATE)

    # NOTE: Not enough precision in signal['dt']
    if shift0 == 79:
        shift0 = 80

    for i in range(0, nsym * nspersym):
        k = shift0 + i
        if k > 0 and k < SAMPLES:
            g[i+nfilt] = iq_dat[k] * r[i].conjugate()

    return g


def read_all_wspr(all_wspr_file):
    """Reads the output of wsprd from an ALL_WSPR.txt format file.

    Sample file contents:
    sprd mple  -9  1.13   0.0014463  ND6P DM04 30            0  0.68  1  1    0  0   0     1   612
    sprd mple -15  0.07   0.0014604  W5BIT EL09 17           0  0.48  1  1    0  0   7     1   500
    sprd mple -23  2.16   0.0014652  G8VDQ IO91 37           0  0.22  3  1   16  0   0  8080  -187
    sprd mple  -6  0.62   0.0014891  WD4LHT EL89 30          0  0.66  1  1    0  0   2     1   539
    sprd mple  -1 -0.79   0.0015028  NM7J DM26 30            0  0.33  1  1    0  0  26     1    45
    sprd mple -21  0.54   0.0015174  KI7CI DM09 37           0  0.36  1  1    0  0   3     1   557
    sprd mple -18 -1.94   0.0015298  DJ6OL JO52 37           0  0.35  1  1    0  0  11     1   364
    sprd mple -11  0.83   0.0015868  W3HH EL89 30            0  0.73  1  1    0  0   0     1   734
    sprd mple -25  0.75   0.0015944  W3BI FN20 30            0  0.14  2  1    0  0  30  1223   -90

    (date, time, snr, dt, freq, message, drift, sync, ipass, blocksize, jitter, decodetype, nhardmin, cycles, metric)

    Sample output:
    [
        {
            'date': 'sprd',
            'time': 'mple',
            'snr': '-9',
            'dt': 1.1,
            'freq': -53.7,
            'drift': '0',
            'callsign': 'ND6P',
            'locator': 'DM04',
            'power': '30'
        },
        {
            'date': 'sprd',
            'time': 'mple',
            'snr': '-15',
            'dt': '0.1',
    <snip>
    """

    decodes = []

    with open(all_wspr_file, 'r') as all_wspr:
        for line in all_wspr:
            parsed_line = line.split()
            if len(parsed_line) < 9:
                raise WSPRSpreadError("Unable to parse ALL_WSPR.txt line: " + line)
            decodes.append(
                {
                    'date'      : parsed_line[0],
                    'time'      : parsed_line[1],
                    'snr'       : int(parsed_line[2]), 
                    'dt'        : float(parsed_line[3]),
                    'freq'      : float(parsed_line[4]) * 1e6 - 1500,
                    'callsign'  : parsed_line[5], 
                    'locator'   : parsed_line[6],
                    'power'     : parsed_line[7],
                    'drift'     : int(parsed_line[8]),
                    'sync'      : float(parsed_line[9]),
                    'ipass'     : int(parsed_line[10]),
                    'blocksize' : int(parsed_line[11]),
                    'jitter'    : int(parsed_line[12]),
                    'decodetype': int(parsed_line[13]),
                    'nhardmin'  : int(parsed_line[14]),
                    'cycles'    : int(parsed_line[15]),
                    'metric'    : int(parsed_line[16]),
                }
            )
    return decodes

parser = argparse.ArgumentParser(description="Calculates the Doppler spread (w50) of signals decoded by wsprd")
parser.add_argument("wav_file", help="WAV file of WSPR audio (same format was wsprd)")
parser.add_argument("all_wspr", help="ALL_WSPR.txt format file with wsprd output (defaults to ALL_WSPR.txt)", default='ALL_WSPR.txt')
parser.add_argument("-d", "--debug", help="Print debugging messages and create FFT graphs", action='store_true')
args = parser.parse_args()

decodes = read_all_wspr(args.all_wspr)
iq_dat = read_wav_file(args.wav_file)

with open("wspr_spread_iq_dat.csv", "w") as f:
    for point in iq_dat:
        f.write(f"{point.real},{point.imag}\n")

if args.debug:
    freqs = np.fft.fftfreq(SAMPLES, 1 / SAMPLE_RATE)
    fft = np.square(np.abs(np.fft.fft(iq_dat)))
    plot_fft('iq_dat_fft.png', "IQ Data FFT (wspr_spread.py)", freqs, fft)

for signal in decodes:
    if args.debug:
        print(signal)

    ref_symbols = wspr.encode(signal['callsign'], signal['locator'], signal['power'])
    r = reference(signal['freq'], ref_symbols, signal['callsign'])

    with open(f"{signal['callsign']}_ref_wspr_spread.csv", "w") as f:
        for point in r:
           f.write(f"{point.real},{point.imag}\n")

    if args.debug:
        freqs = np.fft.fftfreq(SAMPLES, 1 / SAMPLE_RATE)
        fft = np.square(np.abs(np.fft.fft(r)))
        plot_fft(f"{signal['callsign']}_r_fft.png", f"{signal['callsign']} R FFT (wspr_spread.py)", freqs, fft)

#    print("Substituting wsprd's id and qd for iq_dat") 
#    with open(f"wsprd/{signal['callsign']}_id_qd.csv", newline='') as csvfile:
#        for i, row in enumerate(csv.reader(csvfile)):
#            idat = float(row[0])
#            qdat = float(row[1])
#            iq_dat[i] = idat + 1j * qdat

#    print("Substituting wsprd's refi and refq for r")
#    with open(f"wsprd/{signal['callsign']}_refi_refq.csv", newline='') as csvfile:
#        for i, row in enumerate(csv.reader(csvfile)):
#            refi = float(row[0])
#            refq = float(row[1])
#            r[i] = refi + 1j * refq

    g = calc_g(signal, iq_dat, r)

#    print("Substituting wsprd's c(t) for g")
#    with open(f"wsprd/{signal['callsign']}_ci_cq.csv", newline='') as csvfile:
#        for i, row in enumerate(csv.reader(csvfile)):
#            ci = float(row[0])
#            cq = float(row[1])
#            g[i] = ci + 1j* cq

    g_fft = np.fft.fftshift(np.fft.fft(g)) # note the shift from "in-order"

    if args.debug:
        freqs = np.arange(-1 * (BW / 2), (BW / 2),  DF)
        plot_fft(f"{signal['callsign']}_g_fft.png", f"{signal['callsign']} G FFT (wspr_spread.py)", freqs, np.square(np.abs(g_fft)))

    # slice out the region of interest
    # and calculate power: real_component^2 + imag_component^2
    region_of_interest = np.square(np.abs(g_fft[REGION_START:REGION_END]))

    # slice the power region out of the region of interest
    power_region = region_of_interest[POWER_START:POWER_END]

    # scale the whole region of interest by the max power found in the power region
    region_of_interest /= power_region.max()

    # slice out the noise regions from the region of interest
    left_noise_region = region_of_interest[LEFT_NOISE_START:LEFT_NOISE_END]
    right_noise_region = region_of_interest[RIGHT_NOISE_START:RIGHT_NOISE_END]

    # get the average in these regions and pick the lowest to be the noise factor
    # NOTE: this might be a good place to add a warning if the noise is high
    left_noise = np.average(left_noise_region)
    right_noise = np.average(right_noise_region)
    if left_noise < right_noise:
        noise = left_noise
    else:
        noise = right_noise
    if args.debug:
        print(f"noise={noise}")

    # subtract out the noise out of the region of interest
    region_of_interest -= noise

    if args.debug:
        freqs = np.arange(-4, 4,  DF)
        plot_fft(f"{signal['callsign']}_region_of_interest.png", f"{signal['callsign']} Region of Interest (wspr_spread.py)", freqs, region_of_interest)

    # find the index where there's 25% of the total power in the power region
    # and the index where there's 75% of the total power.
    # That frequency range is the doppler spread (w50)
    power_region = region_of_interest[POWER_START:POWER_END]
    total_power = np.sum(power_region)
    curr_sum = 0
    index25 = None
    index75 = None
    for i, value in enumerate(power_region):
        curr_sum += value
        if curr_sum >= 0.25 * total_power and not index25:
            index25 = i
        elif curr_sum >= 0.75 * total_power and not index75:
            index75 = i
    signal['doppler_spread'] = (index75 - index25) * DF

    # here's the C format string from wsprd for reference:
    # "%6s %4s %3.0f %5.2f %11.7f  %-22s %2d %5.2f %2d %2d %4d %2d %3d %5u %5d\n"
    freq = (signal['freq'] + 1500) / 1e6
    message = f"{signal['callsign']} {signal['locator']} {signal['power']}"
    print(f"{signal['date']:>6} " + # %6s
          f"{signal['time']:>4} " + # %4s
          f"{signal['snr']:>3} " + # %3.0f
          f"{signal['dt']:>5.2f} " + # %5.2f
          f"{freq:>11.7f}  " + # %11.7f
          f"{message:<22} " + # %-22s
          f"{signal['drift']:2} " + # %2d
          f"{signal['sync']:5.2f} " + # %5.2f
          f"{signal['ipass']:2} " + # %2d
          f"{signal['blocksize']:2} " + # %d2
          f"{signal['jitter']:4} " + # %4d
          f"{signal['decodetype']:2} " + # %2d
          f"{signal['nhardmin']:3} " + # %3d
          f"{signal['cycles']:5} " + # %5u
          f"{signal['metric']:5} " + # %5d
          f"{signal['doppler_spread']:>6.3f}") # ???
