import matplotlib.pyplot as plt
import numpy as np
import csv

#read   3 -0.8   0.001494  0  W5EEX DM33 37
#read   3 -1.6   0.001502  0  WB6CXC CM88 30

SAMPLES = 46080
SAMPLE_RATE = 375
# divide a freq by DF to get an index
DF = SAMPLE_RATE / SAMPLES
BW = 375

def plot_fft(filename, title, freqs, fft):
    """Utility function to plot things"""

    fig = plt.figure()
    plt.title(title)
    plt.xlabel('Freq (Hz)')
    plt.ylabel('Magnitude')
    plt.plot(freqs, fft)
    plt.savefig(filename)
    plt.close(fig)

iq_dat = np.zeros(46080, dtype=complex)
with open("idat_qdat.csv", newline='') as csvfile:
    for i, row in enumerate(csv.reader(csvfile)):
        iq_dat[i] = complex(float(row[0]), float(row[1]))

fft = np.square(np.abs(np.fft.fft(iq_dat)))
fft = np.fft.fftshift(fft)
freqs = np.arange(1500 - BW / 2, 1500 + BW / 2, BW / SAMPLES)

start_freq = -10
end_freq = 6
start_index = int((SAMPLES / 2) + (start_freq / DF))
end_index = int((SAMPLES / 2) + (end_freq / DF))

freqs = freqs[start_index:end_index]
fft = fft[start_index:end_index]
plot_fft('iq_dat_fft.png', "IQ Data FFT (analysis.py)", freqs, fft)

freqs = []
power = []
with open("W5EEX_roi.csv", newline='') as csvfile:
    for row in csv.reader(csvfile):
        freqs.append(float(row[0]))
        power.append(float(row[1]))
fig = plt.figure()
plt.title("W5EEX Region of Interest (analysis.py)")
plt.xlabel("Freq (Hz)")
plt.ylabel("Power")
plt.plot(freqs, power)
plt.savefig('W5EEX_roi.png')
plt.close(fig)

freqs = []
power = []
with open("WB6CXC_roi.csv", newline='') as csvfile:
    for row in csv.reader(csvfile):
        freqs.append(float(row[0]))
        power.append(float(row[1]))
fig = plt.figure()
plt.title("WB6CXC Region of Interest (analysis.py)")
plt.xlabel("Freq (Hz)")
plt.ylabel("Power")
plt.plot(freqs, power)
plt.savefig('WB6CXC_roi.png')
plt.close(fig)

