import numpy as np
import matplotlib.pyplot as plt
params = {
    'axes.labelsize': 16,
    'legend.fontsize': 14,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'text.usetex': False,
    'lines.markeredgewidth': 2,
    'lines.markersize': 7,
    'figure.figsize': [4, 3.6],
    'xtick.direction': 'in',
    'ytick.direction': 'in'
}
plt.rcParams.update(params)

from scipy.fftpack import fft, ifft

lmp_inputs = './in.GK'

kB = 1.38e-23
ps2s = 1e-12
A2m = 1e-10
e2J = 1.6e-19
convert = e2J * e2J / ps2s / A2m
V = 42.432**3
T = 40

f_max = None


def get_timestep():
    with open(lmp_inputs, 'r') as f:
        for line in f:
            if line.startswith("timestep"):
                timestep = float(line.split()[1])

                return timestep


def get_sampling_info():
    data = np.loadtxt("./seed1" + "/flux.dat", skiprows=2)
    nstep = data.shape[0]
    ndump = data[1, 0] - data[0, 0]
    timestep = get_timestep()
    dT = timestep * ndump
    return dT, nstep


def autocorrelation(x):
    xp = x
    f = fft(xp)
    p = np.absolute(f)**2
    pi = ifft(p)
    return np.real(pi) / len(xp)
    #return np.real(pi)[:int(len(xp)/2)]/(len(xp))


def calc_flux_fft(dirname):
    data = np.loadtxt(dirname + "/flux.dat", skiprows=2)

    # Plot the original heat flux signal
    flux = data[:, 1:4]
    #plt.figure("flux")
    #plt.plot(time, flux[:, 0])

    # autocorrection from fft
    plt.figure("flux-acf")

    flux_acf = np.zeros(flux.shape)
    flux_acf[:, 0] = autocorrelation(flux[:, 0])
    flux_acf[:, 1] = autocorrelation(flux[:, 1])
    flux_acf[:, 2] = autocorrelation(flux[:, 2])
    plt.plot(time[:nstep // 2], flux_acf[:nstep // 2, 0])
    plt.plot(time[:nstep // 2], flux_acf[:nstep // 2, 1])
    plt.plot(time[:nstep // 2], flux_acf[:nstep // 2, 2])

    plt.figure('flux-acf-integral')
    ratio = 1.0 / 3.0 * convert / kB / T / T / V * 1 * dT
    flux_acf_sum = np.cumsum(flux_acf, axis=0).sum(axis=1) * ratio
    global acf_sum_ave
    acf_sum_ave += flux_acf_sum
    plt.plot(time, flux_acf_sum)

    plt.figure("flux-fft")
    ratio = 1.0 / 1.0 * convert / kB / T / T / V * 1 * dT / nstep / 2.0
    f_max = 1.0 / dT
    print("F_max: {:.3f}".format(f_max))
    freq = np.linspace(0, 1, nstep) * f_max
    flux_fft = np.fft.fft(flux, axis=0)
    power = np.abs(flux_fft)**2 * ratio
    global flux_power_spectral
    flux_power_spectral += power.mean(axis=1)
    plt.plot(freq, power[:, 0], 'o', ms=1, markerfacecolor='w')
    plt.plot(freq, power[:, 1], 'o', ms=1, markerfacecolor='w')
    plt.plot(freq, power[:, 2], 'o', ms=1, markerfacecolor='w')
    plt.xlim(freq[1], freq_max * 0.5)
    plt.gca().set_xscale("log")

    flux_fft_0 = (flux_fft[1, :])
    tc = (np.abs(flux_fft_0)**2).mean()
    tc = tc * ratio
    print("TC={:.3e}".format(tc))
    return tc


dT, nstep = get_sampling_info()
time = np.linspace(0, nstep, nstep) * dT
freq_max = 1.0 / dT
freq = np.linspace(0, 1, nstep) * freq_max

flux_power_spectral = np.zeros(nstep)
acf_sum_ave = np.zeros(nstep)

tc = np.zeros(10)
for ii in range(10):
    temp = calc_flux_fft('./seed{}'.format(ii + 1))
    tc[ii] = temp

plt.figure('flux-acf-integral')
acf_sum_ave /= 10
plt.plot(time, acf_sum_ave, "r-", lw=2)
plt.xlim(0, time.max()*0.3)
plt.xlabel("Time (ps)")
plt.ylabel("TC (W/mK)")
plt.tight_layout()
plt.savefig("TC_flux_fft.png", dpi=600)

plt.figure("flux-fft")
flux_power_spectral /= 10
plt.plot(freq, flux_power_spectral, 'ro-', ms=5, markerfacecolor='w', lw=3)
plt.xlabel("Frequency (THz)")
plt.ylabel(r"|$J(\omega)^2/2T_f $| [W/mK] ")
plt.tight_layout()
#plt.savefig("Power_spectral.pdf", dpi = 600)
plt.savefig("Power_spectral.png", dpi=600)
plt.show()
print(tc)
print(tc.mean())
