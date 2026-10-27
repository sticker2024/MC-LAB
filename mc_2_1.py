"""Mobile Communications Mini Project 2.1
Part (a): Friis free-space path loss.  Part (b): Monte Carlo BER of BPSK / QPSK / 16-QAM / 64-QAM
over AWGN, compared with theory and with the Shannon-Hartley limit.
Extras (lab tasks 1-6): QPSK, rate-1/2 convolutional code + Viterbi, 64-QAM, log-distance model,
log-normal shadowing, throughput proxy vs Shannon capacity.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.special import erfc
from mc_common import modulate, demodulate, awgn, bits_per_symbol

np.random.seed(0)                       # reproducible bitstreams (lab step 1)
rng = np.random.default_rng(0)          # reproducible noise


# ============================================================ Part (a): Friis
def friis_rx_power(Ptx_dBm, Gtx_dBi, Grx_dBi, d_km, f_MHz):
    """P_RX(dBm) = PTX + GTX + GRX - 20log10(d_km) - 20log10(f_MHz) - 32.44 (element-wise)."""
    return (Ptx_dBm + Gtx_dBi + Grx_dBi
            - 20 * np.log10(d_km) - 20 * np.log10(f_MHz) - 32.44)


PTX, GTX, GRX = 43.0, 15.0, 0.0
d_km = np.linspace(0.05, 10, 500)
freqs = [900, 1900, 3500]

fig, ax = plt.subplots(figsize=(7, 4))
for f in freqs:
    ax.plot(d_km, friis_rx_power(PTX, GTX, GRX, d_km, f), label=f"{f} MHz")
ax.set_xscale("log"); ax.set_xlabel("distance (km, log scale)"); ax.set_ylabel("received power (dBm)")
ax.set_title("Friis free-space received power (PTX=43 dBm, GTX=15 dBi, GRX=0 dBi)")
ax.grid(True, which="both", alpha=0.3); ax.legend()
fig.tight_layout(); fig.savefig("fig_friis.png", dpi=150); plt.close(fig)

print("=== Friis 6 dB-per-doubling check ===")
for f in freqs:
    for d1 in (0.5, 1.0, 4.0):
        diff = friis_rx_power(PTX, GTX, GRX, d1, f) - friis_rx_power(PTX, GTX, GRX, 2 * d1, f)
        print(f"  f={f:4d} MHz: P({d1} km) - P({2*d1} km) = {diff:.4f} dB   (6.0206 expected)")
for d1 in (1.0, 3.0):
    diff = friis_rx_power(PTX, GTX, GRX, d1, 900) - friis_rx_power(PTX, GTX, GRX, d1, 1800)
    print(f"  d={d1} km: P(900 MHz) - P(1800 MHz) = {diff:.4f} dB   (doubling frequency)")


# ============================================================ Part (b): BER
def ber_theory(M, ebn0_lin):
    """BPSK/QPSK: 0.5*erfc(sqrt(Eb/N0)). Square Gray M-QAM: standard high-SNR approximation."""
    if M in (2, 4):
        return 0.5 * erfc(np.sqrt(ebn0_lin))
    k = bits_per_symbol(M)
    return (2 / k) * (1 - 1 / np.sqrt(M)) * erfc(np.sqrt(3 * k * ebn0_lin / (2 * (M - 1))))


N_BITS = 2_400_000                       # divisible by 1, 2, 4, 6 (so every scheme gets whole symbols)
EBN0_DB = np.arange(0, 17, 1.0)
SCHEMES = {"BPSK": 2, "QPSK": 4, "16-QAM": 16, "64-QAM": 64}

# unit tests of the mapper (lab step 7)
for name, M in SCHEMES.items():
    tb = np.random.randint(0, 2, size=6 * 1000)
    s = modulate(tb, M)
    assert np.array_equal(demodulate(s, M), tb), f"{name}: noiseless round trip failed"
    assert abs(np.mean(np.abs(s) ** 2) - 1) < 0.02, f"{name}: average energy is not 1"

sim = {}
for name, M in SCHEMES.items():
    bits = np.random.randint(0, 2, size=N_BITS)
    sym = modulate(bits, M)
    ber = []
    for e in EBN0_DB:
        rx = demodulate(awgn(sym, e, M, rng), M)
        ber.append(np.mean(rx != bits))
    sim[name] = np.array(ber)


def plot_ber(names, path, title, with_theory_all):
    fig, ax = plt.subplots(figsize=(7, 4.6))
    x = np.linspace(0, 16, 400)
    for name in names:
        M = SCHEMES[name]
        b = sim[name]; ok = b > 0
        ax.semilogy(EBN0_DB[ok], b[ok], "o", ms=4, label=f"{name} simulated")
        if with_theory_all or name == "BPSK":
            ax.semilogy(x, ber_theory(M, 10 ** (x / 10)), "-", lw=1, alpha=0.8, label=f"{name} theory")
    ax.set_ylim(1e-6, 1); ax.set_xlim(0, 16)
    ax.set_xlabel("Eb/N0 (dB)"); ax.set_ylabel("bit error rate"); ax.set_title(title)
    ax.grid(True, which="both", alpha=0.3); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)


plot_ber(["BPSK", "16-QAM"], "fig_ber_bpsk_16qam.png", "BER over AWGN: BPSK vs 16-QAM (BPSK theory overlaid)", False)
plot_ber(list(SCHEMES), "fig_ber_all.png", "BER over AWGN: BPSK, QPSK, 16-QAM, 64-QAM", True)

print("\n=== BER sanity: BPSK simulated vs theory ===")
for e in (0, 4, 6, 8, 9):
    i = int(np.where(EBN0_DB == e)[0][0])
    print(f"  Eb/N0={e:2d} dB: simulated {sim['BPSK'][i]:.3e}   theory {ber_theory(2, 10**(e/10)):.3e}")


def ebn0_at(ber_curve, target):
    """Eb/N0 (dB) where the simulated curve crosses `target`, log-interpolated."""
    ok = ber_curve > 0
    x, y = EBN0_DB[ok], np.log10(ber_curve[ok])
    return float(np.interp(np.log10(target), y[::-1], x[::-1]))


TARGET = 1e-3
print(f"\n=== Eb/N0 needed for BER = {TARGET:g} (simulated) ===")
need = {n: ebn0_at(sim[n], TARGET) for n in SCHEMES}
print(f"  {'scheme':<8}{'k':>3}{'Eb/N0 dB':>10}{'Es/N0 dB':>10}{'Shannon min Eb/N0 dB':>24}{'distance to limit dB':>22}")
for n, M in SCHEMES.items():
    k = bits_per_symbol(M); eta = k                                  # bits/s/Hz, Nyquist signalling
    shannon = 10 * np.log10((2 ** eta - 1) / eta)
    print(f"  {n:<8}{k:>3}{need[n]:>10.2f}{need[n] + 10*np.log10(k):>10.2f}{shannon:>24.2f}{need[n] - shannon:>22.2f}")
print(f"  Eb/N0 gap 16-QAM - BPSK = {need['16-QAM'] - need['BPSK']:.2f} dB")
print(f"  Es/N0 (SNR) gap 16-QAM - BPSK = {need['16-QAM'] + 10*np.log10(4) - need['BPSK']:.2f} dB")
print(f"  QPSK - BPSK Eb/N0 difference = {need['QPSK'] - need['BPSK']:.2f} dB (should be ~0)")
print(f"  64-QAM - 16-QAM Eb/N0 gap = {need['64-QAM'] - need['16-QAM']:.2f} dB")

# ------------------------------------------------ throughput proxy vs Shannon (task 6)
fig, ax = plt.subplots(figsize=(7, 4.4))
snr_db = np.linspace(-2, 28, 300)
ax.plot(snr_db, np.log2(1 + 10 ** (snr_db / 10)), "k--", label="Shannon C/B = log2(1+SNR)")
for n, M in SCHEMES.items():
    k = bits_per_symbol(M)
    ax.plot(EBN0_DB + 10 * np.log10(k), k * (1 - sim[n]), "o-", ms=3, label=f"{n}: k(1-BER)")
ax.set_xlabel("SNR = Es/N0 (dB)"); ax.set_ylabel("bits per symbol (proxy for b/s/Hz)")
ax.set_title("Throughput proxy vs Shannon capacity"); ax.grid(True, alpha=0.3); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig("fig_throughput.png", dpi=150); plt.close(fig)

# ------------------------------------------------ rate-1/2 convolutional code + Viterbi (task 2)
NXT = np.zeros((4, 2), int); OUT = np.zeros((4, 2, 2), int)
for s in range(4):
    u1, u2 = (s >> 1) & 1, s & 1                       # state = (u[t-1], u[t-2])
    for u in (0, 1):
        OUT[s, u] = [u ^ u1 ^ u2, u ^ u2]              # generators 7, 5 (octal), K = 3
        NXT[s, u] = (u << 1) | u1


def conv_encode(u):
    nb, L = u.shape
    s = np.zeros(nb, int); coded = np.empty((nb, L, 2), int)
    for t in range(L):
        coded[:, t, :] = OUT[s, u[:, t]]; s = NXT[s, u[:, t]]
    return coded.reshape(nb, 2 * L)


def viterbi(rx):
    nb = rx.shape[0]; L = rx.shape[1] // 2; r = rx.reshape(nb, L, 2)
    pm = np.full((nb, 4), 1e9); pm[:, 0] = 0
    prev = np.zeros((L, nb, 4), np.int8); inp = np.zeros((L, nb, 4), np.int8)
    for t in range(L):
        new = np.full((nb, 4), 1e9); ps = np.zeros((nb, 4), np.int8); pu = np.zeros((nb, 4), np.int8)
        for s in range(4):
            for u in (0, 1):
                ns = NXT[s, u]
                bm = (r[:, t, 0] != OUT[s, u, 0]).astype(float) + (r[:, t, 1] != OUT[s, u, 1])
                cand = pm[:, s] + bm
                better = cand < new[:, ns]
                new[better, ns] = cand[better]; ps[better, ns] = s; pu[better, ns] = u
        pm = new; prev[t] = ps; inp[t] = pu
    s = np.zeros(nb, int); dec = np.empty((nb, L), int); idx = np.arange(nb)
    for t in range(L - 1, -1, -1):                      # tail bits force the end state to 0
        dec[:, t] = inp[t][idx, s]; s = prev[t][idx, s]
    return dec


NB, LB = 2000, 500                                      # 2000 blocks x 500 info bits = 1e6 info bits
info = np.random.randint(0, 2, size=(NB, LB))
u = np.hstack([info, np.zeros((NB, 2), int)])           # 2 zero tail bits
cb = conv_encode(u)
assert np.array_equal(viterbi(cb)[:, :LB], info), "noiseless Viterbi must be exact"
cbits = cb.reshape(-1)
cs = modulate(cbits, 2)
ber_coded = []
for e in EBN0_DB:
    rx = demodulate(awgn(cs, e, 2, rng, code_rate=0.5), 2).reshape(NB, -1)
    ber_coded.append(np.mean(viterbi(rx)[:, :LB] != info))
ber_coded = np.array(ber_coded)

fig, ax = plt.subplots(figsize=(7, 4.4))
x = np.linspace(0, 16, 400)
ax.semilogy(EBN0_DB[sim["BPSK"] > 0], sim["BPSK"][sim["BPSK"] > 0], "o", ms=4, label="uncoded BPSK (sim)")
ax.semilogy(x, ber_theory(2, 10 ** (x / 10)), "-", lw=1, label="uncoded BPSK (theory)")
okc = ber_coded > 0
ax.semilogy(EBN0_DB[okc], ber_coded[okc], "s-", ms=4, label="rate-1/2, K=3 conv + hard Viterbi")
ax.set_ylim(1e-6, 1); ax.set_xlabel("Eb/N0 (dB)"); ax.set_ylabel("BER"); ax.grid(True, which="both", alpha=0.3)
ax.set_title("Coding gain of rate-1/2 convolutional code"); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig("fig_coded.png", dpi=150); plt.close(fig)
print("\n=== Coded vs uncoded BPSK ===")
for tgt in (1e-3, 1e-4, 1e-5):
    try:
        ec, eu = ebn0_at(ber_coded, tgt), ebn0_at(sim["BPSK"], tgt)
        print(f"  BER {tgt:g}: uncoded {eu:.2f} dB, coded {ec:.2f} dB, coding gain {eu - ec:+.2f} dB")
    except Exception as ex:
        print(f"  BER {tgt:g}: not resolved ({ex})")
print("  coded BER by Eb/N0:", " ".join(f"{e:g}dB:{b:.1e}" for e, b in zip(EBN0_DB[::2], ber_coded[::2])))

# ------------------------------------------------ log-distance + shadowing (tasks 4, 5)
f0, d0 = 1900.0, 0.1                                    # reference 100 m, 1900 MHz
PL_d0 = -(friis_rx_power(0.0, 0.0, 0.0, d0, f0))        # free-space path loss at d0 (dB)
d2 = np.logspace(np.log10(0.1), 1, 300)
fig, axs = plt.subplots(1, 2, figsize=(12, 4.2), sharey=True)
for n in (2, 2.5, 3, 3.5, 4):
    PL = PL_d0 + 10 * n * np.log10(d2 / d0)
    axs[0].plot(d2, PTX + GTX + GRX - PL, label=f"n = {n}")
axs[0].set_title("Log-distance model, exponent n = 2 ... 4")
n_sh, sigma = 3.5, 6.0
base = PTX + GTX + GRX - (PL_d0 + 10 * n_sh * np.log10(d2 / d0))
axs[1].plot(d2, base, "k", lw=2, label=f"mean (n={n_sh})")
for i in range(6):
    axs[1].plot(d2, base + rng.normal(0, sigma, size=d2.shape), lw=0.7, alpha=0.8)
axs[1].set_title(f"Add log-normal shadowing, sigma = {sigma:g} dB (6 realizations)")
for a in axs:
    a.set_xscale("log"); a.set_xlabel("distance (km, log scale)"); a.grid(True, which="both", alpha=0.3); a.legend(fontsize=8)
axs[0].set_ylabel("received power (dBm)")
fig.tight_layout(); fig.savefig("fig_logdist.png", dpi=150); plt.close(fig)
print(f"\nLog-distance reference: PL(d0={d0} km, {f0:g} MHz) = {PL_d0:.2f} dB (free space)")
print("Done: fig_friis.png fig_ber_bpsk_16qam.png fig_ber_all.png fig_throughput.png fig_coded.png fig_logdist.png")
