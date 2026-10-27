"""Mobile Communications Mini Project 2.2 - OFDM with cyclic prefix over a multipath channel.

Chain: bits -> Gray QAM -> N=64 subcarriers -> IFFT -> add cyclic prefix -> multipath FIR channel
       (+ AWGN) -> remove CP -> FFT -> one-tap equalizer (known channel) -> decisions.
The same channel realization is used for the "sufficient CP", "no CP" and "too-short CP" receivers.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mc_common import modulate, demodulate, bits_per_symbol

rng = np.random.default_rng(1)
N = 64                                                    # subcarriers / FFT size

# ---- synthetic multipath channel: 3 taps, delays in samples, complex gains, unit total power
DELAYS = [0, 3, 8]
GAINS = np.array([1.0, 0.5 * np.exp(1j * 0.8), 0.3 * np.exp(-1j * 1.9)])   # |H(k)| >= 1 - 0.5 - 0.3 (before norm.)
GAINS = GAINS / np.sqrt(np.sum(np.abs(GAINS) ** 2))
TAU = max(DELAYS)                                         # delay spread in samples -> CP must be >= TAU
h = np.zeros(TAU + 1, complex); h[DELAYS] = GAINS
H = np.fft.fft(h, N)                                      # channel frequency response on the N subcarriers


def ofdm_modulate(X):
    """X: (nsym, N) frequency-domain symbols -> time-domain, unit average sample power."""
    return np.fft.ifft(X, axis=1) * np.sqrt(N)


def ofdm_demodulate(x):
    return np.fft.fft(x, axis=1) / np.sqrt(N)


def add_cp(x, P):
    return np.hstack([x[:, N - P:], x]) if P > 0 else x   # copy of the last P samples in front


def strip_cp(y, P):
    return y[:, P:]


def run_chain(P, M, snr_db=None, nsym=300, bits=None):
    """Return (X_tx, X_hat, bits_tx, bits_rx). snr_db=None -> noiseless."""
    k = bits_per_symbol(M)
    if bits is None:
        bits = rng.integers(0, 2, size=nsym * N * k)
    X = modulate(bits, M).reshape(-1, N)
    tx = add_cp(ofdm_modulate(X), P)
    stream = tx.reshape(-1)                                # back-to-back symbols -> ISI between them
    y = np.convolve(stream, h)[:len(stream)]
    if snr_db is not None:
        s2 = 10 ** (-snr_db / 10)
        y = y + np.sqrt(s2 / 2) * (rng.standard_normal(y.shape) + 1j * rng.standard_normal(y.shape))
    Y = ofdm_demodulate(strip_cp(y.reshape(-1, N + P), P))
    Xhat = Y / H                                           # one-tap zero-forcing equalizer per subcarrier
    bits_rx = demodulate(Xhat.reshape(-1), M)
    return X, Xhat, bits, bits_rx


def metrics(X, Xhat, bits, bits_rx, M):
    """Skip the first symbol (it has no predecessor, hence no ISI)."""
    k = bits_per_symbol(M)
    X, Xhat = X[1:], Xhat[1:]
    b, br = bits.reshape(-1, N * k)[1:], bits_rx.reshape(-1, N * k)[1:]
    mismatch = (b != br).reshape(-1, k)
    return dict(SER=float(np.mean(mismatch.any(axis=1))), BER=float(np.mean(mismatch)),
                MSE=float(np.mean(np.abs(Xhat - X) ** 2)))


# =========================================================== self-checks
print("=== Self-checks ===")
for M in (4, 16):
    tb = rng.integers(0, 2, size=bits_per_symbol(M) * 5000)
    s = modulate(tb, M)
    assert np.array_equal(demodulate(s, M), tb) and abs(np.mean(np.abs(s) ** 2) - 1) < 0.03
print("  mapper: QPSK and 16-QAM round trip exact, unit average energy: OK")
Xt = modulate(rng.integers(0, 2, size=4 * N * 4), 16).reshape(-1, N)
err = np.max(np.abs(ofdm_demodulate(ofdm_modulate(Xt)) - Xt))
assert err < 1e-12; print(f"  noiseless IFFT -> FFT max error = {err:.2e}: OK")
xs = ofdm_modulate(Xt); P0 = 16
assert np.allclose(add_cp(xs, P0)[:, :P0], xs[:, N - P0:]) and np.allclose(add_cp(xs, P0)[:, P0:], xs)
print(f"  CP of {P0} samples equals the tail of the symbol: OK")
for M in (4, 16):
    X, Xh, b, br = run_chain(TAU, M, snr_db=None)
    m_ok = metrics(X, Xh, b, br, M)
    X, Xh, b, br = run_chain(0, M, snr_db=None)
    m_bad = metrics(X, Xh, b, br, M)
    assert m_ok["MSE"] < 1e-20 and m_ok["SER"] == 0
    print(f"  M={M:2d} noiseless: CP={TAU} (= delay spread) MSE={m_ok['MSE']:.1e}, SER={m_ok['SER']:.3f} | "
          f"CP=0 MSE={m_bad['MSE']:.3e}, SER={m_bad['SER']:.3f}")
print(f"  channel: delays {DELAYS} samples, |gains| {np.round(np.abs(GAINS), 3)}, delay spread = {TAU} samples; "
      f"deepest subcarrier fade |H| = {np.abs(H).min():.3f}")

# =========================================================== constellation comparison (16-QAM, same channel)
SNR = 30.0
M = 16
bits = rng.integers(0, 2, size=400 * N * bits_per_symbol(M))        # identical data for the three receivers
cases = [(16, f"CP = 16 (>= {TAU}): sufficient"), (4, f"CP = 4 (< {TAU}): too short"), (0, "CP = 0: no prefix")]
fig, axs = plt.subplots(1, 3, figsize=(13, 4.4))
print(f"\n=== Constellation comparison, 16-QAM, SNR = {SNR:g} dB, identical channel and data ===")
for ax, (P, title) in zip(axs, cases):
    X, Xh, b, br = run_chain(P, M, SNR, bits=bits)
    m = metrics(X, Xh, b, br, M)
    pts = Xh[1:].reshape(-1)
    ax.plot(pts.real, pts.imag, ".", ms=1.5, alpha=0.5)
    ref = np.unique(np.round(X[1:].reshape(-1), 6))
    ax.plot(ref.real, ref.imag, "r+", ms=8)
    ax.set_xlim(-1.8, 1.8); ax.set_ylim(-1.8, 1.8); ax.set_aspect("equal"); ax.grid(alpha=0.3)
    ax.set_title(f"{title}\nSER={m['SER']:.3f}, MSE={10*np.log10(m['MSE']):.1f} dB", fontsize=9)
    print(f"  CP={P:2d}: SER={m['SER']:.4f}  BER={m['BER']:.4f}  MSE={m['MSE']:.3e} ({10*np.log10(m['MSE']):.1f} dB)")
fig.suptitle("Equalized subcarrier symbols (red + = transmitted points)")
fig.tight_layout(); fig.savefig("fig_constellations.png", dpi=150); plt.close(fig)

# =========================================================== error vs CP length sweep
SWEEP_SNR = 30.0
cps = np.arange(0, 21)
res = {}
for Mx in (4, 16):
    k = bits_per_symbol(Mx)
    data = rng.integers(0, 2, size=1200 * N * k)
    res[Mx] = [metrics(*run_chain(P, Mx, SWEEP_SNR, bits=data), Mx) for P in cps]

fig, axs = plt.subplots(1, 2, figsize=(12, 4.2))
for Mx, st in ((4, "o-"), (16, "s-")):
    ser = np.array([r["SER"] for r in res[Mx]]); floor = 0.5 / (1199 * N)
    axs[0].semilogy(cps, np.maximum(ser, floor), st, ms=4, label=f"{'QPSK' if Mx == 4 else '16-QAM'}")
    axs[1].plot(cps, [10 * np.log10(r["MSE"]) for r in res[Mx]], st, ms=4, label=f"{'QPSK' if Mx == 4 else '16-QAM'}")
for a in axs:
    a.axvline(TAU, color="k", ls="--", lw=1, label=f"delay spread = {TAU} samples")
    a.set_xlabel("cyclic prefix length (samples)"); a.grid(True, which="both", alpha=0.3); a.legend(fontsize=8)
axs[0].set_ylabel("symbol error rate (zeros drawn at floor)"); axs[1].set_ylabel("MSE of equalized symbols (dB)")
axs[0].set_title(f"SER vs CP length (SNR = {SWEEP_SNR:g} dB)"); axs[1].set_title("MSE vs CP length")
fig.tight_layout(); fig.savefig("fig_cp_sweep.png", dpi=150); plt.close(fig)

print(f"\n=== Sweep at SNR = {SWEEP_SNR:g} dB: SER and MSE(dB) vs CP length ===")
print(f"  {'CP':>3} | {'QPSK SER':>10} {'MSE dB':>8} | {'16QAM SER':>10} {'MSE dB':>8}")
for i, P in enumerate(cps):
    q, s = res[4][i], res[16][i]
    print(f"  {P:>3} | {q['SER']:>10.4f} {10*np.log10(q['MSE']):>8.1f} | {s['SER']:>10.4f} {10*np.log10(s['MSE']):>8.1f}")

# =========================================================== time-domain view of CP and ISI
P = 16
X = modulate(rng.integers(0, 2, size=3 * N * 2), 4).reshape(3, N)
tx = add_cp(ofdm_modulate(X), P); stream = tx.reshape(-1); y = np.convolve(stream, h)[:len(stream)]
L = N + P
fig, axs = plt.subplots(2, 1, figsize=(11, 5), sharex=True)
for ax, sig, name in ((axs[0], stream, "transmitted"), (axs[1], y, "after multipath channel")):
    ax.plot(np.abs(sig), lw=0.9)
    for j in range(3):
        ax.axvspan(j * L, j * L + P, color="orange", alpha=0.25, label="cyclic prefix" if j == 0 else None)
        ax.axvline(j * L, color="k", lw=0.6)
    ax.set_ylabel(f"|x|, {name}"); ax.legend(loc="upper right", fontsize=8)
axs[1].set_xlabel("sample index (three consecutive OFDM symbols, CP shaded)")
fig.tight_layout(); fig.savefig("fig_waveform.png", dpi=150); plt.close(fig)
print("\nDone: fig_constellations.png fig_cp_sweep.png fig_waveform.png")
