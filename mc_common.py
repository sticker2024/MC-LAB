"""Shared helpers: Gray-coded square QAM / BPSK mapping and min-distance demapping.

All constellations are normalized to unit average symbol energy (Es = 1).
Bit layout for square M-QAM with k = log2(M) bits/symbol: the first k/2 bits select the
in-phase (I) level, the last k/2 bits select the quadrature (Q) level, each Gray coded.
"""
import numpy as np


def _axis_tables(bits_per_axis):
    """Gray-coded PAM levels. Returns (levels_for_gray_value, gray_for_level_index)."""
    L = 2 ** bits_per_axis
    gray_of_idx = np.array([i ^ (i >> 1) for i in range(L)])      # level index -> Gray label
    level_of_gray = np.empty(L)
    for i in range(L):
        level_of_gray[gray_of_idx[i]] = 2 * i - (L - 1)             # ..., -3, -1, +1, +3, ...
    return level_of_gray, gray_of_idx


def qam_norm(M):
    """Scale factor so average energy is 1: Es_raw = 2(M-1)/3 for odd-integer square QAM."""
    return np.sqrt(2 * (M - 1) / 3)


def modulate(bits, M):
    """bits: 1-D array of 0/1. Returns unit-energy complex symbols (BPSK is real)."""
    bits = np.asarray(bits, dtype=np.int64)
    if M == 2:
        return (2.0 * bits - 1.0).astype(complex)                  # bit 0 -> -1, bit 1 -> +1
    k = int(np.log2(M)); h = k // 2
    b = bits.reshape(-1, k)
    w = 1 << np.arange(h - 1, -1, -1)
    gi, gq = b[:, :h] @ w, b[:, h:] @ w
    lv, _ = _axis_tables(h)
    return (lv[gi] + 1j * lv[gq]) / qam_norm(M)


def demodulate(sym, M):
    """Minimum-distance decision (independent per axis for square QAM). Returns bit array."""
    if M == 2:
        return (sym.real > 0).astype(np.int64)
    k = int(np.log2(M)); h = k // 2; L = 2 ** h
    _, gray_of_idx = _axis_tables(h)
    x = sym * qam_norm(M)
    def axis_bits(v):
        idx = np.clip(np.round((v + (L - 1)) / 2), 0, L - 1).astype(np.int64)
        g = gray_of_idx[idx]
        return ((g[:, None] >> np.arange(h - 1, -1, -1)) & 1)
    return np.hstack([axis_bits(x.real), axis_bits(x.imag)]).reshape(-1)


def bits_per_symbol(M):
    return int(np.log2(M))


def awgn(sym, ebn0_db, M, rng, code_rate=1.0):
    """Add complex AWGN for a given Eb/N0 (dB), Eb = ENERGY PER INFORMATION BIT.
    Es = 1 per channel symbol. Uncoded: Eb = Es/k. With a rate-R code each information bit
    spans 1/R coded bits, so Eb = Es/(k*R) and N0 = Eb/(Eb/N0) = 1/(k*R*ebn0)."""
    k = bits_per_symbol(M)
    ebn0 = 10 ** (ebn0_db / 10)
    n0 = 1.0 / (k * code_rate * ebn0)
    sigma = np.sqrt(n0 / 2)                                         # per real dimension
    noise = sigma * (rng.standard_normal(sym.shape) + 1j * rng.standard_normal(sym.shape))
    return sym + noise
