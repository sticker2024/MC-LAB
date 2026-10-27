# Mobile Communications Mini Projects, Unit 2: written deliverables

All numbers below are copied from the scripts' printed output (`output_2_1.txt`, `output_2_2.txt`, `output_2_3.txt`). Re-run the scripts to reproduce them (seeds are fixed).

---

## Mini Project 2.1: Friis path loss and BER over AWGN

**Friis check.** With P_TX = 43 dBm, G_TX = 15 dBi, G_RX = 0 dBi, doubling the distance (at 900, 1900 and 3500 MHz, from 0.5, 1 and 4 km) lowers P_RX by 6.0206 dB every time, and doubling the frequency at a fixed distance does the same. This is the 20·log10(2) term in the equation.

**Measured Eb/N0 needed for BER = 10^-3** (2.4 million bits per scheme, Gray mapping, minimum-distance decisions):

| scheme | bits/symbol | Eb/N0 (dB) | Es/N0 = SNR (dB) | Shannon minimum Eb/N0 at that spectral efficiency (dB) | distance above the limit (dB) |
|---|---|---|---|---|---|
| BPSK | 1 | 6.78 | 6.78 | 0.00 | 6.78 |
| QPSK | 2 | 6.74 | 9.75 | 1.76 | 4.98 |
| 16-QAM | 4 | 10.52 | 16.54 | 5.74 | 4.78 |
| 64-QAM | 6 | 14.73 | 22.51 | 10.21 | 4.52 |

The simulated BPSK curve follows 0.5·erfc(sqrt(Eb/N0)) (for example 3.67e-5 simulated against 3.36e-5 theory at 9 dB), and QPSK sits on top of BPSK (difference 0.04 dB) even though it carries 2 bits per symbol.

**Stated SNR gap.** 16-QAM needs **3.74 dB more Eb/N0** than BPSK at BER 10^-3, which is **9.76 dB more SNR (Es/N0)**.

**Note on the manual's "8-10 dB" statement.** The manual says 16-QAM needs roughly 8-10 dB more Eb/N0 than BPSK. Our measurement and the closed-form BER agree that the Eb/N0 gap is about 4 dB. The 8-10 dB figure matches the *SNR* (Es/N0) gap, which is Eb/N0 plus 10·log10(bits per symbol) = +6.02 dB for 16-QAM. State both numbers in the report and say which one you mean.

**Paragraph for the report (Shannon-Hartley).** C = B·log2(1 + SNR) means that carrying η bits per second per hertz error-free requires SNR >= 2^η - 1. Going from BPSK (η = 1, SNR >= 0 dB) to 16-QAM (η = 4, SNR >= 15 = 11.76 dB) therefore costs at least 11.76 dB of SNR in the limit, and our uncoded schemes need 6.78 dB and 16.54 dB respectively, a measured gap of 9.76 dB. The measured gap is smaller than 11.76 dB only because BPSK is further from the capacity limit (6.78 dB above it) than 16-QAM (4.78 dB above it); neither scheme beats its own limit. The theorem does not predict the BER curve itself. It describes an ideal capacity-achieving code with arbitrarily long blocks, while uncoded modulation with a symbol-by-symbol minimum-distance decision always stays several dB above the ceiling, and a denser constellation (smaller minimum distance for the same average energy) is hit harder by noise. This assumes Nyquist signalling with one complex symbol per hertz.

**Lab tasks completed.**
- QPSK: Eb/N0 within 0.04 dB of BPSK, as expected.
- 64-QAM: 4.21 dB more Eb/N0 than 16-QAM, widening the gap as the log2(1 + SNR) term predicts.
- Rate-1/2, K = 3 (generators 7, 5 octal) convolutional code with hard-decision Viterbi decoding: coding gain +1.06 dB at BER 10^-3, +1.29 dB at 10^-4, +1.67 dB at 10^-5.
- Log-distance model with exponent n = 2 to 4 (reference 78.02 dB at 100 m, 1900 MHz) and log-normal shadowing with sigma = 6 dB (six realizations, drawn independently at each distance point, so no spatial correlation).
- Throughput proxy k(1 - BER) plotted against Shannon capacity log2(1 + SNR): every scheme stays below the curve.

---

## Mini Project 2.2: OFDM with a cyclic prefix

**Set-up.** N = 64 subcarriers, QPSK and 16-QAM. Channel: three taps at delays 0, 3, 8 samples with magnitudes 0.864, 0.432, 0.259 (powers 0.746, 0.187, 0.067, total 1), so the delay spread is 8 samples. Receiver: strip the CP, FFT, one-tap zero-forcing equalizer using the known channel. Symbols are sent back to back so each symbol is smeared into the next.

**Verification.**
- A noiseless IFFT followed by FFT returns the subcarrier values with maximum error 5e-16.
- The CP equals the last P samples of the symbol.
- Noiseless with CP = 8 the MSE is 1e-31 and the SER is 0 for both QPSK and 16-QAM. Noiseless with CP = 0 the MSE is about 0.09 (16-QAM SER 9.7%).

**Results at 30 dB SNR, 16-QAM, same channel and data (`fig_constellations.png`).**

| CP length | SER | MSE |
|---|---|---|
| 16 (>= 8) | 0.0002 | -26.3 dB |
| 4 (< 8) | 0.0216 | -16.5 dB |
| 0 | 0.1054 | -10.0 dB |

**Error versus CP length (`fig_cp_sweep.png`).** The MSE falls in two steps. It drops from -12.9 dB at CP = 2 to -15.7 dB at CP = 3, where the delay-3 tap's echo becomes absorbed, then falls from -21.3 dB at CP = 7 to -26.2 dB at CP = 8, where the delay-8 tap's echo is absorbed. From CP = 8 to 20 it stays at -26.3 dB. QPSK shows no errors in about 77,000 symbols from CP = 7 onward, while 16-QAM settles at a noise-limited SER of roughly 1e-4 to 3e-4 from CP = 8 onward (noise is amplified on the weakest subcarrier, |H| = 0.196, by the equalizer).

**Explanation for the report.** The cyclic prefix does two jobs. As a **guard interval**, it makes the multipath echoes of the previous symbol land inside the first P samples of the next one, and the receiver discards those samples, so no previous-symbol energy reaches the FFT window as long as P is at least the delay spread (8 samples here). As a **cyclic extension**, it is a copy of the symbol's own tail, so the channel's linear convolution with the transmitted block becomes a circular convolution over the N retained samples. The FFT then gives exactly Y_k = H_k·X_k, which keeps the subcarriers orthogonal and lets a single complex division per subcarrier undo the channel (this is why the noiseless MSE is 1e-31). When P is smaller than the delay spread, the last (8 - P) echo samples of the previous symbol leak into the FFT window (inter-symbol interference), and the broken circularity also makes subcarriers interfere with each other. The error therefore falls step by step as each tap's delay is covered and drops sharply once P reaches 8. Beyond that point nothing more can be gained because only noise remains (MSE flat at -26.3 dB), yet each extra sample costs throughput: the useful fraction N/(N + P) is 88.9% at P = 8 but 76.2% at P = 20. The best CP is therefore just above the channel's delay spread.

---

## Mini Project 2.3: WPAN / WLAN / WMAN comparison

### Data table (from `wireless_dataset.csv`)

| technology | class | range (m) | peak rate (Mbps) | tx power | power class |
|---|---|---|---|---|---|
| Bluetooth 5 (BLE) | WPAN | 40 (indoor worst case) to 200 (outdoor line of sight) | 0.125 to 2 (coded PHY to 2M PHY) | +4 dBm (LE Class 2; Class 1 up to +20 dBm) | low |
| ZigBee (802.15.4) | WPAN | 10 (nominal) to 100 (line of sight) | 0.25 at 2.4 GHz | -3 dBm | very low |
| Wi-Fi 6 (802.11ax) | WLAN | under 50 (indoor) to 300 (outdoor) | 4,803 (4-stream HE160 AP) to 9,600 (theoretical) | +28 dBm (5 GHz AP aggregate; EIRP limit 36 dBm) | high (AP draws 26-32 W) |
| WiMAX (802.16e) | WMAN | 2,000 to 5,000 (cell radius); fixed profile up to 10,000 | 15 (mobile peak) to 63 (downlink per sector, 10 MHz); fixed 70 | +30 dBm (base station, experimental link budget) | high |

Derived quantities (peak rate divided by the stated quantity): rate per metre of range 0.0025 (ZigBee), 0.010 (BLE), 32 (Wi-Fi 6), 0.0126 (WiMAX) Mbps/m; rate per mW of transmit power 0.50, 0.80, 15.2, 0.063 Mbps/mW; rate times range 0.025, 0.40, 2,880, 315 Mbps·km. Charts: `fig_tradeoff.png` (range against rate, marker size by transmit power) and `fig_derived.png`.

### Recommendations (checked in `output_2_3.txt`)

**Scenario 1, smart home.** Requirements: R1 longest link about 25 m through walls; R2 40 battery sensors; R3 2 kbps each (0.08 Mbps total); R4 coin-cell, so transmit power at most +10 dBm; R5 at least 2 years of battery life.
**Recommendation: Bluetooth 5 (BLE).** It is the only option that meets every numeric requirement: 40 m indoor range against 25 m needed, up to 2 Mbps against 0.08 Mbps, +4 dBm (about 2.5 mW) against the +10 dBm ceiling. ZigBee has the lowest transmit power (-3 dBm, 0.5 mW) and a published battery life of 100 to 1,000+ days, but its nominal 10 m single-hop range falls short of 25 m, so it would need mesh routers, and 2 years is only reached at the top of its quoted battery range. Wi-Fi 6 (+28 dBm, about 630 mW) and WiMAX (+30 dBm, 1,000 mW) break the power limit by more than 18 dB. Secondary factors: BLE and ZigBee both use unlicensed 2.4 GHz; BLE talks directly to phones without a hub (general knowledge). Caveat on R5: no source gives BLE battery life in days. The one measured comparison in S1 (5-hour test) shows Bluetooth 5 using 31% (indoor) and 41% (outdoor) of its battery while the 802.15.4 node used 100%, but the setups differ, so treat it as indicative only.

**Scenario 2, campus network.** Requirements: R1 0.49 km² (about 700 m by 700 m); R2 5,000 users, 20% active; R3 1,000 users × 2 Mbps = 2,000 Mbps aggregate; R4 mains power; R5 at least 40 m per link indoors.
**Recommendation: Wi-Fi 6.** Tiling the campus with 50 m cells (ideal, no overlap) needs 63 access points. At the AP product rate of 4,803 Mbps, derated to 10% for real shared-medium throughput (an assumption), that is about 30,259 Mbps of aggregate capacity, 15 times the 2,000 Mbps needed. WiMAX covers the whole campus from one base station but offers 15 to 63 Mbps per sector: 63 Mbps shared by 1,000 active users is 63 kbps each, 3.2% of the demand even at peak. BLE (98 cells, about 1.2 Mbps after derating) and ZigBee (10 m range) are not designed for this. Secondary factors: Wi-Fi 6 uses unlicensed spectrum while WiMAX bands (2.3, 2.5, 3.5 GHz) normally require a licence; 63 APs at 26.1 W is about 1.6 kW of PoE power and a cabling job; the S3 study cites Target Wake Time and spatial reuse (up to 30% throughput gain in dense outdoor deployments). Real deployments need more than 63 APs because circular cells never tile without overlap.

**Scenario 3, metro fixed link.** Requirements: R1 one fixed 6 km link, no fibre; R2 at least 20 Mbps sustained; R3 mains or solar power; R4 no mobility.
**Recommendation: WiMAX (802.16 fixed profile).** It is the only technology whose range reaches 6 km: up to 10,000 m for the fixed profile (1.7× margin) at up to 70 Mbps (3.5× the 20 Mbps requirement). The mobile 802.16e cell radius (2 to 5 km) would not reach 6 km, so the fixed configuration with a directional antenna is needed. Wi-Fi 6 offers the rate but only 300 m outdoors (20× too short); BLE (200 m) and ZigBee (100 m) fail on both range and rate. Secondary factors and caveats: WiMAX needs licensed spectrum in most countries; the S7 experiment reached only about 300 m at 32 dBm EIRP and predicted about 2.1 km at 61 dBm EIRP, so range depends strongly on permitted power, antenna height and line of sight, and a 6 km link needs a proper link-budget study; WiMAX has largely been replaced by LTE and 5G in commercial networks (general knowledge, not from the sources), which affects equipment availability; directional high-gain antennas can extend Wi-Fi links well beyond the 300 m figure, but no source in this dataset quantifies that.

### Source notes and limits of the dataset
- S1 (Bluetooth 5 paper), S3 (5G vs Wi-Fi 6 paper) and S7 (WiMAX experiment) are academic papers; S2 is university course notes (about 2007, so ZigBee figures are old); S4 is a vendor product page (transmit power and consumption of one Wi-Fi 6 AP); S5 and S6 are an educational site and lecture notes (secondary sources).
- Sources disagree on WiMAX peak rate (15 Mbps mobile in S5, 63 Mbps downlink per sector in S6, 70 Mbps fixed in S5), so the table keeps the span.
- Subscriber-station transmit power for WiMAX and client-device power for Wi-Fi 6 were not found; the power column reflects access-point/base-station values, and the "power class" labels are the author's assessment from the stated numbers.
- Range depends heavily on environment; each row records its assumption in `wireless_dataset.csv`.

Sources:
- [Bluetooth 5: a concrete step forward towards the IoT (arXiv:1711.00257)](https://ar5iv.arxiv.org/html/1711.00257)
- [IEEE ZigBee / 802.15.4 notes, UNC Charlotte](https://webpages.charlotte.edu/~jmconrad/ECGR6185-2007-01/notes/IEEE_ZigBee.pdf)
- [Revisiting Wireless Internet Connectivity: 5G vs Wi-Fi 6 (arXiv:2010.11601)](https://arxiv.org/pdf/2010.11601)
- [Aruba AP-577 (802.11ax) product page](https://www.securewirelessworks.com/Aruba-AP-577.asp)
- [Understanding IEEE 802.16 WiMAX standards, RF Wireless World](https://www.rfwireless-world.com/articles/understanding-ieee-802-16-wimax-standards)
- [WiMAX lecture notes, McGill University](https://www.info413.ece.mcgill.ca/lecture%20note/E4%20WiMAX.pdf)
- [WiMAX link budget and experiment, IEEE VTC Fall 2008 (Univ. of Bristol)](https://research-information.bris.ac.uk/ws/files/3011100/tran_IEEE_VTC_fall2008.pdf)
