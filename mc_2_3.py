"""Mobile Communications Mini Project 2.3 - WPAN / WLAN / WMAN trade-off tool.

Dataset values come from the sources listed in SOURCES (retrieved 27 Oct 2026). Range figures vary
strongly with environment; every row records the assumption it uses. Where a quantity is not
published in the sources, it is left as NaN rather than guessed.
"""
import math
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SOURCES = {
    "S1": "Bluetooth 5: a concrete step forward towards the IoT, arXiv:1711.00257 - https://ar5iv.arxiv.org/html/1711.00257",
    "S2": "IEEE ZigBee / 802.15.4 course notes, UNC Charlotte - https://webpages.charlotte.edu/~jmconrad/ECGR6185-2007-01/notes/IEEE_ZigBee.pdf",
    "S3": "Revisiting Wireless Internet Connectivity: 5G vs Wi-Fi 6, arXiv:2010.11601 - https://arxiv.org/pdf/2010.11601",
    "S4": "Aruba AP-577 (802.11ax) product page, Secure Wireless Works - https://www.securewirelessworks.com/Aruba-AP-577.asp",
    "S5": "Understanding IEEE 802.16 WiMAX standards, RF Wireless World - https://www.rfwireless-world.com/articles/understanding-ieee-802-16-wimax-standards",
    "S6": "WiMAX lecture notes, McGill University - https://www.info413.ece.mcgill.ca/lecture%20note/E4%20WiMAX.pdf",
    "S7": "WiMAX link-budget/experiment paper, IEEE VTC Fall 2008 (Univ. of Bristol repository) - https://research-information.bris.ac.uk/ws/files/3011100/tran_IEEE_VTC_fall2008.pdf",
}

NAN = float("nan")
DATA = [
    dict(tech="Bluetooth 5 (BLE)", network="WPAN", band_GHz="2.4",
         range_min_m=40, range_max_m=200, rate_low_Mbps=0.125, rate_high_Mbps=2.0,
         tx_dBm=4, power_class="low",
         fixed_range_max_m=NAN, fixed_rate_Mbps=NAN,
         assumption="range: ~40 m indoor worst case to ~200 m outdoor line of sight; rate: coded PHY 125 kbps to 2M PHY; "
                    "tx: LE Class 2 (+4 dBm), Class 1 allows up to +20 dBm; application throughput ~1.4 Mbps",
         sources="S1"),
    dict(tech="ZigBee (IEEE 802.15.4)", network="WPAN", band_GHz="2.4 (also 868/915 MHz)",
         range_min_m=10, range_max_m=100, rate_low_Mbps=0.25, rate_high_Mbps=0.25,
         tx_dBm=-3, power_class="very low",
         fixed_range_max_m=NAN, fixed_rate_Mbps=NAN,
         assumption="range: nominal 10 m per node, up to ~100 m single-hop line of sight; rate: 250 kbps at 2.4 GHz "
                    "(20/40 kbps at 868/915 MHz); tx: -3 dBm; battery life listed as 100-1000+ days",
         sources="S2"),
    dict(tech="Wi-Fi 6 (802.11ax)", network="WLAN", band_GHz="2.4 / 5 / 6",
         range_min_m=50, range_max_m=300, rate_low_Mbps=4803, rate_high_Mbps=9600,
         tx_dBm=28, power_class="high",
         fixed_range_max_m=NAN, fixed_rate_Mbps=NAN,
         assumption="range: <50 m indoors, up to 300 m outdoors; rate: 4.8 Gbps for a 4-stream HE160 AP product up to "
                    "9.6 Gbps theoretical; tx: AP aggregate +28 dBm on 5 GHz (EIRP limit 36 dBm); AP draws 26-32 W",
         sources="S3, S4"),
    dict(tech="WiMAX (IEEE 802.16e)", network="WMAN", band_GHz="2.3 / 2.5 / 3.5",
         range_min_m=2000, range_max_m=5000, rate_low_Mbps=15, rate_high_Mbps=63,
         tx_dBm=30, power_class="high",
         fixed_range_max_m=10000, fixed_rate_Mbps=70,
         assumption="range: 2-5 km typical cell radius (mobile profile; one experiment reached only ~300 m at 32 dBm EIRP); "
                    "rate: 15 Mbps mobile peak (S5) to 63 Mbps downlink per sector at 10 MHz (S6); fixed profile up to 10 km "
                    "and 70 Mbps; tx: 30 dBm base station in the S7 link budget (permitted EIRP up to 61 dBm); "
                    "subscriber-station power not published in the sources",
         sources="S5, S6, S7"),
]
df = pd.DataFrame(DATA)
df["tx_mW"] = 10 ** (df.tx_dBm / 10)
df["Mbps_per_m"] = df.rate_high_Mbps / df.range_max_m                 # throughput per metre of range
df["Mbps_per_mW"] = df.rate_high_Mbps / df.tx_mW                      # rate per unit transmit power
df["Mbps_km"] = df.rate_high_Mbps * df.range_max_m / 1000             # rate x reach "capacity-distance" product
df["cell_area_km2"] = math.pi * (df.range_max_m / 1000) ** 2

pd.set_option("display.width", 200); pd.set_option("display.max_columns", 30)
print("=== Dataset ===")
print(df[["tech", "network", "range_min_m", "range_max_m", "rate_low_Mbps", "rate_high_Mbps",
          "tx_dBm", "power_class"]].to_string(index=False))
print("\n=== Derived quantities ===")
print(df[["tech", "tx_mW", "Mbps_per_m", "Mbps_per_mW", "Mbps_km"]].round(4).to_string(index=False))
df.to_csv("wireless_dataset.csv", index=False)

# ------------------------------------------------------------------ charts
colors = {"WPAN": "#4C78A8", "WLAN": "#F58518", "WMAN": "#54A24B"}
fig, ax = plt.subplots(figsize=(8.5, 6))
for _, r in df.iterrows():
    xc = math.sqrt(r.range_min_m * r.range_max_m); yc = math.sqrt(r.rate_low_Mbps * r.rate_high_Mbps)
    ax.errorbar(xc, yc, xerr=[[xc - r.range_min_m], [r.range_max_m - xc]],
                yerr=[[yc - r.rate_low_Mbps], [r.rate_high_Mbps - yc]],
                fmt="none", ecolor=colors[r.network], alpha=0.6, lw=1.5, capsize=3)
    ax.scatter(xc, yc, s=(r.tx_dBm + 10) * 14, color=colors[r.network], edgecolor="k", zorder=3)
    off = (-40, -34) if r.tech.startswith("ZigBee") else (12, 12)           # keep the two WPAN labels apart
    ax.annotate(f"{r.tech}\n{r.tx_dBm:+.0f} dBm", (xc, yc), textcoords="offset points", xytext=off, fontsize=8)
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("range (m, log) - bar spans indoor/typical to outdoor/maximum")
ax.set_ylabel("data rate (Mbps, log) - bar spans product/practical to peak")
ax.set_title("Range - throughput - power trade-off (marker area grows with transmit power)")
ax.grid(True, which="both", alpha=0.3)
for n, c in colors.items():
    ax.scatter([], [], color=c, edgecolor="k", label=n)
ax.legend(loc="lower right", title="network class")
ax.set_xlim(5, 2e4); ax.set_ylim(0.05, 3e4)
fig.tight_layout(); fig.savefig("fig_tradeoff.png", dpi=150); plt.close(fig)

fig, axs = plt.subplots(1, 3, figsize=(13, 4))
for ax, col, title in zip(axs, ["Mbps_per_m", "Mbps_per_mW", "Mbps_km"],
                          ["peak rate per metre of range (Mbps/m)", "peak rate per mW of transmit power (Mbps/mW)",
                           "rate x range (Mbps*km)"]):
    ax.bar(range(len(df)), df[col], color=[colors[n] for n in df.network])
    ax.set_yscale("log"); ax.set_xticks(range(len(df)))
    ax.set_xticklabels([t.split(" (")[0] for t in df.tech], rotation=25, ha="right", fontsize=8)
    ax.set_title(title, fontsize=9); ax.grid(True, axis="y", which="both", alpha=0.3)
fig.tight_layout(); fig.savefig("fig_derived.png", dpi=150); plt.close(fig)

# ------------------------------------------------------------------ scenarios
def check(row, req, cells_needed=None):
    """Evaluate one technology against a requirement dict. Returns list of (label, ok, detail)."""
    out = []
    r_use = req.get("range_basis", "indoor")
    reach = row.range_min_m if r_use == "indoor" else (
        row.fixed_range_max_m if (r_use == "fixed" and not math.isnan(row.fixed_range_max_m)) else row.range_max_m)
    out.append(("range >= %g m" % req["link_m"], reach >= req["link_m"], f"offers {reach:g} m"))
    rate = row.fixed_rate_Mbps if (r_use == "fixed" and not math.isnan(row.fixed_rate_Mbps)) else row.rate_high_Mbps
    if "aggregate_Mbps" in req:
        # Ideal tiling of the area with circular cells of the indoor/typical radius (optimistic: no overlap),
        # conservative per-cell rate (product-level, not theoretical peak), derated because real shared-medium
        # throughput is only a fraction of the PHY rate (ASSUMPTION: 10%).
        cells = max(1, math.ceil(req["area_m2"] / (math.pi * row.range_min_m ** 2)))
        derate = req.get("derate", 0.10)
        total = cells * row.rate_low_Mbps * derate
        out.append((f"aggregate >= {req['aggregate_Mbps']:g} Mbps", total >= req["aggregate_Mbps"],
                    f"{cells} cells x {row.rate_low_Mbps:g} Mbps x {derate:.0%} = {total:g} Mbps"))
    else:
        out.append((f"rate >= {req['rate_Mbps']:g} Mbps", rate >= req["rate_Mbps"], f"offers up to {rate:g} Mbps"))
    if "max_tx_dBm" in req:
        out.append((f"tx <= {req['max_tx_dBm']:g} dBm", row.tx_dBm <= req["max_tx_dBm"], f"{row.tx_dBm:+g} dBm"))
    return out


SCENARIOS = {
    "S1 Smart home": dict(
        text=["R1 coverage: 2-storey house, longest link ~25 m through walls", "R2 devices: 40 battery sensors/switches",
              "R3 traffic: 2 kbps per sensor -> 0.08 Mbps aggregate", "R4 power: coin-cell, so transmit power <= +10 dBm",
              "R5 battery life: >= 2 years (730 days)"],
        req=dict(link_m=25, rate_Mbps=0.08, max_tx_dBm=10, range_basis="indoor")),
    "S2 Campus network": dict(
        text=["R1 coverage: 0.49 km2 campus (~700 m x 700 m), indoors and outdoors", "R2 users: 5,000 registered, 20% active at once",
              "R3 demand: 1,000 active users x 2 Mbps = 2,000 Mbps aggregate", "R4 power: mains-powered infrastructure, no battery limit",
              "R5 per-link reach: >= 40 m indoors"],
        req=dict(link_m=40, aggregate_Mbps=2000, area_m2=0.49e6, range_basis="indoor")),
    "S3 Metro fixed link": dict(
        text=["R1 distance: single fixed link of 6 km across a river, no fibre", "R2 throughput: >= 20 Mbps sustained",
              "R3 power: mains/solar, no battery constraint", "R4 mobility: none (both ends fixed)"],
        req=dict(link_m=6000, rate_Mbps=20, range_basis="fixed")),
}
print("\n=== Scenario requirement checks ===")
for name, sc in SCENARIOS.items():
    print(f"\n{name}")
    for t in sc["text"]:
        print("   ", t)
    for _, r in df.iterrows():
        res = check(r, sc["req"])
        verdict = "MEETS ALL" if all(ok for _, ok, _ in res) else "fails: " + "; ".join(l for l, ok, _ in res if not ok)
        print(f"  - {r.tech:<24} {verdict}")
        for l, ok, d in res:
            print(f"        [{'ok' if ok else 'NO'}] {l:<26} {d}")
