"""Maps and charts for a region: python src/figures.py <region> (run analyze.py first)."""
import sys
from pathlib import Path
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
R = sys.argv[1]; OUT = Path(__file__).resolve().parents[1] / "outputs" / R
L = np.load(OUT / "layers_coarse.npz"); CID = L["cid"]; FOOT = L["foot"]
C = [("Outside region / nodata", "#ffffff"), ("Open water", "#476ba1"), ("Developed", "#d1271e"), ("Barren", "#b3ac9f"),
     ("Forest + woody wetland", "#1c5f2c"), ("Shrub/grass", "#ccb879"), ("Agriculture (crops + pasture)", "#dcd939"),
     ("Herbaceous wetland", "#6c9fb8")]
cm = ListedColormap([c[1] for c in C]); title = {"houston": "Harris County (Houston), TX", "gulf": "U.S. Gulf Coast counties (TX-FL)"}[R]
rows = np.where((CID > 0).any(axis=1))[0]; cols = np.where((CID > 0).any(axis=0))[0]
sl = (slice(rows[0], rows[-1] + 1), slice(cols[0], cols[-1] + 1))
note = "1938-1992 = MODELED (USGS backcast, 250 m); 1985-2024 = satellite-derived Annual NLCD (30 m, resampled to 250 m). Unmatched classes = nodata."
PAN = [("Backcast_1938", "1938 (modeled)"), ("Backcast_1970", "1970 (modeled)"), ("Backcast_1992", "1992 (modeled)"), ("NLCD_2024", "2024 (Annual NLCD)")]
if R == "gulf":
    fig, axs = plt.subplots(2, 2, figsize=(16, 7.6)); axs = axs.ravel(); tight = (0, 0.12, 1, 0.95)
else:
    fig, axs = plt.subplots(1, 4, figsize=(20, 5.6)); tight = (0, 0.14, 1, 0.95)
for a, (k, t) in zip(axs, PAN):
    a.imshow(L[k][sl], cmap=cm, vmin=0, vmax=7, interpolation="nearest"); a.set_title(t, fontsize=11); a.axis("off")
fig.legend(handles=[Patch(color=c[1], label=c[0]) for c in C[1:]], loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 0.03))
fig.suptitle(f"{title} - harmonized land cover, 1938-2024", fontsize=14); fig.text(0.5, 0.008, note, ha="center", fontsize=8)
fig.tight_layout(rect=tight); fig.savefig(OUT / "fig1_landcover_maps.png", dpi=130); plt.close(fig)

a, b = L["Backcast_1938"], L["NLCD_2024"]; ch = np.zeros(a.shape, np.uint8)
ch[FOOT & (a == b)] = 1; ch[FOOT & (a != b) & (b == 2)] = 2; ch[FOOT & (a != b) & (a == 2)] = 3
ch[FOOT & (a != b) & (b != 2) & (a != 2)] = 4
fig, ax = plt.subplots(figsize=(14, 6.5) if R == "gulf" else (8, 7))
ax.imshow(ch[sl], cmap=ListedColormap(["#ffffff", "#e6e6e6", "#d1271e", "#7b3294", "#f4a582"]), vmin=0, vmax=4, interpolation="nearest"); ax.axis("off")
ax.set_title(f"{title}: class change 1938 -> 2024 (modeled 1938 vs observed 2024)")
ax.legend(handles=[Patch(color=c, label=l) for c, l in zip(["#e6e6e6", "#d1271e", "#7b3294", "#f4a582"],
          ["No class change", "Converted to Developed", "Developed -> other", "Other class change"])], loc="lower left", fontsize=8)
fig.tight_layout(); fig.savefig(OUT / "fig2_change_map.png", dpi=130); plt.close(fig)

t = pd.read_csv(OUT / "class_area_coarse_km2.csv").set_index("class")
yrs = [int(c.split()[-1]) for c in t.columns]
fig, ax = plt.subplots(figsize=(11, 5.5))
for (n, r), col in zip(t.iterrows(), [c[1] for c in C[1:]]): ax.plot(yrs, r.values, "o-", color=col, label=n, lw=2)
ax.axvspan(1938, 1992, color="#999", alpha=0.12); ax.axvline(1998.5, color="#555", ls="--", lw=1)
ax.text(1940, ax.get_ylim()[1] * 0.97, "MODELED backcast", fontsize=9, va="top"); ax.text(2000, ax.get_ylim()[1] * 0.97, "observed (Annual NLCD)", fontsize=9, va="top")
ax.set_ylabel("Area (km2)"); ax.set_xlabel("Year"); ax.set_title(f"{title}: class area through time\n(dashed line = product switch; level shifts across it are partly method differences)", fontsize=11)
ax.legend(fontsize=8, frameon=False, loc="center left", bbox_to_anchor=(1.01, 0.5)); fig.tight_layout(); fig.savefig(OUT / "fig3_class_area_trajectory.png", dpi=130); plt.close(fig)

if R == "gulf":
    bc = pd.read_csv(OUT / "by_county_coarse_km2.csv")
    bc["lab"] = bc.county + ", " + bc.state; bc = bc.sort_values("developed_NLCD_2024_km2", ascending=False).head(15).iloc[::-1]
    fig, ax = plt.subplots(figsize=(9, 6)); y = np.arange(len(bc))
    ax.barh(y + 0.2, bc.developed_NLCD_2024_km2, 0.4, color="#d1271e", label="2024 (observed)")
    ax.barh(y - 0.2, bc.developed_Backcast_1938_km2, 0.4, color="#e8a39d", label="1938 (modeled)")
    ax.set_yticks(y); ax.set_yticklabels(bc.lab); ax.set_xlabel("Developed area (km2)"); ax.legend(frameon=False)
    ax.set_title("Largest developed footprints, Gulf Coast counties/parishes"); fig.tight_layout(); fig.savefig(OUT / "fig4_county_developed.png", dpi=130); plt.close(fig)
