"""Land use / land cover change, 1938 -> 2024, for a county-defined region.

usage: python src/analyze.py <region> <raw_dir>      (region in regions.REGIONS)

Inputs: USGS modeled backcast 1938-1992 (250 m) + Annual NLCD C1.1 (30 m).
All layers go on the 250 m backcast grid; only classes present in BOTH products
are kept (everything else = nodata 0). Tables use the footprint of cells that are
valid in every timeline layer so classes are compared like-for-like.
"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd, rasterio
from rasterio import features
from rasterio.warp import reproject, Resampling
from rasterio.windows import Window, from_bounds
from pyproj import Transformer
from shapely.geometry import shape
from shapely.ops import transform as shp_transform
sys.path.insert(0, str(Path(__file__).parent))
from regions import REGIONS

REGION, RAW = sys.argv[1], Path(sys.argv[2])
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs" / REGION; OUT.mkdir(parents=True, exist_ok=True)
CELL = 0.25 * 0.25                       # km2 per 250 m cell (EPSG:5070 is equal-area)
BACK_YEARS = (1938, 1950, 1970, 1985, 1992)
NLCD_YEARS = (1985, 2005, 2015, 2024)

DETAIL = {1: "Open water", 2: "Developed", 3: "Barren", 4: "Forest", 5: "Shrubland", 6: "Grassland",
          7: "Hay/Pasture", 8: "Cultivated crops", 9: "Woody wetland", 10: "Herbaceous wetland"}
# backcast legend: 6 Mining / 17 Ice have no NLCD counterpart -> nodata
BACK = {1: 1, 2: 2, 7: 3, 8: 4, 9: 4, 10: 4, 12: 5, 11: 6, 14: 7, 13: 8, 16: 9, 15: 10}
NLCD = {11: 1, 21: 2, 22: 2, 23: 2, 24: 2, 31: 3, 41: 4, 42: 4, 43: 4, 52: 5, 71: 6, 81: 7, 82: 8, 90: 9, 95: 10}
COARSE = {1: "Open water", 2: "Developed", 3: "Barren", 4: "Forest + woody wetland", 5: "Shrub/grass",
          6: "Agriculture (crops + pasture)", 7: "Herbaceous wetland"}
D2C = {1: 1, 2: 2, 3: 3, 4: 4, 9: 4, 5: 5, 6: 5, 7: 6, 8: 6, 10: 7}
def lut(m):
    t = np.zeros(256, np.uint8)
    for k, v in m.items(): t[k] = v
    return t
BL, NL, CL = lut(BACK), lut(NLCD), lut(D2C)

# ---- region polygons -> county id raster on the backcast grid -------------
gj = {f["id"]: f for f in json.load(open(RAW / "counties.json"))["features"]}
fips = [(st, f) for st, fl in REGIONS[REGION].items() for f in fl]
src = rasterio.open(RAW / "CONUS_Backcasting_y1992.tif")
tf = Transformer.from_crs("EPSG:4326", src.crs, always_xy=True).transform
geoms = [shp_transform(tf, shape(gj[f]["geometry"])) for _, f in fips]
names = [gj[f]["properties"]["NAME"] for _, f in fips]
minx, miny, maxx, maxy = (min(g.bounds[0] for g in geoms), min(g.bounds[1] for g in geoms),
                          max(g.bounds[2] for g in geoms), max(g.bounds[3] for g in geoms))
win = from_bounds(minx, miny, maxx, maxy, src.transform).round_offsets().round_lengths()
TR, SHAPE = src.window_transform(win), (win.height, win.width)
CID = features.rasterize([(g, i + 1) for i, g in enumerate(geoms)], SHAPE, transform=TR, dtype="uint8")
MASK = CID > 0
print(REGION, "grid", SHAPE, "counties", len(fips), "cells", int(MASK.sum()), flush=True)

def load_back(y):
    with rasterio.open(RAW / f"CONUS_Backcasting_y{y}.tif") as s:
        w = Window(round((TR.c - s.transform.c) / 250), round((s.transform.f - TR.f) / 250), SHAPE[1], SHAPE[0])
        a = s.read(1, window=w)
    return np.where(MASK, BL[a], 0).astype(np.uint8)

def load_nlcd(y, strip=800):
    out = np.zeros(SHAPE, np.uint8)
    with rasterio.open(next((RAW / "nlcd").glob(f"*LndCov_{y}_CU_C1V1.tif"))) as s:
        for r0 in range(0, SHAPE[0], strip):
            r1 = min(SHAPE[0], r0 + strip)
            if not MASK[r0:r1].any(): continue
            tr = TR * rasterio.Affine.translation(0, r0)
            b = rasterio.transform.array_bounds(r1 - r0, SHAPE[1], tr)
            w = from_bounds(*b, s.transform).round_offsets().round_lengths()
            a = NL[s.read(1, window=w)]
            reproject(a, out[r0:r1], src_transform=s.window_transform(w), src_crs=s.crs,
                      dst_transform=tr, dst_crs=src.crs, resampling=Resampling.mode)
    return np.where(MASK, out, 0).astype(np.uint8)

L = {f"Backcast {y}": load_back(y) for y in BACK_YEARS}
for y in NLCD_YEARS:
    L[f"NLCD {y}"] = load_nlcd(y); print("loaded NLCD", y, flush=True)
LC = {k: CL[v] for k, v in L.items()}                       # coarse versions
TL = ["Backcast 1938", "Backcast 1950", "Backcast 1970", "Backcast 1992", "NLCD 2005", "NLCD 2015", "NLCD 2024"]
FOOT = np.all([LC[k] > 0 for k in TL], axis=0)               # common footprint

def areas(layers, scheme, names_, sel, keys):
    return pd.DataFrame([{"class": n, **{k: round(float(((layers[k] == i) & sel).sum() * CELL), 1) for k in keys}}
                         for i, n in names_.items()])
areas(LC, "c", COARSE, FOOT, TL).to_csv(OUT / "class_area_coarse_km2.csv", index=False)
areas(L, "d", DETAIL, np.all([L[k] > 0 for k in TL], axis=0), TL).to_csv(OUT / "class_area_detailed_km2.csv", index=False)

def trans(a, b, sel, n, nm):
    m = np.zeros((n, n)); np.add.at(m, (a[sel] - 1, b[sel] - 1), CELL)
    return pd.DataFrame(m.round(1), index=[nm[i] for i in range(1, n + 1)], columns=[nm[i] for i in range(1, n + 1)])
trans(LC["Backcast 1938"], LC["NLCD 2024"], FOOT, 7, COARSE).to_csv(OUT / "transition_coarse_1938_to_2024_km2.csv")
trans(LC["Backcast 1992"], LC["NLCD 2024"], FOOT, 7, COARSE).to_csv(OUT / "transition_coarse_1992_to_2024_km2.csv")
trans(L["Backcast 1938"], L["NLCD 2024"], np.all([L[k] > 0 for k in TL], axis=0), 10, DETAIL).to_csv(OUT / "transition_detailed_1938_to_2024_km2.csv")

# ---- 1985 overlap check (both products exist) -----------------------------
ok = (LC["Backcast 1985"] > 0) & (LC["NLCD 1985"] > 0)
ct = pd.crosstab(pd.Series(LC["Backcast 1985"][ok], name="backcast_1985"), pd.Series(LC["NLCD 1985"][ok], name="nlcd_1985")) * CELL
ct.index = [COARSE[i] for i in ct.index]; ct.columns = [COARSE[i] for i in ct.columns]
ct.round(1).to_csv(OUT / "overlap_1985_confusion_coarse_km2.csv")
okd = (L["Backcast 1985"] > 0) & (L["NLCD 1985"] > 0)
rows = []
for i, (st, f) in enumerate(fips):
    m = ok & (CID == i + 1); md = okd & (CID == i + 1); n = int(m.sum())
    if n == 0: continue
    dv = lambda k: round(float(((LC[k] == 2) & m).sum() * CELL), 1)
    rows.append({"state": st, "county": names[i], "fips": f, "land_water_km2": round(n * CELL, 1),
                 "agree_coarse_pct": round(100 * float((LC["Backcast 1985"][m] == LC["NLCD 1985"][m]).mean()), 1),
                 "agree_detailed_pct": round(100 * float((L["Backcast 1985"][md] == L["NLCD 1985"][md]).mean()), 1),
                 "developed_backcast85_km2": dv("Backcast 1985"), "developed_nlcd85_km2": dv("NLCD 1985")})
agr = pd.DataFrame(rows); agr.to_csv(OUT / "overlap_1985_by_county.csv", index=False)

# ---- per-county / per-state developed area and class shares ---------------
rows = []
for i, (st, f) in enumerate(fips):
    m = FOOT & (CID == i + 1)
    if not m.any(): continue
    r = {"state": st, "county": names[i], "fips": f, "footprint_km2": round(float(m.sum() * CELL), 1)}
    for k in TL: r[f"developed_{k.replace(' ', '_')}_km2"] = round(float(((LC[k] == 2) & m).sum() * CELL), 1)
    r["agri_1938_km2"] = round(float(((LC["Backcast 1938"] == 6) & m).sum() * CELL), 1)
    r["agri_2024_km2"] = round(float(((LC["NLCD 2024"] == 6) & m).sum() * CELL), 1)
    rows.append(r)
bycounty = pd.DataFrame(rows); bycounty.to_csv(OUT / "by_county_coarse_km2.csv", index=False)
st_of = np.array([s for s, _ in fips])
bystate = []
for st in REGIONS[REGION]:
    m = FOOT & np.isin(CID, np.where(st_of == st)[0] + 1)
    bystate.append({"state": st, **{k: round(float(((LC[k] == 2) & m).sum() * CELL), 1) for k in TL}})
pd.DataFrame(bystate).to_csv(OUT / "developed_by_state_km2.csv", index=False)

stats = {"region": REGION, "counties": len(fips), "region_cells": int(MASK.sum()),
         "region_km2": round(float(MASK.sum() * CELL), 1), "common_footprint_km2": round(float(FOOT.sum() * CELL), 1),
         "pct_footprint": round(100 * float(FOOT.sum() / MASK.sum()), 1),
         "agree_coarse_1985_pct": round(100 * float((LC["Backcast 1985"][ok] == LC["NLCD 1985"][ok]).mean()), 1),
         "agree_detailed_1985_pct": round(100 * float((L["Backcast 1985"][okd] == L["NLCD 1985"][okd]).mean()), 1),
         "pct_cells_class_changed_1938_2024": round(100 * float((LC["Backcast 1938"][FOOT] != LC["NLCD 2024"][FOOT]).mean()), 1)}
json.dump(stats, open(OUT / "summary.json", "w"), indent=2)
np.savez_compressed(OUT / "layers_coarse.npz", cid=CID, foot=FOOT, **{k.replace(" ", "_"): v for k, v in LC.items()})
json.dump({"transform": list(TR)[:6], "crs": src.crs.to_string(), "counties": [[s, n, f] for (s, f), n in zip(fips, names)]}, open(OUT / "grid.json", "w"))
print(pd.read_csv(OUT / "class_area_coarse_km2.csv").to_string(index=False)); print(stats)
