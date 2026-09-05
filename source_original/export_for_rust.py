#!/usr/bin/env python
"""
Export LEMON EC+EO resting EEG to flat binaries for the Rust dynamic-coherence kernel.

For every subject with both conditions + WST + LPS + age:
  - load EC and EO, set montage, pick the channel set COMMON TO ALL SUBJECTS
    (fixed order → constant n_ch, so dimensionality features aren't confounded by channel count)
  - write raw float32 (row-major n_ch × n_samples) for each condition
  - write one global manifest.csv (sid, wst, lps, age, sfreq, n_ch, n_ec, n_eo)
  - write channels.json (names + 3D positions + pairwise distances)

Output dir: experiments/eeg_coherence/rust_export/
"""
import os, csv, json, warnings
import numpy as np
import mne
from multiprocessing import Pool, cpu_count

warnings.filterwarnings("ignore"); mne.set_log_level("ERROR")

DATA = "/Users/ximon/Documents/interp/data/lemon"
EEG  = os.path.join(DATA, "eeg_preprocessed")
BEH  = os.path.join(DATA, "Behavioural_Data_MPILMBB_LEMON")
OUT  = "/Users/ximon/Documents/interp/experiments/eeg_coherence/rust_export"
os.makedirs(OUT, exist_ok=True)


def load_cog():
    out = {}
    with open(os.path.join(BEH, "Cognitive_Test_Battery_LEMON/WST/WST.csv")) as f:
        for r in csv.reader(f):
            if r[0].startswith("sub-") and len(r) > 3 and r[3].strip():
                out.setdefault(r[0], {})["WST"] = float(r[3])
    with open(os.path.join(BEH, "Cognitive_Test_Battery_LEMON/LPS/LPS.csv")) as f:
        for r in csv.reader(f):
            if r[0].startswith("sub-") and len(r) > 1 and r[1].strip():
                out.setdefault(r[0], {})["LPS"] = float(r[1])
    with open(os.path.join(BEH, "META_File_IDs_Age_Gender_Education_Drug_Smoke_SKID_LEMON.csv")) as f:
        rd = csv.reader(f); next(rd)
        for r in rd:
            if r[0] in out and len(r) > 2 and "-" in str(r[2]):
                a, b = str(r[2]).split("-"); out[r[0]]["age"] = (float(a) + float(b)) / 2
    return {k: v for k, v in out.items() if all(x in v for x in ("WST", "LPS", "age"))}


TARGET_SFREQ = 125.0


def ch_names_of(sid):
    """Return (EC_channels, EO_channels) sets for a subject."""
    out = []
    for cond in ("EC", "EO"):
        try:
            raw = mne.io.read_raw_eeglab(os.path.join(EEG, sid, f"{sid}_{cond}.set"),
                                         preload=False, verbose=False)
            out.append(set(raw.ch_names))
        except Exception:
            out.append(set())
    return tuple(out)


_COMMON = None   # set in each worker via initializer
_COG = None


def _init_worker(common, cog):
    global _COMMON, _COG
    _COMMON = common
    _COG = cog


def export_one(sid):
    """EC is required (primary, within-rest dynamics). EO exported if available (secondary)."""
    try:
        rows = {"EC": None, "EO": None}
        for cond in ("EC", "EO"):
            setf = os.path.join(EEG, sid, f"{sid}_{cond}.set")
            if not os.path.exists(setf):
                continue
            raw = mne.io.read_raw_eeglab(setf, preload=True, verbose=False)
            if not set(_COMMON) <= set(raw.ch_names):
                continue  # this condition lacks a canonical channel; skip it
            raw.pick(_COMMON)
            raw.reorder_channels(_COMMON)
            if raw.info["sfreq"] > TARGET_SFREQ + 1:
                raw.resample(TARGET_SFREQ, verbose=False)
            data = raw.get_data().astype(np.float32)   # (n_ch, n_samples)
            data.tofile(os.path.join(OUT, f"{sid}_{cond}.f32"))
            rows[cond] = (data.shape[1], raw.info["sfreq"])
        if rows["EC"] is None:
            return None  # EC is mandatory
        sf = rows["EC"][1]
        n_ec = rows["EC"][0]
        n_eo = rows["EO"][0] if rows["EO"] else 0
        return (sid, _COG[sid]["WST"], _COG[sid]["LPS"], _COG[sid]["age"],
                sf, len(_COMMON), n_ec, n_eo)
    except Exception as e:
        print(f"  FAIL {sid}: {e}")
        return None


def main():
    cog = load_cog()
    subs = [s for s in sorted(cog)
            if os.path.exists(os.path.join(EEG, s, f"{s}_EC.set"))
            and os.path.exists(os.path.join(EEG, s, f"{s}_EO.set"))]
    print(f"Subjects with EC+EO+cognitive: {len(subs)}")

    # ── Pass 1: canonical channel set = channels present in ≥95% of recordings ──
    print("Finding canonical channel set...")
    with Pool(cpu_count() - 2) as p:
        chsets = p.map(ch_names_of, subs)  # list of (EC_set, EO_set)

    from collections import Counter
    cnt = Counter()
    n_rec = 0
    for ec, eo in chsets:
        if ec: cnt.update(ec); n_rec += 1
        if eo: cnt.update(eo); n_rec += 1
    thresh = 0.95 * n_rec
    canonical = sorted([ch for ch, c in cnt.items() if c >= thresh])
    print(f"Canonical channels (≥95% of {n_rec} recordings): {len(canonical)}")

    # Keep subjects whose EC contains all canonical channels (EC is primary).
    # EO is exported when it also has them (handled per-subject in export_one).
    canon_set = set(canonical)
    keep = [s for s, (ec, eo) in zip(subs, chsets) if canon_set <= ec]
    n_both = sum(1 for s, (ec, eo) in zip(subs, chsets)
                 if canon_set <= ec and canon_set <= eo)
    print(f"Subjects with canonical channels in EC (primary): {len(keep)} of {len(subs)}")
    print(f"  of which also have EO (secondary contrast): {n_both}")
    subs = keep
    common = canonical

    # reference montage positions
    raw0 = mne.io.read_raw_eeglab(os.path.join(EEG, subs[0], f"{subs[0]}_EC.set"),
                                  preload=False, verbose=False)
    raw0.set_montage("standard_1005", match_case=False, on_missing="ignore", verbose=False)
    pos = raw0.get_montage().get_positions()["ch_pos"]
    positions = np.array([pos[c] for c in common])  # (n_ch, 3)
    D = np.linalg.norm(positions[:, None, :] - positions[None, :, :], axis=-1)

    with open(os.path.join(OUT, "channels.json"), "w") as f:
        json.dump({"channels": common,
                   "positions": positions.tolist(),
                   "distances": D.tolist()}, f)

    # ── Pass 2: export each subject ──
    print(f"Exporting {len(subs)} subjects...")
    with Pool(cpu_count() - 2, initializer=_init_worker, initargs=(common, cog)) as p:
        results = [r for r in p.map(export_one, subs) if r is not None]

    with open(os.path.join(OUT, "manifest.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sid", "WST", "LPS", "age", "sfreq", "n_ch", "n_ec", "n_eo"])
        w.writerows(results)

    print(f"Exported {len(results)} subjects to {OUT}")
    print(f"  manifest.csv, channels.json, {len(results)*2} .f32 files")


if __name__ == "__main__":
    main()
