#!/usr/bin/env python3
"""
01_etl.py — Cohort construction for Mg → AKI study (v5: risk-set design)

Reads raw PhysioNet CSVs → outputs analysis-ready files.
No lab summary columns; R computes covariates at match time.

Outputs per database:
  did_all_{db}.csv       — All patients, time-invariant covariates
  did_labs_all_{db}.csv  — All electrolyte/vital measurements with timestamps
  did_cr_all_{db}.csv    — All creatinine measurements with timestamps
  did_consort.csv        — CONSORT numbers (appended)

Usage:
  python 01_etl.py              # both databases
  python 01_etl.py eicu         # eICU only
  python 01_etl.py mimic        # MIMIC only
"""

import os
import sys
from importlib.util import module_from_spec, spec_from_file_location

import duckdb
import numpy as np
import pandas as pd

_spec = spec_from_file_location(
    "config", os.path.join(os.path.dirname(__file__), "00_config.py")
)
cfg = module_from_spec(_spec)
_spec.loader.exec_module(cfg)

# Pull everything from config into local namespace
for _k in dir(cfg):
    if not _k.startswith("_"):
        globals()[_k] = getattr(cfg, _k)


# ═══════════════════════════════════════════════════════════════════
# HELPERS
# ═══════════════════════════════════════════════════════════════════
def _resolve(root, name):
    for ext in [".csv.gz", ".csv"]:
        p = os.path.join(root, name + ext)
        if os.path.exists(p):
            return p
    return None


def load_filtered(name, root, pids=None, pid_col="patientunitstayid", extra_where=None):
    """Load CSV via DuckDB with predicate pushdown."""
    p = _resolve(root, name)
    if p is None:
        print(f"  WARN: {name} not found in {root}")
        return pd.DataFrame()
    con = duckdb.connect()
    try:
        clauses = []
        if pids is not None:
            pid_df = pd.DataFrame({"_pid": pd.array(list(pids), dtype="int64")})
            con.register("_pf", pid_df)
            clauses.append(f'"{pid_col}" IN (SELECT _pid FROM _pf)')
        if extra_where:
            clauses.append(f"({extra_where})")
        where = " WHERE " + " AND ".join(clauses) if clauses else ""
        df = con.execute(
            f"SELECT * FROM read_csv_auto('{p}', header=true, ignore_errors=true){where}"
        ).df()
        df.columns = df.columns.str.lower()
        print(f"  {name}: {len(df):,} rows")
        return df
    finally:
        con.close()


def load_pd(path, **kw):
    """Load CSV with pandas (for small tables)."""
    for ext in ["", ".gz"]:
        fp = path + ext if not path.endswith((".csv", ".csv.gz")) else path
        if os.path.exists(fp):
            df = pd.read_csv(fp, low_memory=False, **kw)
            print(f"  {os.path.basename(fp)}: {len(df):,} rows")
            return df
    return pd.DataFrame()


def gz(p):
    return p if os.path.exists(p) else p.replace(".csv.gz", ".csv")


def matches_any(series, patterns):
    s = series.astype(str).str.lower()
    mask = pd.Series(False, index=series.index)
    for p in patterns:
        mask |= s.str.contains(p.lower(), na=False)
    return mask


def compute_egfr(cr, age, is_female):
    cr = np.asarray(cr, dtype=np.float64)
    age = np.asarray(age, dtype=np.float64)
    fem = np.asarray(is_female, dtype=bool)
    kappa = np.where(fem, 0.7, 0.9)
    alpha = np.where(fem, -0.241, -0.302)
    r = cr / kappa
    return (
        142
        * np.power(np.minimum(r, 1.0), alpha)
        * np.power(np.maximum(r, 1.0), -1.200)
        * np.power(0.9938, age)
        * np.where(fem, 1.012, 1.0)
    )


def matches_icd(dx_df, hadm_ids, code_map):
    sub = dx_df[dx_df.hadm_id.isin(hadm_ids)]
    hits = set()
    for ver, prefixes in code_map.items():
        v = sub[sub.icd_version == ver]
        for p in prefixes:
            hits |= set(v[v.icd_code.str.startswith(p)].hadm_id)
    return hits


def pct(n, total):
    return f"{n:,} ({100*n/max(total,1):.1f}%)"


def desc(s, name):
    v = s.dropna()
    if len(v) == 0:
        print(f"    {name}: all missing")
        return
    print(
        f"    {name}: n={len(v)}, median={v.median():.2f}, "
        f"IQR=[{v.quantile(.25):.2f}–{v.quantile(.75):.2f}]"
    )


def save_all_patients(treated, control, db_tag):
    """Merge treated + control into did_all_{db}.csv."""
    pid_col = (
        "patientunitstayid" if "patientunitstayid" in treated.columns else "stay_id"
    )

    for df in [treated, control]:
        df["surg_cabg"] = (df.surgery_type == "cabg").astype(int)
        df["surg_valve"] = (df.surgery_type == "valve").astype(int)
        df["surg_combined"] = (df.surgery_type == "combined").astype(int)
    treated["treated"] = 1
    control["treated"] = 0

    for c in ["mg_offset_h", "mg_offset_min"]:
        if c not in control.columns:
            control[c] = np.nan

    # Rename pid BEFORE selecting columns
    treated = treated.rename(columns={pid_col: "pid"})
    control = control.rename(columns={pid_col: "pid"})
    keep = [
        c for c in ALL_PATIENTS_COLS if c in treated.columns and c in control.columns
    ]
    all_pts = pd.concat([treated[keep], control[keep]], ignore_index=True)

    tag = db_tag.lower()
    path = os.path.join(RESULTS, f"did_all_{tag}.csv")
    all_pts.to_csv(path, index=False)
    n_trt = int(all_pts.treated.sum())
    print(
        f"\n  ✓ did_all_{tag}.csv: {len(all_pts):,} pts ({n_trt} treated + {len(all_pts)-n_trt} control)"
    )
    return all_pts


# ═══════════════════════════════════════════════════════════════════
#  eICU-CRD
# ═══════════════════════════════════════════════════════════════════
def run_eicu():
    SEP = "=" * 70
    print(f"\n{SEP}\n01_etl.py — eICU-CRD (v5: risk-set design)\n{SEP}")
    consort = {}

    # ── 1. Patient table → cardiac pids ──────────────────────────
    patient = load_pd(os.path.join(EICU_ROOT, "patient.csv.gz"))
    consort["total_icu"] = len(patient)

    mask = matches_any(
        patient.apacheadmissiondx, CARDIAC_DX_PATTERNS
    ) | patient.unittype.isin(CARDIAC_UNIT_TYPES)
    cardiac = patient[mask].copy()
    cardiac["age"] = pd.to_numeric(
        cardiac["age"].astype(str).str.replace(">", ""), errors="coerce"
    ).fillna(90)
    cardiac = cardiac[cardiac.age >= MIN_AGE]
    cardiac = (
        cardiac.sort_values("hospitaladmitoffset")
        .groupby("uniquepid")
        .first()
        .reset_index()
    )
    cardiac["is_female"] = (cardiac.gender.str.lower() == "female").astype(int)
    cardiac["surgery_type"] = cardiac.apacheadmissiondx.apply(
        lambda s: (
            "combined"
            if any(p in str(s).lower() for p in ["cabg", "coronary artery bypass"])
            and any(p in str(s).lower() for p in ["valve"])
            else (
                "cabg"
                if any(p in str(s).lower() for p in ["cabg", "coronary artery bypass"])
                else "valve" if "valve" in str(s).lower() else "other_cardiac"
            )
        )
    )
    pids = set(cardiac.patientunitstayid)
    consort["cardiac_adult_first"] = len(cardiac)
    print(f"  Cardiac surgery, adult, 1st stay: {len(cardiac):,}")

    # ── 2. DuckDB load big tables ────────────────────────────────
    print(f"\n  Loading tables via DuckDB ({len(pids):,} pids)...")
    lab = load_filtered("lab", EICU_ROOT, pids)
    med = load_filtered("medication", EICU_ROOT, pids)
    inf = load_filtered("infusionDrug", EICU_ROOT, pids)
    diag = load_filtered("diagnosis", EICU_ROOT, pids)
    pasthx = load_filtered("pastHistory", EICU_ROOT, pids)
    vital = load_filtered("vitalPeriodic", EICU_ROOT, pids)
    admDrug = load_filtered("admissionDrug", EICU_ROOT, pids)

    # ── 3. ESKD exclusion ────────────────────────────────────────
    eskd = set()
    if len(pasthx) > 0 and "pasthistorypath" in pasthx.columns:
        eskd |= set(
            pasthx[
                pasthx.patientunitstayid.isin(pids)
                & matches_any(pasthx.pasthistorypath, ESKD_PATTERNS)
            ].patientunitstayid
        )
    if len(diag) > 0:
        eskd |= set(
            diag[
                diag.patientunitstayid.isin(pids)
                & matches_any(diag.diagnosisstring, ESKD_PATTERNS)
            ].patientunitstayid
        )
    cardiac = cardiac[~cardiac.patientunitstayid.isin(eskd)]
    pids = set(cardiac.patientunitstayid)
    consort["post_eskd"] = len(cardiac)
    print(f"  ESKD excluded: {len(eskd):,} → {len(cardiac):,}")

    # ── 4. IV Mg identification ──────────────────────────────────
    print("\n  IV Magnesium identification...")
    med_ok = (
        med[
            ~med.get("drugordercancelled", "").astype(str).str.contains("Yes", na=False)
        ]
        if len(med) > 0
        else pd.DataFrame()
    )
    frames = []
    if len(med_ok) > 0:
        mg_m = med_ok[
            matches_any(med_ok.drugname, MG_SUPP_PATTERNS)
            & (med_ok.drugstartoffset >= 0)
        ]
        if "routeadmin" in mg_m.columns:
            mg_m = mg_m[
                mg_m.routeadmin.str.lower().str.contains(
                    "iv|intravenous|inject", na=False
                )
                | mg_m.routeadmin.isna()
            ]
        if len(mg_m) > 0:
            frames.append(
                mg_m[["patientunitstayid", "drugstartoffset"]].rename(
                    columns={"drugstartoffset": "mg_offset_min"}
                )
            )
    if len(inf) > 0 and "drugname" in inf.columns:
        mg_i = inf[
            inf.patientunitstayid.isin(pids)
            & matches_any(inf.drugname, MG_SUPP_PATTERNS)
            & (inf.infusionoffset >= 0)
        ]
        if len(mg_i) > 0:
            frames.append(
                mg_i[["patientunitstayid", "infusionoffset"]].rename(
                    columns={"infusionoffset": "mg_offset_min"}
                )
            )

    if not frames:
        print("  ERROR: no IV Mg found")
        return None

    mg_all = pd.concat(frames, ignore_index=True)
    first_mg = (
        mg_all.sort_values("mg_offset_min")
        .groupby("patientunitstayid")
        .first()
        .reset_index()
    )
    first_mg["mg_offset_h"] = first_mg.mg_offset_min / 60.0

    treated_pids = set(first_mg.patientunitstayid)
    control_pids = pids - treated_pids
    consort["treated_any_ivmg"] = len(treated_pids)
    consort["control_no_ivmg"] = len(control_pids)
    print(f"  Treated (any IV Mg): {pct(len(treated_pids), len(pids))}")
    desc(first_mg.mg_offset_h, "Onset (h)")
    for c in [6, 12, 24, 48]:
        print(f"    ≤{c}h: {pct((first_mg.mg_offset_h<=c).sum(), len(first_mg))}")

    # ── 5. Creatinine (for eGFR + exclusions) ────────────────────
    cr = lab[
        lab.patientunitstayid.isin(pids)
        & lab.labname.str.lower().str.contains("creatinine", na=False)
        & lab.labresult.between(CR_MIN, CR_MAX)
    ].copy()

    # First postop Cr per patient (for eGFR — time-invariant baseline)
    cr_first = (
        cr[cr.labresultoffset >= 0]
        .sort_values("labresultoffset")
        .groupby("patientunitstayid")
        .first()
        .reset_index()
    )

    # Treated: need Cr before Mg for AKI-at-T0 check
    cr_t = cr[cr.patientunitstayid.isin(treated_pids)].merge(
        first_mg[["patientunitstayid", "mg_offset_min"]], on="patientunitstayid"
    )
    cr_pre_t = cr_t[
        (cr_t.labresultoffset >= 0) & (cr_t.labresultoffset < cr_t.mg_offset_min)
    ]
    has_cr_pre = set(cr_pre_t.patientunitstayid.unique())
    consort["treated_has_cr_pre"] = len(has_cr_pre)
    consort["treated_no_cr_pre"] = len(treated_pids) - len(has_cr_pre)

    # Hospital Cr fallback
    hosp_off = cardiac.set_index("patientunitstayid")["hospitaladmitoffset"].to_dict()
    cr_t["hosp_off"] = cr_t.patientunitstayid.map(hosp_off)
    cr_h = cr_t[
        (cr_t.labresultoffset >= cr_t.hosp_off - 360)
        & (cr_t.labresultoffset <= cr_t.hosp_off + 360)
        & (cr_t.labresultoffset < 0)
    ].copy()
    cr_h["_dist"] = (cr_h.labresultoffset - cr_h.hosp_off).abs()
    cr_hosp = (
        cr_h.sort_values(["patientunitstayid", "_dist"])
        .groupby("patientunitstayid")
        .first()
        .reset_index()
    )
    rescue_pids = set(cr_hosp.patientunitstayid) - has_cr_pre
    n_rescue = len(rescue_pids & treated_pids)
    consort["treated_hosp_fallback"] = n_rescue
    has_cr_pre |= rescue_pids
    consort["treated_has_cr_pre_with_fallback"] = len(has_cr_pre & treated_pids)
    print(f"  Cr_pre: ICU={consort['treated_has_cr_pre']}, +fallback={n_rescue}")

    # Exclusions: Cr_pre >= 4.0
    high_cr = set(cr_first[cr_first.labresult >= BASELINE_CR_MAX].patientunitstayid)
    consort["excl_cr_high"] = len(high_cr & has_cr_pre)

    # Prevalent AKI (at time of first Cr — simplified for ETL; R does T0-specific check)
    # ETL excludes obvious prevalent AKI; R does precise T0-based check
    keep_treated = (has_cr_pre & treated_pids) - high_cr
    consort["treated_final"] = len(keep_treated)

    # Controls: need ≥2 Cr
    cr_ctrl = cr[(cr.patientunitstayid.isin(control_pids)) & (cr.labresultoffset >= 0)]
    c2 = set(
        cr_ctrl.groupby("patientunitstayid").size().pipe(lambda s: s[s >= 2]).index
    )
    c_excl = set(
        cr_first[
            cr_first.patientunitstayid.isin(c2)
            & (cr_first.labresult >= BASELINE_CR_MAX)
        ].patientunitstayid
    )
    keep_control = c2 - c_excl
    consort["control_has_2cr"] = len(c2)
    consort["control_final"] = len(keep_control)

    # ── 6. Build treated DataFrame ───────────────────────────────
    print(f"\n  Building treated: {len(keep_treated):,}")
    treated = cardiac[cardiac.patientunitstayid.isin(keep_treated)].copy()
    treated = treated.merge(
        first_mg[["patientunitstayid", "mg_offset_min", "mg_offset_h"]],
        on="patientunitstayid",
    )

    # First Cr → eGFR
    treated = treated.merge(
        cr_first[["patientunitstayid", "labresult"]].rename(
            columns={"labresult": "first_cr"}
        ),
        on="patientunitstayid",
        how="left",
    )
    treated["egfr"] = compute_egfr(treated.first_cr, treated.age, treated.is_female)

    # BMI
    if "admissionweight" in cardiac.columns and "admissionheight" in cardiac.columns:
        wt = cardiac.set_index("patientunitstayid")["admissionweight"].to_dict()
        ht = cardiac.set_index("patientunitstayid")["admissionheight"].to_dict()
        treated["bmi"] = treated.patientunitstayid.apply(
            lambda x: (
                wt.get(x, np.nan) / ((ht.get(x, np.nan) / 100) ** 2)
                if pd.notna(wt.get(x)) and pd.notna(ht.get(x)) and ht.get(x, 0) > 0
                else np.nan
            )
        )
        treated.loc[~treated.bmi.between(10, 80), "bmi"] = np.nan

    # Comorbidities
    for como, pats in EICU_COMORB.items():
        cp = set()
        if len(pasthx) > 0 and "pasthistorypath" in pasthx.columns:
            cp = set(
                pasthx[
                    pasthx.patientunitstayid.isin(keep_treated)
                    & matches_any(pasthx.pasthistorypath, pats)
                ].patientunitstayid
            )
        treated[como] = treated.patientunitstayid.isin(cp).astype(int)

    # Chronic drugs
    if len(admDrug) > 0 and "drugname" in admDrug.columns:
        ad = admDrug[admDrug.patientunitstayid.isin(keep_treated)]
        for col, pats in CHRONIC_DRUG_CLASSES.items():
            hits = (
                set(ad[matches_any(ad.drugname, pats)].patientunitstayid)
                if len(ad) > 0
                else set()
            )
            treated[col] = treated.patientunitstayid.isin(hits).astype(int)
    else:
        for col in CHRONIC_DRUG_CLASSES:
            treated[col] = 0

    # ICU discharge + mortality
    treated["icu_discharge_h"] = treated.patientunitstayid.map(
        (cardiac.set_index("patientunitstayid")["unitdischargeoffset"] / 60).to_dict()
    )
    treated["icu_outcome"] = treated.patientunitstayid.map(
        cardiac.set_index("patientunitstayid")["unitdischargestatus"]
        .str.lower()
        .to_dict()
    )
    treated["hosp_mortality"] = (
        (treated.icu_outcome == "expired").astype(int)
        if "hospitaldischargestatus" not in cardiac.columns
        else treated.patientunitstayid.map(
            cardiac.set_index("patientunitstayid")["hospitaldischargestatus"]
            .str.lower()
            .eq("expired")
            .to_dict()
        )
        .fillna(0)
        .astype(int)
    )

    # ── 7. Build control DataFrame ───────────────────────────────
    print(f"  Building control: {len(keep_control):,}")
    control = cardiac[cardiac.patientunitstayid.isin(keep_control)].copy()

    control = control.merge(
        cr_first[["patientunitstayid", "labresult"]].rename(
            columns={"labresult": "first_cr"}
        ),
        on="patientunitstayid",
        how="left",
    )
    control["egfr"] = compute_egfr(control.first_cr, control.age, control.is_female)

    if "admissionweight" in cardiac.columns:
        control["bmi"] = control.patientunitstayid.apply(
            lambda x: (
                wt.get(x, np.nan) / ((ht.get(x, np.nan) / 100) ** 2)
                if pd.notna(wt.get(x)) and pd.notna(ht.get(x)) and ht.get(x, 0) > 0
                else np.nan
            )
        )
        control.loc[~control.bmi.between(10, 80), "bmi"] = np.nan

    for como, pats in EICU_COMORB.items():
        cp = set()
        if len(pasthx) > 0 and "pasthistorypath" in pasthx.columns:
            cp = set(
                pasthx[
                    pasthx.patientunitstayid.isin(keep_control)
                    & matches_any(pasthx.pasthistorypath, pats)
                ].patientunitstayid
            )
        control[como] = control.patientunitstayid.isin(cp).astype(int)

    if len(admDrug) > 0 and "drugname" in admDrug.columns:
        ad_c = admDrug[admDrug.patientunitstayid.isin(keep_control)]
        for col, pats in CHRONIC_DRUG_CLASSES.items():
            hits = (
                set(ad_c[matches_any(ad_c.drugname, pats)].patientunitstayid)
                if len(ad_c) > 0
                else set()
            )
            control[col] = control.patientunitstayid.isin(hits).astype(int)
    else:
        for col in CHRONIC_DRUG_CLASSES:
            control[col] = 0

    control["icu_discharge_h"] = control.patientunitstayid.map(
        (cardiac.set_index("patientunitstayid")["unitdischargeoffset"] / 60).to_dict()
    )
    control["icu_outcome"] = control.patientunitstayid.map(
        cardiac.set_index("patientunitstayid")["unitdischargestatus"]
        .str.lower()
        .to_dict()
    )
    control["hosp_mortality"] = (
        control.patientunitstayid.map(
            cardiac.set_index("patientunitstayid")
            .get(
                "hospitaldischargestatus",
                cardiac.set_index("patientunitstayid")["unitdischargestatus"],
            )
            .str.lower()
            .eq("expired")
            .to_dict()
        )
        .fillna(0)
        .astype(int)
    )

    # ── 8. Secondary outcomes ────────────────────────────────────
    print("\n  Secondary outcomes...")
    all_p = keep_treated | keep_control

    # POAF
    prior_af = set()
    if len(pasthx) > 0 and "pasthistorypath" in pasthx.columns:
        prior_af = set(
            pasthx[
                pasthx.patientunitstayid.isin(all_p)
                & matches_any(pasthx.pasthistorypath, AF_PRIOR_PATTERNS_EICU)
            ].patientunitstayid
        )
    af_dx = diag[
        diag.patientunitstayid.isin(all_p)
        & matches_any(diag.diagnosisstring, AF_PATTERNS_EICU)
        & (diag.diagnosisoffset > 0)
    ]
    af_pids = set(af_dx.patientunitstayid) - prior_af
    for df in [treated, control]:
        df["poaf"] = df.patientunitstayid.isin(af_pids).astype(int)
    print(f"    POAF: trt={treated.poaf.sum()}, ctl={control.poaf.sum()}")

    # Encephalopathy
    enc = set(
        diag[
            diag.patientunitstayid.isin(all_p)
            & matches_any(diag.diagnosisstring, ENCEPH_PATTERNS_EICU)
            & (diag.diagnosisoffset > 0)
        ].patientunitstayid
    )
    for df in [treated, control]:
        df["encephalopathy"] = df.patientunitstayid.isin(enc).astype(int)

    # Ventricular arrhythmia
    varr = set(
        diag[
            diag.patientunitstayid.isin(all_p)
            & matches_any(diag.diagnosisstring, VARR_PATTERNS_EICU)
            & (diag.diagnosisoffset > 0)
        ].patientunitstayid
    )
    for df in [treated, control]:
        df["vent_arrhythmia"] = df.patientunitstayid.isin(varr).astype(int)

    # ── 9. Export ────────────────────────────────────────────────
    save_all_patients(treated, control, "eicu")

    # did_cr_all: ALL Cr measurements
    all_ids = keep_treated | keep_control
    cr_exp = cr[
        cr.patientunitstayid.isin(all_ids)
        & (cr.labresultoffset >= 0)
        & (cr.labresult <= CR_POST_PLAUSIBLE_MAX)
    ][["patientunitstayid", "labresult", "labresultoffset"]].copy()
    cr_exp["offset_h"] = cr_exp.labresultoffset / 60.0
    cr_exp.to_csv(os.path.join(RESULTS, "did_cr_all_eicu.csv"), index=False)
    print(
        f"  ✓ did_cr_all_eicu.csv: {len(cr_exp):,} Cr for {cr_exp.patientunitstayid.nunique():,} pts"
    )

    # did_labs_all: ALL electrolyte + vital measurements
    lab_rows = []
    for lab_name, patterns in EICU_LAB_PATTERNS.items():
        sub = lab[
            lab.patientunitstayid.isin(all_ids)
            & matches_any(lab.labname, patterns)
            & (lab.labresultoffset >= 0)
        ]
        if len(sub) > 0:
            lab_rows.append(
                pd.DataFrame(
                    {
                        "patientunitstayid": sub.patientunitstayid,
                        "lab_name": lab_name,
                        "value": sub.labresult,
                        "offset_h": sub.labresultoffset / 60.0,
                    }
                )
            )
    # Heart rate from vitalPeriodic
    if len(vital) > 0 and "heartrate" in vital.columns:
        hr = vital[
            vital.patientunitstayid.isin(all_ids)
            & vital.heartrate.notna()
            & vital.heartrate.between(20, 250)
            & (vital.observationoffset >= 0)
        ]
        if len(hr) > 0:
            lab_rows.append(
                pd.DataFrame(
                    {
                        "patientunitstayid": hr.patientunitstayid,
                        "lab_name": "heartrate",
                        "value": hr.heartrate,
                        "offset_h": hr.observationoffset / 60.0,
                    }
                )
            )

    labs_all = pd.concat(lab_rows, ignore_index=True) if lab_rows else pd.DataFrame()
    labs_all.to_csv(os.path.join(RESULTS, "did_labs_all_eicu.csv"), index=False)
    print(f"  ✓ did_labs_all_eicu.csv: {len(labs_all):,} measurements")
    for ln in labs_all.lab_name.unique():
        n = (labs_all.lab_name == ln).sum()
        np_ = labs_all[labs_all.lab_name == ln].patientunitstayid.nunique()
        print(f"    {ln}: {n:,} values across {np_:,} pts")

    # CONSORT
    consort["db"] = "eICU"
    print(f"\n  CONSORT: {consort}")
    return consort


# ═══════════════════════════════════════════════════════════════════
#  MIMIC-IV
# ═══════════════════════════════════════════════════════════════════
def run_mimic():
    SEP = "=" * 70
    print(f"\n{SEP}\n01_etl.py — MIMIC-IV (v5: risk-set design)\n{SEP}")
    consort = {}

    # ── 1. Small tables → cardiac stays ──────────────────────────
    patients = pd.read_csv(gz(f"{MIMIC_HOSP}/patients.csv.gz"))
    admissions = pd.read_csv(gz(f"{MIMIC_HOSP}/admissions.csv.gz"))
    icustays = pd.read_csv(gz(f"{MIMIC_ICU}/icustays.csv.gz"))
    dx = pd.read_csv(gz(f"{MIMIC_HOSP}/diagnoses_icd.csv.gz"))
    px = pd.read_csv(gz(f"{MIMIC_HOSP}/procedures_icd.csv.gz"))
    presc = pd.read_csv(
        gz(f"{MIMIC_HOSP}/prescriptions.csv.gz"),
        usecols=["subject_id", "hadm_id", "drug", "starttime", "stoptime", "route"],
    )
    consort["total_icu"] = len(icustays)

    px["icd_code"] = px.icd_code.astype(str).str.strip()
    dx["icd_code"] = dx.icd_code.astype(str).str.strip()

    icu = icustays.merge(
        patients[["subject_id", "gender", "anchor_age"]], on="subject_id"
    )
    icu = icu.merge(
        admissions[["hadm_id", "admittime", "dischtime", "hospital_expire_flag"]],
        on="hadm_id",
    )
    icu["intime"] = pd.to_datetime(icu.intime)
    icu["outtime"] = pd.to_datetime(icu.outtime)

    cardiac_hadms = set()
    for ver, codes in [(9, CABG_ICD9 + VALVE_ICD9), (10, CABG_ICD10 + VALVE_ICD10)]:
        v = px[px.icd_version == ver]
        for c in codes:
            cardiac_hadms |= set(v[v.icd_code.str.startswith(c)].hadm_id)
    cvicu_stays = set(icu[icu.first_careunit == CVICU].stay_id)
    cardiac_stays = set(icu[icu.hadm_id.isin(cardiac_hadms)].stay_id) | cvicu_stays

    cardiac = icu[icu.stay_id.isin(cardiac_stays)].copy()
    cardiac = cardiac[cardiac.anchor_age >= MIN_AGE]
    cardiac = cardiac.sort_values("intime").groupby("subject_id").first().reset_index()
    cardiac["age"] = cardiac.anchor_age
    cardiac["is_female"] = (cardiac.gender == "F").astype(int)

    def classify_mimic(hadm_id):
        codes = set(px[px.hadm_id == hadm_id].icd_code)
        has_cabg = any(codes & set(c) for c in [CABG_ICD9, CABG_ICD10])
        has_valve = any(codes & set(c) for c in [VALVE_ICD9, VALVE_ICD10])
        if has_cabg and has_valve:
            return "combined"
        if has_cabg:
            return "cabg"
        if has_valve:
            return "valve"
        return "other_cardiac"

    cardiac["surgery_type"] = cardiac.hadm_id.apply(classify_mimic)

    stays = set(cardiac.stay_id)
    hadms = set(cardiac.hadm_id.dropna().astype(int))
    consort["cardiac_adult_first"] = len(cardiac)
    print(f"  Cardiac surgery, adult, 1st stay: {len(cardiac):,}")

    # ESKD
    eskd_hadms = set()
    for pat in ESKD_PATTERNS:
        eskd_hadms |= set(
            dx[
                dx.hadm_id.isin(hadms)
                & dx.icd_code.str.lower().str.contains(pat.replace(" ", ""), na=False)
            ].hadm_id
        )
    for ver, codes in ESKD_ICD.items():
        v = dx[(dx.hadm_id.isin(hadms)) & (dx.icd_version == ver)]
        for c in codes:
            eskd_hadms |= set(v[v.icd_code.str.startswith(c)].hadm_id)
    eskd_stays = set(cardiac[cardiac.hadm_id.isin(eskd_hadms)].stay_id)
    cardiac = cardiac[~cardiac.stay_id.isin(eskd_stays)]
    stays = set(cardiac.stay_id)
    hadms = set(cardiac.hadm_id.dropna().astype(int))
    consort["post_eskd"] = len(cardiac)
    print(f"  ESKD excluded: {len(eskd_stays):,} → {len(cardiac):,}")

    # ── 2. DuckDB big tables ─────────────────────────────────────
    item_list = ",".join(str(i) for i in ALL_LAB_ITEMS_MIMIC)
    print(f"\n  Loading big tables via DuckDB ({len(hadms):,} hadm_ids)...")
    labs = load_filtered(
        "labevents",
        MIMIC_HOSP,
        hadms,
        pid_col="hadm_id",
        extra_where=f"itemid IN ({item_list})",
    )
    ie = load_filtered("inputevents", MIMIC_ICU, stays, pid_col="stay_id")
    hr_items = ",".join(str(i) for i in VITAL_HR_MIMIC)
    ce_hr = load_filtered(
        "chartevents",
        MIMIC_ICU,
        stays,
        pid_col="stay_id",
        extra_where=f"itemid IN ({hr_items})",
    )

    # ── 3. IV Mg identification ──────────────────────────────────
    print("\n  IV Magnesium identification...")
    ie_mg = ie[(ie.stay_id.isin(stays)) & (ie.itemid.isin(MG_SUPP_ITEMS_MIMIC))].copy()
    if "statusdescription" in ie_mg.columns:
        ie_mg = ie_mg[~ie_mg.statusdescription.str.contains("Rewritten", na=False)]
    ie_mg = ie_mg[ie_mg.amount.notna() & (ie_mg.amount > 0)]
    ie_mg["starttime"] = pd.to_datetime(ie_mg.starttime)
    ie_mg = ie_mg.merge(cardiac[["stay_id", "intime"]], on="stay_id")
    ie_mg["mg_offset_h"] = (ie_mg.starttime - ie_mg.intime).dt.total_seconds() / 3600
    ie_mg = ie_mg[ie_mg.mg_offset_h >= 0]
    first_mg = ie_mg.sort_values("mg_offset_h").groupby("stay_id").first().reset_index()
    first_mg["mg_offset_min"] = first_mg.mg_offset_h * 60

    treated_stays = set(first_mg.stay_id)
    control_stays = stays - treated_stays
    consort["treated_any_ivmg"] = len(treated_stays)
    consort["control_no_ivmg"] = len(control_stays)
    print(f"  Treated: {pct(len(treated_stays), len(stays))}")
    desc(first_mg.mg_offset_h, "Onset (h)")
    for c in [6, 12, 24, 48]:
        print(f"    ≤{c}h: {pct((first_mg.mg_offset_h<=c).sum(), len(first_mg))}")

    # ── 4. Creatinine ────────────────────────────────────────────
    cr = labs[
        labs.itemid.isin(LAB_CR_MIMIC) & labs.valuenum.between(CR_MIN, CR_MAX)
    ].copy()
    cr["charttime"] = pd.to_datetime(cr.charttime)
    cr = cr.merge(cardiac[["stay_id", "hadm_id", "intime"]], on="hadm_id")
    cr["offset_h"] = (cr.charttime - cr.intime).dt.total_seconds() / 3600
    cr["offset_min"] = cr.offset_h * 60

    cr_first = (
        cr[cr.offset_h >= 0]
        .sort_values("offset_h")
        .groupby("stay_id")
        .first()
        .reset_index()
    )

    # Treated: Cr before Mg
    cr_t = cr[cr.stay_id.isin(treated_stays)].merge(
        first_mg[["stay_id", "mg_offset_h", "mg_offset_min"]], on="stay_id"
    )
    has_cr = set(
        cr_t[(cr_t.offset_h >= 0) & (cr_t.offset_min < cr_t.mg_offset_min)].stay_id
    )
    consort["treated_has_cr_pre"] = len(has_cr)
    consort["treated_no_cr_pre"] = len(treated_stays) - len(has_cr)

    # Hosp Cr fallback
    adm_times = admissions[["hadm_id", "admittime"]].copy()
    adm_times.admittime = pd.to_datetime(adm_times.admittime)
    cr_hadm = cr.merge(adm_times, on="hadm_id", how="left")
    cr_hadm["pre_adm"] = cr_hadm.charttime < cr_hadm.admittime
    hosp_cr = (
        cr_hadm[cr_hadm.pre_adm & cr_hadm.stay_id.isin(treated_stays - has_cr)]
        .sort_values("charttime", ascending=False)
        .groupby("stay_id")
        .first()
        .reset_index()
    )
    rescue = set(hosp_cr.stay_id)
    consort["treated_hosp_fallback"] = len(rescue)
    has_cr |= rescue
    consort["treated_has_cr_pre_with_fallback"] = len(has_cr)

    high_cr = set(cr_first[cr_first.valuenum >= BASELINE_CR_MAX].stay_id)
    consort["excl_cr_high"] = len(high_cr & has_cr)
    keep_treated = (has_cr & treated_stays) - high_cr
    consort["treated_final"] = len(keep_treated)

    cr_ctrl = cr[(cr.stay_id.isin(control_stays)) & (cr.offset_h >= 0)]
    c2 = set(cr_ctrl.groupby("stay_id").size().pipe(lambda s: s[s >= 2]).index)
    c_excl = set(
        cr_first[
            cr_first.stay_id.isin(c2) & (cr_first.valuenum >= BASELINE_CR_MAX)
        ].stay_id
    )
    keep_control = c2 - c_excl
    consort["control_has_2cr"] = len(c2)
    consort["control_final"] = len(keep_control)

    # ── 5. Build treated ─────────────────────────────────────────
    print(f"\n  Building treated: {len(keep_treated):,}")
    treated = cardiac[cardiac.stay_id.isin(keep_treated)].copy()
    treated = treated.merge(
        first_mg[["stay_id", "mg_offset_min", "mg_offset_h"]], on="stay_id"
    )
    treated["mg_starttime"] = treated.intime + pd.to_timedelta(
        treated.mg_offset_h, unit="h"
    )

    treated = treated.merge(
        cr_first[["stay_id", "valuenum"]].rename(columns={"valuenum": "first_cr"}),
        on="stay_id",
        how="left",
    )
    treated["egfr"] = compute_egfr(treated.first_cr, treated.age, treated.is_female)

    for como, code_map in MIMIC_COMORB_ICD.items():
        treated[como] = treated.hadm_id.isin(
            matches_icd(dx, set(treated.hadm_id.dropna().astype(int)), code_map)
        ).astype(int)

    # Chronic drugs
    presc["starttime"] = pd.to_datetime(presc.starttime, errors="coerce")
    adm_t = admissions[["hadm_id", "admittime"]].copy()
    adm_t.admittime = pd.to_datetime(adm_t.admittime)
    p = presc[presc.hadm_id.isin(set(treated.hadm_id.dropna().astype(int)))].merge(
        adm_t, on="hadm_id"
    )
    p = p.merge(
        treated[["hadm_id", "stay_id", "intime"]].drop_duplicates("hadm_id"),
        on="hadm_id",
        how="left",
    )
    p_chronic = p[
        (p.starttime < p.admittime)
        | (
            (p.starttime < p.intime)
            & p.route.str.lower().str.contains(ORAL_ROUTE_RE, na=False)
        )
    ]
    for col, pats in CHRONIC_DRUG_CLASSES.items():
        hits = set(
            p_chronic[
                p_chronic.drug.str.lower().str.contains(
                    "|".join(x.lower() for x in pats), na=False
                )
            ].stay_id
        )
        treated[col] = treated.stay_id.isin(hits).astype(int)

    # BMI from chartevents
    bmi_items = {226707, 226730}  # height, weight
    ce_bmi = load_filtered(
        "chartevents",
        MIMIC_ICU,
        keep_treated,
        pid_col="stay_id",
        extra_where=f"itemid IN ({','.join(str(i) for i in bmi_items)})",
    )
    if len(ce_bmi) > 0:
        ce_bmi["charttime"] = pd.to_datetime(ce_bmi.charttime)
        ce_bmi = ce_bmi.merge(treated[["stay_id", "intime"]], on="stay_id")
        # Take first height and weight
        for item, col in [(226707, "_ht"), (226730, "_wt")]:
            sub = (
                ce_bmi[ce_bmi.itemid == item]
                .sort_values("charttime")
                .groupby("stay_id")
                .first()
                .reset_index()
            )
            treated = treated.merge(
                sub[["stay_id", "valuenum"]].rename(columns={"valuenum": col}),
                on="stay_id",
                how="left",
            )
        if "_ht" in treated.columns and "_wt" in treated.columns:
            treated["bmi"] = treated._wt / ((treated._ht / 100) ** 2)
            treated.loc[~treated.bmi.between(10, 80), "bmi"] = np.nan
            treated.drop(columns=["_ht", "_wt"], inplace=True, errors="ignore")

    treated["icu_discharge_h"] = (
        treated.outtime - treated.intime
    ).dt.total_seconds() / 3600
    treated["icu_outcome"] = treated.hospital_expire_flag.map(
        {1: "expired", 0: "alive"}
    ).fillna("unknown")
    treated["hosp_mortality"] = treated.hospital_expire_flag.fillna(0).astype(int)

    # ── 6. Build control ─────────────────────────────────────────
    print(f"  Building control: {len(keep_control):,}")
    control = cardiac[cardiac.stay_id.isin(keep_control)].copy()

    control = control.merge(
        cr_first[["stay_id", "valuenum"]].rename(columns={"valuenum": "first_cr"}),
        on="stay_id",
        how="left",
    )
    control["egfr"] = compute_egfr(control.first_cr, control.age, control.is_female)

    for como, code_map in MIMIC_COMORB_ICD.items():
        control[como] = control.hadm_id.isin(
            matches_icd(dx, set(control.hadm_id.dropna().astype(int)), code_map)
        ).astype(int)

    p_c = presc[presc.hadm_id.isin(set(control.hadm_id.dropna().astype(int)))].merge(
        adm_t, on="hadm_id"
    )
    p_c = p_c.merge(
        control[["hadm_id", "stay_id", "intime"]].drop_duplicates("hadm_id"),
        on="hadm_id",
        how="left",
    )
    pc_chr = p_c[
        (p_c.starttime < p_c.admittime)
        | (
            (p_c.starttime < p_c.intime)
            & p_c.route.str.lower().str.contains(ORAL_ROUTE_RE, na=False)
        )
    ]
    for col, pats in CHRONIC_DRUG_CLASSES.items():
        hits = set(
            pc_chr[
                pc_chr.drug.str.lower().str.contains(
                    "|".join(x.lower() for x in pats), na=False
                )
            ].stay_id
        )
        control[col] = control.stay_id.isin(hits).astype(int)

    # BMI for controls
    ce_bmi_c = load_filtered(
        "chartevents",
        MIMIC_ICU,
        keep_control,
        pid_col="stay_id",
        extra_where=f"itemid IN ({','.join(str(i) for i in bmi_items)})",
    )
    if len(ce_bmi_c) > 0:
        ce_bmi_c["charttime"] = pd.to_datetime(ce_bmi_c.charttime)
        for item, col in [(226707, "_ht"), (226730, "_wt")]:
            sub = (
                ce_bmi_c[ce_bmi_c.itemid == item]
                .sort_values("charttime")
                .groupby("stay_id")
                .first()
                .reset_index()
            )
            control = control.merge(
                sub[["stay_id", "valuenum"]].rename(columns={"valuenum": col}),
                on="stay_id",
                how="left",
            )
        if "_ht" in control.columns and "_wt" in control.columns:
            control["bmi"] = control._wt / ((control._ht / 100) ** 2)
            control.loc[~control.bmi.between(10, 80), "bmi"] = np.nan
            control.drop(columns=["_ht", "_wt"], inplace=True, errors="ignore")

    control["icu_discharge_h"] = (
        control.outtime - control.intime
    ).dt.total_seconds() / 3600
    control["icu_outcome"] = control.hospital_expire_flag.map(
        {1: "expired", 0: "alive"}
    ).fillna("unknown")
    control["hosp_mortality"] = control.hospital_expire_flag.fillna(0).astype(int)

    # ── 7. Secondary outcomes ────────────────────────────────────
    print("\n  Secondary outcomes...")
    all_hadms = set(treated.hadm_id.dropna().astype(int)) | set(
        control.hadm_id.dropna().astype(int)
    )
    all_sids = set(treated.subject_id) | set(control.subject_id)

    # POAF
    prior_af_sids = set()
    all_dx = dx[dx.subject_id.isin(all_sids)]
    af_hadms = set()
    for code in AF_ICD9:
        af_hadms |= set(
            all_dx[
                (all_dx.icd_version == 9) & all_dx.icd_code.str.startswith(code)
            ].hadm_id
        )
    af_hadms |= set(
        all_dx[
            (all_dx.icd_version == 10) & all_dx.icd_code.str.startswith(AF_ICD10_PREFIX)
        ].hadm_id
    )
    hadm_sid = dict(
        zip(
            pd.concat(
                [treated[["hadm_id", "subject_id"]], control[["hadm_id", "subject_id"]]]
            ).hadm_id,
            pd.concat(
                [treated[["hadm_id", "subject_id"]], control[["hadm_id", "subject_id"]]]
            ).subject_id,
        )
    )
    for sid in all_sids:
        sid_hadms = set(all_dx[all_dx.subject_id == sid].hadm_id)
        current = {h for h, s in hadm_sid.items() if s == sid}
        if (sid_hadms - current) & af_hadms:
            prior_af_sids.add(sid)
    af_current = matches_icd(dx, all_hadms, {9: AF_ICD9, 10: [AF_ICD10_PREFIX]})
    poaf_hadms = af_current - {h for h, s in hadm_sid.items() if s in prior_af_sids}
    for df in [treated, control]:
        df["poaf"] = df.hadm_id.isin(poaf_hadms).astype(int)
        df["prior_af"] = df.subject_id.isin(prior_af_sids).astype(int)
    print(f"    POAF: trt={treated.poaf.sum()}, ctl={control.poaf.sum()}")

    enc_hadms = matches_icd(dx, all_hadms, {9: ENCEPH_ICD9, 10: ENCEPH_ICD10_PREFIX})
    for df in [treated, control]:
        df["encephalopathy"] = df.hadm_id.isin(enc_hadms).astype(int)

    varr_hadms = matches_icd(dx, all_hadms, {9: VARR_ICD9, 10: VARR_ICD10})
    for df in [treated, control]:
        df["vent_arrhythmia"] = df.hadm_id.isin(varr_hadms).astype(int)

    # ── 8. Export ────────────────────────────────────────────────
    save_all_patients(treated, control, "mimic")

    # did_cr_all
    all_s = keep_treated | keep_control
    cr_exp = cr[
        cr.stay_id.isin(all_s)
        & (cr.offset_h >= 0)
        & (cr.valuenum <= CR_POST_PLAUSIBLE_MAX)
    ]
    cr_exp = cr_exp[["stay_id", "valuenum", "offset_min", "offset_h"]].copy()
    cr_exp = cr_exp.rename(
        columns={"valuenum": "labresult", "offset_min": "labresultoffset"}
    )
    cr_exp.to_csv(os.path.join(RESULTS, "did_cr_all_mimic.csv"), index=False)
    print(
        f"  ✓ did_cr_all_mimic.csv: {len(cr_exp):,} Cr for {cr_exp.stay_id.nunique():,} pts"
    )

    # did_labs_all
    lab_rows = []
    mimic_lab_map = {
        "magnesium": LAB_MG_MIMIC,
        "potassium": LAB_K_MIMIC,
        "calcium": LAB_CA_MIMIC,
        "lactate": LAB_LAC_MIMIC,
    }
    for lab_name, items in mimic_lab_map.items():
        sub = labs[
            labs.itemid.isin(items)
            & labs.hadm_id.isin(
                set(cardiac[cardiac.stay_id.isin(all_s)].hadm_id.dropna().astype(int))
            )
        ]
        if len(sub) > 0:
            sub = sub.copy()
            sub["charttime"] = pd.to_datetime(sub.charttime)
            sub = sub.merge(cardiac[["stay_id", "hadm_id", "intime"]], on="hadm_id")
            sub["offset_h"] = (sub.charttime - sub.intime).dt.total_seconds() / 3600
            sub = sub[sub.offset_h >= 0]
            lab_rows.append(
                pd.DataFrame(
                    {
                        "stay_id": sub.stay_id,
                        "lab_name": lab_name,
                        "value": sub.valuenum,
                        "offset_h": sub.offset_h,
                    }
                )
            )

    # Heart rate
    if len(ce_hr) > 0:
        hr = ce_hr[ce_hr.stay_id.isin(all_s)].copy()
        hr["charttime"] = pd.to_datetime(hr.charttime)
        hr = hr.merge(cardiac[["stay_id", "intime"]], on="stay_id")
        hr["offset_h"] = (hr.charttime - hr.intime).dt.total_seconds() / 3600
        hr = hr[(hr.offset_h >= 0) & hr.valuenum.between(20, 250)]
        if len(hr) > 0:
            lab_rows.append(
                pd.DataFrame(
                    {
                        "stay_id": hr.stay_id,
                        "lab_name": "heartrate",
                        "value": hr.valuenum,
                        "offset_h": hr.offset_h,
                    }
                )
            )

    labs_all = pd.concat(lab_rows, ignore_index=True) if lab_rows else pd.DataFrame()
    labs_all.to_csv(os.path.join(RESULTS, "did_labs_all_mimic.csv"), index=False)
    print(f"  ✓ did_labs_all_mimic.csv: {len(labs_all):,} measurements")
    for ln in labs_all.lab_name.unique():
        n = (labs_all.lab_name == ln).sum()
        np_ = labs_all[labs_all.lab_name == ln].stay_id.nunique()
        print(f"    {ln}: {n:,} values across {np_:,} pts")

    consort["db"] = "MIMIC"
    print(f"\n  CONSORT: {consort}")
    return consort


# ═══════════════════════════════════════════════════════════════════
#  MAIN
# ═══════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    print("=" * 70)
    print("01_etl.py — Cohort Construction (v5: risk-set design)")
    print("  Outputs: did_all, did_labs_all, did_cr_all per database")
    print("  No lab summary columns — R computes at match time")
    print("=" * 70)

    args = [a.lower() for a in sys.argv[1:]]
    consorts = []
    if not args or "eicu" in args:
        c = run_eicu()
        if c:
            consorts.append(c)
    if not args or "mimic" in args:
        c = run_mimic()
        if c:
            consorts.append(c)

    if consorts:
        pd.DataFrame(consorts).to_csv(
            os.path.join(RESULTS, "did_consort.csv"), index=False
        )
        print(f"\n  ✓ did_consort.csv saved")

    print("\n" + "=" * 70)
    print("NEXT: Rscript 02_psm.R eicu / mimic")
    print("=" * 70)
