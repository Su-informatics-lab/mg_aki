#!/usr/bin/env python3
"""
Mg Reserve → Cardiac Surgery AKI  (eICU-CRD)
01_etl.py — cohort construction, exposure, outcome, covariates

Run:  python 01_etl.py                   # uses EICU_DATA env or demo default
      EICU_DATA=/path/to/full python 01_etl.py   # full dataset

Design anchored to:
  - TTE skill: ACNU design for Analysis B (Mg supplementation)
  - Eadon baseline Cr: lowest within 48h (sensitivity)
  - Approach 2 baseline Cr: pre-ICU if available, else first post (primary)
  - Temporal anchoring: AKI onset MUST follow Mg measurement
"""

import os
import warnings

warnings.filterwarnings("ignore")

from importlib.util import module_from_spec, spec_from_file_location

import numpy as np
import pandas as pd

# Load config from same directory
_cfg_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "00_config.py")
_spec = spec_from_file_location("config", _cfg_path)
cfg = module_from_spec(_spec)
_spec.loader.exec_module(cfg)


# =====================================================================
# HELPERS
# =====================================================================
def save(df, name):
    path = os.path.join(cfg.RESULTS, name)
    df.to_csv(path, index=False)
    print(f"  → {path}  ({len(df):,} × {df.shape[1]})")


def load_csv(table_name, usecols=None):
    """Load eICU CSV, trying common naming patterns and casing."""
    candidates = [
        f"{table_name}.csv",
        f"{table_name}.csv.gz",
        f"{table_name.lower()}.csv",
        f"{table_name.lower()}.csv.gz",
    ]
    for pattern in candidates:
        path = os.path.join(cfg.DATA_ROOT, pattern)
        if os.path.exists(path):
            # If usecols given, lowercase them to match actual headers
            uc = [c.lower() for c in usecols] if usecols else None
            df = pd.read_csv(path, low_memory=False, usecols=uc)
            df.columns = df.columns.str.lower()
            print(f"  {table_name}: {len(df):,} rows")
            return df
    print(f"  WARNING: {table_name} not found in {cfg.DATA_ROOT}")
    return pd.DataFrame()


def age_numeric(age_str):
    """Parse eICU age (varchar) → numeric. '> 89' → 90."""
    if pd.isna(age_str):
        return np.nan
    s = str(age_str).strip()
    if s.startswith(">"):
        return cfg.AGE_CAP
    try:
        return float(s)
    except ValueError:
        return np.nan


def matches_any(series, patterns):
    """Case-insensitive check if series contains any pattern."""
    pat = "|".join(patterns)
    return series.str.lower().str.contains(pat, na=False)


def compute_egfr(cr, age, is_female):
    """CKD-EPI 2021 race-free."""
    cr = np.asarray(cr, dtype=np.float64)
    age = np.asarray(age, dtype=np.float64)
    fem = np.asarray(is_female, dtype=bool)
    kappa = np.where(fem, 0.7, 0.9)
    alpha = np.where(fem, -0.241, -0.302)
    ratio = cr / kappa
    return (
        142
        * np.power(np.minimum(ratio, 1.0), alpha)
        * np.power(np.maximum(ratio, 1.0), -1.200)
        * np.power(0.9938, age)
        * np.where(fem, 1.012, 1.0)
    )


# =====================================================================
# STEP 0: Load tables
# =====================================================================
def load_tables():
    print("=" * 70)
    print("STEP 0: Loading eICU tables")
    print("=" * 70)
    print(f"  DATA_ROOT: {cfg.DATA_ROOT}")

    # Check for zip that needs extraction
    zip_path = os.path.join(
        os.path.dirname(cfg.DATA_ROOT),
        "eicu-collaborative-research-database-demo-2.0.1.zip",
    )
    if not os.path.isdir(cfg.DATA_ROOT) and os.path.exists(zip_path):
        print(f"  Extracting {zip_path}...")
        import zipfile

        with zipfile.ZipFile(zip_path, "r") as z:
            z.extractall(os.path.dirname(cfg.DATA_ROOT))

    tables = {}
    tables["patient"] = load_csv("patient")
    tables["lab"] = load_csv("lab")
    tables["medication"] = load_csv(
        "medication",
        usecols=[
            "patientunitstayid",
            "medicationid",
            "drugstartoffset",
            "drugstopoffset",
            "drugname",
            "drugordercancelled",
            "routeadmin",
            "dosage",
        ],
    )
    tables["infusionDrug"] = load_csv("infusionDrug")
    tables["diagnosis"] = load_csv("diagnosis")
    tables["pastHistory"] = load_csv("pastHistory")
    tables["treatment"] = load_csv("treatment")
    tables["admissionDrug"] = load_csv("admissionDrug")
    tables["hospital"] = load_csv("hospital")

    # APACHE tables
    tables["apachePatientResult"] = load_csv("apachePatientResult")
    tables["apacheApsVar"] = load_csv("apacheApsVar")
    tables["apachePredVar"] = load_csv("apachePredVar")
    tables["intakeOutput"] = load_csv("intakeOutput")

    return tables


# =====================================================================
# STEP 1: Cardiac surgery cohort
# =====================================================================
def build_cardiac_cohort(t):
    print("\n" + "=" * 70)
    print("STEP 1: Cardiac Surgery Cohort")
    print("=" * 70)
    consort = {}
    pt = t["patient"].copy()
    consort["total_icu_stays"] = len(pt)
    print(f"  Total ICU stays: {len(pt):,}")

    # Parse age
    pt["age_num"] = pt["age"].apply(age_numeric)
    adults = pt[pt.age_num >= cfg.MIN_AGE].copy()
    consort["adults"] = len(adults)
    print(f"  Adults (≥{cfg.MIN_AGE}): {len(adults):,}")

    # Cardiac surgery: apacheadmissiondx OR unittype
    by_dx = matches_any(adults.apacheadmissiondx, cfg.CARDIAC_DX_PATTERNS)
    by_unit = adults.unittype.isin(cfg.CARDIAC_UNIT_TYPES)
    cardiac = adults[by_dx | by_unit].copy()
    consort["cardiac_surgery"] = len(cardiac)
    print(f"  Cardiac surgery (dx|unittype): {len(cardiac):,}")
    print(f"    by apacheadmissiondx: {by_dx.sum():,}")
    print(f"    by unittype:          {by_unit.sum():,}")

    # First ICU stay per patient
    # Larger hospitaladmitoffset → closer to hospital admit → earlier stay
    cardiac = (
        cardiac.sort_values("hospitaladmitoffset", ascending=False)
        .groupby("uniquepid")
        .first()
        .reset_index()
    )
    consort["first_stay"] = len(cardiac)
    print(f"  First ICU stay per patient: {len(cardiac):,}")

    # Classify surgery type
    def _surgery_type(dx):
        if pd.isna(dx):
            return "unknown"
        dx = str(dx).lower()
        has_cabg = any(k in dx for k in cfg.SURGERY_TYPE_MAP["cabg"])
        has_valve = any(k in dx for k in cfg.SURGERY_TYPE_MAP["valve"])
        if has_cabg and has_valve:
            return "combined"
        if has_cabg:
            return "cabg"
        if has_valve:
            return "valve"
        return "other_cardiac"

    cardiac["surgery_type"] = cardiac.apacheadmissiondx.apply(_surgery_type)
    print(f"  Surgery types: {cardiac.surgery_type.value_counts().to_dict()}")

    return cardiac, consort


# =====================================================================
# STEP 2: Mg exposure (first Mg within window)
# =====================================================================
def extract_mg_exposure(t, cardiac):
    print("\n" + "=" * 70)
    print("STEP 2: Magnesium Exposure")
    print("=" * 70)
    lab = t["lab"]
    pids = set(cardiac.patientunitstayid)

    # All Mg labs for cardiac patients
    mg = lab[
        lab.patientunitstayid.isin(pids)
        & matches_any(lab.labname, cfg.MG_LABNAMES)
        & lab.labresult.between(cfg.MG_PLAUSIBLE_MIN, cfg.MG_PLAUSIBLE_MAX)
    ].copy()
    print(f"  Mg measurements (cardiac patients): {len(mg):,}")
    print(f"  Patients with any Mg: {mg.patientunitstayid.nunique():,}")

    # First Mg within window [0, MG_WINDOW_MIN]
    mg_window = mg[
        (mg.labresultoffset >= 0) & (mg.labresultoffset <= cfg.MG_WINDOW_MIN)
    ]
    first_mg = (
        mg_window.sort_values("labresultoffset")
        .groupby("patientunitstayid")
        .first()
        .reset_index()[["patientunitstayid", "labresult", "labresultoffset"]]
        .rename(
            columns={
                "labresult": "first_mg_value",
                "labresultoffset": "mg_offset",
            }
        )
    )
    first_mg["mg_hours"] = first_mg.mg_offset / 60.0
    print(f"  Patients with Mg within {cfg.MG_WINDOW_HOURS}h: {len(first_mg):,}")
    print(
        f"  Mg: median={first_mg.first_mg_value.median():.1f}, "
        f"IQR=[{first_mg.first_mg_value.quantile(.25):.1f}, "
        f"{first_mg.first_mg_value.quantile(.75):.1f}]"
    )

    # Quartiles
    first_mg["mg_quartile"] = pd.qcut(
        first_mg.first_mg_value, 4, labels=["Q1", "Q2", "Q3", "Q4"]
    )
    # Clinical categories
    first_mg["mg_category"] = pd.cut(
        first_mg.first_mg_value,
        bins=[0, 1.8, 2.3, 999],
        labels=["hypo", "normal", "hyper"],
    )
    print(f"  Mg categories: {first_mg.mg_category.value_counts().to_dict()}")

    return first_mg


# =====================================================================
# STEP 3: Baseline creatinine
# =====================================================================
def compute_baseline_cr(t, eligible_pids, mg_offsets):
    """
    Primary (Approach 2): lowest Cr with offset ∈ [-720, 0].
      Fallback: first Cr with offset > 0 and < mg_offset.
    Sensitivity (Approach 1 / Eadon): lowest Cr within [-60, 2880].

    mg_offsets: dict {patientunitstayid: mg_offset} for temporal anchoring.
    """
    print("\n" + "=" * 70)
    print("STEP 3: Baseline Creatinine")
    print("=" * 70)
    lab = t["lab"]

    cr = lab[
        lab.patientunitstayid.isin(eligible_pids)
        & matches_any(lab.labname, cfg.CR_LABNAMES)
        & lab.labresult.between(cfg.CR_PLAUSIBLE_MIN, cfg.CR_PLAUSIBLE_MAX)
    ].copy()
    print(f"  Cr measurements (eligible): {len(cr):,}")

    # ── Primary: Approach 2 ──────────────────────────────────────
    # Pre-ICU Cr: offset ∈ [-720, 0]
    pre_icu = cr[
        (cr.labresultoffset >= cfg.BASELINE_PRE_ICU_WINDOW_MIN)
        & (cr.labresultoffset <= cfg.BASELINE_PRE_ICU_WINDOW_MAX)
    ]
    pre_icu_bl = (
        pre_icu.sort_values("labresult")
        .groupby("patientunitstayid")
        .first()
        .reset_index()[["patientunitstayid", "labresult"]]
        .rename(columns={"labresult": "baseline_cr"})
    )
    pre_icu_bl["baseline_source"] = "pre_icu"
    n_pre = len(pre_icu_bl)

    # Fallback: first Cr > 0 but BEFORE Mg measurement
    pids_need_fallback = eligible_pids - set(pre_icu_bl.patientunitstayid)
    if pids_need_fallback:
        fb_pool = cr[
            cr.patientunitstayid.isin(pids_need_fallback) & (cr.labresultoffset > 0)
        ].copy()
        # Only keep Cr measured BEFORE the Mg lab for each patient
        fb_rows = []
        for pid, grp in fb_pool.groupby("patientunitstayid"):
            mg_off = mg_offsets.get(pid, cfg.MG_WINDOW_MIN)
            before_mg = grp[grp.labresultoffset < mg_off]
            if len(before_mg) > 0:
                fb_rows.append(before_mg.sort_values("labresultoffset").iloc[0])
        if fb_rows:
            fb_df = pd.DataFrame(fb_rows)[["patientunitstayid", "labresult"]].rename(
                columns={"labresult": "baseline_cr"}
            )
            fb_df["baseline_source"] = "first_post_icu"
        else:
            fb_df = pd.DataFrame(
                columns=["patientunitstayid", "baseline_cr", "baseline_source"]
            )
    else:
        fb_df = pd.DataFrame(
            columns=["patientunitstayid", "baseline_cr", "baseline_source"]
        )
    n_fb = len(fb_df)

    baseline_primary = pd.concat([pre_icu_bl, fb_df], ignore_index=True)
    print(f"  PRIMARY baseline:")
    print(f"    Pre-ICU (lowest in [-12h, 0]): {n_pre:,}")
    print(f"    Fallback (first post-ICU < Mg): {n_fb:,}")
    print(f"    Total with baseline: {len(baseline_primary):,}")

    # ── Sensitivity: Approach 1 (Eadon) ──────────────────────────
    eadon_pool = cr[
        (cr.labresultoffset >= cfg.BASELINE_EADON_WINDOW_MIN)
        & (cr.labresultoffset <= cfg.BASELINE_EADON_WINDOW_MAX)
    ]
    baseline_eadon = (
        eadon_pool.sort_values("labresult")
        .groupby("patientunitstayid")
        .first()
        .reset_index()[["patientunitstayid", "labresult"]]
        .rename(columns={"labresult": "baseline_cr_eadon"})
    )
    print(f"  EADON baseline (lowest in [-1h, 48h]): {len(baseline_eadon):,}")

    return baseline_primary, baseline_eadon, cr


# =====================================================================
# STEP 4: AKI phenotype with temporal anchoring
# =====================================================================
def apply_aki_phenotype(cohort, cr_all, mg_offsets):
    """
    AKI defined as Cr ≥ 1.5× baseline AFTER the Mg measurement time,
    within 7 days of ICU admission. Temporal anchoring ensures
    Mg precedes AKI onset → addresses reverse causation.
    """
    print("\n" + "=" * 70)
    print("STEP 4: AKI Phenotyping (temporally anchored)")
    print("=" * 70)
    cohort = cohort.copy()

    results = []
    for _, row in cohort.iterrows():
        pid = row.patientunitstayid
        bl_cr = row.baseline_cr
        mg_off = mg_offsets.get(pid, 0)

        pt_cr = cr_all[cr_all.patientunitstayid == pid].copy()

        # Follow-up Cr: AFTER Mg measurement, within 7 days
        fu = pt_cr[
            (pt_cr.labresultoffset > mg_off)
            & (pt_cr.labresultoffset <= cfg.AKI_WINDOW_7D_MIN)
        ].sort_values("labresultoffset")

        # Primary: Cr ≥ 1.5× baseline
        aki_15x = 0
        aki_time = np.nan
        if len(fu) > 0 and bl_cr > 0:
            fu_ratio = fu.labresult / bl_cr
            hits = fu[fu_ratio >= cfg.AKI_RATIO_STAGE1]
            if len(hits) > 0:
                aki_15x = 1
                aki_time = hits.labresultoffset.iloc[0]

        # Secondary: delta ≥ 0.3 within 48h sliding window
        aki_delta03 = 0
        if len(fu) >= 2:
            for i in range(len(fu)):
                for j in range(i + 1, len(fu)):
                    off_i = fu.labresultoffset.iloc[i]
                    off_j = fu.labresultoffset.iloc[j]
                    if (off_j - off_i) <= cfg.AKI_WINDOW_48H_MIN:
                        delta = fu.labresult.iloc[j] - fu.labresult.iloc[i]
                        if delta >= cfg.AKI_DELTA_48H:
                            aki_delta03 = 1
                            break
                if aki_delta03:
                    break

        # Stage 2/3
        max_cr = fu.labresult.max() if len(fu) > 0 else np.nan
        max_ratio = max_cr / bl_cr if bl_cr > 0 and not np.isnan(max_cr) else 0
        aki_stage2 = int(max_ratio >= cfg.AKI_RATIO_STAGE2)
        aki_stage3 = int(
            max_ratio >= cfg.AKI_RATIO_STAGE3
            or (max_cr >= cfg.AKI_CR_ABSOLUTE if not np.isnan(max_cr) else False)
        )

        # Washout: Cr ≥ 1.5× baseline AT or BEFORE Mg measurement
        pre_mg = pt_cr[(pt_cr.labresultoffset > 0) & (pt_cr.labresultoffset <= mg_off)]
        prevalent_aki = 0
        if len(pre_mg) > 0 and bl_cr > 0:
            if (pre_mg.labresult / bl_cr).max() >= cfg.AKI_RATIO_STAGE1:
                prevalent_aki = 1

        results.append(
            {
                "patientunitstayid": pid,
                "aki_primary": aki_15x,
                "aki_time_offset": aki_time,
                "aki_delta03": aki_delta03,
                "aki_stage2": aki_stage2,
                "aki_stage3": aki_stage3,
                "max_followup_cr": max_cr,
                "max_cr_ratio": max_ratio,
                "n_followup_cr": len(fu),
                "prevalent_aki": prevalent_aki,
            }
        )

    aki_df = pd.DataFrame(results)
    cohort = cohort.merge(aki_df, on="patientunitstayid")

    # Exclude prevalent AKI (washout)
    pre_washout = len(cohort)
    cohort = cohort[cohort.prevalent_aki == 0].copy()
    n_washed = pre_washout - len(cohort)
    print(f"  Washout (prevalent AKI at Mg time): {n_washed}")

    # Time-to-AKI for survival analysis (minutes → hours)
    cohort["time_to_aki_hours"] = cohort.aki_time_offset / 60.0
    # Censoring time: min(unitdischargeoffset, 7d) — from Mg time
    if "unitdischargeoffset" in cohort.columns:
        cohort["censor_offset"] = cohort[["unitdischargeoffset"]].apply(
            lambda r: min(r.unitdischargeoffset, cfg.AKI_WINDOW_7D_MIN),
            axis=1,
        )
        cohort["time_to_event_hours"] = np.where(
            cohort.aki_primary == 1,
            cohort.time_to_aki_hours,
            (cohort.censor_offset - cohort.mg_offset) / 60.0,
        )
    else:
        cohort["time_to_event_hours"] = np.where(
            cohort.aki_primary == 1,
            cohort.time_to_aki_hours,
            cfg.AKI_WINDOW_7D_MIN / 60.0,
        )

    n_aki = cohort.aki_primary.sum()
    print(f"  Eligible after washout: {len(cohort):,}")
    print(f"  AKI (primary, ≥1.5×, post-Mg): {n_aki} ({n_aki/len(cohort)*100:.1f}%)")
    print(f"  AKI delta≥0.3: {cohort.aki_delta03.sum()}")
    print(f"  AKI stage 2: {cohort.aki_stage2.sum()}")
    print(f"  AKI stage 3: {cohort.aki_stage3.sum()}")

    return cohort


# =====================================================================
# STEP 5: Covariates
# =====================================================================
def build_covariates(t, cohort):
    print("\n" + "=" * 70)
    print("STEP 5: Covariates")
    print("=" * 70)
    pids = set(cohort.patientunitstayid)

    # ── Demographics (already in cohort from patient table) ──────
    cohort["sex"] = cohort.gender.map({"Male": "Male", "Female": "Female"}).fillna(
        "Other"
    )
    cohort["is_female"] = (cohort.sex == "Female").astype(int)

    # eGFR
    cohort["baseline_egfr"] = compute_egfr(
        cohort.baseline_cr, cohort.age_num, cohort.is_female == 1
    )
    print(
        f"  eGFR: median={cohort.baseline_egfr.median():.1f}, "
        f"<60: {(cohort.baseline_egfr < 60).sum()}"
    )

    # BMI
    cohort["bmi"] = np.where(
        (cohort.admissionheight > 0) & (cohort.admissionweight > 0),
        cohort.admissionweight / (cohort.admissionheight / 100) ** 2,
        np.nan,
    )

    # ── Comorbidities from pastHistory ───────────────────────────
    ph = t["pastHistory"]
    ph_elig = ph[ph.patientunitstayid.isin(pids)]
    for cmorb, keywords in cfg.COMORBIDITY_KEYWORDS.items():
        flagged = set(
            ph_elig[
                matches_any(ph_elig.pasthistorypath, keywords)
                | matches_any(ph_elig.pasthistoryvalue, keywords)
            ].patientunitstayid
        )
        cohort[f"hx_{cmorb}"] = cohort.patientunitstayid.isin(flagged).astype(int)
        print(f"    hx_{cmorb}: {cohort[f'hx_{cmorb}'].sum()}")

    # ── APACHE from apachePredVar ────────────────────────────────
    apv = t["apachePredVar"]
    if len(apv) > 0:
        apv_elig = apv[apv.patientunitstayid.isin(pids)].copy()
        apache_cols = ["electivesurgery", "dialysis"]
        avail = [c for c in apache_cols if c in apv_elig.columns]
        if avail:
            cohort = cohort.merge(
                apv_elig[["patientunitstayid"] + avail].drop_duplicates(
                    "patientunitstayid"
                ),
                on="patientunitstayid",
                how="left",
            )
            # ESKD exclusion via dialysis flag
            if "dialysis" in cohort.columns:
                pre_eskd = len(cohort)
                cohort = cohort[cohort.dialysis != 1].copy()
                print(
                    f"  Excluded pre-existing dialysis (APACHE): {pre_eskd - len(cohort)}"
                )

    # ── APACHE IV score ──────────────────────────────────────────
    apr = t["apachePatientResult"]
    if len(apr) > 0:
        apr_elig = apr[apr.patientunitstayid.isin(pids)].copy()
        if "apachescore" in apr_elig.columns:
            apache_score = (
                apr_elig.sort_values("apachescore", ascending=False)
                .groupby("patientunitstayid")
                .first()
                .reset_index()[["patientunitstayid", "apachescore"]]
            )
            cohort = cohort.merge(apache_score, on="patientunitstayid", how="left")
            print(
                f"  APACHE IV: median={cohort.apachescore.median():.0f}, "
                f"available={cohort.apachescore.notna().sum()}"
            )

    # ── Nephrotoxins from admissionDrug + medication ─────────────
    ad = t["admissionDrug"]
    med = t["medication"]
    ad_elig = ad[ad.patientunitstayid.isin(pids)] if len(ad) > 0 else pd.DataFrame()
    med_elig = (
        med[med.patientunitstayid.isin(pids) & (med.drugordercancelled != "Yes")]
        if len(med) > 0
        else pd.DataFrame()
    )

    for drug_class, patterns in cfg.NEPHROTOXIN_CLASSES.items():
        flagged = set()
        if len(ad_elig) > 0 and "drugname" in ad_elig.columns:
            flagged |= set(
                ad_elig[matches_any(ad_elig.drugname, patterns)].patientunitstayid
            )
        if len(med_elig) > 0:
            # Pre-Mg medications only
            for pid in pids:
                mg_off = cohort.loc[cohort.patientunitstayid == pid, "mg_offset"]
                if len(mg_off) == 0:
                    continue
                mg_off = mg_off.iloc[0]
                pt_meds = med_elig[
                    (med_elig.patientunitstayid == pid)
                    & (med_elig.drugstartoffset <= mg_off)
                ]
                if len(pt_meds) > 0 and matches_any(pt_meds.drugname, patterns).any():
                    flagged.add(pid)
        cohort[f"nephrotox_{drug_class}"] = cohort.patientunitstayid.isin(
            flagged
        ).astype(int)
        print(f"    nephrotox_{drug_class}: {cohort[f'nephrotox_{drug_class}'].sum()}")

    # ── Mg supplementation (for TTE Analysis B) ──────────────────
    mg_supp_pids = set()
    if len(med_elig) > 0:
        mg_meds = med_elig[matches_any(med_elig.drugname, cfg.MG_SUPP_DRUG_PATTERNS)]
        for pid in pids:
            mg_off = cohort.loc[cohort.patientunitstayid == pid, "mg_offset"]
            if len(mg_off) == 0:
                continue
            mg_off = mg_off.iloc[0]
            pt_mg = mg_meds[
                (mg_meds.patientunitstayid == pid)
                & (mg_meds.drugstartoffset >= mg_off)
                & (mg_meds.drugstartoffset <= mg_off + cfg.MG_SUPP_GRACE_MIN)
            ]
            if len(pt_mg) > 0:
                mg_supp_pids.add(pid)
    cohort["mg_supplementation"] = cohort.patientunitstayid.isin(mg_supp_pids).astype(
        int
    )
    print(
        f"  Mg supplementation within {cfg.MG_SUPP_GRACE_HOURS}h: "
        f"{cohort.mg_supplementation.sum()}"
    )

    # ── RRT detection ────────────────────────────────────────────
    tx = t["treatment"]
    rrt_pids = set()
    if len(tx) > 0:
        tx_elig = tx[tx.patientunitstayid.isin(pids)]
        rrt_pids = set(
            tx_elig[
                matches_any(
                    tx_elig.treatmentstring,
                    ["dialysis", "crrt", "hemodialysis", "cvvh"],
                )
            ].patientunitstayid
        )
    cohort["rrt_7d"] = cohort.patientunitstayid.isin(rrt_pids).astype(int)
    print(f"  RRT within ICU stay: {cohort.rrt_7d.sum()}")

    # ── Mortality ────────────────────────────────────────────────
    cohort["icu_mortality"] = (
        cohort.unitdischargestatus.str.lower() == "expired"
    ).astype(int)
    cohort["hosp_mortality"] = (
        cohort.hospitaldischargestatus.str.lower() == "expired"
    ).astype(int)
    print(f"  ICU mortality: {cohort.icu_mortality.sum()}")
    print(f"  Hospital mortality: {cohort.hosp_mortality.sum()}")

    return cohort


# =====================================================================
# STEP 6: Build TTE cohort (Analysis B)
# =====================================================================
def build_tte_cohort(cohort):
    """
    TTE Analysis B: Among hypomagnesemic patients (Mg < 2.0),
    does Mg supplementation within 6h reduce AKI?

    ACNU design (Hernán & Robins 2016):
    - Time zero: first Mg lab showing < 2.0
    - Treatment: Mg supplementation within grace period
    - Comparator: No Mg supplementation
    """
    print("\n" + "=" * 70)
    print("STEP 6: TTE Cohort (Analysis B — Mg supplementation)")
    print("=" * 70)

    tte = cohort[cohort.first_mg_value < cfg.MG_HYPO_THRESHOLD].copy()
    print(f"  Hypomagnesemia (Mg < {cfg.MG_HYPO_THRESHOLD}): {len(tte):,}")
    print(f"  Treated (Mg supplementation): {tte.mg_supplementation.sum()}")
    print(f"  Untreated: {(tte.mg_supplementation == 0).sum()}")
    if len(tte) > 0:
        print(
            f"  AKI in treated: "
            f"{tte[tte.mg_supplementation==1].aki_primary.sum()} / "
            f"{tte.mg_supplementation.sum()}"
        )
        print(
            f"  AKI in untreated: "
            f"{tte[tte.mg_supplementation==0].aki_primary.sum()} / "
            f"{(tte.mg_supplementation==0).sum()}"
        )

    save(tte, "06_tte_cohort.csv")
    return tte


# =====================================================================
# MAIN
# =====================================================================
def main():
    print("=" * 70)
    print("Mg Reserve → Cardiac Surgery AKI  (eICU-CRD)")
    print("=" * 70)

    t = load_tables()
    cardiac, consort = build_cardiac_cohort(t)

    # Step 2: Mg exposure
    first_mg = extract_mg_exposure(t, cardiac)
    cohort = cardiac.merge(first_mg, on="patientunitstayid")
    consort["has_mg"] = len(cohort)
    mg_offsets = dict(zip(cohort.patientunitstayid, cohort.mg_offset))

    # Step 3: Baseline Cr
    baseline_primary, baseline_eadon, cr_all = compute_baseline_cr(
        t, set(cohort.patientunitstayid), mg_offsets
    )
    cohort = cohort.merge(baseline_primary, on="patientunitstayid")
    cohort = cohort.merge(baseline_eadon, on="patientunitstayid", how="left")
    consort["has_baseline_cr"] = len(cohort)

    # Exclude baseline Cr ≥ 4.0
    pre_cr_excl = len(cohort)
    cohort = cohort[cohort.baseline_cr < cfg.BASELINE_CR_MAX].copy()
    consort["excluded_high_cr"] = pre_cr_excl - len(cohort)
    print(
        f"  Excluded baseline Cr ≥ {cfg.BASELINE_CR_MAX}: {pre_cr_excl - len(cohort)}"
    )

    # ESKD exclusion from pastHistory
    ph = t["pastHistory"]
    eskd_pids = set()
    if len(ph) > 0:
        eskd_pids = set(
            ph[
                ph.patientunitstayid.isin(cohort.patientunitstayid)
                & (
                    matches_any(ph.pasthistorypath, cfg.ESKD_DX_PATTERNS)
                    | matches_any(ph.pasthistoryvalue, cfg.ESKD_DX_PATTERNS)
                )
            ].patientunitstayid
        )
    pre_eskd = len(cohort)
    cohort = cohort[~cohort.patientunitstayid.isin(eskd_pids)].copy()
    consort["excluded_eskd_hx"] = pre_eskd - len(cohort)
    print(f"  Excluded ESKD (pastHistory): {pre_eskd - len(cohort)}")

    # Step 4: AKI
    cohort = apply_aki_phenotype(cohort, cr_all, mg_offsets)
    consort["eligible_final"] = len(cohort)
    consort["aki_cases"] = int(cohort.aki_primary.sum())
    consort["aki_controls"] = int((cohort.aki_primary == 0).sum())

    # Step 5: Covariates
    cohort = build_covariates(t, cohort)
    consort["final_n"] = len(cohort)

    # ── Save outputs ─────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("Saving outputs")
    print("=" * 70)

    # CONSORT
    consort_df = pd.DataFrame([{"step": k, "n": int(v)} for k, v in consort.items()])
    save(consort_df, "00_consort.csv")

    # Full cohort (Analysis A + covariates)
    save(cohort, "01_analysis_a_cohort.csv")

    # Step 6: TTE cohort (Analysis B)
    tte = build_tte_cohort(cohort)

    # Column summary
    print("\n  CONSORT flowchart:")
    for _, row in consort_df.iterrows():
        print(f"    {row['step']:35s} {row['n']:>8,}")

    print(f"\n  Final cohort: {len(cohort):,} rows, {cohort.shape[1]} cols")
    print(f"  Analysis A (prognostic): {len(cohort):,}")
    print(f"  Analysis B (TTE, Mg<{cfg.MG_HYPO_THRESHOLD}): {len(tte):,}")
    print(f"\n  Next: Rscript 02_psm.R && Rscript 03_models.R")


if __name__ == "__main__":
    main()
