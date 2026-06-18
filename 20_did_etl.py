#!/usr/bin/env python3
"""
20_did_etl.py — DiD cohort construction for Mg → AKI study (v2)

  Design:
    - Time zero = first postop IV Mg administration (patient-specific)
    - Cr_pre  = latest Cr between ICU admission and IV Mg (no hosp fallback)
    - Cr_post = first Cr after IV Mg at three windows
    - Outcome = ΔCr (continuous)
    - Controls = patients who never received postop IV Mg

  Outputs per database:
    20_did_treated_{db}.csv   — treated arm
    20_did_control_{db}.csv   — control arm
    20_did_cr_all_{db}.csv    — all Cr for temporal matching
    20_did_consort.csv        — CONSORT numbers

  Run:  python 20_did_etl.py          # both
        python 20_did_etl.py eicu     # eICU only
        python 20_did_etl.py mimic    # MIMIC only
"""

import os
import sys
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

# ═══════════════════════════════════════════════════════════════════
# PATHS
# ═══════════════════════════════════════════════════════════════════
RESULTS = os.path.expanduser("~/mg_aki/results")
os.makedirs(RESULTS, exist_ok=True)

_FULL = os.path.expanduser("~/mg_aki/eicu-crd-2.0")
_DEMO = os.path.expanduser("~/mg_aki/eicu-collaborative-research-database-demo-2.0.1")
EICU_ROOT = _FULL if os.path.isdir(_FULL) else _DEMO

MIMIC_ROOT = os.path.expanduser("~/mg_aki/mimic-iv-3.1")
MIMIC_HOSP = os.path.join(MIMIC_ROOT, "hosp")
MIMIC_ICU = os.path.join(MIMIC_ROOT, "icu")

# ═══════════════════════════════════════════════════════════════════
# SHARED CONSTANTS
# ═══════════════════════════════════════════════════════════════════
MIN_AGE = 18
CR_MIN, CR_MAX = 0.1, 25.0
BASELINE_CR_MAX = 4.0

# Three Cr_post windows: (min_h, max_h) after IV Mg
CR_POST_WINDOWS = {
    "6_24h": (6, 24),  # primary
    "6_48h": (6, 48),  # sensitivity: wider
    "0_24h": (0, 24),  # sensitivity: no floor
}

# eICU: cardiac surgery identification
CARDIAC_DX_PATTERNS = [
    "cabg",
    "valve",
    "cardiac surgery",
    "open heart",
    "coronary artery bypass",
    "aortic valve",
    "mitral valve",
    "cardiothoracic",
    "aortic surgery",
    "tricuspid",
    "pulmonic valve",
]
CARDIAC_UNIT_TYPES = {"CSICU", "CTICU", "CCU-CTICU"}

# eICU: IV Mg drug patterns
MG_SUPP_PATTERNS = [
    "magnesium",
    "mag sulfate",
    "mgso4",
    "mag oxide",
    "mag gluconate",
    "mag hydroxide",
    "mag chloride",
]

# MIMIC: item IDs
MG_SUPP_ITEMS_MIMIC = [222011, 227523]
LAB_CR_MIMIC = [50912, 52546]
LAB_MG_MIMIC = [50960]
LAB_K_MIMIC = [50971]
LAB_CA_MIMIC = [50893]
LAB_LAC_MIMIC = [50813]
VITAL_HR_MIMIC = [220045]
VASO_ITEMS_MIMIC = [221906, 221289, 222315, 221749, 221662, 221653, 221986]
K_SUPP_ITEMS_MIMIC = [225166, 225168, 222139, 227521, 227522]
BLOOD_ITEMS_MIMIC = [225168, 220970, 225170, 225171, 226368, 226369, 226370, 226371]

# MIMIC: surgery ICD codes
CABG_ICD9 = ["3610", "3611", "3612", "3613", "3614", "3615", "3616", "3617", "3619"]
VALVE_ICD9 = [
    "3521",
    "3522",
    "3523",
    "3524",
    "3525",
    "3526",
    "3527",
    "3528",
    "3511",
    "3512",
    "3513",
    "3514",
]
CABG_ICD10 = ["0210", "0211", "0212", "0213"]
VALVE_ICD10 = ["02RF", "02RG", "02RH", "02RJ", "02QF", "02QG", "02QH", "02QJ"]
CVICU = "Cardiac Vascular Intensive Care Unit (CVICU)"

# eICU comorbidity patterns (pastHistory)
EICU_COMORB = {
    "heart_failure": ["heart failure", "chf", "cardiomyopathy"],
    "hypertension": ["hypertension"],
    "diabetes": ["diabetes"],
    "ckd": ["chronic kidney", "chronic renal", "ckd"],
    "copd": ["copd", "chronic obstructive", "emphysema"],
    "pvd": ["peripheral vascular", "pvd", "claudication"],
    "stroke": ["stroke", "cva", "cerebrovascular"],
    "liver_disease": ["cirrhosis", "hepatitis", "liver disease", "liver failure"],
}

# MIMIC comorbidity ICD codes
MIMIC_COMORB_ICD = {
    "heart_failure": {9: ["4280", "4281", "4289", "428"], 10: ["I50"]},
    "hypertension": {
        9: ["401", "402", "403", "404", "405"],
        10: ["I10", "I11", "I12", "I13", "I15"],
    },
    "diabetes": {9: ["250"], 10: ["E08", "E09", "E10", "E11", "E12", "E13"]},
    "ckd": {9: ["585", "586"], 10: ["N18", "N19"]},
    "copd": {
        9: ["490", "491", "492", "493", "494", "496"],
        10: ["J40", "J41", "J42", "J43", "J44", "J45", "J47"],
    },
    "pvd": {9: ["4431", "4439", "4471"], 10: ["I73"]},
    "stroke": {
        9: ["430", "431", "432", "433", "434", "435", "436"],
        10: ["I60", "I61", "I62", "I63", "I64", "I65", "I66", "G45"],
    },
    "liver_disease": {
        9: ["571"],
        10: ["K70", "K71", "K72", "K73", "K74", "K75", "K76"],
    },
}

# Drug patterns for both databases
DRUG_CLASSES = {
    "loop_diuretics": ["furosemide", "bumetanide", "torsemide", "lasix"],
    "nsaids": ["ibuprofen", "ketorolac", "naproxen", "diclofenac", "celecoxib"],
    "acei_arb": [
        "lisinopril",
        "enalapril",
        "ramipril",
        "captopril",
        "losartan",
        "valsartan",
        "irbesartan",
        "candesartan",
        "olmesartan",
        "telmisartan",
    ],
    "ppi": [
        "omeprazole",
        "pantoprazole",
        "lansoprazole",
        "esomeprazole",
        "rabeprazole",
    ],
    "beta_blockers": [
        "metoprolol",
        "atenolol",
        "propranolol",
        "carvedilol",
        "labetalol",
        "bisoprolol",
        "esmolol",
        "nadolol",
    ],
    "steroids": [
        "methylprednisolone",
        "hydrocortisone",
        "dexamethasone",
        "prednisone",
        "prednisolone",
        "solumedrol",
    ],
    "antiarrhythmics": [
        "amiodarone",
        "lidocaine",
        "procainamide",
        "flecainide",
        "sotalol",
        "dronedarone",
        "digoxin",
    ],
}

ESKD_PATTERNS = [
    "dialysis",
    "esrd",
    "end stage renal",
    "end-stage renal",
    "renal transplant",
    "kidney transplant",
]

VASO_PATTERNS = [
    "norepinephrine",
    "vasopressin",
    "epinephrine",
    "phenylephrine",
    "dopamine",
    "dobutamine",
    "milrinone",
]

TRANSFUSION_PATTERNS = [
    "transfusion",
    "blood product",
    "packed red",
    "prbc",
    "red blood cell",
    "ffp",
    "fresh frozen",
    "platelet",
    "cryoprecipitate",
]


# ═══════════════════════════════════════════════════════════════════
# HELPERS
# ═══════════════════════════════════════════════════════════════════
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


def pct(n, total):
    return f"{n:,} ({100*n/max(total,1):.1f}%)"


def desc(s, name, ind="    "):
    v = s.dropna()
    if len(v) == 0:
        print(f"{ind}{name}: all missing")
        return
    print(
        f"{ind}{name}: n={len(v)}, mean={v.mean():.2f}, "
        f"median={v.median():.2f}, IQR=[{v.quantile(.25):.2f}–"
        f"{v.quantile(.75):.2f}], range=[{v.min():.2f}–{v.max():.2f}]"
    )


def matches_icd(dx_df, hadm_ids, code_map):
    """Return set of hadm_ids matching ICD code prefixes."""
    sub = dx_df[dx_df.hadm_id.isin(hadm_ids)]
    hits = set()
    for ver, prefixes in code_map.items():
        v = sub[sub.icd_version == ver]
        for p in prefixes:
            hits |= set(v[v.icd_code.str.startswith(p)].hadm_id)
    return hits


def classify_surgery_eicu(dx_str):
    s = str(dx_str).lower()
    has_cabg = any(p in s for p in ["cabg", "coronary artery bypass"])
    has_valve = any(
        p in s
        for p in [
            "valve",
            "aortic valve",
            "mitral valve",
            "tricuspid",
            "pulmonic valve",
        ]
    )
    if has_cabg and has_valve:
        return "combined"
    if has_cabg:
        return "cabg"
    if has_valve:
        return "valve"
    return "other_cardiac"


def build_cohort_common(treated, control, db_tag):
    """Shared final steps: surgery dummies, print summary, save."""
    for df in [treated, control]:
        df["surg_cabg"] = (df.surgery_type == "cabg").astype(int)
        df["surg_valve"] = (df.surgery_type == "valve").astype(int)
        df["surg_combined"] = (df.surgery_type == "combined").astype(int)
        if "first_lactate" in df.columns:
            df["lactate_missing"] = df.first_lactate.isna().astype(int)

    treated["treated"] = 1
    control["treated"] = 0

    # KDIGO staging for reference (from primary 6-24h window)
    if "cr_post_6_24h" in treated.columns:
        treated["aki_kdigo1_ref"] = (
            (treated.cr_post_6_24h / treated.cr_pre >= 1.5)
            | ((treated.cr_post_6_24h - treated.cr_pre) >= 0.3)
        ).astype(int)
        treated.loc[treated.cr_post_6_24h.isna(), "aki_kdigo1_ref"] = np.nan

    tag = db_tag.lower()
    treated.to_csv(os.path.join(RESULTS, f"20_did_treated_{tag}.csv"), index=False)
    control.to_csv(os.path.join(RESULTS, f"20_did_control_{tag}.csv"), index=False)
    print(
        f"\n  Saved: 20_did_treated_{tag}.csv ({len(treated)} pts, {treated.shape[1]} cols)"
    )
    print(
        f"  Saved: 20_did_control_{tag}.csv ({len(control)} pts, {control.shape[1]} cols)"
    )


def print_consort(consort, treated, db_tag):
    """CONSORT summary and subgroup counts."""
    SEP = "=" * 70
    print(f"\n{SEP}\n{db_tag} CONSORT SUMMARY\n{SEP}")
    for label, key in [
        ("Total ICU admissions", "total_icu"),
        ("Cardiac surgery, adult, 1st stay", "cardiac_adult_first"),
        ("Post-ESKD exclusion", "post_eskd"),
        ("Received any postop IV Mg", "treated_any_ivmg"),
        ("  Has ICU Cr before IV Mg", "treated_has_cr_pre"),
        ("  No ICU Cr before IV Mg (LOST)", "treated_no_cr_pre"),
        ("  Excl: Cr_pre ≥ 4.0", "excl_cr_high"),
        ("  Excl: prevalent AKI", "excl_prevalent_aki"),
        ("  TREATED FINAL", "treated_final"),
        ("Never received IV Mg", "control_no_ivmg"),
        ("  Has ≥2 postop Cr", "control_has_2cr"),
        ("  CONTROL FINAL", "control_final"),
    ]:
        print(f"  {label:<42s} {consort.get(key, 0):>8,}")

    n_trt = len(treated)
    print(f"\n  Cr_post availability (treated final, n={n_trt}):")
    for wname in CR_POST_WINDOWS:
        col = f"cr_post_{wname}"
        if col in treated.columns:
            n = treated[col].notna().sum()
            tag = " ◀ PRIMARY" if wname == "6_24h" else ""
            print(f"    {wname:>6s}: {pct(n, n_trt)}{tag}")

    for wname in CR_POST_WINDOWS:
        col = f"cr_post_{wname}"
        if col in treated.columns:
            n_ready = treated[col].notna().sum()
            print(f"\n  DiD-ready (Cr_pre + {wname}): {pct(n_ready, n_trt)}")
            if n_ready > 0:
                desc(
                    treated.loc[treated[col].notna(), f"delta_cr_{wname}"],
                    f"ΔCr ({wname})",
                )

    print(f"\n  Surgery type (treated):")
    for st in ["cabg", "valve", "combined", "other_cardiac"]:
        print(f"    {st:<15s}: {pct((treated.surgery_type==st).sum(), n_trt)}")

    if "first_mg_value" in treated.columns:
        print(f"\n  Serum Mg strata (treated):")
        for lo, hi, lbl in [
            (0, 1.8, "<1.8"),
            (1.8, 2.0, "1.8-2.0"),
            (2.0, 2.3, "2.0-2.3"),
            (2.3, 99, ">2.3"),
        ]:
            n = treated.first_mg_value.between(lo, hi, inclusive="left").sum()
            print(f"    {lbl:<10s}: {pct(n, n_trt)}")


# ═══════════════════════════════════════════════════════════════════
#  eICU-CRD
# ═══════════════════════════════════════════════════════════════════
def run_eicu():
    SEP = "=" * 70
    print(f"\n{SEP}\neICU-CRD: DiD Cohort Construction\n  Data: {EICU_ROOT}\n{SEP}")
    consort = {}

    # ── load tables ──────────────────────────────────────────────
    def ld(name, **kw):
        for ext in [".csv.gz", ".csv"]:
            p = os.path.join(EICU_ROOT, name + ext)
            if os.path.exists(p):
                df = pd.read_csv(p, low_memory=False, **kw)
                df.columns = df.columns.str.lower()
                print(f"  Loaded {name}: {len(df):,}")
                return df
        print(f"  WARNING: {name} not found")
        return pd.DataFrame()

    patient = ld("patient")
    lab = ld("lab")
    med = ld(
        "medication",
        usecols=[
            "patientunitstayid",
            "drugstartoffset",
            "drugstopoffset",
            "drugname",
            "drugordercancelled",
            "routeadmin",
            "dosage",
        ],
    )
    inf = ld("infusionDrug")
    diag = ld("diagnosis")
    pasthx = ld("pastHistory")
    tx_tbl = ld("treatment")
    vital = ld(
        "vitalPeriodic", usecols=["patientunitstayid", "observationoffset", "heartrate"]
    )

    consort["total_icu"] = len(patient)
    print(f"\n  Total ICU admissions: {len(patient):,}")

    # ── cardiac surgery, adult, first stay ───────────────────────
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
    pids = set(cardiac.patientunitstayid)
    consort["cardiac_adult_first"] = len(cardiac)
    print(f"  Cardiac surgery, adult, 1st stay: {len(cardiac):,}")

    # ── ESKD ─────────────────────────────────────────────────────
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
    print(f"  ESKD excluded: {len(eskd):,} → remaining: {len(cardiac):,}")

    # surgery type
    cardiac["surgery_type"] = cardiac.apacheadmissiondx.apply(classify_surgery_eicu)

    # ── IV Mg: any postop ────────────────────────────────────────
    print("\n── IV Magnesium identification ──")
    med_ok = (
        med[med.patientunitstayid.isin(pids) & (med.drugordercancelled != "Yes")]
        if len(med) > 0
        else pd.DataFrame()
    )

    frames = []
    if len(med_ok) > 0:
        mg_m = med_ok[
            matches_any(med_ok.drugname, MG_SUPP_PATTERNS)
            & (med_ok.drugstartoffset >= 0)
        ].copy()
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
    print(f"  IV Mg (any postop): {pct(len(treated_pids), len(pids))}")
    print(f"  No IV Mg (controls): {len(control_pids):,}")
    desc(first_mg.mg_offset_h, "IV Mg timing (h from ICU)")
    for c in [6, 12, 24, 48]:
        print(
            f"    Within {c:2d}h: {pct((first_mg.mg_offset_h<=c).sum(), len(first_mg))}"
        )

    # ── Creatinine ───────────────────────────────────────────────
    print("\n── Creatinine ──")
    cr = lab[
        lab.patientunitstayid.isin(pids)
        & lab.labname.str.lower().str.contains("creatinine", na=False)
        & lab.labresult.between(CR_MIN, CR_MAX)
    ].copy()
    print(
        f"  Cr measurements: {len(cr):,} across {cr.patientunitstayid.nunique():,} pts"
    )

    # Cr_pre: latest ICU Cr before IV Mg
    print("\n── Cr_pre (latest ICU Cr before IV Mg) ──")
    cr_t = cr[cr.patientunitstayid.isin(treated_pids)].merge(
        first_mg[["patientunitstayid", "mg_offset_min"]], on="patientunitstayid"
    )
    cr_pre_cand = cr_t[
        (cr_t.labresultoffset >= 0) & (cr_t.labresultoffset < cr_t.mg_offset_min)
    ]
    cr_pre = (
        cr_pre_cand.sort_values("labresultoffset", ascending=False)
        .groupby("patientunitstayid")
        .first()
        .reset_index()
        .rename(columns={"labresult": "cr_pre", "labresultoffset": "cr_pre_offset_min"})
    )
    cr_pre["cr_pre_offset_h"] = cr_pre.cr_pre_offset_min / 60.0
    cr_pre["gap_to_ivmg_h"] = (cr_pre.mg_offset_min - cr_pre.cr_pre_offset_min) / 60.0

    n_has = len(cr_pre)
    n_no = len(treated_pids) - n_has
    consort["treated_has_cr_pre"] = n_has
    consort["treated_no_cr_pre"] = n_no
    print(f"  Has ICU Cr before IV Mg: {pct(n_has, len(treated_pids))}")
    print(f"  No ICU Cr before IV Mg:  {pct(n_no, len(treated_pids))}")
    desc(cr_pre.cr_pre, "Cr_pre (mg/dL)")
    desc(cr_pre.gap_to_ivmg_h, "Gap Cr_pre → IV Mg (h)")

    # Comparison: ICU Cr_pre vs hosp Cr
    print("\n── ICU Cr_pre vs hospitalization Cr ──")
    hosp_off = cardiac.set_index("patientunitstayid")["hospitaladmitoffset"].to_dict()
    cr_t["hosp_off"] = cr_t.patientunitstayid.map(hosp_off)
    cr_h = cr_t[
        (cr_t.labresultoffset >= cr_t.hosp_off - 360)
        & (cr_t.labresultoffset <= cr_t.hosp_off + 360)
        & (cr_t.labresultoffset < 0)
    ].copy()
    cr_h["dist"] = (cr_h.labresultoffset - cr_h.hosp_off).abs()
    cr_hosp = (
        cr_h.sort_values("dist").groupby("patientunitstayid").first().reset_index()
    )
    both = cr_pre[["patientunitstayid", "cr_pre"]].merge(
        cr_hosp[["patientunitstayid", "labresult"]].rename(
            columns={"labresult": "cr_hosp"}
        ),
        on="patientunitstayid",
        how="left",
    )
    nb = both.cr_hosp.notna().sum()
    if nb > 0:
        bv = both.dropna(subset=["cr_hosp"])
        d = bv.cr_pre - bv.cr_hosp
        print(f"  Both available: {nb}")
        desc(bv.cr_hosp, "Hosp Cr")
        desc(bv.cr_pre, "ICU Cr_pre")
        desc(d, "Δ (ICU − hosp)")
        print(
            f"    r = {bv.cr_pre.corr(bv.cr_hosp):.3f}, "
            f"|Δ|>0.3: {pct((d.abs()>0.3).sum(), nb)}"
        )
    n_recov = cr_hosp.patientunitstayid.isin(
        treated_pids - set(cr_pre.patientunitstayid)
    ).sum()
    print(f"  Lost patients recoverable via hosp Cr: {n_recov}")

    # Cr_post windows
    print("\n── Cr_post windows ──")
    cr_post_all = cr_t[cr_t.labresultoffset > cr_t.mg_offset_min].copy()
    cr_post_all["post_h"] = (
        cr_post_all.labresultoffset - cr_post_all.mg_offset_min
    ) / 60.0
    for wname, (lo, hi) in CR_POST_WINDOWS.items():
        cand = cr_post_all[(cr_post_all.post_h >= lo) & (cr_post_all.post_h <= hi)]
        n = cand.groupby("patientunitstayid").ngroups
        tag = " ◀ PRIMARY" if wname == "6_24h" else ""
        print(f"  {wname:>6s}: {pct(n, n_has)}{tag}")

    # Exclusions
    cr_pre_ok = cr_pre[cr_pre.cr_pre < BASELINE_CR_MAX]
    consort["excl_cr_high"] = n_has - len(cr_pre_ok)
    print(f"\n  Excl Cr_pre ≥ {BASELINE_CR_MAX}: {n_has - len(cr_pre_ok)}")

    # Prevalent AKI
    cr_earliest = (
        cr_pre_cand.sort_values("labresultoffset")
        .groupby("patientunitstayid")
        .first()
        .reset_index()
        .rename(columns={"labresult": "cr_earliest"})
    )
    prev = cr_pre_ok.merge(
        cr_earliest[["patientunitstayid", "cr_earliest"]],
        on="patientunitstayid",
        how="left",
    )
    prev["prevalent"] = (prev.cr_pre / prev.cr_earliest.clip(lower=0.1) >= 1.5).astype(
        int
    )
    n_prev = prev.prevalent.sum()
    consort["excl_prevalent_aki"] = n_prev
    print(f"  Prevalent AKI: {pct(n_prev, len(prev))}")

    # ── TREATED ──────────────────────────────────────────────────
    print(f"\n{'─'*50}\nBuilding treated cohort...")
    keep_pids = set(cr_pre_ok.patientunitstayid) - set(
        prev[prev.prevalent == 1].patientunitstayid
    )
    treated = cardiac[cardiac.patientunitstayid.isin(keep_pids)].copy()
    treated = treated.merge(
        first_mg[["patientunitstayid", "mg_offset_min", "mg_offset_h"]],
        on="patientunitstayid",
    )
    treated = treated.merge(
        cr_pre[
            [
                "patientunitstayid",
                "cr_pre",
                "cr_pre_offset_min",
                "cr_pre_offset_h",
                "gap_to_ivmg_h",
            ]
        ],
        on="patientunitstayid",
    )
    treated["egfr"] = compute_egfr(treated.cr_pre, treated.age, treated.is_female)

    # Cr_post at each window
    for wname, (lo, hi) in CR_POST_WINDOWS.items():
        cand = cr_post_all[
            cr_post_all.patientunitstayid.isin(set(treated.patientunitstayid))
            & (cr_post_all.post_h >= lo)
            & (cr_post_all.post_h <= hi)
        ]
        fp = (
            cand.sort_values("labresultoffset")
            .groupby("patientunitstayid")
            .first()
            .reset_index()[["patientunitstayid", "labresult", "labresultoffset"]]
            .rename(
                columns={
                    "labresult": f"cr_post_{wname}",
                    "labresultoffset": f"cr_post_offset_{wname}",
                }
            )
        )
        treated = treated.merge(fp, on="patientunitstayid", how="left")
        treated[f"delta_cr_{wname}"] = treated[f"cr_post_{wname}"] - treated.cr_pre

    consort["treated_final"] = len(treated)

    # Comorbidities
    for como, pats in EICU_COMORB.items():
        cp = set()
        if len(pasthx) > 0 and "pasthistorypath" in pasthx.columns:
            cp = set(
                pasthx[
                    pasthx.patientunitstayid.isin(set(treated.patientunitstayid))
                    & matches_any(pasthx.pasthistorypath, pats)
                ].patientunitstayid
            )
        treated[como] = treated.patientunitstayid.isin(cp).astype(int)

    # Drugs before t0 (patient-specific)
    print("  Covariates (pre-t0)...")
    mg_off_d = dict(zip(treated.patientunitstayid, treated.mg_offset_min))
    tpids = set(treated.patientunitstayid)
    for dc, pats in DRUG_CLASSES.items():
        flagged = set()
        if len(med_ok) > 0:
            for pid in tpids:
                mo = mg_off_d.get(pid, 0)
                pt = med_ok[
                    (med_ok.patientunitstayid == pid) & (med_ok.drugstartoffset <= mo)
                ]
                if len(pt) > 0 and matches_any(pt.drugname, pats).any():
                    flagged.add(pid)
        treated[dc] = treated.patientunitstayid.isin(flagged).astype(int)

    # Labs before t0
    for col, lab_pats in [
        ("first_mg_value", ["magnesium"]),
        ("first_potassium", ["potassium"]),
        ("first_calcium", ["calcium"]),
        ("first_lactate", ["lactate"]),
    ]:
        lsub = lab[
            lab.patientunitstayid.isin(tpids) & matches_any(lab.labname, lab_pats)
        ]
        vals = {}
        for pid in tpids:
            mo = mg_off_d.get(pid, 0)
            pt = lsub[
                (lsub.patientunitstayid == pid)
                & (lsub.labresultoffset >= 0)
                & (lsub.labresultoffset < mo)
            ]
            if len(pt) > 0:
                vals[pid] = pt.sort_values("labresultoffset").iloc[0].labresult
        treated[col] = treated.patientunitstayid.map(vals)
        print(f"    {col}: {pct(treated[col].notna().sum(), len(treated))}")

    # ◆ Heart rate from vitalPeriodic — before t0
    if len(vital) > 0 and "heartrate" in vital.columns:
        vt = vital[
            vital.patientunitstayid.isin(tpids)
            & vital.heartrate.notna()
            & vital.heartrate.between(20, 250)
        ]
        hr_vals = {}
        for pid in tpids:
            mo = mg_off_d.get(pid, 0)
            pt = vt[
                (vt.patientunitstayid == pid)
                & (vt.observationoffset >= 0)
                & (vt.observationoffset < mo)
            ]
            if len(pt) > 0:
                hr_vals[pid] = pt.sort_values("observationoffset").iloc[0].heartrate
        treated["first_heartrate"] = treated.patientunitstayid.map(hr_vals)
        print(
            f"    first_heartrate: {pct(treated.first_heartrate.notna().sum(), len(treated))}"
        )
    else:
        treated["first_heartrate"] = np.nan
        print(f"    first_heartrate: vitalPeriodic not available")

    # Vasopressor & transfusion before t0
    if len(med_ok) > 0:
        vp = set()
        for pid in tpids:
            mo = mg_off_d.get(pid, 0)
            pt = med_ok[
                (med_ok.patientunitstayid == pid)
                & (med_ok.drugstartoffset >= 0)
                & (med_ok.drugstartoffset <= mo)
            ]
            if len(pt) > 0 and matches_any(pt.drugname, VASO_PATTERNS).any():
                vp.add(pid)
        treated["vasopressor_pre_t0"] = treated.patientunitstayid.isin(vp).astype(int)
        print(
            f"    vasopressor_pre_t0: {pct(treated.vasopressor_pre_t0.sum(), len(treated))}"
        )

    if len(tx_tbl) > 0:
        tp = set()
        tx_m = tx_tbl[
            tx_tbl.patientunitstayid.isin(tpids)
            & matches_any(tx_tbl.treatmentstring, TRANSFUSION_PATTERNS)
        ]
        for pid in tpids:
            mo = mg_off_d.get(pid, 0)
            pt = tx_m[
                (tx_m.patientunitstayid == pid)
                & (tx_m.treatmentoffset >= 0)
                & (tx_m.treatmentoffset <= mo)
            ]
            if len(pt) > 0:
                tp.add(pid)
        treated["transfusion_pre_t0"] = treated.patientunitstayid.isin(tp).astype(int)
        print(
            f"    transfusion_pre_t0: {pct(treated.transfusion_pre_t0.sum(), len(treated))}"
        )

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

    # Mortality
    mort_map_icu = (
        cardiac.set_index("patientunitstayid")["unitdischargestatus"]
        .str.lower()
        .eq("expired")
    ).to_dict()
    mort_map_hosp = (
        cardiac.set_index("patientunitstayid")["hospitaldischargestatus"]
        .str.lower()
        .eq("expired")
    ).to_dict()
    treated["icu_mortality"] = (
        treated.patientunitstayid.map(mort_map_icu).fillna(False).astype(int)
    )
    treated["hosp_mortality"] = (
        treated.patientunitstayid.map(mort_map_hosp).fillna(False).astype(int)
    )

    # ── CONTROL ──────────────────────────────────────────────────
    print(f"\n{'─'*50}\nBuilding control cohort...")
    cr_ctrl = cr[(cr.patientunitstayid.isin(control_pids)) & (cr.labresultoffset >= 0)]
    c2 = set(
        cr_ctrl.groupby("patientunitstayid").size().pipe(lambda s: s[s >= 2]).index
    )
    control = cardiac[cardiac.patientunitstayid.isin(c2)].copy()
    consort["control_has_2cr"] = len(control)
    print(f"  Controls with ≥2 postop Cr: {pct(len(c2), len(control_pids))}")

    c1 = (
        cr_ctrl[cr_ctrl.patientunitstayid.isin(c2)]
        .sort_values("labresultoffset")
        .groupby("patientunitstayid")
        .first()
        .reset_index()
    )
    excl = set(c1[c1.labresult >= BASELINE_CR_MAX].patientunitstayid)
    control = control[~control.patientunitstayid.isin(excl)].copy()
    print(f"  Excl first Cr ≥ {BASELINE_CR_MAX}: {len(excl)}")

    # eGFR
    c1_ok = c1[c1.patientunitstayid.isin(set(control.patientunitstayid))]
    control = control.merge(
        c1_ok[["patientunitstayid", "labresult", "labresultoffset"]].rename(
            columns={
                "labresult": "first_postop_cr",
                "labresultoffset": "first_cr_offset_min",
            }
        ),
        on="patientunitstayid",
        how="left",
    )
    control["egfr"] = compute_egfr(
        control.first_postop_cr, control.age, control.is_female
    )

    # Comorbidities
    for como, pats in EICU_COMORB.items():
        cp = set()
        if len(pasthx) > 0 and "pasthistorypath" in pasthx.columns:
            cp = set(
                pasthx[
                    pasthx.patientunitstayid.isin(set(control.patientunitstayid))
                    & matches_any(pasthx.pasthistorypath, pats)
                ].patientunitstayid
            )
        control[como] = control.patientunitstayid.isin(cp).astype(int)

    # Drugs: first 6h (common anchor for controls)
    cpids = set(control.patientunitstayid)
    for dc, pats in DRUG_CLASSES.items():
        f = set()
        if len(med_ok) > 0:
            early = med_ok[
                (med_ok.patientunitstayid.isin(cpids))
                & (med_ok.drugstartoffset >= 0)
                & (med_ok.drugstartoffset <= 360)
            ]
            if len(early) > 0:
                f = set(early[matches_any(early.drugname, pats)].patientunitstayid)
        control[dc] = control.patientunitstayid.isin(f).astype(int)

    # Labs: first 6h
    for col, lab_pats in [
        ("first_mg_value", ["magnesium"]),
        ("first_potassium", ["potassium"]),
        ("first_calcium", ["calcium"]),
        ("first_lactate", ["lactate"]),
    ]:
        lsub = lab[
            lab.patientunitstayid.isin(cpids)
            & matches_any(lab.labname, lab_pats)
            & (lab.labresultoffset >= 0)
            & (lab.labresultoffset <= 360)
        ]
        fv = (
            lsub.sort_values("labresultoffset")
            .groupby("patientunitstayid")
            .first()
            .reset_index()[["patientunitstayid", "labresult"]]
            .rename(columns={"labresult": col})
        )
        control = control.merge(fv, on="patientunitstayid", how="left")

    # Heart rate: first 6h
    if len(vital) > 0 and "heartrate" in vital.columns:
        vt_c = vital[
            vital.patientunitstayid.isin(cpids)
            & vital.heartrate.notna()
            & vital.heartrate.between(20, 250)
            & (vital.observationoffset >= 0)
            & (vital.observationoffset <= 360)
        ]
        hr_c = (
            vt_c.sort_values("observationoffset")
            .groupby("patientunitstayid")
            .first()
            .reset_index()[["patientunitstayid", "heartrate"]]
            .rename(columns={"heartrate": "first_heartrate"})
        )
        control = control.merge(hr_c, on="patientunitstayid", how="left")
    else:
        control["first_heartrate"] = np.nan

    # Vasopressor/transfusion: first 6h
    if len(med_ok) > 0:
        ve = med_ok[
            (med_ok.patientunitstayid.isin(cpids))
            & (med_ok.drugstartoffset >= 0)
            & (med_ok.drugstartoffset <= 360)
        ]
        vp = (
            set(ve[matches_any(ve.drugname, VASO_PATTERNS)].patientunitstayid)
            if len(ve) > 0
            else set()
        )
        control["vasopressor_pre_t0"] = control.patientunitstayid.isin(vp).astype(int)

    if len(tx_tbl) > 0:
        te = tx_tbl[
            tx_tbl.patientunitstayid.isin(cpids)
            & matches_any(tx_tbl.treatmentstring, TRANSFUSION_PATTERNS)
            & (tx_tbl.treatmentoffset >= 0)
            & (tx_tbl.treatmentoffset <= 360)
        ]
        control["transfusion_pre_t0"] = control.patientunitstayid.isin(
            set(te.patientunitstayid)
        ).astype(int)

    # BMI
    if "admissionweight" in cardiac.columns and "admissionheight" in cardiac.columns:
        control["bmi"] = control.patientunitstayid.apply(
            lambda x: (
                wt.get(x, np.nan) / ((ht.get(x, np.nan) / 100) ** 2)
                if pd.notna(wt.get(x)) and pd.notna(ht.get(x)) and ht.get(x, 0) > 0
                else np.nan
            )
        )
        control.loc[~control.bmi.between(10, 80), "bmi"] = np.nan

    # Mortality
    control["icu_mortality"] = (
        control.patientunitstayid.map(mort_map_icu).fillna(False).astype(int)
    )
    control["hosp_mortality"] = (
        control.patientunitstayid.map(mort_map_hosp).fillna(False).astype(int)
    )

    consort["control_final"] = len(control)

    # ── All Cr for temporal matching ─────────────────────────────
    all_pids = set(treated.patientunitstayid) | set(control.patientunitstayid)
    cr_exp = cr[cr.patientunitstayid.isin(all_pids) & (cr.labresultoffset >= 0)][
        ["patientunitstayid", "labresult", "labresultoffset"]
    ].copy()
    cr_exp["cr_offset_h"] = cr_exp.labresultoffset / 60.0
    cr_exp.to_csv(os.path.join(RESULTS, "20_did_cr_all_eicu.csv"), index=False)
    print(
        f"\n  Exported {len(cr_exp):,} Cr for {cr_exp.patientunitstayid.nunique()} pts"
    )

    # ── SMD check ────────────────────────────────────────────────
    print(f"\n{'─'*50}\nCOVARIATE TIMING ASSESSMENT")
    for v in ["vasopressor_pre_t0", "transfusion_pre_t0"]:
        if v in treated.columns and v in control.columns:
            tr = treated[v].mean()
            cr_ = control[v].mean()
            smd = abs(tr - cr_) / np.sqrt((tr * (1 - tr) + cr_ * (1 - cr_)) / 2 + 1e-10)
            warn = " ⚠ consider dropping" if smd > 0.25 else ""
            print(f"  {v}: trt={100*tr:.1f}% ctrl={100*cr_:.1f}% SMD={smd:.3f}{warn}")

    # ── save & report ────────────────────────────────────────────
    build_cohort_common(treated, control, "eicu")
    print_consort(consort, treated, "eICU")
    return consort


# ═══════════════════════════════════════════════════════════════════
#  MIMIC-IV
# ═══════════════════════════════════════════════════════════════════
def run_mimic():
    SEP = "=" * 70
    print(f"\n{SEP}\nMIMIC-IV: DiD Cohort Construction\n  Data: {MIMIC_ROOT}\n{SEP}")
    consort = {}

    # ── load ─────────────────────────────────────────────────────
    patients = pd.read_csv(gz(f"{MIMIC_HOSP}/patients.csv.gz"))
    admissions = pd.read_csv(gz(f"{MIMIC_HOSP}/admissions.csv.gz"))
    icustays = pd.read_csv(gz(f"{MIMIC_ICU}/icustays.csv.gz"))
    dx = pd.read_csv(gz(f"{MIMIC_HOSP}/diagnoses_icd.csv.gz"))
    px = pd.read_csv(gz(f"{MIMIC_HOSP}/procedures_icd.csv.gz"))
    presc = pd.read_csv(
        gz(f"{MIMIC_HOSP}/prescriptions.csv.gz"),
        usecols=["subject_id", "hadm_id", "drug", "starttime", "stoptime", "route"],
    )

    print("  Loading labevents (chunked)...")
    needed = set(
        LAB_CR_MIMIC + LAB_MG_MIMIC + LAB_K_MIMIC + LAB_CA_MIMIC + LAB_LAC_MIMIC
    )
    lc = []
    for ch in pd.read_csv(
        gz(f"{MIMIC_HOSP}/labevents.csv.gz"),
        usecols=["subject_id", "hadm_id", "itemid", "charttime", "valuenum"],
        dtype={"subject_id": int, "hadm_id": "Int64", "itemid": int},
        chunksize=5_000_000,
    ):
        lc.append(ch[ch.itemid.isin(needed)])
    labs = pd.concat(lc, ignore_index=True)
    print(f"  labevents (filtered): {len(labs):,}")

    print("  Loading inputevents...")
    _ie_want = {
        "subject_id",
        "stay_id",
        "hadm_id",
        "itemid",
        "starttime",
        "amount",
        "amountuom",
        "statusdescription",
    }
    ie = pd.read_csv(
        gz(f"{MIMIC_ICU}/inputevents.csv.gz"),
        usecols=lambda c: c in _ie_want,
        low_memory=False,
    )
    print(f"  inputevents: {len(ie):,}")

    print("  Loading chartevents (HR only)...")
    ce_chunks = []
    for ch in pd.read_csv(
        gz(f"{MIMIC_ICU}/chartevents.csv.gz"),
        usecols=["stay_id", "itemid", "charttime", "valuenum"],
        chunksize=10_000_000,
    ):
        ce_chunks.append(ch[ch.itemid.isin(VITAL_HR_MIMIC)])
    ce_hr = pd.concat(ce_chunks, ignore_index=True) if ce_chunks else pd.DataFrame()
    print(f"  chartevents HR: {len(ce_hr):,}")

    consort["total_icu"] = len(icustays)
    print(f"\n  Total ICU stays: {len(icustays):,}")

    # ── cardiac surgery identification ───────────────────────────
    px["icd_code"] = px["icd_code"].astype(str).str.strip()
    dx["icd_code"] = dx["icd_code"].astype(str).str.strip()

    def classify_mimic(hadm_id):
        sub = px[px.hadm_id == hadm_id]
        codes = set(sub.icd_code)
        has_cabg = any(codes & set(c) for c in [CABG_ICD9, CABG_ICD10])
        has_valve = any(codes & set(c) for c in [VALVE_ICD9, VALVE_ICD10])
        if has_cabg and has_valve:
            return "combined"
        if has_cabg:
            return "cabg"
        if has_valve:
            return "valve"
        return "other_cardiac"

    # CVICU patients or those with cardiac surgery ICD codes
    icu = icustays.merge(
        patients[["subject_id", "gender", "anchor_age"]], on="subject_id"
    )
    icu = icu.merge(
        admissions[["hadm_id", "admittime", "dischtime", "hospital_expire_flag"]],
        on="hadm_id",
    )
    icu["intime"] = pd.to_datetime(icu["intime"])
    icu["outtime"] = pd.to_datetime(icu["outtime"])

    # Cardiac: CVICU or ICD procedure codes
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
    cardiac["age"] = cardiac["anchor_age"]
    cardiac["is_female"] = (cardiac.gender == "F").astype(int)
    cardiac["surgery_type"] = cardiac.hadm_id.apply(classify_mimic)
    stays = set(cardiac.stay_id)
    hadms = set(cardiac.hadm_id.dropna().astype(int))
    consort["cardiac_adult_first"] = len(cardiac)
    print(f"  Cardiac surgery, adult, 1st stay: {len(cardiac):,}")

    # ── ESKD ─────────────────────────────────────────────────────
    eskd_hadms = set()
    for pat in ESKD_PATTERNS:
        eskd_hadms |= set(
            dx[
                dx.hadm_id.isin(hadms)
                & dx.icd_code.str.lower().str.contains(pat.replace(" ", ""), na=False)
            ].hadm_id
        )
    # Also check ICD codes for ESKD
    for ver, codes in [
        (9, ["5856", "V4511", "V560", "V561", "V562"]),
        (10, ["N186", "Z491", "Z492", "Z9911", "Z940"]),
    ]:
        v = dx[(dx.hadm_id.isin(hadms)) & (dx.icd_version == ver)]
        for c in codes:
            eskd_hadms |= set(v[v.icd_code.str.startswith(c)].hadm_id)
    eskd_stays = set(cardiac[cardiac.hadm_id.isin(eskd_hadms)].stay_id)
    cardiac = cardiac[~cardiac.stay_id.isin(eskd_stays)]
    stays = set(cardiac.stay_id)
    hadms = set(cardiac.hadm_id.dropna().astype(int))
    consort["post_eskd"] = len(cardiac)
    print(f"  ESKD excluded: {len(eskd_stays):,} → remaining: {len(cardiac):,}")

    # ── IV Mg from inputevents ───────────────────────────────────
    print("\n── IV Magnesium identification ──")
    ie_mg = ie[(ie.stay_id.isin(stays)) & (ie.itemid.isin(MG_SUPP_ITEMS_MIMIC))].copy()
    # Filter out cancelled/rewritten orders
    if "statusdescription" in ie_mg.columns:
        ie_mg = ie_mg[~ie_mg.statusdescription.str.contains("Rewritten", na=False)]
    ie_mg = ie_mg[ie_mg.amount.notna() & (ie_mg.amount > 0)]
    ie_mg["starttime"] = pd.to_datetime(ie_mg["starttime"])
    ie_mg = ie_mg.merge(cardiac[["stay_id", "intime"]], on="stay_id")
    ie_mg["mg_offset_h"] = (ie_mg.starttime - ie_mg.intime).dt.total_seconds() / 3600
    ie_mg = ie_mg[ie_mg.mg_offset_h >= 0]  # postop only

    first_mg = ie_mg.sort_values("mg_offset_h").groupby("stay_id").first().reset_index()
    first_mg["mg_offset_min"] = first_mg.mg_offset_h * 60

    treated_stays = set(first_mg.stay_id)
    control_stays = stays - treated_stays
    consort["treated_any_ivmg"] = len(treated_stays)
    consort["control_no_ivmg"] = len(control_stays)
    print(f"  IV Mg (any postop): {pct(len(treated_stays), len(stays))}")
    print(f"  No IV Mg: {len(control_stays):,}")
    desc(first_mg.mg_offset_h, "IV Mg timing (h from ICU)")
    for c in [6, 12, 24, 48]:
        print(
            f"    Within {c:2d}h: {pct((first_mg.mg_offset_h<=c).sum(), len(first_mg))}"
        )

    # ── Creatinine ───────────────────────────────────────────────
    print("\n── Creatinine ──")
    cr_labs = labs[
        labs.itemid.isin(LAB_CR_MIMIC)
        & labs.hadm_id.isin(hadms)
        & labs.valuenum.between(CR_MIN, CR_MAX)
    ].copy()
    cr_labs["charttime"] = pd.to_datetime(cr_labs["charttime"])
    cr_labs = cr_labs.merge(cardiac[["stay_id", "hadm_id", "intime"]], on="hadm_id")
    cr_labs["offset_h"] = (cr_labs.charttime - cr_labs.intime).dt.total_seconds() / 3600
    cr_labs["offset_min"] = cr_labs.offset_h * 60
    print(
        f"  Cr measurements: {len(cr_labs):,} across {cr_labs.stay_id.nunique():,} pts"
    )

    # Cr_pre: latest ICU Cr before IV Mg
    print("\n── Cr_pre (latest ICU Cr before IV Mg) ──")
    cr_t = cr_labs[cr_labs.stay_id.isin(treated_stays)].merge(
        first_mg[["stay_id", "mg_offset_h", "mg_offset_min"]], on="stay_id"
    )
    cr_pre_cand = cr_t[(cr_t.offset_h >= 0) & (cr_t.offset_min < cr_t.mg_offset_min)]
    cr_pre = (
        cr_pre_cand.sort_values("offset_min", ascending=False)
        .groupby("stay_id")
        .first()
        .reset_index()
        .rename(
            columns={
                "valuenum": "cr_pre",
                "offset_min": "cr_pre_offset_min",
                "offset_h": "cr_pre_offset_h",
            }
        )
    )
    cr_pre["gap_to_ivmg_h"] = (cr_pre.mg_offset_min - cr_pre.cr_pre_offset_min) / 60.0

    n_has = len(cr_pre)
    n_no = len(treated_stays) - n_has
    consort["treated_has_cr_pre"] = n_has
    consort["treated_no_cr_pre"] = n_no
    print(f"  Has ICU Cr before IV Mg: {pct(n_has, len(treated_stays))}")
    print(f"  No ICU Cr before IV Mg:  {pct(n_no, len(treated_stays))}")
    desc(cr_pre.cr_pre, "Cr_pre (mg/dL)")
    desc(cr_pre.gap_to_ivmg_h, "Gap Cr_pre → IV Mg (h)")

    # Hosp Cr comparison
    print("\n── ICU Cr_pre vs hospitalization Cr ──")
    adm_times = admissions[["hadm_id", "admittime"]].copy()
    adm_times["admittime"] = pd.to_datetime(adm_times["admittime"])
    cr_t2 = cr_t.merge(adm_times, on="hadm_id", how="left")
    cr_t2["h_from_admit"] = (
        cr_t2.charttime - cr_t2.admittime
    ).dt.total_seconds() / 3600
    cr_h = cr_t2[
        (cr_t2.h_from_admit >= -24) & (cr_t2.h_from_admit <= 24) & (cr_t2.offset_h < 0)
    ]
    if len(cr_h) > 0:
        cr_h["dist"] = cr_h.h_from_admit.abs()
        cr_hosp = cr_h.sort_values("dist").groupby("stay_id").first().reset_index()
        both = cr_pre[["stay_id", "cr_pre"]].merge(
            cr_hosp[["stay_id", "valuenum"]].rename(columns={"valuenum": "cr_hosp"}),
            on="stay_id",
            how="left",
        )
        nb = both.cr_hosp.notna().sum()
        if nb > 0:
            bv = both.dropna(subset=["cr_hosp"])
            d = bv.cr_pre - bv.cr_hosp
            print(f"  Both available: {nb}")
            desc(bv.cr_hosp, "Hosp Cr")
            desc(bv.cr_pre, "ICU Cr_pre")
            desc(d, "Δ (ICU − hosp)")
            print(
                f"    r = {bv.cr_pre.corr(bv.cr_hosp):.3f}, "
                f"|Δ|>0.3: {pct((d.abs()>0.3).sum(), nb)}"
            )
    else:
        print(f"  No pre-ICU Cr found for comparison")

    # Cr_post windows
    print("\n── Cr_post windows ──")
    cr_post_all = cr_t[cr_t.offset_min > cr_t.mg_offset_min].copy()
    cr_post_all["post_h"] = (cr_post_all.offset_min - cr_post_all.mg_offset_min) / 60.0
    for wname, (lo, hi) in CR_POST_WINDOWS.items():
        cand = cr_post_all[(cr_post_all.post_h >= lo) & (cr_post_all.post_h <= hi)]
        n = cand.groupby("stay_id").ngroups
        tag = " ◀ PRIMARY" if wname == "6_24h" else ""
        print(f"  {wname:>6s}: {pct(n, n_has)}{tag}")

    # Exclusions
    cr_pre_ok = cr_pre[cr_pre.cr_pre < BASELINE_CR_MAX]
    consort["excl_cr_high"] = n_has - len(cr_pre_ok)
    print(f"\n  Excl Cr_pre ≥ {BASELINE_CR_MAX}: {n_has - len(cr_pre_ok)}")

    # Prevalent AKI
    cr_earliest = (
        cr_pre_cand.sort_values("offset_min")
        .groupby("stay_id")
        .first()
        .reset_index()
        .rename(columns={"valuenum": "cr_earliest"})
    )
    prev = cr_pre_ok.merge(
        cr_earliest[["stay_id", "cr_earliest"]], on="stay_id", how="left"
    )
    prev["prevalent"] = (prev.cr_pre / prev.cr_earliest.clip(lower=0.1) >= 1.5).astype(
        int
    )
    n_prev = prev.prevalent.sum()
    consort["excl_prevalent_aki"] = n_prev
    print(f"  Prevalent AKI: {pct(n_prev, len(prev))}")

    # ── TREATED ──────────────────────────────────────────────────
    print(f"\n{'─'*50}\nBuilding treated cohort...")
    keep = set(cr_pre_ok.stay_id) - set(prev[prev.prevalent == 1].stay_id)
    treated = cardiac[cardiac.stay_id.isin(keep)].copy()
    treated = treated.merge(
        first_mg[["stay_id", "mg_offset_min", "mg_offset_h", "starttime"]].rename(
            columns={"starttime": "mg_starttime"}
        ),
        on="stay_id",
    )
    treated = treated.merge(
        cr_pre[
            [
                "stay_id",
                "cr_pre",
                "cr_pre_offset_min",
                "cr_pre_offset_h",
                "gap_to_ivmg_h",
            ]
        ],
        on="stay_id",
    )
    treated["egfr"] = compute_egfr(treated.cr_pre, treated.age, treated.is_female)

    for wname, (lo, hi) in CR_POST_WINDOWS.items():
        cand = cr_post_all[
            cr_post_all.stay_id.isin(set(treated.stay_id))
            & (cr_post_all.post_h >= lo)
            & (cr_post_all.post_h <= hi)
        ]
        fp = (
            cand.sort_values("offset_min")
            .groupby("stay_id")
            .first()
            .reset_index()[["stay_id", "valuenum", "offset_min"]]
            .rename(
                columns={
                    "valuenum": f"cr_post_{wname}",
                    "offset_min": f"cr_post_offset_{wname}",
                }
            )
        )
        treated = treated.merge(fp, on="stay_id", how="left")
        treated[f"delta_cr_{wname}"] = treated[f"cr_post_{wname}"] - treated.cr_pre

    consort["treated_final"] = len(treated)

    # Comorbidities
    for como, code_map in MIMIC_COMORB_ICD.items():
        treated[como] = treated.hadm_id.isin(
            matches_icd(dx, set(treated.hadm_id.dropna().astype(int)), code_map)
        ).astype(int)

    # Drugs before t0
    print("  Covariates (pre-t0)...")
    presc["starttime"] = pd.to_datetime(presc["starttime"], errors="coerce")
    presc_t = presc[
        presc.hadm_id.isin(set(treated.hadm_id.dropna().astype(int)))
    ].merge(treated[["hadm_id", "stay_id", "intime", "mg_starttime"]], on="hadm_id")
    presc_t["ptime"] = pd.to_datetime(presc_t["starttime"])
    presc_pre = presc_t[presc_t.ptime <= presc_t.mg_starttime]

    for dc, pats in DRUG_CLASSES.items():
        f = set(
            presc_pre[
                presc_pre.drug.str.lower().str.contains(
                    "|".join(p.lower() for p in pats), na=False
                )
            ].stay_id
        )
        treated[dc] = treated.stay_id.isin(f).astype(int)

    # Labs before t0
    for col, items in [
        ("first_mg_value", LAB_MG_MIMIC),
        ("first_potassium", LAB_K_MIMIC),
        ("first_calcium", LAB_CA_MIMIC),
        ("first_lactate", LAB_LAC_MIMIC),
    ]:
        lsub = cr_labs  # reuse merged labs with offset
        # Need fresh merge for non-Cr labs
        ll = labs[labs.itemid.isin(items) & labs.hadm_id.isin(hadms)].copy()
        ll["charttime"] = pd.to_datetime(ll["charttime"])
        ll = ll.merge(
            treated[["stay_id", "hadm_id", "intime", "mg_offset_min"]], on="hadm_id"
        )
        ll["offset_min"] = (ll.charttime - ll.intime).dt.total_seconds() / 60
        ll = ll[(ll.offset_min >= 0) & (ll.offset_min < ll.mg_offset_min)]
        fv = (
            ll.sort_values("offset_min")
            .groupby("stay_id")
            .first()
            .reset_index()[["stay_id", "valuenum"]]
            .rename(columns={"valuenum": col})
        )
        treated = treated.merge(fv, on="stay_id", how="left")
        print(f"    {col}: {pct(treated[col].notna().sum(), len(treated))}")

    # Heart rate before t0
    if len(ce_hr) > 0:
        ce = ce_hr[ce_hr.stay_id.isin(set(treated.stay_id))].copy()
        ce["charttime"] = pd.to_datetime(ce["charttime"])
        ce = ce.merge(treated[["stay_id", "intime", "mg_offset_min"]], on="stay_id")
        ce["offset_min"] = (ce.charttime - ce.intime).dt.total_seconds() / 60
        ce = ce[
            (ce.offset_min >= 0)
            & (ce.offset_min < ce.mg_offset_min)
            & ce.valuenum.between(20, 250)
        ]
        hr = (
            ce.sort_values("offset_min")
            .groupby("stay_id")
            .first()
            .reset_index()[["stay_id", "valuenum"]]
            .rename(columns={"valuenum": "first_heartrate"})
        )
        treated = treated.merge(hr, on="stay_id", how="left")
        print(
            f"    first_heartrate: {pct(treated.first_heartrate.notna().sum(), len(treated))}"
        )
    else:
        treated["first_heartrate"] = np.nan

    # Vasopressor before t0
    ie_t = ie[ie.stay_id.isin(set(treated.stay_id))].copy()
    ie_t["starttime"] = pd.to_datetime(ie_t["starttime"])
    ie_t = ie_t.merge(treated[["stay_id", "intime", "mg_starttime"]], on="stay_id")
    ie_pre = ie_t[ie_t.starttime <= ie_t.mg_starttime]
    vp = set(ie_pre[ie_pre.itemid.isin(VASO_ITEMS_MIMIC)].stay_id)
    treated["vasopressor_pre_t0"] = treated.stay_id.isin(vp).astype(int)
    print(
        f"    vasopressor_pre_t0: {pct(treated.vasopressor_pre_t0.sum(), len(treated))}"
    )

    # Transfusion before t0
    bp = set(ie_pre[ie_pre.itemid.isin(BLOOD_ITEMS_MIMIC)].stay_id)
    treated["transfusion_pre_t0"] = treated.stay_id.isin(bp).astype(int)
    print(
        f"    transfusion_pre_t0: {pct(treated.transfusion_pre_t0.sum(), len(treated))}"
    )

    # Mortality
    treated["hosp_mortality"] = treated.hospital_expire_flag.fillna(0).astype(int)

    # ── CONTROL ──────────────────────────────────────────────────
    print(f"\n{'─'*50}\nBuilding control cohort...")
    LANDMARK_H = 6  # common anchor for controls
    cr_ctrl = cr_labs[cr_labs.stay_id.isin(control_stays) & (cr_labs.offset_h >= 0)]
    c2 = set(cr_ctrl.groupby("stay_id").size().pipe(lambda s: s[s >= 2]).index)
    control = cardiac[cardiac.stay_id.isin(c2)].copy()
    consort["control_has_2cr"] = len(control)
    print(f"  Controls with ≥2 postop Cr: {pct(len(c2), len(control_stays))}")

    c1 = (
        cr_ctrl[cr_ctrl.stay_id.isin(c2)]
        .sort_values("offset_h")
        .groupby("stay_id")
        .first()
        .reset_index()
    )
    excl = set(c1[c1.valuenum >= BASELINE_CR_MAX].stay_id)
    control = control[~control.stay_id.isin(excl)].copy()
    print(f"  Excl first Cr ≥ {BASELINE_CR_MAX}: {len(excl)}")

    c1_ok = c1[c1.stay_id.isin(set(control.stay_id))]
    control = control.merge(
        c1_ok[["stay_id", "valuenum", "offset_min"]].rename(
            columns={"valuenum": "first_postop_cr", "offset_min": "first_cr_offset_min"}
        ),
        on="stay_id",
        how="left",
    )
    control["egfr"] = compute_egfr(
        control.first_postop_cr, control.age, control.is_female
    )

    # Comorbidities
    for como, code_map in MIMIC_COMORB_ICD.items():
        control[como] = control.hadm_id.isin(
            matches_icd(dx, set(control.hadm_id.dropna().astype(int)), code_map)
        ).astype(int)

    # Drugs: first 6h
    ctrl_hadms = set(control.hadm_id.dropna().astype(int))
    presc_c = presc[presc.hadm_id.isin(ctrl_hadms)].merge(
        control[["hadm_id", "stay_id", "intime"]], on="hadm_id"
    )
    presc_c["ptime"] = pd.to_datetime(presc_c["starttime"], errors="coerce")
    presc_c["off_h"] = (presc_c.ptime - presc_c.intime).dt.total_seconds() / 3600
    presc_early = presc_c[(presc_c.off_h >= 0) & (presc_c.off_h <= LANDMARK_H)]

    for dc, pats in DRUG_CLASSES.items():
        f = set(
            presc_early[
                presc_early.drug.str.lower().str.contains(
                    "|".join(p.lower() for p in pats), na=False
                )
            ].stay_id
        )
        control[dc] = control.stay_id.isin(f).astype(int)

    # Labs: first 6h
    for col, items in [
        ("first_mg_value", LAB_MG_MIMIC),
        ("first_potassium", LAB_K_MIMIC),
        ("first_calcium", LAB_CA_MIMIC),
        ("first_lactate", LAB_LAC_MIMIC),
    ]:
        ll = labs[labs.itemid.isin(items) & labs.hadm_id.isin(ctrl_hadms)].copy()
        ll["charttime"] = pd.to_datetime(ll["charttime"])
        ll = ll.merge(control[["stay_id", "hadm_id", "intime"]], on="hadm_id")
        ll["offset_h"] = (ll.charttime - ll.intime).dt.total_seconds() / 3600
        ll = ll[(ll.offset_h >= 0) & (ll.offset_h <= LANDMARK_H)]
        fv = (
            ll.sort_values("offset_h")
            .groupby("stay_id")
            .first()
            .reset_index()[["stay_id", "valuenum"]]
            .rename(columns={"valuenum": col})
        )
        control = control.merge(fv, on="stay_id", how="left")

    # HR: first 6h
    if len(ce_hr) > 0:
        ce_c = ce_hr[ce_hr.stay_id.isin(set(control.stay_id))].copy()
        ce_c["charttime"] = pd.to_datetime(ce_c["charttime"])
        ce_c = ce_c.merge(control[["stay_id", "intime"]], on="stay_id")
        ce_c["off_h"] = (ce_c.charttime - ce_c.intime).dt.total_seconds() / 3600
        ce_c = ce_c[
            (ce_c.off_h >= 0)
            & (ce_c.off_h <= LANDMARK_H)
            & ce_c.valuenum.between(20, 250)
        ]
        hr_c = (
            ce_c.sort_values("off_h")
            .groupby("stay_id")
            .first()
            .reset_index()[["stay_id", "valuenum"]]
            .rename(columns={"valuenum": "first_heartrate"})
        )
        control = control.merge(hr_c, on="stay_id", how="left")
    else:
        control["first_heartrate"] = np.nan

    # Vasopressor/transfusion: first 6h
    ie_c = ie[ie.stay_id.isin(set(control.stay_id))].copy()
    ie_c["starttime"] = pd.to_datetime(ie_c["starttime"])
    ie_c = ie_c.merge(control[["stay_id", "intime"]], on="stay_id")
    ie_c["off_h"] = (ie_c.starttime - ie_c.intime).dt.total_seconds() / 3600
    ie_early = ie_c[(ie_c.off_h >= 0) & (ie_c.off_h <= LANDMARK_H)]
    control["vasopressor_pre_t0"] = control.stay_id.isin(
        set(ie_early[ie_early.itemid.isin(VASO_ITEMS_MIMIC)].stay_id)
    ).astype(int)
    control["transfusion_pre_t0"] = control.stay_id.isin(
        set(ie_early[ie_early.itemid.isin(BLOOD_ITEMS_MIMIC)].stay_id)
    ).astype(int)

    control["hosp_mortality"] = control.hospital_expire_flag.fillna(0).astype(int)
    consort["control_final"] = len(control)

    # BMI — single chartevents pass for weight + height (both cohorts)
    print("\n  BMI extraction (chartevents)...")
    try:
        _bmi_items = {226512, 226730}
        bmi_chunks = []
        for ch in pd.read_csv(
            gz(f"{MIMIC_ICU}/chartevents.csv.gz"),
            usecols=["stay_id", "itemid", "valuenum"],
            chunksize=10_000_000,
        ):
            bmi_chunks.append(ch[ch.itemid.isin(_bmi_items) & ch.valuenum.notna()])
        bmi_df = pd.concat(bmi_chunks, ignore_index=True)
        all_bmi_stays = set(treated.stay_id) | set(control.stay_id)
        bmi_df = bmi_df[bmi_df.stay_id.isin(all_bmi_stays)]
        wt_first = (
            bmi_df[bmi_df.itemid == 226512]
            .groupby("stay_id")["valuenum"]
            .first()
            .to_dict()
        )
        ht_first = (
            bmi_df[bmi_df.itemid == 226730]
            .groupby("stay_id")["valuenum"]
            .first()
            .to_dict()
        )
        for df in [treated, control]:
            df["bmi"] = df.stay_id.apply(
                lambda x: (
                    wt_first.get(x, np.nan) / ((ht_first.get(x, np.nan) / 100) ** 2)
                    if pd.notna(wt_first.get(x))
                    and pd.notna(ht_first.get(x))
                    and ht_first.get(x, 0) > 0
                    else np.nan
                )
            )
            df.loc[~df.bmi.between(10, 80), "bmi"] = np.nan
        print(f"    Treated BMI: {pct(treated.bmi.notna().sum(), len(treated))}")
        print(f"    Control BMI: {pct(control.bmi.notna().sum(), len(control))}")
    except Exception as e:
        print(f"    BMI failed: {e}")
        treated["bmi"] = np.nan
        control["bmi"] = np.nan

    # All Cr for temporal matching
    all_s = set(treated.stay_id) | set(control.stay_id)
    cr_exp = cr_labs[cr_labs.stay_id.isin(all_s) & (cr_labs.offset_h >= 0)][
        ["stay_id", "valuenum", "offset_min", "offset_h"]
    ].copy()
    cr_exp = cr_exp.rename(
        columns={"valuenum": "labresult", "offset_min": "labresultoffset"}
    )
    cr_exp.to_csv(os.path.join(RESULTS, "20_did_cr_all_mimic.csv"), index=False)
    print(f"\n  Exported {len(cr_exp):,} Cr for {cr_exp.stay_id.nunique()} pts")

    # SMD check
    print(f"\n{'─'*50}\nCOVARIATE TIMING ASSESSMENT")
    for v in ["vasopressor_pre_t0", "transfusion_pre_t0"]:
        if v in treated.columns and v in control.columns:
            tr = treated[v].mean()
            cr_ = control[v].mean()
            smd = abs(tr - cr_) / np.sqrt((tr * (1 - tr) + cr_ * (1 - cr_)) / 2 + 1e-10)
            warn = " ⚠ consider dropping" if smd > 0.25 else ""
            print(f"  {v}: trt={100*tr:.1f}% ctrl={100*cr_:.1f}% SMD={smd:.3f}{warn}")

    build_cohort_common(treated, control, "mimic")
    print_consort(consort, treated, "MIMIC-IV")
    return consort


# ═══════════════════════════════════════════════════════════════════
#  MAIN
# ═══════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    print("=" * 70)
    print("20_did_etl.py — DiD Cohort Construction (v2: +MIMIC, +HR fix)")
    print("  Windows: 6-24h (primary), 6-48h, 0-24h")
    print("=" * 70)

    args = [a.lower() for a in sys.argv[1:]]
    run_all = len(args) == 0
    rows = []

    if run_all or "eicu" in args:
        c = run_eicu()
        if c:
            c["db"] = "eICU"
            rows.append(c)

    if run_all or "mimic" in args:
        c = run_mimic()
        if c:
            c["db"] = "MIMIC"
            rows.append(c)

    if rows:
        pd.DataFrame(rows).to_csv(
            os.path.join(RESULTS, "20_did_consort.csv"), index=False
        )
        print(f"\n  Saved: 20_did_consort.csv")

    print(f"\n{'='*70}")
    print("NEXT STEPS")
    print("=" * 70)
    print("  1. Review CONSORT numbers — check attrition")
    print("  2. Run 21_did_matching.R (1:4 PSM + temporal alignment)")
    print("  3. Run 22_did_analysis.R (DiD estimation)")
    print("=" * 70)
