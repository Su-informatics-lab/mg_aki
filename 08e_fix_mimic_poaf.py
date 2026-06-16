#!/usr/bin/env python3
"""
08e_fix_mimic_poaf.py — Fix MIMIC POAF phenotype

Problem: 01b_mimic_etl.py uses the same ICD codes on the same
diagnoses table for both pre-existing AF and current-admission AF,
yielding 0 POAF events.

Fix: Pre-existing AF = AF ICD codes in ANY PRIOR admission
     New-onset POAF  = AF ICD codes in CURRENT admission
                       AND no AF in any prior admission

Also adds: postop antiarrhythmic initiation as a secondary POAF signal.

Run: python 08e_fix_mimic_poaf.py
"""

import os
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

MIMIC = os.path.expanduser("~/mg_aki/mimic-iv-3.1")
RESULTS = os.path.expanduser("~/mg_aki/results")
HOSP = os.path.join(MIMIC, "hosp")
ICU = os.path.join(MIMIC, "icu")


def gz(path):
    return path if os.path.exists(path) else path.replace(".csv.gz", ".csv")


# AF ICD codes
AF_ICD9 = ["42731"]
AF_ICD10 = ["I48"]  # I480, I481, I482, I4891, etc.

# Postop antiarrhythmic drugs (new initiation = POAF signal)
ANTIARR_DRUGS = [
    "amiodarone",
    "sotalol",
    "flecainide",
    "propafenone",
    "ibutilide",
    "dofetilide",
    "dronedarone",
]
# Amiodarone in inputevents
AMIO_ITEMS = [221347, 228339, 229654, 230034]


def main():
    print("=" * 60)
    print("08e: Fix MIMIC POAF Phenotype")
    print("=" * 60)

    # ── Load cohort ──────────────────────────────────────────────
    cohort = pd.read_csv(os.path.join(RESULTS, "04_mimic_cohort.csv"))
    cohort["intime"] = pd.to_datetime(cohort["intime"])
    print(f"  Cohort: {len(cohort)} patients")

    hadms = set(cohort.hadm_id.dropna().astype(int))
    subjects = set(cohort.subject_id)

    # ── Load diagnoses + admissions ──────────────────────────────
    dx = pd.read_csv(gz(f"{HOSP}/diagnoses_icd.csv.gz"))
    dx["icd_code"] = dx["icd_code"].astype(str).str.strip()
    admissions = pd.read_csv(gz(f"{HOSP}/admissions.csv.gz"))
    admissions["admittime"] = pd.to_datetime(admissions["admittime"])

    print(f"  Diagnoses: {len(dx):,}")

    # ── Identify AF codes in ALL admissions for these subjects ───
    def has_af_codes(dx_sub):
        """Return set of hadm_ids with any AF ICD code."""
        af_hadm = set()
        for ver, prefixes in [(9, AF_ICD9), (10, AF_ICD10)]:
            v = dx_sub[dx_sub.icd_version == ver]
            for p in prefixes:
                af_hadm |= set(v[v.icd_code.str.startswith(p)].hadm_id)
        return af_hadm

    # All admissions for cohort subjects (not just the cardiac surgery one)
    all_subject_admissions = admissions[admissions.subject_id.isin(subjects)].copy()
    all_subject_dx = dx[dx.hadm_id.isin(set(all_subject_admissions.hadm_id))]

    # AF codes across ALL admissions
    all_af_hadms = has_af_codes(all_subject_dx)
    print(f"  Admissions with ANY AF code: {len(all_af_hadms)}")

    # ── For each cohort patient: check prior admissions for AF ───
    print("\n  Building prior-admission AF flags...")

    # Map: subject_id → current hadm_id + admittime
    cohort_admit = cohort[["subject_id", "hadm_id"]].merge(
        admissions[["hadm_id", "admittime"]], on="hadm_id"
    )

    preexist_af_subjects = set()
    current_af_hadms = set()

    for _, row in cohort_admit.iterrows():
        sid = row.subject_id
        current_hadm = int(row.hadm_id)
        current_admit = row.admittime

        # Prior admissions for this subject (before current admission)
        prior_adm = all_subject_admissions[
            (all_subject_admissions.subject_id == sid)
            & (all_subject_admissions.admittime < current_admit)
            & (all_subject_admissions.hadm_id != current_hadm)
        ]
        prior_hadms = set(prior_adm.hadm_id)

        # Does any prior admission have AF codes?
        if prior_hadms & all_af_hadms:
            preexist_af_subjects.add(sid)

        # Does current admission have AF codes?
        if current_hadm in all_af_hadms:
            current_af_hadms.add(current_hadm)

    print(f"  Pre-existing AF (from prior admissions): {len(preexist_af_subjects)}")
    print(f"  Current-admission AF codes: {len(current_af_hadms)}")

    # ── New-onset POAF = current AF AND NOT pre-existing ─────────
    poaf_hadms = current_af_hadms - {
        int(row.hadm_id)
        for _, row in cohort[cohort.subject_id.isin(preexist_af_subjects)].iterrows()
    }

    cohort["preexisting_af_fixed"] = cohort.subject_id.isin(
        preexist_af_subjects
    ).astype(int)
    cohort["poaf_fixed"] = cohort.hadm_id.astype(int).isin(poaf_hadms).astype(int)
    # Exclude pre-existing AF from POAF
    cohort.loc[cohort.preexisting_af_fixed == 1, "poaf_fixed"] = np.nan

    n_preaf = cohort.preexisting_af_fixed.sum()
    n_poaf = cohort.poaf_fixed.sum()
    n_eligible = (cohort.preexisting_af_fixed == 0).sum()
    print(f"\n  FIXED phenotype:")
    print(f"    Pre-existing AF: {n_preaf}")
    print(f"    POAF-eligible: {n_eligible}")
    print(f"    New-onset POAF: {int(n_poaf)} ({100*n_poaf/n_eligible:.1f}%)")

    # ── Secondary: postop antiarrhythmic initiation ──────────────
    print("\n  Secondary POAF signal: postop antiarrhythmic initiation...")
    stays = set(cohort.stay_id)

    # Check prescriptions for new antiarrhythmics
    try:
        presc = pd.read_csv(
            gz(f"{HOSP}/prescriptions.csv.gz"),
            usecols=["subject_id", "hadm_id", "starttime", "drug"],
        )
        presc["starttime"] = pd.to_datetime(presc["starttime"], errors="coerce")
        presc_cohort = presc[presc.hadm_id.isin(hadms)].merge(
            cohort[["hadm_id", "stay_id", "intime"]], on="hadm_id"
        )
        presc_cohort["offset_h"] = (
            presc_cohort.starttime - presc_cohort.intime
        ).dt.total_seconds() / 3600

        # Postop antiarrhythmics (0-168h = 7 days)
        postop_aa = presc_cohort[
            (presc_cohort.offset_h > 0) & (presc_cohort.offset_h <= 168)
        ]
        aa_mask = (
            postop_aa.drug.str.lower()
            .fillna("")
            .apply(lambda d: any(a in d for a in ANTIARR_DRUGS))
        )
        postop_aa_stays = set(postop_aa[aa_mask].stay_id)

        # Pre-existing antiarrhythmic (before surgery)
        preop_aa = presc_cohort[presc_cohort.offset_h <= 0]
        preop_aa_mask = (
            preop_aa.drug.str.lower()
            .fillna("")
            .apply(lambda d: any(a in d for a in ANTIARR_DRUGS))
        )
        preop_aa_stays = set(preop_aa[preop_aa_mask].stay_id)

        # New antiarrhythmic = postop AND NOT preop
        new_aa_stays = postop_aa_stays - preop_aa_stays
        # Exclude pre-existing AF
        new_aa_stays -= set(cohort[cohort.preexisting_af_fixed == 1].stay_id)

        cohort["poaf_antiarr"] = cohort.stay_id.isin(new_aa_stays).astype(int)
        cohort.loc[cohort.preexisting_af_fixed == 1, "poaf_antiarr"] = np.nan
        n_aa = cohort.poaf_antiarr.sum()
        print(
            f"    New postop antiarrhythmic: {int(n_aa)} ({100*n_aa/n_eligible:.1f}%)"
        )
    except Exception as e:
        print(f"    Antiarrhythmic detection failed: {e}")
        cohort["poaf_antiarr"] = np.nan

    # Check inputevents for amiodarone drips
    try:
        ie = pd.read_csv(
            gz(f"{ICU}/inputevents.csv.gz"),
            usecols=["stay_id", "itemid", "starttime"],
        )
        ie["starttime"] = pd.to_datetime(ie["starttime"], errors="coerce")
        ie_cohort = ie[ie.stay_id.isin(stays)].merge(
            cohort[["stay_id", "intime"]], on="stay_id"
        )
        ie_cohort["offset_h"] = (
            ie_cohort.starttime - ie_cohort.intime
        ).dt.total_seconds() / 3600

        amio_postop = ie_cohort[
            ie_cohort.itemid.isin(AMIO_ITEMS)
            & (ie_cohort.offset_h > 0)
            & (ie_cohort.offset_h <= 168)
        ]
        amio_stays = set(amio_postop.stay_id) - set(
            cohort[cohort.preexisting_af_fixed == 1].stay_id
        )
        print(f"    Postop amiodarone drip: {len(amio_stays)}")
    except Exception as e:
        print(f"    Amiodarone detection failed: {e}")
        amio_stays = set()

    # ── Composite POAF: ICD OR new antiarrhythmic OR amiodarone ──
    poaf_composite_stays = (
        set(cohort[cohort.poaf_fixed == 1].stay_id) | new_aa_stays | amio_stays
    )
    cohort["poaf_composite_fixed"] = cohort.stay_id.isin(poaf_composite_stays).astype(
        int
    )
    cohort.loc[cohort.preexisting_af_fixed == 1, "poaf_composite_fixed"] = np.nan
    n_comp = cohort.poaf_composite_fixed.sum()
    print(
        f"\n    Composite POAF (ICD + antiarrhythmic + amiodarone): "
        f"{int(n_comp)} ({100*n_comp/n_eligible:.1f}%)"
    )

    # ── Compare old vs new ───────────────────────────────────────
    print(f"\n  ── Old vs New ──")
    print(
        f"    Old preexisting_af: {cohort.preexisting_af.sum() if 'preexisting_af' in cohort.columns else 'N/A'}"
    )
    print(f"    New preexisting_af: {n_preaf}")
    print(f"    Old POAF: {cohort.poaf.sum() if 'poaf' in cohort.columns else 'N/A'}")
    print(f"    New POAF (ICD): {int(n_poaf)}")
    print(f"    New POAF (composite): {int(n_comp)}")

    # ── By treatment group ───────────────────────────────────────
    print(f"\n  ── POAF by Treatment Group ──")
    mg_supp_col = (
        "mg_supplementation" if "mg_supplementation" in cohort.columns else "mg_supp"
    )
    if mg_supp_col not in cohort.columns:
        for c in cohort.columns:
            if "mg_supp" in c.lower():
                mg_supp_col = c
                break

    elig = cohort[cohort.preexisting_af_fixed == 0]
    for pc in ["poaf_fixed", "poaf_composite_fixed"]:
        if pc in elig.columns:
            trt = elig[elig[mg_supp_col] == 1]
            ctrl = elig[elig[mg_supp_col] == 0]
            print(f"    {pc}:")
            print(
                f"      Supplemented: {trt[pc].sum()}/{len(trt)} ({100*trt[pc].mean():.1f}%)"
            )
            print(
                f"      Not supplemented: {ctrl[pc].sum()}/{len(ctrl)} ({100*ctrl[pc].mean():.1f}%)"
            )

    # ── Save ─────────────────────────────────────────────────────
    # Update the cohort CSV with fixed columns
    out_cols = [
        "stay_id",
        "preexisting_af_fixed",
        "poaf_fixed",
        "poaf_antiarr",
        "poaf_composite_fixed",
    ]
    patch = cohort[out_cols]
    patch.to_csv(os.path.join(RESULTS, "08e_mimic_poaf_patch.csv"), index=False)
    print(f"\n  Saved: 08e_mimic_poaf_patch.csv ({len(patch)} rows)")

    # Also save the updated full cohort
    # Replace old columns
    if "preexisting_af" in cohort.columns:
        cohort["preexisting_af_old"] = cohort["preexisting_af"]
        cohort["preexisting_af"] = cohort["preexisting_af_fixed"]
    if "poaf" in cohort.columns:
        cohort["poaf_old"] = cohort["poaf"]
        cohort["poaf"] = cohort["poaf_fixed"]
    cohort["poaf_composite"] = cohort["poaf_composite_fixed"]

    outpath = os.path.join(RESULTS, "04_mimic_cohort.csv")
    cohort.to_csv(outpath, index=False)
    print(f"  Updated: {outpath}")

    print(f"\n  Next: re-run Rscript 08d_poaf.R with fixed phenotype")


if __name__ == "__main__":
    main()
