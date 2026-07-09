#!/usr/bin/env python3
"""
probe_como_pret0.py — Compare pre-T0 comorbidity extraction approaches (MIMIC-IV)

Purpose: Determine how to get strictly pre-T0 comorbidities for the PS model.
Tests 5 approaches for 8 comorbidities (HF, HTN, DM, CKD, COPD, PVD, stroke, liver):

  A. stay_wide   — current approach: any ICD code in diagnoses_icd for this hadm_id
  B. seq_gt1     — only seq_num > 1 (exclude primary dx)
  C. prior_adm   — ICD codes from PRIOR admissions (same subject_id, earlier admittime)
  D. union_bc    — B union C (prior + current non-primary)
  E. pmh_note    — regex from "Past Medical History" section of discharge notes

Output: prevalence table + per-patient concordance vs. stay-wide (A)

Usage:  module load Python/3.10.8-GCCcore-12.2.0
        source ~/alcrx/.venv/bin/activate
        python probe_como_pret0.py
"""

import csv
import gzip
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# ================================================================
# PATHS
# ================================================================
MIMIC_HOSP = os.path.expanduser("~/mg_aki/mimic-iv-3.1/hosp")
NOTE_PATH = (
    Path.home() / "mg_aki/physionet.org/files/mimic-iv-note/2.2/note/discharge.csv.gz"
)
RESULTS = os.path.expanduser("~/mg_aki/results")

# ================================================================
# COMORBIDITY ICD DEFINITIONS (same as 00_config.py)
# ================================================================
COMORB_ICD = {
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

# ================================================================
# PMH REGEX PATTERNS (for discharge note Past Medical History section)
# ================================================================
PMH_PATTERNS = {
    "heart_failure": [
        r"(?i)\b(?:CHF|congestive\s+heart\s+failure|heart\s+failure|HFrEF|HFpEF|"
        r"cardiomyopathy|systolic\s+dysfunction|diastolic\s+dysfunction|"
        r"reduced\s+ejection\s+fraction|EF\s*(?:of\s*)?\d{1,2}\s*%)\b"
    ],
    "hypertension": [
        r"(?i)\b(?:hypertension|HTN|high\s+blood\s+pressure|elevated\s+BP)\b"
    ],
    "diabetes": [
        r"(?i)\b(?:diabetes\s+mellitus|DM\s*(?:type|II|2|1)|IDDM|NIDDM|"
        r"insulin[\s-]dependent|type\s*[12]\s*(?:DM|diabetes)|T[12]DM|"
        r"diabetes(?:\s+type\s*[12])?)\b"
    ],
    "ckd": [
        r"(?i)\b(?:chronic\s+kidney\s+disease|CKD|chronic\s+renal\s+(?:insufficiency|failure|disease)|"
        r"CRI|stage\s*[2345]\s*(?:CKD|kidney)|GFR\s*(?:of\s*)?\d{1,2}\b)\b"
    ],
    "copd": [r"(?i)\b(?:COPD|chronic\s+obstructive|emphysema|chronic\s+bronchitis)\b"],
    "pvd": [
        r"(?i)\b(?:peripheral\s+(?:vascular|arterial)\s+disease|PVD|PAD|"
        r"claudication|arterial\s+insufficiency)\b"
    ],
    "stroke": [
        r"(?i)\b(?:(?:cerebrovascular|cerebral)\s+(?:accident|event)|CVA|stroke|TIA|"
        r"transient\s+ischemic|brain\s+(?:infarct|hemorrhage))\b"
    ],
    "liver_disease": [
        r"(?i)\b(?:cirrhosis|hepatitis\s*[BC]?|liver\s+(?:disease|failure|cirrhosis)|"
        r"hepatic\s+(?:failure|dysfunction|steatosis)|NASH|NAFLD|"
        r"fatty\s+liver|portal\s+hypertension|varic(?:es|eal))\b"
    ],
}

# Section extractor
SECTION_PAT = re.compile(
    r"^((?:Chief Complaint|Major Surgical or Invasive Procedure|"
    r"History of Present Illness|Past Medical History|"
    r"Social History|Family History|Physical Exam|"
    r"Pertinent Results|Brief Hospital Course|"
    r"Medications on Admission|Discharge Medications|"
    r"Discharge Disposition|Discharge Diagnosis|"
    r"Discharge Condition|Discharge Instructions|"
    r"Followup Instructions|Allergies|Service)\s*):?\s*$",
    re.MULTILINE | re.IGNORECASE,
)


def extract_sections(text):
    """Split discharge text into {section_name: section_text}."""
    headers = list(SECTION_PAT.finditer(text))
    sections = {}
    for i, m in enumerate(headers):
        name = m.group(1).strip().lower()
        start = m.end()
        end = headers[i + 1].start() if i + 1 < len(headers) else len(text)
        sections[name] = text[start:end].strip()
    return sections


def matches_icd(dx_sub, code_map):
    """Return set of hadm_ids matching any prefix in code_map."""
    hits = set()
    for ver, prefixes in code_map.items():
        v = dx_sub[dx_sub.icd_version == ver]
        for p in prefixes:
            hits |= set(v[v.icd_code.str.startswith(p)].hadm_id)
    return hits


# ================================================================
# MAIN
# ================================================================
def main():
    SEP = "=" * 70
    print(f"{SEP}\nprobe_como_pret0.py -- MIMIC-IV comorbidity timing audit\n{SEP}\n")

    # Load cohort (from existing pipeline output)
    cohort_path = os.path.join(RESULTS, "did_all_mimic.csv")
    if not os.path.exists(cohort_path):
        print(f"ERROR: {cohort_path} not found. Run 01_etl.py mimic first.")
        sys.exit(1)
    cohort = pd.read_csv(cohort_path)
    cohort_hadms = set(cohort.hadm_id.dropna().astype(int))
    # Map hadm_id -> subject_id for prior-admission lookup
    adm = pd.read_csv(os.path.join(MIMIC_HOSP, "admissions.csv.gz"))
    adm["admittime"] = pd.to_datetime(adm.admittime)
    # Build hadm -> subject mapping from cohort
    cohort_subj = cohort[["pid", "hadm_id"]].dropna(subset=["hadm_id"])
    cohort_subj["hadm_id"] = cohort_subj.hadm_id.astype(int)
    # pid in did_all_mimic IS stay_id; need subject_id
    icu = pd.read_csv(
        os.path.join(os.path.expanduser("~/mg_aki/mimic-iv-3.1/icu"), "icustays.csv.gz")
    )
    sid_map = icu.set_index("stay_id")["subject_id"].to_dict()
    cohort_subj["subject_id"] = cohort_subj.pid.map(sid_map)
    hadm_to_subj = cohort_subj.set_index("hadm_id")["subject_id"].to_dict()
    subject_ids = set(cohort_subj.subject_id.dropna().astype(int))

    N = len(cohort_hadms)
    print(f"Cohort: {N:,} hadm_ids, {len(subject_ids):,} subjects\n")

    # -- Load diagnoses_icd ----------------------------------------
    print("Loading diagnoses_icd...")
    dx = pd.read_csv(os.path.join(MIMIC_HOSP, "diagnoses_icd.csv.gz"))
    dx["icd_code"] = dx.icd_code.astype(str).str.strip()
    print(f"  Total rows: {len(dx):,}")

    # Diagnoses for COHORT admissions
    dx_cohort = dx[dx.hadm_id.isin(cohort_hadms)].copy()
    print(f"  Cohort rows: {len(dx_cohort):,}")
    print(f"  seq_num range: {dx_cohort.seq_num.min()} - {dx_cohort.seq_num.max()}")
    print(f"  seq_num distribution (cohort):")
    for s in [1, 2, 3, 5, 10, 20]:
        n = (dx_cohort.seq_num <= s).sum()
        print(f"    seq_num <= {s}: {n:,} ({100*n/len(dx_cohort):.1f}%)")

    # PRIOR admissions: all admissions for cohort subjects BEFORE the cohort admission
    print("\nFinding prior admissions...")
    adm_subj = adm[adm.subject_id.isin(subject_ids)].copy()
    # For each cohort hadm, find the admittime
    cohort_admit = adm[adm.hadm_id.isin(cohort_hadms)][
        ["hadm_id", "subject_id", "admittime"]
    ]
    # Prior = same subject, earlier admittime
    prior_hadms = set()
    hadm_to_prior = {}  # cohort_hadm -> set of prior hadm_ids
    for _, row in cohort_admit.iterrows():
        subj = row.subject_id
        this_time = row.admittime
        priors = set(
            adm_subj[
                (adm_subj.subject_id == subj) & (adm_subj.admittime < this_time)
            ].hadm_id
        )
        hadm_to_prior[row.hadm_id] = priors
        prior_hadms |= priors
    n_with_prior = sum(1 for v in hadm_to_prior.values() if len(v) > 0)
    print(
        f"  Cohort patients with >= 1 prior admission: {n_with_prior:,} / {N:,} ({100*n_with_prior/N:.1f}%)"
    )
    print(f"  Total prior hadm_ids: {len(prior_hadms):,}")

    dx_prior = dx[dx.hadm_id.isin(prior_hadms)].copy()
    print(f"  Prior admission dx rows: {len(dx_prior):,}")

    # -- Approach A: stay-wide (current) ---------------------------
    print(f"\n{'-'*70}")
    print("Computing 5 approaches for 8 comorbidities...")
    print(f"{'-'*70}")

    results = {}  # como -> {approach -> set of hadm_ids}

    for como, code_map in COMORB_ICD.items():
        r = {}
        # A. stay-wide (current implementation)
        r["A_stay_wide"] = matches_icd(dx_cohort, code_map)

        # B. seq_num > 1 only (exclude primary dx)
        dx_seq = dx_cohort[dx_cohort.seq_num > 1]
        r["B_seq_gt1"] = matches_icd(dx_seq, code_map)

        # C. prior admissions only
        hits_prior = set()
        for ver, prefixes in code_map.items():
            v = dx_prior[dx_prior.icd_version == ver]
            for p in prefixes:
                prior_hits = set(v[v.icd_code.str.startswith(p)].hadm_id)
                # Map prior hadm_ids back to cohort hadm_ids
                for ch, priors in hadm_to_prior.items():
                    if priors & prior_hits:
                        hits_prior.add(ch)
        r["C_prior_adm"] = hits_prior

        # D. union of B and C
        r["D_union_bc"] = r["B_seq_gt1"] | r["C_prior_adm"]

        results[como] = r

    # -- Approach E: PMH from discharge notes ----------------------
    print(f"\nExtracting PMH from discharge notes...")
    if not NOTE_PATH.exists():
        print(f"  WARN: {NOTE_PATH} not found. Skipping PMH approach.")
        for como in COMORB_ICD:
            results[como]["E_pmh_note"] = set()
    else:
        cohort_hadm_strs = {str(h) for h in cohort_hadms}
        pmh_hits = {como: set() for como in COMORB_ICD}
        n_notes = 0
        n_has_pmh = 0
        pmh_lengths = []

        with gzip.open(NOTE_PATH, "rt", encoding="utf-8", errors="replace") as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                if row["hadm_id"] not in cohort_hadm_strs:
                    continue
                n_notes += 1
                hadm = int(row["hadm_id"])
                sections = extract_sections(row["text"])
                pmh_text = sections.get("past medical history", "")
                if not pmh_text:
                    continue
                n_has_pmh += 1
                pmh_lengths.append(len(pmh_text))

                for como, patterns in PMH_PATTERNS.items():
                    for pat in patterns:
                        if re.search(pat, pmh_text):
                            pmh_hits[como].add(hadm)
                            break

                if n_notes % 2000 == 0:
                    print(f"    Processed {n_notes:,} cohort notes...")

        print(f"  Cohort notes found: {n_notes:,}")
        print(
            f"  Notes with PMH section: {n_has_pmh:,} ({100*n_has_pmh/max(n_notes,1):.1f}%)"
        )
        if pmh_lengths:
            print(
                f"  PMH length: median={int(np.median(pmh_lengths))}, "
                f"IQR=[{int(np.percentile(pmh_lengths, 25))}-{int(np.percentile(pmh_lengths, 75))}]"
            )

        for como in COMORB_ICD:
            results[como]["E_pmh_note"] = pmh_hits[como]

    # -- Print results ---------------------------------------------
    print(f"\n{'='*90}")
    print(f"  PREVALENCE COMPARISON (N = {N:,} cohort admissions)")
    print(f"{'='*90}")
    print(
        f"  {'Comorbidity':<20} {'A:stay_wide':>12} {'B:seq>1':>12} {'C:prior':>12} {'D:B+C':>12} {'E:PMH_note':>12}"
    )
    print(f"  {'-'*20} {'-'*12} {'-'*12} {'-'*12} {'-'*12} {'-'*12}")
    for como in COMORB_ICD:
        r = results[como]
        row = f"  {como:<20}"
        for app in [
            "A_stay_wide",
            "B_seq_gt1",
            "C_prior_adm",
            "D_union_bc",
            "E_pmh_note",
        ]:
            n = len(r[app])
            row += f" {n:>6} ({100*n/N:4.1f}%)"
        print(row)

    # -- Concordance: A vs D (potential contamination) -------------
    print(f"\n{'='*90}")
    print(f"  CONTAMINATION AUDIT: patients in A (stay-wide) but NOT in D (pre-T0)")
    print(
        f"  These are ICD codes that might be POST-operative complications, not pre-existing"
    )
    print(f"{'='*90}")
    print(
        f"  {'Comorbidity':<20} {'A_only':>10} {'D_only':>10} {'Both':>10} {'A_only/A':>12}"
    )
    print(f"  {'-'*20} {'-'*10} {'-'*10} {'-'*10} {'-'*12}")
    for como in COMORB_ICD:
        a = results[como]["A_stay_wide"]
        d = results[como]["D_union_bc"]
        a_only = len(a - d)
        d_only = len(d - a)
        both = len(a & d)
        frac = 100 * a_only / max(len(a), 1)
        print(f"  {como:<20} {a_only:>10} {d_only:>10} {both:>10} {frac:>10.1f}%")

    # -- STROKE deep dive ------------------------------------------
    print(f"\n{'='*90}")
    print(f"  STROKE DEEP DIVE: seq_num distribution for stroke ICD codes")
    print(f"{'='*90}")
    stroke_codes = COMORB_ICD["stroke"]
    stroke_dx = pd.DataFrame()
    for ver, prefixes in stroke_codes.items():
        v = dx_cohort[dx_cohort.icd_version == ver]
        for p in prefixes:
            stroke_dx = pd.concat([stroke_dx, v[v.icd_code.str.startswith(p)]])
    if len(stroke_dx) > 0:
        print(f"  Total stroke ICD rows in cohort: {len(stroke_dx):,}")
        print(f"  Unique hadm_ids with stroke code: {stroke_dx.hadm_id.nunique():,}")
        print(f"  seq_num distribution:")
        for s in sorted(stroke_dx.seq_num.unique()):
            n = (stroke_dx.seq_num == s).sum()
            print(f"    seq_num={s}: {n:,} rows")
        # Primary diagnosis stroke (seq_num=1) -- likely post-op stroke
        prim_stroke = set(stroke_dx[stroke_dx.seq_num == 1].hadm_id)
        prior_stroke = results["stroke"]["C_prior_adm"]
        prim_only = prim_stroke - prior_stroke
        print(f"\n  Primary dx stroke (seq=1): {len(prim_stroke):,} hadm_ids")
        print(
            f"  Of these, NOT in prior admissions: {len(prim_only):,} (likely post-op stroke)"
        )
    else:
        print("  No stroke ICD codes found in cohort")

    # -- Recommendation --------------------------------------------
    print(f"\n{'='*90}")
    print("  RECOMMENDATION")
    print(f"{'='*90}")
    print("""
  Compare approaches D (prior_adm + seq>1) vs E (PMH notes) vs A (current).

  Key question per comorbidity:
    - If A_only is small (<5% of A), stay-wide is acceptable (low contamination)
    - If A_only is large, use D or E for that comorbidity
    - For STROKE specifically, seq_num=1 stroke codes are likely post-op -> must exclude
    - PMH note approach (E) is cleanest but limited to ~76% note coverage

  Hybrid recommendation:
    - Use D (prior + seq>1) as primary for all 8 comorbidities
    - Cross-validate against E (PMH notes) for the ~76% with notes
    - Report concordance in supplement
""")

    # Save detailed results
    out_path = os.path.join(RESULTS, "probe_como_pret0.csv")
    rows = []
    for como in COMORB_ICD:
        for app_name, app_set in results[como].items():
            rows.append(
                {
                    "comorbidity": como,
                    "approach": app_name,
                    "n_positive": len(app_set),
                    "prevalence_pct": round(100 * len(app_set) / N, 2),
                }
            )
    pd.DataFrame(rows).to_csv(out_path, index=False)
    print(f"  Saved: {out_path}")


if __name__ == "__main__":
    main()
