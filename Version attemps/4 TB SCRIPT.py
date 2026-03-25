#Usage
"""
Description:
This script analyzes tuberculosis sequencing data to identify drug resistance
variants and antimicrobial resistance (AMR) genes.

It performs:
1. TB-Profiler analysis to detect Mycobacterium tuberculosis drug resistance variants.
2. SpoTyping to get spoligotype code.
3. Read assembly using MEGAHIT followed by AMR/virulence gene detection using
   ABRicate with the VFDB database.
4. Parse JSON output to extract detailed lineage SNPs, QC metrics, and resistance data.
5. Generate structured report matching client format (INDEX, SUMMARY, SPOLIGOTYPE,
   DRUG RESISTANCE, MUTATIONS, QC METRICS, VIRULENCE, PIPELINE, CITATIONS).

===================================================================================

Usage
-----
python3 TB_Script_Final.py <sample_dir> <sample_name> <output_prefix> <threads>

Arguments
---------
sample_dir     : Directory containing FASTQ files for the sample
sample_name    : Sample identifier
output_prefix  : Prefix used for naming final output txt file
threads        : No. of threads to use
"""

#import libraries
import os
import subprocess
import shutil
import sys
import re
import json
import pandas as pd

#############################
#Conda env for activation
#############################
conda_base = "/media/carme/Idaa/conda/anaconda3/etc/profile.d/conda.sh"
conda_env_tbprofiler = "tb_profiler_env"

abr_env_conda = "/media/carme/ganesh/miniconda3/etc/profile.d/conda.sh"
conda_env_abricate = "abricate_env"

##################
sample_dir = sys.argv[1]
sample = sys.argv[2]
output_prefix = sys.argv[3]
threads = sys.argv[4]

fastq_file = os.path.join(sample_dir, f"{sample}.fastq")

################
#Output paths###
################
tb_profiler_dir = os.path.join(sample_dir, f"{sample}_tbprofler")
megahit_out_dir = f"{sample_dir}/Megahit/"
abricate_result = f"{sample_dir}/{sample}_abricate_vfdb.tab"
vfdb_annot = "/media/carme/Research_Projects/BugSpeaks_MicrobiomeAnalysis/Tuberculosis_Firoz/VFDB/VFs_VFDB.csv"
spoligotype_output = f"{sample_dir}/{sample}_spoligotype"


######################################################################
# Lineage family abbreviation lookup
######################################################################
FAMILY_ABBREV = {
    'Indo-Oceanic': 'EAI',
    'East-Asian': 'Beijing',
    'East-African-Indian': 'CAS',
    'Euro-American': 'Euro-American',
    'West-African 1': 'Cameroon',
    'West-African 2': 'West African 2',
    'Ethiopia': 'Ethiopia',
    'La1': 'La1',
    'M. bovis': 'Bovis',
    'M. caprae': 'Caprae',
}


######################################################################
# Step 1: Run tb-profiler
######################################################################
def run_tb_profiler(conda_path, env_name):
    cmd = (
        f"source {conda_path} && "
        f"conda activate {env_name} && "
        f"tb-profiler profile --read1 {fastq_file} --platform nanopore --mapper minimap2 "
        f"--caller freebayes --prefix {sample} --dir {tb_profiler_dir} --call_whole_genome "
        f"--depth 5 --af 0.10 --suspect --txt --csv --json"
    )
    subprocess.run(cmd, shell=True, executable="/bin/bash", check=True)


######################################################################
# Step 2: Run SpoTyping for spoligotype
######################################################################
def run_spoltyping(conda_path, env_name):
    """
    Runs SpoTyping to calculate spoligotype (TB fingerprint).
    tb-profiler cannot do this - that's why JSON shows "spoligotype": null
    """
    cmd = (
        f"source {conda_path} && "
        f"conda activate {env_name} && "
        f"SpoTyping.py {fastq_file} -o {spoligotype_output}"
    )
    try:
        subprocess.run(cmd, shell=True, executable="/bin/bash", check=True)
        print("SpoTyping completed successfully")
    except subprocess.CalledProcessError:
        print("Warning: SpoTyping failed. Continuing without spoligotype.")


######################################################################
# Step 3: Parse JSON file
######################################################################
def parse_tb_json():
    """
    Reads the JSON file and extracts all detailed data.
    """
    json_file = f"{tb_profiler_dir}/results/{sample}.results.json"

    if not os.path.exists(json_file):
        print(f"Warning: JSON file not found at {json_file}")
        return None

    with open(json_file, 'r') as f:
        data = json.load(f)

    return data


######################################################################
# Step 4: Parse spoligotype result
######################################################################
def parse_spoligotype():
    """
    Reads SpoTyping output file and extracts spoligotype codes.
    """
    spol_file = f"{spoligotype_output}.txt"

    # Try different possible file names
    possible_files = [
        spol_file,
        spoligotype_output,
        f"{spoligotype_output}/spoligotype.txt",
        f"{sample_dir}/{sample}_spoligotype.txt"
    ]

    for f in possible_files:
        if os.path.exists(f):
            spol_file = f
            break
    else:
        print("Spoligotype file not found")
        return {'octal': 'N/A', 'binary': 'N/A', 'SB_number': 'N/A'}

    try:
        with open(spol_file, 'r') as f:
            lines = f.readlines()

        result = {'octal': 'N/A', 'binary': 'N/A', 'SB_number': 'N/A'}
        for line in lines:
            if line.strip() and not line.startswith('#'):
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    # SpoTyping output: col0=filename, col1=binary(43chars), col2=octal(15digits)
                    result['binary'] = parts[1] if len(parts) > 1 else 'N/A'
                    result['octal'] = parts[2] if len(parts) > 2 else 'N/A'
                    break
        return result
    except Exception as e:
        print(f"Warning: Could not parse spoligotype: {e}")
        return {'octal': 'N/A', 'binary': 'N/A', 'SB_number': 'N/A'}


######################################################################
# Step 5: Megahit (UNCHANGED)
######################################################################
def run_megahit():
    subprocess.run(f"megahit -t {threads} -r {fastq_file} -o {megahit_out_dir}", shell=True, check=True)
    shutil.move(f"{megahit_out_dir}/final.contigs.fa", f"{sample}.fa")
    shutil.rmtree(megahit_out_dir)


######################################################################
# Step 6: Abricate (UNCHANGED)
######################################################################
def run_abricate(conda_path, env_name):
    cmd = (
        f"source {conda_path} && "
        f"conda activate {env_name} && "
        f"abricate --db vfdb {sample}.fa > {abricate_result}"
    )
    subprocess.run(cmd, shell=True, executable="/bin/bash", check=True)


######################################################################
# Step 7: Merge VFDB results (UNCHANGED)
######################################################################
def extract_and_merge_vfdb_results():
    abricate_df = pd.read_csv(abricate_result, sep="\t")

    if abricate_df.empty:
        print("Warning: No virulence factors found in abricate output.")
        return pd.DataFrame()

    vfdb_df = pd.read_csv(vfdb_annot)

    abricate_df['VFID'] = abricate_df['PRODUCT'].apply(
        lambda x: re.search(r'\(VF[0-9]+\)', str(x)).group(0)[1:-1]
        if isinstance(x, str) and re.search(r'\(VF[0-9]+\)', x) else None
    )

    merged_result = pd.merge(abricate_df, vfdb_df, how='left', left_on='VFID', right_on='VFID').fillna("NA")
    return merged_result


######################################################################
# ===================================================================
# REPORT GENERATION — CLIENT FORMAT
# ===================================================================
# The following functions build the report from scratch using JSON
# data, spoligotype data, and VFDB results. No TB-Profiler .txt
# file is used as a base. This eliminates the duplication bug.
# ===================================================================
######################################################################


######################################################################
# Helper: Build position-to-variant lookup for lineage table
######################################################################
def build_variant_position_lookup(data):
    """
    Maps chromosomal positions to gene_name, ref, alt from all variant arrays.
    Used by the lineage table to get Locus/WT/Call columns.
    """
    lookup = {}
    for var in data.get('dr_variants', []) + data.get('other_variants', []):
        pos = var.get('pos')
        if pos and pos not in lookup:
            lookup[pos] = {
                'gene_name': var.get('gene_name', '-'),
                'ref': var.get('ref', '-'),
                'alt': var.get('alt', '-'),
            }
    return lookup


######################################################################
# Helper: Categorize drugs into Resistant/Susceptible/Unknown
######################################################################
def categorize_drugs(data):
    """
    Separates drugs into three categories based on confidence level.
    Returns drug-confidence pairs for the client report format.

    - Resistant: "Assoc w R" or "Assoc w R - Interim" confidence
    - Unknown: "Uncertain significance" or "Indeterminate"
    - Susceptible: Everything else (shown as "No mutation")
    """
    all_drugs = {
        'rifampicin', 'isoniazid', 'ethambutol', 'pyrazinamide',
        'moxifloxacin', 'levofloxacin', 'ofloxacin',
        'bedaquiline', 'delamanid', 'linezolid', 'streptomycin',
        'amikacin', 'kanamycin', 'capreomycin', 'clofazimine',
        'ethionamide', 'para-aminosalicylic_acid', 'cycloserine'
    }

    # Track the best (highest-ranked) confidence per drug
    confidence_rank = {
        'not assoc w r': 0, 'not assoc w r - interim': 0,
        'uncertain significance': 1, 'indeterminate': 1,
        'assoc w r - interim': 2, 'assoc w r': 3,
    }

    drug_best = {}  # drug -> (rank, confidence_string)

    # Check dr_variants
    for var in data.get('dr_variants', []):
        for annot in var.get('annotation', []):
            drug = annot.get('drug', '').lower()
            conf = annot.get('confidence', '')
            if not drug:
                continue
            rank = confidence_rank.get(conf.lower(), 1)
            if drug not in drug_best or rank > drug_best[drug][0]:
                drug_best[drug] = (rank, conf)

    # Check other_variants for unknown
    for var in data.get('other_variants', []):
        for annot in var.get('annotation', []):
            drug = annot.get('drug', '').lower()
            conf = annot.get('confidence', '')
            if not drug:
                continue
            rank = confidence_rank.get(conf.lower(), 1)
            if drug not in drug_best or rank > drug_best[drug][0]:
                drug_best[drug] = (rank, conf)

    resistant = []
    unknown = []
    susceptible = []

    for drug in sorted(all_drugs):
        if drug in drug_best:
            rank, conf = drug_best[drug]
            if rank >= 2:  # Assoc w R or Assoc w R - Interim
                resistant.append((drug, conf))
            elif rank == 1:  # Uncertain / Indeterminate
                unknown.append((drug, conf))
            else:  # Not assoc w R
                susceptible.append((drug, 'No mutation'))
        else:
            susceptible.append((drug, 'No mutation'))

    return {
        'resistant': resistant,
        'susceptible': susceptible,
        'unknown': unknown
    }


######################################################################
# Section builders
######################################################################

def build_section_header(title):
    """Standard section header matching client format."""
    return "\n" + "=" * 60 + "\n" + title + "\n" + "=" * 60 + "\n\n"


def build_index_section():
    """Builds INDEX / Table of Contents."""
    s = build_section_header("INDEX")
    sections = [
        "Summary",
        "Spoligotype",
        "Drug resistance summary",
        "Mutations - Drug Resistance",
        "Mutations - Lineage",
        "Quality Metrics - Coverage",
        "Significant Virulence Factors",
        "Significant Virulence Factors - Description",
        "Other Virulence Factors",
        "Analysis pipeline specifications and target region",
        "Citations",
    ]
    s += "Section\n"
    s += "-" * 50 + "\n"
    for i, section in enumerate(sections, 1):
        s += f"{i}.\t{section}\n"
    return s


def build_summary_section(json_data, spoligotype):
    """
    Builds SUMMARY section matching client format.
    12 fields in a two-column table.
    """
    s = build_section_header("SUMMARY")

    qc = json_data.get('qc', {})
    median_depth = qc.get('target_median_depth', 0) or 0
    pct_mapped = qc.get('percent_reads_mapped', 'N/A')

    # Get family and RD from deepest lineage entry
    lineage_list = json_data.get('lineage', [])
    family = 'N/A'
    rd = 'N/A'
    if lineage_list:
        deepest = lineage_list[-1]
        raw_family = deepest.get('family', 'N/A')
        abbrev = FAMILY_ABBREV.get(raw_family, '')
        if abbrev and abbrev != raw_family:
            family = f"{raw_family} ({abbrev})"
        else:
            family = raw_family
        rd = deepest.get('rd', 'N/A') or 'N/A'

    # Notes
    notes_list = json_data.get('notes', [])
    notes = '; '.join(notes_list) if notes_list else '-'

    # Format as tab-separated two-column table
    s += "Field\tValue\n"
    s += f"Sample name\t{json_data.get('id', sample)}\n"
    s += f"Species\tMycobacterium tuberculosis\n"
    s += f"Phylogenetic lineage\t{json_data.get('sub_lineage', 'N/A')}\n"
    s += f"Family\t{family}\n"
    s += f"Spoligotype Octal\t{spoligotype.get('octal', 'N/A')}\n"
    s += f"Spoligotype Binary\t{spoligotype.get('binary', 'N/A')}\n"
    s += f"Median Depth\t{median_depth}\n"
    s += f"% reads mapping\t{pct_mapped}\n"
    s += f"Drug resistance\t{json_data.get('drtype', 'N/A')}\n"
    s += f"Notes\t{notes}\n"
    s += f"Sublineage\t{json_data.get('sub_lineage', 'N/A')}\n"
    s += f"RD\t{rd}\n"

    # QC warning if sample failed
    sample_failed = median_depth < 5
    if sample_failed:
        s += "\n*** QC WARNING: SAMPLE FAILED ***\n"
        s += f"Median depth ({median_depth}) is below minimum threshold (5).\n"
        s += "Lineage, resistance, and variant calls are UNRELIABLE.\n"
        s += "Resistance calls (especially rrs/rrl) may be FALSE POSITIVES.\n"

    return s


def build_spoligotype_section(spoligotype):
    """Builds SPOLIGOTYPE section."""
    s = build_section_header("SPOLIGOTYPE")
    s += "Field\tValue\n"
    s += f"Spoligotype_octal_code\t{spoligotype.get('octal', 'N/A')}\n"
    s += f"Spoligotype_binary\t{spoligotype.get('binary', 'N/A')}\n"
    s += f"Spoligotype SB Number\t{spoligotype.get('SB_number', 'N/A')}\n"
    return s


def build_drug_resistance_summary_section(json_data):
    """
    Builds DRUG RESISTANCE SUMMARY with three sub-tables.
    Each drug on its own row with its confidence level.
    """
    s = build_section_header("DRUG RESISTANCE SUMMARY")
    categories = categorize_drugs(json_data)

    # Sub-table 1: Predicted Drug Resistance
    s += "Predicted Drug Resistance\n"
    s += "Drug\tConfidence\n"
    if categories['resistant']:
        for drug, conf in categories['resistant']:
            s += f"{drug}\t{conf}\n"
    else:
        s += "None\n"

    # Sub-table 2: Predicted Drug Susceptibility
    s += "\nPredicted Drug Susceptibility\n"
    s += "Drug\tConfidence\n"
    if categories['susceptible']:
        for drug, conf in categories['susceptible']:
            s += f"{drug}\t{conf}\n"
    else:
        s += "None\n"

    # Sub-table 3: Other / Unknown
    s += "\nOther / Unknown\n"
    s += "Drug\tConfidence\n"
    if categories['unknown']:
        for drug, conf in categories['unknown']:
            s += f"{drug}\t{conf}\n"
    else:
        s += "None\n"

    return s


def build_mutations_resistance_section(json_data):
    """
    Builds MUTATIONS - DRUG RESISTANCE table.
    One row per drug per variant (not one row per variant).
    Columns match client format exactly.
    """
    s = build_section_header("MUTATIONS - DRUG RESISTANCE")

    rows = []
    for var in json_data.get('dr_variants', []):
        nucleotide_change = var.get('nucleotide_change', '') or var.get('change', '')
        protein_change = var.get('protein_change', '') or ''
        freq = var.get('freq', 0)
        rel_cov = round(freq * 100, 1) if isinstance(freq, (int, float)) else freq

        for annot in var.get('annotation', []):
            drug = annot.get('drug', '')
            if not drug:
                continue
            rows.append({
                'Drug': drug,
                'Locus': var.get('gene_name', ''),
                'Locus Tag': var.get('locus_tag', ''),
                'Variant Nucleotide change': nucleotide_change,
                'Variant Amino Acid Change': protein_change,
                'Confidence': annot.get('confidence', ''),
                'Coverage': var.get('depth', ''),
                'Relative coverage (%)': rel_cov,
                'Notes': annot.get('comment', '') or ''
            })

    if rows:
        df = pd.DataFrame(rows)
        s += df.to_csv(sep="\t", index=False)
    else:
        s += "No resistance mutations detected.\n"

    return s


def build_mutations_lineage_section(json_data, pos_lookup):
    """
    Builds MUTATIONS - LINEAGE table.
    Columns: Lineage | Locus | Position | WT | Call | Coverage | Relative coverage (%)
    Uses pos_lookup to get Locus/WT/Call from variant arrays.
    """
    s = build_section_header("MUTATIONS - LINEAGE")

    rows = []
    for lineage_info in json_data.get('lineage', []):
        lineage_name = lineage_info.get('lineage', '')
        family_raw = lineage_info.get('family', '')
        abbrev = FAMILY_ABBREV.get(family_raw, '')
        if abbrev and abbrev != family_raw:
            display_lineage = f"{lineage_name} ({family_raw} / {abbrev})"
        else:
            display_lineage = f"{lineage_name} ({family_raw})" if family_raw else lineage_name

        for snp in lineage_info.get('support', []):
            pos = snp.get('pos', '')
            var_info = pos_lookup.get(pos, {})
            rows.append({
                'Lineage': display_lineage,
                'Locus': var_info.get('gene_name', '-'),
                'Position': pos,
                'WT': var_info.get('ref', '-'),
                'Call': var_info.get('alt', '-'),
                'Coverage': snp.get('target_allele_count', ''),
                'Relative coverage (%)': snp.get('target_allele_percent', '')
            })

    if rows:
        df = pd.DataFrame(rows)
        s += df.to_csv(sep="\t", index=False)
    else:
        s += "No lineage mutation data available.\n"

    return s


def build_quality_metrics_section(json_data):
    """
    Builds QUALITY METRICS - COVERAGE section.
    Two sub-tables matching client format:
      Table 1: Locus | QC | Percentage depth > 5
      Table 2: Locus | QC | Median
    """
    s = build_section_header("QUALITY METRICS - COVERAGE")

    qc_data = json_data.get('qc', {})
    target_qc = qc_data.get('target_qc', [])

    if not target_qc:
        s += "No QC data available.\n"
        return s

    # Table 1: Percentage depth > 5
    s += "Quality Metrics - Percentage depth > 5\n\n"
    s += "Locus\tQC\tPercentage depth > 5\n"
    for item in target_qc:
        pct = item.get('percent_depth_pass', 0)
        qc_status = 'PASS' if pct >= 90 else 'FAIL'
        s += f"{item.get('target', '')}\t{qc_status}\t{pct}\n"

    # Table 2: Median depth
    s += "\n\nQuality Metrics - Median depth\n\n"
    s += "Locus\tQC\tMedian\n"
    for item in target_qc:
        med = item.get('median_depth', 0)
        qc_status = 'PASS' if med >= 5 else 'FAIL'
        s += f"{item.get('target', '')}\t{qc_status}\t{med}\n"

    return s


def build_virulence_sections(has_vfdb, df1, df2, df_all):
    """
    Builds all three virulence factor sections.
    """
    s = ""

    # Significant Virulence Factors
    s += build_section_header("SIGNIFICANT VIRULENCE FACTORS")
    if has_vfdb and not df1.empty:
        s += df1.to_csv(sep="\t", index=False)
    else:
        s += "No significant virulence factors detected.\n"

    # Significant VF - Description
    s += build_section_header("SIGNIFICANT VIRULENCE FACTORS - DESCRIPTION")
    if has_vfdb and not df2.empty:
        s += df2.to_csv(sep="\t", index=False)
    else:
        s += "No descriptions available.\n"

    # Other Virulence Factors
    s += build_section_header("OTHER VIRULENCE FACTORS")
    if has_vfdb and not df_all.empty:
        s += df_all.to_csv(sep="\t", index=False)
    else:
        s += "No other virulence factors detected.\n"

    return s


def build_pipeline_section(json_data):
    """
    Builds ANALYSIS PIPELINE SPECIFICATIONS from JSON pipeline object.
    No longer parsed from .txt file.
    """
    s = build_section_header("ANALYSIS PIPELINE SPECIFICATIONS")

    s += "Task\tTool\tVersion\n"

    pipeline = json_data.get('pipeline', {})

    # TB-Profiler main entry
    tb_version = pipeline.get('software_version', 'N/A')
    s += f"resistance_prediction\ttb-profiler\t{tb_version}\n"

    # Sub-tools from pipeline.software[]
    for sw in pipeline.get('software', []):
        s += f"{sw.get('process', '')}\t{sw.get('software', '')}\t{sw.get('version', '')}\n"

    # Our additions
    s += "virulence_prediction\tabricate\t1.0.1\n"
    s += "spoligotyping\tSpoTyping\t2.1\n"

    # DB version info
    db_info = pipeline.get('db_version', {})
    if db_info:
        s += f"\nDatabase: {db_info.get('name', 'tbdb')}"
        s += f" (commit: {db_info.get('commit', 'N/A')},"
        s += f" date: {db_info.get('date', 'N/A')})\n"

    return s


def build_target_regions_section(json_data):
    """
    Builds TARGET REGIONS list from QC target names.
    """
    s = build_section_header("TARGET REGIONS")

    target_qc = json_data.get('qc', {}).get('target_qc', [])
    targets = [item.get('target', '') for item in target_qc if item.get('target')]

    if targets:
        s += f"Target genes ({len(targets)} regions):\n\n"
        # Display as comma-separated list, wrapping every 8 genes
        for i in range(0, len(targets), 8):
            chunk = targets[i:i+8]
            s += ", ".join(chunk) + "\n"
    else:
        s += "No target region data available.\n"

    return s


def build_citations_section():
    """Builds CITATIONS section with all tool references."""
    s = build_section_header("CITATIONS")

    s += ("Phelan JE, O'Sullivan DM, Machado D, et al. Integrating informatics tools and "
          "portable sequencing technology for rapid detection of resistance among tuberculosis "
          "bacteria. Genome Medicine 11, 41. 2019.\n\n")

    s += "https://github.com/tseemann/abricate\n\n"

    s += ("Zhou S, Liu B, Zheng D, Chen L, Yang J. VFDB 2025: an integrated resource for "
          "exploring anti-virulence compounds. Nucleic Acids Res. 2025.\n\n")

    s += ("Xia E, Teo YY, Ong RT. SpoTyping: fast and accurate in silico Mycobacterium "
          "spoligotyping from sequence reads. Genome Med. 2016;8(1):19.\n\n")

    s += ("Li H. Minimap2: pairwise alignment for nucleotide sequences. "
          "Bioinformatics. 2018;34(18):3094-3100.\n\n")

    s += ("Li D, Liu CM, Luo R, Sadakane K, Lam TW. MEGAHIT: an ultra-fast single-node "
          "solution for large and complex metagenomics assembly via succinct de Bruijn graph. "
          "Bioinformatics. 2015;31(10):1674-1676.\n")

    return s


######################################################################
# Step 8: Generate Final Report (CLIENT FORMAT)
######################################################################
def generate_report(merged_result, txt_dir, json_data, spoligotype):
    """
    Generates complete report matching client format.

    DOES NOT use TB-Profiler .txt file. All data comes from JSON,
    spoligotype, and VFDB results. This eliminates the duplication bug
    permanently — the report is built from scratch every time.
    """
    merged_df = merged_result

    # --- Process VFDB results ---
    has_vfdb = not merged_df.empty
    df1 = pd.DataFrame()
    df2 = pd.DataFrame()
    df_all = pd.DataFrame()

    if has_vfdb:
        # Drop columns not needed in the report
        cols_to_drop = [c for c in ['#FILE', 'RESISTANCE', 'Reference', 'COVERAGE_MAP', 'DATABASE']
                        if c in merged_df.columns]
        merged_df = merged_df.drop(columns=cols_to_drop)
        merged_df = merged_df.sort_values(by=["%COVERAGE", "%IDENTITY"], ascending=[False, False])

        # Nanopore-adjusted thresholds (80% instead of 95% for Illumina)
        sig_cov_thresh = 80
        sig_id_thresh = 80

        df1_sig = merged_df[
            (merged_df["%COVERAGE"] >= sig_cov_thresh) &
            (merged_df["%IDENTITY"] >= sig_id_thresh)
        ]
        # VF data columns (first 12)
        df1 = df1_sig.iloc[:, 0:12]
        # VF description columns (from column 11 onward, deduplicated)
        df2 = df1_sig.iloc[:, 11:].drop_duplicates()
        # Other (below threshold) — using OR logic
        df1_notsig = merged_df[
            (merged_df["%COVERAGE"] < sig_cov_thresh) |
            (merged_df["%IDENTITY"] < sig_id_thresh)
        ]
        df_all = df1_notsig.iloc[:, 0:12]

    # --- Build position lookup for lineage table ---
    pos_lookup = build_variant_position_lookup(json_data)

    # --- Output file path ---
    output_filename = f"{output_prefix}_{sample}.results.txt"
    txt_file_path = os.path.join(txt_dir, output_filename)

    # Remove previous report if it exists
    if os.path.exists(txt_file_path):
        os.remove(txt_file_path)

    # ===================================================================
    # BUILD COMPLETE REPORT FROM SCRATCH
    # ===================================================================
    report = ""

    # 1. INDEX
    report += build_index_section()

    # 2. SUMMARY
    report += build_summary_section(json_data, spoligotype)

    # 3. SPOLIGOTYPE
    report += build_spoligotype_section(spoligotype)

    # 4. DRUG RESISTANCE SUMMARY
    report += build_drug_resistance_summary_section(json_data)

    # 5. MUTATIONS - DRUG RESISTANCE
    report += build_mutations_resistance_section(json_data)

    # 6. MUTATIONS - LINEAGE
    report += build_mutations_lineage_section(json_data, pos_lookup)

    # 7. QUALITY METRICS - COVERAGE
    report += build_quality_metrics_section(json_data)

    # 8-10. VIRULENCE FACTORS
    report += build_virulence_sections(has_vfdb, df1, df2, df_all)

    # 11. ANALYSIS PIPELINE SPECIFICATIONS
    report += build_pipeline_section(json_data)

    # 12. TARGET REGIONS
    report += build_target_regions_section(json_data)

    # 13. CITATIONS
    report += build_citations_section()

    # --- Write output ---
    with open(txt_file_path, "w") as f:
        f.write(report)

    print(f"\nReport generated: {txt_file_path}")
    print(f"Output file: {txt_file_path}")


######################################################################
# MAIN FUNCTION
######################################################################
def main():
    print("\n" + "=" * 60)
    print(f"TB ANALYSIS PIPELINE - Sample: {sample}")
    print("=" * 60)

    # Step 1: Run tb-profiler (with --json)
    print("\n[1/7] Running tb-profiler...")
    run_tb_profiler(conda_base, conda_env_tbprofiler)

    # Step 2: Run SpoTyping
    print("\n[2/7] Running SpoTyping for spoligotype...")
    run_spoltyping(conda_base, conda_env_tbprofiler)

    # Step 3: Parse JSON output
    print("\n[3/7] Parsing tb-profiler JSON output...")
    json_data = parse_tb_json()

    # Step 4: Parse spoligotype
    print("\n[4/7] Parsing spoligotype results...")
    spoligotype = parse_spoligotype()

    # Step 5: Run Megahit
    print("\n[5/7] Running Megahit assembly...")
    run_megahit()

    # Step 6: Run Abricate
    print("\n[6/7] Running Abricate (VFDB)...")
    run_abricate(abr_env_conda, conda_env_abricate)

    # Step 7: Generate final report
    print("\n[7/7] Generating final report...")
    vfdb_result = extract_and_merge_vfdb_results()
    txt_dir = f"{tb_profiler_dir}/results/"
    generate_report(vfdb_result, txt_dir, json_data, spoligotype)

    print("\n" + "=" * 60)
    print("DONE!")
    print("=" * 60 + "\n")


if __name__ == "__main__":
    main()
