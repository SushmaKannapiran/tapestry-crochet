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
python3 TB_Script_Final.py <sample_dir> <sample_name> <threads>

Arguments
---------
sample_dir     : Directory containing FASTQ files for the sample
sample_name    : Sample identifier
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
threads = sys.argv[3]
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
        f"python2 /media/carme/Idaa/Scripts/Tuberculosis_Report/Jan_2026/SpoTyping-v2.1-commandLine/SpoTyping.py {fastq_file} -O {spoligotype_output} -o {sample}"
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
    # Build fastq base name (without .fastq extension) for SpoTyping default output
    fastq_base = os.path.splitext(os.path.basename(fastq_file))[0]

    # Try different possible file names/locations
    possible_files = [
        f"{spoligotype_output}/{sample}",
        f"{spoligotype_output}/{fastq_base}",
        f"{spoligotype_output}.txt",
        f"{sample_dir}/{sample}_spoligotype.txt",
        f"{sample_dir}/{fastq_base}",
        fastq_base,
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
# H37Rv GENE REFERENCE DATA (NC_000962.3)
# ===================================================================
# Coordinates from the M. tuberculosis H37Rv reference genome GFF.
# Used for gene size, start/stop positions in QC and variant tables.
# ===================================================================
######################################################################

H37RV_GENES = {
    "PPE35": {"start": 2167649, "stop": 2170612, "size": 2964, "strand": "-", "locus_tag": "Rv1918c"},
    "Rv0010c": {"start": 13133, "stop": 13558, "size": 426, "strand": "-", "locus_tag": "Rv0010c"},
    "Rv0565c": {"start": 656010, "stop": 657470, "size": 1461, "strand": "-", "locus_tag": "Rv0565c"},
    "Rv1129c": {"start": 1253074, "stop": 1254534, "size": 1461, "strand": "-", "locus_tag": "Rv1129c"},
    "Rv1258c": {"start": 1406081, "stop": 1407340, "size": 1260, "strand": "-", "locus_tag": "Rv1258c"},
    "Rv1979c": {"start": 2221719, "stop": 2223164, "size": 1446, "strand": "-", "locus_tag": "Rv1979c"},
    "Rv2477c": {"start": 2782366, "stop": 2784042, "size": 1677, "strand": "-", "locus_tag": "Rv2477c"},
    "Rv2680": {"start": 2996105, "stop": 2996737, "size": 633, "strand": "+", "locus_tag": "Rv2680"},
    "Rv2681": {"start": 2996739, "stop": 2998055, "size": 1317, "strand": "+", "locus_tag": "Rv2681"},
    "Rv2752c": {"start": 3064515, "stop": 3066191, "size": 1677, "strand": "-", "locus_tag": "Rv2752c"},
    "Rv3083": {"start": 3448504, "stop": 3449991, "size": 1488, "strand": "+", "locus_tag": "Rv3083"},
    "Rv3236c": {"start": 3611959, "stop": 3613116, "size": 1158, "strand": "-", "locus_tag": "Rv3236c"},
    "aftB": {"start": 4266953, "stop": 4268836, "size": 1884, "strand": "-", "locus_tag": "Rv3805c"},
    "ahpC": {"start": 2726193, "stop": 2726780, "size": 588, "strand": "+", "locus_tag": "Rv2428"},
    "ald": {"start": 3086820, "stop": 3087935, "size": 1116, "strand": "+", "locus_tag": "Rv2780"},
    "alr": {"start": 3840194, "stop": 3841420, "size": 1227, "strand": "-", "locus_tag": "Rv3423c"},
    "atpE": {"start": 1461045, "stop": 1461290, "size": 246, "strand": "+", "locus_tag": "Rv1305"},
    "bacA": {"start": 2062809, "stop": 2064728, "size": 1920, "strand": "-", "locus_tag": "Rv1819c"},
    "ccsA": {"start": 619891, "stop": 620865, "size": 975, "strand": "+", "locus_tag": "Rv0529"},
    "clpC1": {"start": 4038158, "stop": 4040704, "size": 2547, "strand": "-", "locus_tag": "Rv3596c"},
    "ddn": {"start": 3986844, "stop": 3987299, "size": 456, "strand": "+", "locus_tag": "Rv3547"},
    "dnaA": {"start": 1, "stop": 1524, "size": 1524, "strand": "+", "locus_tag": "Rv0001"},
    "eis": {"start": 2714124, "stop": 2715332, "size": 1209, "strand": "-", "locus_tag": "Rv2416c"},
    "embA": {"start": 4243233, "stop": 4246517, "size": 3285, "strand": "+", "locus_tag": "Rv3794"},
    "embB": {"start": 4246514, "stop": 4249810, "size": 3297, "strand": "+", "locus_tag": "Rv3795"},
    "embC": {"start": 4239863, "stop": 4243147, "size": 3285, "strand": "+", "locus_tag": "Rv3793"},
    "embR": {"start": 1416181, "stop": 1417347, "size": 1167, "strand": "-", "locus_tag": "Rv1267c"},
    "ethA": {"start": 4326004, "stop": 4327473, "size": 1470, "strand": "-", "locus_tag": "Rv3854c"},
    "ethR": {"start": 4327549, "stop": 4328199, "size": 651, "strand": "+", "locus_tag": "Rv3855"},
    "fbiA": {"start": 3640543, "stop": 3641538, "size": 996, "strand": "+", "locus_tag": "Rv3261"},
    "fbiB": {"start": 3641535, "stop": 3642881, "size": 1347, "strand": "+", "locus_tag": "Rv3262"},
    "fbiC": {"start": 1302931, "stop": 1305501, "size": 2571, "strand": "+", "locus_tag": "Rv1173"},
    "fbiD": {"start": 3339118, "stop": 3339762, "size": 645, "strand": "+", "locus_tag": "Rv2983"},
    "fgd1": {"start": 490783, "stop": 491793, "size": 1011, "strand": "+", "locus_tag": "Rv0407"},
    "folC": {"start": 2746135, "stop": 2747598, "size": 1464, "strand": "-", "locus_tag": "Rv2447c"},
    "gid": {"start": 4407528, "stop": 4408202, "size": 675, "strand": "-", "locus_tag": "Rv3919c"},
    "glpK": {"start": 4138202, "stop": 4139755, "size": 1554, "strand": "-", "locus_tag": "Rv3696c"},
    "gyrA": {"start": 7302, "stop": 9818, "size": 2517, "strand": "+", "locus_tag": "Rv0006"},
    "gyrB": {"start": 5240, "stop": 7267, "size": 2028, "strand": "+", "locus_tag": "Rv0005"},
    "hadA": {"start": 731930, "stop": 732406, "size": 477, "strand": "+", "locus_tag": "Rv0635"},
    "inhA": {"start": 1674202, "stop": 1675011, "size": 810, "strand": "+", "locus_tag": "Rv1484"},
    "kasA": {"start": 2518115, "stop": 2519365, "size": 1251, "strand": "+", "locus_tag": "Rv2245"},
    "katG": {"start": 2153889, "stop": 2156111, "size": 2223, "strand": "-", "locus_tag": "Rv1908c"},
    "lpqB": {"start": 3623159, "stop": 3624910, "size": 1752, "strand": "-", "locus_tag": "Rv3244c"},
    "mmaA3": {"start": 737268, "stop": 738149, "size": 882, "strand": "-", "locus_tag": "Rv0643c"},
    "mmpL5": {"start": 775586, "stop": 778480, "size": 2895, "strand": "-", "locus_tag": "Rv0676c"},
    "mmpR5": {"start": 778990, "stop": 779487, "size": 498, "strand": "+", "locus_tag": "Rv0678"},
    "mmpS5": {"start": 778477, "stop": 778905, "size": 429, "strand": "-", "locus_tag": "Rv0677c"},
    "mshA": {"start": 575348, "stop": 576790, "size": 1443, "strand": "+", "locus_tag": "Rv0486"},
    "mtrA": {"start": 3626663, "stop": 3627349, "size": 687, "strand": "-", "locus_tag": "Rv3246c"},
    "mtrB": {"start": 3624910, "stop": 3626613, "size": 1704, "strand": "-", "locus_tag": "Rv3245c"},
    "ndh": {"start": 2101651, "stop": 2103042, "size": 1392, "strand": "-", "locus_tag": "Rv1854c"},
    "nusG": {"start": 734254, "stop": 734970, "size": 717, "strand": "+", "locus_tag": "Rv0639"},
    "panD": {"start": 4043862, "stop": 4044281, "size": 420, "strand": "-", "locus_tag": "Rv3601c"},
    "pepQ": {"start": 2859300, "stop": 2860418, "size": 1119, "strand": "-", "locus_tag": "Rv2535c"},
    "pncA": {"start": 2288681, "stop": 2289241, "size": 561, "strand": "-", "locus_tag": "Rv2043c"},
    "ribD": {"start": 2986839, "stop": 2987615, "size": 777, "strand": "+", "locus_tag": "Rv2671"},
    "rplC": {"start": 800809, "stop": 801462, "size": 654, "strand": "+", "locus_tag": "Rv0701"},
    "rpoA": {"start": 3877464, "stop": 3878507, "size": 1044, "strand": "-", "locus_tag": "Rv3457c"},
    "rpoB": {"start": 759807, "stop": 763325, "size": 3519, "strand": "+", "locus_tag": "Rv0667"},
    "rpoC": {"start": 763370, "stop": 767320, "size": 3951, "strand": "+", "locus_tag": "Rv0668"},
    "rpsA": {"start": 1833542, "stop": 1834987, "size": 1446, "strand": "+", "locus_tag": "Rv1630"},
    "rpsL": {"start": 781560, "stop": 781934, "size": 375, "strand": "+", "locus_tag": "Rv0682"},
    "rrl": {"start": 1473658, "stop": 1476795, "size": 3138, "strand": "+", "locus_tag": "EBG00000313339"},
    "rrs": {"start": 1471846, "stop": 1473382, "size": 1537, "strand": "+", "locus_tag": "EBG00000313325"},
    "sigE": {"start": 1364413, "stop": 1365186, "size": 774, "strand": "+", "locus_tag": "Rv1221"},
    "thyA": {"start": 3073680, "stop": 3074471, "size": 792, "strand": "-", "locus_tag": "Rv2764c"},
    "thyX": {"start": 3067193, "stop": 3067945, "size": 753, "strand": "-", "locus_tag": "Rv2754c"},
    "tlyA": {"start": 1917940, "stop": 1918746, "size": 807, "strand": "+", "locus_tag": "Rv1694"},
    "tsnR": {"start": 1853606, "stop": 1854388, "size": 783, "strand": "+", "locus_tag": "Rv1644"},
    "ubiA": {"start": 4268925, "stop": 4269833, "size": 909, "strand": "-", "locus_tag": "Rv3806c"},
    "whiB6": {"start": 4338171, "stop": 4338521, "size": 351, "strand": "-", "locus_tag": "Rv3862c"},
    "whiB7": {"start": 3568401, "stop": 3568679, "size": 279, "strand": "-", "locus_tag": "Rv3197A"},
}

# Map locus_tag → gene_name for reverse lookup
LOCUS_TAG_TO_GENE = {v['locus_tag']: k for k, v in H37RV_GENES.items()}


######################################################################
# ===================================================================
# REPORT GENERATION — CLIENT FORMAT (V5.0)
# ===================================================================
######################################################################


def build_section_header(title):
    """Standard section header."""
    return "\n" + "=" * 60 + "\n" + title + "\n" + "=" * 60 + "\n\n"


def build_variant_position_lookup(data):
    """Maps chromosomal positions to variant details from all variant arrays."""
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


def get_variant_type(var):
    """Determines variant type label matching client format."""
    vtype = var.get('type', '')
    ref = var.get('ref', '')
    alt = var.get('alt', '')

    if var.get('sv', False):
        return 'Deletion' if (var.get('sv_len') or 0) < 0 else 'Structural Variant'
    if len(ref) > len(alt):
        return 'Deletion'
    if len(alt) > len(ref):
        return 'Insertion'
    return 'Mutation'


def get_position_in_gene(var):
    """Extracts position within gene from nucleotide_change (e.g., c.873G>C → 873)."""
    nc = var.get('nucleotide_change', '')
    if not nc:
        return ''
    # Match patterns like c.873G>C, c.-777C>T, c.1234delA, c.1234_1235insA
    m = re.search(r'c\.(-?\d+)', nc)
    if m:
        return m.group(1)
    return ''


def get_strand_info(var):
    """Determines strand from forward/reverse read counts."""
    fwd = var.get('forward_reads', 0)
    rev = var.get('reverse_reads', 0)
    if fwd > 0 and rev > 0:
        return '+/-'
    elif fwd > 0:
        return '+'
    elif rev > 0:
        return '-'
    return 'N/A'


def get_gene_drugs(gene_name):
    """Gets drug associations for a gene from H37RV_GENES."""
    # Map gene names to their associated drugs
    GENE_DRUGS = {
        'rpoB': 'rifampicin', 'rpoC': 'rifampicin',
        'katG': 'isoniazid', 'inhA': 'isoniazid, ethionamide', 'ahpC': 'isoniazid',
        'kasA': 'isoniazid', 'ndh': 'isoniazid',
        'embB': 'ethambutol', 'embA': 'ethambutol', 'embC': 'ethambutol', 'embR': 'ethambutol',
        'pncA': 'pyrazinamide', 'panD': 'pyrazinamide', 'rpsA': 'pyrazinamide', 'clpC1': 'pyrazinamide',
        'gyrA': 'fluoroquinolones', 'gyrB': 'fluoroquinolones',
        'rrs': 'aminoglycosides', 'rrl': 'linezolid, capreomycin',
        'rpsL': 'streptomycin', 'gid': 'streptomycin',
        'eis': 'kanamycin, amikacin',
        'rplC': 'linezolid',
        'ethA': 'ethionamide', 'ethR': 'ethionamide',
        'mmpR5': 'bedaquiline, clofazimine', 'mmpL5': 'bedaquiline, clofazimine',
        'mmpS5': 'bedaquiline, clofazimine',
        'atpE': 'bedaquiline',
        'pepQ': 'bedaquiline, clofazimine',
        'fbiA': 'delamanid', 'fbiB': 'delamanid', 'fbiC': 'delamanid',
        'fbiD': 'delamanid', 'fgd1': 'delamanid, clofazimine', 'ddn': 'delamanid',
        'folC': 'para-aminosalicylic_acid', 'ribD': 'para-aminosalicylic_acid',
        'thyA': 'para-aminosalicylic_acid', 'thyX': 'para-aminosalicylic_acid',
        'alr': 'cycloserine',
        'tlyA': 'capreomycin',
        'sigE': 'multiple', 'whiB6': 'multiple', 'whiB7': 'multiple',
        'ubiA': 'ethambutol',
    }
    return GENE_DRUGS.get(gene_name, '-')


def categorize_drugs(data):
    """
    Separates drugs into three categories with confidence levels.
    Returns drug-confidence pairs for the client report format.
    """
    all_drugs = {
        'rifampicin', 'isoniazid', 'ethambutol', 'pyrazinamide',
        'moxifloxacin', 'levofloxacin', 'ofloxacin',
        'bedaquiline', 'delamanid', 'linezolid', 'streptomycin',
        'amikacin', 'kanamycin', 'capreomycin', 'clofazimine',
        'ethionamide', 'para-aminosalicylic_acid', 'cycloserine'
    }

    confidence_rank = {
        'not assoc w r': 0, 'not assoc w r - interim': 0,
        'uncertain significance': 1, 'indeterminate': 1,
        'assoc w r - interim': 2, 'assoc w r': 3,
    }

    drug_best = {}

    for var in data.get('dr_variants', []):
        for annot in var.get('annotation', []):
            drug = annot.get('drug', '').lower()
            conf = annot.get('confidence', '')
            if not drug:
                continue
            rank = confidence_rank.get(conf.lower(), 1)
            if drug not in drug_best or rank > drug_best[drug][0]:
                drug_best[drug] = (rank, conf)

    for var in data.get('other_variants', []):
        for annot in var.get('annotation', []):
            drug = annot.get('drug', '').lower()
            conf = annot.get('confidence', '')
            if not drug:
                continue
            rank = confidence_rank.get(conf.lower(), 1)
            if drug not in drug_best or rank > drug_best[drug][0]:
                drug_best[drug] = (rank, conf)

    resistant, unknown, susceptible = [], [], []

    for drug in sorted(all_drugs):
        if drug in drug_best:
            rank, conf = drug_best[drug]
            if rank >= 2:
                resistant.append((drug, conf))
            elif rank == 1:
                unknown.append((drug, conf))
            else:
                susceptible.append((drug, 'No mutation'))
        else:
            susceptible.append((drug, 'No mutation'))

    return {'resistant': resistant, 'susceptible': susceptible, 'unknown': unknown}


######################################################################
# Section builders
######################################################################

def build_index_section():
    """Builds INDEX / Table of Contents."""
    s = build_section_header("INDEX")
    sections = [
        "Summary",
        "Spoligotype",
        "Drug resistance summary",
        "Genomic variants resistance",
        "Mutations - Lineage",
        "Quality metrics I coverage",
        "Quality metrics II resistance",
        "Quality metrics III resistance",
        "Significant Virulence Factors",
        "Significant Virulence Factors - Description",
        "Other Virulence Factors",
        "Analysis pipeline specifications and target region",
        "Info",
        "Citations",
    ]
    s += "Section\n"
    s += "-" * 50 + "\n"
    for i, section in enumerate(sections, 1):
        s += f"{i}.\t{section}\n"
    return s


def build_summary_section(json_data, spoligotype):
    """Builds SUMMARY section — 12 fields matching client format."""
    s = build_section_header("SUMMARY")

    qc = json_data.get('qc', {})
    median_depth = qc.get('target_median_depth', 0) or 0
    pct_mapped = qc.get('percent_reads_mapped', 'N/A')

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

    notes_list = json_data.get('notes', [])
    notes = '; '.join(notes_list) if notes_list else '-'

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

    if median_depth < 5:
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
    """Builds DRUG RESISTANCE SUMMARY — three sub-tables with Drug | Confidence."""
    s = build_section_header("DRUG RESISTANCE SUMMARY")
    categories = categorize_drugs(json_data)

    s += "Predicted Drug Resistance\n"
    s += "Drug\tConfidence\n"
    if categories['resistant']:
        for drug, conf in categories['resistant']:
            s += f"{drug}\t{conf}\n"
    else:
        s += "None\n"

    s += "\nPredicted Drug Susceptibility\n"
    s += "Drug\tConfidence\n"
    if categories['susceptible']:
        for drug, conf in categories['susceptible']:
            s += f"{drug}\t{conf}\n"
    else:
        s += "None\n"

    s += "\nOther / Unknown\n"
    s += "Drug\tConfidence\n"
    if categories['unknown']:
        for drug, conf in categories['unknown']:
            s += f"{drug}\t{conf}\n"
    else:
        s += "None\n"

    return s


def build_genomic_variants_section(json_data):
    """
    Builds GENOMIC VARIANTS RESISTANCE table — ALL variants
    (dr_variants + other_variants), matching client format.
    Columns: Antibiotic | Locus | Gene | Size | Position | Genomic Position |
             Type | WT | Call | Relative coverage (%) | Strand
    """
    s = build_section_header("GENOMIC VARIANTS RESISTANCE")

    rows = []
    all_variants = json_data.get('dr_variants', []) + json_data.get('other_variants', [])

    for var in all_variants:
        gene_name = var.get('gene_name', '')
        gene_info = H37RV_GENES.get(gene_name, {})
        gene_size = gene_info.get('size', '-')

        pos_in_gene = get_position_in_gene(var)
        genomic_pos = var.get('pos', '')
        var_type = get_variant_type(var)
        strand = get_strand_info(var)

        freq = var.get('freq', 0)
        rel_cov = round(freq * 100, 1) if isinstance(freq, (int, float)) else freq

        # Get drug from annotation or from gene association
        drugs = set()
        for annot in var.get('annotation', []):
            d = annot.get('drug', '')
            if d:
                drugs.add(d)
        if not drugs:
            drug_str = get_gene_drugs(gene_name)
        else:
            drug_str = ', '.join(sorted(drugs))

        rows.append({
            'Antibiotic': drug_str,
            'Locus': gene_info.get('locus_tag', var.get('locus_tag', '')),
            'Gene': gene_name,
            'Size': gene_size,
            'Position': pos_in_gene,
            'Genomic Position': genomic_pos,
            'Type': var_type,
            'WT': var.get('ref', ''),
            'Call': var.get('alt', ''),
            'Relative coverage (%)': rel_cov,
            'Strand': strand,
        })

    if rows:
        df = pd.DataFrame(rows)
        # Sort by gene name then genomic position
        df = df.sort_values(by=['Gene', 'Genomic Position'])
        s += df.to_csv(sep="\t", index=False)
    else:
        s += "No genomic variants detected.\n"

    return s


def build_mutations_lineage_section(json_data, pos_lookup):
    """
    Builds MUTATIONS - LINEAGE table.
    Only shows the DEEPEST lineage entry (not all parent levels).
    """
    s = build_section_header("MUTATIONS - LINEAGE")

    lineage_list = json_data.get('lineage', [])
    if not lineage_list:
        s += "No lineage mutation data available.\n"
        return s

    # Only use the deepest lineage entry
    deepest = lineage_list[-1]
    lineage_name = deepest.get('lineage', '')
    family_raw = deepest.get('family', '')
    abbrev = FAMILY_ABBREV.get(family_raw, '')
    if abbrev and abbrev != family_raw:
        display_lineage = f"{lineage_name} ({family_raw} / {abbrev})"
    else:
        display_lineage = f"{lineage_name} ({family_raw})" if family_raw else lineage_name

    rows = []
    for snp in deepest.get('support', []):
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
    Builds QUALITY METRICS — THREE tables matching client format.

    Table I (Coverage): Antibiotic | Locus | Gene | QC | Size | Percentage depth > 5 |
                        Median Depth | Start | Stop
    Table II (Resistance): Antibiotic | Locus | Gene | QC | NA mutations | AA mutations
    Table III (Resistance): Antibiotic | Locus | Gene | QC | Number of deletions |
                            Number of insertions
    """
    target_qc = json_data.get('qc', {}).get('target_qc', [])

    if not target_qc:
        s = build_section_header("QUALITY METRICS I COVERAGE")
        s += "No QC data available.\n"
        return s

    # Pre-compute per-gene variant counts for Tables II and III
    gene_variants = {}  # gene -> counts dict
    all_variants = json_data.get('dr_variants', []) + json_data.get('other_variants', [])

    unknown_confidences = {'uncertain significance', 'indeterminate', 'not assoc w r', 'not assoc w r - interim'}

    for var in all_variants:
        gene = var.get('gene_name', '')
        if not gene:
            continue
        if gene not in gene_variants:
            gene_variants[gene] = {'na_mut': 0, 'aa_mut': 0, 'del': 0, 'ins': 0,
                                   'unk_mut': 0, 'unk_del': 0, 'unk_ins': 0}

        ref = var.get('ref', '')
        alt = var.get('alt', '')
        nc = var.get('nucleotide_change', '')
        pc = var.get('protein_change', '')

        # Determine if this variant is "unknown" based on confidence
        is_unknown = False
        for annot in var.get('annotation', []):
            conf = annot.get('confidence', '').lower()
            if conf in unknown_confidences:
                is_unknown = True
                break

        # Count nucleotide mutations
        if nc:
            gene_variants[gene]['na_mut'] += 1
        # Count amino acid mutations
        if pc:
            gene_variants[gene]['aa_mut'] += 1
        # Count unknown mutations
        if is_unknown and (nc or pc):
            gene_variants[gene]['unk_mut'] += 1
        # Count deletions and insertions
        if len(ref) > len(alt):
            gene_variants[gene]['del'] += 1
            if is_unknown:
                gene_variants[gene]['unk_del'] += 1
        elif len(alt) > len(ref):
            gene_variants[gene]['ins'] += 1
            if is_unknown:
                gene_variants[gene]['unk_ins'] += 1

    # --- Table I: Coverage ---
    s = build_section_header("QUALITY METRICS I COVERAGE")
    s += "Antibiotic\tLocus\tGene\tQC\tSize\tPercentage depth > 5\tMedian Depth\tStart\tStop\n"
    for item in target_qc:
        gene = item.get('target', '')
        gene_info = H37RV_GENES.get(gene, {})
        pct = item.get('percent_depth_pass', 0)
        med = item.get('median_depth', 0)
        qc_status = 'PASS' if pct >= 90 else 'FAIL'
        drug = get_gene_drugs(gene)
        locus = gene_info.get('locus_tag', '-')
        size = gene_info.get('size', '-')
        start = gene_info.get('start', '-')
        stop = gene_info.get('stop', '-')
        s += f"{drug}\t{locus}\t{gene}\t{qc_status}\t{size}\t{pct}\t{med}\t{start}\t{stop}\n"

    # --- Table II: Resistance mutations ---
    s += build_section_header("QUALITY METRICS II RESISTANCE")
    s += "Antibiotic\tLocus\tGene\tQC\tNA mutations\tAA mutations\n"
    for item in target_qc:
        gene = item.get('target', '')
        gene_info = H37RV_GENES.get(gene, {})
        pct = item.get('percent_depth_pass', 0)
        qc_status = 'PASS' if pct >= 90 else 'FAIL'
        drug = get_gene_drugs(gene)
        locus = gene_info.get('locus_tag', '-')
        gv = gene_variants.get(gene, {'na_mut': 0, 'aa_mut': 0})
        s += f"{drug}\t{locus}\t{gene}\t{qc_status}\t{gv['na_mut']}\t{gv['aa_mut']}\n"

    # --- Table III: Deletions/Insertions (matching client: 8 data columns) ---
    s += build_section_header("QUALITY METRICS III RESISTANCE")
    s += "Antibiotic\tLocus\tGene\tQC\tNA mutations\tAA mutations\tUnknown mutations\tNumber of deletions\tUnknown deletions\tNumber of insertions\tUnknown insertions\n"
    for item in target_qc:
        gene = item.get('target', '')
        gene_info = H37RV_GENES.get(gene, {})
        pct = item.get('percent_depth_pass', 0)
        qc_status = 'PASS' if pct >= 90 else 'FAIL'
        drug = get_gene_drugs(gene)
        locus = gene_info.get('locus_tag', '-')
        gv = gene_variants.get(gene, {'na_mut': 0, 'aa_mut': 0, 'del': 0, 'ins': 0, 'unk_mut': 0, 'unk_del': 0, 'unk_ins': 0})
        s += f"{drug}\t{locus}\t{gene}\t{qc_status}\t{gv['na_mut']}\t{gv['aa_mut']}\t{gv.get('unk_mut', 0)}\t{gv['del']}\t{gv.get('unk_del', 0)}\t{gv['ins']}\t{gv.get('unk_ins', 0)}\n"

    return s


def build_virulence_sections(has_vfdb, df1, df2, df_all):
    """Builds all three virulence factor sections."""
    s = ""

    s += build_section_header("SIGNIFICANT VIRULENCE FACTORS")
    if has_vfdb and not df1.empty:
        s += df1.to_csv(sep="\t", index=False)
    else:
        s += "No significant virulence factors detected.\n"

    s += build_section_header("SIGNIFICANT VIRULENCE FACTORS - DESCRIPTION")
    if has_vfdb and not df2.empty:
        s += df2.to_csv(sep="\t", index=False)
    else:
        s += "No descriptions available.\n"

    s += build_section_header("OTHER VIRULENCE FACTORS")
    if has_vfdb and not df_all.empty:
        s += df_all.to_csv(sep="\t", index=False)
    else:
        s += "No other virulence factors detected.\n"

    return s


def build_pipeline_section(json_data):
    """Builds ANALYSIS PIPELINE SPECIFICATIONS from JSON pipeline object."""
    s = build_section_header("ANALYSIS PIPELINE SPECIFICATIONS")

    s += "Task\tTool\tVersion\n"
    pipeline = json_data.get('pipeline', {})
    tb_version = pipeline.get('software_version', 'N/A')
    s += f"resistance_prediction\ttb-profiler\t{tb_version}\n"

    for sw in pipeline.get('software', []):
        s += f"{sw.get('process', '')}\t{sw.get('software', '')}\t{sw.get('version', '')}\n"

    s += "virulence_prediction\tabricate\t1.0.1\n"
    s += "spoligotyping\tSpoTyping\t2.1\n"
    s += "assembly\tMEGAHIT\t1.2.9\n"

    db_info = pipeline.get('db_version', {})
    if db_info:
        s += f"\nDatabase: {db_info.get('name', 'tbdb')}"
        s += f" (commit: {db_info.get('commit', 'N/A')},"
        s += f" date: {db_info.get('date', 'N/A')})\n"

    return s


def build_target_regions_section(json_data):
    """Builds TARGET REGIONS list from QC target names."""
    s = build_section_header("TARGET REGIONS")

    target_qc = json_data.get('qc', {}).get('target_qc', [])
    targets = [item.get('target', '') for item in target_qc if item.get('target')]

    if targets:
        s += f"Target genes ({len(targets)} regions):\n\n"
        for i in range(0, len(targets), 8):
            chunk = targets[i:i+8]
            s += ", ".join(chunk) + "\n"
    else:
        s += "No target region data available.\n"

    return s


def build_info_section(json_data, spoligotype):
    """
    Builds INFO section with General, Lineage, Resistance metadata
    and Disclaimer — matching client format.
    """
    s = build_section_header("INFO")

    # General
    s += "General\n"
    s += "-" * 40 + "\n"
    s += f"Sample ID:\t{json_data.get('id', sample)}\n"
    s += f"Species:\tMycobacterium tuberculosis\n"
    s += f"Analysis date:\t{json_data.get('timestamp', 'N/A')}\n"
    pipeline = json_data.get('pipeline', {})
    s += f"TB-Profiler version:\t{pipeline.get('software_version', 'N/A')}\n"
    db_info = pipeline.get('db_version', {})
    s += f"Database:\t{db_info.get('name', 'tbdb')} ({db_info.get('commit', 'N/A')})\n"
    qc = json_data.get('qc', {})
    s += f"Reads mapped:\t{qc.get('num_reads_mapped', 'N/A')}\n"
    s += f"% reads mapped:\t{qc.get('percent_reads_mapped', 'N/A')}\n"
    s += f"Median depth:\t{qc.get('target_median_depth', 'N/A')}\n"

    # Lineage
    s += "\nLineage\n"
    s += "-" * 40 + "\n"
    s += f"Main lineage:\t{json_data.get('main_lineage', 'N/A')}\n"
    s += f"Sub lineage:\t{json_data.get('sub_lineage', 'N/A')}\n"
    lineage_list = json_data.get('lineage', [])
    if lineage_list:
        deepest = lineage_list[-1]
        s += f"Family:\t{deepest.get('family', 'N/A')}\n"
        s += f"RD:\t{deepest.get('rd', 'N/A')}\n"
    s += f"Spoligotype octal:\t{spoligotype.get('octal', 'N/A')}\n"
    s += f"Spoligotype binary:\t{spoligotype.get('binary', 'N/A')}\n"

    # Resistance
    s += "\nResistance\n"
    s += "-" * 40 + "\n"
    s += f"Drug resistance type:\t{json_data.get('drtype', 'N/A')}\n"
    dr_vars = json_data.get('dr_variants', [])
    s += f"Number of resistance variants:\t{len(dr_vars)}\n"
    other_vars = json_data.get('other_variants', [])
    s += f"Number of other variants:\t{len(other_vars)}\n"

    # Disclaimer
    s += "\nDisclaimer\n"
    s += "-" * 40 + "\n"
    s += ("This report is generated for research purposes. Drug resistance predictions "
          "are based on the WHO catalogue of mutations and should be confirmed by "
          "phenotypic drug susceptibility testing (DST) before clinical decisions are made. "
          "The analysis uses the TB-Profiler bioinformatics pipeline with the tbdb mutation "
          "database. Results should be interpreted by qualified personnel in the context of "
          "clinical and epidemiological data.\n")

    return s


def build_citations_section():
    """Builds CITATIONS section."""
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
# Step 8: Generate Final Report (CLIENT FORMAT V5.0)
######################################################################
def generate_report(merged_result, txt_dir, json_data, spoligotype):
    """
    Generates complete report matching client format.
    Builds entire report from scratch using JSON data.
    No TB-Profiler .txt file is used.
    """
    merged_df = merged_result

    # --- Process VFDB results ---
    has_vfdb = not merged_df.empty
    df1 = pd.DataFrame()
    df2 = pd.DataFrame()
    df_all = pd.DataFrame()

    if has_vfdb:
        cols_to_drop = [c for c in ['#FILE', 'RESISTANCE', 'Reference', 'COVERAGE_MAP', 'DATABASE']
                        if c in merged_df.columns]
        merged_df = merged_df.drop(columns=cols_to_drop)
        merged_df = merged_df.sort_values(by=["%COVERAGE", "%IDENTITY"], ascending=[False, False])

        sig_cov_thresh = 80
        sig_id_thresh = 80

        df1_sig = merged_df[
            (merged_df["%COVERAGE"] >= sig_cov_thresh) &
            (merged_df["%IDENTITY"] >= sig_id_thresh)
        ]
        df1 = df1_sig.iloc[:, 0:12]
        df2 = df1_sig.iloc[:, 11:].drop_duplicates()
        df1_notsig = merged_df[
            (merged_df["%COVERAGE"] < sig_cov_thresh) |
            (merged_df["%IDENTITY"] < sig_id_thresh)
        ]
        df_all = df1_notsig.iloc[:, 0:12]

    # --- Build position lookup for lineage table ---
    pos_lookup = build_variant_position_lookup(json_data)

    # --- Output file path ---
    output_filename = f"{sample}_Updated.results.txt"
    txt_file_path = os.path.join(txt_dir, output_filename)

    if os.path.exists(txt_file_path):
        os.remove(txt_file_path)

    # ===================================================================
    # BUILD COMPLETE REPORT
    # ===================================================================
    report = ""

    report += build_index_section()
    report += build_summary_section(json_data, spoligotype)
    report += build_spoligotype_section(spoligotype)
    report += build_drug_resistance_summary_section(json_data)
    report += build_genomic_variants_section(json_data)
    report += build_mutations_lineage_section(json_data, pos_lookup)
    report += build_quality_metrics_section(json_data)
    report += build_virulence_sections(has_vfdb, df1, df2, df_all)
    report += build_pipeline_section(json_data)
    report += build_target_regions_section(json_data)
    report += build_info_section(json_data, spoligotype)
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
    #run_tb_profiler(conda_base, conda_env_tbprofiler)

    # Step 2: Run SpoTyping
    print("\n[2/7] Running SpoTyping for spoligotype...")
    #run_spoltyping(conda_base, conda_env_tbprofiler)

    # Step 3: Parse JSON output
    print("\n[3/7] Parsing tb-profiler JSON output...")
    json_data = parse_tb_json()

    # Step 4: Parse spoligotype
    print("\n[4/7] Parsing spoligotype results...")
    spoligotype = parse_spoligotype()

    # Step 5: Run Megahit
    print("\n[5/7] Running Megahit assembly...")
    #run_megahit()

    # Step 6: Run Abricate
    print("\n[6/7] Running Abricate (VFDB)...")
    #run_abricate(abr_env_conda, conda_env_abricate)

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
