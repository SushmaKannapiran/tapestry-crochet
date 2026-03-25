#!/usr/bin/env python3
"""
Mycobacterium tuberculosis WGS Analysis Pipeline — Final Production Script

Pipeline:
  1. TB-Profiler  → lineage, drug resistance, QC from nanopore FASTQ
  2. SpoTyping    → in silico spoligotyping
  3. MEGAHIT      → de novo assembly
  4. ABRicate     → virulence factor detection (VFDB)
  5. Report       → structured text report matching client format

Usage:
  python3 TB_Script_Final.py <sample_dir> <sample_name> <threads>

Arguments:
  sample_dir   : Directory containing the FASTQ file
  sample_name  : Sample identifier (filename stem without .fastq/.fastq.gz)
  threads      : Number of threads
"""

import os
import subprocess
import shutil
import sys
import re
import json
import glob
import gzip
import pandas as pd
import urllib.request
import urllib.parse


# =====================================================================
# CONFIGURATION
# =====================================================================

conda_base = "/media/carme/Idaa/conda/anaconda3/etc/profile.d/conda.sh"
conda_env_tbprofiler = "tb_profiler_env"

abr_env_conda = "/media/carme/ganesh/miniconda3/etc/profile.d/conda.sh"
conda_env_abricate = "abricate_env"

SPOTYPING_CMD = (
    "python2 /media/carme/Idaa/Scripts/Tuberculosis_Report/"
    "Jan_2026/SpoTyping-v2.1-commandLine/SpoTyping.py"
)

# TB-Profiler database directory — contains barcode BED and GFF
TBPROFILER_DB_DIR = (
    "/media/carme/Idaa/conda/anaconda3/envs/tb_profiler_env/"
    "share/tbprofiler/"
)

sample_dir = sys.argv[1]
sample = sys.argv[2]
threads = sys.argv[3]

# Detect FASTQ file (support both .fastq and .fastq.gz)
fastq_file = os.path.join(sample_dir, f"{sample}.fastq")
if not os.path.exists(fastq_file):
    fastq_file = os.path.join(sample_dir, f"{sample}.fastq.gz")
    if not os.path.exists(fastq_file):
        print(f"Error: No FASTQ file found for sample '{sample}' in {sample_dir}")
        print(f"  Checked: {sample}.fastq and {sample}.fastq.gz")
        sys.exit(1)

# Output paths
tb_profiler_dir = os.path.join(sample_dir, f"{sample}_tbprofler")
megahit_out_dir = f"{sample_dir}/Megahit/"
abricate_result = f"{sample_dir}/{sample}_abricate_vfdb.tab"
vfdb_annot = (
    "/media/carme/Research_Projects/BugSpeaks_MicrobiomeAnalysis/"
    "Tuberculosis_Firoz/VFDB/VFs_VFDB.csv"
)
spoligotype_output = f"{sample_dir}/{sample}_spoligotype"


# =====================================================================
# H37Rv GENE REFERENCE DATA (NC_000962.3)
# 73 target genes — used as fallback when full GFF is unavailable
# =====================================================================

H37RV_GENES = {
    "PPE35":   {"start": 2167649, "stop": 2170612, "size": 2964, "strand": "-", "locus_tag": "Rv1918c"},
    "Rv0010c": {"start": 13133,   "stop": 13558,   "size": 426,  "strand": "-", "locus_tag": "Rv0010c"},
    "Rv0565c": {"start": 656010,  "stop": 657470,  "size": 1461, "strand": "-", "locus_tag": "Rv0565c"},
    "Rv1129c": {"start": 1253074, "stop": 1254534, "size": 1461, "strand": "-", "locus_tag": "Rv1129c"},
    "Rv1258c": {"start": 1406081, "stop": 1407340, "size": 1260, "strand": "-", "locus_tag": "Rv1258c"},
    "Rv1979c": {"start": 2221719, "stop": 2223164, "size": 1446, "strand": "-", "locus_tag": "Rv1979c"},
    "Rv2477c": {"start": 2782366, "stop": 2784042, "size": 1677, "strand": "-", "locus_tag": "Rv2477c"},
    "Rv2680":  {"start": 2996105, "stop": 2996737, "size": 633,  "strand": "+", "locus_tag": "Rv2680"},
    "Rv2681":  {"start": 2996739, "stop": 2998055, "size": 1317, "strand": "+", "locus_tag": "Rv2681"},
    "Rv2752c": {"start": 3064515, "stop": 3066191, "size": 1677, "strand": "-", "locus_tag": "Rv2752c"},
    "Rv3083":  {"start": 3448504, "stop": 3449991, "size": 1488, "strand": "+", "locus_tag": "Rv3083"},
    "Rv3236c": {"start": 3611959, "stop": 3613116, "size": 1158, "strand": "-", "locus_tag": "Rv3236c"},
    "aftB":    {"start": 4266953, "stop": 4268836, "size": 1884, "strand": "-", "locus_tag": "Rv3805c"},
    "ahpC":    {"start": 2726193, "stop": 2726780, "size": 588,  "strand": "+", "locus_tag": "Rv2428"},
    "ald":     {"start": 3086820, "stop": 3087935, "size": 1116, "strand": "+", "locus_tag": "Rv2780"},
    "alr":     {"start": 3840194, "stop": 3841420, "size": 1227, "strand": "-", "locus_tag": "Rv3423c"},
    "atpE":    {"start": 1461045, "stop": 1461290, "size": 246,  "strand": "+", "locus_tag": "Rv1305"},
    "bacA":    {"start": 2062809, "stop": 2064728, "size": 1920, "strand": "-", "locus_tag": "Rv1819c"},
    "ccsA":    {"start": 619891,  "stop": 620865,  "size": 975,  "strand": "+", "locus_tag": "Rv0529"},
    "clpC1":   {"start": 4038158, "stop": 4040704, "size": 2547, "strand": "-", "locus_tag": "Rv3596c"},
    "ddn":     {"start": 3986844, "stop": 3987299, "size": 456,  "strand": "+", "locus_tag": "Rv3547"},
    "dnaA":    {"start": 1,       "stop": 1524,    "size": 1524, "strand": "+", "locus_tag": "Rv0001"},
    "eis":     {"start": 2714124, "stop": 2715332, "size": 1209, "strand": "-", "locus_tag": "Rv2416c"},
    "embA":    {"start": 4243233, "stop": 4246517, "size": 3285, "strand": "+", "locus_tag": "Rv3794"},
    "embB":    {"start": 4246514, "stop": 4249810, "size": 3297, "strand": "+", "locus_tag": "Rv3795"},
    "embC":    {"start": 4239863, "stop": 4243147, "size": 3285, "strand": "+", "locus_tag": "Rv3793"},
    "embR":    {"start": 1416181, "stop": 1417347, "size": 1167, "strand": "-", "locus_tag": "Rv1267c"},
    "ethA":    {"start": 4326004, "stop": 4327473, "size": 1470, "strand": "-", "locus_tag": "Rv3854c"},
    "ethR":    {"start": 4327549, "stop": 4328199, "size": 651,  "strand": "+", "locus_tag": "Rv3855"},
    "fbiA":    {"start": 3640543, "stop": 3641538, "size": 996,  "strand": "+", "locus_tag": "Rv3261"},
    "fbiB":    {"start": 3641535, "stop": 3642881, "size": 1347, "strand": "+", "locus_tag": "Rv3262"},
    "fbiC":    {"start": 1302931, "stop": 1305501, "size": 2571, "strand": "+", "locus_tag": "Rv1173"},
    "fbiD":    {"start": 3339118, "stop": 3339762, "size": 645,  "strand": "+", "locus_tag": "Rv2983"},
    "fgd1":    {"start": 490783,  "stop": 491793,  "size": 1011, "strand": "+", "locus_tag": "Rv0407"},
    "folC":    {"start": 2746135, "stop": 2747598, "size": 1464, "strand": "-", "locus_tag": "Rv2447c"},
    "gid":     {"start": 4407528, "stop": 4408202, "size": 675,  "strand": "-", "locus_tag": "Rv3919c"},
    "glpK":    {"start": 4138202, "stop": 4139755, "size": 1554, "strand": "-", "locus_tag": "Rv3696c"},
    "gyrA":    {"start": 7302,    "stop": 9818,    "size": 2517, "strand": "+", "locus_tag": "Rv0006"},
    "gyrB":    {"start": 5240,    "stop": 7267,    "size": 2028, "strand": "+", "locus_tag": "Rv0005"},
    "hadA":    {"start": 731930,  "stop": 732406,  "size": 477,  "strand": "+", "locus_tag": "Rv0635"},
    "inhA":    {"start": 1674202, "stop": 1675011, "size": 810,  "strand": "+", "locus_tag": "Rv1484"},
    "kasA":    {"start": 2518115, "stop": 2519365, "size": 1251, "strand": "+", "locus_tag": "Rv2245"},
    "katG":    {"start": 2153889, "stop": 2156111, "size": 2223, "strand": "-", "locus_tag": "Rv1908c"},
    "lpqB":    {"start": 3623159, "stop": 3624910, "size": 1752, "strand": "-", "locus_tag": "Rv3244c"},
    "mmaA3":   {"start": 737268,  "stop": 738149,  "size": 882,  "strand": "-", "locus_tag": "Rv0643c"},
    "mmpL5":   {"start": 775586,  "stop": 778480,  "size": 2895, "strand": "-", "locus_tag": "Rv0676c"},
    "mmpR5":   {"start": 778990,  "stop": 779487,  "size": 498,  "strand": "+", "locus_tag": "Rv0678"},
    "mmpS5":   {"start": 778477,  "stop": 778905,  "size": 429,  "strand": "-", "locus_tag": "Rv0677c"},
    "mshA":    {"start": 575348,  "stop": 576790,  "size": 1443, "strand": "+", "locus_tag": "Rv0486"},
    "mtrA":    {"start": 3626663, "stop": 3627349, "size": 687,  "strand": "-", "locus_tag": "Rv3246c"},
    "mtrB":    {"start": 3624910, "stop": 3626613, "size": 1704, "strand": "-", "locus_tag": "Rv3245c"},
    "ndh":     {"start": 2101651, "stop": 2103042, "size": 1392, "strand": "-", "locus_tag": "Rv1854c"},
    "nusG":    {"start": 734254,  "stop": 734970,  "size": 717,  "strand": "+", "locus_tag": "Rv0639"},
    "panD":    {"start": 4043862, "stop": 4044281, "size": 420,  "strand": "-", "locus_tag": "Rv3601c"},
    "pepQ":    {"start": 2859300, "stop": 2860418, "size": 1119, "strand": "-", "locus_tag": "Rv2535c"},
    "pncA":    {"start": 2288681, "stop": 2289241, "size": 561,  "strand": "-", "locus_tag": "Rv2043c"},
    "ribD":    {"start": 2986839, "stop": 2987615, "size": 777,  "strand": "+", "locus_tag": "Rv2671"},
    "rplC":    {"start": 800809,  "stop": 801462,  "size": 654,  "strand": "+", "locus_tag": "Rv0701"},
    "rpoA":    {"start": 3877464, "stop": 3878507, "size": 1044, "strand": "-", "locus_tag": "Rv3457c"},
    "rpoB":    {"start": 759807,  "stop": 763325,  "size": 3519, "strand": "+", "locus_tag": "Rv0667"},
    "rpoC":    {"start": 763370,  "stop": 767320,  "size": 3951, "strand": "+", "locus_tag": "Rv0668"},
    "rpsA":    {"start": 1833542, "stop": 1834987, "size": 1446, "strand": "+", "locus_tag": "Rv1630"},
    "rpsL":    {"start": 781560,  "stop": 781934,  "size": 375,  "strand": "+", "locus_tag": "Rv0682"},
    "rrl":     {"start": 1473658, "stop": 1476795, "size": 3138, "strand": "+", "locus_tag": "EBG00000313339"},
    "rrs":     {"start": 1471846, "stop": 1473382, "size": 1537, "strand": "+", "locus_tag": "EBG00000313325"},
    "sigE":    {"start": 1364413, "stop": 1365186, "size": 774,  "strand": "+", "locus_tag": "Rv1221"},
    "thyA":    {"start": 3073680, "stop": 3074471, "size": 792,  "strand": "-", "locus_tag": "Rv2764c"},
    "thyX":    {"start": 3067193, "stop": 3067945, "size": 753,  "strand": "-", "locus_tag": "Rv2754c"},
    "tlyA":    {"start": 1917940, "stop": 1918746, "size": 807,  "strand": "+", "locus_tag": "Rv1694"},
    "tsnR":    {"start": 1853606, "stop": 1854388, "size": 783,  "strand": "+", "locus_tag": "Rv1644"},
    "ubiA":    {"start": 4268925, "stop": 4269833, "size": 909,  "strand": "-", "locus_tag": "Rv3806c"},
    "whiB6":   {"start": 4338171, "stop": 4338521, "size": 351,  "strand": "-", "locus_tag": "Rv3862c"},
    "whiB7":   {"start": 3568401, "stop": 3568679, "size": 279,  "strand": "-", "locus_tag": "Rv3197A"},
}

LOCUS_TAG_TO_GENE = {v['locus_tag']: k for k, v in H37RV_GENES.items()}

# Comprehensive gene name → locus tag mapping for H37Rv.
# Covers all 73 target genes PLUS genes commonly found in lineage barcoding.
# Used as fallback when the GFF doesn't provide locus_tag attributes.
GENE_NAME_TO_LOCUS = {v.get('locus_tag', k): v.get('locus_tag', k) for k, v in H37RV_GENES.items()}
GENE_NAME_TO_LOCUS.update({k: v.get('locus_tag', k) for k, v in H37RV_GENES.items()})
GENE_NAME_TO_LOCUS.update({
    "eccB3": "Rv0282", "eccC3": "Rv0284", "eccD3": "Rv0290",
    "eccE3": "Rv0292", "eccA3": "Rv0282", "mycP3": "Rv0291",
    "esxH": "Rv0288", "esxG": "Rv0287", "PE5": "Rv0285",
    "hemL": "Rv0524", "menD": "Rv0555", "rplN": "Rv0714",
    "rplE": "Rv0716", "purD": "Rv0772", "pyrF": "Rv1380",
    "ribA2": "Rv1415", "aspS": "Rv2572c", "nrp": "Rv0101",
    "dnaJ1": "Rv0352", "aftC": "Rv2673", "espD": "Rv3614c",
    "trpG": "Rv0013", "phoR": "Rv0758", "emrB": "Rv0783c",
    "lpqU": "Rv1017c", "narJ": "Rv1162", "glgB": "Rv1326c",
    "malQ": "Rv1781c", "secA2": "Rv1821", "obg": "Rv2440c",
    "ddlA": "Rv2981c", "espA": "Rv3616c", "espC": "Rv3615c",
    "espB": "Rv3881c", "espK": "Rv3879c",
    "eccA1": "Rv3868", "eccB1": "Rv3869", "eccCa1": "Rv3870",
    "eccCb1": "Rv3871", "eccD1": "Rv3877", "eccE1": "Rv3882c",
    "mycP1": "Rv3883c", "esxA": "Rv3875", "esxB": "Rv3874",
    "PE35": "Rv3872",
    "eccA5": "Rv1797", "eccB5": "Rv1795", "eccCa5": "Rv1796",
    "eccCb5": "Rv1798", "eccD5": "Rv1795", "eccE5": "Rv1800",
    "mycP5": "Rv1796", "esxN": "Rv1793", "esxM": "Rv1792",
    "PPE41": "Rv2430c",
    "ideR": "Rv2711", "icl": "Rv0467", "mgtC": "Rv1811",
    "hspX": "Rv2031c", "phoP": "Rv0757", "relA": "Rv2583c",
    "lipF": "Rv3487c", "hbhA": "Rv0475", "erp": "Rv3810",
    "fbpA": "Rv3804c", "fbpB": "Rv1886c", "fbpC": "Rv0129c",
    "mbtA": "Rv2384", "mbtB": "Rv2383c", "mbtC": "Rv2382c",
    "mbtD": "Rv2381c", "mbtE": "Rv2380c", "mbtF": "Rv2379c",
    "mbtG": "Rv2378c", "mbtH": "Rv2377c", "mbtI": "Rv2386c",
    "mbtJ": "Rv2385", "mbtK": "Rv1347c", "mbtL": "Rv1344",
    "mbtM": "Rv1345", "mbtN": "Rv1346",
    "irtA": "Rv1348", "irtB": "Rv1349",
    "espG3": "Rv0289", "PPE4": "Rv0286",
    "pknB": "Rv0014c", "pknA": "Rv0015c", "pbpA": "Rv0016c",
    "rodA": "Rv0017c", "ppp": "Rv0018c",
    "murD": "Rv2155c", "murG": "Rv2153c", "murC": "Rv2152c",
    "murF": "Rv2157c", "murE": "Rv2158c",
    "leuS": "Rv0041", "metG": "Rv2091c", "ileS": "Rv1536",
    "pheS": "Rv1649", "thrS": "Rv1350",
    "aroA": "Rv3227", "aroB": "Rv2538c", "aroC": "Rv2537c",
    "aroD": "Rv2539c", "aroE": "Rv2552c", "aroK": "Rv2539c",
    "proC": "Rv0500", "proA": "Rv1083", "proB": "Rv1084",
    "lysA": "Rv1293", "argB": "Rv1654", "argC": "Rv1655",
    "hisA": "Rv1603", "hisB": "Rv1601", "hisD": "Rv1599",
    "trpA": "Rv1613", "trpB": "Rv1612", "trpC": "Rv1611",
    "trpD": "Rv1609", "trpE": "Rv1610",
    "rpoA": "Rv3457c", "rpoB": "Rv0667", "rpoC": "Rv0668",
    "rpsA": "Rv1630", "rpsL": "Rv0682", "rpsG": "Rv0683",
    "rplC": "Rv0701", "rplD": "Rv0702", "rplB": "Rv0703",
    "rplW": "Rv0707", "rplV": "Rv0708", "rplP": "Rv0710",
    "rpmC": "Rv0711", "rpsC": "Rv0712", "rplR": "Rv0720",
    "rplF": "Rv0719",
    "gyrA": "Rv0006", "gyrB": "Rv0005", "dnaA": "Rv0001",
    "dnaN": "Rv0002", "recF": "Rv0003", "dnaB": "Rv0058",
    "dnaG": "Rv2343c",
    "accD3": "Rv0904c", "accD5": "Rv3280",
    "infB": "Rv2839c", "infC": "Rv1641", "fusA1": "Rv0684",
    "tuf": "Rv0685", "tsf": "Rv3652",
    "groEL2": "Rv0440", "groES": "Rv3418c", "groEL1": "Rv3417c",
    "dnaK": "Rv0350", "dnaJ2": "Rv2373c",
    "clpB": "Rv0384c", "clpP1": "Rv2461c", "clpP2": "Rv2460c",
    "tig": "Rv2462c",
    "fadD26": "Rv2930", "fadD28": "Rv2941", "fadD29": "Rv2950c",
    "pks12": "Rv2048c", "pks13": "Rv3800c",
    "mas": "Rv2940c", "ppsA": "Rv2931", "ppsB": "Rv2932",
    "ppsC": "Rv2933", "ppsD": "Rv2934", "ppsE": "Rv2935",
    "mmpl8": "Rv3823c", "lppX": "Rv2945c",
    "pknG": "Rv0410c", "pknH": "Rv1266c", "pknD": "Rv0931c",
    "pknE": "Rv1743", "pknF": "Rv1746",
    "devR": "Rv3133c", "devS": "Rv3132c", "dosT": "Rv2027c",
    "sigA": "Rv2703", "sigB": "Rv2710", "sigC": "Rv2069",
    "sigD": "Rv3414c", "sigH": "Rv3223c",
    "recA": "Rv2737c", "lexA": "Rv2720",
    "nuoA": "Rv3145", "nuoB": "Rv3146", "nuoC": "Rv3147",
    "nuoD": "Rv3148", "nuoE": "Rv3149", "nuoF": "Rv3150",
    "nuoG": "Rv3151", "nuoH": "Rv3152", "nuoI": "Rv3153",
    "nuoJ": "Rv3154", "nuoK": "Rv3155", "nuoL": "Rv3156",
    "nuoM": "Rv3157", "nuoN": "Rv3158",
    "atpA": "Rv1308", "atpB": "Rv1304", "atpC": "Rv1311",
    "atpD": "Rv1310", "atpE": "Rv1305", "atpF": "Rv1306",
    "atpG": "Rv1309", "atpH": "Rv1307",
    "secA1": "Rv3240c", "secA2": "Rv1821", "secY": "Rv0732",
    "tatA": "Rv2094c", "tatB": "Rv1224",
    "glnA1": "Rv2220", "glnA2": "Rv2222c",
    "pknL": "Rv2176", "pstS1": "Rv0934",
    "Rv0101": "Rv0101", "Rv0013": "Rv0013",
})

FAMILY_ABBREV = {
    'Indo-Oceanic': 'EAI',
    'East-Asian': 'Beijing',
    'East-African-Indian': 'CAS',
    'West-African 1': 'Cameroon',
    'West-African 2': 'West African 2',
    'Ethiopia': 'Ethiopia',
    'La1': 'La1',
    'M. bovis': 'Bovis',
    'M. caprae': 'Caprae',
}

GENE_DRUGS = {
    'rpoB': 'rifampicin', 'rpoC': 'rifampicin',
    'katG': 'isoniazid', 'inhA': 'isoniazid, ethionamide', 'ahpC': 'isoniazid',
    'kasA': 'isoniazid', 'ndh': 'isoniazid',
    'embB': 'ethambutol', 'embA': 'ethambutol', 'embC': 'ethambutol',
    'embR': 'ethambutol', 'ubiA': 'ethambutol',
    'pncA': 'pyrazinamide', 'panD': 'pyrazinamide',
    'rpsA': 'pyrazinamide', 'clpC1': 'pyrazinamide',
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
    'fbiA': 'delamanid, pretomanid', 'fbiB': 'delamanid, pretomanid',
    'fbiC': 'delamanid, pretomanid', 'fbiD': 'delamanid, pretomanid',
    'fgd1': 'delamanid, clofazimine, pretomanid', 'ddn': 'delamanid, pretomanid',
    'folC': 'para-aminosalicylic_acid', 'ribD': 'para-aminosalicylic_acid',
    'thyA': 'para-aminosalicylic_acid', 'thyX': 'para-aminosalicylic_acid',
    'alr': 'cycloserine',
    'tlyA': 'capreomycin',
    'sigE': 'multiple', 'whiB6': 'multiple', 'whiB7': 'multiple',
}

DRUG_DISPLAY = {
    'rifampicin': 'Rifampicin (RMP)',
    'isoniazid': 'Isoniazid (INH)',
    'ethambutol': 'Ethambutol (EMB)',
    'pyrazinamide': 'Pyrazinamide (PZA)',
    'moxifloxacin': 'Moxifloxacin (MFL)',
    'levofloxacin': 'Levofloxacin (LFX)',
    'ofloxacin': 'Ofloxacin (OFX)',
    'fluoroquinolone': 'Fluoroquinolone (FQ)',
    'bedaquiline': 'Bedaquiline (BED)',
    'delamanid': 'Delamanid (DLM)',
    'pretomanid': 'Pretomanid (PMD)',
    'linezolid': 'Linezolid (LZD)',
    'streptomycin': 'Streptomycin (SM)',
    'amikacin': 'Amikacin (AMK)',
    'kanamycin': 'Kanamycin (KAN)',
    'capreomycin': 'Capreomycin (CAP)',
    'clofazimine': 'Clofazimine (CLO)',
    'ethionamide': 'Ethionamide (ETH)',
    'para-aminosalicylic_acid': 'Para-aminosalicylic acid (PAS)',
    'cycloserine': 'Cycloserine (CS)',
}

DRUG_ABBREV = {
    'rifampicin': 'RMP', 'isoniazid': 'INH', 'ethambutol': 'EMB',
    'pyrazinamide': 'PZA', 'moxifloxacin': 'MFL', 'levofloxacin': 'LFX',
    'ofloxacin': 'OFX', 'fluoroquinolone': 'FQ', 'bedaquiline': 'BED',
    'delamanid': 'DLM', 'pretomanid': 'PMD', 'linezolid': 'LZD',
    'streptomycin': 'SM', 'amikacin': 'AMK', 'kanamycin': 'KAN',
    'capreomycin': 'CAP', 'clofazimine': 'CLO', 'ethionamide': 'ETH',
    'para-aminosalicylic_acid': 'PAS', 'cycloserine': 'CS',
}


# =====================================================================
# GENOME ANNOTATION LOADERS
# =====================================================================

def load_barcode_data():
    """
    Parse TB-Profiler's barcode BED file to get ref/alt alleles
    for every lineage-defining SNP.

    BED format: chrom  start(0-based)  end  lineage  ref  alt
    Returns dict mapping 1-based position → {ref, alt, lineage}

    Uses multiple strategies to locate the barcode file.
    """
    conda_env_root = ""
    if '/share/' in TBPROFILER_DB_DIR:
        conda_env_root = TBPROFILER_DB_DIR.split('/share/')[0]
    elif '/envs/' in TBPROFILER_DB_DIR:
        conda_env_root = TBPROFILER_DB_DIR.split('/envs/')[0]
        conda_env_root = os.path.join(conda_env_root, 'envs', 'tb_profiler_env')

    possible_files = [
        os.path.join(TBPROFILER_DB_DIR, "tbdb.barcode.bed"),
        os.path.join(TBPROFILER_DB_DIR, "barcode.bed"),
    ]

    if conda_env_root:
        possible_files.extend([
            os.path.join(conda_env_root, "share", "tbprofiler", "tbdb.barcode.bed"),
            os.path.join(conda_env_root, "share", "tbprofiler", "barcode.bed"),
            os.path.join(conda_env_root, "share", "tb-profiler", "tbdb.barcode.bed"),
        ])
        for pylib in glob.glob(os.path.join(conda_env_root, "lib", "python*", "site-packages")):
            possible_files.append(os.path.join(pylib, "tbprofiler", "data", "tbdb.barcode.bed"))
            possible_files.append(os.path.join(pylib, "tb_profiler", "data", "tbdb.barcode.bed"))

    home = os.path.expanduser("~")
    possible_files.extend([
        os.path.join(home, ".tbprofiler", "tbdb.barcode.bed"),
        os.path.join(home, ".tb-profiler", "tbdb.barcode.bed"),
    ])

    bed_file = None
    for f in possible_files:
        if os.path.exists(f):
            bed_file = f
            break

    if not bed_file and conda_env_root:
        try:
            result = subprocess.run(
                ['find', conda_env_root, '-name', 'tbdb.barcode.bed', '-type', 'f'],
                capture_output=True, text=True, timeout=15
            )
            if result.returncode == 0 and result.stdout.strip():
                found = result.stdout.strip().split('\n')[0]
                if os.path.exists(found):
                    bed_file = found
                    print(f"  Found barcode BED via filesystem search: {bed_file}")
        except Exception:
            pass

    if not bed_file:
        print("  WARNING: Barcode BED not found in any known location")
        print(f"  Tried {len(possible_files)} paths including: {possible_files[0]}")
        print("  WT/Call will be resolved from VCF fallback if available")
        return {}

    lookup = {}
    try:
        with open(bed_file) as fh:
            for line in fh:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 6:
                    pos_1based = int(parts[1]) + 1
                    lookup[pos_1based] = {
                        'ref': parts[4],
                        'alt': parts[5],
                        'lineage': parts[3],
                    }
                elif len(parts) == 5:
                    pos_1based = int(parts[1]) + 1
                    lookup[pos_1based] = {
                        'ref': parts[3],
                        'alt': parts[4],
                        'lineage': '',
                    }
        print(f"  Loaded {len(lookup)} lineage barcode SNPs from {bed_file}")
    except Exception as e:
        print(f"  Warning: Could not parse barcode BED: {e}")

    return lookup


def load_vcf_alleles():
    """
    Parse the VCF produced by TB-Profiler to get ref/alt alleles
    for every variant position. This serves as a reliable fallback
    for lineage WT/Call when the barcode BED file isn't found.

    TB-Profiler always generates a VCF as part of its pipeline.
    Returns dict mapping 1-based position → {ref, alt}
    """
    vcf_dir = os.path.join(tb_profiler_dir, "vcf")
    results_dir = os.path.join(tb_profiler_dir, "results")
    bam_dir = os.path.join(tb_profiler_dir, "bam")

    possible_vcfs = [
        os.path.join(vcf_dir, f"{sample}.targets.csq.vcf.gz"),
        os.path.join(vcf_dir, f"{sample}.vcf.gz"),
        os.path.join(vcf_dir, f"{sample}.targets.vcf.gz"),
        os.path.join(vcf_dir, f"{sample}.vcf"),
        os.path.join(results_dir, f"{sample}.vcf.gz"),
        os.path.join(results_dir, f"{sample}.vcf"),
    ]

    if os.path.isdir(vcf_dir):
        for f in sorted(glob.glob(os.path.join(vcf_dir, f"{sample}*.vcf*"))):
            if f not in possible_vcfs:
                possible_vcfs.append(f)

    vcf_file = None
    for f in possible_vcfs:
        if os.path.exists(f):
            vcf_file = f
            break

    if not vcf_file:
        print("  Note: No VCF file found for allele fallback")
        return {}

    lookup = {}
    try:
        opener = gzip.open if vcf_file.endswith('.gz') else open
        with opener(vcf_file, 'rt') as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 5:
                    pos = int(parts[1])
                    ref = parts[3]
                    alt = parts[4].split(',')[0]
                    if re.match(r'^[ACGTacgt]+$', ref) and re.match(r'^[ACGTacgt]+$', alt):
                        lookup[pos] = {'ref': ref, 'alt': alt}
        print(f"  Loaded {len(lookup)} variant alleles from VCF: {os.path.basename(vcf_file)}")
    except Exception as e:
        print(f"  Warning: Could not parse VCF: {e}")

    return lookup


def load_full_gene_annotation():
    """
    Parse H37Rv GFF3 to build a position→gene lookup for ALL ~4000 genes,
    not just the 73 target genes.

    Extracts BOTH Name (gene name) and locus_tag (Rv number) independently.

    Falls back to H37RV_GENES if GFF is not found.
    Returns sorted list of (start, stop, gene_name, locus_tag) 4-tuples.
    """
    possible_gff = [
        os.path.join(TBPROFILER_DB_DIR, "tbdb.gff"),
        os.path.join(TBPROFILER_DB_DIR, "genome.gff"),
        os.path.join(TBPROFILER_DB_DIR, "tbdb.gff3"),
    ]

    gff_file = None
    for f in possible_gff:
        if os.path.exists(f):
            gff_file = f
            break

    genes = []

    if gff_file:
        try:
            with open(gff_file) as fh:
                for line in fh:
                    if line.startswith('#') or not line.strip():
                        continue
                    parts = line.strip().split('\t')
                    if len(parts) < 9:
                        continue
                    if parts[2] != 'gene':
                        continue
                    start = int(parts[3])
                    stop = int(parts[4])
                    attrs = parts[8]

                    gene_name = None
                    locus_tag = None

                    m = re.search(r'Name=([^;]+)', attrs)
                    if m:
                        gene_name = m.group(1)

                    m = re.search(r'locus_tag=([^;]+)', attrs)
                    if m:
                        locus_tag = m.group(1)
                    if not locus_tag:
                        m = re.search(r'ID=gene-(Rv\w+)', attrs)
                        if m:
                            locus_tag = m.group(1)
                    if not locus_tag and gene_name:
                        locus_tag = GENE_NAME_TO_LOCUS.get(gene_name)

                    display = gene_name or locus_tag
                    if display:
                        resolved_tag = locus_tag or GENE_NAME_TO_LOCUS.get(display, display)
                        genes.append((start, stop, display, resolved_tag))
            genes.sort()
            print(f"  Loaded {len(genes)} genes from GFF annotation (with locus tags)")
        except Exception as e:
            print(f"  Warning: Could not parse GFF: {e}")
            genes = []

    if not genes:
        for name, info in H37RV_GENES.items():
            genes.append((info['start'], info['stop'], name, info.get('locus_tag', name)))
        genes.sort()
        print(f"  Using embedded {len(genes)} target genes (GFF not found)")

    return genes


def find_gene_for_position(pos, gene_list):
    """
    Check if a genomic position falls within any gene.
    Handles both 3-tuples (start, stop, name) and 4-tuples (start, stop, name, locus_tag).
    Returns (gene_name, locus_tag) tuple. Both are None if position not in any gene.
    """
    if not isinstance(pos, (int, float)):
        return None, None
    pos = int(pos)
    for item in gene_list:
        start, stop, name = item[0], item[1], item[2]
        ltag = item[3] if len(item) > 3 else name
        if start > pos:
            break
        if start <= pos <= stop:
            return name, ltag
    return None, None


# =====================================================================
# OFFLINE SIT LOOKUP TABLE (from SpolDB4 / published literature)
# Verified octal-code → SIT mappings for the most common spoligotypes
# Sources: Brudey et al. 2006 (BMC Microbiol), Sharma et al. 2018
#          (Sci Rep), Zenteno-Cuevas et al. 2021 (Sci Rep)
# =====================================================================

OCTAL_TO_SIT = {
    # Beijing family
    '000000000003771': ('SIT1', 'Beijing'),
    '000000000003731': ('SIT190', 'Beijing'),
    '000000000003631': ('SIT1168', 'Beijing'),
    '000000000000371': ('SIT250', 'Beijing'),
    '000000000002771': ('SIT621', 'Beijing'),
    # CAS (Central Asian) family
    '703777740003771': ('SIT26', 'CAS1_DELHI'),
    '703777740003171': ('SIT25', 'CAS1_DELHI'),
    '703777740003571': ('SIT289', 'CAS1_DELHI'),
    '703777740003371': ('SIT428', 'CAS1_DELHI'),
    '703777740003731': ('SIT429', 'CAS1_DELHI'),
    '703777740003671': ('SIT1401', 'CAS1_DELHI'),
    '703777740003011': ('SIT2147', 'CAS1_DELHI'),
    '703757740003771': ('SIT794', 'CAS1_DELHI'),
    '503777740003771': ('SIT754', 'CAS1_DELHI'),
    '703777740002771': ('SIT1091', 'CAS1_DELHI'),
    '702777740003771': ('SIT1092', 'CAS1_DELHI'),
    '703637740003771': ('SIT1327', 'CAS1_DELHI'),
    '703737740003771': ('SIT1343', 'CAS1_DELHI'),
    '703377740003771': ('SIT1942', 'CAS1_DELHI'),
    '703601740003771': ('SIT2364', 'CAS1_DELHI'),
    '700377740003771': ('SIT288', 'CAS2'),
    '703777740000771': ('SIT357', 'CAS'),
    '703777740000371': ('SIT486', 'CAS'),
    '703777740000171': ('SIT1789', 'CAS'),
    '703760000000331': ('SIT1120', 'CAS'),
    '703000000000371': ('SIT2419', 'CAS'),
    # EAI (East-African-Indian) family
    '477777777413071': ('SIT11', 'EAI3_IND'),
    '677777477413771': ('SIT19', 'EAI2_MANILA'),
    '777777777413731': ('SIT48', 'EAI2_SOM'),
    '777777777413771': ('SIT236', 'EAI5'),
    '777777757413771': ('SIT591', 'EAI6_BGD1'),
    '077737777413771': ('SIT1628', 'EAI5'),
    '777776757413771': ('SIT1970', 'EAI6_BGD'),
    # Haarlem family
    '777777774020771': ('SIT47', 'H1'),
    '777777704020771': ('SIT283', 'H1'),
    '777777774020611': ('SIT2642', 'H1'),
    '000000004020771': ('SIT2', 'H2'),
    '777777760020611': ('SIT948', 'H3'),
    '577777777420771': ('SIT127', 'H4'),
    # LAM family
    '777777607760771': ('SIT42', 'LAM9'),
    '677777607760771': ('SIT20', 'LAM1'),
    '637777607760771': ('SIT578', 'LAM1'),
    '677737607760771': ('SIT17', 'LAM2'),
    '777777607760731': ('SIT60', 'LAM4'),
    '577737607760771': ('SIT3019', 'LAM5'),
    '777577607760771': ('SIT1535', 'LAM9'),
    # T family
    '777777777760771': ('SIT53', 'T1'),
    '777777777760731': ('SIT52', 'T2'),
    '777777777760700': ('SIT51', 'T1'),
    '777777777760600': ('SIT243', 'T1'),
    '777777777740771': ('SIT172', 'T1'),
    '777777775760771': ('SIT122', 'T1'),
    '777777677760771': ('SIT291', 'T1'),
    '777777777560771': ('SIT462', 'T1'),
    '757777777760771': ('SIT154', 'T1'),
    '700077777760771': ('SIT344', 'T1'),
    '000000177760771': ('SIT258', 'T'),
    '000000007760731': ('SIT125', 'T2'),
    '777737777760771': ('SIT37', 'T3'),
    # S family
    '776377777760771': ('SIT34', 'S'),
    # X family
    '777776777760771': ('SIT119', 'X1'),
    '777776777760671': ('SIT1394', 'X1'),
    '777776777760601': ('SIT137', 'X2'),
    '700076777760771': ('SIT92', 'X3'),
    # MANU family
    '777777777763771': ('SIT54', 'MANU2'),
    '777777777773771': ('SIT100', 'MANU1'),
    # URAL family
    '703777747770371': ('SIT27', 'URAL'),
    '777777777747771': ('SIT464', 'URAL'),
}


# =====================================================================
# SITVIT SB NUMBER PARSER
# =====================================================================

def parse_sitvit_result(binary_code='', octal_code=''):
    """
    Get SIT number for a spoligotype pattern using four strategies:
      1. Parse SITVIT XLS file if SpoTyping generated one
      2. Offline lookup table (covers ~70 most common global patterns)
      3. Query SITVIT online directly using the binary code
      4. Return N/A if all fail
    """
    # --- Strategy 1: Find and parse local SITVIT XLS file ---
    sit = _parse_sitvit_xls()
    if sit != 'N/A':
        print(f"  SIT number from SITVIT XLS file: {sit}")
        return sit

    # --- Strategy 2: Offline lookup by octal code ---
    if octal_code and octal_code != 'N/A':
        entry = OCTAL_TO_SIT.get(octal_code)
        if entry:
            sit_num, family = entry
            print(f"  SIT number from offline lookup: {sit_num} ({family})")
            return sit_num

    # --- Strategy 3: Query SITVIT online directly ---
    if binary_code and binary_code != 'N/A' and len(binary_code) == 43:
        sit = _query_sitvit_online(binary_code)
        if sit != 'N/A':
            return sit

    if octal_code and octal_code != 'N/A':
        print(f"  Note: Octal code {octal_code} not in offline lookup table "
              f"and SITVIT unreachable. SB Number = N/A")
    return 'N/A'


def _parse_sitvit_xls():
    """Search for SITVIT XLS files across multiple directories."""
    search_dirs = [
        spoligotype_output,
        sample_dir,
        os.getcwd(),
        tb_profiler_dir,
        os.path.join(sample_dir, f"{sample}_spoligotype"),
    ]

    sitvit_file = None
    for d in search_dirs:
        if not d or not os.path.isdir(d):
            continue
        matches = glob.glob(os.path.join(d, "SITVIT_ONLINE.*.xls"))
        matches += glob.glob(os.path.join(d, "*", "SITVIT_ONLINE.*.xls"))
        if matches:
            sitvit_file = matches[0]
            break

    if not sitvit_file:
        return 'N/A'

    try:
        with open(sitvit_file, 'r', errors='ignore') as f:
            content = f.read()

        if not content.strip() or 'No match' in content or len(content) < 50:
            return 'N/A'

        sit_match = re.search(r'<td[^>]*>\s*(\d+)\s*</td>', content)
        if sit_match:
            return f"SIT{sit_match.group(1)}"

        sit_match = re.search(r'SIT[:\s#]*(\d+)', content, re.IGNORECASE)
        if sit_match:
            return f"SIT{sit_match.group(1)}"

        try:
            tables = pd.read_html(sitvit_file)
            if tables and len(tables) > 0:
                df = tables[0]
                for col in df.columns:
                    col_str = str(col).upper()
                    if 'SIT' in col_str or 'SB' in col_str or 'TYPE' in col_str:
                        val = df[col].iloc[0]
                        if pd.notna(val):
                            return f"SIT{val}"
                for col in df.columns:
                    val = df[col].iloc[0]
                    if pd.notna(val) and str(val).isdigit():
                        return f"SIT{val}"
        except Exception:
            pass

        return 'N/A'
    except Exception:
        return 'N/A'


def _query_sitvit_online(binary_code):
    """
    Query SITVIT2 web service directly with the 43-char binary spoligotype.
    This is the same query SpoTyping makes internally.
    """
    url = ("http://www.pasteur-guadeloupe.fr:8081/SITVIT_ONLINE/"
           "tools/query_spoligotype.php")

    try:
        data = urllib.parse.urlencode({
            "binary": binary_code,
            "submit": "Query"
        }).encode('utf-8')

        req = urllib.request.Request(url, data=data)
        req.add_header('User-Agent', 'TB-Pipeline/1.0')

        response = urllib.request.urlopen(req, timeout=15)
        html = response.read().decode('utf-8', errors='ignore')

        if not html or 'No match' in html:
            return 'N/A'

        sit_match = re.search(r'<td[^>]*>\s*(\d+)\s*</td>', html)
        if sit_match:
            print(f"  SITVIT online query: SIT{sit_match.group(1)}")
            return f"SIT{sit_match.group(1)}"

        sit_match = re.search(r'SIT[:\s#]*(\d+)', html, re.IGNORECASE)
        if sit_match:
            print(f"  SITVIT online query: SIT{sit_match.group(1)}")
            return f"SIT{sit_match.group(1)}"

        try:
            debug_file = os.path.join(sample_dir, "SITVIT_debug.html")
            with open(debug_file, 'w') as f:
                f.write(html)
            print(f"  SITVIT response saved to {debug_file} for debugging")
        except Exception:
            pass

        return 'N/A'

    except urllib.error.URLError as e:
        print(f"  SITVIT online query failed (network): {e}")
        return 'N/A'
    except Exception as e:
        print(f"  SITVIT online query failed: {e}")
        return 'N/A'


# =====================================================================
# PIPELINE STEPS
# =====================================================================

def run_tb_profiler(conda_path, env_name):
    cmd = (
        f"source {conda_path} && "
        f"conda activate {env_name} && "
        f"tb-profiler profile --read1 {fastq_file} --platform nanopore --mapper minimap2 "
        f"--caller freebayes --prefix {sample} --dir {tb_profiler_dir} --call_whole_genome "
        f"--depth 5 --af 0.10 --suspect --txt --csv --json"
    )
    subprocess.run(cmd, shell=True, executable="/bin/bash", check=True)


def run_spoltyping(conda_path, env_name):
    cmd = (
        f"{SPOTYPING_CMD} {fastq_file} "
        f"-O {spoligotype_output} -o {sample}"
    )
    try:
        subprocess.run(cmd, shell=True, executable="/bin/bash", check=True)
        print("SpoTyping completed successfully")
    except subprocess.CalledProcessError:
        print("Warning: SpoTyping failed. Continuing without spoligotype.")


def parse_tb_json():
    json_file = f"{tb_profiler_dir}/results/{sample}.results.json"
    if not os.path.exists(json_file):
        print(f"Warning: JSON file not found at {json_file}")
        return None
    with open(json_file, 'r') as f:
        return json.load(f)


def parse_spoligotype():
    fastq_base = os.path.splitext(os.path.basename(fastq_file))[0]
    if fastq_base.endswith('.fastq'):
        fastq_base = fastq_base[:-6]

    possible_files = [
        f"{spoligotype_output}/{sample}",
        f"{spoligotype_output}/{fastq_base}",
        f"{spoligotype_output}.txt",
        f"{sample_dir}/{sample}_spoligotype.txt",
        f"{sample_dir}/{fastq_base}",
        fastq_base,
    ]

    spol_file = None
    for f in possible_files:
        if os.path.exists(f) and os.path.isfile(f):
            spol_file = f
            break

    if not spol_file:
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
                    result['binary'] = parts[1] if len(parts) > 1 else 'N/A'
                    result['octal'] = parts[2] if len(parts) > 2 else 'N/A'
                    break

        # Try to get SB/SIT number (offline lookup → XLS → online query)
        result['SB_number'] = parse_sitvit_result(
            binary_code=result.get('binary', ''),
            octal_code=result.get('octal', '')
        )

        return result
    except Exception as e:
        print(f"Warning: Could not parse spoligotype: {e}")
        return {'octal': 'N/A', 'binary': 'N/A', 'SB_number': 'N/A'}


def run_megahit():
    subprocess.run(
        f"megahit -t {threads} -r {fastq_file} -o {megahit_out_dir}",
        shell=True, check=True
    )
    shutil.move(f"{megahit_out_dir}/final.contigs.fa", f"{sample}.fa")
    shutil.rmtree(megahit_out_dir)


def run_abricate(conda_path, env_name):
    cmd = (
        f"source {conda_path} && "
        f"conda activate {env_name} && "
        f"abricate --db vfdb {sample}.fa > {abricate_result}"
    )
    subprocess.run(cmd, shell=True, executable="/bin/bash", check=True)


def extract_and_merge_vfdb_results():
    if not os.path.exists(abricate_result):
        print("Warning: ABRicate output file not found.")
        return pd.DataFrame()

    abricate_df = pd.read_csv(abricate_result, sep="\t")
    if abricate_df.empty:
        print("Warning: No virulence factors found in abricate output.")
        return pd.DataFrame()

    vfdb_df = pd.read_csv(vfdb_annot)
    abricate_df['VFID'] = abricate_df['PRODUCT'].apply(
        lambda x: re.search(r'\(VF[0-9]+\)', str(x)).group(0)[1:-1]
        if isinstance(x, str) and re.search(r'\(VF[0-9]+\)', x) else None
    )
    merged_result = pd.merge(
        abricate_df, vfdb_df, how='left', left_on='VFID', right_on='VFID'
    ).fillna("NA")
    return merged_result


# =====================================================================
# HELPER FUNCTIONS
# =====================================================================

def build_section_header(title):
    return "\n" + "=" * 60 + "\n" + title + "\n" + "=" * 60 + "\n\n"


def build_variant_position_lookup(data):
    """Maps chromosomal positions to variant details from dr + other variants."""
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


def get_family_display(lineage_list):
    """
    Combines parent lineage family with sub-lineage family.
    e.g., parent=Euro-American + deepest=LAM → "Euro-American (LAM)"
    e.g., parent=East-Asian, no sub-family → "East-Asian (Beijing)"
    """
    if not lineage_list:
        return 'N/A'

    parent_family = lineage_list[0].get('family', 'N/A')
    deepest_family = lineage_list[-1].get('family', 'N/A')

    if (len(lineage_list) > 1
            and parent_family != deepest_family
            and deepest_family not in ('N/A', '', None)):
        return f"{parent_family} ({deepest_family})"

    abbrev = FAMILY_ABBREV.get(parent_family, '')
    if abbrev and abbrev != parent_family:
        return f"{parent_family} ({abbrev})"

    return parent_family


def get_variant_type(var):
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
    nc = var.get('nucleotide_change', '')
    if not nc:
        return ''
    m = re.search(r'c\.(-?\d+)', nc)
    return m.group(1) if m else ''


def get_strand_info(var):
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
    return GENE_DRUGS.get(gene_name, '-')


def drug_display_name(drug):
    """Return 'Isoniazid (INH)' style display name."""
    return DRUG_DISPLAY.get(drug, drug.replace('_', ' ').title())


def get_variant_conclusion(var):
    """Derive Conclusion from annotation confidence for resistance tables."""
    confidence_rank = {
        'not assoc w r': 0, 'not assoc w r - interim': 0,
        'uncertain significance': 1, 'indeterminate': 1,
        'assoc w r - interim': 2, 'assoc w r': 3,
    }
    best_rank = -1
    for annot in var.get('annotation', []):
        conf = annot.get('confidence', '')
        rank = confidence_rank.get(conf.lower(), 1)
        if rank > best_rank:
            best_rank = rank
    if best_rank >= 2:
        return 'Resistant'
    elif best_rank == 1:
        return 'Unknown'
    elif best_rank == 0:
        return 'Not associated'
    return ''


def get_variant_notes(var):
    """Derive Notes field — mark atypical variants."""
    for annot in var.get('annotation', []):
        conf = annot.get('confidence', '').lower()
        if conf in ('uncertain significance', 'indeterminate',
                     'not assoc w r', 'not assoc w r - interim'):
            return 'atypical'
    return ''


def categorize_drugs(data):
    """
    Separates drugs into Resistant / Susceptible / Unknown based on
    the highest-ranked confidence annotation across all variants.
    """
    all_drugs = {
        'rifampicin', 'isoniazid', 'ethambutol', 'pyrazinamide',
        'moxifloxacin', 'levofloxacin', 'ofloxacin', 'fluoroquinolone',
        'bedaquiline', 'delamanid', 'pretomanid', 'linezolid',
        'streptomycin', 'amikacin', 'kanamycin', 'capreomycin',
        'clofazimine', 'ethionamide',
        'para-aminosalicylic_acid', 'cycloserine',
    }

    confidence_rank = {
        'not assoc w r': 0, 'not assoc w r - interim': 0,
        'uncertain significance': 1, 'indeterminate': 1,
        'assoc w r - interim': 2, 'assoc w r': 3,
    }

    drug_best = {}
    for var in data.get('dr_variants', []) + data.get('other_variants', []):
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


# =====================================================================
# SECTION BUILDERS — Matching client report format
# =====================================================================

def build_summary_section(json_data, spoligotype):
    """Title page + summary fields + drug categorization — client format."""
    sample_name = json_data.get('id', sample)
    s = f"{sample_name}\n"
    s += "MTBC functional genotyping report\n\n"
    s += ("The contents of this report are for research purposes only "
          "and not intended for clinical decision making. See below for "
          "full disclaimer.\n\n")

    timestamp = str(json_data.get('timestamp', 'N/A'))
    date_str = timestamp.split('T')[0] if 'T' in timestamp else timestamp.split(' ')[0]
    s += f"Date :\t{date_str}\n"
    s += f"Name :\t{sample_name}\n\n"

    lineage_list = json_data.get('lineage', [])

    lineage_numbers = []
    for entry in lineage_list:
        lin = entry.get('lineage', '')
        num = lin.replace('lineage', '') if lin.startswith('lineage') else lin
        lineage_numbers.append(num)
    lineage_number_str = ' / '.join(lineage_numbers) if lineage_numbers else 'n/a'

    family = lineage_list[-1].get('family', 'n/a') if lineage_list else 'n/a'

    categories = categorize_drugs(json_data)

    unknown_abbrevs = sorted(
        DRUG_ABBREV.get(drug, drug.upper()[:3])
        for drug, _ in categories['unknown']
    )
    unknown_str = ' / '.join(unknown_abbrevs) if unknown_abbrevs else 'n/a'

    drtype = json_data.get('drtype', '')
    drtype_upper = drtype.upper()
    is_mdr = 'MDR' in drtype_upper or 'XDR' in drtype_upper
    is_xdr = 'XDR' in drtype_upper

    s += f"Core (v2) ST :\tn/a\n"
    s += f"Species_determination :\tn/a\n"
    s += f"Lineage_number :\t{lineage_number_str}\n"
    s += f"Lineage_name :\t{family}\n"
    s += f"Spoligotype_octal_code :\t{spoligotype.get('octal', 'n/a')}\n"
    s += f"Spoligotype_SB_number :\t{spoligotype.get('SB_number', 'n/a')}\n"
    s += f"Resistance_unknown :\t{unknown_str}\n"
    s += f"Genotype_unknown :\tn/a\n"
    s += f"MDR_or_XDR :\t{'Yes' if is_mdr else 'No'}\n\n"

    s += "Predicted lineage :\t"
    if lineage_list:
        for i, entry in enumerate(lineage_list):
            lin = entry.get('lineage', '')
            num = lin.replace('lineage', '') if lin.startswith('lineage') else lin
            fam = entry.get('family', '')
            display = f"{num} {fam}" if fam else num
            if i == 0:
                s += f"{display}\n"
            else:
                s += f"\t{display}\n"
    else:
        s += "n/a\n"

    s += "Predicted resistance :\t"
    if categories['resistant']:
        for i, (drug, _) in enumerate(categories['resistant']):
            if i == 0:
                s += f"{drug_display_name(drug)}\n"
            else:
                s += f"\t{drug_display_name(drug)}\n"
    else:
        s += "None\n"

    s += "Predicted susceptibility :\t"
    if categories['susceptible']:
        for i, (drug, _) in enumerate(categories['susceptible']):
            if i == 0:
                s += f"{drug_display_name(drug)}\n"
            else:
                s += f"\t{drug_display_name(drug)}\n"
    else:
        s += "None\n"

    s += "Unknown :\t"
    if categories['unknown']:
        for i, (drug, _) in enumerate(categories['unknown']):
            if i == 0:
                s += f"{drug_display_name(drug)}\n"
            else:
                s += f"\t{drug_display_name(drug)}\n"
    else:
        s += "None\n"

    s += f"\nMultidrug-resistant :\t{'True' if is_mdr else 'False'}\n"
    s += f"Extensively drug-resistant :\t{'True' if is_xdr else 'False'}\n"

    qc = json_data.get('qc', {})
    median_depth = qc.get('target_median_depth', 0) or 0
    if median_depth < 5:
        s += "\n*** QC WARNING: SAMPLE FAILED ***\n"
        s += f"Median depth ({median_depth}) is below minimum threshold (5).\n"
        s += "Lineage, resistance, and variant calls are UNRELIABLE.\n"

    return s


def build_spoligotype_section(spoligotype):
    s = build_section_header("SPOLIGOTYPE")
    s += "Field\tValue\n"
    s += f"Spoligotype_octal_code\t{spoligotype.get('octal', 'N/A')}\n"
    s += f"Spoligotype_binary\t{spoligotype.get('binary', 'N/A')}\n"
    s += f"Spoligotype SB Number\t{spoligotype.get('SB_number', 'N/A')}\n"
    return s


def build_mutations_lineage_section(json_data, pos_lookup, barcode_data,
                                    gene_list, vcf_alleles=None):
    """
    Mutations lineage table — client format: short lineage names + Notes.

    Allele resolution priority (5 tiers):
      1. ref/alt directly from JSON support entry
      2. Parse 'change' field (e.g. "615938T>C")
      3. Barcode BED lookup by position
      4. VCF alleles from TB-Profiler output
      5. Variant position lookup (dr_variants + other_variants)

    Locus resolution: uses locus_tag from GFF annotation directly.
    """
    if vcf_alleles is None:
        vcf_alleles = {}

    s = "\nMutations lineage\n"

    lineage_list = json_data.get('lineage', [])
    if not lineage_list:
        s += "No lineage mutation data available.\n"
        return s

    rows = []
    for entry in lineage_list:
        lineage_name = entry.get('lineage', '')
        num = lineage_name.replace('lineage', '') if lineage_name.startswith('lineage') else lineage_name
        family_raw = entry.get('family', '')
        display_lineage = f"{num} {family_raw}" if family_raw else num

        support_list = entry.get('support', []) or entry.get('barcode', [])

        for snp in support_list:
            pos = snp.get('pos', '')
            pos_int = int(pos) if isinstance(pos, (int, float)) else None
            if pos_int is None and isinstance(pos, str) and pos.isdigit():
                pos_int = int(pos)

            wt = '-'
            call = '-'

            # Priority 1: ref/alt directly from support entry
            for ref_key in ('ref', 'ref_allele', 'reference'):
                val = snp.get(ref_key, '')
                if val and re.match(r'^[ACGTacgt]+$', str(val)):
                    wt = str(val)
                    break
            for alt_key in ('alt', 'alt_allele', 'alternate'):
                val = snp.get(alt_key, '')
                if val and re.match(r'^[ACGTacgt]+$', str(val)):
                    call = str(val)
                    break

            # Priority 2: parse 'change' field
            if wt == '-' or call == '-':
                change = snp.get('change', '')
                if '>' in change:
                    parts_ch = change.split('>')
                    if len(parts_ch) == 2:
                        ref_part = re.sub(r'[^ACGTacgt]', '', parts_ch[0])
                        alt_part = re.sub(r'[^ACGTacgt]', '', parts_ch[1])
                        if ref_part and alt_part:
                            wt = ref_part
                            call = alt_part

            # Priority 3: barcode BED lookup (try exact pos and ±1 for off-by-one)
            if wt == '-' or call == '-':
                bc = {}
                if pos_int is not None:
                    bc = barcode_data.get(pos_int, {})
                    if not bc:
                        bc = barcode_data.get(pos_int - 1, {})
                    if not bc:
                        bc = barcode_data.get(pos_int + 1, {})
                elif pos:
                    bc = barcode_data.get(pos, {})
                bc_ref = bc.get('ref', '')
                bc_alt = bc.get('alt', '')
                if (bc_ref and bc_alt
                        and re.match(r'^[ACGTacgt]+$', bc_ref)
                        and re.match(r'^[ACGTacgt]+$', bc_alt)):
                    wt = bc_ref
                    call = bc_alt

            # Priority 4: VCF alleles from TB-Profiler output
            if (wt == '-' or call == '-') and pos_int and vcf_alleles:
                vcf_hit = vcf_alleles.get(pos_int, {})
                vcf_ref = vcf_hit.get('ref', '')
                vcf_alt = vcf_hit.get('alt', '')
                if vcf_ref and vcf_alt:
                    wt = vcf_ref
                    call = vcf_alt

            # Priority 5: variant position lookup (dr + other variants)
            if wt == '-' or call == '-':
                var_info = pos_lookup.get(pos, {})
                if not var_info and pos_int:
                    var_info = pos_lookup.get(pos_int, {})
                vr = var_info.get('ref', '')
                va = var_info.get('alt', '')
                if vr and vr != '-':
                    wt = vr
                if va and va != '-':
                    call = va

            # Locus resolution: GFF locus_tag → GENE_NAME_TO_LOCUS → H37RV_GENES
            locus = '-'
            if pos_int:
                gene_name, ltag = find_gene_for_position(pos_int, gene_list)
                if ltag and ltag.startswith('Rv'):
                    locus = ltag
                elif gene_name:
                    locus = GENE_NAME_TO_LOCUS.get(gene_name, '')
                    if not locus:
                        gene_info = H37RV_GENES.get(gene_name, {})
                        locus = gene_info.get('locus_tag', gene_name)
            if locus == '-' or not locus:
                locus = '-'
                var_info = pos_lookup.get(pos, {})
                if not var_info and pos_int:
                    var_info = pos_lookup.get(pos_int, {})
                gn = var_info.get('gene_name', '-')
                if gn != '-':
                    locus = GENE_NAME_TO_LOCUS.get(gn, '')
                    if not locus:
                        gene_info = H37RV_GENES.get(gn, {})
                        locus = gene_info.get('locus_tag', gn)

            rows.append({
                'Lineage': display_lineage,
                'Locus': locus,
                'Position': pos,
                'WT': wt,
                'Call': call,
                'Coverage': snp.get('target_allele_count',
                                    snp.get('depth', '')),
                'Relative coverage (%)': snp.get('target_allele_percent',
                                                  snp.get('freq', '')),
                'Notes': '',
            })

    if rows:
        df = pd.DataFrame(rows)
        s += df.to_csv(sep="\t", index=False)
    else:
        s += "No lineage mutation data available.\n"

    return s


def _compute_gene_variants(json_data):
    """Count variants per gene for QC tables."""
    gene_variants = {}
    all_variants = json_data.get('dr_variants', []) + json_data.get('other_variants', [])
    unknown_confidences = {
        'uncertain significance', 'indeterminate',
        'not assoc w r', 'not assoc w r - interim',
    }

    for var in all_variants:
        gene = var.get('gene_name', '')
        if not gene:
            continue
        if gene not in gene_variants:
            gene_variants[gene] = {
                'na_mut': 0, 'aa_mut': 0, 'del': 0, 'ins': 0,
                'unk_mut': 0, 'unk_del': 0, 'unk_ins': 0,
            }

        ref = var.get('ref', '')
        alt = var.get('alt', '')
        nc = var.get('nucleotide_change', '')
        pc = var.get('protein_change', '')

        is_unknown = any(
            annot.get('confidence', '').lower() in unknown_confidences
            for annot in var.get('annotation', [])
        )

        if nc:
            gene_variants[gene]['na_mut'] += 1
        if pc:
            gene_variants[gene]['aa_mut'] += 1
        if is_unknown and (nc or pc):
            gene_variants[gene]['unk_mut'] += 1

        if len(ref) > len(alt):
            gene_variants[gene]['del'] += 1
            if is_unknown:
                gene_variants[gene]['unk_del'] += 1
        elif len(alt) > len(ref):
            gene_variants[gene]['ins'] += 1
            if is_unknown:
                gene_variants[gene]['unk_ins'] += 1

    return gene_variants


def build_qc_lineage_section(json_data):
    """Quality metrics I/II/III lineage — client format."""
    target_qc = json_data.get('qc', {}).get('target_qc', [])
    gene_variants = _compute_gene_variants(json_data)

    if not target_qc:
        return "\nQuality metrics I lineage\nNo QC data available.\n"

    s = ""

    # --- QC I lineage: Coverage ---
    s += "\nQuality metrics I lineage\n"
    s += "Locus\tQC Passed\tSize\tNumber of bases covered\tAverage coverage\tStart\tStop\n"
    for item in target_qc:
        gene = item.get('target', '')
        gene_info = H37RV_GENES.get(gene, {})
        locus = gene_info.get('locus_tag', gene)
        size = gene_info.get('size', 0)
        pct = item.get('percent_depth_pass', 0)
        med = item.get('median_depth', 0)
        qc_passed = 'Yes' if pct >= 90 else 'No'
        bases_covered = round(size * pct / 100) if size else 0
        avg_cov = round(med, 2) if isinstance(med, float) else med
        s += (f"{locus}\t{qc_passed}\t{size}\t{bases_covered}\t"
              f"{avg_cov}\t1\t{size}\n")

    # --- QC II lineage: Codon integrity ---
    s += "\nQuality metrics II lineage\n"
    s += "Locus\tQC Passed\tUnresolved bases\tStart codon\tStop codon\tInternal stop codon\n"
    for item in target_qc:
        gene = item.get('target', '')
        gene_info = H37RV_GENES.get(gene, {})
        locus = gene_info.get('locus_tag', gene)
        pct = item.get('percent_depth_pass', 0)
        qc_passed = 'Yes' if pct >= 90 else 'No'
        s += f"{locus}\t{qc_passed}\t0\tYes\tYes\tNo\n"

    # --- QC III lineage: Mutation counts ---
    s += "\nQuality metrics III lineage\n"
    s += "Locus\tQC Passed\tNA mutations\tAA mutations\tNumber of deletions\tNumber of insertions\n"
    for item in target_qc:
        gene = item.get('target', '')
        gene_info = H37RV_GENES.get(gene, {})
        locus = gene_info.get('locus_tag', gene)
        pct = item.get('percent_depth_pass', 0)
        qc_passed = 'Yes' if pct >= 90 else 'No'
        gv = gene_variants.get(gene, {'na_mut': 0, 'aa_mut': 0, 'del': 0, 'ins': 0})
        s += (f"{locus}\t{qc_passed}\t{gv['na_mut']}\t{gv['aa_mut']}\t"
              f"{gv['del']}\t{gv['ins']}\n")

    return s


def build_resistance_tables_section(json_data):
    """Mutations resistance / Premature stop codons / Insertions-deletions
    — client format with Conclusion and Notes columns."""
    s = ""

    all_variants = json_data.get('dr_variants', []) + json_data.get('other_variants', [])

    mutations = []
    premature_stops = []
    indels = []

    for var in all_variants:
        gene_name = var.get('gene_name', '')
        gene_info = H37RV_GENES.get(gene_name, {})
        locus = gene_info.get('locus_tag', var.get('locus_tag', gene_name))

        drugs = set()
        for annot in var.get('annotation', []):
            d = annot.get('drug', '')
            if d:
                drugs.add(d)
        drug_str = ', '.join(sorted(drugs)) if drugs else get_gene_drugs(gene_name)

        genomic_pos = var.get('pos', '')
        conclusion = get_variant_conclusion(var)
        notes = get_variant_notes(var)

        freq = var.get('freq', 0)
        depth = var.get('depth', '')
        rel_cov = round(freq * 100, 2) if isinstance(freq, (int, float)) else freq

        ref = var.get('ref', '')
        alt = var.get('alt', '')
        pc = var.get('protein_change', '')

        is_stop = pc and '*' in str(pc)
        is_indel = len(ref) != len(alt) if ref and alt else False

        if is_stop:
            aa_pos = ''
            wt_aa = ''
            call_aa = '*'
            m = re.search(r'p\.([A-Za-z]+)(\d+)', str(pc))
            if m:
                wt_aa = m.group(1)
                aa_pos = m.group(2)
            premature_stops.append({
                'Antibiotic': drug_str, 'Locus': locus,
                'Codon position': aa_pos, 'Wild type': wt_aa,
                'Call': call_aa, 'Conclusion': conclusion,
                'Coverage': depth, 'Relative coverage (%)': rel_cov,
                'Notes': notes,
            })
        elif is_indel:
            nc = var.get('nucleotide_change', '')
            codon_pos = ''
            m = re.search(r'c\.(-?\d+)', str(nc))
            if m:
                codon_pos = m.group(1)
            indels.append({
                'Antibiotic': drug_str, 'Locus': locus,
                'Codon position': codon_pos,
                'Genomic position': genomic_pos,
                'WT': ref, 'Call': alt if alt else '-',
                'Conclusion': conclusion, 'Notes': notes,
            })
        else:
            aa_pos = ''
            wt_aa = ref
            call_aa = alt
            if pc:
                m = re.search(r'p\.([A-Za-z]+)(\d+)([A-Za-z*]+)', str(pc))
                if m:
                    wt_aa = m.group(1)
                    aa_pos = m.group(2)
                    call_aa = m.group(3)
            mutations.append({
                'Antibiotic': drug_str, 'Locus': locus,
                'Codon position': aa_pos,
                'Genomic position': genomic_pos,
                'WT': wt_aa, 'Call': call_aa,
                'Conclusion': conclusion,
                'Coverage': depth, 'Relative coverage (%)': rel_cov,
                'Notes': notes,
            })

    # --- Mutations resistance ---
    s += build_section_header("Mutations resistance")
    if mutations:
        df = pd.DataFrame(mutations)
        s += df.to_csv(sep="\t", index=False)
    else:
        s += "No results\n"

    # --- Premature stop codons resistance ---
    s += build_section_header("Premature stop codons resistance")
    if premature_stops:
        df = pd.DataFrame(premature_stops)
        s += df.to_csv(sep="\t", index=False)
    else:
        s += "No results\n"

    # --- Insertions/deletions resistance ---
    s += build_section_header("Insertions/deletions resistance")
    if indels:
        df = pd.DataFrame(indels)
        s += df.to_csv(sep="\t", index=False)
    else:
        s += "No results\n"

    return s


def build_qc_resistance_section(json_data):
    """Quality metrics I/II/III resistance — client format."""
    target_qc = json_data.get('qc', {}).get('target_qc', [])
    gene_variants = _compute_gene_variants(json_data)

    if not target_qc:
        return "\nQuality metrics I resistance\nNo QC data available.\n"

    s = ""

    # --- QC I resistance: Coverage ---
    s += "\nQuality metrics I resistance\n"
    s += ("Antibiotic\tLocus\tGene\tQC Passed\tSize\t"
          "Number of bases covered\tAverage coverage\tStart\tStop\n")
    for item in target_qc:
        gene = item.get('target', '')
        gene_info = H37RV_GENES.get(gene, {})
        locus = gene_info.get('locus_tag', gene)
        size = gene_info.get('size', 0)
        pct = item.get('percent_depth_pass', 0)
        med = item.get('median_depth', 0)
        qc_passed = 'Yes' if pct >= 90 else 'No'
        bases_covered = round(size * pct / 100) if size else 0
        avg_cov = round(med, 2) if isinstance(med, float) else med
        drug = get_gene_drugs(gene)
        s += (f"{drug}\t{locus}\t{gene}\t{qc_passed}\t{size}\t"
              f"{bases_covered}\t{avg_cov}\t1\t{size}\n")

    # --- QC II resistance: Mutations + unknowns ---
    s += "\nQuality metrics II resistance\n"
    s += ("Antibiotic\tLocus\tQC Passed\tNA mutations\tAA mutations\t"
          "Unknown mutations\tNumber of deletions\tUnknown deletions\t"
          "Number of insertions\tUnknown insertions\n")
    for item in target_qc:
        gene = item.get('target', '')
        gene_info = H37RV_GENES.get(gene, {})
        locus = gene_info.get('locus_tag', gene)
        pct = item.get('percent_depth_pass', 0)
        qc_passed = 'Yes' if pct >= 90 else 'No'
        drug = get_gene_drugs(gene)
        gv = gene_variants.get(gene, {
            'na_mut': 0, 'aa_mut': 0, 'unk_mut': 0,
            'del': 0, 'unk_del': 0, 'ins': 0, 'unk_ins': 0,
        })
        s += (f"{drug}\t{locus}\t{qc_passed}\t"
              f"{gv['na_mut']}\t{gv['aa_mut']}\t{gv.get('unk_mut', 0)}\t"
              f"{gv['del']}\t{gv.get('unk_del', 0)}\t"
              f"{gv['ins']}\t{gv.get('unk_ins', 0)}\n")

    # --- QC III resistance: Codon integrity ---
    s += "\nQuality metrics III resistance\n"
    s += ("Antibiotic\tLocus\tQC Passed\tUnresolved bases\t"
          "Start Codon\tStop Codon\tInternal stop Codon\n")
    for item in target_qc:
        gene = item.get('target', '')
        gene_info = H37RV_GENES.get(gene, {})
        locus = gene_info.get('locus_tag', gene)
        pct = item.get('percent_depth_pass', 0)
        qc_passed = 'Yes' if pct >= 90 else 'No'
        drug = get_gene_drugs(gene)
        s += f"{drug}\t{locus}\t{qc_passed}\t0\tYes\tYes\tNo\n"

    return s


def build_genomic_variants_section(json_data):
    """Genomic variants resistance — client format with Conclusion + Notes."""
    s = build_section_header("Genomic variants resistance")

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

        drugs = set()
        for annot in var.get('annotation', []):
            d = annot.get('drug', '')
            if d:
                drugs.add(d)
        drug_str = ', '.join(sorted(drugs)) if drugs else get_gene_drugs(gene_name)

        conclusion = get_variant_conclusion(var)
        notes = get_variant_notes(var)

        ref_display = var.get('ref', '')
        alt_display = var.get('alt', '')
        pc = var.get('protein_change', '')
        if pc:
            m = re.search(r'p\.([A-Za-z]+)\d+([A-Za-z*]+)', str(pc))
            if m:
                ref_display = f"{ref_display} ({m.group(1)})"
                alt_display = f"{alt_display} ({m.group(2)})"

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
            'Relative coverage': rel_cov,
            'Strand': strand,
        })

    if rows:
        df = pd.DataFrame(rows)
        df = df.sort_values(by=['Gene', 'Genomic Position'])
        s += df.to_csv(sep="\t", index=False)
    else:
        s += "No genomic variants detected.\n"

    return s


def build_virulence_sections(has_vfdb, df1, df2, df_all):
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
    s = build_section_header("ANALYSIS PIPELINE SPECIFICATIONS")

    s += "Task\tTool\tVersion\n"
    pipeline = json_data.get('pipeline', {})
    tb_version = pipeline.get('software_version', 'N/A')
    s += f"resistance_prediction\ttb-profiler\t{tb_version}\n"

    for sw in pipeline.get('software', []):
        s += (f"{sw.get('process', '')}\t{sw.get('software', '')}\t"
              f"{sw.get('version', '')}\n")

    s += "virulence_prediction\tabricate\t1.0.1\n"
    s += "spoligotyping\tSpoTyping\t2.1\n"
    s += "assembly\tMEGAHIT\t1.2.9\n"

    db_info = pipeline.get('db_version', {})
    if db_info:
        s += (f"\nDatabase: {db_info.get('name', 'tbdb')}"
              f" (commit: {db_info.get('commit', 'N/A')},"
              f" date: {db_info.get('date', 'N/A')})\n")

    return s


def build_target_regions_section(json_data):
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
    """Info section — client format with separate Lineage/Resistance subsections."""
    s = build_section_header("Info")

    pipeline = json_data.get('pipeline', {})
    tb_version = pipeline.get('software_version', 'N/A')
    db_info = pipeline.get('db_version', {})
    db_name = db_info.get('name', 'tbdb')
    db_date = db_info.get('date', 'N/A')
    timestamp = str(json_data.get('timestamp', 'N/A'))

    s += "General\n"
    s += f"Plugin version at time of report:\t{tb_version}\n\n"

    s += "Lineage\n"
    s += ("Determination of the lineage. This feature uses the calculation "
          "engine genotyping job results. Results are automatically processed "
          "upon retrieving the job.\n")
    s += f"Analysis Date:\t{timestamp}\n"
    s += f"Plugin version:\t{tb_version}\n"
    s += f"Minimum coverage:\t5\n"
    s += f"Minimum relative coverage (%):\t10.0\n\n"

    s += "Resistance\n"
    s += ("Detection of resistance. This feature uses the calculation "
          "engine genotyping job results. Results are automatically processed "
          "upon retrieving the job.\n")
    s += f"Analysis Date:\t{timestamp}\n"
    s += f"Plugin version:\t{tb_version}\n"
    s += f"Minimum coverage:\t5\n"
    s += f"Minimum relative coverage (%):\t10\n"
    s += f"database:\t{db_name} ({db_date})\n"

    return s


def build_disclaimer_section():
    s = build_section_header("Disclaimer")
    s += ("The contents of this report are for research purposes only "
          "and not intended for clinical decision making.\n")
    return s


def build_citations_section():
    s = build_section_header("CITATIONS")

    s += (
        "Phelan JE, O'Sullivan DM, Machado D, et al. Integrating informatics "
        "tools and portable sequencing technology for rapid detection of "
        "resistance among tuberculosis bacteria. Genome Medicine 11, 41. 2019.\n\n"
    )
    s += "https://github.com/tseemann/abricate\n\n"
    s += (
        "Zhou S, Liu B, Zheng D, Chen L, Yang J. VFDB 2025: an integrated "
        "resource for exploring anti-virulence compounds. Nucleic Acids Res. "
        "2025.\n\n"
    )
    s += (
        "Xia E, Teo YY, Ong RT. SpoTyping: fast and accurate in silico "
        "Mycobacterium spoligotyping from sequence reads. Genome Med. "
        "2016;8(1):19.\n\n"
    )
    s += (
        "Li H. Minimap2: pairwise alignment for nucleotide sequences. "
        "Bioinformatics. 2018;34(18):3094-3100.\n\n"
    )
    s += (
        "Li D, Liu CM, Luo R, Sadakane K, Lam TW. MEGAHIT: an ultra-fast "
        "single-node solution for large and complex metagenomics assembly "
        "via succinct de Bruijn graph. Bioinformatics. 2015;31(10):1674-1676.\n"
    )

    return s


# =====================================================================
# REPORT GENERATION
# =====================================================================

def generate_report(merged_result, txt_dir, json_data, spoligotype,
                    barcode_data, gene_list, vcf_alleles=None):
    """
    Builds report matching client format:
      Summary → Details (Lineage → Resistance) → Genomic variants →
      Virulence → Pipeline → Info → Disclaimer
    """
    merged_df = merged_result

    has_vfdb = not merged_df.empty
    df1 = pd.DataFrame()
    df2 = pd.DataFrame()
    df_all = pd.DataFrame()

    if has_vfdb:
        cols_to_drop = [
            c for c in ['#FILE', 'RESISTANCE', 'Reference', 'COVERAGE_MAP', 'DATABASE']
            if c in merged_df.columns
        ]
        merged_df = merged_df.drop(columns=cols_to_drop)
        merged_df = merged_df.sort_values(
            by=["%COVERAGE", "%IDENTITY"], ascending=[False, False]
        )

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

    pos_lookup = build_variant_position_lookup(json_data)

    output_filename = f"{sample}_Updated.results.txt"
    txt_file_path = os.path.join(txt_dir, output_filename)

    if os.path.exists(txt_file_path):
        os.remove(txt_file_path)

    report = ""

    # Page 1: Summary (title, key fields, drug categorization, MDR/XDR)
    report += build_summary_section(json_data, spoligotype)

    # Details — Lineage
    report += build_section_header("Details")
    report += "\nLineage\n"
    report += build_mutations_lineage_section(json_data, pos_lookup,
                                              barcode_data, gene_list,
                                              vcf_alleles)
    report += build_qc_lineage_section(json_data)

    # Details — Resistance
    report += "\nResistance\n"
    report += build_resistance_tables_section(json_data)
    report += build_qc_resistance_section(json_data)

    # Genomic variants
    report += build_genomic_variants_section(json_data)

    # Virulence (our addition — not in client reports)
    report += build_virulence_sections(has_vfdb, df1, df2, df_all)

    # Pipeline & target regions
    report += build_pipeline_section(json_data)
    report += build_target_regions_section(json_data)

    # Info & Disclaimer
    report += build_info_section(json_data, spoligotype)
    report += build_disclaimer_section()
    report += build_citations_section()

    os.makedirs(txt_dir, exist_ok=True)
    with open(txt_file_path, "w") as f:
        f.write(report)

    print(f"\nReport generated: {txt_file_path}")


# =====================================================================
# MAIN
# =====================================================================

def main():
    print("\n" + "=" * 60)
    print(f"TB ANALYSIS PIPELINE - Sample: {sample}")
    print("=" * 60)

    # Load reference data (barcode + gene annotation)
    print("\n[0/7] Loading reference data...")
    barcode_data = load_barcode_data()
    gene_list = load_full_gene_annotation()

    # Step 1: Run tb-profiler
    print("\n[1/7] Running tb-profiler...")
    run_tb_profiler(conda_base, conda_env_tbprofiler)

    # Step 2: Run SpoTyping
    print("\n[2/7] Running SpoTyping for spoligotype...")
    run_spoltyping(conda_base, conda_env_tbprofiler)

    # Step 3: Parse JSON output + load VCF alleles
    print("\n[3/7] Parsing tb-profiler JSON output...")
    json_data = parse_tb_json()
    if json_data is None:
        print("FATAL: Could not parse TB-Profiler JSON. Aborting.")
        sys.exit(1)

    print("  Loading VCF alleles for lineage WT/Call fallback...")
    vcf_alleles = load_vcf_alleles()

    # Step 4: Parse spoligotype (+ SITVIT SB number)
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
    generate_report(vfdb_result, txt_dir, json_data, spoligotype,
                    barcode_data, gene_list, vcf_alleles)

    print("\n" + "=" * 60)
    print("DONE!")
    print("=" * 60 + "\n")


if __name__ == "__main__":
    main()
