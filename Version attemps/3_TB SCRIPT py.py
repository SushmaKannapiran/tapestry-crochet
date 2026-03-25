#Usage 
"""
Description:
This script analyzes tuberculosis sequencing data to identify drug resistance
variants and antimicrobial resistance (AMR) genes.

It performs:
1. TB-Profiler analysis to detect Mycobacterium tuberculosis drug resistance variants.
2. SpoTyping to get spoligotype code (NEW).
3. Read assembly using MEGAHIT followed by AMR/virulence gene detection using
   ABRicate with the VFDB database.
4. Parse JSON output to extract detailed lineage SNPs, QC metrics, and resistance data (NEW).

===================================================================================
WHAT WAS CHANGED:
===================================================================================

CHANGE 1: Added --json flag to tb-profiler command (line 79)
    WHY: To get detailed output file with all the data we need

CHANGE 2: Added SpoTyping function (lines 85-97)
    WHY: tb-profiler cannot calculate spoligotype. SpoTyping does this.
    INSTALL: conda install -c bioconda spoltyping -y

CHANGE 3: Added JSON parsing function (lines 103-182)
    WHY: The JSON file has all the detailed data (lineage SNPs, resistance 
         mutations with confidence, QC metrics) but old script didn't read it.

CHANGE 4: Added function to create lineage mutations table (lines 188-214)
    WHY: Client report shows SNPs that define the lineage. This data is in 
         JSON under "lineage" -> "support"

CHANGE 5: Added function to create QC tables (lines 220-250)
    WHY: Client report has Quality Metrics tables. This data is in JSON 
         under "qc" -> "target_qc"

CHANGE 6: Added function to categorize drugs (lines 256-290)
    WHY: Client report separates drugs into Resistant/Susceptible/Unknown.
         We use "confidence" field from JSON to categorize.

CHANGE 7: Updated report generation (lines 350-500)
    WHY: Add all new sections to final report in client format.

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
# Step 1: Run tb-profiler
# CHANGED: Added --json flag
######################################################################
def run_tb_profiler(conda_path, env_name):
    cmd = (
        f"source {conda_path} && "
        f"conda activate {env_name} && "
        f"tb-profiler profile --read1 {fastq_file} --platform nanopore --mapper minimap2 "
        f"--caller freebayes --prefix {sample} --dir {tb_profiler_dir} --call_whole_genome "
        f"--depth 5 --af 0.10 --suspect --txt --csv --json"  # ADDED --json
    )
    subprocess.run(cmd, shell=True, executable="/bin/bash", check=True)


######################################################################
# Step 2: NEW - Run SpoTyping for spoligotype
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
# Step 3: NEW - Parse JSON file
######################################################################
def parse_tb_json():
    """
    Reads the JSON file and extracts all detailed data.
    
    The JSON structure (from TB-Profiler v6.6.5):
    - "main_lineage": "lineage1"
    - "sub_lineage": "lineage1.1.3.3"
    - "drtype": "HR-TB"
    - "lineage": [...] contains SNPs that define lineage
    - "dr_variants": [...] contains resistance mutations
    - "other_variants": [...] contains other mutations
    - "qc": {"target_qc": [...]} contains coverage data
    """
    json_file = f"{tb_profiler_dir}/results/{sample}.results.json"
    
    if not os.path.exists(json_file):
        print(f"Warning: JSON file not found at {json_file}")
        return None
    
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    return data


######################################################################
# Step 4: NEW - Parse spoligotype result
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
# Step 5: NEW - Create Lineage Mutations Table
######################################################################
def create_lineage_table(data):
    """
    Creates table showing SNPs that define the lineage.
    
    Client report format:
    Lineage | Locus | Position | WT | Call | Coverage | Relative coverage %
    
    Data comes from JSON: "lineage" -> "support"
    """
    rows = []
    
    for lineage_info in data.get('lineage', []):
        lineage_name = lineage_info.get('lineage', '')
        family = lineage_info.get('family', '')
        
        for snp in lineage_info.get('support', []):
            rows.append({
                'Lineage': f"{lineage_name} ({family})",
                'Position': snp.get('pos', ''),
                'Coverage': snp.get('target_allele_count', ''),
                'Relative_coverage_%': snp.get('target_allele_percent', '')
            })
    
    return pd.DataFrame(rows) if rows else pd.DataFrame()


######################################################################
# Step 6: NEW - Create QC Tables
######################################################################
def create_qc_table(data):
    """
    Creates Quality Metrics table.
    
    Client report format:
    Locus | QC Passed | Percent Depth Pass | Median Depth
    
    Data comes from JSON: "qc" -> "target_qc"
    """
    rows = []
    
    qc_data = data.get('qc', {})
    target_qc = qc_data.get('target_qc', [])
    
    for item in target_qc:
        rows.append({
            'Locus': item.get('target', ''),
            'QC_Passed': 'Yes' if item.get('percent_depth_pass', 0) >= 90 else 'No',
            'Percent_Depth_Pass': item.get('percent_depth_pass', ''),
            'Median_Depth': item.get('median_depth', '')
        })
    
    return pd.DataFrame(rows) if rows else pd.DataFrame()


######################################################################
# Step 7: NEW - Create Resistance Mutations Table
######################################################################
def create_resistance_table(data):
    """
    Creates detailed resistance mutations table.
    
    Client report format:
    Drug | Locus | Position | WT | Call | Confidence | Coverage | Notes
    
    Data comes from JSON: "dr_variants"
    """
    rows = []
    
    # Confidence ranking: higher index = stronger evidence
    confidence_rank = {
        'not assoc w r': 0,
        'not assoc w r - interim': 0,
        'uncertain significance': 1,
        'indeterminate': 1,
        'assoc w r - interim': 2,
        'assoc w r': 3,
    }

    for var in data.get('dr_variants', []):
        # Get drugs from annotation, keeping the HIGHEST confidence
        drugs = []
        best_confidence = ''
        best_rank = -1
        comment = ''

        for annot in var.get('annotation', []):
            if annot.get('drug'):
                drugs.append(annot.get('drug'))
                conf = annot.get('confidence', '')
                rank = confidence_rank.get(conf.lower(), 1)
                if rank > best_rank:
                    best_rank = rank
                    best_confidence = conf
                if annot.get('comment', ''):
                    comment = annot.get('comment', '')

        rows.append({
            'Drug': ', '.join(set(drugs)),
            'Locus': var.get('gene_name', ''),
            'Locus_Tag': var.get('locus_tag', ''),
            'Position': var.get('pos', ''),
            'WT': var.get('ref', ''),
            'Call': var.get('alt', ''),
            'Change': var.get('change', ''),
            'Confidence': best_confidence,
            'Coverage': var.get('depth', ''),
            'Frequency': var.get('freq', ''),
            'Notes': comment
        })
    
    return pd.DataFrame(rows) if rows else pd.DataFrame()


######################################################################
# Step 8: NEW - Categorize drugs into Resistant/Susceptible/Unknown
######################################################################
def categorize_drugs(data):
    """
    Separates drugs into three categories based on confidence level.
    
    - Resistant: "Assoc w R" confidence
    - Unknown: "Uncertain significance" or low confidence
    - Susceptible: Everything else
    """
    all_drugs = {
        'rifampicin', 'isoniazid', 'ethambutol', 'pyrazinamide',
        'moxifloxacin', 'levofloxacin', 'ofloxacin',
        'bedaquiline', 'delamanid', 'linezolid', 'streptomycin',
        'amikacin', 'kanamycin', 'capreomycin', 'clofazimine',
        'ethionamide', 'para-aminosalicylic_acid', 'cycloserine'
    }
    
    resistant = set()
    unknown = set()
    
    # Check dr_variants
    for var in data.get('dr_variants', []):
        for annot in var.get('annotation', []):
            drug = annot.get('drug', '').lower()
            confidence = annot.get('confidence', '').lower()
            
            if 'assoc w r' in confidence:
                resistant.add(drug)
            elif 'uncertain' in confidence or 'indeterminate' in confidence:
                unknown.add(drug)
    
    # Check other_variants for unknown
    for var in data.get('other_variants', []):
        for annot in var.get('annotation', []):
            drug = annot.get('drug', '').lower()
            confidence = annot.get('confidence', '').lower()
            
            if 'uncertain' in confidence:
                unknown.add(drug)
    
    # If a drug is in resistant, remove it from unknown
    unknown = unknown - resistant
    susceptible = all_drugs - resistant - unknown

    return {
        'resistant': sorted(list(resistant)),
        'susceptible': sorted(list(susceptible)),
        'unknown': sorted(list(unknown))
    }


######################################################################
# Step 9: Megahit (UNCHANGED)
######################################################################
def run_megahit():
    subprocess.run(f"megahit -t {threads} -r {fastq_file} -o {megahit_out_dir}", shell=True, check=True)
    shutil.move(f"{megahit_out_dir}/final.contigs.fa", f"{sample}.fa")
    shutil.rmtree(megahit_out_dir)


######################################################################
# Step 10: Abricate (UNCHANGED)
######################################################################
def run_abricate(conda_path, env_name):
    cmd = (
        f"source {conda_path} && "
        f"conda activate {env_name} && "
        f"abricate --db vfdb {sample}.fa > {abricate_result}"
    )
    subprocess.run(cmd, shell=True, executable="/bin/bash", check=True)


######################################################################
# Step 11: Merge VFDB results (UNCHANGED)
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
# Step 12: Generate Final Report (UPDATED)
######################################################################
def generate_report(merged_result, txt_dir, json_data, spoligotype):
    """
    Generates final report with all new sections matching client format.
    """
    merged_df = merged_result

    # Process VFDB results (handle empty abricate output)
    has_vfdb = not merged_df.empty
    if has_vfdb:
        merged_df = merged_df.drop(columns=['#FILE', 'RESISTANCE', 'Reference', 'COVERAGE_MAP', 'DATABASE'])
        merged_df = merged_df.sort_values(by=["%COVERAGE", "%IDENTITY"], ascending=[False, False])

        # Nanopore-adjusted thresholds: 80% coverage, 80% identity
        # (Nanopore has higher per-read error rates than Illumina,
        # so 95% identity is too strict and almost nothing passes)
        sig_cov_thresh = 80
        sig_id_thresh = 80

        df1_sig = merged_df[(merged_df["%COVERAGE"] >= sig_cov_thresh) & (merged_df["%IDENTITY"] >= sig_id_thresh)]
        df1 = df1_sig.iloc[:, 0:12]
        df2 = df1_sig.iloc[:, 11:].drop_duplicates()
        df1_notsig = merged_df[(merged_df["%COVERAGE"] < sig_cov_thresh) | (merged_df["%IDENTITY"] < sig_id_thresh)]
        df_all = df1_notsig.iloc[:, 0:12]
    
    # Find and copy txt file — pick the ORIGINAL tb-profiler output only
    # (exclude any previously generated report files with output_prefix)
    txt_files = [f for f in os.listdir(txt_dir) if f.endswith('.txt') and not f.startswith(output_prefix)]
    if not txt_files:
        # Fallback: pick any txt file
        txt_files = [f for f in os.listdir(txt_dir) if f.endswith('.txt')]
    txt_file = txt_files[0]
    src = os.path.join(txt_dir, txt_file)
    copy = os.path.join(txt_dir, f"{output_prefix}_{txt_file}")
    shutil.copy(src, copy)
    txt_file_path = copy
    
    with open(txt_file_path, "r") as f:
        content = f.read()
    
    # Find pipeline section
    pipeline_start = content.find("Analysis pipeline specifications")
    if pipeline_start == -1:
        raise ValueError("Pipeline block not found.")
    
    pipeline_end = content.find("Citation", pipeline_start)
    citation_start = pipeline_end
    citation_end_match = re.search(r"Genome Medicine 11, 41\. 2019", content[citation_start:], re.MULTILINE)
    citation_end = citation_start + citation_end_match.end()
    
    pipeline_block = content[pipeline_start:citation_end]
    before_pipeline = content[:pipeline_start]
    # Discard anything after the citation — prevents duplication if the
    # txt file was from a previous run that already had appended sections
    after_pipeline = ""
    
    # Update pipeline block (count=1 to prevent duplicate insertion)
    pipeline_block = re.sub(
        r"^(depth_calculation.*)$",
        r"\1\nvirulence_prediction\tabricate\t1.0.1\nspoligotyping\tSpoTyping\t2.1",
        pipeline_block,
        count=1,
        flags=re.MULTILINE
    )
    
    # Add citations
    pipeline_block += (
        "\n\nhttps://github.com/tseemann/abricate\n\n"
        "Zhou S, Liu B, Zheng D, Chen L, Yang J. VFDB 2025: an integrated resource for "
        "exploring anti-virulence compounds. Nucleic Acids Res. 2025.\n\n"
        "Xia E, Teo YY, Ong RT. SpoTyping: fast and accurate in silico Mycobacterium "
        "spoligotyping from sequence reads. Genome Med. 2016;8(1):19.\n"
    )
    
    # =====================================================================
    # BUILD NEW SECTIONS (MATCHING CLIENT FORMAT)
    # =====================================================================
    
    new_sections = ""
    
    # --- SUMMARY SECTION ---
    new_sections += "\n" + "="*60 + "\n"
    new_sections += "SUMMARY\n"
    new_sections += "="*60 + "\n\n"
    
    median_depth = json_data.get('qc', {}).get('target_median_depth', 0) or 0
    sample_failed = median_depth < 5

    new_sections += f"Sample Name:\t\t{json_data.get('id', sample)}\n"
    new_sections += f"Main Lineage:\t\t{json_data.get('main_lineage', 'N/A')}\n"
    new_sections += f"Sub Lineage:\t\t{json_data.get('sub_lineage', 'N/A')}\n"
    new_sections += f"Drug Resistance Type:\t{json_data.get('drtype', 'N/A')}\n"
    new_sections += f"Median Depth:\t\t{median_depth}\n"

    if sample_failed:
        new_sections += "\n*** QC WARNING: SAMPLE FAILED ***\n"
        new_sections += f"Median depth ({median_depth}) is below minimum threshold (5).\n"
        new_sections += "Lineage, resistance, and variant calls are UNRELIABLE.\n"
        new_sections += "Resistance calls (especially rrs/rrl) may be FALSE POSITIVES.\n"

    # --- SPOLIGOTYPE SECTION ---
    new_sections += "\n" + "="*60 + "\n"
    new_sections += "SPOLIGOTYPE\n"
    new_sections += "="*60 + "\n\n"
    
    new_sections += f"Spoligotype_octal_code:\t\t{spoligotype.get('octal', 'N/A')}\n"
    new_sections += f"Spoligotype_binary:\t\t{spoligotype.get('binary', 'N/A')}\n"
    new_sections += f"Spoligotype_SB_number:\t\t{spoligotype.get('SB_number', 'N/A')}\n"
    
    # --- DRUG RESISTANCE SUMMARY ---
    drug_categories = categorize_drugs(json_data)
    
    new_sections += "\n" + "="*60 + "\n"
    new_sections += "DRUG RESISTANCE SUMMARY\n"
    new_sections += "="*60 + "\n\n"
    
    new_sections += f"Predicted Resistance:\t\t{', '.join(drug_categories['resistant']) if drug_categories['resistant'] else 'None'}\n"
    new_sections += f"Predicted Susceptibility:\t{', '.join(drug_categories['susceptible']) if drug_categories['susceptible'] else 'None'}\n"
    new_sections += f"Unknown:\t\t\t{', '.join(drug_categories['unknown']) if drug_categories['unknown'] else 'None'}\n"
    
    # --- MUTATIONS LINEAGE TABLE ---
    new_sections += "\n" + "="*60 + "\n"
    new_sections += "MUTATIONS LINEAGE\n"
    new_sections += "="*60 + "\n\n"
    
    lineage_df = create_lineage_table(json_data)
    if not lineage_df.empty:
        new_sections += lineage_df.to_csv(sep="\t", index=False)
    else:
        new_sections += "No lineage mutation data available.\n"
    
    # --- MUTATIONS RESISTANCE TABLE ---
    new_sections += "\n" + "="*60 + "\n"
    new_sections += "MUTATIONS RESISTANCE\n"
    new_sections += "="*60 + "\n\n"
    
    resistance_df = create_resistance_table(json_data)
    if not resistance_df.empty:
        new_sections += resistance_df.to_csv(sep="\t", index=False)
    else:
        new_sections += "No resistance mutations detected.\n"
    
    # --- QUALITY METRICS TABLE ---
    new_sections += "\n" + "="*60 + "\n"
    new_sections += "QUALITY METRICS - COVERAGE\n"
    new_sections += "="*60 + "\n\n"
    
    qc_df = create_qc_table(json_data)
    if not qc_df.empty:
        new_sections += qc_df.to_csv(sep="\t", index=False)
    else:
        new_sections += "No QC data available.\n"
    
    # --- VIRULENCE SECTIONS ---
    virulence_sections = ""

    virulence_sections += "\n" + "="*60 + "\n"
    virulence_sections += "SIGNIFICANT VIRULENCE FACTORS\n"
    virulence_sections += "="*60 + "\n\n"
    if has_vfdb and not df1.empty:
        virulence_sections += df1.to_csv(sep="\t", index=False)
    else:
        virulence_sections += "No significant virulence factors detected.\n"

    virulence_sections += "\n" + "="*60 + "\n"
    virulence_sections += "SIGNIFICANT VIRULENCE FACTORS - DESCRIPTION\n"
    virulence_sections += "="*60 + "\n\n"
    if has_vfdb and not df2.empty:
        virulence_sections += df2.to_csv(sep="\t", index=False)
    else:
        virulence_sections += "No descriptions available.\n"

    virulence_sections += "\n" + "="*60 + "\n"
    virulence_sections += "OTHER VIRULENCE FACTORS\n"
    virulence_sections += "="*60 + "\n\n"
    if has_vfdb and not df_all.empty:
        virulence_sections += df_all.to_csv(sep="\t", index=False)
    else:
        virulence_sections += "No other virulence factors detected.\n"
    
    # --- FINAL ASSEMBLY ---
    final_output = before_pipeline + new_sections + virulence_sections + "\n\n" + pipeline_block + after_pipeline
    
    with open(txt_file_path, "w") as f:
        f.write(final_output)
    
    print(f"\nReport generated: {txt_file_path}")


######################################################################
# MAIN FUNCTION
######################################################################
def main():
    print("\n" + "="*60)
    print(f"TB ANALYSIS PIPELINE - Sample: {sample}")
    print("="*60)
    
    # Step 1: Run tb-profiler (now with --json)
    print("\n[1/7] Running tb-profiler...")
    run_tb_profiler(conda_base, conda_env_tbprofiler)
    
    # Step 2: Run SpoTyping (NEW)
    print("\n[2/7] Running SpoTyping for spoligotype...")
    run_spoltyping(conda_base, conda_env_tbprofiler)
    
    # Step 3: Parse JSON output (NEW)
    print("\n[3/7] Parsing tb-profiler JSON output...")
    json_data = parse_tb_json()
    
    # Step 4: Parse spoligotype (NEW)
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
    
    print("\n" + "="*60)
    print("DONE!")
    print("="*60 + "\n")


if __name__ == "__main__":
    main()
