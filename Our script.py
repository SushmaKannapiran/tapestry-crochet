#Usage 
"""
Description:
This script analyzes tuberculosis sequencing data to identify drug resistance
variants and antimicrobial resistance (AMR) genes.

It performs:
1. TB-Profiler analysis to detect Mycobacterium tuberculosis drug resistance variants.
2. Read assembly using MEGAHIT followed by AMR/virulence gene detection using
   ABRicate with the VFDB database.

Usage
-----
python3 Tuberculosis_ARG_Oman.py <sample_dir> <sample_name> <output_prefix> <threads>

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
import pandas as pd

#############################
#Conda env for activation
#############################
conda_base = "/media/carme/Idaa/conda/anaconda3/etc/profile.d/conda.sh"
conda_env_tbprofiler = "tb_profiler_env"

abr_env_conda = "/media/carme/ganesh/miniconda3/etc/profile.d/conda.sh"
conda_env_abricate = "abricate_env"
##################
sample_dir=sys.argv[1]
sample=sys.argv[2]
output_prefix=sys.argv[3]
threads =sys.argv[4]

fastq_file = os.path.join(sample_dir, f"{sample}.fastq")


################
#Output paths###
################

tb_profiler_dir = os.path.join(sample_dir, f"{sample}_tbprofler")
megahit_out_dir = f"{sample_dir}/Megahit/"
abricate_result = f"{sample_dir}/{sample}_abricate_vfdb.tab"
vfdb_annot = "/media/carme/Research_Projects/BugSpeaks_MicrobiomeAnalysis/Tuberculosis_Firoz/VFDB/VFs_VFDB.csv"

######################################################################
# Step 1: Activate tb-profiler environment and run tb-profiler profile

def run_tb_profiler(conda_path,env_name):
    cmd = (
        f"source {conda_path} && "
        f"conda activate {env_name} && "
        f"tb-profiler profile --read1 {fastq_file} --platform nanopore --mapper minimap2 "
        f"--caller freebayes --prefix {sample} --dir {tb_profiler_dir} --call_whole_genome "
        f"--depth 5 --af 0.10 --suspect --txt --csv"
    )

    subprocess.run(cmd, shell=True, executable="/bin/bash", check=True)


# Step 2: Megahit 
def run_megahit():
    subprocess.run(f"megahit -t {threads} -r {fastq_file} -o {megahit_out_dir}", shell=True, check=True)
    shutil.move(f"{megahit_out_dir}/final.contigs.fa", f"{sample}.fa")
    shutil.rmtree(megahit_out_dir)  # Clean up the megahit directory

# Step 3: Abricate [VFDB]
def run_abricate(conda_path,env_name):
    cmd = (
        f"source {conda_path} && "
        f"conda activate {env_name} && "
        f"abricate --db vfdb {sample}.fa > {abricate_result}"
    )

    subprocess.run(cmd, shell=True, executable="/bin/bash", check=True)


# Step 4: Merge VFDB profiler results (extracting and merging with Python)
def extract_and_merge_vfdb_results():
    import pandas as pd
    abricate_df = pd.read_csv(abricate_result, sep="\t")
    vfdb_df = pd.read_csv(vfdb_annot)

    # Extract VFID from the 'PRODUCT' column in abricate result (13th column)
    abricate_df['VFID'] = abricate_df['PRODUCT'].apply(lambda x: re.search(r'\(VF[0-9]+\)', str(x)).group(0)[1:-1] if isinstance(x, str) and re.search(r'\(VF[0-9]+\)', x) else None)
    print(abricate_df.head())

    # Merge the two dataframes based on the VFID column
    merged_result = pd.merge(abricate_df, vfdb_df, how='left', left_on='VFID', right_on='VFID').fillna("NA")
    print("MergedResult",merged_result.shape)
    return merged_result 

# Step5 : Append VFDB results into Tb_profiler results 
def append_to_txt_file(merged_result, txt_dir):
    merged_df = merged_result

    # Drop unwanted columns: '#FILE', 'RESISTANCE', 'Reference', and 'COVERAGE_MAP'
    merged_df = merged_df.drop(columns=['#FILE', 'RESISTANCE', 'Reference', 'COVERAGE_MAP','DATABASE'])
    merged_df = merged_df.sort_values(by=["%COVERAGE", "%IDENTITY"], ascending=[False, False])
    
    df1_sig = merged_df[(merged_df["%COVERAGE"] > 95) & (merged_df["%IDENTITY"] > 95)]

    df1 = df1_sig.iloc[:,0:12]
    df2 = df1_sig.iloc[:,11:].drop_duplicates()

    #Factors with < 95% coverage , identity
    df1_notsig = merged_df[(merged_df["%COVERAGE"] < 95) & (merged_df["%IDENTITY"] < 95)]


    df_all = df1_notsig.iloc[:,0:12]

    # Convert the dataframe to a tab-separated string (without the index)
    #merged_str = merged_df.to_csv(sep='\t', index=False, header=False)

    # Define the path for the txt file (appending to the first .txt file found in the directory)
    txt_file = next(f for f in os.listdir(txt_dir) if f.endswith('.txt'))

    src = os.path.join(txt_dir, txt_file)
    copy = os.path.join(txt_dir, f"{output_prefix}_{txt_file}")
    shutil.copy(src, copy)

    txt_file_path = copy

    with open(txt_file_path, "r") as f:
        content = f.read()

    # Identify pipeline block
    pipeline_start = content.find("Analysis pipeline specifications")
    if pipeline_start == -1:
        raise ValueError("Pipeline block not found.")

    # Extract: before, pipeline, after
    pipeline_end = content.find("Citation", pipeline_start)  # Start of citation
    citation_start = pipeline_end

    # Find end of citation section
    citation_end_match = re.search(
        r"Genome Medicine 11, 41\. 2019", content[citation_start:], re.MULTILINE
    )
    citation_end = citation_start + citation_end_match.end()

    pipeline_block = content[pipeline_start:citation_end]
    before_pipeline = content[:pipeline_start]
    after_pipeline = content[citation_end:]

    # --- Modify pipeline: add the abricate line BEFORE "Citation" ---
    pipeline_block = re.sub(
    r"^(depth_calculation.*)$",
    r"\1\nvirulence_prediction\tabricate\t1.0.1",
    pipeline_block,
    flags=re.MULTILINE
    )


    # --- Append extra citation items ---
    pipeline_block += (
        "\n\nhttps://github.com/tseemann/abricate\n\n"
        "Zhou S, Liu B, Zheng D, Chen L, Yang J. VFDB 2025: an integrated resource for "
        "exploring anti-virulence compounds. Nucleic Acids Res. 2025 Jan 6;53(D1):D871-D877. "
        "doi: 10.1093/nar/gkae968. PMID: 39470738; PMCID: PMC11701737.\n"
    )

    # --- Build virulence section ---
    virulence_sections = ""
    virulence_sections += "\n\nSignificant Virulence Factors\n" + "-"*50 + "\n"
    virulence_sections += df1.to_csv(sep="\t", index=False)

    virulence_sections += "\n\nSignificant Virulence Factors - Description\n" + "-"*50 + "\n"
    virulence_sections += df2.to_csv(sep="\t", index=False)

    virulence_sections += "\n\nOther Virulence Factors\n" + "-"*50 + "\n"
    virulence_sections += df_all.to_csv(sep="\t", index=False)

    # --- Final Assembly ---
    final_output = before_pipeline + virulence_sections + "\n\n" + pipeline_block + after_pipeline

    # Write back
    with open(txt_file_path, "w") as f:
        f.write(final_output)

    print("Update complete!")


############
merged_vfdb_results_csv = f'{sample_dir}/{sample}_merged_vfdb_results.csv'

# Example usage:
#merged_vfdb_results_csv = f'{sample_dir}/{sample}_merged_vfdb_results.csv'
tb_dir = f"{sample_dir}/{sample}_tbprofler/"
txt_dir = f"{tb_dir}/results/" 

############################################################################


def main():
    # Step 1: Activate tb-profiler environment and run profiling
    run_tb_profiler(conda_base,conda_env_tbprofiler)

    # Step 2: Run Megahit
    run_megahit()

    # Step 3: Activate abricate environment and run abricate
    run_abricate(abr_env_conda,conda_env_abricate)   

    # Step 4: Extract and merge VFDB results
    result=extract_and_merge_vfdb_results()

    #Step 5: Append the VFDB results to tb_profiler result 
    append_to_txt_file(result, txt_dir)

if __name__ == "__main__":
    main()
