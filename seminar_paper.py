# Umutcan Unaldi
# 7025677
# Saarland University
# Bioinformatics M.Sc.
# Seminar Paper

import os
import pandas as pd
import openpyxl
import warnings
import subprocess
import requests
from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaWriter
import re
from collections import Counter

# Get rid of unknown extension warning
warnings.simplefilter("ignore", UserWarning)

# mmc3.xlsx consists of all the interactors and isoforms
df_mmc3 = pd.read_excel("mmc3.xlsx", sheet_name="2B-Isoform PPIs")

# Extract interactors
interactors = df_mmc3.iloc[2:, 6].unique()

# Unconventional naming by the authors is fixed
interactors[interactors == "LOC100288797"] = "TMEM239"

def search(gene_name):
    ''' Searches for interactors in UniProt database
    and returns information like accession number, gene name,
    synonyms of the all found isoforms of that interactor [1].'''

    # Filters: Human, and Swiss-Prot reviewed
    query = f'(gene_exact:{gene_name})+(organism_id:9606+AND+reviewed:true)'
    # Output format option is set to JSON with included isoform information
    url = (f"https://rest.uniprot.org/uniprotkb/search?query={query}"
           f"&format=json&fields=accession,gene_names,protein_name&includeIsoform=true")
    response = requests.get(url)

    return response.json().get("results", [])

def parse(entry):
    ''' Parses the entry retrieved from search function [2].'''
    all_names = []
    gene_data = entry.get('genes', [])
    for gene in gene_data:
        if 'geneName' in gene:
            all_names.append(gene['geneName']['value'])
        for syn in gene.get('synonyms', []):
            all_names.append(syn['value'])

    # Get the recommended protein name
    protein_name = (entry.get('proteinDescription', {}).
                    get('recommendedName', {}).
                    get('fullName', {}).
                    get('value', ''))
    return all_names, protein_name

def json(accession, fasta):
    ''' Fetches the results in JSON format for an interactors
    and also returns the FASTA files for an interactors.'''
    if fasta == False:
        url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
    else:
        url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text if fasta else response.json()
    else:
        print(f"Couldn't fetch {'FASTA' if fasta else 'JSON'} for {accession}")
        return None

all_data = list()
# Cache is to prevent duplicates
cache = dict()
# Keeps record of the non-retrieved interactors
failed_genes = list()

if not os.path.exists("fasta_interactors"):
    os.makedirs("fasta_interactors")
    # Handles the interactor names used in research, and retrieves the information
    for research_name in interactors:
        # Removes isoform information suffices (_ORF1, _ORF2, etc.)
        prefix_interactor = research_name.split('_')[0]
        # Uses it from cache if it is already processed
        if prefix_interactor in cache:
            entries = cache[prefix_interactor]
        # Adds to cache
        else:
            entries = search(prefix_interactor)
            cache[prefix_interactor] = entries

        # Each isoform has different accession numbers (Q15645-1, Q15645-2, etc.)
        accession_numbers = list()
        for entry in entries:
            accession_numbers.append(entry['primaryAccession'])
        accession_numbers.sort()

        # Targets isoforms
        if '_ORF' in research_name:
            # It is failed if it doesn't find multiple accession numbers from the UniProt
            if len(accession_numbers) == 1 or len(accession_numbers) == 0:
                failed_genes.append(research_name)
                print(f"{research_name} failed because there is no isoforms in UniProt")
                continue

            # Extracts the last letter, which indicates the
            # isoform number (1 for _ORF1, 2 for _ORF2, etc.)
            isoform_num = int(research_name[-1])
            # Assigns the suffices for isoforms.
            # Q15645 and Q15645-1 always returns the first isoform.
            if isoform_num == 1:
                accession = accession_numbers[int(isoform_num - 1)] + str(-1)
            else:
                accession = accession_numbers[int(isoform_num - 1)]
        # It is failed if UniProt search didn't return any accession numbers.
        else:
            if len(accession_numbers) == 0:
                failed_genes.append(research_name)
                continue
            accession = accession_numbers[0]
        print(f"Selected UniProt entry for {research_name}: {accession}")

        full_entry = json(accession, False)
        fasta = json(accession, True)

        synonyms, protein_name = parse(full_entry)

        if fasta:
            lines = fasta.strip().splitlines()

            # Checks if first line is a proper FASTA header
            if not lines or not lines[0].startswith(">"):
                print(f"{research_name} doesn't have a proper FASTA header")
                continue

            # Separates the header and the sequence
            header = lines[0]
            # Also, removes line breaks and spaces
            sequence = ''.join(lines[1:]).replace(" ", "")

            # Formats the header
            synonym_str = "_".join(synonyms)
            prot_name_clean = protein_name.replace(" ", "_")
            # New header includes all essential information
            customised_header = (f">{research_name}:"
                                 f"{prefix_interactor}:"
                                 f"{accession}:"
                                 f"Synonyms-{synonym_str}:"
                                 f"{prot_name_clean}")

            # Ensures not to repeat writing the files for each run of the code (debugging).


            # Writes to files for each isoform that are retrieved
            output_path = os.path.join("fasta_interactors", f"{research_name}.fasta")
            print(f"Output path is: {output_path}")
            with open(output_path, "w") as f:
                f.write(customised_header + "\n")
                f.write(sequence + "\n")

        else:
            print(f"{research_name}, no FASTA information")

if failed_genes:
    print(f"Failed genes: {failed_genes}")

def write_isoform_fasta(output_dir):
    '''This function writes fasta files for isoforms
    that are retrieved from the research paper.'''

    # mmc2.xlsx has all isoform fasta information
    path = "mmc2.xlsx"
    sheet_name = "1B-Ref+Alt ORFs"
    df = pd.read_excel(path, sheet_name=sheet_name)

    # Extracts ORF column from the Excel and zip it with Isoform ID.
    isoform_sequence_map = dict(zip(df["Isoform_ID"],
                                    df["Isoform__Open_Reading_Frame_Sequence"]))

    # Ensures not to repeat writing the files for each run of the code (debugging).
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

        # Creates the fasta files for each isoform
        for isoform_name, sequence in isoform_sequence_map.items():
            filename = os.path.join(output_dir, f"{isoform_name}.fasta")
            with open(filename, "w") as fasta_file:
                fasta_file.write(f">{isoform_name}\n{sequence}\n")

write_isoform_fasta("fasta_isoforms")

input_fasta_isoforms = "fasta_isoforms"
output_translated_fasta = "translated_isoforms"

# Ensures not to repeat writing the files for each run of the code (debugging).
if not os.path.exists(output_translated_fasta):
    os.makedirs(output_translated_fasta)

    # Need to extract each frame
    frames = [1, 2, 3, -1, -2, -3]

    # Retrieves all the fasta files and their path
    for fasta_file in os.listdir(input_fasta_isoforms):
        if fasta_file.endswith(".fasta"):
            input_path = os.path.join(input_fasta_isoforms, fasta_file)
            output_path = os.path.join(output_translated_fasta, fasta_file)

            with open(output_path, "w") as outfile:
                for frame in frames:
                    # Translates in each frame with SeqKit [3]
                    cmd = f"seqkit translate --frame {frame} --line-width 0 {input_path}"
                    # Processs the code line and waits until each run finishes (subprocess).
                    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

                    # Tag headers with frame information
                    frame_tag = f"_frame{frame}"
                    for line in result.stdout.splitlines():
                        if line.startswith(">"):
                            base_header = line[1:]
                            new_header = f">{base_header}{frame_tag}"
                            outfile.write(new_header + "\n")
                        else:
                            outfile.write(line + "\n")

output_longest_translated_fasta = "translated_isoforms_longest"

def retrieve_sensible_fasta(sequence):
    ''' Returns the sequence before
    it hits the stop codon.'''
    if '*' in sequence:
        return sequence.split('*')[0]
    return sequence

# Ensures not to repeat writing the files for each run of the code (debugging).
if not os.path.exists(output_longest_translated_fasta):
    os.makedirs(output_longest_translated_fasta)

    for tr_fasta in os.listdir(output_translated_fasta):
        if tr_fasta.endswith(".fasta"):
            # Creates the paths for fasta files
            input_path = os.path.join(output_translated_fasta, tr_fasta)
            output_path = os.path.join(output_longest_translated_fasta, tr_fasta)

            # Iterate the frames into a list [4] [5]
            frames = list(SeqIO.parse(input_path, "fasta"))
            best_frame = None
            longest_seq = -1

            for frame in frames:
                shortened_seq = retrieve_sensible_fasta(str(frame.seq))
                # Iterates over all frames and keeps the longest seq
                if len(shortened_seq) > longest_seq:
                    longest_seq = len(shortened_seq)
                    best_frame = frame
                    best_frame.seq = best_frame.seq[:len(shortened_seq)]

            # Writes the longest frame
            if best_frame:
                with open(output_path, "w") as out_f:
                    writer = FastaWriter(out_f, wrap=None)
                    writer.write_file([best_frame])

# Assigns isoform id, interactor id and interactors with isoform information
# First two rows are excluded because they are nonsense information
iso_id = df_mmc3.iloc[2:, 2]
int_id = df_mmc3.iloc[2:, 6]
int_info = df_mmc3.iloc[2:, 7]

int_id[int_id == "LOC100288797"] = "TMEM239"

# Maps all three information for pairs
map_iso_int_excel = [[iso, interactor, info] for iso, interactor, info in
                     zip(iso_id, int_id, int_info)]

# Creates the dataframe for the pair map
df_pairs = pd.DataFrame(map_iso_int_excel,
                        columns=['Isoform_ID', 'Interactor_ID', 'Interaction_Found'])
print("Interaction info in dataframe:")
print(df_pairs)


# Saves only the interaction data
df_pairs.to_csv("df_strctmn.csv", index=False)

# Gets rid of the missing values
df_pairs = df_pairs.dropna(subset=['Interaction_Found'])

# Adds a Base_Isoform column to make it easy to select
df_pairs['Base_Isoform'] = df_pairs['Isoform_ID'].str.split('_').str[0]

print(
    f"\nDataframe after NAN values and missing interactors "
    f"are removed and Base_Isoform column is added:\n {df_pairs}\n")

# Groups by Base_Isoform name and Interactor_ID to get
# cases where an isoform interacted the same with all
# interactor proteins [6]
conflicting_interactions = df_pairs.groupby(['Base_Isoform', 'Interactor_ID'])[
                               'Interaction_Found'].nunique() > 1

# Filters based on conflicting interactions
df_pairs_filtered = df_pairs[df_pairs.set_index(['Base_Isoform', 'Interactor_ID']).index.isin(
    conflicting_interactions[conflicting_interactions].index)]

print(f"Filtered dataframe:\n {df_pairs_filtered}\n")

# Some statistics:
print(f"\nUnique Isoform Number: {df_pairs_filtered['Isoform_ID'].nunique()}")
print(f"Unique Gene Number: {df_pairs_filtered['Base_Isoform'].nunique()}")
print(f"Unique Interactor Protein Number: {df_pairs_filtered['Interactor_ID'].nunique()}\n")

# Merging and writing the pairs step
interactor_dir = "fasta_interactors"
output_pairs_dir = "fasta_pairs"
isoform_filtered_dir = "fasta_isoforms_filtered"
interactor_filtered_dir = "fasta_interactors_filtered"

# Ensures not to repeat writing the files for each run of the code (debugging).
for path in [output_pairs_dir, isoform_filtered_dir, interactor_filtered_dir]:
    os.makedirs(path, exist_ok=True)

for _, row in df_pairs.iterrows():
    isoform_id = row["Isoform_ID"]
    interactor_id = row["Interactor_ID"]

    # File paths for the output files
    isoform_fasta = os.path.join(output_longest_translated_fasta, f"{isoform_id}.fasta")
    interactor_fasta = os.path.join(interactor_dir, f"{interactor_id}.fasta")
    output_fasta = os.path.join(output_pairs_dir, f"{isoform_id}__{interactor_id}.fasta")
    isoform_out = os.path.join(isoform_filtered_dir, f"{isoform_id}.fasta")
    interactor_out = os.path.join(interactor_filtered_dir, f"{interactor_id}.fasta")

    if os.path.exists(isoform_fasta) and os.path.exists(interactor_fasta):
        with open(output_fasta, "w") as outfile:
            # Writes isoform FASTA
            with open(isoform_fasta, "r") as infile:
                iso_data = infile.read().strip()
                outfile.write(iso_data + "\n")
            # Writes interactor FASTA
            with open(interactor_fasta, "r") as infile:
                interactor_data = infile.read().strip()
                outfile.write(interactor_data + "\n")

        # Also writes individual isoform FASTA files
        if not os.path.exists(isoform_out):
            with open(isoform_out, "w") as f:
                f.write(iso_data + "\n")

        # And, writes individual interactor FASTA files
        if not os.path.exists(interactor_out):
            with open(interactor_out, "w") as f:
                f.write(interactor_data + "\n")

    else:
        print(f"{isoform_id} and {interactor_id} are missing")

all_fasta = list()
for dir in [isoform_filtered_dir, interactor_filtered_dir]:
    for file in os.listdir(dir):
        if file.endswith(".fasta"):
            all_fasta.append(os.path.join(dir, file))

with open("structman_input.fasta", "w") as fw:
    for path in all_fasta:
        with open(path, "r") as fr:
            fw.write(fr.read())

structman_output_dir = "structman_output"
def swap_colon(pdb):
    '''Swaps chains in Structure Recommendation
    column in .protein_protein_interactions'''
    if ":" in pdb:
        l, r = pdb.rsplit(":", 1)
        if l and r:
            l, r = l[:-1] + r[0], l[-1] + r[1:]
        return l + ":" + r
    return pdb

# Needs the output file
if os.path.exists(structman_output_dir):
    temp = pd.read_csv(
        f'{structman_output_dir}/structman_input.protein_protein_interactions'
                                '.tsv', sep='\t')

    # Selects only the relevant columns
    df_strctmn = temp[['Input Protein ID A', 'Input Protein ID B', 'Structure Recommendation',
                            'Interaction Score']].copy()

    # Changes the column names
    df_strctmn.columns= ['Isoform_ID', 'Interactor_ID', 'Structure_Recommendation',
                               'Interaction_Score']

    # Mask that looks for cloned isoforms on interactor column [7]
    swap_mask = df_strctmn['Interactor_ID'].str.contains('frame')

    # Swaps the values when it finds an isoform on interactor column
    df_strctmn.loc[swap_mask, ['Isoform_ID', 'Interactor_ID']] = (
        df_strctmn.loc[swap_mask, ['Interactor_ID', 'Isoform_ID']].values)

    # Change the Structure Recommendation structure (7ANK:B:C would be changed to
    # 7ANK:C:B if the column were swapped.)
    df_strctmn.loc[swap_mask, 'Structure_Recommendation'] = \
        df_strctmn.loc[swap_mask, 'Structure_Recommendation'].apply(swap_colon)

    # Gets rid of the _frame* suffices [8]
    df_strctmn['Isoform_ID'] = df_strctmn['Isoform_ID'].astype(str).str.replace(
        r'_frame\d+$', '', regex=True)
    df_strctmn['Interactor_ID'] = df_strctmn['Interactor_ID'].astype(str).str.replace(
        r'_frame\d+$', '', regex=True)

    # Filters by the selected df_pairs [9]
    filter_pairs = set(
        df_pairs_filtered.apply(lambda row: frozenset([row['Isoform_ID'], row['Interactor_ID']]),
                                axis=1))
    ppi_filtered = df_strctmn[df_strctmn.apply(lambda row: frozenset(
        [row['Isoform_ID'], row['Interactor_ID']]) in filter_pairs,
                                                 axis=1)]

    # Merges the Y2H information
    sorted_ppi = ppi_filtered.merge(
        df_pairs_filtered[['Isoform_ID', 'Interactor_ID', 'Interaction_Found']],
        on=['Isoform_ID', 'Interactor_ID'],
        how='left'
    )

    # Reformats the order of the columns
    ppi_condensed = sorted_ppi[["Isoform_ID", "Interactor_ID", "Interaction_Found",
                                "Interaction_Score",
                                "Structure_Recommendation"]]

    ppi_condensed = ppi_condensed.sort_values(by='Isoform_ID')

    ppi_condensed = ppi_condensed.drop_duplicates(subset=['Isoform_ID', 'Interactor_ID',
                                                    'Interaction_Score',
                                                    'Structure_Recommendation'])

    isoform_id = ppi_condensed['Isoform_ID'].tolist()

    # Gets the prefixes
    prefix_iso = [x[:-2] for x in isoform_id]

    # Removes incomparable ones
    c = Counter(prefix_iso)
    to_remove = {x for x, count in c.items() if count == 1}
    final_df = ppi_condensed.copy()
    for i, iso in ppi_condensed['Isoform_ID'].items():
        if iso[:-2] in to_remove:
            final_df = final_df.drop(index = i)

    print(f"Final result:\n{final_df}")
    final_df.to_csv(f"{structman_output_dir}/structman_result.tsv", sep="\t", index=False)





'''
######### REFERENCES ######### 
1 - https://www.uniprot.org/help/api_queries
2 - https://requests.readthedocs.io/en/latest/user/quickstart/
3 - https://bioinf.shenwei.me/seqkit/usage/
4 - https://stackoverflow.com/questions/66923711/how-to-parse-fasta-string-using-seq-io-from-biopython
5 - https://biopython.org/wiki/SeqIO
6 - https://stackoverflow.com/questions/72651199/finding-non-unique-rows-in-pandas-dataframe
7 - https://stackoverflow.com/questions/69663149/how-to-change-values-of-masked-column-in-pandas
8 - https://stackoverflow.com/questions/70079432/remove-suffix-if-string-matches-regular-expression-in-pandas
9 - https://stackoverflow.com/questions/55425324/pandas-drop-duplicates-based-on-2-columns-sometimes-reversed
'''
