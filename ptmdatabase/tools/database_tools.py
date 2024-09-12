from Bio import SeqIO
import os
import pandas as pd
import concurrent.futures
import re

# Database library: 5 original PTMs databases (Phospho, N- and O-linked Glyco, Acetylation, and Ubiquitination) were generated using 
# 1. The PTM text file from (https://awi.cuhk.edu.cn/dbPTM/download.php).
# 2. The uniprot database (https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz). 
# 3. While using the code in PTM_Database\ptmdatabase\Database_library\Original_database_generation.ipynb -> add the PTM information in the entry and add annotation for the PTM site in the protein sequence.

# Workflow:
# 1. The code extract the protein ID and peptide sequence from the data file. 
# 2. It will used the protein ID to find the corresponding global protein entries in the uniprot fasta file and pasted those entries to the generated database.
# 3. It will then used both the protein ID and the peptide sequence to determine the exact PTM sites in the corresponding global protein sequence.
# 4. After obtaining a list of the PTM sites, the code will match the protein IDs and the PTM sites to the corresponding entries in the specific original PTM database, which can be found in PTM_Database\ptmdatabase\Database_library.
# 5. The code will then proceed to paste those entries to the generated database, which already contained the global protein entries in there.
# 6. For the PTM sites that do not exist in the original PTM databases, the code will automatically create new entries for those PTM sites using the corresponding Global protein entries.
# 7. All of the unmatched protein ID (Proteins that are listed in the matrix file but cannot be found in the UniProt database), peptide sequence (Peptides that are identified in the matrix file but cannot be found within the corresponding protein sequence in the UniProt database), and PTM sites (Modifications that are identified in the matrix file but cannot be found in the PTM-specific library) are recorded in the Excel list located in the same directory of the generated database. 


def parse_matrix_file(file_path):
    if file_path.endswith('.xlsx'):
        df = pd.read_excel(file_path)
    elif file_path.endswith('.tsv'):
        df = pd.read_csv(file_path, sep='\t')
    else:
        raise ValueError("Unsupported file format. Only .xlsx and .tsv are supported.")
    return df

def format_fasta_sequence(sequence, line_length=60):
    return '\n'.join([sequence[i:i+line_length] for i in range(0, len(sequence), line_length)])

def load_uniprot_sequences(fasta_file):
    uniprot_sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        protein_id = record.id.split('|')[1]
        description = record.description
        sequence = str(record.seq)
        uniprot_sequences[protein_id] = {'header': description, 'sequence': sequence}
    return uniprot_sequences

def load_ptm_sequences(fasta_file):
    ptm_sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        description = record.description
        sequence = str(record.seq)
        key = '|'.join(description.split('|')[:3]) + '|'
        ptm_sequences[key] = {'header': description, 'sequence': sequence}
    return ptm_sequences

# Global processing 
def extract_modifications(peptide, ptm_type):
    modifications = []
    clean_peptide = ""
    i = 0
    while i < len(peptide):
        if peptide[i] == '[':
            end = peptide.find(']', i)
            if end != -1:
                mod_annotation = peptide[i+1:end]
                mod_residue = clean_peptide[-1]
                relative_position = len(clean_peptide) - 1
                if ptm_type == 'Phosphorylation' and mod_residue in "STY":
                    if mod_annotation == 'P' or re.match(r'79(\.\d+)?', mod_annotation):
                        modifications.append((mod_residue, f"{mod_residue}{relative_position + 1}P", relative_position))
                elif ptm_type == 'Acetylation' and mod_residue == 'K':
                    if mod_annotation == 'A' or re.match(r'42(\.\d+)?', mod_annotation):
                        modifications.append((mod_residue, f"{mod_residue}{relative_position + 1}A", relative_position))
                elif ptm_type == 'Ubiquitination' and mod_residue == 'K':
                    if mod_annotation == 'U' or re.match(r'114(\.\d+)?', mod_annotation):
                        modifications.append((mod_residue, f"{mod_residue}{relative_position + 1}U", relative_position))
                i = end + 1
            else:
                clean_peptide += peptide[i]
                i += 1
        else:
            clean_peptide += peptide[i]
            i += 1
    return clean_peptide, modifications

def generate_ptm_entries(peptide_list, uniprot_sequences, ptm_type):
    ptm_entries = []
    missing_peptides = []
    inferred_protein_ids = set()

    # Step 1: Build peptide-to-protein mapping and protein-to-peptide mapping
    peptide_to_proteins = {}
    protein_to_peptides = {}

    for peptide in peptide_list:
        peptide_sequence, modifications = extract_modifications(peptide, ptm_type)
        if not modifications:  # Skip if no modifications are found
            continue

        found_protein = False
        potential_proteins = []

        for protein_id, protein_data in uniprot_sequences.items():
            protein_sequence = protein_data['sequence']
            peptide_start = protein_sequence.find(peptide_sequence)

            if peptide_start != -1:
                found_protein = True
                potential_proteins.append(protein_id)

                # Track which peptides are covered by this protein
                if protein_id not in protein_to_peptides:
                    protein_to_peptides[protein_id] = []
                protein_to_peptides[protein_id].append(peptide_sequence)

        if found_protein:
            peptide_to_proteins[peptide_sequence] = potential_proteins
        else:
            missing_peptides.append(peptide)

    # Step 2: Assign unique peptides to their corresponding proteins
    unique_peptides = {p: ps[0] for p, ps in peptide_to_proteins.items() if len(ps) == 1}

    for peptide, protein_id in unique_peptides.items():
        inferred_protein_ids.add(protein_id)

        # Generate the PTM entry for unique peptides
        mod_descriptions = process_modifications(peptide, protein_id, uniprot_sequences, modifications)
        ptm_entries.append(mod_descriptions)

        # Remove unique peptides from further processing
        peptide_to_proteins.pop(peptide)

    # Step 3: Greedily assign shared peptides
    while peptide_to_proteins:
        # Find the protein that covers the most unassigned peptides
        best_protein = max(protein_to_peptides, key=lambda p: len(set(protein_to_peptides[p]) & set(peptide_to_proteins.keys())))
        inferred_protein_ids.add(best_protein)

        # Assign all peptides covered by this protein
        for peptide in protein_to_peptides[best_protein]:
            if peptide in peptide_to_proteins:
                mod_descriptions = process_modifications(peptide, best_protein, uniprot_sequences, modifications)
                ptm_entries.append(mod_descriptions)
                peptide_to_proteins.pop(peptide)

    return ptm_entries, missing_peptides, inferred_protein_ids


# Helper function to process modifications
def process_modifications(peptide_sequence, protein_id, uniprot_sequences, modifications):
    protein_data = uniprot_sequences[protein_id]
    protein_sequence = protein_data['sequence']
    peptide_start = protein_sequence.find(peptide_sequence)

    mod_descriptions = []
    
    # Loop through modifications and adjust them relative to the protein sequence
    for mod in modifications:
        if len(mod) == 3:  # Ensure it follows the (residue, description, position) format
            mod_residue, mod_desc, relative_position = mod
            site_position = peptide_start + relative_position + 1
            mod_descriptions.append(f"{mod_residue}{site_position}P")
        else:
            raise ValueError(f"Modification format is incorrect: {mod}. Expected (residue, description, position).")

    # Construct the new header
    mod_description = '_'.join(mod_descriptions)
    new_header = f"sp|{protein_id}|{mod_description}|{protein_data['header'].split('|', 2)[2]}"

    # Annotate the protein sequence with modification annotations
    modified_protein_sequence = list(protein_sequence)
    for mod in modifications:
        if len(mod) == 3:
            mod_residue, mod_desc, relative_position = mod
            site_position = peptide_start + relative_position + 1
            modified_protein_sequence[site_position - 1] += f"[{mod_desc[-1]}]"

    modified_protein_sequence = ''.join(modified_protein_sequence)
    return (new_header, modified_protein_sequence)


# Glyco processing
def extract_glyco_modifications(peptide):
    modifications = []
    clean_peptide = ""
    i = 0
    glyco_pattern = re.compile(r'([HNFSG]\d+)+')
    while i < len(peptide):
        if peptide[i] == '[':
            # Find the closing parenthesis
            end = peptide.find(']', i)
            if end != -1:
                mod_annotation = peptide[i+1:end]
                mod_residue = clean_peptide[-1]
                relative_position = len(clean_peptide) - 1
                if glyco_pattern.fullmatch(mod_annotation):
                    modifications.append((mod_residue, mod_annotation, relative_position))
                i = end + 1
            else:
                clean_peptide += peptide[i]
                i += 1
        else:
            clean_peptide += peptide[i]
            i += 1
    return clean_peptide, modifications

def generate_ptm_entries_glyco(peptide_list, uniprot_sequences, ptm_type):
    ptm_entries = []
    missing_peptides = []
    inferred_protein_ids = set()

    for peptide in peptide_list:
        peptide_sequence, modifications = extract_glyco_modifications(peptide)
        found_protein = False

        if not modifications:
            continue  # Skip if there are no modifications in the peptide

        for protein_id, protein_data in uniprot_sequences.items():
            protein_sequence = protein_data['sequence']
            peptide_start = protein_sequence.find(peptide_sequence)
            if peptide_start != -1:
                found_protein = True
                inferred_protein_ids.add(protein_id)

                for mod in modifications:
                    mod_residue, mod_annotation, relative_position = mod
                    # Calculate the protein-level position
                    site_position = peptide_start + relative_position + 1
                    
                    if ptm_type == 'N-linked Glycosylation':
                        mod_description = f"N{site_position}[{mod_annotation}]"
                    elif ptm_type == 'O-linked Glycosylation':
                        mod_description = f"{mod_residue}{site_position}[{mod_annotation}]"

                    # Update header and sequence
                    new_header = f"sp|{protein_id}|{mod_description}|{protein_data['header'].split('|', 2)[2]}"
                    modified_protein_sequence = list(protein_sequence)
                    modified_protein_sequence[site_position - 1] += f"[{mod_annotation}]"
                    modified_protein_sequence = ''.join(modified_protein_sequence)
                    ptm_entries.append((new_header, modified_protein_sequence))

                break

        if not found_protein:
            missing_peptides.append(('', peptide))

    return ptm_entries, missing_peptides, inferred_protein_ids

def write_fasta(output_file, uniprot_sequences, ptm_entries, inferred_protein_ids, include_global_protein_entries=False):
    with open(output_file, 'w') as file:
        written_entries = set()
        write_count = 0
        
        for header, sequence in ptm_entries:
            formatted_sequence = format_fasta_sequence(sequence)
            entry = (header, formatted_sequence)
            if entry not in written_entries:
                file.write(f">{header}\n{formatted_sequence}\n")
                written_entries.add(entry)
                write_count += 1
        
        if include_global_protein_entries:
            for protein_id in inferred_protein_ids:
                if protein_id in uniprot_sequences:
                    data = uniprot_sequences[protein_id]
                    header = data['header']
                    sequence = format_fasta_sequence(data['sequence'])
                    entry = (header, sequence)
                    if entry not in written_entries:
                        file.write(f">{header}\n{sequence}\n")
                        written_entries.add(entry)
                        write_count += 1

        print(f"Total unique entries written: {write_count}")

def write_missing_info(output_file_dir, missing_peptides):
    # Convert the missing peptides list into a DataFrame and remove duplicates
    missing_peptides_df = pd.DataFrame(missing_peptides, columns=['Peptide Sequence']).drop_duplicates()

    # Set the output file path with the new name
    output_file = os.path.join(output_file_dir, 'missing_peptides.xlsx')

    # Write the missing peptides to an Excel file
    with pd.ExcelWriter(output_file, engine='xlsxwriter') as writer:
        missing_peptides_df.to_excel(writer, sheet_name='Missing Peptides', index=False)

    print(f"Missing peptides have been recorded in {output_file}")

def count_entries_in_fasta(fasta_file):
    entries = set()
    protein_ids = set()
    for record in SeqIO.parse(fasta_file, "fasta"):
        header = record.description
        sequence = str(record.seq)
        entries.add((header, sequence))
        protein_id = record.id.split('|')[1]
        protein_ids.add(protein_id)
    return len(entries), len(protein_ids)


# MsPycloud code:
# def generate_ptm_entries(df, uniprot_sequences, ptm_sequences, ptm_type):
#     ptm_entries = []
#     missing_ptms = []
#     missing_proteins = []
#     missing_peptides = []

#     for index, row in df.iterrows():
#         accession_entries = row['Protein.Group.Accessions'].split(';')
#         found_protein = False
#         for accession in accession_entries:
#             protein_id = accession.split('|')[1] if '|' in accession else accession
#             if protein_id in uniprot_sequences:
#                 found_protein = True
#                 protein_data = uniprot_sequences[protein_id]
#                 protein_sequence = protein_data['sequence']
#                 peptide_sequence = row['Sequence']
#                 mod_str = row['Modifications']

#                 if ptm_type == 'Phosphorylation':
#                     peptide_start = protein_sequence.find(peptide_sequence)
#                     if peptide_start == -1:
#                         missing_peptides.append((protein_id, peptide_sequence))
#                         continue

#                     modifications = extract_modifications(peptide_sequence, mod_str, ptm_type)
#                     if not modifications:
#                         continue

#                     for mod in modifications:
#                         site_position = peptide_start + peptide_sequence.find(mod[0]) + 1
#                         mod_description = f"{mod[0]}{site_position}P"
#                         ptm_header = f"sp|{protein_id}|{mod_description}|"
#                         if ptm_header in ptm_sequences:
#                             ptm_data = ptm_sequences[ptm_header]
#                             ptm_entries.append((ptm_data['header'], ptm_data['sequence']))
#                         else:
#                             missing_ptms.append((protein_id, site_position, mod[0]))

#                 elif ptm_type in ['N-linked Glycosylation', 'O-linked Glycosylation']:
#                     clean_peptide = peptide_sequence.split('-')[0]
#                     glyco_sites = identify_glycosylation_sites(clean_peptide, ptm_type)

#                     for site in glyco_sites:
#                         glyco_position = protein_sequence.find(clean_peptide) + clean_peptide.find(site[0]) + 1
#                         if ptm_type == 'N-linked Glycosylation':
#                             glyco_description = f"N{glyco_position}nG"
#                         elif ptm_type == 'O-linked Glycosylation':
#                             glyco_description = f"{site[0]}{glyco_position}oG"
                        
#                         ptm_header = f"sp|{protein_id}|{glyco_description}|"
#                         if ptm_header in ptm_sequences:
#                             ptm_data = ptm_sequences[ptm_header]
#                             ptm_entries.append((ptm_data['header'], ptm_data['sequence']))
#                         else:
#                             missing_ptms.append((protein_id, glyco_position, site[0]))

#             else:
#                 missing_proteins.append(protein_id)
#         if not found_protein:
#             missing_proteins.append(accession.split('|')[1] if '|' in accession else accession)
    
#     new_ptm_entries = create_missing_ptm_entries(missing_ptms, uniprot_sequences, ptm_type)
#     ptm_entries.extend(new_ptm_entries)

#     return ptm_entries, missing_ptms, missing_proteins, missing_peptides

# def extract_modifications(peptide, mod_str, ptm_type):
#     modifications = []
#     current_position = 0
#     while '[' in mod_str:
#         start = mod_str.find('[')
#         end = mod_str.find(']')
#         mod_residue = mod_str[start-1]
#         mod_position = current_position + start

#         if ptm_type == 'Phosphorylation' and mod_residue in "STY":
#             mod_description = f"{mod_residue}{mod_position + 1}P"
#             modifications.append(mod_description)
#         elif ptm_type == 'N-linked Glycosylation' and mod_residue == 'N':
#             mod_description = f"{mod_residue}{mod_position + 1}nG"
#             modifications.append(mod_description)
#         elif ptm_type == 'O-linked Glycosylation' and mod_residue in "ST":
#             mod_description = f"{mod_residue}{mod_position + 1}oG"
#             modifications.append(mod_description)

#         mod_str = mod_str[end+1:]
#         current_position += start
#     return modifications

# def identify_glycosylation_sites(peptide, ptm_type):
#     sites = []
#     if ptm_type == 'N-linked Glycosylation':
#         for i in range(len(peptide) - 2):
#             if peptide[i] == 'N' and peptide[i+2] in 'ST' and peptide[i+1] != 'P':
#                 sites.append((peptide[i], i))
#     elif ptm_type == 'O-linked Glycosylation':
#         for i in range(len(peptide)):
#             if peptide[i] in 'ST':
#                 sites.append((peptide[i], i))
#     return sites

# def create_missing_ptm_entries(missing_ptms, uniprot_sequences, ptm_type):
#     new_ptm_entries = []

#     for protein_id, site_position, mod_residue in missing_ptms:
#         if protein_id in uniprot_sequences:
#             protein_data = uniprot_sequences[protein_id]
#             protein_sequence = protein_data['sequence']
#             if site_position <= len(protein_sequence):
#                 modified_protein_sequence = list(protein_sequence)
#                 mod_position_in_protein = site_position - 1

#                 if ptm_type == 'Phosphorylation':
#                     modification = "[P]"
#                     mod_annotation = f"{mod_residue}{site_position}P"
#                 elif ptm_type == 'Acetylation':
#                     modification = "[A]"
#                     mod_annotation = f"N{site_position}A"
#                 elif ptm_type == 'Ubiquitination':
#                     modification = "[U]"
#                     mod_annotation = f"{mod_residue}{site_position}U"

#                 modified_protein_sequence[mod_position_in_protein] += modification
#                 modified_protein_sequence = ''.join(modified_protein_sequence)

#                 original_header = protein_data['header']
#                 parts = original_header.split('|')
#                 new_header = f"sp|{parts[1]}|{mod_annotation}|{parts[2]}"

#                 new_ptm_entries.append((new_header, modified_protein_sequence))

#     return new_ptm_entries

# def generate_all_ptm_entries(protein_ids, ptm_sequences):
#     additional_ptm_entries = []
#     additional_info = []
#     total_proteins = len(protein_ids)

#     def process_protein(protein_id):
#         temp_entries = []
#         temp_info = []
#         for ptm_type, sequences in ptm_sequences.items():
#             for ptm_header, ptm_data in sequences.items():
#                 if protein_id in ptm_header:
#                     temp_entries.append((ptm_data['header'], ptm_data['sequence']))
#                     temp_info.append((protein_id, ptm_data['header'], ptm_type))
#         return temp_entries, temp_info

#     with concurrent.futures.ThreadPoolExecutor() as executor:
#         futures = {executor.submit(process_protein, protein_id): protein_id for protein_id in protein_ids}
#         for future in concurrent.futures.as_completed(futures):
#             ptm_entries, ptm_info = future.result()
#             additional_ptm_entries.extend(ptm_entries)
#             additional_info.extend(ptm_info)

#     return additional_ptm_entries, additional_info