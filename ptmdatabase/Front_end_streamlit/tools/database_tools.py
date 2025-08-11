from Bio import SeqIO
import pandas as pd
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


def parse_matrix_file(file):
    if file.name.endswith('.xlsx'):
        df = pd.read_excel(file)
    elif file.name.endswith('.tsv'):
        df = pd.read_csv(file, sep='\t')
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

def extract_modifications_multi(peptide, ptm_types):
    """
    Extract modifications for a given peptide based on multiple PTM types.
    Returns a tuple of (clean_peptide, modifications) where modifications is a list
    of tuples (ptm_type, residue, relative_position, formatted_annotation).
    """
    peptide = re.sub(r'(?:n)\[\d+(?:\.\d+)?\]', '', peptide)
    modifications = []
    clean_peptide = ""
    i = 0

    # Glyco pattern
    glyco_pattern = re.compile(r'N\d+H\d+F\d+S\d+G\d+')

    # Phospho numeric annotations to recognize (string form inside brackets)
    phospho_numeric = r'(?:79|181|167|243)(?:\.\d+)?$'

    while i < len(peptide):
        ch = peptide[i]

        # ----- Bracketed annotation branch -----
        if ch == '[':
            end = peptide.find(']', i)
            if end != -1:
                mod_annotation = peptide[i+1:end]
                if not clean_peptide:
                    # No residue yet to attach this mod to -> skip this block safely
                    i = end + 1
                    continue

                mod_residue = clean_peptide[-1]
                relative_position = len(clean_peptide) - 1

                # Check each requested PTM type
                for ptm in ptm_types:
                    if ptm == 'Phosphorylation' and mod_residue in "STY":
                        # Accept 'P' or numeric masses 79/181/167/243 (optionally with decimals)
                        if (mod_annotation == 'P') or re.match(r'^' + phospho_numeric, mod_annotation):
                            modifications.append((
                                ptm,
                                mod_residue,
                                relative_position,
                                f"{mod_residue}{relative_position+1}[phospho]"
                            ))

                    elif ptm == 'Acetylation':
                        if (mod_annotation == 'A') or re.match(r'^42(?:\.\d+)?$', mod_annotation):
                            modifications.append((
                                ptm,
                                mod_residue,
                                relative_position,
                                f"{mod_residue}{relative_position+1}[ac]"
                            ))

                    elif ptm == 'Ubiquitination':
                        if (mod_annotation == 'U') or re.match(r'^114(?:\.\d+)?$', mod_annotation):
                            modifications.append((
                                ptm,
                                mod_residue,
                                relative_position,
                                f"{mod_residue}{relative_position+1}[ub]"
                            ))

                    elif ptm in ['N-linked Glycosylation', 'O-linked Glycosylation']:
                        if glyco_pattern.match(mod_annotation):
                            if ptm == 'N-linked Glycosylation':
                                modifications.append((
                                    ptm,
                                    mod_residue,
                                    relative_position,
                                    f"N{relative_position+1}[{mod_annotation}]"
                                ))
                            else:
                                modifications.append((
                                    ptm,
                                    mod_residue,
                                    relative_position,
                                    f"{mod_residue}{relative_position+1}[{mod_annotation}]"
                                ))
                i = end + 1
            else:
                # unmatched '[', treat it as a literal char
                clean_peptide += ch
                i += 1

        else:
            # ----- Non-bracket branch -----
            # Lowercase s/t/y are considered phospho (if requested)
            # ----- Non-bracket branch -----
            for ptm in ptm_types:
                if ptm == 'Phosphorylation' and ch in 'sty':
                    # normalize to uppercase in clean peptide
                    upper = ch.upper()
                    relative_position = len(clean_peptide)
                    clean_peptide += upper
                    modifications.append((
                        'Phosphorylation',
                        upper,
                        relative_position,
                        f"{upper}{relative_position+1}[phospho]"
                    ))
                    break
            else:
                # if no PTM condition matched, just add the character as-is
                clean_peptide += ch.upper() if ch.isalpha() else ch
            i += 1

    return clean_peptide, modifications

def process_modifications_multi(peptide_sequence, protein_id, uniprot_sequences, modifications):
    """
    Processes a peptide's modifications by mapping them to the protein sequence.
    Groups modifications by protein site so that modifications on the same residue 
    (even from different PTM types) are combined.
    Returns a tuple (new_header, annotated_protein_sequence).
    """
    protein_data = uniprot_sequences[protein_id]
    protein_sequence = protein_data['sequence']
    peptide_start = protein_sequence.find(peptide_sequence)
    if peptide_start == -1:
        raise ValueError(f"Peptide sequence not found in protein {protein_id}.")

    # Group modifications by site position (1-based index)
    site_mods = {}
    for mod in modifications:
        ptm_type, mod_residue, relative_position, mod_annotation = mod
        site_position = peptide_start + relative_position + 1
        # Extract just the inner content (e.g., "phospho" or "acetyl")
        mod_content = mod_annotation[mod_annotation.find('[') + 1:-1]
        site_mods.setdefault(site_position, []).append(mod_content)

    combined_mod_descriptions = []
    modified_protein_sequence = list(protein_sequence)

    # Combine modifications for each site
    for site_position in sorted(site_mods.keys()):
        residue = protein_sequence[site_position - 1]
        # Combine the modification contents (remove duplicates)
        mods_combined = ",".join(sorted(set(site_mods[site_position])))
        combined_annotation = f"{residue}{site_position}[{mods_combined}]"
        combined_mod_descriptions.append(combined_annotation)
        annotation_tag = f"[{mods_combined}]"
        if annotation_tag not in modified_protein_sequence[site_position - 1]:
            modified_protein_sequence[site_position - 1] += annotation_tag

    # Build header by joining separate residue annotations with underscores
    mod_description_str = '_'.join(combined_mod_descriptions)
    new_header = f"sp|{protein_id}|Mod:{mod_description_str}|{protein_data['header'].split('|', 2)[2]}"
    return new_header, ''.join(modified_protein_sequence)

def generate_ptm_entries_multi(peptide_list, uniprot_sequences, ptm_types):
    """
    Processes a list of peptides with multiple PTM types.
    Ensures that modifications for each peptide are extracted together.
    Returns (ptm_entries, missing_peptides, inferred_protein_ids).
    """
    ptm_entries = []
    missing_peptides = []
    inferred_protein_ids = set()
    peptide_to_proteins = {}
    protein_to_peptides = {}

    for peptide in peptide_list:
        # Extract all modifications for the given PTM types at once:
        peptide_sequence, modifications = extract_modifications_multi(peptide, ptm_types)
        if not modifications:
            continue

        found_protein = False
        potential_proteins = []
        for protein_id, protein_data in uniprot_sequences.items():
            if protein_data['sequence'].find(peptide_sequence) != -1:
                found_protein = True
                potential_proteins.append(protein_id)
                protein_to_peptides.setdefault(protein_id, []).append(peptide_sequence)

        if found_protein:
            peptide_to_proteins[peptide_sequence] = {
                'proteins': potential_proteins,
                'modifications': modifications
            }
        else:
            missing_peptides.append(peptide)

    # Assign unique peptides first
    unique_peptides = {p: data for p, data in peptide_to_proteins.items() if len(data['proteins']) == 1}
    for peptide_sequence, data in unique_peptides.items():
        protein_id = data['proteins'][0]
        inferred_protein_ids.add(protein_id)
        header, annotated_seq = process_modifications_multi(peptide_sequence, protein_id, uniprot_sequences, data['modifications'])
        ptm_entries.append((header, annotated_seq))
        peptide_to_proteins.pop(peptide_sequence)

    # Greedily assign remaining shared peptides
    while peptide_to_proteins:
        best_protein = max(protein_to_peptides,
                           key=lambda p: len(set(protein_to_peptides[p]) & set(peptide_to_proteins.keys())))
        inferred_protein_ids.add(best_protein)
        for peptide_sequence in protein_to_peptides[best_protein]:
            if peptide_sequence in peptide_to_proteins:
                data = peptide_to_proteins[peptide_sequence]
                header, annotated_seq = process_modifications_multi(peptide_sequence, best_protein, uniprot_sequences, data['modifications'])
                ptm_entries.append((header, annotated_seq))
                peptide_to_proteins.pop(peptide_sequence)

    return ptm_entries, missing_peptides, inferred_protein_ids

def write_fasta(output_file, uniprot_sequences, ptm_entries, inferred_protein_ids, include_global_protein_entries=False):
    written_entries = set()
    write_count = 0
    
    for header, sequence in ptm_entries:
        formatted_sequence = format_fasta_sequence(sequence)
        entry = (header, formatted_sequence)
        if entry not in written_entries:
            output_file.write(f">{header}\n{formatted_sequence}\n")  # Writing to StringIO
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
                    output_file.write(f">{header}\n{sequence}\n")  # Writing to StringIO
                    written_entries.add(entry)
                    write_count += 1

    print(f"Total unique entries written: {write_count}")

def write_missing_info(output_file, missing_peptides):
    # Convert the missing peptides list into a DataFrame and remove duplicates
    missing_peptides_df = pd.DataFrame(missing_peptides, columns=['Peptide Sequence']).drop_duplicates()
    output_file.write(missing_peptides_df.to_csv(index=False))  

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
