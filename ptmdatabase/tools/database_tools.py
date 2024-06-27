import re
from Bio import SeqIO
import os
import pandas as pd
import streamlit as st

#Workflow: 5 PTMs database 

def parse_matrix_file(matrix_file):
    return pd.read_csv(matrix_file, delimiter='\t')

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

def extract_modifications(peptide, mod_str, ptm_type):
    modifications = []
    current_position = 0
    while '[' in mod_str:
        start = mod_str.find('[')
        end = mod_str.find(']')
        mod_residue = mod_str[start-1]
        mod_position = current_position + start

        if ptm_type == 'Phosphorylation' and mod_residue in "STY":
            mod_description = f"{mod_residue}{mod_position + 1}P"
            modifications.append(mod_description)
        elif ptm_type == 'N-linked Glycosylation' and mod_residue == 'N':
            mod_description = f"{mod_residue}{mod_position + 1}nG"
            modifications.append(mod_description)

        mod_str = mod_str[end+1:]
        current_position += start
    return modifications

def identify_glycosylation_sites(peptide):
    sites = []
    for i in range(len(peptide) - 2):
        if peptide[i] == 'N' and peptide[i+2] in 'ST' and peptide[i+1] != 'P':
            sites.append((peptide[i], i))
    return sites

def generate_ptm_entries(df, uniprot_sequences, ptm_sequences, ptm_type):
    ptm_entries = []
    missing_ptms = []
    missing_proteins = []
    missing_peptides = []

    for index, row in df.iterrows():
        accession_entries = row['Protein.Group.Accessions'].split(';')
        found_protein = False
        for accession in accession_entries:
            protein_id = accession.split('|')[1] if '|' in accession else accession
            if protein_id in uniprot_sequences:
                found_protein = True
                protein_data = uniprot_sequences[protein_id]
                protein_sequence = protein_data['sequence']
                peptide_sequence = row['Sequence']
                mod_str = row['Modifications']

                if ptm_type == 'Phosphorylation':
                    peptide_start = protein_sequence.find(peptide_sequence)
                    if peptide_start == -1:
                        missing_peptides.append((protein_id, peptide_sequence))
                        continue

                    modifications = extract_modifications(peptide_sequence, mod_str, ptm_type)
                    if not modifications:
                        continue

                    for mod in modifications:
                        site_position = peptide_start + peptide_sequence.find(mod[0]) + 1
                        mod_description = f"{mod[0]}{site_position}P"
                        ptm_header = f"sp|{protein_id}|{mod_description}|"
                        if ptm_header in ptm_sequences:
                            ptm_data = ptm_sequences[ptm_header]
                            ptm_entries.append((ptm_data['header'], ptm_data['sequence']))
                        else:
                            missing_ptms.append((protein_id, site_position, mod[0]))

                elif ptm_type == 'N-linked Glycosylation':
                    clean_peptide = peptide_sequence.split('-')[0]
                    glyco_sites = identify_glycosylation_sites(clean_peptide)

                    for site in glyco_sites:
                        glyco_position = protein_sequence.find(clean_peptide) + clean_peptide.find(site[0]) + 1
                        glyco_description = f"N{glyco_position}nG"
                        ptm_header = f"sp|{protein_id}|{glyco_description}|"
                        if ptm_header in ptm_sequences:
                            ptm_data = ptm_sequences[ptm_header]
                            ptm_entries.append((ptm_data['header'], ptm_data['sequence']))
                        else:
                            missing_ptms.append((protein_id, glyco_position, 'N'))

            else:
                missing_proteins.append(protein_id)
        if not found_protein:
            missing_proteins.append(accession.split('|')[1] if '|' in accession else accession)
    
    new_ptm_entries = create_missing_ptm_entries(missing_ptms, uniprot_sequences, ptm_type)
    ptm_entries.extend(new_ptm_entries)

    return ptm_entries, missing_ptms, missing_proteins, missing_peptides

def create_missing_ptm_entries(missing_ptms, uniprot_sequences, ptm_type):
    new_ptm_entries = []

    for protein_id, site_position, mod_residue in missing_ptms:
        if protein_id in uniprot_sequences:
            protein_data = uniprot_sequences[protein_id]
            protein_sequence = protein_data['sequence']
            if site_position <= len(protein_sequence):
                modified_protein_sequence = list(protein_sequence)
                mod_position_in_protein = site_position - 1

                if ptm_type == 'Phosphorylation':
                    modification = "(P)"
                    mod_annotation = f"{mod_residue}{site_position}P"
                elif ptm_type == 'N-linked Glycosylation':
                    modification = "(nG)"
                    mod_annotation = f"N{site_position}nG"

                modified_protein_sequence[mod_position_in_protein] += modification
                modified_protein_sequence = ''.join(modified_protein_sequence)

                original_header = protein_data['header']
                parts = original_header.split('|')
                new_header = f"sp|{parts[1]}|{mod_annotation}|{parts[2]}"

                new_ptm_entries.append((new_header, modified_protein_sequence))

    return new_ptm_entries

def write_fasta(output_file, uniprot_sequences, protein_ids, ptm_entries=[]):
    with open(output_file, 'w') as file:
        written_entries = set()
        for protein_id in protein_ids:
            if protein_id in uniprot_sequences:
                data = uniprot_sequences[protein_id]
                header = data['header']
                sequence = format_fasta_sequence(data['sequence'])
                entry = (header, sequence)
                if entry not in written_entries:
                    file.write(f">{header}\n{sequence}\n")
                    written_entries.add(entry)
        
        for header, sequence in ptm_entries:
            formatted_sequence = format_fasta_sequence(sequence)
            entry = (header, formatted_sequence)
            if entry not in written_entries:
                file.write(f">{header}\n{formatted_sequence}\n")
                written_entries.add(entry)

def write_missing_info(output_file_dir, missing_ptms, missing_proteins, missing_peptides):
    missing_ptms_df = pd.DataFrame(missing_ptms, columns=['Protein ID', 'Site Position', 'Modification Residue']).drop_duplicates()
    missing_proteins_df = pd.DataFrame(missing_proteins, columns=['Missing Protein IDs']).drop_duplicates()
    missing_peptides_df = pd.DataFrame(missing_peptides, columns=['Protein ID', 'Peptide Sequence']).drop_duplicates()

    output_file = os.path.join(output_file_dir, 'missing_protein_and_sites.xlsx')
    
    with pd.ExcelWriter(output_file, engine='xlsxwriter') as writer:
        missing_ptms_df.to_excel(writer, sheet_name='Missing PTMs', index=False)
        missing_proteins_df.to_excel(writer, sheet_name='Missing Proteins', index=False)
        missing_peptides_df.to_excel(writer, sheet_name='Missing Peptides', index=False)

def count_entries_in_fasta(fasta_file):
    entries = set()
    protein_ids = set()
    for record in SeqIO.parse(fasta_file, "fasta"):
        entries.add(record.description)
        protein_id = record.id.split('|')[1]
        protein_ids.add(protein_id)
    return len(entries), len(protein_ids)
