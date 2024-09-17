import requests 
import streamlit as st
import pandas as pd
import os
import sys
from Bio import SeqIO
import time
from tqdm import tqdm
from multiprocessing import Pool, cpu_count
from pathlib import Path
from tools.database_tools import (
    parse_matrix_file, 
    load_uniprot_sequences, 
    load_ptm_sequences, 
    generate_ptm_entries, 
    write_fasta, 
    write_missing_info, 
    count_entries_in_fasta,
    generate_ptm_entries_glyco,
)

def send_fasta_to_backend(fasta_file_path, username):
    with open(fasta_file_path, "rb") as file:
        files = {'file': file}
        data = {'username': username}  # Include the email in the form data
        response = requests.post("http://127.0.0.1:8000/upload-fasta/", files=files, data=data)
        return response.json()

def initialize_session_state():
    base_dir = Path(__file__).resolve().parent.parent
    database_dir = base_dir / 'Database_library'

    if 'work_dir' not in st.session_state:
        st.session_state['work_dir'] = ""
    if 'original_fasta_dir' not in st.session_state:
        st.session_state['original_fasta_dir'] = os.path.join(database_dir, 'uniprotkb_proteome_UP000005640_AND_revi_2024_07_23.fasta')
    if 'new_db_dir' not in st.session_state:
        st.session_state['new_db_dir'] = ""
    if 'missing_info_file' not in st.session_state:
        st.session_state['missing_info_file'] = ""

def process_peptide_phosphorylation(args):
    chunk, uniprot_sequences = args
    ptm_entries, missing_peptides, inferred_protein_ids = [], [], set()
    for peptide in chunk:
        entries, peptides, inferred_ids = generate_ptm_entries([peptide], uniprot_sequences, 'Phosphorylation')
        ptm_entries.extend(entries)
        missing_peptides.extend(peptides)
        inferred_protein_ids.update(inferred_ids)
    return ptm_entries, missing_peptides, inferred_protein_ids
    
    return ptm_entries, missing_peptides, inferred_protein_ids

def process_peptide_acetylation(args):
    chunk, uniprot_sequences = args
    ptm_entries, missing_peptides, inferred_protein_ids = [], [], set()
    for peptide in chunk:
        entries, peptides, inferred_ids = generate_ptm_entries([peptide], uniprot_sequences, 'Acetylation')
        ptm_entries.extend(entries)
        missing_peptides.extend(peptides)
        inferred_protein_ids.update(inferred_ids)
    return ptm_entries, missing_peptides, inferred_protein_ids

def process_peptide_ubiquitination(args):
    chunk, uniprot_sequences = args
    ptm_entries, missing_peptides, inferred_protein_ids = [], [], set()
    for peptide in chunk:
        entries, peptides, inferred_ids = generate_ptm_entries([peptide], uniprot_sequences, 'Ubiquitination')
        ptm_entries.extend(entries)
        missing_peptides.extend(peptides)
        inferred_protein_ids.update(inferred_ids)
    return ptm_entries, missing_peptides, inferred_protein_ids


def process_peptide_glycosylation(args):
    chunk, uniprot_sequences, ptm_type = args
    ptm_entries, missing_peptides, inferred_protein_ids = [], [], set()
    for peptide in chunk:
        entries, peptides, inferred_ids = generate_ptm_entries_glyco([peptide], uniprot_sequences, ptm_type)
        ptm_entries.extend(entries)
        missing_peptides.extend(peptides)
        inferred_protein_ids.update(inferred_ids)
    return ptm_entries, missing_peptides, inferred_protein_ids

def chunk_list(lst, n):
    """Divide list lst into n chunks."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def main():

    st.markdown(
        """
        <style>
            .title {
                font-family: "Arial", sans-serif;
                color: #008080;
                font-size: 42px;
                font-weight: bold;
            }
            .subtitle {
                font-family: "Arial", sans-serif;
                color: #800000;
                font-size: 22px;
            }
            .main-header {
                font-family: "Arial", sans-serif; 
                font-size: 28px;
                font-weight: bold;
                text-decoration: underline;
            }
        </style>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("<div class='title'>Database Generation</div>", unsafe_allow_html=True)
    st.write("\n")

    initialize_session_state()

    with st.form(key='database_generation_form', clear_on_submit=False):
        matrix_file = st.file_uploader('Peptide List (xlsx or tsv):', type=['xlsx', 'tsv'])

        new_db_dir = st.text_input('Directory to Store Generated Database:', value=st.session_state['new_db_dir'])
        st.session_state['new_db_dir'] = new_db_dir

        modification_types = st.multiselect(
            'Select PTM Types to Process',
            ['Phosphorylation', 'Acetylation', 'Ubiquitination', 'N-linked Glycosylation', 'O-linked Glycosylation']
        )

        include_global_protein_entries = st.checkbox('Include Global Protein Entries', value=False)

        submit_button = st.form_submit_button(label='Generate Database')

        if submit_button:
            if not matrix_file:
                st.error("No file uploaded. Please upload a valid .xlsx or .tsv file.")
                return

            st.session_state['work_dir'] = matrix_file
            output_file = new_db_dir

            # Ensure the directory exists
            if not Path(output_file).exists():
                st.error(f"The directory '{output_file}' does not exist. Please provide a valid directory.")
                return

            missing_info_file = os.path.dirname(output_file)
            st.session_state['missing_info_file'] = missing_info_file

            try:
                df = parse_matrix_file(matrix_file)
                uniprot_sequences = load_uniprot_sequences(st.session_state['original_fasta_dir'])

                peptide_list = df.iloc[:, 0].tolist()

                ptm_entries = []
                missing_peptides = []
                inferred_protein_ids = set()

                start_time = time.time()
                num_cpus = cpu_count()
                chunked_peptide_list = list(chunk_list(peptide_list, max(1, len(peptide_list) // num_cpus)))      

                with Pool(num_cpus) as pool:
                    if 'Phosphorylation' in modification_types:
                        args = [(chunk, uniprot_sequences) for chunk in chunked_peptide_list]
                        results = list(tqdm(pool.imap(process_peptide_phosphorylation, args), total=len(chunked_peptide_list), desc="Processing Phosphorylation"))
                        for result in results:
                            phospho_ptm_entries, phospho_missing_peptides, phospho_inferred_protein_ids = result
                            ptm_entries.extend(phospho_ptm_entries)
                            missing_peptides.extend(phospho_missing_peptides)
                            inferred_protein_ids.update(phospho_inferred_protein_ids)

                    if 'Acetylation' in modification_types:
                        args = [(chunk, uniprot_sequences) for chunk in chunked_peptide_list]
                        results = list(tqdm(pool.imap(process_peptide_acetylation, args), total=len(chunked_peptide_list), desc="Processing Acetylation"))
                        for result in results:
                            acetyl_ptm_entries, acetyl_missing_peptides, acetyl_inferred_protein_ids = result
                            ptm_entries.extend(acetyl_ptm_entries)
                            missing_peptides.extend(acetyl_missing_peptides)
                            inferred_protein_ids.update(acetyl_inferred_protein_ids)

                    if 'Ubiquitination' in modification_types:
                        args = [(chunk, uniprot_sequences) for chunk in chunked_peptide_list]
                        results = list(tqdm(pool.imap(process_peptide_ubiquitination, args), total=len(chunked_peptide_list), desc="Processing Ubiquitination"))
                        for result in results:
                            ubiquitin_ptm_entries, ubiquitin_missing_peptides, ubiquitin_inferred_protein_ids = result
                            ptm_entries.extend(ubiquitin_ptm_entries)
                            missing_peptides.extend(ubiquitin_missing_peptides)
                            inferred_protein_ids.update(ubiquitin_inferred_protein_ids)

                    if 'N-linked Glycosylation' in modification_types:
                        args = [(chunk, uniprot_sequences, 'N-linked Glycosylation') for chunk in chunked_peptide_list]
                        results = list(tqdm(pool.imap(process_peptide_glycosylation, args), total=len(chunked_peptide_list), desc="Processing N-linked Glycosylation"))
                        for result in results:
                            nlinked_ptm_entries, nlinked_missing_peptides, nlinked_inferred_protein_ids = result
                            ptm_entries.extend(nlinked_ptm_entries)
                            missing_peptides.extend(nlinked_missing_peptides)
                            inferred_protein_ids.update(nlinked_inferred_protein_ids)

                    if 'O-linked Glycosylation' in modification_types:
                        args = [(chunk, uniprot_sequences, 'O-linked Glycosylation') for chunk in chunked_peptide_list]
                        results = list(tqdm(pool.imap(process_peptide_glycosylation, args), total=len(chunked_peptide_list), desc="Processing O-linked Glycosylation"))
                        for result in results:
                            olinked_ptm_entries, olinked_missing_peptides, olinked_inferred_protein_ids = result
                            ptm_entries.extend(olinked_ptm_entries)
                            missing_peptides.extend(olinked_missing_peptides)
                            inferred_protein_ids.update(olinked_inferred_protein_ids)

                write_fasta(output_file, uniprot_sequences, ptm_entries, inferred_protein_ids, include_global_protein_entries)
                
                total_entries, unique_protein_ids = count_entries_in_fasta(output_file)
                st.write(f"Total entries in generated database: {total_entries}")
                st.write(f"Unique protein IDs in generated database: {unique_protein_ids}")

                elapsed_time = time.time() - start_time
                st.write(f"Elapsed time: {elapsed_time:.2f} seconds")

                st.success("FASTA database has been successfully created with protein and PTM entries.")
                write_missing_info(missing_info_file, missing_peptides)

                # ***Send the generated FASTA file to the backend***
                username = st.session_state.get('username', 'anonymous')  # Get the username stored in session
                response = send_fasta_to_backend(output_file, username)

            except Exception as e:
                st.error(f"An error occurred: {e}")



if __name__ == '__main__':
    if len(sys.argv) > 1:
        work_dir = sys.argv[1]
        if not st.session_state.get('work_dir'):
            st.session_state['work_dir'] = work_dir

    main()