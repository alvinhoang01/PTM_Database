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

def send_fasta_to_backend(fasta_file_path, email):
    with open(fasta_file_path, "rb") as file:
        files = {'file': file}
        data = {'email': email}  # Include the email in the form data
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
        matrix_file = st.text_input('Peptide List (xlsx or tsv):', value=st.session_state['work_dir'])
        new_db_dir = st.text_input('Directory to Store Generated Database:', value=st.session_state['new_db_dir'])
        st.session_state['new_db_dir'] = new_db_dir

        modification_types = st.multiselect(
            'Select PTM Types to Process',
            ['Phosphorylation', 'Acetylation', 'Ubiquitination', 'N-linked Glycosylation', 'O-linked Glycosylation']
        )

        include_global_protein_entries = st.checkbox('Include Global Protein Entries', value=False)

        submit_button = st.form_submit_button(label='Generate Database')
        if submit_button:
            st.session_state['work_dir'] = matrix_file
            output_file = new_db_dir
            missing_info_file = os.path.dirname(output_file)
            st.session_state['missing_info_file'] = missing_info_file
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

                # (Add other modification processing here)

            write_fasta(output_file, uniprot_sequences, ptm_entries, inferred_protein_ids, include_global_protein_entries)
            
            total_entries, unique_protein_ids = count_entries_in_fasta(output_file)
            st.write(f"Total entries in generated database: {total_entries}")
            st.write(f"Unique protein IDs in generated database: {unique_protein_ids}")

            elapsed_time = time.time() - start_time
            st.write(f"Elapsed time: {elapsed_time:.2f} seconds")

            st.success("FASTA database has been successfully created with protein and PTM entries.")
            write_missing_info(missing_info_file, missing_peptides)

            # **Add this line to send the generated FASTA file to the backend**
            response = send_fasta_to_backend(output_file, st.session_state['email'])

if __name__ == '__main__':
    if len(sys.argv) > 1:
        work_dir = sys.argv[1]
        if not st.session_state.get('work_dir'):
            st.session_state['work_dir'] = work_dir

    main()
