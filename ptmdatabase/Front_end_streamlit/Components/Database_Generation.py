import requests
import io
from io import StringIO
import tempfile
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

def send_fasta_to_backend(fasta_data, input_filename, username):
    # Derive the FASTA file name from the input file name
    fasta_file_name = os.path.splitext(input_filename)[0] + '.fasta'

    # Send the in-memory fasta data as a file-like object with the correct file name
    files = {'file': (fasta_file_name, fasta_data, 'text/plain')}
    data = {'username': username, 'filename': fasta_file_name}  # Send the file name to the backend
    response = requests.post("http://3.91.75.15:8000/upload-fasta/", files=files, data=data)

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

def process_peptide_phospho_acetyl_ubiquitin(args, progress_bar, status_text, total_peptides):
    chunk, uniprot_sequences, ptm_types = args  # Handle multiple PTM types
    ptm_entries, missing_peptides, inferred_protein_ids = [], [], set()

    for ptm_type in ptm_types:  # Loop through each PTM type
        for i, peptide in enumerate(chunk):
            # Generate PTM entries for each type
            entries, peptides, inferred_ids = generate_ptm_entries([peptide], uniprot_sequences, ptm_type)
            ptm_entries.extend(entries)
            missing_peptides.extend(peptides)
            inferred_protein_ids.update(inferred_ids)

            # Update progress
            progress = (i + 1) / total_peptides
            progress_bar.progress(progress)
            status_text.text(f"Processing {ptm_type} peptide {i+1}/{total_peptides} ({progress * 100:.2f}%)")

    return ptm_entries, missing_peptides, inferred_protein_ids

def process_peptide_glycosylation(args, progress_bar, status_text, total_peptides):
    chunk, uniprot_sequences, ptm_types = args
    ptm_entries, missing_peptides, inferred_protein_ids = [], [], set()

    for ptm_type in ptm_types: 
        for i, peptide in enumerate(chunk):
            entries, peptides, inferred_ids = generate_ptm_entries_glyco([peptide], uniprot_sequences, ptm_type)
            ptm_entries.extend(entries)
            missing_peptides.extend(peptides)
            inferred_protein_ids.update(inferred_ids)

            # Update progress
            progress = (i + 1) / total_peptides
            progress_bar.progress(progress)
            status_text.text(f"Processing glycosylation {i+1}/{total_peptides} ({progress * 100:.2f}%)")

    return ptm_entries, missing_peptides, inferred_protein_ids


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
        fasta_file = st.file_uploader('Upload FASTA file', type=['fasta'])
        modification_types = st.multiselect(
            'Select PTM Types to Process',
            ['Phosphorylation', 'Acetylation', 'Ubiquitination', 'N-linked Glycosylation', 'O-linked Glycosylation']
        )

        include_global_protein_entries = st.checkbox('Include Global Protein Entries', value=False)
        use_existing_db = st.checkbox("Use existing human FASTA database", key='use_existing_db')

        submit_button = st.form_submit_button(label='Generate Database')

    if submit_button:
        if not matrix_file:
            st.error("No file uploaded. Please upload a valid .xlsx or .tsv file.")
            return

        try:
            # In-memory buffers for fasta and missing info
            fasta_buffer = io.StringIO()
            missing_info_buffer = io.StringIO()

            # Parse the matrix file
            df = parse_matrix_file(matrix_file)

            # Use either the uploaded FASTA file or the existing human FASTA database
            if use_existing_db:
                uniprot_sequences = load_uniprot_sequences(st.session_state['original_fasta_dir'])
            else:
                if fasta_file is not None:
                    # Use `load_uniprot_sequences` to parse the uploaded FASTA file
                    fasta_content = io.StringIO(fasta_file.getvalue().decode("utf-8"))
                    uniprot_sequences = load_uniprot_sequences(fasta_content)
                else:
                    st.error("Please upload a FASTA file or select the existing human database.")
                    st.stop()

            peptide_list = df.iloc[:, 0].tolist()

            ptm_entries = []
            missing_peptides = []
            inferred_protein_ids = set()

            start_time = time.time()

            # Add progress bar and status text
            progress_bar = st.progress(0)
            status_text = st.empty()

            total_peptides = len(peptide_list)

            # Process phosphorylation, acetylation, ubiquitination
            if any(ptm in modification_types for ptm in ['Phosphorylation', 'Acetylation', 'Ubiquitination']):
                result = process_peptide_phospho_acetyl_ubiquitin((peptide_list, uniprot_sequences, modification_types), progress_bar, status_text, total_peptides)
                phospho_ptm_entries, phospho_missing_peptides, phospho_inferred_protein_ids = result
                ptm_entries.extend(phospho_ptm_entries)
                missing_peptides.extend(phospho_missing_peptides)
                inferred_protein_ids.update(phospho_inferred_protein_ids)

            # Process glycosylation
            if any(ptm in modification_types for ptm in ['N-linked Glycosylation', 'O-linked Glycosylation']):
                result = process_peptide_glycosylation((peptide_list, uniprot_sequences, modification_types), progress_bar, status_text, total_peptides)
                glyco_ptm_entries, glyco_missing_peptides, glyco_inferred_protein_ids = result
                ptm_entries.extend(glyco_ptm_entries)
                missing_peptides.extend(glyco_missing_peptides)
                inferred_protein_ids.update(glyco_inferred_protein_ids)

            # Write the FASTA file in-memory
            write_fasta(fasta_buffer, uniprot_sequences, ptm_entries, inferred_protein_ids, include_global_protein_entries)
            fasta_buffer.seek(0)  # Move to the start of the buffer

            # Write missing info file in-memory
            write_missing_info(missing_info_buffer, missing_peptides)
            missing_info_buffer.seek(0)  # Move to the start of the buffer

            # Display information
            total_entries, unique_protein_ids = count_entries_in_fasta(fasta_buffer)
            st.write(f"Total entries in generated database: {total_entries}")
            st.write(f"Unique protein IDs in generated database: {unique_protein_ids}")

            st.success("FASTA database has been successfully created with protein and PTM entries.")

            # Get the base name of the uploaded file (without the extension)
            input_filename = os.path.splitext(matrix_file.name)[0]

            # Create the output file name by appending .fasta to the input file name
            output_filename = f"{input_filename}.fasta"

            # Move download buttons outside the form
            fasta_buffer.seek(0)  # Move cursor to the start of the buffer
            st.download_button(
                label="Download FASTA file",
                data=fasta_buffer.getvalue(),
                file_name=output_filename,  # Use the derived file name
                mime="text/plain"
            )

            missing_info_buffer.seek(0)  # Move cursor to the start of the buffer
            st.download_button(
                label="Download Missing Info",
                data=missing_info_buffer.getvalue(),
                file_name="missing_info.csv",
                mime="text/csv"
            )

            # Send the FASTA file to the backend
            username = st.session_state.get('username', 'anonymous')
            send_fasta_to_backend(fasta_buffer.getvalue(), input_filename, username)

        except Exception as e:
            st.error(f"An error occurred: {e}")
            st.write(f"Error details: {str(e)}")  # Print detailed error message

if __name__ == '__main__':
    main()
