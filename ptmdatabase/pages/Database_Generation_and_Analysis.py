import streamlit as st
import pandas as pd
import os
import sys
from tools.database_tools import (
    parse_matrix_file, 
    load_uniprot_sequences, 
    load_ptm_sequences, 
    generate_ptm_entries, 
    write_fasta, 
    write_missing_info, 
    count_entries_in_fasta
)


def initialize_session_state():
    if 'work_dir' not in st.session_state:
        st.session_state['work_dir'] = ""
    if 'original_fasta_dir' not in st.session_state:
        st.session_state['original_fasta_dir'] = ""
    if 'phosphosite_fasta_dir' not in st.session_state:
        st.session_state['phosphosite_fasta_dir'] = ""
    if 'new_db_dir' not in st.session_state:
        st.session_state['new_db_dir'] = ""

def main():
    st.set_page_config(
        page_title="Database Generation and Analysis",
        page_icon="🧬",
        layout="wide",
        initial_sidebar_state="expanded",
    )

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

    st.markdown("<div class='title'>Database Generation and Analysis</div>", unsafe_allow_html=True)
    st.write("\n")

    st.sidebar.title("Toolbox")
    page = st.sidebar.radio("Select Page", ["Database Generation", "Database Analysis"])

    initialize_session_state()

    if page == "Database Generation":
        st.header("Database Generation")

        with st.form(key='database_generation_form', clear_on_submit=False):
            matrix_file = st.text_input('Frag-Pipe or MS-Pycloud Search Data File:', value=st.session_state['work_dir'])
            
            # Path to the global (pure protein) database
            st.session_state['original_fasta_dir'] = 'C:\\Users\\maitr\\Documents\\JHU\\Database_code\\qcmspycloud\\qcmspycloud\\Database_library\\uniprot_sprot.fasta'
            # Path to the phosphosite database
            st.session_state['phosphosite_fasta_dir'] = 'C:\\Users\\maitr\\Documents\\JHU\\Database_code\\qcmspycloud\\qcmspycloud\\Database_library\\Phosphosite.fasta'
            # Path to the N-linked glycosylation database
            st.session_state['nlinked_fasta_dir'] = 'C:\\Users\\maitr\\Documents\\JHU\\Database_code\\qcmspycloud\\qcmspycloud\\Database_library\\N-linked_Glycosite.fasta'
            
            new_db_dir = st.text_input('Directory to Store Generated Database:', value=st.session_state['new_db_dir'])
            st.session_state['new_db_dir'] = new_db_dir

            include_phosphosites = st.checkbox('Include Phosphosite Entries', value=False)
            include_nlinked_glycosylation = st.checkbox('Include N-linked Glycosylation Entries', value=False)

            submit_button = st.form_submit_button(label='Generate Database')
            if submit_button:
                st.session_state['work_dir'] = matrix_file
                output_file = new_db_dir
                missing_info_file = os.path.dirname(output_file)
                df = parse_matrix_file(matrix_file)
                uniprot_sequences = load_uniprot_sequences(st.session_state['original_fasta_dir'])

                # Extract unique protein IDs from the matrix file
                protein_ids = set()
                for accession_list in df['Protein.Group.Accessions']:
                    for accession in accession_list.split(';'):
                        protein_id = accession.split('|')[1] if '|' in accession else accession
                        protein_ids.add(protein_id)

                ptm_entries = []
                if include_phosphosites:
                    phosphosite_sequences = load_ptm_sequences(st.session_state['phosphosite_fasta_dir'])
                    phospho_ptm_entries, missing_phosphosites, missing_proteins, missing_peptides = generate_ptm_entries(df, uniprot_sequences, phosphosite_sequences, 'Phosphorylation')
                    ptm_entries.extend(phospho_ptm_entries)
                    write_missing_info(missing_info_file, missing_phosphosites, missing_proteins, missing_peptides)

                if include_nlinked_glycosylation:
                    nlinked_sequences = load_ptm_sequences(st.session_state['nlinked_fasta_dir'])
                    nlinked_ptm_entries, missing_nlinked, missing_proteins, missing_peptides = generate_ptm_entries(df, uniprot_sequences, nlinked_sequences, 'N-linked Glycosylation')
                    ptm_entries.extend(nlinked_ptm_entries)
                    write_missing_info(missing_info_file, missing_nlinked, missing_proteins, missing_peptides)

                # Write protein and PTM entries to the output file
                write_fasta(output_file, uniprot_sequences, protein_ids, ptm_entries)
                
                # Count the entries and unique protein IDs in the generated FASTA file
                total_entries, unique_protein_ids = count_entries_in_fasta(output_file)
                
                # Display the number of entries and unique protein IDs
                st.write(f"Total entries in generated database: {total_entries}")
                st.write(f"Unique protein IDs in generated database: {unique_protein_ids}")

                st.success("FASTA database has been successfully created with protein and PTM entries.")

if __name__ == '__main__':
    if len(sys.argv) > 1:
        work_dir = sys.argv[1]
        if not st.session_state.get('work_dir'):
            st.session_state['work_dir'] = work_dir

    main()