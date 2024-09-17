import streamlit as st
import os
import sys
import time

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

    st.markdown("<div class='title'>Welcome to Proteoform Database Generation App</div>", unsafe_allow_html=True)
    st.write("\n")

    st.header("Database Upgrade")
    st.write("""
    Upgrade your FASTA database by adding the modification version of the peptides from your search results.
    The application offers entries for the following modifications with their annotations in the protein sequence.
    - Phosphorylation [P]
    - Acetylation [A]
    - Ubiquitination [U]
    - O-linked and N-linked [NxHxFxSxGx]
    """)
    
    # Example of an exit button
    left_column2, right_column2 = st.columns([18, 2])
    with right_column2:
        st.write("\n")
        st.write("\n")
        if st.button('Exit'):
            st.error("The app has been stopped. Please close the browser window manually.")
            time.sleep(5)
            os._exit(0)

# This file will be imported by `app.py` to render the home page
