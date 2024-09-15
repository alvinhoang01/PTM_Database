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

    st.header("Overview")
    st.write("""
    Different from traditional databases, this method integrates phosphosites directly into the entry name and protein sequence. 
    Additionally, the application offers the option to customize the phosphosites in the protein sequence to a letter or symbol of choice...
    """)

    st.header("Database Upgrade")
    st.write("""
    Upgrade your database by adding new proteins and phosphorylation sites based on analysis results.
    This process involves:
    1. Analyzing the current database for performance.
    2. Identifying missing proteins or phosphorylation sites.
    3. Integrating new information into the existing database.
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
