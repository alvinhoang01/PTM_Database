import streamlit as st
import os
import sys
import time

def main():
    st.set_page_config(
        page_title="Phospho Proteomic Database Generation App",
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

    st.markdown("<div class='title'>Welcome to Phosphoproteomic Database Generation App</div>", unsafe_allow_html=True)
    st.write("\n")


    st.header("Overview")
    st.write("""
    Different from traditional database, this method integrate the phosphosites directly to the entry name and the protein sequence. 
    Additionally, the application also offer the function to customize the phosphosites in the protein sequence to the letter or symbol as their desire. 
    This database, even though only compatible with MSFragger at the moment since this is the only search engine that offers four letter words (B, X, Z) with customizable dynamic modification mass
    This idea aims to reduce the search space and searching time, while also achieve a similar or better performance than traditional database where only the dynamic modification is specified in the search engine parameter.        
    """)

    st.header("Database Upgrade")
    st.write("""
    Upgrade your database by adding new proteins and phosphorylation sites based on analysis results.
    This process involves:
    1. Analyzing the current database for performance.
    2. Identifying missing proteins or phosphorylation sites.
    3. Integrating the new information into the existing database.
    This ensures that your database is comprehensive and up-to-date.
    """)


    left_column2, right_column2 = st.columns([18, 2])
    with right_column2:
        st.write("\n")
        st.write("\n")
        if st.button('Exit'):
            st.error("The app has been stopped. Please close the browser window manually.")
            time.sleep(5)
            os._exit(0)

if __name__ == '__main__':
    if len(sys.argv) > 1:
        work_dir = sys.argv[1]
        if not st.session_state.get('work_dir'):
            st.session_state['work_dir'] = work_dir

    main()
