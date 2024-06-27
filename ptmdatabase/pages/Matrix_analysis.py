import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3

# Function to plot bar chart
def plot_bar_chart(data, indices, legend=False):
    df = pd.DataFrame(data, index=indices)
    ax = df.plot(kind='bar', stacked=True, figsize=(10, 6))
    plt.xlabel('Database')
    plt.ylabel('Counts')
    plt.xticks(rotation=0)
    if legend:
        plt.legend(title="Modification Sites", loc='upper right')
    else:
        ax.get_legend().remove()
    st.pyplot(plt)

# Function to plot Venn diagram for two sets
def plot_venn_diagram(set1, set2, label1, label2):
    fig, ax = plt.subplots()
    v = venn2([set1, set2], set_labels=(label1, label2))
    for label_id in ['10', '01', '11']:
        if v.get_label_by_id(label_id):
            v.get_label_by_id(label_id).set_text(len(eval(f"set1 & set2" if label_id == '11' else f"set1 - set2" if label_id == '10' else "set2 - set1")))
    st.pyplot(fig)

# Function to plot Venn diagram for three sets
def plot_venn_diagram_three_sets(set1, set2, set3, label1, label2, label3):
    fig, ax = plt.subplots()
    v = venn3([set1, set2, set3], set_labels=(label1, label2, label3))
    
    # Define the labels and the corresponding sets
    labels = {
        '100': set1 - set2 - set3,
        '010': set2 - set1 - set3,
        '001': set3 - set1 - set2,
        '110': (set1 & set2) - set3,
        '101': (set1 & set3) - set2,
        '011': (set2 & set3) - set1,
        '111': set1 & set2 & set3
    }
    
    # Print the contents of each intersection for debugging
    for label_id, subset in labels.items():
        print(f"Label {label_id}: {len(subset)} elements")
    
    # Iterate over the labels to set the text
    for label_id, subset in labels.items():
        if v.get_label_by_id(label_id):
            v.get_label_by_id(label_id).set_text(len(subset))
        else:
            # Create a text label for regions that have no elements and no existing label
            label_positions = {'100': (-0.4, 0.2), '010': (0.4, 0.2), '001': (0.0, -0.4),
                               '110': (0.0, 0.4), '101': (-0.2, -0.2), '011': (0.2, -0.2), '111': (0.0, 0.0)}
            x, y = label_positions[label_id]
            ax.text(x, y, str(len(subset)), ha='center', va='center')

    st.pyplot(fig)

# Load and preprocess data
def load_and_preprocess_data(file_path):
    df = pd.read_csv(file_path, sep='\t')
    df.columns = df.columns.str.strip()
    return df

# Extract core protein ID
def extract_core_protein_id(df):
    df['Core Protein ID'] = df['Protein ID'].str.strip()
    return df

# Replace modified amino acids in the modified dataset
def replace_modified_amino_acids(df):
    df['Peptide'] = df['Peptide'].str.replace('B', 'S').str.replace('Z', 'T').str.replace('X', 'Y')
    return df

# Extract phosphorylation sites
def extract_phosphorylation_sites(mod_str, is_modified=False):
    if pd.isna(mod_str) or mod_str == '':
        return set()
    if is_modified:
        mod_str = mod_str.replace('B', 'S').replace('Z', 'T').replace('X', 'Y')
    return set(m.split('(')[0] for m in mod_str.replace(')', '').split(', ') if '79.9663' in m or '181.0160' in m or '166.9960' in m or '243.0260' in m)

# Extract detailed phosphorylation site counts
def extract_detailed_phosphorylation_sites(mod_str, is_modified=False):
    if pd.isna(mod_str) or mod_str == '':
        return {'S': 0, 'T': 0, 'Y': 0}
    if is_modified:
        mod_str = mod_str.replace('B', 'S').replace('Z', 'T').replace('X', 'Y')
    counts = {'S': 0, 'T': 0, 'Y': 0}
    for m in mod_str.replace(')', '').split(', '):
        if '79.9663' in m or '181.0160' in m or '166.9960' in m or '243.0260' in m:
            if 'S' in m:
                counts['S'] += 1
            elif 'T' in m:
                counts['T'] += 1
            elif 'Y' in m:
                counts['Y'] += 1
    return counts

# Calculate similar and unique counts
def calculate_counts(original_set, modified_set):
    similar_count = len(original_set & modified_set)
    unique_original_count = len(original_set - modified_set)
    unique_modified_count = len(modified_set - original_set)
    return similar_count, unique_original_count, unique_modified_count

def calculate_counts_three_sets(set1, set2, set3):
    similar_all = set1.intersection(set2).intersection(set3)
    unique_set1 = set1.difference(set2).difference(set3)
    unique_set2 = set2.difference(set1).difference(set3)
    unique_set3 = set3.difference(set1).difference(set2)
    return len(similar_all), len(unique_set1), len(unique_set2), len(unique_set3)

# Main Streamlit UI
st.title("Matrix Data Analysis")

original_path = st.text_input('Enter the directory for the original data matrix:', '')
modified_path = st.text_input('Enter the directory for the modified data matrix:', '')
modified_path_v2 = st.text_input('Enter the directory for the modified_v2 data matrix:', '')

if st.button('Analyze'):
    if original_path or modified_path or modified_path_v2:
        original_df, modified_df, modified_df_v2 = None, None, None
        if original_path:
            original_df = load_and_preprocess_data(original_path)
            original_df = extract_core_protein_id(original_df)
        if modified_path:
            modified_df = load_and_preprocess_data(modified_path)
            modified_df = extract_core_protein_id(modified_df)
            modified_df = replace_modified_amino_acids(modified_df)
        if modified_path_v2:
            modified_df_v2 = load_and_preprocess_data(modified_path_v2)
            modified_df_v2 = extract_core_protein_id(modified_df_v2)
            modified_df_v2 = replace_modified_amino_acids(modified_df_v2)

        if original_df is not None:
            original_df['Assigned Modifications'] = original_df['Assigned Modifications'].fillna('')
        if modified_df is not None:
            modified_df['Assigned Modifications'] = modified_df['Assigned Modifications'].fillna('')
        if modified_df_v2 is not None:
            modified_df_v2['Assigned Modifications'] = modified_df_v2['Assigned Modifications'].fillna('')

        # Protein sets
        original_protein_set = set(original_df['Core Protein ID']) if original_df is not None else set()
        modified_protein_set = set(modified_df['Core Protein ID']) if modified_df is not None else set()
        modified_protein_set_v2 = set(modified_df_v2['Core Protein ID']) if modified_df_v2 is not None else set()
        
        # Peptide sets
        original_peptide_set = set(original_df['Peptide']) if original_df is not None else set()
        modified_peptide_set = set(modified_df['Peptide']) if modified_df is not None else set()
        modified_peptide_set_v2 = set(modified_df_v2['Peptide']) if modified_df_v2 is not None else set()

        # Phosphorylation sites
        original_sites, modified_sites, modified_sites_v2 = {}, {}, {}
        phospho_counts_original, phospho_counts_modified, phospho_counts_modified_v2 = {'S': 0, 'T': 0, 'Y': 0}, {'S': 0, 'T': 0, 'Y': 0}, {'S': 0, 'T': 0, 'Y': 0}
        if original_df is not None:
            for _, row in original_df.iterrows():
                peptide = row['Peptide']
                sites = extract_phosphorylation_sites(row['Assigned Modifications'])
                counts = extract_detailed_phosphorylation_sites(row['Assigned Modifications'])
                phospho_counts_original['S'] += counts['S']
                phospho_counts_original['T'] += counts['T']
                phospho_counts_original['Y'] += counts['Y']
                if peptide in original_sites:
                    original_sites[peptide].update(sites)
                else:
                    original_sites[peptide] = sites

        if modified_df is not None:
            for _, row in modified_df.iterrows():
                peptide = row['Peptide']
                sites = extract_phosphorylation_sites(row['Assigned Modifications'], is_modified=True)
                counts = extract_detailed_phosphorylation_sites(row['Assigned Modifications'], is_modified=True)
                phospho_counts_modified['S'] += counts['S']
                phospho_counts_modified['T'] += counts['T']
                phospho_counts_modified['Y'] += counts['Y']
                if peptide in modified_sites:
                    modified_sites[peptide].update(sites)
                else:
                    modified_sites[peptide] = sites

        if modified_df_v2 is not None:
            for _, row in modified_df_v2.iterrows():
                peptide = row['Peptide']
                sites = extract_phosphorylation_sites(row['Assigned Modifications'], is_modified=True)
                counts = extract_detailed_phosphorylation_sites(row['Assigned Modifications'], is_modified=True)
                phospho_counts_modified_v2['S'] += counts['S']
                phospho_counts_modified_v2['T'] += counts['T']
                phospho_counts_modified_v2['Y'] += counts['Y']
                if peptide in modified_sites_v2:
                    modified_sites_v2[peptide].update(sites)
                else:
                    modified_sites_v2[peptide] = sites

        all_orig_sites = set()
        all_mod_sites = set()
        all_mod_sites_v2 = set()

        if original_df is not None:
            for peptide, orig_sites in original_sites.items():
                all_orig_sites.update(f"{peptide}_{site}" for site in orig_sites)
        if modified_df is not None:
            for peptide, mod_sites in modified_sites.items():
                all_mod_sites.update(f"{peptide}_{site}" for site in mod_sites)
        if modified_df_v2 is not None:
            for peptide, mod_sites_v2 in modified_sites_v2.items():
                all_mod_sites_v2.update(f"{peptide}_{site}" for site in mod_sites_v2)

        # Plotting and visualization
        if original_df is not None and modified_df is not None and modified_df_v2 is not None:
            st.write("### Protein counts")
            col1, col2 = st.columns([11.9, 7.5])
            with col1:
                plot_bar_chart({'Protein Counts': [len(original_protein_set), len(modified_protein_set), len(modified_protein_set_v2)]}, ['Original', 'Modified v1', 'Modified v2'])
            with col2:
                plot_venn_diagram_three_sets(original_protein_set, modified_protein_set, modified_protein_set_v2, 'Original', 'Modified_v1', 'Modified_v2')

            st.write("### Peptide Counts")
            col3, col4 = st.columns([11.9, 7.5])
            with col3:
                plot_bar_chart({'Peptide Counts': [len(original_peptide_set), len(modified_peptide_set), len(modified_peptide_set_v2)]}, ['Original', 'Modified v1', 'Modified v2'])
            with col4:
                plot_venn_diagram_three_sets(original_peptide_set, modified_peptide_set, modified_peptide_set_v2, 'Original', 'Modified_v1', 'Modified_v2')

            st.write("### Phosphorylation Site Counts")
            col5, col6 = st.columns([11.9, 7.7])
            with col5:
                plot_bar_chart({'Phospho S': [phospho_counts_original['S'], phospho_counts_modified['S'], phospho_counts_modified_v2['S']],
                                'Phospho T': [phospho_counts_original['T'], phospho_counts_modified['T'], phospho_counts_modified_v2['T']],
                                'Phospho Y': [phospho_counts_original['Y'], phospho_counts_modified['Y'], phospho_counts_modified_v2['Y']]}, 
                               ['Original', 'Modified v1', 'Modified v2'], legend=True)
            with col6:
                plot_venn_diagram_three_sets(all_orig_sites, all_mod_sites, all_mod_sites_v2, 'Original', 'Modified_v1', 'Modified_v2')

            
            # Corrected Summary Table for Original, Modified v1, and Modified v2
            st.write("### Summary Table for Original, Modified v1, and Modified v2")
            summary_data = {
                'Metric': [
                    '**Proteins**', 'Total Count', 'Similar Count', 'Unique Count',
                    '**Peptides**', 'Total Count', 'Similar Count', 'Unique Count',
                    '**Phosphorylation Sites**', 'Total Count', 'Similar Count', 'Unique Count', 
                    '**Phospho S**', 'Total Count', 'Similar Count', 'Unique Count',
                    '**Phospho T**', 'Total Count', 'Similar Count', 'Unique Count',
                    '**Phospho Y**', 'Total Count', 'Similar Count', 'Unique Count'
                ],
                'Original': [
                    '', len(original_protein_set), len(original_protein_set & modified_protein_set & modified_protein_set_v2), len(original_protein_set - modified_protein_set - modified_protein_set_v2),
                    '', len(original_peptide_set), len(original_peptide_set & modified_peptide_set & modified_peptide_set_v2), len(original_peptide_set - modified_peptide_set - modified_peptide_set_v2),
                    '', len(all_orig_sites), len(all_orig_sites & all_mod_sites & all_mod_sites_v2), len(all_orig_sites - all_mod_sites - all_mod_sites_v2),
                    '', phospho_counts_original['S'], min(phospho_counts_original['S'], phospho_counts_modified['S'], phospho_counts_modified_v2['S']), max(phospho_counts_original['S'] - min(phospho_counts_modified['S'], phospho_counts_modified_v2['S']), 0),
                    '', phospho_counts_original['T'], min(phospho_counts_original['T'], phospho_counts_modified['T'], phospho_counts_modified_v2['T']), max(phospho_counts_original['T'] - min(phospho_counts_modified['T'], phospho_counts_modified_v2['T']), 0),
                    '', phospho_counts_original['Y'], min(phospho_counts_original['Y'], phospho_counts_modified['Y'], phospho_counts_modified_v2['Y']), max(phospho_counts_original['Y'] - min(phospho_counts_modified['Y'], phospho_counts_modified_v2['Y']), 0)
                ],
                'Modified v1': [
                    '', len(modified_protein_set), len(original_protein_set & modified_protein_set & modified_protein_set_v2), len(modified_protein_set - original_protein_set - modified_protein_set_v2),
                    '', len(modified_peptide_set), len(original_peptide_set & modified_peptide_set & modified_peptide_set_v2), len(modified_peptide_set - original_peptide_set - modified_peptide_set_v2),
                    '', len(all_mod_sites), len(all_orig_sites & all_mod_sites & all_mod_sites_v2), len(all_mod_sites - all_orig_sites - all_mod_sites_v2),
                    '', phospho_counts_modified['S'], min(phospho_counts_original['S'], phospho_counts_modified['S'], phospho_counts_modified_v2['S']), max(phospho_counts_modified['S'] - min(phospho_counts_original['S'], phospho_counts_modified_v2['S']), 0),
                    '', phospho_counts_modified['T'], min(phospho_counts_original['T'], phospho_counts_modified['T'], phospho_counts_modified_v2['T']), max(phospho_counts_modified['T'] - min(phospho_counts_original['T'], phospho_counts_modified_v2['T']), 0),
                    '', phospho_counts_modified['Y'], min(phospho_counts_original['Y'], phospho_counts_modified['Y'], phospho_counts_modified_v2['Y']), max(phospho_counts_modified['Y'] - min(phospho_counts_original['Y'], phospho_counts_modified_v2['Y']), 0)
                ],
                'Modified v2': [
                    '', len(modified_protein_set_v2), len(original_protein_set & modified_protein_set & modified_protein_set_v2), len(modified_protein_set_v2 - original_protein_set - modified_protein_set),
                    '', len(modified_peptide_set_v2), len(original_peptide_set & modified_peptide_set & modified_peptide_set_v2), len(modified_peptide_set_v2 - original_peptide_set - modified_peptide_set),
                    '', len(all_mod_sites_v2), len(all_orig_sites & all_mod_sites & all_mod_sites_v2), len(all_mod_sites_v2 - all_orig_sites - all_mod_sites),
                    '', phospho_counts_modified_v2['S'], min(phospho_counts_original['S'], phospho_counts_modified['S'], phospho_counts_modified_v2['S']), max(phospho_counts_modified_v2['S'] - min(phospho_counts_original['S'], phospho_counts_modified['S']), 0),
                    '', phospho_counts_modified_v2['T'], min(phospho_counts_original['T'], phospho_counts_modified['T'], phospho_counts_modified_v2['T']), max(phospho_counts_modified_v2['T'] - min(phospho_counts_original['T'], phospho_counts_modified['T']), 0),
                    '', phospho_counts_modified_v2['Y'], min(phospho_counts_original['Y'], phospho_counts_modified['Y'], phospho_counts_modified_v2['Y']), max(phospho_counts_modified_v2['Y'] - min(phospho_counts_original['Y'], phospho_counts_modified['Y']), 0)
                ]
            }

            summary_df = pd.DataFrame(summary_data)
            st.table(summary_df.astype(str))

        elif original_df is not None and modified_df is not None:
            st.write("### Comparison between Original and Modified v1 Databases")
            col1, col2 = st.columns([11.9, 7.5])
            with col1:
                plot_bar_chart({'Counts': [len(original_protein_set), len(modified_protein_set)]}, ['Original', 'Modified v1'])
            with col2:
                plot_venn_diagram(original_protein_set, modified_protein_set, 'Original Proteins', 'Modified v1 Proteins')

            st.write("### Peptide Counts")
            col3, col4 = st.columns([11.9, 7.5])
            with col3:
                plot_bar_chart({'Counts': [len(original_peptide_set), len(modified_peptide_set)]}, ['Original', 'Modified v1'])
            with col4:
                plot_venn_diagram(original_peptide_set, modified_peptide_set, 'Original Peptides', 'Modified v1 Peptides')

            st.write("### Phosphorylation Site Counts")
            col5, col6 = st.columns([11.9, 7.7])
            with col5:
                plot_bar_chart({'Phospho S': [phospho_counts_original['S'], phospho_counts_modified['S']],
                                'Phospho T': [phospho_counts_original['T'], phospho_counts_modified['T']],
                                'Phospho Y': [phospho_counts_original['Y'], phospho_counts_modified['Y']]}, 
                               ['Original', 'Modified v1'], legend=True)
            with col6:
                plot_venn_diagram(all_orig_sites, all_mod_sites, 'Original Phosphosites', 'Modified v1 Phosphosites')


            # Summary Table for Original and Modified v1
            st.write("### Summary Table for Original and Modified v1")
            summary_data = {
                'Metric': [
                    '**Proteins**', 'Total Count', 'Similar Count', 'Unique Count',
                    '**Peptides**', 'Total Count', 'Similar Count', 'Unique Count',
                    '**Phosphorylation Sites**', 'Total Count', 'Similar Count', 'Unique Count',
                    '**Phospho S**', 'Total Count', 'Similar Count', 'Unique Count',
                    '**Phospho T**', 'Total Count', 'Similar Count', 'Unique Count',
                    '**Phospho Y**', 'Total Count', 'Similar Count', 'Unique Count'
                ],
                'Original': [
                    '', len(original_protein_set), len(original_protein_set & modified_protein_set), len(original_protein_set - modified_protein_set),
                    '', len(original_peptide_set), len(original_peptide_set & modified_peptide_set), len(original_peptide_set - modified_peptide_set),
                    '', len(all_orig_sites), len(all_orig_sites & all_mod_sites), len(all_orig_sites - all_mod_sites),
                    '', phospho_counts_original['S'], min(phospho_counts_original['S'], phospho_counts_modified['S']), phospho_counts_original['S'] - phospho_counts_modified['S'],
                    '', phospho_counts_original['T'], min(phospho_counts_original['T'], phospho_counts_modified['T']), phospho_counts_original['T'] - phospho_counts_modified['T'],
                    '', phospho_counts_original['Y'], min(phospho_counts_original['Y'], phospho_counts_modified['Y']), phospho_counts_original['Y'] - phospho_counts_modified['Y']
                ],
                'Modified v1': [
                    '', len(modified_protein_set), len(original_protein_set & modified_protein_set), len(modified_protein_set - original_protein_set),
                    '', len(modified_peptide_set), len(original_peptide_set & modified_peptide_set), len(modified_peptide_set - original_peptide_set),
                    '', len(all_mod_sites), len(all_orig_sites & all_mod_sites), len(all_mod_sites - all_orig_sites),
                    '', phospho_counts_modified['S'], min(phospho_counts_original['S'], phospho_counts_modified['S']), phospho_counts_modified['S'] - phospho_counts_original['S'],
                    '', phospho_counts_modified['T'], min(phospho_counts_original['T'], phospho_counts_modified['T']), phospho_counts_modified['T'] - phospho_counts_original['T'],
                    '', phospho_counts_modified['Y'], min(phospho_counts_original['Y'], phospho_counts_modified['Y']), phospho_counts_modified['Y'] - phospho_counts_original['Y']
                ]
            }
            summary_df = pd.DataFrame(summary_data)
            st.table(summary_df.astype(str))

        elif original_df is not None and modified_df_v2 is not None:
            st.write("### Comparison between Original and Modified v2 Databases")
            col1, col2 = st.columns([11.9, 7.5])
            with col1:
                plot_bar_chart({'Counts': [len(original_protein_set), len(modified_protein_set_v2)]}, ['Original', 'Modified v2'])
            with col2:
                plot_venn_diagram(original_protein_set, modified_protein_set_v2, 'Original Proteins', 'Modified v2 Proteins')

            st.write("### Peptide Counts")
            col3, col4 = st.columns([11.9, 7.5])
            with col3:
                plot_bar_chart({'Counts': [len(original_peptide_set), len(modified_peptide_set_v2)]}, ['Original', 'Modified v2'])
            with col4:
                plot_venn_diagram(original_peptide_set, modified_peptide_set_v2, 'Original Peptides', 'Modified v2 Peptides')

            st.write("### Phosphorylation Site Counts")
            col5, col6 = st.columns([11.9, 7.7])
            with col5:
                plot_bar_chart({'Phospho S': [phospho_counts_original['S'], phospho_counts_modified_v2['S']],
                                'Phospho T': [phospho_counts_original['T'], phospho_counts_modified_v2['T']],
                                'Phospho Y': [phospho_counts_original['Y'], phospho_counts_modified_v2['Y']]}, 
                               ['Original', 'Modified v2'], legend=True)
            with col6:
                plot_venn_diagram(all_orig_sites, all_mod_sites_v2, 'Original Phosphosites', 'Modified v2 Phosphosites')

            # Summary Table for Original and Modified v2
            st.write("### Summary Table for Original and Modified v2")
            summary_data = {
                'Metric': [
                    '**Proteins**', 'Total Count', 'Similar Count', 'Unique Count',
                    '**Peptides**', 'Total Count', 'Similar Count', 'Unique Count',
                    '**Phosphorylation Sites**', 'Total Count', 'Similar Count', 'Unique Count',
                    '**Phospho S**', 'Total Count', 'Similar Count', 'Unique Count',
                    '**Phospho T**', 'Total Count', 'Similar Count', 'Unique Count',
                    '**Phospho Y**', 'Total Count', 'Similar Count', 'Unique Count'
                ],
                'Original': [
                    '', len(original_protein_set), len(original_protein_set & modified_protein_set_v2), len(original_protein_set - modified_protein_set_v2),
                    '', len(original_peptide_set), len(original_peptide_set & modified_peptide_set_v2), len(original_peptide_set - modified_peptide_set_v2),
                    '', len(all_orig_sites), len(all_orig_sites & all_mod_sites_v2), len(all_orig_sites - all_mod_sites_v2),
                    '', phospho_counts_original['S'], min(phospho_counts_original['S'], phospho_counts_modified_v2['S']), phospho_counts_original['S'] - phospho_counts_modified_v2['S'],
                    '', phospho_counts_original['T'], min(phospho_counts_original['T'], phospho_counts_modified_v2['T']), phospho_counts_original['T'] - phospho_counts_modified_v2['T'],
                    '', phospho_counts_original['Y'], min(phospho_counts_original['Y'], phospho_counts_modified_v2['Y']), phospho_counts_original['Y'] - phospho_counts_modified_v2['Y']
                ],
                'Modified v2': [
                    '', len(modified_protein_set_v2), len(original_protein_set & modified_protein_set_v2), len(modified_protein_set_v2 - original_protein_set),
                    '', len(modified_peptide_set_v2), len(original_peptide_set & modified_peptide_set_v2), len(modified_peptide_set_v2 - original_peptide_set),
                    '', len(all_mod_sites_v2), len(all_orig_sites & all_mod_sites_v2), len(all_mod_sites_v2 - all_orig_sites),
                    '', phospho_counts_modified_v2['S'], min(phospho_counts_original['S'], phospho_counts_modified_v2['S']), phospho_counts_modified_v2['S'] - phospho_counts_original['S'],
                    '', phospho_counts_modified_v2['T'], min(phospho_counts_original['T'], phospho_counts_modified_v2['T']), phospho_counts_modified_v2['T'] - phospho_counts_original['T'],
                    '', phospho_counts_modified_v2['Y'], min(phospho_counts_original['Y'], phospho_counts_modified_v2['Y']), phospho_counts_modified_v2['Y'] - phospho_counts_original['Y']
                ]
            }
            summary_df = pd.DataFrame(summary_data)
            st.table(summary_df.astype(str))

        elif modified_df is not None and modified_df_v2 is not None:
            st.write("### Comparison between Modified v1 and Modified v2 Databases")
            col1, col2 = st.columns([11.9, 7.5])
            with col1:
                plot_bar_chart({'Counts': [len(modified_protein_set), len(modified_protein_set_v2)]}, ['Modified v1', 'Modified v2'])
            with col2:
                plot_venn_diagram(modified_protein_set, modified_protein_set_v2, 'Modified v1 Proteins', 'Modified v2 Proteins')

            st.write("### Peptide Counts")
            col3, col4 = st.columns([11.9, 7.5])
            with col3:
                plot_bar_chart({'Counts': [len(modified_peptide_set), len(modified_peptide_set_v2)]}, ['Modified v1', 'Modified v2'])
            with col4:
                plot_venn_diagram(modified_peptide_set, modified_peptide_set_v2, 'Modified v1 Peptides', 'Modified v2 Peptides')

            st.write("### Phosphorylation Site Counts")
            col5, col6 = st.columns([11.9, 7.7])
            with col5:
                plot_bar_chart({'Phospho S': [phospho_counts_modified['S'], phospho_counts_modified_v2['S']],
                                'Phospho T': [phospho_counts_modified['T'], phospho_counts_modified_v2['T']],
                                'Phospho Y': [phospho_counts_modified['Y'], phospho_counts_modified_v2['Y']]}, 
                               ['Modified v1', 'Modified v2'], legend=True)
            with col6:
                plot_venn_diagram(all_mod_sites, all_mod_sites_v2, 'Modified v1 Phosphosites', 'Modified v2 Phosphosites')


            # Summary Table for Modified v1 and Modified v2
            st.write("### Summary Table for Modified v1 and Modified v2")
            summary_data = {
                'Metric': [
                    '**Proteins**', 'Total Count', 'Similar Count', 'Unique Count',
                    '**Peptides**', 'Total Count', 'Similar Count', 'Unique Count',
                    '**Phosphorylation Sites**', 'Total Count', 'Similar Count', 'Unique Count',
                    '**Phospho S**', 'Total Count', 'Similar Count', 'Unique Count',
                    '**Phospho T**', 'Total Count', 'Similar Count', 'Unique Count',
                    '**Phospho Y**', 'Total Count', 'Similar Count', 'Unique Count'
                ],
                'Modified v1': [
                    '', len(modified_protein_set), len(modified_protein_set & modified_protein_set_v2), len(modified_protein_set - modified_protein_set_v2),
                    '', len(modified_peptide_set), len(modified_peptide_set & modified_peptide_set_v2), len(modified_peptide_set - modified_peptide_set_v2),
                    '', len(all_mod_sites), len(all_mod_sites & all_mod_sites_v2), len(all_mod_sites - all_mod_sites_v2),
                    '', phospho_counts_modified['S'], min(phospho_counts_modified['S'], phospho_counts_modified_v2['S']), phospho_counts_modified['S'] - phospho_counts_modified_v2['S'],
                    '', phospho_counts_modified['T'], min(phospho_counts_modified['T'], phospho_counts_modified_v2['T']), phospho_counts_modified['T'] - phospho_counts_modified_v2['T'],
                    '', phospho_counts_modified['Y'], min(phospho_counts_modified['Y'], phospho_counts_modified_v2['Y']), phospho_counts_modified['Y'] - phospho_counts_modified_v2['Y']
                ],
                'Modified v2': [
                    '', len(modified_protein_set_v2), len(modified_protein_set & modified_protein_set_v2), len(modified_protein_set_v2 - modified_protein_set),
                    '', len(modified_peptide_set_v2), len(modified_peptide_set & modified_peptide_set_v2), len(modified_peptide_set_v2 - modified_peptide_set),
                    '', len(all_mod_sites_v2), len(all_mod_sites & all_mod_sites_v2), len(all_mod_sites_v2 - all_mod_sites),
                    '', phospho_counts_modified_v2['S'], min(phospho_counts_modified['S'], phospho_counts_modified_v2['S']), phospho_counts_modified_v2['S'] - phospho_counts_modified['S'],
                    '', phospho_counts_modified_v2['T'], min(phospho_counts_modified['T'], phospho_counts_modified_v2['T']), phospho_counts_modified_v2['T'] - phospho_counts_modified['T'],
                    '', phospho_counts_modified_v2['Y'], min(phospho_counts_modified['Y'], phospho_counts_modified_v2['Y']), phospho_counts_modified_v2['Y'] - phospho_counts_modified['Y']
                ]
            }
            summary_df = pd.DataFrame(summary_data)
            st.table(summary_df.astype(str))

    else:
        st.error('Please enter at least one directory.')
