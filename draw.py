import streamlit as st
from streamlit_ketcher import st_ketcher
import pubchempy as pcp

def get_smiles_from_name(name):
    try:
        compound = pcp.get_compounds(name, 'name')[0]
        return compound.isomeric_smiles
    except:
        return None

def main():
    st.title("Chemical Structure Editor with User Prompts")

    # User input for chemical name
    chemical_name = st.text_input("Enter a chemical name or SMILES string:")

    if chemical_name:
        # Attempt to convert name to SMILES, otherwise use the input as-is
        smiles = get_smiles_from_name(chemical_name) or chemical_name

        st.write("Draw the structure for: ", smiles)

    # Ketcher component
    drawn_smiles, molblock = st_ketcher(height=400)

    if drawn_smiles:
        st.text("Drawn SMILES string:")
        st.write(drawn_smiles)

    if molblock:
        st.text("MolBlock:")
        st.write(molblock)

if __name__ == "__main__":
    main()
