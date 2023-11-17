import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw

def name_to_smiles(name):
    try:
        compound = pcp.get_compounds(name, 'name')[0]
        return compound.isomeric_smiles
    except:
        return None

def draw_molecule(smiles):
    molecule = Chem.MolFromSmiles(smiles)
    return Draw.MolToImage(molecule)

def main():
    st.title('Chemical Structure Drawer from Name')

    chemical_name = st.text_input("Enter a chemical name:", "")

    if chemical_name:
        smiles_string = name_to_smiles(chemical_name)
        if smiles_string:
            image = draw_molecule(smiles_string)
            st.image(image)
        else:
            st.error("Could not find the chemical. Please try a different name.")

if __name__ == "__main__":
    main()
