import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw

def draw_molecule(smiles):
    molecule = Chem.MolFromSmiles(smiles)
    return Draw.MolToImage(molecule)

def main():
    st.title('Chemical Structure Drawer')

    smiles_string = st.text_input("Enter a SMILES string:", "")
    
    if smiles_string:
        try:
            image = draw_molecule(smiles_string)
            st.image(image)
        except:
            st.error("Invalid SMILES string. Please try again.")

if __name__ == "__main__":
    main()
