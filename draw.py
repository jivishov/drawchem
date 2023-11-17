import streamlit as st
import pubchempy as pcp
import requests
from PIL import Image
from io import BytesIO

def name_to_smiles(name):
    try:
        compound = pcp.get_compounds(name, 'name')[0]
        return compound.isomeric_smiles
    except:
        return None

def smiles_to_image(smiles, width=300, height=300):
    url = f"https://cactus.nci.nih.gov/chemical/structure/{smiles}/image"
    params = {'width': width, 'height': height}
    response = requests.get(url, params=params)
    if response.status_code == 200:
        return Image.open(BytesIO(response.content))
    return None

def main():
    st.title('Chemical Structure Drawer from Name')

    chemical_name = st.text_input("Enter a chemical name:", "")

    if chemical_name:
        smiles_string = name_to_smiles(chemical_name)
        if smiles_string:
            image = smiles_to_image(smiles_string)
            if image:
                st.image(image, use_column_width=True)
            else:
                st.error("Could not generate the image for this chemical.")
        else:
            st.error("Could not find the chemical. Please try a different name.")

if __name__ == "__main__":
    main()
