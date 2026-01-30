import streamlit as st
import py3dmol
from Bio.Seq import Seq
import numpy as np

# --- НАСТРОЙКИ ---
st.set_page_config(page_title="Sfiral Protein Lab", layout="wide", page_icon="🧬")

# --- СТИЛИ ---
st.markdown("""
<style>
    .stApp {background-color: #0e1117; color: #fff;}
    h1 {color: #00CCFF;}
</style>
""", unsafe_allow_html=True)

st.title("🧬 Protein-Sfiral: Time-Genetics Folding")
st.caption("Testing the Kushelev Hypothesis: Same Amino Acids -> Different Geometry (CDS-driven)")

# --- ЛОГИКА ВРЕМЕНИ ---
TIME_GENETICS_MAP = {
    'AAA': {'aa': 'K', 'phi': -65, 'psi': -40, 'delay': 1.0, 'note': 'Fast (Pi-Helix)'},
    'AAG': {'aa': 'K', 'phi': -57, 'psi': -47, 'delay': 1.5, 'note': 'Slow (Alpha-Helix)'},
    'DEFAULT': {'aa': '?', 'phi': -60, 'psi': -45, 'delay': 1.0}
}

def get_codon_params(codon):
    return TIME_GENETICS_MAP.get(codon, TIME_GENETICS_MAP['DEFAULT'])

# --- ИНТЕРФЕЙС ---
col1, col2 = st.columns([1, 2])

with col1:
    st.subheader("📥 Ввод (CDS)")
    dna_input = st.text_area("Введите последовательность кодонов:", height=150, placeholder="AAA AAA AAG AAG...")
    
    sequence = ""
    if dna_input:
        sequence = dna_input.replace("\n", "").replace(" ", "").upper()
        st.success(f"Кодонов: {len(sequence)//3}")

# --- ГЕНЕРАЦИЯ PDB ---
def generate_structure(dna_seq):
    codons = [dna_seq[i:i+3] for i in range(0, len(dna_seq), 3)]
    pdb_str = ""
    atom_id = 1
    res_id = 1
    x, y, z = 0.0, 0.0, 0.0
    
    for codon in codons:
        params = get_codon_params(codon)
        phi = params['phi']
        
        # Геометрия сворачивания
        x += 1.5 * np.cos(np.radians(phi))
        y += 1.5 * np.sin(np.radians(phi))
        z += 0.8 
        
        # Формируем атом
        pdb_str += f"ATOM  {atom_id:5d}  CA  LYS A{res_id:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 {params['delay']:5.2f}           C\n"
        atom_id += 1
        res_id += 1
    return pdb_str

# --- ВИЗУАЛИЗАЦИЯ ---
with col2:
    if sequence and len(sequence) % 3 == 0:
        pdb = generate_structure(sequence)
        
        view = py3dmol.view(width=600, height=400)
        view.addModel(pdb, 'pdb')
        view.setStyle({'stick': {'radius': 0.2}, 'sphere': {'scale': 0.3}})
        view.setStyle({'cartoon': {'color': 'spectrum'}})
        view.zoomTo()
        
        output = view._make_html()
        st.components.v1.html(output, width=600, height=400)
