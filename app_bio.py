import streamlit as st
import numpy as np
import py3dmol
from stmol import showmol  # <--- Профессиональная библиотека для рендера
from Bio.Seq import Seq

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

# --- ЛОГИКА ВРЕМЯГЕНЕТИКИ ---
TIME_GENETICS_MAP = {
    'AAA': {'aa': 'K', 'phi': -65, 'psi': -40, 'delay': 1.0, 'note': 'Fast (Pi-Helix)'},
    'AAG': {'aa': 'K', 'phi': -57, 'psi': -47, 'delay': 1.5, 'note': 'Slow (Alpha-Helix)'},
    'DEFAULT': {'aa': '?', 'phi': -60, 'psi': -45, 'delay': 1.0}
}

def get_codon_params(codon):
    return TIME_GENETICS_MAP.get(codon, TIME_GENETICS_MAP['DEFAULT'])

# --- ИНТЕРФЕЙС ВВОДА ---
col1, col2 = st.columns([1, 2])

with col1:
    st.subheader("📥 Ввод (CDS)")
    dna_input = st.text_area("Введите последовательность кодонов:", height=150, placeholder="AAA AAA AAG AAG...")
    
    sequence = ""
    if dna_input:
        sequence = dna_input.replace("\n", "").replace(" ", "").upper()
        # Проверка на валидность
        if len(sequence) % 3 == 0 and len(sequence) > 0:
            st.success(f"✅ Кодонов загружено: {len(sequence)//3}")
        else:
            if len(sequence) > 0:
                st.error("⚠ Длина последовательности должна делиться на 3")

# --- ГЕНЕРАЦИЯ СТРУКТУРЫ ---
def generate_structure(dna_seq):
    codons = [dna_seq[i:i+3] for i in range(0, len(dna_seq), 3)]
    pdb_str = ""
    atom_id = 1
    res_id = 1
    x, y, z = 0.0, 0.0, 0.0
    
    for codon in codons:
        params = get_codon_params(codon)
        phi = params['phi']
        
        # Геометрия сворачивания (Time-Geometry)
        x += 1.5 * np.cos(np.radians(phi))
        y += 1.5 * np.sin(np.radians(phi))
        z += 0.8 
        
        # Формируем атом CA (Carbon Alpha)
        pdb_str += f"ATOM  {atom_id:5d}  CA  LYS A{res_id:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 {params['delay']:5.2f}           C\n"
        atom_id += 1
        res_id += 1
    
    # Соединяем атомы связями (CONECT) для визуализации
    for i in range(1, atom_id - 1):
        pdb_str += f"CONECT{i:5d}{i+1:5d}\n"
        
    return pdb_str

# --- ВИЗУАЛИЗАЦИЯ (STMOL) ---
with col2:
    if sequence and len(sequence) % 3 == 0:
        pdb = generate_structure(sequence)
        
        # Настройка вьюера
        view = py3dmol.view(width=600, height=400)
        view.addModel(pdb, 'pdb')
        
        # Объединенные стили (исправляем ошибку перезаписи)
        view.setStyle({
            'stick': {'radius': 0.15, 'color': 'lightgrey'},
            'sphere': {'scale': 0.25},
            'cartoon': {'color': 'spectrum'}
        })
        
        view.zoomTo()
        
        # Рендер через stmol (надежнее, чем raw html)
        showmol(view, height=400, width=600)
    else:
        st.info("Ожидание данных... Введите кодоны слева.")