import streamlit as st
import subprocess
import sys
import time

# --- 1. НАСТРОЙКИ (Должны быть в самом верху) ---
st.set_page_config(page_title="Sfiral Protein Lab", layout="wide", page_icon="🧬")

# --- 2. АВТО-УСТАНОВЩИК (AGRESSIVE INSTALLER) ---
try:
    import py3dmol
    from Bio.Seq import Seq
except ImportError:
    st.warning("⚙️ Обнаружен сбой облака. Запускаю принудительную установку библиотек... (10-20 сек)")
    # Устанавливаем тихо, чтобы не пугать пользователя
    subprocess.check_call([sys.executable, "-m", "pip", "install", "py3dmol", "biopython", "numpy"])
    st.success("✅ Готово! Перезагружаю систему...")
    time.sleep(1)
    st.rerun() # <--- ВОТ ЭТА КОМАНДА СПАСЕТ СИТУАЦИЮ

# Если мы здесь, значит библиотеки установлены. Подключаем их.
import py3dmol
from Bio.Seq import Seq
import numpy as np

# --- 3. СТИЛИ И ИНТЕРФЕЙС ---
st.markdown("""
<style>
    .stApp {background-color: #0e1117; color: #fff;}
    h1 {color: #00CCFF;}
    .report {background: #161b22; padding: 15px; border-radius: 10px; border: 1px solid #30363d;}
</style>
""", unsafe_allow_html=True)

st.title("🧬 Protein-Sfiral: Time-Genetics Folding")
st.caption("Testing the Kushelev Hypothesis: Same Amino Acids -> Different Geometry (CDS-driven)")

# --- 4. ЛОГИКА ---
TIME_GENETICS_MAP = {
    'AAA': {'aa': 'K', 'phi': -65, 'psi': -40, 'delay': 1.0, 'note': 'Fast (Pi-Helix)'},
    'AAG': {'aa': 'K', 'phi': -57, 'psi': -47, 'delay': 1.5, 'note': 'Slow (Alpha-Helix)'},
    'DEFAULT': {'aa': '?', 'phi': -60, 'psi': -45, 'delay': 1.0}
}

def get_codon_params(codon):
    return TIME_GENETICS_MAP.get(codon, TIME_GENETICS_MAP['DEFAULT'])

col1, col2 = st.columns([1, 2])

with col1:
    st.subheader("📥 Ввод Последовательности (CDS)")
    st.info("Загрузите ДНК-последовательность (Нуклеотиды: A, T, G, C)")
    
    uploaded_file = st.file_uploader("Перетащите файл сюда (.txt, .fasta)", type=["txt", "fasta"])
    dna_input = st.text_area("Или введите вручную:", height=150, placeholder="Например: AAA AAA AAA AAG AAG AAG...")

    sequence = ""
    if uploaded_file is not None:
        stringio = uploaded_file.getvalue().decode("utf-8")
        sequence = stringio.replace("\n", "").replace(" ", "").upper()
    elif dna_input:
        sequence = dna_input.replace("\n", "").replace(" ", "").upper()

    if sequence:
        if len(sequence) % 3 != 0:
            st.error(f"⚠ Длина ДНК ({len(sequence)}) не делится на 3!")
        else:
            st.success(f"✅ Загружено {len(sequence)//3} кодонов.")

def generate_structure_from_time(dna_seq):
    codons = [dna_seq[i:i+3] for i in range(0, len(dna_seq), 3)]
    pdb_str = ""
    atom_id = 1
    res_id = 1
    x, y, z = 0.0, 0.0, 0.0
    
    for codon in codons:
        params = get_codon_params(codon)
        aa_name = "LYS" if params['aa'] == 'K' else "ALA"
        phi = params['phi']
        x += 1.5 * np.cos(np.radians(phi))
        y += 1.5 * np.sin(np.radians(phi))
        z += 0.8 
        pdb_str += f"ATOM  {atom_id:5d}  CA  {aa_name} A{res_id:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 {params['delay']:5.2f}           C\n"
        atom_id += 1
        res_id += 1
    return pdb_str, codons

with col2:
    st.subheader("🧬 3D-Симуляция Структуры")
    if sequence and len(sequence) % 3 == 0:
        pdb_data, parsed_codons = generate_structure_from_time(sequence)
        aaa_count = parsed_codons.count('AAA')
        aag_count = parsed_codons.count('AAG')
        st.write(f"**Анализ состава:** AAA (Fast): {aaa_count} | AAG (Slow): {aag_count}")
        if aag_count > 0 and aaa_count > 0:
            st.warning("Обнаружена программная аллотропия! Белок имеет разные фазовые состояния.")
        
        view = py3dmol.view(width=600, height=400)
        view.addModel(pdb_data, 'pdb')
        view.setStyle({'stick': {'radius': 0.2}, 'sphere': {'scale': 0.3}})
        view.setStyle({'cartoon': {'colorscheme': {'prop': 'b', 'gradient': 'roygb', 'min': 1, 'max': 1.5}}})
        view.zoomTo()
        output = view._make_html()
        st.components.v1.html(output, width=600, height=400)
    else:
        st.info("Ожидание данных...")
