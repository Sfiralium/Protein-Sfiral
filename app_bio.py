import streamlit as st
import numpy as np
import py3Dmol
from stmol import showmol
import io
import re

st.set_page_config(page_title="Helix Engine PRO", layout="wide", page_icon="🧬")

# --- СТИЛИ ---
st.markdown("""
<style>
    .stApp {background-color: #000000; color: #fff;}
    div.stButton > button:first-child {background-color: #00FF00; color: black; border: none; height: 50px; width: 100%; font-size: 20px;}
    .stSidebar {background-color: #111;}
</style>
""", unsafe_allow_html=True)

if 'key' not in st.session_state: st.session_state.key = 0
def reset(): st.session_state.key += 1; st.rerun()

# --- МАТРИЦЫ ---
def rot_x(deg):
    r = np.radians(deg)
    c, s = np.cos(r), np.sin(r)
    return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])
def rot_y(deg):
    r = np.radians(deg)
    c, s = np.cos(r), np.sin(r)
    return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])
def rot_z(deg):
    r = np.radians(deg)
    c, s = np.cos(r), np.sin(r)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])

# --- ТИПЫ СПИРАЛЕЙ ---
def get_type(codon):
    c = codon.upper().replace("U", "T")
    if c.endswith('T'): return 'pi'
    if c.endswith('C'): return 'alpha'
    if c.endswith('A'): return 'pi' if c in ['TTA','CAA','AAA','GAA'] else 'beta'
    if c.endswith('G'): return 'alpha' if c in ['TTG','ATG','TGG'] else '310'
    return 'alpha'

# --- ИНТЕРФЕЙС ---
with st.sidebar:
    st.title("🧬 Настройка 4-х форм")
    
    # БАЗОВЫЕ НАСТРОЙКИ (СКРЫТЫ, НО РАБОТАЮТ)
    base_z = 97.8
    base_y = -55.97
    base_x = -29.7
    move_vec = np.array([-0.8, 2.15, -1.37])

    states = {}
    
    # 🔴 ALPHA: Установлены ваши "Золотые углы" (-45, 150, -16)
    with st.expander("🔴 ALPHA (Red / Standard)", expanded=True):
        states['alpha'] = {
            'x': st.slider("Alpha X", -180, 180, -45, key="ax"),
            'y': st.slider("Alpha Y", -180, 180, 150, key="ay"),
            'z': st.slider("Alpha Z", -180, 180, -16, key="az"),
            'color': '#FF0000'
        }
        
    # 🔵 PI: 
    with st.expander("🔵 PI (Blue / Wide)"):
        states['pi'] = {
            'x': st.slider("Pi X", -180, 180, 10, key="px"),
            'y': st.slider("Pi Y", -180, 180, 10, key="py"),
            'z': st.slider("Pi Z", -180, 180, 0, key="pz"),
            'color': '#00FFFF'
        }
        
    # 🟢 BETA: Установлены ваши углы (114, -21, -17)
    with st.expander("🟢 BETA (Green / Flat)", expanded=True):
        states['beta'] = {
            'x': st.slider("Beta X", -180, 180, 114, key="bx"),
            'y': st.slider("Beta Y", -180, 180, -21, key="by"),
            'z': st.slider("Beta Z", -180, 180, -17, key="bz"),
            'color': '#00FF00'
        }
        
    # 🟡 3-10:
    with st.expander("🟡 3-10 (Yellow / Tight)"):
        states['310'] = {
            'x': st.slider("3-10 X", -180, 180, 45, key="tx"),
            'y': st.slider("3-10 Y", -180, 180, 0, key="ty"),
            'z': st.slider("3-10 Z", -180, 180, 0, key="tz"),
            'color': '#FFFF00'
        }

st.title("🧬 Helix Engine PRO")
st.caption("Mode: Real Gene Analysis")

c1, c2 = st.columns([1, 4])
with c1:
    if st.button("🗑 СБРОС"): reset()

# Гемоглобин по умолчанию (для демонстрации)
hemoglobin_cds = "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCACTGTGACAAGCTGCACGTGGATCCTGAGAACTTCAGGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGCAAAGAATTCACCCCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAATGCCCTGGCCCACAAGTATCACTAA"

uploaded_file = st.file_uploader("Файл гена (.fasta)", key=f"up_{st.session_state.key}")
text_input = st.text_area("Или код CDS:", value=hemoglobin_cds, height=150, key=f"txt_{st.session_state.key}")

# ЛОГИКА
raw = ""
if uploaded_file: raw = io.StringIO(uploaded_file.getvalue().decode("utf-8")).read()
elif text_input: raw = text_input
if ">" in raw: raw = re.sub(r'^>.*\n', '', raw, flags=re.MULTILINE)
clean = re.sub(r'[^a-zA-Z]', '', raw).upper().replace("U", "T")

if len(clean) > 0:
    st.success(f"Кодонов: {len(clean)//3}")

if st.button("🚀 ПОСТРОИТЬ ГЕН"):
    codons = [clean[i:i+3] for i in range(0, len(clean), 3)]
    
    current_pos = np.array([0.0, 0.0, 0.0])
    current_rot = np.identity(3)
    
    points = []
    pdb = ""
    atom_id = 1
    
    R_base = rot_z(base_z) @ rot_y(base_y) @ rot_x(base_x)
    
    for i, codon in enumerate(codons):
        t = get_type(codon)
        s = states[t]
        
        R_codon = rot_y(s['y']) @ rot_x(s['x']) @ rot_z(s['z'])
        R_step = R_base @ R_codon
        current_rot = current_rot @ R_step
        step = current_rot @ move_vec
        current_pos = current_pos + step
        
        points.append({'x': current_pos[0], 'y': current_pos[1], 'z': current_pos[2], 'color': s['color']})
        pdb += f"ATOM  {atom_id:5d}  CA  UNK A{i+1:4d}    {current_pos[0]:8.3f}{current_pos[1]:8.3f}{current_pos[2]:8.3f}  1.00 50.00           C\n"
        atom_id += 1

    view = py3Dmol.view(width=1000, height=800)
    prev = None
    for pt in points:
        view.addSphere({'center': pt, 'radius': 0.8, 'color': pt['color']})
        if prev: view.addCylinder({'start': prev, 'end': pt, 'radius': 0.3, 'color': 'white'})
        prev = pt
    
    view.setBackgroundColor('#000000')
    view.zoomTo()
    showmol(view, height=800, width=1000)
    
    c_dw, _ = st.columns([1, 4])
    with c_dw:
        st.download_button("💾 Скачать PDB", pdb, "helix_structure.pdb")