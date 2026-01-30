import streamlit as st
import numpy as np
import py3Dmol
from stmol import showmol
import io
import re

# --- НАСТРОЙКИ ---
st.set_page_config(page_title="Helix Engine", layout="wide")
st.markdown("""
<style>
    .stApp {background-color: #000000; color: #fff;}
    div.stButton > button:first-child {
        background-color: #00FF00; color: black; font-weight: bold; border: none; height: 50px; width: 100%;
    }
</style>
""", unsafe_allow_html=True)

# --- СБРОС ---
if 'key' not in st.session_state: st.session_state.key = 0
def reset(): st.session_state.key += 1; st.rerun()

# --- ПАРАМЕТРЫ (ГЕОМЕТРИЯ) ---
KUSHELEV_TABLE = {
    'AAA': {'phi': -65, 'radius': 2.0, 'color': '#00FFFF'}, # Узкий (Cyan)
    'AAG': {'phi': -45, 'radius': 6.0, 'color': '#FF00FF'}, # Широкий (Magenta)
    'DEFAULT': {'phi': -60, 'radius': 4.0, 'color': '#FFFFFF'}
}
def get_params(c): return KUSHELEV_TABLE.get(c, KUSHELEV_TABLE['DEFAULT'])

# --- ИНТЕРФЕЙС ---
c1, c2 = st.columns([1, 3])

with c1:
    st.title("🧬 Helix Engine")
    st.caption("Visualizing Geometric Allotropy")
    if st.button("🗑 СБРОС"): reset()
    
    # ПРИВЯЗЫВАЕМ КЛЮЧИ (KEY) К СЧЕТЧИКУ СБРОСА
    uploaded_file = st.file_uploader("Файл", key=f"up_{st.session_state.key}")
    
    # ИСПРАВЛЕНО: Теперь текстовое поле тоже имеет динамический ключ
    text_input = st.text_area("Или текст", height=150, key=f"text_{st.session_state.key}")
    
    # Чтение и Очистка
    raw = ""
    if uploaded_file:
        raw = io.StringIO(uploaded_file.getvalue().decode("utf-8")).read()
    elif text_input:
        raw = text_input

    if ">" in raw: raw = re.sub(r'^>.*\n', '', raw, flags=re.MULTILINE)
    clean = re.sub(r'[^a-zA-Z]', '', raw).upper()
    
    # Валидация
    final_seq = clean
    if len(clean) > 0:
        rem = len(clean) % 3
        if rem != 0:
            final_seq = clean[:-rem]
            st.warning(f"✂️ Отрезано {rem} симв.")
        st.success(f"Кодонов: {len(final_seq)//3}")
        
    run = st.button("🚀 ЗАПУСТИТЬ")

# --- ГЕНЕРАЦИЯ ---
def get_coords(seq):
    codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
    points = []
    
    x, y, z = 0.0, 0.0, 0.0
    current_angle = 0.0 
    z_step = 0.8
    
    pdb = ""
    atom_id = 1
    
    for i, codon in enumerate(codons):
        p = get_params(codon)
        
        current_angle += p['phi']
        rad = np.radians(current_angle)
        r = p['radius']
        
        x = r * np.cos(rad)
        y = r * np.sin(rad)
        z += z_step
        
        points.append({'x': x, 'y': y, 'z': z, 'color': p['color']})
        
        pdb += f"ATOM  {atom_id:5d}  CA  LYS A{i+1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 50.00           C\n"
        atom_id += 1
        
    return points, pdb

# --- ВИЗУАЛИЗАЦИЯ ---
with c2:
    if run and len(final_seq) >= 3:
        points, pdb_data = get_coords(final_seq)
        
        view = py3Dmol.view(width=800, height=700)
        
        prev = None
        for pt in points:
            view.addSphere({'center': pt, 'radius': 0.8, 'color': pt['color']})
            if prev:
                view.addCylinder({'start': prev, 'end': pt, 'radius': 0.3, 'color': pt['color']})
            prev = pt
            
        view.setBackgroundColor('#111111')
        view.zoomTo()
        showmol(view, height=700, width=800)
        st.download_button("💾 PDB", pdb_data, "helix.pdb")