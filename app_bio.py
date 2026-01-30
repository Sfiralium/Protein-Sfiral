import streamlit as st
import numpy as np
import py3Dmol  # Важно: Большая D
from stmol import showmol
import io
import re  # Для жесткой очистки мусора

# --- НАСТРОЙКИ СТРАНИЦЫ ---
st.set_page_config(page_title="Pico-Technology: Kushelev Lab", layout="wide", page_icon="🧬")

# --- СТИЛИ ---
st.markdown("""
<style>
    .stApp {background-color: #000000; color: #fff;}
    h1 {color: #00CCFF;}
    /* Скрываем лишние отступы */
    .block-container {padding-top: 2rem;}
</style>
""", unsafe_allow_html=True)

# --- ИНИЦИАЛИЗАЦИЯ (Reset) ---
if 'uploader_key' not in st.session_state: st.session_state.uploader_key = 0
def reset_app():
    st.session_state.uploader_key += 1
    st.experimental_rerun()

# --- ШАПКА ---
c1, c2 = st.columns([4, 1])
with c1:
    st.title("🧬 Pico-Technology")
    st.caption("Лаборатория Кушелева: Визуализация геометрической аллотропии")
with c2:
    if st.button("🗑 Сбросить"): reset_app()

# --- ПАРАМЕТРЫ КУШЕЛЕВА ---
KUSHELEV_TABLE = {
    'AAA': {'phi': -65, 'color': 'blue', 'bfactor': 0},   # Синий (Fast)
    'AAG': {'phi': -57, 'color': 'red', 'bfactor': 100},  # Красный (Slow)
    'DEFAULT': {'phi': -60, 'color': 'grey', 'bfactor': 50}
}
def get_params(codon): return KUSHELEV_TABLE.get(codon, KUSHELEV_TABLE['DEFAULT'])

# --- ИНТЕРФЕЙС ---
col_input, col_view = st.columns([1, 2])

with col_input:
    st.info("📂 Данные")
    
    # 1. Загрузка
    uploaded_file = st.file_uploader("Файл .fasta / .txt", type=["txt", "fasta"], key=f"up_{st.session_state.uploader_key}")
    dna_text = st.text_area("Или текст:", height=100, placeholder="AAA AAA AAG AAG...")

    # 2. Обработка (Считывание)
    raw_data = ""
    if uploaded_file:
        stringio = io.StringIO(uploaded_file.getvalue().decode("utf-8"))
        raw_data = stringio.read()
    elif dna_text:
        raw_data = dna_text

    # 3. ЖЕСТКАЯ ОЧИСТКА (RegEx)
    # Удаляем всё, что НЕ буквы A-Z (удаляются пробелы, \n, \r, цифры, >заголовки)
    # Если это FASTA, сначала уберем первую строку
    if ">" in raw_data:
        raw_data = re.sub(r'^>.*\n', '', raw_data, flags=re.MULTILINE)
    
    # Оставляем только буквы
    clean_seq = re.sub(r'[^a-zA-Z]', '', raw_data).upper()

    # 4. Обрезка под кратность 3
    final_seq = clean_seq
    if len(clean_seq) > 0:
        rem = len(clean_seq) % 3
        if rem != 0:
            final_seq = clean_seq[:-rem] # Молча отрезаем лишнее
            st.warning(f"✂️ Отрезано {rem} лишних символов (было {len(clean_seq)}).")
        else:
            st.success(f"✅ Данные приняты. Кодонов: {len(final_seq)//3}")

# --- ГЕНЕРАЦИЯ PDB ---
def make_pdb(seq):
    codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
    pdb = ""
    atom_id = 1
    res_id = 1
    x, y, z = 0.0, 0.0, 0.0
    
    for codon in codons:
        p = get_params(codon)
        rad = np.radians(p['phi'])
        
        # Спираль
        x += 2.0 * np.cos(rad)
        y += 2.0 * np.sin(rad)
        z += 0.8
        
        # Строгое форматирование PDB (Column-aligned)
        # ATOM      1  CA  LYS A   1      12.345  12.345  12.345  1.00 50.00           C
        pdb += f"ATOM  {atom_id:5d}  CA  LYS A{res_id:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 {p['bfactor']:5.2f}           C\n"
        
        atom_id += 1
        res_id += 1
        
    for i in range(1, atom_id - 1):
        pdb += f"CONECT{i:5d}{i+1:5d}\n"
        
    return pdb

# --- ВИЗУАЛИЗАЦИЯ ---
with col_view:
    if len(final_seq) >= 3:
        st.subheader("🧬 3D Модель")
        
        # Считаем статистику
        c1, c2 = st.columns(2)
        c1.metric("AAA (Синий)", final_seq.count('AAA'))
        c2.metric("AAG (Красный)", final_seq.count('AAG'))
        
        # Генерируем PDB
        pdb_data = make_pdb(final_seq)
        
        # Рисуем
        view = py3Dmol.view(width=800, height=600)
        view.addModel(pdb_data, 'pdb')
        
        # СТИЛЬ: Stick + Sphere (Самый надежный)
        view.setStyle({'stick': {'radius': 0.2, 'colorscheme': {'prop': 'b', 'gradient': 'bwr', 'min': 0, 'max': 100}}})
        view.addStyle({'sphere': {'scale': 0.4, 'colorscheme': {'prop': 'b', 'gradient': 'bwr', 'min': 0, 'max': 100}}})
        
        view.setBackgroundColor('#000000')
        view.zoomTo()
        showmol(view, height=600, width=800)
        
    else:
        st.info("👈 Загрузите файл, чтобы увидеть магию.")