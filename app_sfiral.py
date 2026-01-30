import streamlit as st
import numpy as np
import py3Dmol
from stmol import showmol
import io
import re

# --- КОНФИГУРАЦИЯ СФИРАЛИ ---
st.set_page_config(page_title="Sfiral Topology Core", layout="wide", page_icon="☯️")

st.markdown("""
<style>
    .stApp {background-color: #0e1117; color: #fff;}
    div.stButton > button:first-child {background-color: #FFD700; color: black; border: none; height: 50px; font-weight: bold;}
    .metric-card {background-color: #222; padding: 10px; border-radius: 5px; border: 1px solid #444;}
</style>
""", unsafe_allow_html=True)

if 'key' not in st.session_state: st.session_state.key = 0
def reset(): st.session_state.key += 1; st.rerun()

# --- МАТЕМАТИЧЕСКОЕ ЯДРО (ЗОЛОТЫЕ УГЛЫ) ---
# Мы используем те углы, которые доказали эффективность на Инсулине
BASE_ANGLES = {'z': 97.8, 'y': -55.97, 'x': -29.7}
MOVE_VEC = np.array([-0.8, 2.15, -1.37])

# Настройки состояний (Alpha, Beta, Pi, 3-10)
STATES = {
    'alpha': {'x': -44, 'y': 147, 'z': -16}, # Красный (Виток)
    'pi':    {'x': 10,  'y': 10,  'z': 0},   # Синий (Виток)
    'beta':  {'x': 114, 'y': -21, 'z': -17}, # Зеленый (S-ПЕТЛЯ / ИНВЕРСИЯ)
    '310':   {'x': 45,  'y': 0,   'z': 0}
}

# Матрицы поворота
def rot_x(deg):
    r = np.radians(deg); c, s = np.cos(r), np.sin(r)
    return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])
def rot_y(deg):
    r = np.radians(deg); c, s = np.cos(r), np.sin(r)
    return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])
def rot_z(deg):
    r = np.radians(deg); c, s = np.cos(r), np.sin(r)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])

# Типизация кодонов
def get_type(codon):
    c = codon.upper().replace("U", "T")
    if c.endswith('T'): return 'pi'
    if c.endswith('C'): return 'alpha'
    if c.endswith('A'): return 'pi' if c in ['TTA','CAA','AAA','GAA'] else 'beta'
    if c.endswith('G'): return 'alpha' if c in ['TTG','ATG','TGG'] else '310'
    return 'alpha'

# --- ИНТЕРФЕЙС ---
c1, c2 = st.columns([1, 3])

with c1:
    st.title("☯️ SFIRAL CORE")
    st.caption("Топологический анализ Свития")
    
    # Инсулин по умолчанию
    insulin_cds = "ATGGCCCTGTGGATGCGCCTCCTGCCCCTGCTGGCGCTGCTGGCCCTCTGGGGACCTGACCCAGCCGCAGCCTTTGTGAACCAACACCTGTGCGGCTCACACCTGGTGGAAGCTCTCTACCTAGTGTGCGGGGAACGAGGCTTCTTCTACACACCCAAGACCCGCCGGGAGGCAGAGGACCTGCAGGTGGGGCAGGTGGAGCTGGGCGGGGGCCCTGGTGCAGGCAGCCTGCAGCCCTTGGCCCTGGAGGGGTCCCTGCAGAAGCGTGGCATTGTGGAACAATGCTGTACCAGCATCTGCTCCCTCTACCAGCTGGAGAACTACTGCAACTAG"
    text_input = st.text_area("Ввод CDS:", value=insulin_cds, height=150)
    
    clean = re.sub(r'[^a-zA-Z]', '', text_input).upper().replace("U", "T")
    
    st.markdown("---")
    st.markdown("**Легенда Сфирали:**")
    st.markdown("🟡 **Правый виток (Yang):** Накопление")
    st.markdown("⚪ **S-Петля (Dao):** Точка Инверсии")
    st.markdown("🟣 **Левый виток (Yin):** Компенсация")
    
    run = st.button("🧬 АНАЛИЗ СВИТИЯ")

# --- ЛОГИКА СФИРАЛИ ---
def analyze_sfiral(seq):
    codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
    
    current_pos = np.array([0.0, 0.0, 0.0])
    current_rot = np.identity(3)
    
    points = []
    
    # Базовая матрица
    R_base = rot_z(BASE_ANGLES['z']) @ rot_y(BASE_ANGLES['y']) @ rot_x(BASE_ANGLES['x'])
    
    # Вектор начала (для расчета компенсации)
    start_pos = np.array([0.0, 0.0, 0.0])
    
    for i, codon in enumerate(codons):
        t = get_type(codon)
        s = STATES[t]
        
        # Расчет поворота
        R_codon = rot_y(s['y']) @ rot_x(s['x']) @ rot_z(s['z'])
        R_step = R_base @ R_codon
        current_rot = current_rot @ R_step
        step = current_rot @ MOVE_VEC
        current_pos = current_pos + step
        
        # --- ОПРЕДЕЛЕНИЕ ФАЗЫ СФИРАЛИ ---
        # Мы считаем, что 'beta' (зеленые углы) - это те самые S-петли (Инверторы)
        # Остальные (alpha, pi) - это Витки (накопители)
        
        sfiral_phase = "turn"
        color = "#FFD700" # Gold (Правый/Базовый)
        
        if t == 'beta':
            sfiral_phase = "s-loop"
            color = "#FFFFFF" # White (Инверсия)
        elif t == 'alpha':
            color = "#FFD700" # Gold
        elif t == 'pi':
            color = "#00FFFF" # Cyan (Можно трактовать как другой виток)
        else:
            color = "#888888"

        points.append({
            'x': current_pos[0], 'y': current_pos[1], 'z': current_pos[2],
            'color': color,
            'type': t,
            'phase': sfiral_phase
        })
        
    return points, current_pos

# --- ВИЗУАЛИЗАЦИЯ ---
with c2:
    if run and len(clean) > 0:
        points, end_pos = analyze_sfiral(clean)
        
        # Расчет "Коэффициента Свития" (Compensation Index)
        # Насколько конец близок к началу относительно длины цепи?
        total_length = len(points) * np.linalg.norm(MOVE_VEC)
        displacement = np.linalg.norm(end_pos) # Расстояние от (0,0,0) до конца
        
        # Коэффициент Сфиральности (100% = идеальное кольцо/нуль)
        sfiral_score = max(0, 100 - (displacement / total_length * 100))
        
        # МЕТРИКИ
        m1, m2, m3 = st.columns(3)
        m1.metric("Длина цепи", f"{len(points)} звеньев")
        m2.metric("Смещение (Displacement)", f"{displacement:.1f}")
        m3.metric("КОЭФФИЦИЕНТ СВИТИЯ", f"{sfiral_score:.1f}%")
        
        if sfiral_score > 80:
            st.success("✅ ВЫСОКАЯ СТЕПЕНЬ КОМПЕНСАЦИИ (СФИРАЛЬ ЗАМКНУТА)")
        elif sfiral_score > 50:
            st.warning("⚠️ ЧАСТИЧНАЯ КОМПЕНСАЦИЯ")
        else:
            st.error("❌ РАЗОМКНУТАЯ СПИРАЛЬ (НЕТ СВИТИЯ)")
            
        # 3D Viewer
        view = py3Dmol.view(width=900, height=700)
        
        prev = None
        for pt in points:
            # S-петли рисуем КРУПНЕЕ и БЕЛЫМ
            radius = 1.0 if pt['phase'] == 's-loop' else 0.6
            
            view.addSphere({'center': pt, 'radius': radius, 'color': pt['color']})
            if prev:
                # Цвет связи - градиент или цвет текущего
                view.addCylinder({'start': prev, 'end': pt, 'radius': 0.2, 'color': pt['color']})
            prev = pt
            
        # Рисуем вектор компенсации (от конца к началу)
        view.addArrow({
            'start': {'x': end_pos[0], 'y': end_pos[1], 'z': end_pos[2]},
            'end': {'x': 0, 'y': 0, 'z': 0},
            'color': 'red',
            'radius': 0.3
        })
        
        view.setBackgroundColor('#0e1117')
        view.zoomTo()
        showmol(view, height=700, width=900)
        
        st.info("ℹ️ Красная стрелка показывает вектор, необходимый для полного обнуления (возврата в Великий Предел).")