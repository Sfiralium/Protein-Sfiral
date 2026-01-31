import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import py3Dmol
from stmol import showmol
import re

# --- КОНФИГУРАЦИЯ ---
st.set_page_config(page_title="Sfiral Lab v1.0", layout="wide", page_icon="☯️")

st.markdown("""
<style>
    .stApp {background-color: #0e1117; color: #ddd;}
    div.stButton > button:first-child {background-color: #D4AF37; color: black; border: none; height: 50px; font-weight: bold; font-size: 18px;}
    h1, h2, h3 {color: #D4AF37;}
    .stTabs [data-baseweb="tab-list"] {gap: 24px;}
    .stTabs [data-baseweb="tab"] {height: 50px; white-space: pre-wrap; background-color: #222; border-radius: 4px 4px 0 0; color: #fff;}
    .stTabs [aria-selected="true"] {background-color: #D4AF37; color: #000;}
</style>
""", unsafe_allow_html=True)

# --- ГЕОМЕТРИЧЕСКОЕ ЯДРО (CONSTANTS) ---
BASE_ANGLES = {'z': 97.8, 'y': -55.97, 'x': -29.7}
MOVE_VEC = np.array([-0.8, 2.15, -1.37])

# "Золотые Углы" Сфирали (доказанные на Инсулине)
STATES = {
    'alpha': {'x': -44, 'y': 147, 'z': -16}, # Ян (Накопление)
    'pi':    {'x': 10,  'y': 10,  'z': 0},   # Ян (Вариация)
    'beta':  {'x': 114, 'y': -21, 'z': -17}, # Дао (S-Петля / Переход)
    '310':   {'x': 45,  'y': 0,   'z': 0}    # Инь (Сжатие)
}

# --- МАТЕМАТИКА ---
def rot_x(deg): r=np.radians(deg); c,s=np.cos(r),np.sin(r); return np.array([[1,0,0],[0,c,-s],[0,s,c]])
def rot_y(deg): r=np.radians(deg); c,s=np.cos(r),np.sin(r); return np.array([[c,0,s],[0,1,0],[-s,0,c]])
def rot_z(deg): r=np.radians(deg); c,s=np.cos(r),np.sin(r); return np.array([[c,-s,0],[s,c,0],[0,0,1]])

def get_type(codon):
    c = codon.upper().replace("U", "T")
    if c.endswith('T'): return 'pi'
    if c.endswith('C'): return 'alpha'
    if c.endswith('A'): return 'pi' if c in ['TTA','CAA','AAA','GAA'] else 'beta'
    if c.endswith('G'): return 'alpha' if c in ['TTG','ATG','TGG'] else '310'
    return 'alpha'

def build_sfiral_model(seq):
    codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
    pos = np.array([0.0, 0.0, 0.0])
    rot = np.identity(3)
    
    points = []
    R_base = rot_z(BASE_ANGLES['z']) @ rot_y(BASE_ANGLES['y']) @ rot_x(BASE_ANGLES['x'])
    
    atom_id = 1
    pdb_str = ""
    
    # 1. Генерация точек
    for i, codon in enumerate(codons):
        t = get_type(codon)
        s = STATES[t]
        
        rot = rot @ (R_base @ (rot_y(s['y']) @ rot_x(s['x']) @ rot_z(s['z'])))
        pos = pos + (rot @ MOVE_VEC)
        
        # Определение фазы
        phase = "coil"
        color = "#FFD700" # Золото (по умолчанию)
        radius = 0.6
        pdb_factor = 1.0
        
        if t == 'beta':
            phase = "s-loop"
            color = "#FFFFFF" # Белый (Инверсия)
            radius = 1.1
            pdb_factor = 99.0
        elif t == 'alpha':
            color = "#FFD700"
        elif t == 'pi':
            color = "#00FFFF" # Циан
        elif t == '310':
            color = "#FF00FF"
            
        points.append({
            'pos': pos, 'color': color, 'r': radius, 'phase': phase, 'idx': i
        })
        
        pdb_str += f"ATOM  {atom_id:5d}  CA  UNK A{i+1:4d}    {pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}  1.00{pdb_factor:6.2f}           C\n"
        atom_id += 1
        
    return points, pdb_str

# --- ИНТЕРФЕЙС ---
st.title("☯️ SFIRAL LABORATORY")
st.caption("Единая система анализа Времени и Пространства")

# Сайдбар (Ввод данных)
with st.sidebar:
    st.header("🧬 Ввод Данных")
    insulin_def = "ATGGCCCTGTGGATGCGCCTCCTGCCCCTGCTGGCGCTGCTGGCCCTCTGGGGACCTGACCCAGCCGCAGCCTTTGTGAACCAACACCTGTGCGGCTCACACCTGGTGGAAGCTCTCTACCTAGTGTGCGGGGAACGAGGCTTCTTCTACACACCCAAGACCCGCCGGGAGGCAGAGGACCTGCAGGTGGGGCAGGTGGAGCTGGGCGGGGGCCCTGGTGCAGGCAGCCTGCAGCCCTTGGCCCTGGAGGGGTCCCTGCAGAAGCGTGGCATTGTGGAACAATGCTGTACCAGCATCTGCTCCCTCTACCAGCTGGAGAACTACTGCAACTAG"
    txt = st.text_area("CDS Код:", value=insulin_def, height=200)
    clean_seq = re.sub(r'[^a-zA-Z]', '', txt).upper().replace("U", "T")
    st.info(f"Длина: {len(clean_seq)//3} кодонов")
    
    run = st.button("ЗАПУСТИТЬ АНАЛИЗ")

# --- ОСНОВНАЯ ЛОГИКА ---
if run and len(clean_seq) > 0:
    points, pdb_data = build_sfiral_model(clean_seq)
    
    # Вычисления для графиков
    coords = np.array([p['pos'] for p in points])
    center = np.mean(coords, axis=0)
    distances = [np.linalg.norm(p['pos'] - center) for p in points]
    s_loops = [p['idx'] for p in points if p['phase'] == 's-loop']
    
    # ТАБЫ
    tab1, tab2 = st.tabs(["🧬 Сфиральное Свитие (3D)", "⏳ График Ритма (2D)"])
    
    # --- ТАБ 1: 3D СВИТИЕ ---
    with tab1:
        c1, c2 = st.columns([3, 1])
        with c1:
            view = py3Dmol.view(width=800, height=600)
            prev = None
            for p in points:
                view.addSphere({'center':{'x':p['pos'][0],'y':p['pos'][1],'z':p['pos'][2]}, 'radius':p['r'], 'color':p['color']})
                if prev:
                    lnk_col = "#FFF" if (p['phase']=='s-loop' or prev['phase']=='s-loop') else p['color']
                    view.addCylinder({'start':{'x':prev['pos'][0],'y':prev['pos'][1],'z':prev['pos'][2]}, 
                                      'end':{'x':p['pos'][0],'y':p['pos'][1],'z':p['pos'][2]}, 'radius':0.2, 'color':lnk_col})
                prev = p
            
            # Вектор ошибки (Красный)
            end_pos = points[-1]['pos']
            view.addArrow({'start': {'x':end_pos[0], 'y':end_pos[1], 'z':end_pos[2]}, 
                           'end': {'x':0, 'y':0, 'z':0}, 'color': 'red', 'radius': 0.3})
            
            view.setBackgroundColor('#0e1117')
            view.zoomTo()
            showmol(view, height=600, width=800)
            
        with c2:
            displacement = np.linalg.norm(end_pos)
            st.metric("Вектор Ошибки", f"{displacement:.1f} Å")
            st.metric("Кол-во S-петель", len(s_loops))
            
            st.download_button("💾 Скачать PDB (Паспорт)", pdb_data, "sfiral_passport.pdb")
            st.markdown("""
            **Легенда:**
            * ⚪ **Белый:** Точка S-петли (Инверсия)
            * 🟡 **Золотой:** Виток Накопления
            * 🔴 **Красная стрелка:** Остаточный вектор
            """)

    # --- ТАБ 2: ГРАФИК РИТМА ---
    with tab2:
        fig, ax = plt.subplots(figsize=(12, 6))
        fig.patch.set_facecolor('#0e1117')
        ax.set_facecolor('#0e1117')
        
        # Линия
        ax.plot(distances, color='#FFD700', linewidth=2, label='Амплитуда')
        
        # Маркеры S-петель
        for l in s_loops:
            ax.axvline(x=l, color='white', linestyle='--', alpha=0.5)
            
        ax.spines['bottom'].set_color('#fff')
        ax.spines['left'].set_color('#fff')
        ax.tick_params(colors='white')
        ax.set_xlabel('Время (Такты)', color='white')
        ax.set_ylabel('Удаление от Центра', color='white')
        ax.set_title(f'Хронограмма Сфирали: {len(s_loops)} переходов', color='#D4AF37')
        
        st.pyplot(fig)
        
        if len(s_loops) > 1:
            avg_period = np.mean(np.diff(s_loops))
            st.info(f"Средний Период Ритма: **{avg_period:.1f}** тактов")
        else:
            st.warning("Ритм линейный (нет или мало точек перехода)")