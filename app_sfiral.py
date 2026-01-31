import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import py3Dmol
from stmol import showmol
import re
import os
from PIL import Image

# --- НАСТРОЙКИ СТРАНИЦЫ ---
st.set_page_config(page_title="Sfiral Laboratory: Time-Genetics", layout="wide", page_icon="Sfiralium_Stop.ico")

# --- CSS СТИЛИЗАЦИЯ (PREMIUM GOLD) ---
st.markdown("""
<style>
    /* Основной фон */
    .stApp {background-color: #000000; color: #e0e0e0; font-family: 'Arial', sans-serif;}
    
    /* Увеличение шрифтов */
    p, label, span, div {font-size: 16px !important;}
    
    /* ОГРОМНЫЙ ЗАГОЛОВОК */
    .big-title {
        font-size: 60px !important;
        font-weight: bold;
        color: #D4AF37;
        margin-bottom: 0px;
        text-transform: uppercase;
        letter-spacing: 2px;
    }
    .subtitle {
        font-size: 24px !important;
        color: #888;
        margin-top: -10px;
        margin-bottom: 30px;
    }
    
    /* Золотые заголовки блоков */
    h1, h2, h3 {color: #D4AF37 !important;}
    
    /* Кнопка запуска */
    div.stButton > button:first-child {
        background-color: #D4AF37; 
        color: #000000; 
        border: none; 
        height: 65px; 
        width: 100%;
        font-weight: bold; 
        font-size: 24px !important;
        border-radius: 8px;
        transition: 0.3s;
        margin-top: 20px;
        text-transform: uppercase;
    }
    div.stButton > button:first-child:hover {
        background-color: #FFD700;
        box-shadow: 0 0 20px #D4AF37;
    }

    /* Табы */
    .stTabs [data-baseweb="tab-list"] {gap: 15px;}
    .stTabs [data-baseweb="tab"] {
        background-color: #1a1a1a; 
        border-radius: 5px 5px 0 0; 
        color: #aaa;
        border: 1px solid #333;
        font-size: 18px !important;
        padding: 15px 30px;
    }
    .stTabs [aria-selected="true"] {
        background-color: #D4AF37 !important; 
        color: #000 !important;
        font-weight: bold;
        border-bottom: none;
    }
    
    /* Инфо-боксы (Разъяснения) */
    .info-box {
        background-color: #0a0a0a; 
        border-left: 6px solid #D4AF37; 
        padding: 20px; 
        margin-bottom: 25px;
        border-radius: 5px;
        border: 1px solid #222;
        font-size: 18px;
        line-height: 1.6;
        box-shadow: 0 4px 6px rgba(0,0,0,0.3);
    }
    .info-title {
        color: #D4AF37;
        font-weight: bold;
        font-size: 20px;
        margin-bottom: 10px;
        display: block;
    }
</style>
""", unsafe_allow_html=True)

# --- ГЕОМЕТРИЯ ---
BASE_ANGLES = {'z': 97.8, 'y': -55.97, 'x': -29.7}
MOVE_VEC = np.array([-0.8, 2.15, -1.37])
PHI = 1.61803398875
STATES = {
    'alpha': {'x': -44, 'y': 147, 'z': -16}, 
    'pi':    {'x': 10,  'y': 10,  'z': 0},   
    'beta':  {'x': 114, 'y': -21, 'z': -17}, # S-Петля
    '310':   {'x': 45,  'y': 0,   'z': 0}
}

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
    current_polarity = 1 
    atom_id = 1
    pdb_str = ""
    
    for i, codon in enumerate(codons):
        t = get_type(codon)
        s = STATES[t]
        rot = rot @ (R_base @ (rot_y(s['y']) @ rot_x(s['x']) @ rot_z(s['z'])))
        pos = pos + (rot @ MOVE_VEC)
        
        phase = "coil"
        radius = 0.6
        pdb_bfactor = 1.0 
        color = "#FFD700" if current_polarity > 0 else "#00FFFF"
        
        if t == 'beta':
            phase = "s-loop"
            color = "#FFFFFF"
            radius = 1.2
            pdb_bfactor = 99.0 
            current_polarity *= -1 
            
        points.append({'pos': pos, 'color': color, 'r': radius, 'phase': phase, 'idx': i, 'polarity': current_polarity})
        pdb_str += f"ATOM  {atom_id:5d}  CA  UNK A{i+1:4d}    {pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}  1.00{pdb_bfactor:6.2f}           C\n"
        atom_id += 1
        
    return points, pdb_str

# --- ШАПКА: ЛОГО И НАЗВАНИЕ ---
c_logo, c_title = st.columns([1, 5])
with c_logo:
    # Ищем иконку
    if os.path.exists("Sfiralium_Stop.ico"):
        st.image("Sfiralium_Stop.ico", width=150)
    else:
        st.markdown("# ☯️")

with c_title:
    st.markdown('<div class="big-title">SFIRAL LABORATORY</div>', unsafe_allow_html=True)
    st.markdown('<div class="subtitle">Time-Genetics Analysis System</div>', unsafe_allow_html=True)

# --- ОБЩАЯ КОНЦЕПЦИЯ ---
st.markdown("""
<div class="info-box">
<span class="info-title">📚 О Концепции</span>
Система основана на фундаментальных трудах <b>О.С. Басаргина</b> (<i>«Сфираль времени», «Антисимметрия времени», «Волновой геном»</i>).<br><br>
<b>Времягенетика (Time-Genetics)</b> постулирует первичность Ритма (Времени) над Материей. Белок рассматривается не как случайный клубок химии, а как <b>«Сфираль»</b> — устройство, компенсирующее напряжение времени через зеркальную антисимметрию и S-образные переходы (Великий Предел).
</div>
""", unsafe_allow_html=True)

# --- САЙДБАР: ВВОД ДАННЫХ ---
with st.sidebar:
    st.header("🧬 Ввод Данных")
    
    # Загрузка
    uploaded_file = st.file_uploader("Загрузить файл CDS (.txt, .fasta)", type=['txt', 'fasta'])
    
    st.markdown("---")
    preset = st.selectbox("Или выбрать пример:", ["Инсулин (Ритм 7 - Эталон)", "Гемоглобин (Ритм 90 - Структура)", "Коллаген (Ритм 0 - Линия)", "Ввести вручную"])
    
    seq_input = ""
    insulin_seq = "ATGGCCCTGTGGATGCGCCTCCTGCCCCTGCTGGCGCTGCTGGCCCTCTGGGGACCTGACCCAGCCGCAGCCTTTGTGAACCAACACCTGTGCGGCTCACACCTGGTGGAAGCTCTCTACCTAGTGTGCGGGGAACGAGGCTTCTTCTACACACCCAAGACCCGCCGGGAGGCAGAGGACCTGCAGGTGGGGCAGGTGGAGCTGGGCGGGGGCCCTGGTGCAGGCAGCCTGCAGCCCTTGGCCCTGGAGGGGTCCCTGCAGAAGCGTGGCATTGTGGAACAATGCTGTACCAGCATCTGCTCCCTCTACCAGCTGGAGAACTACTGCAACTAG"
    hbb_seq = "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCACTGTGACAAGCTGCACGTGGATCCTGAGAACTTCAGGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGCAAAGAATTCACCCCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAATGCCCTGGCCCACAAGTATCACTAA"
    
    if uploaded_file:
        seq_input = uploaded_file.getvalue().decode("utf-8")
    elif preset == "Инсулин (Ритм 7 - Эталон)": seq_input = insulin_seq
    elif preset == "Гемоглобин (Ритм 90 - Структура)": seq_input = hbb_seq
    
    txt = st.text_area("CDS Код:", value=seq_input, height=150)
    clean_seq = re.sub(r'[^a-zA-Z]', '', txt).upper().replace("U", "T")
    
    if len(clean_seq) > 0:
        st.success(f"Длина цепи: {len(clean_seq)//3} кодонов")
        
    run_btn = st.button("ЗАПУСТИТЬ АНАЛИЗ")

# --- ЛОГИКА АНАЛИЗА ---
if run_btn and len(clean_seq) > 0:
    points, pdb_data = build_sfiral_model(clean_seq)
    
    # Расчеты метрик
    center = np.mean([p['pos'] for p in points], axis=0)
    dists = [np.linalg.norm(p['pos'] - center) for p in points]
    signed_dists = [d * p['polarity'] for d, p in zip(dists, points)]
    s_loops = [p['idx'] for p in points if p['phase'] == 's-loop']
    
    intervals = []
    if len(s_loops) > 0:
        prev = 0
        for l in s_loops:
            intervals.append(l - prev)
            prev = l
            
    # Вкладки
    t1, t2, t3 = st.tabs(["🧬 3D СВИТИЕ (ФОРМА)", "📉 АНТИСИММЕТРИЯ (РИТМ)", "🌀 ФРАКТАЛ (ЗОЛОТОЕ СЕЧЕНИЕ)"])
    
    # --- ТАБ 1: 3D ---
    with t1:
        st.markdown("""
        <div class="info-box">
        <span class="info-title">1. Геометрия Сфирали</span>
        Здесь мы видим, как Линейное Время (код ДНК) сворачивается в Пространственную Форму.
        <ul>
        <li><b>Белые узлы:</b> S-петли (Точки перехода).</li>
        <li><b>Красная стрелка:</b> Вектор ошибки (насколько система не скомпенсирована). В идеальной Сфирали начало и конец сходятся в ноль.</li>
        </ul>
        </div>
        """, unsafe_allow_html=True)
        
        c1, c2 = st.columns([3, 1])
        with c1:
            view = py3Dmol.view(width=900, height=700)
            prev = None
            for p in points:
                view.addSphere({'center':{'x':p['pos'][0],'y':p['pos'][1],'z':p['pos'][2]}, 'radius':p['r'], 'color':p['color']})
                if prev:
                    lnk = "#FFF" if (p['phase']=='s-loop' or prev['phase']=='s-loop') else p['color']
                    view.addCylinder({'start':{'x':prev['pos'][0],'y':prev['pos'][1],'z':prev['pos'][2]}, 'end':{'x':p['pos'][0],'y':p['pos'][1],'z':p['pos'][2]}, 'radius':0.2, 'color':lnk})
                prev = p
            # Стрелка ошибки
            end_pos = points[-1]['pos']
            view.addArrow({'start':{'x':end_pos[0],'y':end_pos[1],'z':end_pos[2]}, 'end':{'x':0,'y':0,'z':0}, 'color':'#FF0000', 'radius':0.3})
            view.setBackgroundColor('#000')
            view.zoomTo()
            showmol(view, height=700, width=900)
        with c2:
            st.metric("Кол-во S-петель", len(s_loops))
            err_vec = np.linalg.norm(end_pos)
            st.metric("Вектор Ошибки (Å)", f"{err_vec:.1f}")
            st.download_button("💾 СКАЧАТЬ PDB ПАСПОРТ", pdb_data, "sfiral_model.pdb")

    # --- ТАБ 2: АНТИСИММЕТРИЯ ---
    with t2:
        st.markdown("""
        <div class="info-box">
        <span class="info-title">2. Дыхание Антисимметрии</span>
        Сфираль — это не статика, а пульсация. График показывает, как система меняет полярность.
        <ul>
        <li><b>Золотая зона (Ян):</b> Фаза расширения.</li>
        <li><b>Голубая зона (Инь):</b> Фаза сжатия.</li>
        <li><b>Красная линия (0):</b> Великий Предел. Пересечение этой линии возможно только через S-петлю.</li>
        </ul>
        </div>
        """, unsafe_allow_html=True)
        
        fig, ax = plt.subplots(figsize=(12, 6))
        fig.patch.set_facecolor('#050505')
        ax.set_facecolor('#050505')
        x_vals = range(len(points))
        
        ax.fill_between(x_vals, signed_dists, 0, where=[d>=0 for d in signed_dists], facecolor='#FFD700', alpha=0.5, label='Ян (+)')
        ax.fill_between(x_vals, signed_dists, 0, where=[d<0 for d in signed_dists], facecolor='#00FFFF', alpha=0.5, label='Инь (-)')
        
        ax.plot(x_vals, signed_dists, color='white', linewidth=1.5)
        ax.axhline(0, color='red', linewidth=2, label='Великий Предел')
        
        ax.tick_params(colors='white', labelsize=12)
        ax.legend(facecolor='#222', labelcolor='white', fontsize=12)
        st.pyplot(fig)

    # --- ТАБ 3: ФРАКТАЛ ---
    with t3:
        st.markdown("""
        <div class="info-box">
        <span class="info-title">3. Фрактальность и Золотое Сечение</span>
        Мы проверяем, насколько ритм S-петель соответствует математике живой природы (Числам Фибоначчи и числу Φ ≈ 1.618).
        <ul>
        <li>Если ритм стремится к <b>1.618</b> — структура живая, развивающаяся (Сфираль).</li>
        <li>Если ритм стремится к <b>1.0</b> — структура кристаллическая, цикличная.</li>
        </ul>
        </div>
        """, unsafe_allow_html=True)
        
        if len(intervals) < 2:
            st.warning("⚠️ Недостаточно S-петель для анализа фрактальности.")
        else:
            c1, c2 = st.columns(2)
            with c1:
                st.write("**Длины тактов (Интервалы):**", intervals)
                ratios = [intervals[i+1]/intervals[i] for i in range(len(intervals)-1) if intervals[i]!=0]
                
            with c2:
                # ГРАФИК ОТКЛОНЕНИЯ ОТ ФИ
                fig2, ax2 = plt.subplots(figsize=(8, 5))
                fig2.patch.set_facecolor('#050505')
                ax2.set_facecolor('#050505')
                
                ax2.plot(ratios, marker='o', color='#00FFFF', linewidth=2, label='Ритм Гена')
                ax2.axhline(PHI, color='#FFD700', linestyle='-', linewidth=3, label='Золотое Сечение (1.618)')
                ax2.axhline(1.0, color='#666', linestyle='--', label='Кристалл (1.0)')
                
                ax2.tick_params(colors='white', labelsize=12)
                ax2.legend(facecolor='#222', labelcolor='white')
                ax2.set_title("Поиск Золотой Пропорции", color='white', fontsize=16)
                
                st.pyplot(fig2)
            
            st.download_button("💾 СКАЧАТЬ ДАННЫЕ АНАЛИЗА", str(intervals), "analysis.txt")

elif run_btn:
    st.error("Ошибка: Введите данные для анализа.")