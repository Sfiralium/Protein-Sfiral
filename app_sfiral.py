import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import py3Dmol
from stmol import showmol
import re

# --- КОНФИГУРАЦИЯ ЛАБОРАТОРИИ ---
st.set_page_config(page_title="Sfiral Laboratory: Universal", layout="wide", page_icon="☯️")

st.markdown("""
<style>
    .stApp {background-color: #0e1117; color: #e0e0e0;}
    /* Золотые кнопки */
    div.stButton > button:first-child {
        background-color: #D4AF37; 
        color: black; 
        border: none; 
        height: 50px; 
        font-weight: bold; 
        font-size: 16px;
    }
    h1, h2, h3 {color: #D4AF37 !important;}
    .stTabs [aria-selected="true"] {background-color: #D4AF37; color: #000;}
    .metric-box {border: 1px solid #333; padding: 10px; border-radius: 5px; background: #1a1a1a;}
</style>
""", unsafe_allow_html=True)

# --- КОНСТАНТЫ СФИРАЛИ ---
BASE_ANGLES = {'z': 97.8, 'y': -55.97, 'x': -29.7}
MOVE_VEC = np.array([-0.8, 2.15, -1.37])
PHI = 1.61803398875

# Золотые Углы (Инсулиновый Стандарт)
STATES = {
    'alpha': {'x': -44, 'y': 147, 'z': -16}, # Ян (Накопление)
    'pi':    {'x': 10,  'y': 10,  'z': 0},   # Ян (Вариация)
    'beta':  {'x': 114, 'y': -21, 'z': -17}, # Дао (S-Петля / Инвертор)
    '310':   {'x': 45,  'y': 0,   'z': 0}    # Инь (Сжатие)
}

# --- МАТЕМАТИЧЕСКОЕ ЯДРО ---
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

# --- ГЕНЕРАТОР МОДЕЛИ ---
def build_sfiral_model(seq):
    codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
    pos = np.array([0.0, 0.0, 0.0])
    rot = np.identity(3)
    points = []
    
    # Базовая ротация
    R_base = rot_z(BASE_ANGLES['z']) @ rot_y(BASE_ANGLES['y']) @ rot_x(BASE_ANGLES['x'])
    
    current_polarity = 1 # 1 = Ян, -1 = Инь
    atom_id = 1
    pdb_str = ""
    
    for i, codon in enumerate(codons):
        t = get_type(codon)
        s = STATES[t]
        
        # Вращение и шаг
        rot = rot @ (R_base @ (rot_y(s['y']) @ rot_x(s['x']) @ rot_z(s['z'])))
        pos = pos + (rot @ MOVE_VEC)
        
        # Свойства точки
        phase = "coil"
        radius = 0.6
        pdb_bfactor = 1.0
        
        # Логика цвета зависит от Полярности
        color = "#FFD700" if current_polarity > 0 else "#00FFFF" # Золото vs Циан
        
        if t == 'beta':
            phase = "s-loop"
            color = "#FFFFFF" # Белый (Точка Перехода)
            radius = 1.2
            pdb_bfactor = 99.0
            # ГЛАВНОЕ: S-петля переключает полярность Вселенной белка
            current_polarity *= -1 
            
        points.append({
            'pos': pos, 'color': color, 'r': radius, 
            'phase': phase, 'idx': i, 'polarity': current_polarity
        })
        
        # Запись в PDB (сохраняем фазу в B-factor)
        pdb_str += f"ATOM  {atom_id:5d}  CA  UNK A{i+1:4d}    {pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}  1.00{pdb_bfactor:6.2f}           C\n"
        atom_id += 1
        
    return points, pdb_str

# --- ИНТЕРФЕЙС ---
st.title("☯️ SFIRAL LABORATORY: UNIVERSAL")
st.caption("Time-Genetics Analysis System | Version 2.0 Final")

# САЙДБАР
with st.sidebar:
    st.header("Ввод Генетического Кода")
    
    # Пресеты
    preset = st.selectbox("Быстрая загрузка:", ["Инсулин (Ритм 7)", "Гемоглобин (Ритм 90)", "Коллаген (Ритм 0)", "Свой код"])
    
    insulin_seq = "ATGGCCCTGTGGATGCGCCTCCTGCCCCTGCTGGCGCTGCTGGCCCTCTGGGGACCTGACCCAGCCGCAGCCTTTGTGAACCAACACCTGTGCGGCTCACACCTGGTGGAAGCTCTCTACCTAGTGTGCGGGGAACGAGGCTTCTTCTACACACCCAAGACCCGCCGGGAGGCAGAGGACCTGCAGGTGGGGCAGGTGGAGCTGGGCGGGGGCCCTGGTGCAGGCAGCCTGCAGCCCTTGGCCCTGGAGGGGTCCCTGCAGAAGCGTGGCATTGTGGAACAATGCTGTACCAGCATCTGCTCCCTCTACCAGCTGGAGAACTACTGCAACTAG"
    hbb_seq = "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCACTGTGACAAGCTGCACGTGGATCCTGAGAACTTCAGGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGCAAAGAATTCACCCCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAATGCCCTGGCCCACAAGTATCACTAA"
    col_seq = "GGTCCTCGTGGTCTCCCTGGCCCCCCTGGT" * 5
    
    default_txt = insulin_seq
    if preset == "Гемоглобин (Ритм 90)": default_txt = hbb_seq
    elif preset == "Коллаген (Ритм 0)": default_txt = col_seq
    elif preset == "Свой код": default_txt = ""
    
    txt = st.text_area("CDS Sequence:", value=default_txt, height=150)
    clean_seq = re.sub(r'[^a-zA-Z]', '', txt).upper().replace("U", "T")
    
    st.info(f"Кодонов: {len(clean_seq)//3}")
    run = st.button("🚀 ЗАПУСТИТЬ АНАЛИЗ")

# --- ГЛАВНАЯ СЦЕНА ---
if run and len(clean_seq) > 0:
    points, pdb_data = build_sfiral_model(clean_seq)
    
    # Подготовка данных
    center = np.mean([p['pos'] for p in points], axis=0)
    dists = [np.linalg.norm(p['pos'] - center) for p in points]
    signed_dists = [d * p['polarity'] for d, p in zip(dists, points)]
    s_loops = [p['idx'] for p in points if p['phase'] == 's-loop']
    
    # Вектор ошибки
    end_pos = points[-1]['pos']
    err_vec = np.linalg.norm(end_pos)
    
    # ТАБЫ
    t1, t2, t3 = st.tabs(["🧬 3D ГЕОМЕТРИЯ", "📉 АНТИСИММЕТРИЯ (ДЫХАНИЕ)", "🌀 ФРАКТАЛ И ФИБОНАЧЧИ"])
    
    # --- ТАБ 1: 3D ---
    with t1:
        c1, c2 = st.columns([3, 1])
        with c1:
            view = py3Dmol.view(width=800, height=600)
            prev = None
            for p in points:
                view.addSphere({'center':{'x':p['pos'][0],'y':p['pos'][1],'z':p['pos'][2]}, 'radius':p['r'], 'color':p['color']})
                if prev:
                    # Связь белая, если рядом S-петля
                    lnk = "#FFF" if (p['phase']=='s-loop' or prev['phase']=='s-loop') else p['color']
                    view.addCylinder({'start':{'x':prev['pos'][0],'y':prev['pos'][1],'z':prev['pos'][2]}, 
                                      'end':{'x':p['pos'][0],'y':p['pos'][1],'z':p['pos'][2]}, 'radius':0.2, 'color':lnk})
                prev = p
            
            # Красная стрелка (Вектор Ошибки)
            view.addArrow({'start':{'x':end_pos[0],'y':end_pos[1],'z':end_pos[2]}, 'end':{'x':0,'y':0,'z':0}, 'color':'red', 'radius':0.3})
            view.setBackgroundColor('#0e1117')
            view.zoomTo()
            showmol(view, height=600, width=800)
            
        with c2:
            st.markdown("### Паспорт Сфирали")
            st.metric("Вектор Смещения", f"{err_vec:.1f} Å")
            st.metric("Кол-во S-петель", len(s_loops))
            
            st.download_button("💾 СКАЧАТЬ PDB", pdb_data, "sfiral_model.pdb", mime="chemical/x-pdb")
            st.info("В файле PDB в колонке B-factor записана фаза (1.0 = Виток, 99.0 = S-петля).")

    # --- ТАБ 2: АНТИСИММЕТРИЯ ---
    with t2:
        st.markdown("### Фазовое Дыхание (+/-)")
        fig, ax = plt.subplots(figsize=(12, 5))
        fig.patch.set_facecolor('#0e1117')
        ax.set_facecolor('#0e1117')
        
        x_vals = range(len(points))
        # Заливка Ян (Золото) и Инь (Циан)
        ax.fill_between(x_vals, signed_dists, 0, where=[d>=0 for d in signed_dists], facecolor='#FFD700', alpha=0.4, label='Ян (Расширение)')
        ax.fill_between(x_vals, signed_dists, 0, where=[d<0 for d in signed_dists], facecolor='#00FFFF', alpha=0.4, label='Инь (Сжатие)')
        
        ax.plot(x_vals, signed_dists, color='white', linewidth=1)
        ax.axhline(0, color='red', linestyle='-', linewidth=0.8, label='Великий Предел (0)')
        
        # Точки переходов
        for l in s_loops:
            ax.axvline(x=l, color='white', linestyle='--', alpha=0.3)
            
        ax.tick_params(colors='white')
        ax.legend(facecolor='#222', labelcolor='white')
        st.pyplot(fig)
        st.caption("График показывает, как структура пересекает нулевую отметку в моменты S-петель, меняя хиральность/полярность.")

    # --- ТАБ 3: ФРАКТАЛ ---
    with t3:
        st.markdown("### Поиск Золотого Сечения (Φ = 1.618)")
        
        intervals = []
        if len(s_loops) > 0:
            prev = 0
            for l in s_loops:
                intervals.append(l - prev)
                prev = l
        
        if len(intervals) < 2:
            st.warning("Слишком мало переходов для анализа ритма.")
        else:
            c1, c2 = st.columns(2)
            
            with c1:
                st.write("**Длины фаз (такты):**")
                st.write(intervals)
                
                # Фибоначчи тест
                fibs = [1,2,3,5,8,13,21,34,55,89,144]
                matches = sum(1 for v in intervals if any(abs(v-f)<=1 for f in fibs))
                st.metric("Совпадение с Фибоначчи", f"{(matches/len(intervals))*100:.0f}%")
                
            with c2:
                # Коэффициенты роста
                ratios = [intervals[i+1]/intervals[i] for i in range(len(intervals)-1) if intervals[i]!=0]
                avg_ratio = np.mean(ratios) if ratios else 0
                st.metric("Средний Рост Фазы", f"{avg_ratio:.3f}")
                st.metric("Цель (Φ)", "1.618")
                
            # График отклонений
            fig2, ax2 = plt.subplots(figsize=(10, 3))
            fig2.patch.set_facecolor('#0e1117')
            ax2.set_facecolor('#0e1117')
            ax2.plot(ratios, marker='o', color='#00FFFF', label='Ритм Гена')
            ax2.axhline(PHI, color='#FFD700', linestyle='-', label='Золотое Сечение')
            ax2.tick_params(colors='white')
            ax2.legend()
            st.pyplot(fig2)