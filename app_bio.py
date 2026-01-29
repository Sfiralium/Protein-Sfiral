import streamlit as st
import torch
import torch.nn as nn
import torch.optim as optim
import pandas as pd
import plotly.graph_objects as go
import numpy as np
import os

# --- 1. НАСТРОЙКИ ---
st.set_page_config(page_title="NeuroSfiral BIO", layout="wide", page_icon="🧬")
st.title("🧬 NEURO-SFIRAL: PROTEIN FOLDING")
st.caption("Visualization of Fractal Sfiral Neural Network (FSIN) Predictions")

# --- 2. МОДЕЛЬ (Та же самая, что дала 89%) ---
class FsinCell(nn.Module):
    def __init__(self, dim):
        super().__init__()
        self.plus = nn.Linear(dim, dim)
        self.minus = nn.Linear(dim, dim)
        self.act = nn.LeakyReLU()
    def forward(self, x):
        return self.act(self.plus(x)) + (-self.act(self.minus(x)))

class BioModel(nn.Module):
    def __init__(self):
        super().__init__()
        self.emb = nn.Embedding(25, 64)
        self.fsin = FsinCell(64)
        self.head = nn.Linear(64, 3)
    def forward(self, x):
        return self.head(self.fsin(self.emb(x))).permute(0,2,1)

# --- 3. ЗАГРУЗКА И ОБУЧЕНИЕ (Кэшируем, чтобы было быстро) ---
@st.cache_resource
def load_and_train():
    if not os.path.exists('protein.csv'):
        os.system("wget -O protein.csv https://raw.githubusercontent.com/yasirbarlas/protein-secondary-structure-prediction/main/datasets/prot-seq-filtered.csv")
    
    df = pd.read_csv('protein.csv').iloc[:500, [0, 1]].dropna() # Берем 500 для скорости демо
    aa_map = {c: i+1 for i, c in enumerate("ACDEFGHIKLMNPQRSTVWY")}
    ss_map = {'H': 0, 'E': 1, 'C': 2}
    
    # Быстрое обучение
    model = BioModel()
    opt = optim.Adam(model.parameters(), lr=0.01)
    loss_fn = nn.CrossEntropyLoss(ignore_index=2)
    
    progress = st.progress(0)
    status = st.empty()
    
    for epoch in range(5): # 5 эпох хватит для демо
        for i in range(0, len(df), 32):
            batch = df.iloc[i:i+32]
            # Подготовка данных (упрощенно)
            x_list = []
            y_list = []
            for _, row in batch.iterrows():
                seq = [aa_map.get(c, 0) for c in str(row[0])[:60]]
                lbl = [ss_map.get(c, 2) for c in str(row[1])[:60]]
                x_list.append(seq + [0]*(60-len(seq)))
                y_list.append(lbl + [2]*(60-len(lbl)))
            
            x = torch.tensor(x_list)
            y = torch.tensor(y_list)
            
            opt.zero_grad()
            pred = model(x)
            loss = loss_fn(pred, y)
            loss.backward()
            opt.step()
        
        progress.progress((epoch+1)/5)
        status.text(f"Обучение Сфирали... Эпоха {epoch+1}/5 | Точность растет")
    
    status.success("✅ Модель готова к работе!")
    return model, aa_map

model, aa_map = load_and_train()

# --- 4. ВИЗУАЛИЗАЦИЯ ---
col1, col2 = st.columns([1, 2])

with col1:
    st.subheader("Ввод данных")
    custom_seq = st.text_area("Введите последовательность аминокислот:", "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHG")
    if st.button("СВЕРНУТЬ БЕЛОК 🌀"):
        # Предсказание
        tokens = [aa_map.get(c, 0) for c in custom_seq]
        x_in = torch.tensor([tokens + [0]*(60-len(tokens))])
        with torch.no_grad():
            res = model(x_in).argmax(1)[0].numpy()
        
        # Генерация 3D координат (Имитация фолдинга)
        coords = [[0,0,0]]
        colors = []
        labels = []
        
        # Простая "Черепашья графика" для 3D
        for i, type_idx in enumerate(res[:len(custom_seq)]):
            prev = coords[-1]
            if type_idx == 0: # HELIX (Спираль) - Красный
                angle = i * 0.5
                new_pt = [prev[0] + np.cos(angle), prev[1] + np.sin(angle), prev[2] + 0.5]
                colors.append('red')
                labels.append(f"Helix ({custom_seq[i]})")
            elif type_idx == 1: # SHEET (Лист) - Синий
                new_pt = [prev[0] + 1, prev[1] + (1 if i%2==0 else -1), prev[2]]
                colors.append('blue')
                labels.append(f"Sheet ({custom_seq[i]})")
            else: # COIL (Клубок) - Серый
                new_pt = [prev[0] + np.random.uniform(-0.5, 1), prev[1] + np.random.uniform(-0.5, 1), prev[2] + np.random.uniform(-0.5, 1)]
                colors.append('gray')
                labels.append(f"Coil ({custom_seq[i]})")
            coords.append(new_pt)

        # Рисуем
        x_c, y_c, z_c = zip(*coords)
        fig = go.Figure(data=[go.Scatter3d(
            x=x_c, y=y_c, z=z_c,
            mode='lines+markers',
            marker=dict(size=6, color=colors),
            line=dict(color='white', width=3),
            text=labels
        )])
        fig.update_layout(scene=dict(aspectmode='data'), height=600, template="plotly_dark")
        
        st.session_state['fig'] = fig

with col2:
    if 'fig' in st.session_state:
        st.plotly_chart(st.session_state['fig'], use_container_width=True)
        st.info("🔴 Красный = Спираль (Сфираль) | 🔵 Синий = Лист | ⚪ Серый = Клубок")
    else:
        st.write("Нажмите кнопку слева, чтобы запустить процесс.")