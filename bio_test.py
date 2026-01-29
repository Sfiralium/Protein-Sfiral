import os
import torch
import torch.nn as nn
import torch.optim as optim
import pandas as pd
from torch.utils.data import Dataset, DataLoader

# --- 1. АВТО-ЗАГРУЗКА ДАННЫХ ---
print("📥 Скачиваю базу данных белков...")
if not os.path.exists('protein.csv'):
    os.system("wget -O protein.csv https://raw.githubusercontent.com/yasirbarlas/protein-secondary-structure-prediction/main/datasets/prot-seq-filtered.csv")

# --- 2. ПОДГОТОВКА ДАННЫХ ---
class ProteinData(Dataset):
    def __init__(self):
        try:
            df = pd.read_csv('protein.csv')
            # Берем колонки: seq (буквы) и sst3 (форма)
            self.data = df.iloc[:, [0, 1]].dropna().iloc[:3000] # 3000 белков
            self.aa_map = {c: i+1 for i, c in enumerate("ACDEFGHIKLMNPQRSTVWY")}
            self.ss_map = {'H': 0, 'E': 1, 'C': 2} # Спираль, Лист, Клубок
            print(f"✅ Данные готовы: {len(self.data)} образцов.")
        except:
            print("⚠️ Ошибка данных. Создаю тестовый набор.")
            self.data = pd.DataFrame({'seq':['A']*100, 'sst3':['H']*100})
            self.aa_map = {'A':1}
            self.ss_map = {'H':0}

    def __len__(self): return len(self.data)

    def __getitem__(self, i):
        # Превращаем буквы в цифры
        seq = [self.aa_map.get(c, 0) for c in str(self.data.iloc[i, 0])[:60]]
        lbl = [self.ss_map.get(c, 2) for c in str(self.data.iloc[i, 1])[:60]]
        # Выравниваем длину до 60 (padding)
        seq += [0]*(60-len(seq))
        lbl += [2]*(60-len(lbl))
        return torch.tensor(seq), torch.tensor(lbl)

# --- 3. НЕЙРОСФИРАЛЬ (FSIN) ---
class FsinCell(nn.Module):
    def __init__(self, dim):
        super().__init__()
        self.plus = nn.Linear(dim, dim)
        self.minus = nn.Linear(dim, dim)
        self.act = nn.LeakyReLU()
    def forward(self, x):
        # Зеркальная антисимметрия: (V+) + (-V-)
        return self.act(self.plus(x)) + (-self.act(self.minus(x)))

class BioModel(nn.Module):
    def __init__(self):
        super().__init__()
        self.emb = nn.Embedding(25, 64)
        self.fsin = FsinCell(64)          # Сфиральный слой
        self.head = nn.Linear(64, 3)      # Выход (3 класса формы)
    def forward(self, x):
        return self.head(self.fsin(self.emb(x))).permute(0,2,1)

# --- 4. ЗАПУСК ---
if __name__ == "__main__":
    print("🧬 ЗАПУСК НЕЙРОСФИРАЛИ...")
    loader = DataLoader(ProteinData(), batch_size=32, shuffle=True)
    model = BioModel()
    opt = optim.Adam(model.parameters(), lr=0.005)
    loss_fn = nn.CrossEntropyLoss(ignore_index=2)

    print("\n🚀 ОБУЧЕНИЕ (Наблюдайте за ростом точности):")
    for epoch in range(1, 16): # 15 эпох
        correct = 0
        total = 0
        for x, y in loader:
            opt.zero_grad()
            pred = model(x)
            loss = loss_fn(pred, y)
            loss.backward()
            opt.step()
            
            # Считаем совпадения
            choice = pred.argmax(1)
            mask = (y != 2) # Не считаем пустые места
            correct += (choice[mask] == y[mask]).sum().item()
            total += mask.sum().item()

        acc = correct/total * 100 if total > 0 else 0
        bar = "▓" * int(acc/5) + "░" * (20 - int(acc/5))
        print(f"Эпоха {epoch:02d} | {bar} | Точность: {acc:.1f}%")