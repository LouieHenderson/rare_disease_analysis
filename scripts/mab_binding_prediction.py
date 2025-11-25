from torch.utils.data import Dataset
from torch.utils.data import DataLoader
import torch.optim as optim
import torch.nn as nn
import matplotlib.pyplot as plt

def loadEmbd(dataframe, outfold, sumdf):
    embdic= {}
    for i, row in dataframe.iterrows():
        #print(f'{row["pdb_id"]} {sumdf.loc[sumdf["pdb"] == row["pdb_id"]]["affinity"].values[0]}')
        #print(f'{outfold}{row["pdb_id"]}_vh.pt, {row["pdb_id"]}_vl.pt, {row["pdb_id"]}_ag.pt,')
        gc.collect()
        torch.cuda.empty_cache()
        embd_vh = torch.load(f'{outfold}{row["pdb_id"]}_vh.pt')
        embd_vl = torch.load(f'{outfold}{row["pdb_id"]}_vl.pt')
        embd_ag = torch.load(f'{outfold}{row["pdb_id"]}_ag.pt')
        if "affin" in sumdf:
            #display(sumdf)
            affin = sumdf.loc[sumdf["pdb_id"] == row["pdb_id"]]["affin"].values[0]
        else:
            affin = sumdf.loc[sumdf["pdb"] == row["pdb_id"]]["affinity"].values[0]
        #break
        if np.isnan(affin) == False:
            #print(affin, np.isnan(affin))
            embdic.setdefault(row["pdb_id"], [embd_vh, embd_vl, embd_ag, affin])

    return embdic

def preprocess_dict(raw_dict):
    processed = {}
    for pdb_id, (vh, vl, ag, aff) in raw_dict.items():
        vh = vh.squeeze(0)  # remove batch dim → (L_vh, d)
        vl = vl.squeeze(0)
        ag = ag.squeeze(0)
        processed[pdb_id] = (vh, vl, ag, float(aff))
    return processed

def mean_pool(t):
    return t.mean(dim=0)  # (L, d) → (d,)

def makeRecord(fastaobj, row, thresh, noemp):
    goodrec  = False
    name     = row["light"]+"_"+row["heavy"]+"_"+row["antigen"]
    light, heavy, antig = None, None, None
    affinity = row["delta_g"]
    try:
        light    = str(fastaobj[row["light"]].seq)
        heavy    = str(fastaobj[row["heavy"]].seq)
        antig    = str(fastaobj[row["antigen"]].seq)
        goodrec = True
    except:
        #print("YAYAY")
        light    = row["light"]
        heavy    = row["heavy"]
        antig    = row["antigen"]
      #  print("mistake = ", light, heavy, antig)

    record = {"pdb_id": name, 
              "vh_seq": heavy, 
              "vl_seq": light, 
              "antigen": antig, 
              "affin": affinity}
    
    if noemp == False:
        if len(light) > thresh or len(heavy) > thresh or len(antig) > thresh:
            goodrec = False
    else:
        if len(light) > thresh or len(heavy) > thresh or len(antig) > thresh or "EMPTY" in name:
            goodrec = False

    return record, goodrec

class AntibodyAffinityDataset(Dataset):
    def __init__(self, data_dict, pool_fn=mean_pool):
        self.items = list(data_dict.items())
        self.pool_fn = pool_fn

    def __len__(self):
        return len(self.items)

    def __getitem__(self, idx):
        pdb_id, (vh, vl, ag, y) = self.items[idx]
        vh_vec = self.pool_fn(vh)
        vl_vec = self.pool_fn(vl)
        ag_vec = self.pool_fn(ag)
        y_val = torch.tensor(np.log10(y), dtype=torch.float32)
        return vh_vec, vl_vec, ag_vec, y_val

class AntibodyAffinityDataset2(Dataset):
    def __init__(self, data_dict, pool_fn=mean_pool):
        self.items = list(data_dict.items())
        self.pool_fn = pool_fn

    def __len__(self):
        return len(self.items)

    def __getitem__(self, idx):
        pdb_id, (vh, vl, ag, y) = self.items[idx]
        vh_vec = self.pool_fn(vh)
        vl_vec = self.pool_fn(vl)
        ag_vec = self.pool_fn(ag)
        y_val = torch.tensor(y, dtype=torch.float32)
        return vh_vec, vl_vec, ag_vec, y_val

class AffiniPredMulti(nn.Module):
    def __init__(self, emb_dim, hidden_dim=256, fusion_dim=128):
        super().__init__()
        # Indiv branch analysis
        self.seq = nn.Sequential(
            nn.Linear(emb_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(hidden_dim, fusion_dim),
            nn.ReLU()
        )
        # Fusion & regression head
        self.reg_head = nn.Sequential(
            nn.Linear(fusion_dim * 3, 256),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(256, 1)
        )

    def forward(self, vh, vl, ag):
        vh_latent = self.seq(vh)
        vl_latent = self.seq(vl)
        ag_latent = self.seq(ag)

        fused = torch.cat([vh_latent, vl_latent, ag_latent], dim=1)
        out = self.reg_head(fused)
        return out

def SanityCheck(model, test_loader, label, test_loss, train_loss, miny, maxy):
    model.eval()
    y_true_list = []
    y_pred_list = []
    
    with torch.no_grad():
        for vh, vl, ag, y_true in test_loader:
            y_pred = model(vh, vl, ag)
            y_true_list.append(y_true.numpy())
            y_pred_list.append(y_pred.squeeze(1).numpy())
    
    import numpy as np
    y_true_all = np.concatenate(y_true_list)
    y_pred_all = np.concatenate(y_pred_list)
    
    from sklearn.metrics import mean_squared_error, r2_score
    from scipy.stats import pearsonr, spearmanr
    
    mse = mean_squared_error(y_true_all, y_pred_all)
    rmse = np.sqrt(mse)
    r2 = r2_score(y_true_all, y_pred_all)
    pearson_r, _ = pearsonr(y_true_all, y_pred_all)
    spearman_r, _ = spearmanr(y_true_all, y_pred_all)
    
    print(f"Test MSE: {mse:.4f}")
    print(f"Test RMSE: {rmse:.4f}")
    print(f"Test R²: {r2:.4f}")
    print(f"Test Pearson r: {pearson_r:.4f}")
    print(f"Test Spearman r: {spearman_r:.4f}")


    fig, plots = plt.subplots(nrows=1, ncols=2)
    
    #Affinity prediction plot
    plots[0].scatter(y_true_all, y_pred_all, alpha=0.4)
    plots[0].plot([y_true_all.min(), y_true_all.max()],
             [y_true_all.min(), y_true_all.max()],
             'r--')
    plots[0].set_xlabel(f'True {label} (affinity)')
    plots[0].set_ylabel(f"Predicted {label} (affinity)")
    plots[0].set_title("Affinity Prediction on Test")
    plots[0].grid(False)

    #Gradient descent loss plot
    plots[1].plot(train_loss, label="Train loss", linestyle='--')
    plots[1].plot(test_loss, label="Test loss", linestyle='--')
    plots[1].set_xlabel(f'Epoch')
    plots[1].set_ylabel(f"Averaged loss")
    plots[1].set_title("Training validation")
    plots[1].grid(False)
    plots[1].legend()
    plots[1].set_ylim(miny, maxy)
    fig.show()

def trainModel(model, train_loader, test_loader, epoch_num):
    criterion = nn.MSELoss()
    optimizer = optim.Adam(model.parameters(), lr=1e-4)
    train_all = []
    test_all  = []
    
    for epoch in range(epoch_num):
        model.train()
        total_loss = 0.0
        for vh, vl, ag, y in train_loader:
            optimizer.zero_grad()
            pred = model(vh, vl, ag)
            loss = criterion(pred, y.unsqueeze(1))
            loss.backward()
            optimizer.step()
            total_loss += loss.item() * vh.size(0)
    
        avg_loss = total_loss / len(train_loader.dataset)
    
        model.eval()
        test_loss = 0.0
        with torch.no_grad():
            for vh, vl, ag, y in test_loader:
                pred = model(vh, vl, ag)
                loss = criterion(pred, y.unsqueeze(1))
                test_loss += loss.item() * vh.size(0)
    
        test_loss /= len(test_loader.dataset)
        train_all.append(avg_loss)
        test_all.append(test_loss)
        print(f"Epoch {epoch+1:02d} | Train Loss: {avg_loss:.4f} | Test Loss: {test_loss:.4f}")

    return model, train_all, test_all