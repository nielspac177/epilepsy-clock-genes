import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, nibabel as nib
from nilearn.image import resample_to_img
from scipy.stats import spearmanr
from scipy.spatial.distance import cdist
from sklearn.decomposition import PCA
import abagen
from brainsmash.mapgen.base import Base
SEED=20260724

# my DK atlas + LGS + AHBA
aimg=nib.load(str(abagen.fetch_desikan_killiany(native=False)["image"]))
info=pd.read_csv("results/spatial/ahba_dk_info.csv").set_index("id")
expr=pd.read_csv("results/spatial/ahba_dk_expression.csv",index_col=0)
cent=pd.read_csv("results/spatial/ahba_dk_centroids.csv").set_index("id")
m=resample_to_img(nib.load("EEGfMRI_GPFA_NOGAONLY_NOTHRESH_TSCORES_MNI2009b.nii.gz"),aimg,interpolation="continuous")
md=np.asarray(m.dataobj); ad=np.asarray(aimg.dataobj)
lgs=pd.Series({int(l):float(np.nanmean(md[ad==l])) for l in np.unique(ad) if l!=0})
# key each DK id by hemi_label to match ENIGMA
info["key"]=info["hemisphere"]+"_"+info["label"]
key2id={k:i for i,k in info["key"].items()}
CLOCK=["ARNTL","ARNTL2","CLOCK","NPAS2","PER1","PER2","PER3","CRY1","CRY2","NR1D1","NR1D2","RORA","RORB","RORC","CSNK1D","FBXL3","DBP","NFIL3","TIMELESS","BHLHE40","BHLHE41","CIART"]
clk=[g for g in CLOCK if g in expr.columns]
Z=(expr-expr.mean())/expr.std()
pc1=pd.Series(PCA(3).fit_transform(Z.values)[:,0],index=expr.index)
clockmap=Z[clk].mean(1)

def load_enigma(f):
    d=pd.read_csv(f"data/raw/enigma/{f}.csv")
    d["id"]=d["Structure"].map(key2id)
    return d.dropna(subset=["id"]).set_index("id")["d_icv"]

for name in ["gge","allepi","tlemtsl"]:
    e=load_enigma(f"{name}_case-controls_CortThick")
    reg=[i for i in e.index if i in expr.index and i in lgs.index and pd.notna(lgs[i])]
    ev=e.loc[reg].to_numpy()
    D=cdist(cent.loc[reg,["x","y","z"]].to_numpy(),cent.loc[reg,["x","y","z"]].to_numpy())
    np.random.seed(SEED)
    surr=Base(x=ev,D=D,resample=True)(n=1000)
    def sp(vec):
        r=spearmanr(vec,ev).statistic; rn=np.array([spearmanr(vec,s).statistic for s in surr])
        return r,(np.sum(np.abs(rn)>=abs(r))+1)/1001
    r_lgs,p_lgs=sp(lgs.loc[reg].to_numpy())
    r_clk,p_clk=sp(clockmap.loc[reg].to_numpy())
    r_pc1,p_pc1=sp(pc1.loc[reg].to_numpy())
    r_per1,p_per1=sp(expr.loc[reg,"PER1"].to_numpy())
    print(f"[{name}] n={len(reg)}  vs LGS: r={r_lgs:+.3f} p={p_lgs:.3f} | vs clock: r={r_clk:+.3f} p={p_clk:.3f} | vs PC1: r={r_pc1:+.3f} p={p_pc1:.3f} | vs PER1: r={r_per1:+.3f} p={p_per1:.3f}")
