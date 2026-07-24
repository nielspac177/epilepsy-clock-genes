import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, nibabel as nib
from nilearn.image import resample_to_img
from scipy.stats import spearmanr
from scipy.spatial.distance import cdist
from sklearn.decomposition import PCA
import abagen
from brainsmash.mapgen.base import Base

expr = pd.read_csv("results/spatial/ahba_dk_expression.csv", index_col=0)
cent = pd.read_csv("results/spatial/ahba_dk_centroids.csv").set_index("id")
aimg = nib.load(str(abagen.fetch_desikan_killiany(native=False)["image"]))
m = resample_to_img(nib.load("EEGfMRI_GPFA_NOGAONLY_NOTHRESH_TSCORES_MNI2009b.nii.gz"), aimg, interpolation="continuous")
md=np.asarray(m.dataobj); ad=np.asarray(aimg.dataobj)
lgs=pd.Series({int(l):float(np.nanmean(md[ad==l])) for l in np.unique(ad) if l!=0}).reindex(expr.index)
ok=expr.notna().all(1)&lgs.notna(); reg=expr.index[ok]
E=expr.loc[reg]; lgs_v=lgs.loc[reg].to_numpy(); Z=(E-E.mean())/E.std()
pca=PCA(n_components=3).fit(Z.values); pc1=pca.transform(Z.values)[:,0]
xyz=cent.loc[reg,["x","y","z"]].to_numpy(); D=cdist(xyz,xyz)
surr=Base(x=lgs_v,D=D,resample=True)(n=1000)
def sp(vec):
    r=spearmanr(vec,lgs_v).statistic; rn=np.array([spearmanr(vec,s).statistic for s in surr])
    return r,(np.sum(np.abs(rn)>=abs(r))+1)/1001
CLOCK=["ARNTL","ARNTL2","CLOCK","NPAS2","PER1","PER2","PER3","CRY1","CRY2","NR1D1","NR1D2","RORA","RORB","RORC","CSNK1D","FBXL3","DBP","NFIL3","TIMELESS","BHLHE40","BHLHE41","CIART"]
clk=[g for g in CLOCK if g in E.columns]; clockmap=Z[clk].mean(1).to_numpy()
r_pc1,p_pc1=sp(pc1); r_cm,p_cm=sp(clockmap); r_per1,p_per1=sp(E["PER1"].to_numpy())
print(f"PC1_var={pca.explained_variance_ratio_[0]*100:.0f}%")
print(f"LGS vs AHBA-PC1(transcriptional gradient): r={r_pc1:.3f} p_spin={p_pc1:.3f}")
print(f"clock-set expr vs PC1: r={spearmanr(clockmap,pc1).statistic:.3f}")
print(f"clock-set expr vs LGS: r={r_cm:.3f} p_spin={p_cm:.3f}")
print(f"PER1 vs PC1: r={spearmanr(E['PER1'],pc1).statistic:.3f} | PER1 vs LGS: r={r_per1:.3f} p={p_per1:.3f}")
print("DONE")
