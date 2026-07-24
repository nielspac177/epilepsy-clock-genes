import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, nibabel as nib
from nilearn.image import resample_to_img
from scipy.stats import spearmanr
import abagen
expr = pd.read_csv("results/spatial/ahba_dk_expression.csv", index_col=0)
info = pd.read_csv("results/spatial/ahba_dk_info.csv").set_index("id")
aimg = nib.load(str(abagen.fetch_desikan_killiany(native=False)["image"]))
m = resample_to_img(nib.load("EEGfMRI_GPFA_NOGAONLY_NOTHRESH_TSCORES_MNI2009b.nii.gz"), aimg, interpolation="continuous")
md=np.asarray(m.dataobj); ad=np.asarray(aimg.dataobj)
lgs=pd.Series({int(l):float(np.nanmean(md[ad==l])) for l in np.unique(ad) if l!=0}).reindex(expr.index)
CLOCK=["ARNTL","ARNTL2","CLOCK","NPAS2","PER1","PER2","PER3","CRY1","CRY2","NR1D1","NR1D2","RORA","RORB","RORC","CSNK1D","FBXL3","DBP","NFIL3","TIMELESS","BHLHE40","BHLHE41","CIART"]
clk=[g for g in CLOCK if g in expr.columns]
Z=(expr-expr.mean())/expr.std(); clockmap=Z[clk].mean(1)
df=pd.DataFrame({"label":info["label"],"struct":info["structure"],"lgs":lgs,"clock":clockmap,"PER1":expr["PER1"]}).dropna()
print("=== Top 12 LGS network regions (where the network peaks) ===")
top=df.sort_values("lgs",ascending=False).head(12)
for i,r in top.iterrows(): print(f"  {r['label']:28s} {r['struct']:20s} LGS_t={r['lgs']:+.2f}")
print(f"\nsubcortical fraction of top-12 LGS: {(top['struct']!='cortex').sum()}/12")
print("\n=== clock-set / PER1 vs LGS, split by tissue ===")
for name,sub in [("whole-brain",df),("cortex-only",df[df['struct']=='cortex']),("subcortex/brainstem",df[df['struct']!='cortex'])]:
    rc=spearmanr(sub['clock'],sub['lgs']).statistic; rp=spearmanr(sub['PER1'],sub['lgs']).statistic
    print(f"  {name:20s} (n={len(sub):2d}): clock r={rc:+.3f}  PER1 r={rp:+.3f}")
print("\n=== thalamus specifically ===")
th=df[df['label'].str.contains('thalamus',case=False)]
for i,r in th.iterrows(): print(f"  {r['label']:20s} LGS_t={r['lgs']:+.2f}  clock_z={r['clock']:+.2f}  PER1={r['PER1']:.2f}")
print(f"  clock expr rank of thalamus among 83 regions: {[f'{r.label}:{(df.clock<r.clock).sum()+1}/83' for _,r in th.iterrows()]}")
