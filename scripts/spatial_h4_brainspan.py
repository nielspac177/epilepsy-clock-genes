import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, nibabel as nib
from nilearn.image import resample_to_img
from scipy.stats import spearmanr
import abagen

# LGS-t per DK region (reuse parcellation)
aimg=nib.load(str(abagen.fetch_desikan_killiany(native=False)["image"]))
info=pd.read_csv("results/spatial/ahba_dk_info.csv").set_index("id")
m=resample_to_img(nib.load("EEGfMRI_GPFA_NOGAONLY_NOTHRESH_TSCORES_MNI2009b.nii.gz"),aimg,interpolation="continuous")
md=np.asarray(m.dataobj); ad=np.asarray(aimg.dataobj)
lgs={int(l):float(np.nanmean(md[ad==l])) for l in np.unique(ad) if l!=0}
lab2t={}
for i,r in info.iterrows():
    lab2t.setdefault(r["label"],[]).append(lgs.get(i,np.nan))
lab2t={k:np.nanmean(v) for k,v in lab2t.items()}
# BrainSpan cortical region -> DK labels
MAP={"DFC":["rostralmiddlefrontal","caudalmiddlefrontal","superiorfrontal"],
 "VFC":["parsopercularis","parstriangularis","parsorbitalis"],
 "MFC":["medialorbitofrontal","rostralanteriorcingulate","caudalanteriorcingulate"],
 "OFC":["lateralorbitofrontal","medialorbitofrontal"],"M1C":["precentral"],"S1C":["postcentral"],
 "IPC":["inferiorparietal","supramarginal"],"A1C":["transversetemporal"],"STC":["superiortemporal"],
 "ITC":["inferiortemporal"],"V1C":["pericalcarine","lateraloccipital","cuneus"]}
bs_lgs={r:np.nanmean([lab2t[l] for l in ls if l in lab2t]) for r,ls in MAP.items()}

# BrainSpan expression: clock gene rows only
rows=pd.read_csv("data/raw/brainspan/rows_metadata.csv")
cols=pd.read_csv("data/raw/brainspan/columns_metadata.csv")
CLOCK=["ARNTL","ARNTL2","CLOCK","NPAS2","PER1","PER2","PER3","CRY1","CRY2","NR1D1","NR1D2","RORA","RORB","RORC","CSNK1D","CSNK1E","FBXL3","DBP","NFIL3","TIMELESS","BHLHE40","BHLHE41","CIART"]
clk_rows=rows[rows["gene_symbol"].isin(CLOCK)]
idx=clk_rows["row_num"].values-1
expr=pd.read_csv("data/raw/brainspan/expression_matrix.csv",header=None,index_col=0)
E=np.log2(expr.values[idx]+1)  # clock genes x samples
Egenes=clk_rows["gene_symbol"].values
def agebin(a):
    if "pcw" in a: return "1_prenatal"
    v=float(a.split()[0]); u=a.split()[1]
    mos=v if u=="mos" else v*12
    if mos<12: return "2_infancy"
    yr=mos/12
    if yr<=5: return "3_earlychild(LGS-onset)"
    if yr<19: return "4_child_adol"
    return "5_adult"
cols["bin"]=cols["age"].apply(agebin)
print("=== H4: clock-set expression vs LGS map across 11 cortical regions, by developmental stage ===")
for b in sorted(cols["bin"].unique()):
    sel=cols[(cols["bin"]==b)&(cols["structure_acronym"].isin(MAP))]
    if len(sel)<8: 
        print(f"  {b:26s} (too few samples)"); continue
    # mean clock expr per region (z across clock genes then mean)
    perreg={}
    for reg in MAP:
        s=sel[sel["structure_acronym"]==reg]
        if len(s)==0: continue
        sub=E[:,s.index.values]  # subtract 1-index? columns_metadata row order = sample col order
        perreg[reg]=np.nanmean(sub)
    regs=[r for r in perreg if np.isfinite(bs_lgs.get(r,np.nan))]
    x=[perreg[r] for r in regs]; y=[bs_lgs[r] for r in regs]
    r=spearmanr(x,y).statistic
    print(f"  {b:26s} n_reg={len(regs):2d} n_samp={len(sel):3d}  clock-vs-LGS r={r:+.3f}")
