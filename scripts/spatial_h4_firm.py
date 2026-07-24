import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, nibabel as nib
from nilearn.image import resample_to_img
from scipy.stats import spearmanr
import abagen
rng=np.random.default_rng(20260724)
aimg=nib.load(str(abagen.fetch_desikan_killiany(native=False)["image"]))
info=pd.read_csv("results/spatial/ahba_dk_info.csv").set_index("id")
m=resample_to_img(nib.load("EEGfMRI_GPFA_NOGAONLY_NOTHRESH_TSCORES_MNI2009b.nii.gz"),aimg,interpolation="continuous")
md=np.asarray(m.dataobj); ad=np.asarray(aimg.dataobj)
lgs={int(l):float(np.nanmean(md[ad==l])) for l in np.unique(ad) if l!=0}
lab2t={}
for i,r in info.iterrows(): lab2t.setdefault(r["label"],[]).append(lgs.get(i,np.nan))
lab2t={k:np.nanmean(v) for k,v in lab2t.items()}
# BrainSpan region -> DK labels (cortical + SUBCORTICAL, incl mediodorsal thalamus)
MAP={"DFC":["rostralmiddlefrontal","caudalmiddlefrontal","superiorfrontal"],"VFC":["parsopercularis","parstriangularis","parsorbitalis"],
 "MFC":["medialorbitofrontal","rostralanteriorcingulate","caudalanteriorcingulate"],"OFC":["lateralorbitofrontal","medialorbitofrontal"],
 "M1C":["precentral"],"S1C":["postcentral"],"IPC":["inferiorparietal","supramarginal"],"A1C":["transversetemporal"],
 "STC":["superiortemporal"],"ITC":["inferiortemporal"],"V1C":["pericalcarine","lateraloccipital","cuneus"],
 "MD":["thalamusproper"],"DTH":["thalamusproper"],"STR":["caudate","putamen","accumbensarea"],"HIP":["hippocampus"],"AMY":["amygdala"]}
bs_lgs={r:np.nanmean([lab2t[l] for l in ls if l in lab2t]) for r,ls in MAP.items()}
rows=pd.read_csv("data/raw/brainspan/rows_metadata.csv"); cols=pd.read_csv("data/raw/brainspan/columns_metadata.csv")
CLOCK=["ARNTL","ARNTL2","CLOCK","NPAS2","PER1","PER2","PER3","CRY1","CRY2","NR1D1","NR1D2","RORA","RORB","RORC","CSNK1D","CSNK1E","FBXL3","DBP","NFIL3","TIMELESS","BHLHE40","BHLHE41","CIART"]
clk=rows[rows["gene_symbol"].isin(CLOCK)]; idx=clk["row_num"].values-1
E=np.log2(pd.read_csv("data/raw/brainspan/expression_matrix.csv",header=None,index_col=0).values[idx]+1)
def ab(a):
    if "pcw" in a: return "1_prenatal"
    v=float(a.split()[0]); mos=v if a.split()[1]=="mos" else v*12
    if mos<12: return "2_infancy"
    yr=mos/12
    return "3_earlychild" if yr<=5 else ("4_child_adol" if yr<19 else "5_adult")
cols["bin"]=cols["age"].apply(ab)
print("=== H4 firm-up: clock-LGS r by stage (cortex+subcortex, 15 regions), bootstrap 95% CI ===")
stage_r=[]
for b in sorted(cols["bin"].unique()):
    sel=cols[(cols["bin"]==b)&(cols["structure_acronym"].isin(MAP))]
    perreg={r:np.nanmean(E[:,sel[sel["structure_acronym"]==r].index.values]) for r in MAP if len(sel[sel["structure_acronym"]==r])>0}
    regs=[r for r in perreg if np.isfinite(bs_lgs.get(r,np.nan))]
    x=np.array([perreg[r] for r in regs]); y=np.array([bs_lgs[r] for r in regs])
    r0=spearmanr(x,y).statistic
    boot=[spearmanr(x[i],y[i]).statistic for i in [rng.integers(0,len(regs),len(regs)) for _ in range(2000)]]
    lo,hi=np.nanpercentile(boot,[2.5,97.5]); stage_r.append(r0)
    print(f"  {b:14s} n={len(regs):2d}  r={r0:+.3f}  95%CI[{lo:+.2f},{hi:+.2f}]")
# monotonic developmental trend test
tr=spearmanr(range(len(stage_r)),stage_r)
print(f"\nmonotonic trend across 5 stages: Spearman rho={tr.statistic:+.3f}, p={tr.pvalue:.3f}")
