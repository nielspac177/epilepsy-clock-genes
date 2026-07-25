import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, nibabel as nib
from nilearn.image import resample_to_img
from scipy.stats import spearmanr
from scipy.spatial.distance import cdist
import abagen
from brainsmash.mapgen.base import Base
from neuromaps import datasets
SEED=20260724
# curated MNI152 PET/metabolism maps (source,desc,label) — epilepsy/E-I relevant
MAPS=[("norgaard2021","flumazenil","GABAa (flumazenil)"),("dukart2018","flumazenil","GABAa/BZ (flumazenil)"),
 ("hillmer2016","flubatine","nAChR a4b2"),("tuominen","feobv","AChVAChT (FEOBV)"),
 ("dukart2018","fdopa","dopamine synth") ,("kaller2017","sch23390","D1 (SCH23390)"),("sandiego2015","flb457","D2 (FLB457)"),
 ("savli2012","dasb","SERT (DASB)"),("savli2012","way100635","5-HT1A"),("savli2012","altanserin","5-HT2A"),
 ("normandin2015","omar","CB1 (OMAR)"),("kantonen2020","carfentanil","mu-opioid"),
 ("castrillon2023","cmrglc","glucose metab"),("satterthwaite2014","meancbf","cerebral blood flow")]
aimg=nib.load(str(abagen.fetch_desikan_killiany(native=False)["image"]))
ad=np.asarray(aimg.dataobj); labels=[int(l) for l in np.unique(ad) if l!=0]
cent=pd.read_csv("results/spatial/ahba_dk_centroids.csv").set_index("id")
def parc(img):
    r=resample_to_img(img,aimg,interpolation="continuous"); d=np.asarray(r.dataobj)
    return pd.Series({l:float(np.nanmean(d[ad==l])) for l in labels})
lgs=parc(nib.load("EEGfMRI_GPFA_NOGAONLY_NOTHRESH_TSCORES_MNI2009b.nii.gz"))
reg=[l for l in labels if np.isfinite(lgs[l])]
lv=lgs[reg].to_numpy()
D=cdist(cent.loc[reg,["x","y","z"]].to_numpy(),cent.loc[reg,["x","y","z"]].to_numpy())
np.random.seed(SEED); surr=Base(x=lv,D=D,resample=True)(n=1000)
def sp(vec):
    ok=np.isfinite(vec); r=spearmanr(vec[ok],lv[ok]).statistic
    rn=np.array([spearmanr(vec[ok],s[ok]).statistic for s in surr])
    return r,(np.sum(np.abs(rn)>=abs(r))+1)/1001
rows=[]
for src,desc,lab in MAPS:
    try:
        f=datasets.fetch_annotation(source=src,desc=desc,space="MNI152")
        v=parc(nib.load(f) if isinstance(f,str) else f).reindex(reg).to_numpy()
        r,p=sp(v); rows.append((lab,r,p)); print(f"  {lab:22s} r={r:+.3f} p={p:.3f}")
    except Exception as e:
        print(f"  {lab:22s} SKIP ({str(e)[:50]})")
df=pd.DataFrame(rows,columns=["map","r","p_spin"])
ps=df["p_spin"].values; n=len(ps); o=np.argsort(ps); q=np.empty(n); q[o]=np.minimum.accumulate((ps[o]*n/(np.arange(n)+1))[::-1])[::-1]; df["p_fdr"]=np.clip(q,0,1)
df.to_csv("results/spatial/lgs_receptor_fingerprint.tsv",sep="\t",index=False)
print("\n"+df.sort_values("p_spin").to_string(index=False)); print("DONE")
