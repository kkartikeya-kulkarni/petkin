"""Plot saved results, keeping failed boundary fits visible in the annotation."""
from pathlib import Path
import json,sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
folder=Path(sys.argv[1]);rows=json.loads((folder/'pbr28_results.json').read_text())
fig,axes=plt.subplots(1,2,figsize=(11,4))
valid=[r for r in rows if not r['full_fit']['at_bound']]
for r in rows:
    if not r['full_fit']['at_bound']:axes[0].scatter(r['full_fit'].get('raw_VT',r['full_fit']['VT']),r['logan'][1]['VT'])
    axes[1].plot([x['tstar'] for x in r['logan']],[x['VT'] for x in r['logan']],alpha=.65)
limits=[min(r['full_fit'].get('raw_VT',r['full_fit']['VT']) for r in valid),max(r['full_fit'].get('raw_VT',r['full_fit']['VT']) for r in valid)]
axes[0].plot(limits,limits,'k--')
axes[0].set(xlabel='2-tissue VT (mL/mL)',ylabel='Logan VT, t*=40 min (mL/mL)',title='19 interior fits; not a ground-truth comparison')
axes[0].text(.02,.98,'rbqc_1 at bound: 2-tissue VT=10124, Logan VT=2.62\nBoundary fit omitted from plotted range; retained in results',transform=axes[0].transAxes,va='top',fontsize=8,color='red')
axes[1].set(xlabel='Logan t* (min)',ylabel='VT (mL/mL)',title='All 20 real PBR28 scans: t* sensitivity')
fig.tight_layout();fig.savefig(folder/'pbr28_validation.png',dpi=150)
