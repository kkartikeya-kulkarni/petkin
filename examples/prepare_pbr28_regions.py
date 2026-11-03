"""Convert the pinned kinfitr regional PBR28 export into workflow inputs.

Run fetch_pbr28.py and export_pbr28.R first. These are development data with
source-defined anatomical curves, not a new untouched validation cohort.
"""
import argparse,csv,json
from pathlib import Path
import numpy as np

NAMES={'FC':'Frontal cortex','TC':'Temporal cortex','STR':'Striatum','THA':'Thalamus','WB':'Whole brain','CBL':'Cerebellum'}

def convert(data_dir,output_dir):
    data_dir,output_dir=Path(data_dir),Path(output_dir)
    if output_dir.exists() and any(output_dir.iterdir()):raise ValueError('Use an empty output directory')
    output_dir.mkdir(parents=True,exist_ok=True)
    scans=list(csv.DictReader((data_dir/'scans.csv').open()))
    for item in scans:
        scan=item['scan'];folder=output_dir/scan;folder.mkdir();raw=data_dir/scan
        tac=np.atleast_1d(np.genfromtxt(raw/'tacs.csv',names=True,delimiter=','));tac=tac[tac['Duration']>0]
        source=np.genfromtxt(raw/'input.csv',names=True,delimiter=',');time=source['Time'];cp=source['AIF'];cb=source['Blood'];end=float((tac['StartTime']+tac['Duration'])[-1]/60);extension=max(0,end-time[-1])
        if extension:time=np.r_[time,end];cp=np.r_[cp,cp[-1]];cb=np.r_[cb,cb[-1]]
        np.savetxt(folder/'input.csv',np.c_[time,cp,cb],delimiter=',',header='time,plasma,blood',comments='')
        regions=[]
        for key,name in NAMES.items():
            np.savetxt(folder/(key+'.csv'),np.c_[tac['StartTime']/60,(tac['StartTime']+tac['Duration'])/60,tac[key],tac['Weights']],delimiter=',',header='start,end,activity,weight',comments='')
            regions.append(dict(id=key,name=name,path=key+'.csv'))
        (folder/'regions.json').write_text(json.dumps(dict(schema_version=1,scan=scan,regions=regions,source='Pinned kinfitr PBR28 regional export',role='Exposed development data',weighting='Source weights; zero-weight frames retained in curves and excluded from compartment fitting',terminal_input_hold_minutes=extension,anatomy='Source-defined regional TACs; no subject imaging truth claim'),indent=2))
    (output_dir/'scans.json').write_text(json.dumps(scans,indent=2));return scans

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('data_dir');p.add_argument('output_dir');a=p.parse_args();convert(a.data_dir,a.output_dir)
