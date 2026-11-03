"""Auditable PET kinetic analysis; all CSV time columns use minutes."""
import argparse,json
from pathlib import Path
import numpy as np
from . import Frames,InputFunction,fit_compartment,logan,ma1,patlak,profile_likelihood,parametric_bootstrap,profile_macroparameter,analyze_research_tac
from .reporting import file_manifest,write_qc_html


def _csv(tac_path,input_path):
    tac=np.atleast_1d(np.genfromtxt(tac_path,delimiter=',',names=True))
    data=np.atleast_1d(np.genfromtxt(input_path,delimiter=',',names=True))
    f=Frames(tac['start'],tac['end']);inp=InputFunction(data['time'],data['plasma'],data['blood'] if 'blood' in data.dtype.names else None)
    return f,inp,tac


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    for name in ['fit','research-fit','logan','ma1','patlak','profile','bootstrap']:
        q=sub.add_parser(name);q.add_argument('tac');q.add_argument('input');q.add_argument('output')
        q.add_argument('--report',help='Write portable HTML quality summary')
        if name in ['profile','bootstrap']:
            q.add_argument('--fit-result',required=True)
            if name=='profile':
                q.add_argument('--parameter',required=True);q.add_argument('--minimum',type=float,required=True)
                q.add_argument('--maximum',type=float,required=True);q.add_argument('--points',type=int,default=41)
                q.add_argument('--spacing',choices=['log','linear'],default='log')
            else:q.add_argument('--samples',type=int,default=100);q.add_argument('--seed',type=int,default=0)
        else:
            q.add_argument('--vb',type=float,default=0.);q.add_argument('--delay',type=float,default=0.)
            q.add_argument('--input-tail',choices=['error','hold'],default='error')
            if name in ('fit','research-fit'):
                q.add_argument('--model',choices=['1t','2t','2t_irreversible'],default='2t')
                q.add_argument('--fit-vb',action='store_true');q.add_argument('--fit-delay',action='store_true')
                q.add_argument('--starts',type=int,default=8);q.add_argument('--seed',type=int,default=0)
                q.add_argument('--absolute-sigma',action='store_true');q.add_argument('--dt',type=float,default=.025)
                if name=='research-fit':q.add_argument('--assumptions',required=True,help='JSON study records confirming required input assumptions')
            else:q.add_argument('--tstar',type=float,required=True)
    q=sub.add_parser('preprocess');q.add_argument('pet');q.add_argument('mri');q.add_argument('timing');q.add_argument('output_dir')
    q.add_argument('--reference',type=int);q.add_argument('--preserve-blank',action='store_true')
    args=parser.parse_args()
    if args.command=='preprocess':
        import SimpleITK as sitk
        from .imaging import preprocess_dynamic
        from .bids import frames_from_bids
        out=Path(args.output_dir);out.mkdir(parents=True,exist_ok=True)
        pet,mri=sitk.ReadImage(args.pet),sitk.ReadImage(args.mri);metadata=json.loads(Path(args.timing).read_text())
        frames=frames_from_bids(metadata)
        if len(frames.start)!=pet.GetSize()[3]:raise ValueError('Image frame count does not match BIDS timing')
        result=preprocess_dynamic(pet,mri,frames.duration,args.reference,'identity' if args.preserve_blank else 'raise')
        sitk.WriteImage(result['image'],str(out/'pet_in_mri.nii.gz'))
        sitk.WriteImage(result['motion_corrected'],str(out/'motion_corrected.nii.gz'))
        for i,tx in enumerate(result['frame_transforms']):sitk.WriteTransform(tx,str(out/f'reference_to_frame_{i:03d}.tfm'))
        for i,tx in enumerate(result['composed_transforms']):
            tx.FlattenTransform();sitk.WriteTransform(tx,str(out/f'mri_to_original_frame_{i:03d}.h5'))
        sitk.WriteTransform(result['mri_to_reference'],str(out/'mri_to_pet_reference.tfm'))
        qc={k:v for k,v in result.items() if k in ['motion_qc','registration_metric','registration_stop','registration_candidates','final_interpolations_per_frame','limitations']}
        qc['provenance']=file_manifest([args.pet,args.mri,args.timing],vars(args))
        (out/'qc.json').write_text(json.dumps(qc,indent=2,allow_nan=False));return
    f,inp,tac=_csv(args.tac,args.input)
    sigma=tac['sigma'] if 'sigma' in tac.dtype.names else None
    if args.command in ('fit','research-fit'):
        options={}
        function=fit_compartment
        if args.command=='research-fit':
            function=analyze_research_tac;options['assumptions']=json.loads(Path(args.assumptions).read_text())
        result=function(f,inp,tac['activity'],args.model,sigma,vb=args.vb,delay=args.delay,
                               input_tail=args.input_tail,fit_vb=args.fit_vb,fit_delay=args.fit_delay,
                               n_starts=args.starts,seed=args.seed,absolute_sigma=args.absolute_sigma,dt=args.dt,**options)
    elif args.command in ['profile','bootstrap']:
        fit=json.loads(Path(args.fit_result).read_text())
        if args.command=='profile':
            if args.maximum<=args.minimum or args.points<3 or (args.spacing=='log' and args.minimum<=0):raise ValueError('Ordered profile grid required; logarithmic grids must be positive')
            grid=(np.geomspace if args.spacing=='log' else np.linspace)(args.minimum,args.maximum,args.points)
            result=profile_macroparameter(f,inp,tac['activity'],fit,grid,sigma=sigma) if args.parameter in ('VT','Ki') else profile_likelihood(f,inp,tac['activity'],fit,args.parameter,grid,sigma=sigma)
        else:result=parametric_bootstrap(f,inp,tac['activity'],fit,sigma=sigma,n=args.samples,seed=args.seed)
    else:
        method={'logan':logan,'ma1':ma1,'patlak':patlak}[args.command]
        result=method(f,tac['activity'],inp,args.tstar,vb=args.vb,delay=args.delay,input_tail=args.input_tail)
    paths=[args.tac,args.input]+([args.assumptions] if args.command=='research-fit' else [])+([args.fit_result] if args.command in ['profile','bootstrap'] else [])
    result['provenance']=file_manifest(paths,vars(args))
    Path(args.output).write_text(json.dumps(result,indent=2,allow_nan=False))
    if args.report:write_qc_html(result,args.report)

if __name__=='__main__':main()
