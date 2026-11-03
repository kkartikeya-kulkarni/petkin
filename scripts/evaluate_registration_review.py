"""Calculate TRE from independently supplied review evidence; never invent truth.

python scripts/evaluate_registration_review.py PROTOCOL.json LANDMARKS.csv CASES.json OUT.json
Reviewer declarations are documentary attestations, not verified identities.
"""
import argparse,csv,hashlib,json
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import SimpleITK as sitk
from petkin.imaging import landmark_registration_error

def evaluate(protocol_path,landmarks_path,cases_path):
 protocol_path,landmarks_path,cases_path=map(Path,(protocol_path,landmarks_path,cases_path))
 protocol_bytes=protocol_path.read_bytes();protocol=json.loads(protocol_bytes)
 required=['protocol_id','reviewed_by','locked_at','reviewer_statement']
 if any(not isinstance(protocol.get(k),str) or not protocol[k].strip() for k in required):raise ValueError('A reviewed, locked protocol is required')
 try:locked=datetime.fromisoformat(protocol['locked_at'].replace('Z','+00:00'))
 except ValueError as e:raise ValueError('Protocol lock time must be an ISO timestamp') from e
 if locked.tzinfo is None or locked>datetime.now(timezone.utc):raise ValueError('Protocol lock time must include timezone and precede evaluation')
 if protocol.get('coordinate_convention')!='LPS_mm' or protocol.get('independent_of_algorithm_outputs') is not True:raise ValueError('Independent LPS-mm landmark declaration required')
 limits=protocol.get('limits',{});thresholds=[]
 for key in ['mean_of_pair_mean_TRE_mm','maximum_landmark_TRE_mm']:
  value=limits.get(key)
  if isinstance(value,bool) or not isinstance(value,(float,int)) or not np.isfinite(value) or value<=0:raise ValueError('Reviewed positive finite error limits are required')
  thresholds.append(float(value))
 manifest=json.loads(cases_path.read_text());cases=manifest.get('cases',[])
 if manifest.get('protocol_sha256')!=hashlib.sha256(protocol_bytes).hexdigest():raise ValueError('Cases must bind to the locked protocol hash')
 if len(cases)<10 or len({x['scan'] for x in cases})!=len(cases):raise ValueError('At least ten unique measured pairs required')
 with landmarks_path.open(newline='') as file:rows=list(csv.DictReader(file))
 groups={}
 for row in rows:
  if row['scan'] not in {x['scan'] for x in cases}:raise ValueError('Unlisted scan in landmarks')
  groups.setdefault(row['scan'],[]).append(row)
 results=[]
 for case in cases:
  scan=case['scan'];points=groups.get(scan,[])
  if len(points)<3:raise ValueError('At least three landmarks per scan required: '+scan)
  if len({r['landmark_id'] for r in points})!=len(points):raise ValueError('Duplicate landmark IDs: '+scan)
  if any(not r.get('rater_id','').strip() or not r.get('anatomical_definition','').strip() for r in points):raise ValueError('Rater and landmark definitions required')
  if case.get('measured_PET_and_MRI') is not True or not case.get('subject'):raise ValueError('Measured image pair and subject identity declarations required')
  try:
   fixed=np.array([[float(r['MRI_LPS_'+axis+'_mm']) for axis in 'xyz'] for r in points]);moving=np.array([[float(r['PET_LPS_'+axis+'_mm']) for axis in 'xyz'] for r in points])
  except (KeyError,ValueError) as e:raise ValueError('Blank or invalid independent coordinates: '+scan) from e
  if not np.isfinite(fixed).all() or not np.isfinite(moving).all():raise ValueError('Finite coordinates required')
  if min(np.linalg.matrix_rank(p-p.mean(axis=0)) for p in [fixed,moving])<2:raise ValueError('Collinear or collapsed landmarks: '+scan)
  transform_path=(cases_path.parent/case['transform_path']).resolve()
  if hashlib.sha256(transform_path.read_bytes()).hexdigest()!=case['transform_sha256']:raise ValueError('Transform checksum mismatch: '+scan)
  metric=landmark_registration_error(fixed,moving,sitk.ReadTransform(str(transform_path)))
  results.append(dict(scan=scan,subject=case['subject'],**metric))
 mean=float(np.mean([x['mean_error_mm'] for x in results]));maximum=float(max(x['maximum_error_mm'] for x in results))
 return dict(protocol_id=protocol['protocol_id'],protocol_sha256=hashlib.sha256(protocol_bytes).hexdigest(),landmarks_sha256=hashlib.sha256(landmarks_path.read_bytes()).hexdigest(),cases_manifest_sha256=hashlib.sha256(cases_path.read_bytes()).hexdigest(),pairs=len(results),subjects=len({x['subject'] for x in results}),mean_of_pair_mean_TRE_mm=mean,maximum_landmark_TRE_mm=maximum,geometric_criterion_passed=mean<=thresholds[0] and maximum<=thresholds[1],results=results,toolbox_qualified=False,independent_review_still_required=True,limitation='Scores are conditional on the supplied landmarks. Code cannot certify anatomical correctness, source declarations, reviewer independence or the historical timing of protocol locking.')

def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('protocol');parser.add_argument('landmarks');parser.add_argument('cases');parser.add_argument('output');a=parser.parse_args()
 try:result=evaluate(a.protocol,a.landmarks,a.cases)
 except (ValueError,KeyError,OSError,RuntimeError) as e:parser.exit(2,'Evidence rejected: '+str(e)+'\n')
 output=Path(a.output)
 if output.exists():parser.exit(2,'Refusing to overwrite an existing evaluation.\n')
 output.write_text(json.dumps(result,indent=2,allow_nan=False));print('Geometric criterion:',result['geometric_criterion_passed'],'; independent review still required.')
if __name__=='__main__':main()
