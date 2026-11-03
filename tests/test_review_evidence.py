"""Synthetic tests for review scoring. These are not human validation evidence."""
import csv,hashlib,importlib.util,json
from pathlib import Path
import pytest
import SimpleITK as sitk
spec=importlib.util.spec_from_file_location('review',Path(__file__).parents[1]/'scripts/evaluate_registration_review.py');review=importlib.util.module_from_spec(spec);spec.loader.exec_module(review)

def setup(tmp_path):
 p=tmp_path/'protocol.json';p.write_text(json.dumps(dict(protocol_id='SYNTHETIC_TEST_ONLY',reviewed_by='test fixture',locked_at='2026-01-01T00:00:00Z',reviewer_statement='synthetic evaluator test',coordinate_convention='LPS_mm',independent_of_algorithm_outputs=True,limits=dict(mean_of_pair_mean_TRE_mm=1,maximum_landmark_TRE_mm=2))))
 tx=tmp_path/'transform.tfm';transform=sitk.TranslationTransform(3,[1.,2.,3.]);sitk.WriteTransform(transform,str(tx));sha=hashlib.sha256(tx.read_bytes()).hexdigest()
 cases=[dict(scan=f'case{i}',subject=f'person{i}',measured_PET_and_MRI=True,transform_path=tx.name,transform_sha256=sha) for i in range(10)]
 m=tmp_path/'cases.json';m.write_text(json.dumps(dict(protocol_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),cases=cases)))
 points=tmp_path/'landmarks.csv'
 with points.open('w',newline='') as file:
  w=csv.writer(file);w.writerow(['scan','landmark_id','MRI_LPS_x_mm','MRI_LPS_y_mm','MRI_LPS_z_mm','PET_LPS_x_mm','PET_LPS_y_mm','PET_LPS_z_mm','rater_id','anatomical_definition'])
  for case in cases:
   for i,(a,b) in enumerate([([0,0,0],[1,2,3]),([10,0,0],[11,2,3]),([0,10,0],[1,12,3])]):w.writerow([case['scan'],str(i),*a,*b,'synthetic rater','synthetic landmark'])
 return p,points,m

def test_known_translation_scored_without_certifying_toolbox(tmp_path):
 result=review.evaluate(*setup(tmp_path));assert result['pairs']==10 and result['mean_of_pair_mean_TRE_mm']==0
 assert result['geometric_criterion_passed'] and not result['toolbox_qualified'] and result['independent_review_still_required']

def test_inverted_transform_fails_geometric_criterion(tmp_path):
 p,points,m=setup(tmp_path);tx=tmp_path/'transform.tfm';sitk.WriteTransform(sitk.TranslationTransform(3,[-1.,-2.,-3.]),str(tx));manifest=json.loads(m.read_text())
 for case in manifest['cases']:case['transform_sha256']=hashlib.sha256(tx.read_bytes()).hexdigest()
 m.write_text(json.dumps(manifest));result=review.evaluate(p,points,m);assert not result['geometric_criterion_passed']

@pytest.mark.parametrize('problem',['blank','unreviewed','hash','small','collinear','future','bad_timestamp'])
def test_incomplete_or_invalid_evidence_rejected(tmp_path,problem):
 p,points,m=setup(tmp_path)
 if problem=='blank':points.write_text(points.read_text().replace('0,0,0,1,2,3',',0,0,1,2,3',1))
 if problem=='unreviewed':
  value=json.loads(p.read_text());value['limits']['maximum_landmark_TRE_mm']=None;p.write_text(json.dumps(value))
 if problem in ['hash','small']:
  value=json.loads(m.read_text())
  if problem=='hash':value['cases'][0]['transform_sha256']='0'*64
  else:value['cases']=value['cases'][:9]
  m.write_text(json.dumps(value))
 if problem in ['future','bad_timestamp']:
  value=json.loads(p.read_text());value['locked_at']='2099-01-01T00:00:00Z' if problem=='future' else 'not a date';p.write_text(json.dumps(value))
 if problem=='collinear':points.write_text(points.read_text().replace('0,10,0,1,12,3','20,0,0,21,2,3'))
 with pytest.raises(ValueError):review.evaluate(p,points,m)
