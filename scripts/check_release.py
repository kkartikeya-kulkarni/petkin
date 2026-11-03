"""Verify historical source hashes and the consistency of published evidence."""
import hashlib,json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
def run():
 base=ROOT/'validation/v03';freeze=json.loads((base/'IMPLEMENTATION_FREEZE.json').read_text());amendments={}
 for filename,key in [('REGISTRATION_GENERATOR_AMENDMENT.json','sha256'),('VERSION_METADATA_AMENDMENT.json','current_sha256')]:
  entry=json.loads((base/filename).read_text());amendments[entry['file']]=entry[key]
 addendum=ROOT/'validation/workflow_v04/IMPLEMENTATION_ADDENDUM.json'
 if addendum.exists():
  for item in json.loads(addendum.read_text())['files']:
   if item['path'] in ['src/petkin/cli.py','src/petkin/__init__.py']:amendments[item['path']]=item['sha256']
 interop=ROOT/'validation/bids_interop_v04/IMPLEMENTATION_ADDENDUM.json'
 if interop.exists():
  for item in json.loads(interop.read_text())['files']:
   actual=hashlib.sha256((ROOT/item['path']).read_bytes()).hexdigest()
   if actual!=item['sha256']:raise ValueError('Interop implementation differs: '+item['path'])
   amendments[item['path']]=item['sha256']
 for entry in freeze['files']:
  actual=hashlib.sha256((ROOT/entry['path']).read_bytes()).hexdigest()
  if actual not in {entry['sha256'],amendments.get(entry['path'])}:raise ValueError('Frozen implementation changed: '+entry['path'])
 ext=ROOT/'validation/external_ps13';summary=json.loads((ext/'summary.json').read_text());results=json.loads((ext/'results.json').read_text())
 if len(results)!=20 or any(r['status']!='completed' for r in results):raise ValueError('Expected the full 20-scan cohort')
 for mode in ['fixed_delay','fitted_delay']:
  errors=np.array([r['cross_validation'][mode]['held_frame_nrmse'] for r in results])
  for key,actual in [('median_nrmse',np.median(errors)),('q95_nrmse',np.quantile(errors,.95)),('maximum_nrmse',np.max(errors))]:
   if not np.isclose(summary[mode][key],actual,rtol=1e-12):raise ValueError('Summary inconsistent: '+mode+'/'+key)
 if summary['research_released']!=0 or any(r['full_research_fit']['usable_for_quantification'] for r in results):raise ValueError('Unexpected research release')
 qualification=json.loads((ext/'current_qualification.json').read_text())
 if qualification['qualified_for_declared_research_scope'] or set(qualification['incomplete_gates'])!={'EXT-01','REG-02'}:raise ValueError('Unresolved qualification claims changed')
 print('Numerical engine hashes and documented workflow amendments verified. Historical external summaries are consistent.')
 print('Full research qualification: INCOMPLETE (EXT-01, REG-02). CI success does not certify validation.')
if __name__=='__main__':run()
