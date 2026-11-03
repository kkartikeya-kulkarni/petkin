import urllib.request,json,hashlib
from pathlib import Path
import sys
out=Path(sys.argv[1] if len(sys.argv)>1 else 'public-pet'); out.mkdir(parents=True,exist_ok=True); manifest=[]
for sub in ['s028','s029']:
 prefix=f'sub-{sub}/ses-open/'
 for rel in [f'pet/sub-{sub}_ses-open_task-rest_pet.json',f'pet/sub-{sub}_ses-open_task-rest_recording-autosampler_blood.json',f'pet/sub-{sub}_ses-open_task-rest_recording-autosampler_blood.tsv',f'pet/sub-{sub}_ses-open_task-rest_pet.nii.gz',f'anat/sub-{sub}_ses-open_T1w.nii.gz']:
  path=prefix+rel; url='https://s3.amazonaws.com/openneuro.org/ds004513/'+path
  dest=out/rel.split('/')[-1]
  print('Fetching',url,flush=True)
  if not dest.exists():
   temporary=dest.with_suffix(dest.suffix+'.part')
   urllib.request.urlretrieve(url,temporary)
   temporary.replace(dest)
  manifest.append(dict(path=path,url=url,bytes=dest.stat().st_size,sha256=hashlib.sha256(dest.read_bytes()).hexdigest()))
(out/'public_manifest.json').write_text(json.dumps(manifest,indent=2))
