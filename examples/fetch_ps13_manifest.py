"""Download exact external study inputs; fail if public source bytes have changed."""
import sys,json,hashlib,urllib.request
from pathlib import Path
manifest,target=Path(sys.argv[1]),Path(sys.argv[2]);target.mkdir(parents=True,exist_ok=True)
for item in json.loads(manifest.read_text())['files']:
 path=target/item['path'];path.parent.mkdir(parents=True,exist_ok=True)
 if not path.exists():
  temporary=path.with_suffix(path.suffix+'.part');urllib.request.urlretrieve(item['url'],temporary);temporary.replace(path)
 digest=hashlib.sha256()
 with path.open('rb') as file:
  for chunk in iter(lambda:file.read(2**20),b''):digest.update(chunk)
 if path.stat().st_size!=item['bytes'] or digest.hexdigest()!=item['sha256']:raise ValueError(f'Source changed or corrupted: {path}')
 print('Verified',item['path'],flush=True)
