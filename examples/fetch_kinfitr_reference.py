"""Retrieve the exact upstream functions used for independent point comparisons."""
from pathlib import Path
import hashlib,json,sys,urllib.request
manifest=json.loads(Path(__file__).with_name('kinfitr_reference_sources.json').read_text())
out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=True)
for record in manifest['files']:
    data=urllib.request.urlopen(record['url'],timeout=120).read()
    if hashlib.sha256(data).hexdigest()!=record['sha256']:raise RuntimeError('Reference source checksum mismatch')
    (out/record['path']).write_bytes(data)
print('Retrieved unmodified sources at commit',manifest['commit'])
