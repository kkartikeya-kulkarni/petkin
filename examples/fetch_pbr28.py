"""Fetch the exact kinfitr example dataset used in the real-data report."""
from pathlib import Path
import hashlib,json,sys,urllib.request
source=json.loads(Path(__file__).with_name('pbr28_source.json').read_text())
out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=True)
path=out/'pbr28.RData';temporary=out/'pbr28.RData.part'
urllib.request.urlretrieve(source['url'],temporary)
if hashlib.sha256(temporary.read_bytes()).hexdigest()!=source['sha256']: raise RuntimeError('Checksum mismatch')
temporary.replace(path)
(out/'source.json').write_text(json.dumps(source,indent=2))
