"""Download an additional subject; verify pinned dataset checksums before use."""
from pathlib import Path
import concurrent.futures,hashlib,json,sys,urllib.request,time

def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''): h.update(block)
    return h.hexdigest()

def run(destination):
    out=Path(destination); out.mkdir(parents=True,exist_ok=True)
    sources=json.loads(Path(__file__).with_name('extended_sources.json').read_text())
    for record in sources:
        dest=out/Path(record['path']).name
        print('Download:',dest.name,flush=True)
        if dest.exists() and dest.stat().st_size==record['bytes'] and digest(dest)==record['sha256']: continue
        temporary=dest.with_suffix(dest.suffix+'.part')
        with temporary.open('wb') as f: f.truncate(record['bytes'])
        def chunk(start):
            end=min(start+16*1024*1024,record['bytes'])-1
            for attempt in range(3):
                try:
                    req=urllib.request.Request(record['url'],headers={'Range':f'bytes={start}-{end}'})
                    with urllib.request.urlopen(req,timeout=120) as r:
                        if r.status!=206: raise RuntimeError('Range request was ignored')
                        block=r.read()
                    if len(block)!=end-start+1: raise RuntimeError('Truncated range')
                    with temporary.open('r+b') as f: f.seek(start); f.write(block)
                    return
                except Exception:
                    if attempt==2: raise
                    time.sleep(1)
        with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
            list(pool.map(chunk,range(0,record['bytes'],16*1024*1024)))
        if digest(temporary)!=record['sha256']: raise RuntimeError(f'Checksum mismatch: {dest.name}')
        temporary.replace(dest)
    (out/'extended_manifest.json').write_text(json.dumps(sources,indent=2))
if __name__=='__main__': run(sys.argv[1])
