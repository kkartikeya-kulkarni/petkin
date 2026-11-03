"""Lock and retrieve one immutable OpenNeuro Git snapshot, including metadata."""
import hashlib,json,re,shutil,urllib.parse,urllib.request,xml.etree.ElementTree as ET
from pathlib import Path,PurePosixPath


def _read(url):
    request=urllib.request.Request(url,headers={'User-Agent':'petkin-snapshot'})
    with urllib.request.urlopen(request,timeout=90) as response:return response.read()


def _digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as file:
        for chunk in iter(lambda:file.read(2**20),b''):h.update(chunk)
    return h.hexdigest()


def _safe_path(name):
    p=PurePosixPath(name)
    if p.is_absolute() or '..' in p.parts or not p.parts or '\\' in name:raise ValueError('Unsafe snapshot path')
    return p


def lock_snapshot(dataset,revision,output_path,participant=None):
    """Resolve an explicit revision once; all text and annex hashes use that commit."""
    if not re.fullmatch(r'ds\d{6}',dataset):raise ValueError('Expected an OpenNeuro ID such as ds004868')
    if not revision or not re.fullmatch(r'[A-Za-z0-9._/-]+',revision):raise ValueError('An explicit dataset revision is required')
    if participant is not None:
        participant=participant.removeprefix('sub-')
        if not re.fullmatch(r'[A-Za-z0-9]+',participant):raise ValueError('Invalid participant label')
    out=Path(output_path)
    if out.exists():raise ValueError('Snapshot lock already exists; use a new output file')
    api=f'https://api.github.com/repos/OpenNeuroDatasets/{dataset}'
    commit=json.loads(_read(api+'/commits/'+urllib.parse.quote(revision,safe='')))['sha']
    tree=json.loads(_read(api+'/git/trees/'+commit+'?recursive=1'))
    if tree.get('truncated'):raise ValueError('Git tree is truncated; cannot create a complete lock')
    base=f'https://raw.githubusercontent.com/OpenNeuroDatasets/{dataset}/{commit}/';files=[]
    for item in tree['tree']:
        name=item['path'];_safe_path(name)
        if item['type']!='blob' or name.startswith('.') or '/.' in name:continue
        if participant and name.startswith('sub-') and not name.startswith('sub-'+participant+'/'):continue
        entry=dict(path=name,git_blob_sha=item['sha'],url=base+urllib.parse.quote(name,safe='/'))
        if item['mode']=='120000':
            link=_read(entry['url']);actual=hashlib.sha1(b'blob '+str(len(link)).encode()+b'\0'+link).hexdigest()
            if actual!=entry['git_blob_sha']:raise ValueError('Git symlink hash mismatch')
            match=re.search(rb'SHA256E-s(\d+)--([a-f0-9]{64})',link)
            if match is None:raise ValueError('Unsupported annex object: '+name)
            entry.update(kind='annex',bytes=int(match[1]),sha256=match[2].decode())
        else:entry.update(kind='git',bytes=item.get('size'))
        files.append(entry)
    if not files or participant and not any(x['path'].startswith('sub-'+participant+'/') for x in files):raise ValueError('No matching participant files')
    manifest=dict(schema_version=1,dataset=dataset,requested_revision=revision,commit=commit,participant=participant,files=files,
        policy='Every metadata/blood file comes from this Git commit; any renamed image storage object must match the pinned annex content hash')
    out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(manifest,indent=2));return manifest


def _storage_objects(dataset):
    ns={'s':'http://s3.amazonaws.com/doc/2006-03-01/'};items=[];token=None
    while True:
        query={'list-type':'2','prefix':dataset+'/','max-keys':'1000'}
        if token:query['continuation-token']=token
        root=ET.fromstring(_read('https://s3.amazonaws.com/openneuro.org?'+urllib.parse.urlencode(query)))
        items.extend(dict(key=x.findtext('s:Key',namespaces=ns),bytes=int(x.findtext('s:Size',namespaces=ns))) for x in root.findall('s:Contents',ns))
        if root.findtext('s:IsTruncated',namespaces=ns)!='true':break
        token=root.findtext('s:NextContinuationToken',namespaces=ns)
        if not token:raise ValueError('Truncated S3 listing without continuation token')
    return items


def fetch_snapshot(lock_path,output_dir,cache_dir=None):
    """Retrieve locked content; mutable S3 JSON/blood is never substituted.

    A renamed image can be read from local cache or S3 only after its SHA256 is
    verified. The lock's participant names and metadata always remain together.
    """
    lock_path=Path(lock_path);lock=json.loads(lock_path.read_text());out=Path(output_dir)
    if out.exists() and any(out.iterdir()):raise ValueError('Snapshot output directory must be empty')
    names=[x['path'] for x in lock['files']]
    if len(set(names))!=len(names):raise ValueError('Duplicate paths in snapshot lock')
    for name in names:_safe_path(name)
    if not re.fullmatch(r'ds\d{6}',lock['dataset']) or not re.fullmatch(r'[a-f0-9]{40}',lock['commit']):raise ValueError('Invalid snapshot identity')
    out.mkdir(parents=True,exist_ok=True);records=[];storage=None;cache={}
    if cache_dir:
        for p in Path(cache_dir).rglob('*'):
            if p.is_file() and p.name.endswith(('.nii','.nii.gz')):cache.setdefault(p.stat().st_size,[]).append(p)
    base=f"https://raw.githubusercontent.com/OpenNeuroDatasets/{lock['dataset']}/{lock['commit']}/"
    for item in lock['files']:
        destination=out/item['path'];destination.parent.mkdir(parents=True,exist_ok=True);temporary=destination.with_name(destination.name+'.part');source=None
        try:
            if item['kind']=='git':
                url=base+urllib.parse.quote(item['path'],safe='/');data=_read(url)
                if hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()!=item['git_blob_sha']:raise ValueError('Pinned Git blob differs: '+item['path'])
                temporary.write_bytes(data);source=url
            elif item['kind']=='annex':
                if not re.fullmatch(r'[a-f0-9]{64}',item['sha256']):raise ValueError('Invalid annex hash')
                for candidate in cache.get(item['bytes'],[]):
                    if _digest(candidate)==item['sha256']:shutil.copyfile(candidate,temporary);source='SHA256-verified local cache';break
                if source is None:
                    if storage is None:storage=_storage_objects(lock['dataset'])
                    candidates=[x for x in storage if x['bytes']==item['bytes']]
                    # Content, not participant spelling or file size alone, decides equivalence.
                    for candidate in candidates:
                        url='https://s3.amazonaws.com/openneuro.org/'+urllib.parse.quote(candidate['key'],safe='/')
                        with urllib.request.urlopen(url,timeout=90) as incoming,temporary.open('wb') as file:shutil.copyfileobj(incoming,file)
                        if _digest(temporary)==item['sha256']:source=url;break
                    if source is None:raise ValueError('No accessible object matches pinned content: '+item['path'])
                if temporary.stat().st_size!=item['bytes'] or _digest(temporary)!=item['sha256']:raise ValueError('Annex integrity failure')
            else:raise ValueError('Unknown snapshot object kind')
            temporary.replace(destination);records.append(dict(path=item['path'],source=source,sha256=_digest(destination),bytes=destination.stat().st_size))
        finally:
            if temporary.exists():temporary.unlink()
    result=dict(dataset=lock['dataset'],commit=lock['commit'],requested_revision=lock['requested_revision'],lock_sha256=_digest(lock_path),files=records,all_files_verified=True)
    (out/'SNAPSHOT_PROVENANCE.json').write_text(json.dumps(result,indent=2));return result
