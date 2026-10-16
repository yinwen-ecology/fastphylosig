"""Read-only verification of the frozen release; Python 3 standard library only."""
from pathlib import Path, PurePosixPath
import argparse, csv, hashlib, json, subprocess, tarfile, zipfile

SOURCE_SHA256 = '8B1A8940101997492D4F3E2AFE696717BFE63CC9905FE91E225A3534754B997D'
base = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--workspace', type=Path, help='Also verify selected working-tree file hashes')
parser.add_argument('--repo', type=Path, help='Git repository for --git-ref')
parser.add_argument('--git-ref', help='Read-only comparison with the future release commit')
args = parser.parse_args()

def sha(data):
    return hashlib.sha256(data).hexdigest().upper()

def safe_name(name):
    p = PurePosixPath(name)
    return not p.is_absolute() and '..' not in p.parts and '\\' not in name and ':' not in name

manifest = json.loads((base/'FREEZE_MANIFEST.json').read_text(encoding='utf-8'))
assert manifest['version'] == '0.3.4'
for e in manifest['files']:
    assert safe_name(e['path']), e['path']
    p = base/e['path']
    assert p.is_file() and not p.is_symlink(), f'Missing/unsafe file: {p}'
    assert len(p.read_bytes()) == e['bytes'] and sha(p.read_bytes()) == e['sha256'], f'Changed: {p}'
listed = {e['path'] for e in manifest['files']}
actual = {p.relative_to(base).as_posix() for p in base.rglob('*') if p.is_file()
          and '__pycache__' not in p.parts and p.name != 'FREEZE_MANIFEST.json'}
assert listed == actual, f'Unexpected/missing files: {listed ^ actual}'
archive = base/'fastphylosig_0.3.4.tar.gz'
assert sha(archive.read_bytes()) == SOURCE_SHA256
with tarfile.open(archive, 'r:gz') as tf:
    items = tf.getmembers()
    names = [m.name for m in items]
    assert len(names) == len(set(names)), 'Duplicate archive names'
    assert all(safe_name(m.name) and (m.isfile() or m.isdir()) for m in items)
    assert all(m.name.startswith('fastphylosig/') or (m.name == 'fastphylosig' and m.isdir()) for m in items)
    desc = tf.extractfile('fastphylosig/DESCRIPTION').read().decode('utf-8')
    assert 'Version: 0.3.4' in desc
    assert not any(PurePosixPath(n).suffix.lower() in {'.dll','.so','.o','.exe'} for n in names)
    with (base/'manifests/package-files.csv').open(encoding='utf-8',newline='') as f:
        for e in csv.DictReader(f):
            assert sha(tf.extractfile(e['path']).read()) == e['sha256'], e['path']
scope = json.loads((base/'STATISTICAL_SCOPE.json').read_text(encoding='utf-8'))
assert scope['heterogeneity_test_default'] is False
assert scope['methods']['lambda']['approximate'] is True
assert scope['methods']['D']['heterogeneity_status'] == 'experimental'
assert scope['methods']['Delta']['heterogeneity_status'] == 'experimental_not_recommended'
assert not scope['methods']['Delta']['recommended_heterogeneity_test']
with zipfile.ZipFile(base/'fastphylosig-0.3.4-release-materials.zip') as z:
    assert z.testzip() is None
    assert len(z.namelist()) == len(set(z.namelist()))
    assert all(safe_name(n) for n in z.namelist())
    assert sha(z.read('fastphylosig_0.3.4.tar.gz')) == SOURCE_SHA256
    assert z.read('STATISTICAL_SCOPE.json') == (base/'STATISTICAL_SCOPE.json').read_bytes()
for line in (base/'SHA256SUMS.txt').read_text(encoding='utf-8').splitlines():
    digest, name = line.split('  ', 1)
    assert sha((base/name).read_bytes()) == digest
with (base/'manifests/selected-source-files.csv').open(encoding='utf-8',newline='') as f:
    source_entries = list(csv.DictReader(f))
if args.workspace:
    for e in source_entries:
        assert sha((args.workspace/e['path']).read_bytes()) == e['sha256'], e['path']
if args.git_ref:
    repo = args.repo or args.workspace or Path.cwd()
    def normalize(data):
        if b'\x00' in data:
            return data
        try:
            return data.decode('utf-8').replace('\r\n','\n').encode('utf-8')
        except UnicodeDecodeError:
            return data
    with tarfile.open(base/'provenance/selected-source-snapshot.tar.gz','r:gz') as tf:
        for e in source_entries:
            expected = tf.extractfile('fastphylosig/'+e['path']).read()
            got = subprocess.run(['git','-C',str(repo),'show',args.git_ref+':'+e['path']],
                                 capture_output=True)
            assert got.returncode == 0, f'Not in candidate commit: {e["path"]}'
            assert normalize(got.stdout) == normalize(expected), f'Commit differs: {e["path"]}'
print(f'PASS: 0.3.4 freeze, {len(manifest["files"])} files, archive/manifests/scope/ZIP/checksums')
if args.workspace:
    print(f'PASS: {len(source_entries)} selected working-tree source files')
if args.git_ref:
    print('PASS: candidate Git commit matches the selected frozen source')
