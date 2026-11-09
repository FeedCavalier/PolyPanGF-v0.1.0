"""End-to-end release checks. Run on the configured Linux analysis server."""
from collections import Counter
from pathlib import Path
import json
import os
import subprocess
import sys
import tempfile

from polypangf.core import read_tsv, save_json, sha256
from polypangf.cli import FUNCTIONS

root = Path(__file__).resolve().parents[1]
python = root / '.venv/bin/python'
environment = os.environ.copy()
environment['PATH'] = str(root / '.venv/bin') + os.pathsep + environment['PATH']
qa = root / 'validation/release_qa';qa.mkdir(exist_ok=True)


def run(arguments, name):
    with (qa / (name + '.stdout')).open('w') as output, (qa / (name + '.stderr')).open('w') as error:
        subprocess.run([str(arg) for arg in arguments], cwd=root, env=environment, stdout=output, stderr=error, check=True)
    return (qa / (name + '.stdout')).read_text()


run([python, '-m', 'unittest', 'discover', '-s', 'tests', '-v'], 'unit_tests')
run([root / '.venv/bin/ruff', 'check', 'src', 'tests', '--select', 'F'], 'lint')
run([root / '.venv/bin/ruff', 'format', '--check', 'src', 'tests'], 'format')
for name in ['init', 'doctor', 'validate', 'status', 'export', 'run'] + list(FUNCTIONS):
    run([python, '-m', 'polypangf', name, '--help'], 'help_' + name)
run([python, '-m', 'polypangf', '--version'], 'version')
plan_target = root / 'demo/release_dry_run'
assert not plan_target.exists()
run([python, '-m', 'polypangf', 'run', '-c', root / 'demo/config.json', '-o', 'release_dry_run', '--dry-run'], 'dry_run')
assert not plan_target.exists()

configs = {'demo': root/'demo/config.json', 'strawberry_replay': root/'validation/strawberry_replay/config.json',
           'native_family': root/'validation/native_family/config.json', 'psonic_native': root/'validation/psonic_native/config.json',
           'kaks_real': root/'validation/kaks_real/config.json'}
psonic = json.loads(configs['psonic_native'].read_text())
prepared = dict(psonic);prepared['output'] = 'results_prepared'
prepared['orthology'] = dict(psonic['orthology'], prefix=str(configs['psonic_native'].parent/'results/synteny/genome'),
    groups=str(configs['psonic_native'].parent/'results/synteny/Orthogroups.tsv'),
    sequence_ids=str(configs['psonic_native'].parent/'results/synteny/SequenceIDs.txt'),
    species_ids=str(configs['psonic_native'].parent/'results/synteny/SpeciesIDs.txt'),
    ploidy=str(configs['psonic_native'].parent/'results/synteny/relative_ploidy.tsv'))
prepared_path=configs['psonic_native'].parent/'config_prepared.json'
save_json(prepared_path,prepared);configs['psonic_prepared']=prepared_path
summary={}
for name, config in configs.items():
    print('Checking',name,flush=True)
    run([python,'-m','polypangf','run','-c',config,'--resume'],name+'_run')
    output=config.parent/json.loads(config.read_text()).get('output','results')
    before={str(path.relative_to(output)):sha256(path) for path in output.glob('*/state.json')}
    run([python,'-m','polypangf','run','-c',config,'--resume'],name+'_resume')
    after={str(path.relative_to(output)):sha256(path) for path in output.glob('*/state.json')}
    assert before==after,(name,'unexpected recomputation on unchanged input')
    status=json.loads(run([python,'-m','polypangf','status','-c',config,'--verify'],name+'_status'))
    assert all(value.get('outputs_valid',True) for key,value in status.items() if key!='failures')
    assert not status['failures']
    summary[name]={'pan':json.loads((output/'pan/summary.json').read_text()),'resume_verified':True,'output_hashes_verified':True}
    if (output/'enrichment/summary.json').exists():summary[name]['enrichment']=json.loads((output/'enrichment/summary.json').read_text())

replay=summary['strawberry_replay']['pan']
assert (replay['groups'],replay['assigned_loci'],replay['unassigned_loci'])==(122,1494,51)
assert replay['classes']=={'Core':5,'Softcore':10,'Shell':102,'Specific':5}
assert summary['strawberry_replay']['enrichment']['significant']==19
enrich=read_tsv(root/'validation/strawberry_replay/results/enrichment/significant.tsv')
assert Counter(row['scope'] for row in enrich)=={'all_family':5,'coverage_classes':14}
marker_state=json.loads((root/'validation/strawberry_replay/results/markers/state.json').read_text())
assert marker_state['summary']=={'matched_records':145,'unmatched_records':0,'marker_genes':28}
native=read_tsv(root/'validation/native_family/results/identify/candidates.tsv')
assert Counter(row['dataset'] for row in native)=={'F006':57,'F017':85}
for name in ['psonic_native','psonic_prepared']:
    assert summary[name]['pan']['groups']==6 and summary[name]['pan']['assigned_loci']==42
    assert summary[name]['pan']['classes']=={'Core':6}
expected={row['orthogroup']:row for row in json.loads((root/'validation/kaks_real/expected_pairs.json').read_text())}
observed=read_tsv(root/'validation/kaks_real/results/kaks/pairs.tsv')
assert len(observed)==3 and all(row['status']=='pass' for row in observed)
for row in observed:
    for current,old in [('Ka','Ka'),('Ks','Ks'),('omega','Omega')]:
        assert abs(float(row[current])-float(expected[row['group']][old]))<.00011,(row,old)
with tempfile.TemporaryDirectory() as temporary:
    archive=Path(temporary)/'demo.zip'
    exported=json.loads(run([python,'-m','polypangf','export','-c',configs['demo'],'-o',archive],'export'))
    assert archive.is_file() and exported['files']>20
summary['official_psonic_commit']='bcebe9af8d6bf3cfee0779e6e891a82b10ffe486'
summary['official_psonic_script_sha256']=sha256(root/'.external/pSONIC.py')
summary['unit_tests']=15
summary['cli_subcommands']=len(FUNCTIONS)+6
summary['scope']='Functional and numerical regression checks; no comparative algorithm-performance benchmark.'
save_json(qa/'SUMMARY.json',summary)
print(json.dumps(summary,indent=2))
