"""Create a portable replay fixture from this project's frozen 13-species release."""
import argparse
from collections import defaultdict
from pathlib import Path
import json

from polypangf.core import read_tsv, write_tsv
from polypangf.inputs import fasta, write_fasta

parser = argparse.ArgumentParser()
parser.add_argument('--legacy-root', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()
root = args.legacy_root
dest = args.output.resolve();dest.mkdir(parents=True, exist_ok=True)
release = root / 'pgfat_nc85_20261107/PGFAT_reference_release_20261108'
rows = read_tsv(release / 'analysis_evidence/FROZEN_GENE_EVIDENCE.tsv')
master = {row['query_id']: row for row in read_tsv(root / 'pgfat_nc85_20261107/reference_complete_20261108/source_data/CURRENT_MASTER_WITH_SEPARATE_CULTIVATED_SCOPE.tsv')}
sequences = fasta(root / 'pgfat_nc85_20261107/inputs/ALL1545_MADS.faa', protein=True)
by_dataset = defaultdict(list)
for row in rows:by_dataset[row['dataset_id']].append(row)
def locus_id(row):
    return row['query_id'] if row['locus_id'] == 'UNMAPPED_PROTEIN_ONLY' else row['locus_id']
samples = [];supplied = [];groups = [];duplicate = []
for dataset, members in sorted(by_dataset.items()):
    write_fasta(dest / (dataset + '.faa'), [(row['query_id'], sequences[row['query_id']]) for row in members])
    write_tsv(dest / (dataset + '.models.tsv'), [{'protein_id': row['query_id'], 'gene_id': locus_id(row), 'unit_type': 'protein' if row['locus_id'] == 'UNMAPPED_PROTEIN_ONLY' else 'gene'} for row in members], ['protein_id', 'gene_id', 'unit_type'])
    first = members[0]
    samples.append({'dataset': dataset, 'species': master[first['query_id']]['species'], 'ploidy': first['biological_ploidy'],
                    'represented_sets': first['represented_complements'], 'proteins': dataset + '.faa', 'model_map': dataset + '.models.tsv'})
    for row in members:
        supplied.append({'dataset': dataset, 'protein_id': row['query_id'], 'architecture': master[row['query_id']]['itak_architecture']})
        if row['MSA_external_HOG'].startswith('N1.HOG'):
            groups.append({'group': row['MSA_external_HOG'], 'dataset': dataset, 'gene_id': locus_id(row)})
        kind = row['existing_duplication_class']
        duplicate.append({'dataset': dataset, 'gene_id': locus_id(row), 'duplication_class': 'wgd_segmental' if kind == 'wgd' else kind if kind in ['tandem', 'proximal', 'transposed', 'dispersed'] else 'unsupported'})
write_tsv(dest / 'datasets.tsv', samples, ['dataset', 'species', 'ploidy', 'represented_sets', 'proteins', 'model_map'])
write_tsv(dest / 'candidates.tsv', supplied, ['dataset', 'protein_id', 'architecture'])
write_tsv(dest / 'groups.tsv', groups, ['group', 'dataset', 'gene_id'])
write_tsv(dest / 'duplication.tsv', duplicate, ['dataset', 'gene_id', 'duplication_class'])
cfg = {'schema_version': 1, 'family': 'MADS frozen 13-species replay', 'manifest': 'datasets.tsv', 'output': 'results', 'threads': 4, 'seed': 42,
       'identification': {'candidates': 'candidates.tsv'}, 'orthology': {'backend': 'import', 'groups': 'groups.tsv'},
       'pan': {'unit': 'species', 'softcore_fraction': .9}, 'duplication': {'table': 'duplication.tsv', 'source_method': 'previously audited DupGen_finder results'}}
(dest / 'config.json').write_text(json.dumps(cfg, indent=2) + '\n')
(dest / 'expected.json').write_text(json.dumps({'candidates': 1545, 'assigned': 1494, 'unassigned': 51, 'groups': 122,
                                              'classes': {'Core': 5, 'Softcore': 10, 'Shell': 102, 'Specific': 5}}, indent=2) + '\n')
print(dest)
