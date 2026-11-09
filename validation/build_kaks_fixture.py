"""Use three previously audited real one-to-one pairs to verify the native codon adapter."""
import argparse,json
from pathlib import Path
from polypangf.core import read_tsv,write_tsv
from polypangf.inputs import fasta,write_fasta
parser=argparse.ArgumentParser();parser.add_argument('--legacy-root',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args();p=args.legacy_root;d=args.output.resolve();d.mkdir(parents=True,exist_ok=True)
table=read_tsv(p/'pgfat_nc85_20261107/kaks/ALL_PAIRWISE_KAKS.tsv')
pairs=[row for row in table if row['status']=='PASS_FILTERED_YN00' and {row['query1'][:4],row['query2'][:4]}=={'F006','F017'}][:3]
assert len(pairs)==3
proteins=fasta(p/'pgfat_nc85_20261107/inputs/ALL1545_MADS.faa',protein=True)
samples=[];supplied=[];groups=[]
for dataset,species in [('F006','Fragaria_daltoniana'),('F017','Fragaria_vesca')]:
    source=p/'pgfat_nc85_20261107/validated_cds'/dataset/'VALIDATED_CDS.fna'
    if not source.exists():source=p/'pgfat_integration/extended_analysis_20261104/validated_cds'/dataset/'VALIDATED_CDS.fna'
    coding=fasta(source)
    ids=[next(uid for uid in [row['query1'],row['query2']] if uid.startswith(dataset)) for row in pairs]
    write_fasta(d/(dataset+'.faa'),[(uid,proteins[uid]) for uid in ids]);write_fasta(d/(dataset+'.cds.fa'),[(uid,coding[uid]) for uid in ids])
    samples.append({'dataset':dataset,'species':species,'ploidy':2,'represented_sets':1,'proteins':dataset+'.faa','cds':dataset+'.cds.fa'})
    supplied += [{'dataset':dataset,'protein_id':uid} for uid in ids]
    groups += [{'group':row['orthogroup'],'dataset':dataset,'gene_id':uid} for row,uid in zip(pairs,ids)]
write_tsv(d/'datasets.tsv',samples,['dataset','species','ploidy','represented_sets','proteins','cds'])
write_tsv(d/'candidates.tsv',supplied,['dataset','protein_id']);write_tsv(d/'groups.tsv',groups,['group','dataset','gene_id'])
(d/'expected_pairs.json').write_text(json.dumps(pairs,indent=2)+'\n')
cfg={'schema_version':1,'family':'Real KaKs adapter validation','manifest':'datasets.tsv','output':'results','threads':2,'seed':42,
     'identification':{'candidates':'candidates.tsv'},'orthology':{'backend':'import','groups':'groups.tsv'},'kaks':{'anchor':'F017'}}
(d/'config.json').write_text(json.dumps(cfg,indent=2)+'\n');print(d)
