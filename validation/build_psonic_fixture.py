"""Synthetic mixed-ploidy interface test for native whole-genome pSONIC integration."""
from pathlib import Path
import argparse,json,random
from polypangf.core import write_tsv
from polypangf.inputs import write_fasta
parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);parser.add_argument('--script',type=Path,required=True)
args=parser.parse_args();d=args.output.resolve();d.mkdir(parents=True,exist_ok=True)
rng=random.Random(71)
proteins=['M'+''.join(rng.choice('ACDEFGHIKLMNPQRSTVWY') for _ in range(119+13*group)) for group in range(12)]
samples=[];candidates=[]
for dataset,ploidy,copies in [('DIP',2,1),('TET',4,2),('OCT',8,4)]:
    records=[];gff=[]
    for copy in range(copies):
        for group,seq in enumerate(proteins):
            uid=f'{dataset}_h{copy}_g{group:02d}';start=1000+group*1000;end=start+599
            records.append((uid,seq))
            gff += [f'chr{copy+1}\tfixture\tgene\t{start}\t{end}\t.\t+\t.\tID={uid}',
                    f'chr{copy+1}\tfixture\tmRNA\t{start}\t{end}\t.\t+\t.\tID={uid}.t;Parent={uid};protein_id={uid}']
            if group<6:candidates.append({'dataset':dataset,'protein_id':uid,'architecture':'synthetic'})
    write_fasta(d/(dataset+'.faa'),records)
    (d/(dataset+'.gff3')).write_text('##gff-version 3\n'+'\n'.join(gff)+'\n')
    samples.append({'dataset':dataset,'species':'Synthetic_'+dataset,'ploidy':ploidy,'represented_sets':copies,'relative_ploidy':copies,
                    'proteins':dataset+'.faa','gff':dataset+'.gff3'})
write_tsv(d/'datasets.tsv',samples,['dataset','species','ploidy','represented_sets','relative_ploidy','proteins','gff'])
write_tsv(d/'candidates.tsv',candidates,['dataset','protein_id','architecture'])
cfg={'schema_version':1,'family':'Synthetic native pSONIC interface validation','manifest':'datasets.tsv','output':'results','threads':2,'seed':42,
     'identification':{'candidates':'candidates.tsv'},'orthology':{'backend':'psonic','script':str(args.script.resolve())},'pan':{'unit':'species'}}
(d/'config.json').write_text(json.dumps(cfg,indent=2)+'\n')
print(d)
