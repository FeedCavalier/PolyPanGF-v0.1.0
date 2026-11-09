"""Real full-proteome HMM scan plus native family-tool smoke test."""
from pathlib import Path
import argparse,json
import re,gzip
from Bio import SeqIO
from polypangf.core import write_tsv
from polypangf.inputs import fasta,write_fasta
parser=argparse.ArgumentParser();parser.add_argument('--legacy-root',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args();p=args.legacy_root;d=args.output.resolve();d.mkdir(parents=True,exist_ok=True)
full=p/'inputs/F017'
reference=p/'pgfat_nc85_20261107/inputs/ALL1545_MADS.faa'
seqs=fasta(reference,protein=True)
write_fasta(d/'F006.faa',[(uid,seq) for uid,seq in seqs.items() if uid.startswith('F006_')])
with gzip.open(full/'Fragaria_vesca_v4.0.a2.proteins.fa.gz','rt') as handle:
    protein_ids=[record.id for record in SeqIO.parse(handle,'fasta')]
write_tsv(d/'F017.models.tsv',[{'protein_id':uid,'gene_id':re.sub(r'\.t\d+$','',uid),'transcript_id':uid,'cds_id':uid} for uid in protein_ids],['protein_id','gene_id','transcript_id','cds_id'])
samples=[{'dataset':'F017','species':'Fragaria_vesca','ploidy':2,'represented_sets':1,'proteins':str(full/'Fragaria_vesca_v4.0.a2.proteins.fa.gz'),'gff':str(full/'Fragaria_vesca_v4.0.a2.genes.gff3.gz'),'cds':str(full/'Fragaria_vesca_v4.0.a2.cds.fa.gz')},
         {'dataset':'F006','species':'Fragaria_daltoniana','ploidy':2,'represented_sets':1,'proteins':'F006.faa','gff':'','cds':''}]
samples[0]['gff']='';samples[0]['model_map']='F017.models.tsv';samples[1]['model_map']=''
write_tsv(d/'datasets.tsv',samples,['dataset','species','ploidy','represented_sets','proteins','gff','cds','model_map'])
cfg={'schema_version':1,'family':'MADS native tool validation','manifest':'datasets.tsv','output':'results','threads':4,'seed':42,
     'preparation':{'invalid_protein':'mask'},
     'identification':{'profiles':[{'label':'MADS','path':str(p/'pgfat_nc85_20261107/structure_domains/PF00319.hmm'),'primary':True},{'label':'K_box','path':str(p/'pgfat_nc85_20261107/structure_domains/PF01486.hmm'),'primary':False}],'architecture_rules':[{'label':'MIKC_like','requires':['MADS','K_box']},{'label':'MADS_domain_only','requires':['MADS']}]},
     'orthology':{'backend':'orthofinder'},'pan':{'unit':'species','softcore_fraction':.9},'tree':{},'motifs':{'number':2,'min_width':6,'max_width':12}}
(d/'config.json').write_text(json.dumps(cfg,indent=2)+'\n');print(d)
