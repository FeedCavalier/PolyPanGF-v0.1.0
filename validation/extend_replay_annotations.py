"""Add version-matched real GO and author marker inputs to the 13-species replay."""
import argparse,json
from pathlib import Path
from polypangf.core import read_tsv,write_tsv
parser=argparse.ArgumentParser();parser.add_argument('--legacy-root',type=Path,required=True);parser.add_argument('--project',type=Path,required=True)
args=parser.parse_args();p=args.legacy_root;d=args.project.resolve()
frozen=read_tsv(p/'pgfat_nc85_20261107/PGFAT_reference_release_20261108/analysis_evidence/FROZEN_GENE_EVIDENCE.tsv')
names={row['query_id']:row['locus_id'] for row in frozen if row['dataset_id']=='F017'}
terms=read_tsv(p/'bioconductor_GO_KEGG_20261109/GO_TERM2GENE.tsv')
labels={row['term']:row['name'] for row in read_tsv(p/'bioconductor_GO_KEGG_20261109/GO_TERM2NAME.tsv')}
write_tsv(d/'GO.tsv',[{'gene_id':names.get(row['gene'],row['gene']),'term':row['term'],'description':labels.get(row['term'],row['term'])} for row in terms],['gene_id','term','description'])
markers=read_tsv(p/'NC_advanced_GO_20261108/release/source_data/NC_ALL145_MARKER_RECORDS_WITH_FLAGS.tsv')
write_tsv(d/'markers.tsv',[{'dataset':'F017','gene_id':row['v4_gene_id'],'assay':row['assay'],'stage':row['stage'],'cluster':row['cluster'],'padj':row['p_val_adj']} for row in markers],['dataset','gene_id','assay','stage','cluster','padj'])
cfg=json.loads((d/'config.json').read_text());cfg['enrichment']={'dataset':'F017','database':'GO','term2gene':'GO.tsv','min_term_size':2,'fdr':.05};cfg['markers']={'table':'markers.tsv'}
(d/'config.json').write_text(json.dumps(cfg,indent=2)+'\n');print(d)
