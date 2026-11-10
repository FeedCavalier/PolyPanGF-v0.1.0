# PolyPanGF

面向混合倍性基因组的泛基因家族分析工具包，当前版本为 0.1.0。

输入蛋白组、材料清单及家族识别模型，完成候选鉴定、同源群整理、泛家族统计、可选下游分析和标准图表。物种存在记录、材料拷贝数、生物学倍性及组装代表套数分别保存。<img width="2161" height="2881" alt="image" src="https://github.com/user-attachments/assets/a247fe92-1203-4169-970c-55d37cea1a72" />

## 安装和快速运行

支持 Linux 和 Python 3.9 以上版本。所有实际分析在运行工具包的主机执行。

```bash
python -m pip install .
polypangf --version
polypangf init demo
polypangf validate -c demo/config.json
polypangf doctor -c demo/config.json
polypangf run -c demo/config.json
polypangf run -c demo/config.json --resume
polypangf status -c demo/config.json
polypangf export -c demo/config.json -o demo_results.zip
```

`init` 生成明确标注的合成示例，用于检查安装、文件格式和运行入口。示例结果不用于生物学验证。

## 命令

| 命令 | 功能 |
|---|---|
| `init` | 创建示例项目和配置 |
| `validate` | 校验清单、FASTA、GFF3和配置输入 |
| `doctor` | 检查Python依赖与所选模块所需的外部程序 |
| `prepare` | 标准化ID、建立基因模型对应并核验CDS翻译 |
| `identify` | HMMER鉴定或导入候选，在有结构域支持的异构体中选择代表 |
| `synteny` | 全蛋白组OrthoFinder与MCScanX，准备原生pSONIC输入 |
| `orthology` | 导入分组、运行家族OrthoFinder或调用官方pSONIC |
| `pan` | 覆盖度分类、存在矩阵、拷贝数及精确期望累积曲线 |
| `properties` | 蛋白长度、分子量、理论pI与GRAVY |
| `tree` | MAFFT与FastTree家族树 |
| `motifs` | 去重蛋白上的MEME基序识别及全部候选的FIMO扫描 |
| `expression` | 时期表达、Tau、群内多拷贝相关与置换检验 |
| `enrichment` | 自定义GO、KEGG或其他注释的Fisher富集与BH校正 |
| `te` | TE区间合并、方向化覆盖曲线及同染色体长度匹配背景 |
| `kaks` | MAFFT、PAL2NAL与PAML YN00单拷贝配对分析 |
| `duplication` | 导入DupGen_finder等工具提供的复制类型 |
| `markers` | 导入单细胞或空间marker清单，匹配成员并统计覆盖广度 |
| `plot` | PDF、SVG和300 dpi PNG标准图表 |
| `report` | 无外部网络资源的HTML报告及Markdown摘要 |
| `run` | 按配置运行工作流 |
| `status` | 完成状态及失败目录 |
| `export` | 带SHA256清单的ZIP结果包 |

执行 `polypangf COMMAND --help` 查看选项。分析命令支持以下公共参数。

```text
-c, --config       JSON或YAML配置文件
-o, --output       覆盖输出目录，相对路径以配置文件所在目录为准
-t, --threads      外部工具线程数
--seed             随机种子
--resume           校验输入、代码及输出散列后复用已完成步骤
--force            重算目标步骤及其依赖
--dry-run          只显示依赖计划，不创建分析目录
-v, --verbose      调试日志
-q, --quiet        只显示错误及最终JSON结果
```

`run --stages identify,pan` 指定目标阶段，依赖自动补齐。`--resume` 在输入或输出改变后重算对应步骤。`validate --quick` 跳过完整序列扫描，`status --verify` 验证输出散列。每阶段使用独立工作目录，失败目录及日志保留。已有完整结果在新步骤失败时继续保留。项目级文件锁防止并发写入。外部程序通过参数列表调用，路径可以包含空格及中文。单材料可独立运行鉴定和蛋白属性等模块，泛家族比较至少需要两个独立统计单位。

## 材料清单

清单为制表符分隔的TSV，路径相对清单文件。必需列如下。

| 列 | 含义 |
|---|---|
| `dataset` | 唯一材料或组装代码，字母开头，使用字母、数字、下划线和连字符 |
| `species` | 物种标识，同种不同组装填写相同值 |
| `ploidy` | 生物学倍性，例如2、4、8 |
| `represented_sets` | 输入所代表的基本染色体组或分相套数 |
| `proteins` | 蛋白FASTA，支持gzip |

可选列为 `gff`、`cds`、`model_map`、`role` 和 `relative_ploidy`。`role` 可取 `analysis` 或 `outgroup`，外群参与同源关系推断，泛家族覆盖统计仅使用analysis材料。自动全基因组pSONIC流程要求显式提供 `relative_ploidy`，其数值按被比较谱系中的相对拷贝期望设定。

`model_map` 必须包含 `protein_id`、`gene_id`，可附 `transcript_id` 和 `cds_id`。推荐提供与所用注释版本对应的映射。没有GFF3或映射时按蛋白ID登记单位，并明确标记 `protein_id_assumption`。未定位模型应使用各自唯一ID并填写 `unit_type: protein`，避免多个成员共用占位基因ID。不同材料中相同的原始ID通过内部命名空间隔离。

GFF3使用1-based闭区间。TE表使用0-based半开区间。CDS仅在长度为3的倍数且翻译与相应蛋白一致时进入Ka/Ks模块。

非标准蛋白残基默认报错。可显式设置 `preparation.invalid_protein` 为 `mask`，将非法残基置为X，或设为 `skip` 排除对应记录。处理记录写入 `prepare/sequence_audit.tsv`。同一FASTA中的重复ID、冲突坐标及共享占位基因ID会报错，便于在分析前修正来源数据。

## 配置示例

```yaml
schema_version: 1
family: MADS
manifest: datasets.tsv
output: results
threads: 8
seed: 42
identification:
  profiles:
    - label: MADS
      path: profiles/PF00319.hmm
      primary: true
    - label: K_box
      path: profiles/PF01486.hmm
      primary: false
  domain_evalue: 1.0e-5
  min_domain_coverage: 0.0
  architecture_rules:
    - label: MIKC_like
      requires: [MADS, K_box]
    - label: MADS_domain_only
      requires: [MADS]
orthology:
  backend: orthofinder
pan:
  unit: species
  softcore_fraction: 0.9
```

架构规则由用户定义并按顺序匹配。简单结构域组合提供架构标签，具体亚家族分类应结合专业分类工具或参考序列。已有候选可以用 `identification.candidates: candidates.tsv` 替代HMM扫描，表头为 `dataset`、`protein_id`，可选 `architecture`。

### 分组后端

1. `backend: import` 配合 `groups`，接受长表 `group,dataset,gene_id` 或OrthoFinder宽表。导入未知成员默认报错，全基因组结果可显式设置 `allow_nonfamily: true`。重复成员及多重归群报错。没有进入任何群的候选保留为Unassigned。
2. `backend: orthofinder` 运行家族候选的OrthoFinder分组，保留家族分析范围。
3. `backend: psonic` 配合官方 `script` 路径，可自动运行全蛋白组OrthoFinder与MCScanX。该模式要求所有全基因组代表模型可定位，并显式声明相对倍性。官方脚本及其igraph依赖需单独安装，项目不捆绑第三方源码。

```yaml
orthology:
  backend: psonic
  script: software/pSONIC.py
  python: /path/to/python-with-igraph
```

也可直接使用已准备好的全基因组输入，增加以下字段即可跳过自动synteny阶段。

```yaml
  prefix: genome_inputs/project
  groups: genome_inputs/Orthogroups.tsv
  sequence_ids: genome_inputs/SequenceIDs.txt
  species_ids: genome_inputs/SpeciesIDs.txt
  ploidy: genome_inputs/relative_ploidy.tsv
```

`prefix` 对应 `.gff`、`.tandem`、`.collinearity` 文件，基因编码须与OrthoFinder一致。pSONIC相对倍性文件为一行数字，按OrthoFinder物种列顺序排列。工具会为官方无表头的输出恢复物种列名，再提取目标家族。调用原生pSONIC和导入已有分组会记录不同的来源。

部分MCScanX发行版不输出串联边文件。自动准备流程在这种情况下根据同染色体相邻基因的BLAST命中生成串联边，阈值为E≤1×10⁻⁵，并在状态文件登记 `tandem_source`。已有外部串联边文件可通过预计算输入模式完整保留。

### 统计单位

`pan.unit` 可选 `species` 或 `dataset`。同种多个组装在species覆盖率中贡献一次存在记录，拷贝数始终按dataset分别输出。Core为全部单位检出，Specific为单个单位检出，Softcore为其余达到指定覆盖阈值的群，Shell为其他多单位群。两个单位的面板中不会出现Shell。未分组候选不并入Specific。

`copies_per_represented_set` 为观察拷贝数除以声明的代表套数，是描述性指标。其解释取决于组装与注释完整性。数据缺失、基因丢失及祖先拷贝身份需要更多证据判定。累积曲线使用所有等大小子集的精确期望，独立保存统计源表。

### 可选分析

```yaml
tree: {}
motifs:
  number: 10
  min_width: 6
  max_width: 50
  pvalue: 0.0001
expression:
  dataset: OCT
  matrix: expression.tsv
  samples: samples.tsv
  permutations: 19999
  min_max_expression: 1
enrichment:
  dataset: DIP
  database: GO
  term2gene: genome_GO.tsv
  fdr: 0.05
  min_term_size: 2
te:
  intervals: TE_intervals.tsv
  flank: 2000
  bins_per_region: 30
kaks:
  anchor: DIP
  min_codons: 50
  min_coverage: 0.7
  min_ks: 0.01
  max_ks: 2
  max_omega: 99
  timeout: 300
duplication:
  table: duplication_types.tsv
  source_method: DupGen_finder
markers:
  table: marker_records.tsv
  fdr: 0.05
```

- 表达矩阵为 `gene_id` 加样本数值列。元数据为 `sample,stage`，时期顺序沿用元数据的首次出现顺序。多重复按时期取均值。缺失行单独登记，原始数值要求非负且有限。群内相关基于log(1+value)，各可评价多拷贝群等权，单侧置换使用加一校正。本版提供未分层置换，架构分层置换和表达聚类关联仍沿用项目专用脚本。
- 功能注释为 `gene_id,term`，可加 `description`。背景为该版本注释中全部基因，可另用 `background` 指定含gene_id列的表。`valid_terms` 指定含term列的物种适用通路清单。全家族与覆盖类别分开校正，类别内所有条目及零命中项联合BH校正，仅检验具有注释成员的类别。Unassigned单独列出。本版每次配置一个注释体系，GO和KEGG可分别运行至独立输出目录。
- TE表为 `dataset,chrom,start,end`，要求预先筛选为TE区间。同染色体上长度最接近的非候选基因为背景，允许背景复用并在原表保留锚点关系。缺少TE染色体记录的候选不当作零覆盖。`chrom_sizes` 可提供 `dataset,chrom,length` 表以裁剪染色体边界。未提供长度时，上游负坐标截断到0，末端范围按输入坐标计算。
- 复制类型表为 `dataset,gene_id,duplication_class`。类型可取 `wgd_segmental,tandem,proximal,transposed,dispersed,unsupported`。本版导入已有分类结果，未提供记录的候选标记not_evaluated。
- marker表为 `dataset,gene_id,assay,stage,cluster,padj`。完整匹配行、显著清单广度及未匹配记录分别保存。本版使用处理结果和清单，完整单细胞与空间矩阵的聚类分析需另行运行对应流程。

## 输出与追溯

每阶段一个目录，包含 `state.json`、表格及外部命令日志。状态记录配置、代码、输入和输出SHA256及运行耗时。错误保存在 `failures`。`plot` 保存PDF、SVG和PNG，`report/report.html` 可直接在浏览器打开。报告旁的Markdown和JSON便于后续文稿与工作流整合。

配置与材料路径建议固定保存。结果导出包含标准化序列和结果，外部原始输入仍通过状态文件中的路径及SHA256追溯。

## 开发和测试

```bash
python -m pip install -e .
python -m unittest discover -s tests -v
polypangf run -c demo/config.json --dry-run
```

本版重点提供运行入口、数据单位约束及可复用分析模块。与pSONIC等方法的性能比较需要独立基准测试。项目发布许可证尚未指定，第三方软件及数据各自遵循原有许可证。

## 方法来源

- pSONIC https://github.com/conJUSTover/pSONIC ，DOI 10.1093/g3journal/jkab170
- OrthoFinder https://github.com/OrthoFinder/OrthoFinder ，DOI 10.1186/s13059-019-1832-y
- MCScanX DOI 10.1093/nar/gkr1293
- HMMER DOI 10.1371/journal.pcbi.1002195
- MAFFT DOI 10.1093/molbev/mst010
- FastTree DOI 10.1371/journal.pone.0009490
- PAL2NAL DOI 10.1093/nar/gkl315
- YN00 DOI 10.1093/oxfordjournals.molbev.a026236

