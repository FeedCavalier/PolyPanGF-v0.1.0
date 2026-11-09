"""Build a source-and-wheel release with portable examples and checksums."""
from pathlib import Path
import hashlib
import importlib.metadata
import json
import platform
import shutil
import subprocess
import sys
import zipfile
from PIL import Image, ImageDraw

from polypangf.core import save_json, sha256

root=Path(__file__).resolve().parents[1]
release=root/'release';release.mkdir(exist_ok=True)
if not list((release/'dist').glob('*.whl')):
    subprocess.run([sys.executable,'-m','pip','wheel','--no-deps','--no-build-isolation','.', '-w',str(release/'dist')],cwd=root,check=True)
summary=json.loads((root/'validation/release_qa/SUMMARY.json').read_text())
versions={name:importlib.metadata.version(name) for name in ['polypangf','biopython','numpy','scipy','matplotlib','PyYAML','igraph']}
versions.update(python=platform.python_version(),platform=platform.platform())
save_json(release/'RUNTIME_VERSIONS.json',versions)
text='''# PolyPanGF 0.1.0 验证记录

本版完成15项回归测试、22个子命令的帮助检查、无写入运行计划，以及六种端到端场景的运行、重复运行和输出散列核验。构建后的wheel已独立安装并运行完整示例，包内源码与审核版本逐文件一致。

## 真实数据回放

13个草莓物种的1545个冻结候选保留完整，1494个成员进入122群，51个未分组候选独立保留。覆盖度类别为5个Core、10个Softcore、102个Shell和5个Specific。旧数据中共用占位ID的未定位蛋白通过各自唯一ID导入，并保留蛋白单位状态。

同版本GO背景包含18662个已注释单位，85个森林草莓候选均匹配，恢复5项全家族显著条目及14项类别比较。145条单核空间marker记录全部匹配，涉及28个基因。

## 原生外部程序

完整森林草莓蛋白组采用显式转录本映射，HMMER识别85个候选，与另一个草莓物种的57个成员一起运行家族OrthoFinder、MAFFT、FastTree和MEME/FIMO。原始蛋白组中8条含非法残基的记录按显式mask策略处理并登记，坐标冲突保留为输入错误，测试通过模型映射开展序列分析。

原生pSONIC使用固定官方提交bcebe9af8d6bf3cfee0779e6e891a82b10ffe486。完整自动准备和预计算输入两条路径均通过合成用例，二、四、八倍体对应1、2、4套代表，目标家族42个成员恢复为6个Core群。自动准备使用84个全基因组代表模型。当前MCScanX发行版不导出串联边，适配层按相邻基因及BLAST E≤1×10⁻⁵生成对应输入，来源写入状态文件。

三对已审核的真实草莓基因通过MAFFT、PAL2NAL和YN00计算，Ka、Ks及omega与原有估计在0.00011容差内一致。

## 当前范围

上述验证检查功能、数据单位、数值与重试行为。pSONIC真实多物种全基因组的性能基准，以及跨家族方法比较仍需单独开展。复制类型通过标准表导入。表达模块提供未分层群成员置换，项目中的架构分层置换和表达簇关联尚未通用化。完整单细胞与空间矩阵聚类由相应专门流程执行。

原生pSONIC脚本为单独依赖，发布包不捆绑其源码。软件与数据各自遵循来源许可，项目自身发布许可证尚未指定。

## 运行示例

```bash
python -m pip install dist/polypangf-0.1.0-py3-none-any.whl
polypangf run -c examples/demo/config.json
polypangf run -c examples/demo/config.json --resume
polypangf run -c examples/strawberry_replay/config.json
polypangf status -c examples/strawberry_replay/config.json --verify
```

examples中的配置和输入使用相对路径。原生程序测试脚本保留在validation目录，需要自行提供已安装的外部程序及相应原始参考数据。
'''
(release/'VALIDATION.md').write_text(text,encoding='utf-8')
(release/'SERVER_QUICKSTART.md').write_text('''# 服务器使用

源码目录为 /run/media/dd/存储盘/work/PolyPanGF 。

```bash
source /run/media/dd/存储盘/work/PolyPanGF/.venv/bin/activate
polypangf --help
polypangf doctor -c /run/media/dd/存储盘/work/PolyPanGF/demo/config.json
polypangf run -c /run/media/dd/存储盘/work/PolyPanGF/demo/config.json --resume
polypangf status -c /run/media/dd/存储盘/work/PolyPanGF/validation/strawberry_replay/config.json --verify
```

本机专用虚拟环境复用已安装的科学计算依赖，外部程序通过专用bin目录提供。分发到其他机器时按README及environment.yml安装依赖。
''',encoding='utf-8')
members={}
def add(path,name):members[str(name)]=Path(path)
for filename in ['pyproject.toml','README.md','CHANGELOG.md','environment.yml','.gitignore']:
    add(root/filename,filename)
for folder in ['src/polypangf','tests']:
    for path in (root/folder).glob('*.py'):add(path,path.relative_to(root))
for filename in ['build_strawberry_replay.py','build_native_family_test.py','build_psonic_fixture.py','build_kaks_fixture.py','extend_replay_annotations.py','release_check.py','package_release.py','smoke_wheel.py']:
    add(root/'validation'/filename,Path('validation')/filename)
for case,source in [('demo',root/'demo'),('strawberry_replay',root/'validation/strawberry_replay')]:
    for path in source.iterdir():
        if path.is_file() and path.suffix in ['.faa','.tsv','.json','.txt'] and not path.name.startswith('run.'):
            add(path,Path('examples')/case/path.name)
for path in (release/'dist').glob('*.whl'):add(path,Path('dist')/path.name)
for name in ['VALIDATION.md','SERVER_QUICKSTART.md','RUNTIME_VERSIONS.json']:
    add(release/name,name)
add(root/'validation/release_qa/SUMMARY.json','validation/RESULTS.json')
for path in (root/'validation/release_qa').glob('unit_tests.*'):add(path,Path('validation')/path.name)
for case in ['demo','strawberry_replay','native_family','psonic_native','kaks_real']:
    source=root/'demo/results' if case=='demo' else root/'validation'/case/'results'
    for path in source.glob('*/summary.json'):add(path,Path('validation/summaries')/case/path.parent.name/path.name)
for key in ['strawberry_replay_status.stdout','psonic_native_status.stdout','psonic_prepared_status.stdout','native_family_status.stdout','kaks_real_status.stdout']:
    add(root/'validation/release_qa'/key,Path('validation')/key)
manifest=''.join(sha256(path)+'  '+name+'\n' for name,path in sorted(members.items()))
(release/'SHA256SUMS.txt').write_text(manifest)
add(release/'SHA256SUMS.txt','SHA256SUMS.txt')
archive=release/'PolyPanGF-0.1.0.zip'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as output:
    for name,path in sorted(members.items()):output.write(path,'PolyPanGF-0.1.0/'+name)
with zipfile.ZipFile(archive) as output:
    assert output.testzip() is None
    for line in manifest.splitlines():
        digest,name=line.split('  ',1)
        assert hashlib.sha256(output.read('PolyPanGF-0.1.0/'+name)).hexdigest()==digest
preview_paths=list((root/'demo/results/plot').glob('*.png'))+list((root/'validation/strawberry_replay/results/plot').glob('*.png'))
sheet=Image.new('RGB',(1800,450*((len(preview_paths)+2)//3)),'#eeeeee')
for i,path in enumerate(preview_paths):
    im=Image.open(path).convert('RGB');im.thumbnail((585,415));x=(i%3)*600;y=(i//3)*450
    sheet.paste(im,(x,y+25));ImageDraw.Draw(sheet).text((x+5,y+5),path.stem,fill='black')
sheet.save(release/'plot_review.jpg')
save_json(release/'RELEASE_MANIFEST.json',{'archive':archive.name,'sha256':sha256(archive),'files':len(members),'size_bytes':archive.stat().st_size,'validation_passed':True})
print(json.dumps({'archive':str(archive),'files':len(members),'sha256':sha256(archive)},indent=2))
