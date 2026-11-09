# PolyPanGF 0.1.0 验证记录

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
