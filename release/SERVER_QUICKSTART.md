# 服务器使用

源码目录为 /run/media/dd/存储盘/work/PolyPanGF 。

```bash
source /run/media/dd/存储盘/work/PolyPanGF/.venv/bin/activate
polypangf --help
polypangf doctor -c /run/media/dd/存储盘/work/PolyPanGF/demo/config.json
polypangf run -c /run/media/dd/存储盘/work/PolyPanGF/demo/config.json --resume
polypangf status -c /run/media/dd/存储盘/work/PolyPanGF/validation/strawberry_replay/config.json --verify
```

本机专用虚拟环境复用已安装的科学计算依赖，外部程序通过专用bin目录提供。分发到其他机器时按README及environment.yml安装依赖。
