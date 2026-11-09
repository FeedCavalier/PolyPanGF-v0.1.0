from pathlib import Path
import json
import os
import subprocess
import sys
import tempfile
import zipfile

root=Path(__file__).resolve().parents[1]
wheel=next((root/'release/dist').glob('*.whl'))
with zipfile.ZipFile(wheel) as archive:
    for path in (root/'src/polypangf').glob('*.py'):
        assert archive.read('polypangf/'+path.name)==path.read_bytes()
with tempfile.TemporaryDirectory(prefix='polypangf_wheel_') as temporary:
    directory=Path(temporary);target=directory/'installed'
    subprocess.run([sys.executable,'-m','pip','install','--no-deps','--target',str(target),str(wheel)],check=True)
    environment=os.environ.copy();environment['PYTHONPATH']=str(target)
    environment['PATH']=str(root/'.venv/bin')+os.pathsep+environment['PATH']
    location=subprocess.check_output([sys.executable,'-c','import polypangf;print(polypangf.__file__)'],env=environment,text=True).strip()
    assert Path(location).is_relative_to(target)
    subprocess.run([sys.executable,'-m','polypangf','init',str(directory/'demo')],env=environment,check=True,stdout=subprocess.DEVNULL)
    subprocess.run([sys.executable,'-m','polypangf','run','-c',str(directory/'demo/config.json')],env=environment,check=True,stdout=subprocess.DEVNULL)
    assert (directory/'demo/results/report/report.html').is_file()
summary_path=root/'validation/release_qa/SUMMARY.json'
summary=json.loads(summary_path.read_text());summary['wheel_install_verified']=True
summary_path.write_text(json.dumps(summary,indent=2)+'\n')
print('Wheel installation, source identity and end-to-end example verified.')
