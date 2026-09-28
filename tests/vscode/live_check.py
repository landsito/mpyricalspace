"""Opt-in local acceptance check: load each real artifact, never run a build."""
import json
from pathlib import Path
import subprocess
import sys

root = Path(__file__).resolve().parents[2]
helper = root / 'mpyricalspace/vscode_ext/probe.py'


def snapshot():
    return {str(p): (p.stat().st_size, p.stat().st_mtime_ns)
            for p in (root / 'build').rglob('*') if p.is_file()}


def run(*args):
    process = subprocess.run([sys.executable, '-I', '-S', '-B', str(helper), *args],
                             capture_output=True, text=True, timeout=20)
    if process.returncode:
        raise RuntimeError(process.stderr)
    return json.loads(process.stdout.rsplit('MPYRICALSPACE_PROBE=', 1)[1])


before = snapshot()
modules = run('discover', '--root', str(root))
results = {name: run('check', name, info['path']) if info['path'] else
           {'state': 'unknown', 'error': info['error']} for name, info in modules.items()}
changed = before != snapshot()
print(json.dumps({'modules': results, 'build_files_changed': changed}, indent=2))
assert not changed, 'Probe changed build files'
assert all(info['state'] == 'ok' for info in results.values()), 'Some wrappers did not load'
