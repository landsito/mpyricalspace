"""Run directly with python -I -S -B; no package import or pytest conftest."""
import ast
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
import zipfile

ROOT = Path(__file__).resolve().parents[2]
HELPER = ROOT / 'mpyricalspace/vscode_ext/probe.py'
spec = importlib.util.spec_from_file_location('observer_probe', HELPER)
probe = importlib.util.module_from_spec(spec)
spec.loader.exec_module(probe)


class ProbeTests(unittest.TestCase):
    def test_registry_matches_models_without_import(self):
        tree = ast.parse((ROOT / 'mpyricalspace/_build.py').read_text())
        models = next(ast.literal_eval(n.value) for n in tree.body if isinstance(n, ast.Assign)
                      and any(isinstance(t, ast.Name) and t.id == 'EXTENSIONS' for t in n.targets))
        registry = json.loads(HELPER.with_name('modules.json').read_text())
        self.assertEqual(registry, {k: {'label': v[0], 'subdir': v[1]} for k, v in models.items()})

    def plan(self, root, build='cp-test'):
        directory = root / 'build' / build
        filename = directory / ('rocsatf' + probe.importlib.machinery.EXTENSION_SUFFIXES[0])
        plan = directory / 'meson-info/intro-install_plan.json'
        plan.parent.mkdir(parents=True)
        plan.write_text(json.dumps({'targets': {str(filename): {'destination': '{py_platlib}/mpyricalspace/' + filename.name}}}))
        return filename

    def test_editable_discovery_reads_plan_without_build(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td); expected = self.plan(root)
            with patch.object(probe, 'dependency_paths', side_effect=AssertionError('must not search installed package')):
                self.assertEqual(probe.discover(root)['rocsatf']['path'], str(expected.resolve()))

    def test_multiple_builds_require_selection(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td); first = self.plan(root, 'first'); self.plan(root, 'second')
            self.assertIsNone(probe.discover(root)['rocsatf']['path'])
            self.assertEqual(probe.discover(root, str(first.parent))['rocsatf']['path'], str(first.resolve()))

    def test_normal_install_discovery_never_executes_init(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td); package = root / 'site-packages/mpyricalspace'; package.mkdir(parents=True)
            (package / '__init__.py').write_text('raise RuntimeError("DO NOT IMPORT")')
            artifact = package / ('rocsatf' + probe.importlib.machinery.EXTENSION_SUFFIXES[0]); artifact.touch()
            with patch.object(probe, 'dependency_paths', return_value=[str(package.parent)]):
                self.assertEqual(probe.discover(root)['rocsatf']['path'], str(artifact.resolve()))

    def run_check(self, name, path):
        result = subprocess.run([sys.executable, '-I', '-S', '-B', str(HELPER), 'check', name, str(path)],
                                capture_output=True, text=True, timeout=15)
        self.assertEqual(result.returncode, 0, result.stderr)
        return json.loads(result.stdout.rsplit('MPYRICALSPACE_PROBE=', 1)[1])

    def test_missing_and_broken_artifacts(self):
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / ('rocsatf' + probe.importlib.machinery.EXTENSION_SUFFIXES[0])
            self.assertEqual(self.run_check('rocsatf', path)['state'], 'unknown')
            path.write_bytes(b'not a native module')
            self.assertEqual(self.run_check('rocsatf', path)['state'], 'failed')
            self.assertEqual(self.run_check('other', path)['state'], 'unknown')

    def test_startup_hooks_and_pythonpath_are_not_executed(self):
        import os
        with tempfile.TemporaryDirectory() as td:
            root = Path(td); marker = root / 'executed'
            (root / 'sitecustomize.py').write_text(f'open({str(marker)!r}, "w").write("bad")')
            env = dict(os.environ, PYTHONPATH=str(root))
            r = subprocess.run([sys.executable, '-I', '-S', '-B', str(HELPER), 'discover', '--root', td],
                               capture_output=True, text=True, env=env, timeout=15)
            self.assertEqual(r.returncode, 0); self.assertFalse(marker.exists())

    def test_packaged_extension_contains_standalone_helpers(self):
        spec = importlib.util.spec_from_file_location('standalone_packager', ROOT / 'mpyricalspace/_vscode.py')
        packager = importlib.util.module_from_spec(spec); spec.loader.exec_module(packager)
        with tempfile.TemporaryDirectory() as td:
            vsix = packager.build_vsix(td)
            with zipfile.ZipFile(vsix) as archive:
                for file in ('observer.js', 'probe.py', 'modules.json', 'extension.js'):
                    self.assertIn('extension/' + file, archive.namelist())


if __name__ == '__main__':
    unittest.main()
