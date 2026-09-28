"""Standalone observer. Run with Python -I -S -B; never import mpyricalspace."""
import argparse
import importlib.abc
import importlib.machinery
import importlib.util
import json
import os
from pathlib import Path
import sys
import sysconfig


def dependency_paths():
    """Expose dependencies without executing site, .pth files or startup hooks."""
    paths = [sysconfig.get_path(key) for key in ('platlib', 'purelib')]
    allow_user = True
    # Before Python 3.14, -S does not set sys.prefix to a venv's directory.
    executable = Path(sys.executable).absolute()
    for prefix in (executable.parent, executable.parent.parent):
        config = prefix / 'pyvenv.cfg'
        if config.is_file():
            venv_site = (prefix / 'Lib/site-packages' if os.name == 'nt' else
                         prefix / f'lib/python{sys.version_info.major}.{sys.version_info.minor}/site-packages')
            settings = dict(line.lower().split('=', 1) for line in config.read_text().splitlines() if '=' in line)
            system_site = any(key.strip() == 'include-system-site-packages' and value.strip() == 'true'
                              for key, value in settings.items())
            paths = [str(venv_site)] + (paths if system_site else [])
            allow_user = system_site
            break
    # Reading site helpers is safe under -S: site.main() is not run.
    if allow_user:
        import site
        paths.append(site.getusersitepackages())
    return list(dict.fromkeys(p for p in paths if p and Path(p).is_dir()))


def artifact_name(filename):
    for suffix in importlib.machinery.EXTENSION_SUFFIXES:
        if filename.endswith(suffix):
            name = filename[:-len(suffix)]
            if name.isidentifier():
                return name
    return None


def discover(root, build_dir='', package_dir=''):
    registry = json.loads(Path(__file__).with_name('modules.json').read_text())
    candidates = {name: set() for name in registry}
    root = Path(root).resolve()
    plans = ([Path(build_dir).resolve() / 'meson-info/intro-install_plan.json'] if build_dir else
             sorted((root / 'build').glob('**/meson-info/intro-install_plan.json')))
    if not package_dir:
        for plan in plans:
            data = json.loads(plan.read_text())
            for source, info in data.get('targets', {}).items():
                destination = info.get('destination', '').replace('\\', '/')
                if '/mpyricalspace/' not in destination:
                    continue
                name = artifact_name(Path(source).name)
                if name in candidates:
                    candidates[name].add(str(Path(source).resolve()))
    has_build = any(candidates.values()) or bool(build_dir)
    if package_dir or not has_build:
        directories = [Path(package_dir).resolve()] if package_dir else [Path(p) / 'mpyricalspace' for p in dependency_paths()]
        for directory in directories:
            if directory.is_dir():
                for source in directory.iterdir():
                    name = artifact_name(source.name)
                    if name in candidates:
                        candidates[name].add(str(source.resolve()))
    results = {}
    for name, info in registry.items():
        matches = sorted(candidates[name])
        results[name] = dict(info, path=matches[0] if len(matches) == 1 else None,
                             error=('Multiple artifacts; select buildDirectory or packageDirectory.' if len(matches) > 1
                                    else 'No compatible artifact found.' if not matches else None))
    return results


def fingerprint(path):
    st = Path(path).stat()
    return [st.st_dev, st.st_ino, st.st_size, st.st_mtime_ns, st.st_ctime_ns]


def check(name, filename):
    path = Path(filename).resolve()
    if not path.is_file():
        return {'state': 'unknown', 'error': 'Compiled file is missing.'}
    if artifact_name(path.name) != name:
        return {'state': 'unknown', 'error': 'Filename does not match this interpreter or module.'}
    before = fingerprint(path)
    sys.path.extend(p for p in dependency_paths() if p not in sys.path)
    # Do not let dependency imports accidentally enter this package or launch tools.
    class NoPackage(importlib.abc.MetaPathFinder):
        def find_spec(self, fullname, path=None, target=None):
            if fullname == 'mpyricalspace' or fullname.startswith('mpyricalspace.'):
                raise ImportError('Package imports are disabled in the build-status probe.')
    sys.meta_path.insert(0, NoPackage())
    def no_commands(event, args):
        if event in ('subprocess.Popen', 'os.system', 'os.posix_spawn', 'os.exec', 'os.fork', 'os.forkpty'):
            raise RuntimeError('Subprocess creation is disabled in the build-status probe.')
    sys.addaudithook(no_commands)
    try:
        spec = importlib.util.spec_from_file_location(name, path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        result = {'state': 'ok', 'error': None}
    except Exception as exc:
        result = {'state': 'failed', 'error': f'{type(exc).__name__}: {exc}'}
    try:
        if fingerprint(path) != before:
            return {'state': 'unknown', 'error': 'File changed during verification; refresh after the build finishes.'}
    except OSError:
        return {'state': 'unknown', 'error': 'File disappeared during verification.'}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    d = sub.add_parser('discover')
    d.add_argument('--root', required=True)
    d.add_argument('--build-dir', default='')
    d.add_argument('--package-dir', default='')
    c = sub.add_parser('check')
    c.add_argument('name')
    c.add_argument('path')
    args = parser.parse_args()
    try:
        result = (discover(args.root, args.build_dir, args.package_dir) if args.command == 'discover'
                  else check(args.name, args.path))
    except Exception as exc:
        result = {'error': f'{type(exc).__name__}: {exc}'}
    # Native libraries may write to stdout; use a framed result at the end.
    print('MPYRICALSPACE_PROBE=' + json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
