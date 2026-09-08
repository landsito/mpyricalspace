'''
python -m mpyricalspace doctor                                  which compiled extensions loaded
python -m mpyricalspace status                                  local-store coverage
python -m mpyricalspace fetch  --from 2024-01-01 --to 2024-06-01 [--kpap-source ...] [--ae-source ...]
python -m mpyricalspace update [--since 2024-01-01]             extend every series toward today
python -m mpyricalspace config --mongo-uri ... --data-dir ...
python -m mpyricalspace vscode [install|build|status]           companion VS Code extension
'''
import argparse
from datetime import datetime

p   = argparse.ArgumentParser(prog='mpyricalspace')
sub = p.add_subparsers(dest='cmd', required=True)
d = sub.add_parser('doctor'); d.add_argument('--json', action='store_true')
sub.add_parser('status')
u = sub.add_parser('update'); u.add_argument('--since')
f = sub.add_parser('fetch')
f.add_argument('--from', dest='d0', required=True)
f.add_argument('--to',   dest='dn', required=True)
f.add_argument('--kpap-source', default='potsdam')
f.add_argument('--ae-source',   default='auto')
c = sub.add_parser('config')
for opt in ('--mongo-uri', '--mongo-user', '--mongo-password', '--data-dir'):
    c.add_argument(opt)
c.add_argument('--mongo-writes', action='store_true', help='allow fetch/update to populate the mongo backend')
v = sub.add_parser('vscode', help='build / install the companion VS Code extension')
v.add_argument('action', nargs='?', default='install', choices=['install', 'build', 'status'])
v.add_argument('out', nargs='?', default='.', help="output dir for 'build' (default: .)")
a = p.parse_args()

if a.cmd == 'vscode':
    from mpyricalspace import _vscode
    if a.action == 'status':
        print(_vscode.status_text())
    elif a.action == 'build':
        print(_vscode.build_vsix(a.out))
    else:
        raise SystemExit(0 if _vscode.install() else 1)
    raise SystemExit

# first-run, best-effort, silent -- offer the VS Code extension (see mpyricalspace/_vscode.py)
try:
    from mpyricalspace._vscode import maybe_autoinstall
    maybe_autoinstall()
except Exception:
    pass

if a.cmd == 'doctor':
    from mpyricalspace._build import build_report, build_json
    print(build_json() if a.json else build_report())
    raise SystemExit

from mpyricalspace.DataManager import DataManager, configure

if a.cmd == 'config':
    print("wrote %s" % configure(mongo_uri=a.mongo_uri, mongo_user=a.mongo_user, mongo_password=a.mongo_password,
                                 data_dir=a.data_dir, mongo_writes=a.mongo_writes or None))
    raise SystemExit

dm = DataManager()
if a.cmd == 'status':
    dm.status()
elif a.cmd == 'fetch':
    dm.fetch(a.d0, a.dn, kpap_source=a.kpap_source, ae_source=a.ae_source)
    dm.status()
elif a.cmd == 'update':
    dm.update(since=datetime.fromisoformat(a.since) if a.since else None)
    try:
        from mpyricalspace.models import refresh_iri_indices
        refresh_iri_indices()
        print("refreshed IRI ig_rz.dat / apf107.dat")
    except Exception as e:
        print("IRI index refresh skipped: %s" % e)
    dm.status()
