# mpyricalspace build status (VS Code)

An optional observer for compiled wrappers. Version **0.2.0** replaces the old
`doctor` subprocess with a standalone, per-file load check. It never imports
`mpyricalspace`, invokes a build, installs dependencies, or repairs the package.

| Explorer indicator | Meaning |
|---|---|
| Green ✓ | The selected Python interpreter loaded this compiled file successfully. |
| Red ✗ | Loading failed, the probe exited unexpectedly, or it timed out. Hover for the error. |
| No check | Missing, ambiguous, pending, or not yet verified. Hover for details. |

These checks do not run model functions or validate scientific results. A check
also does not prove that an artifact is up to date with its source. It describes
the existing file that was loaded. Ancestor folders show an aggregate count;
the status bar shows the number of wrappers that loaded.

## How it works

1. Read a bundled module-to-source-folder registry and existing Meson install
   plans, or locate an installed package directory without importing it.
2. Start a fresh Python process for **one artifact**, loading it by absolute path.
3. Report the result and exit. There is only one active probe per extension
   instance; repeated events are grouped into a pending check.

Python runs with `-I -S -B`: no `PYTHONPATH`, startup customization, `.pth` files,
editable import hooks, or bytecode writes. The helper adds dependency directories
without executing their startup hooks, blocks imports of `mpyricalspace`, and
rejects Python subprocess-launch requests. Native code is still executed when a
wrapper loads: process isolation is not an OS sandbox or a promise of zero CPU,
memory, or native-code side effects. Native crashes stay in the probe process.

Only changed artifacts are rechecked. Their old check is invalidated immediately;
a result is discarded if the file changes during verification. New artifacts or
changed install plans trigger rediscovery. Manual refresh and interpreter/settings
changes recheck the selection. There is no polling, automatic repair, or refresh
on unrelated task completion. A timeout terminates the probe; on POSIX the entire
probe process group is terminated. Deactivation cancels timers and queued work
and terminates the active probe.

The extension targets the first workspace folder containing both
`mpyricalspace/_build.py` and `src/meson.build`. It decorates that source tree even
when checking a normal installed package. It requires a trusted workspace.

## Artifact discovery and settings

| Setting | Default | Meaning |
|---|---|---|
| `mpyricalspace.pythonPath` | `""` | Python executable. Empty uses the Python extension's selected environment, falling back to `python3`. |
| `mpyricalspace.buildDirectory` | `""` | Existing Meson build directory, e.g. `build/cp312`; relative to the source root or absolute. |
| `mpyricalspace.packageDirectory` | `""` | Directory containing installed wrappers; overrides build discovery. Useful for custom installation layouts. |
| `mpyricalspace.probeTimeout` | `15` | Maximum seconds per probe (1–120). |
| `mpyricalspace.decorateBadges` | `true` | Show ✓/✗ badges as well as folder colors. |

By default, compatible artifacts in `build/**/meson-info/intro-install_plan.json`
are preferred. If none are mapped, the selected interpreter's standard package
directories are searched. Multiple candidates remain unverified: select a build
or package directory explicitly. Incompatible filename tags are ignored. A missing
artifact in a selected build does not silently fall back to an older installed copy.
Custom dependency paths supplied only by `.pth` files are deliberately not executed;
use an interpreter with the needed dependencies in its standard package directories.

## Build and install without importing the library

From the repository root, using the desired Python executable:

```bash
python -I -S -B mpyricalspace/_vscode.py build
code --install-extension ./mpyricalspace-build-status-0.2.0.vsix --force
```

The build command writes the VSIX to the current directory. Optionally pass an
existing output directory after `build`. After installation, the VSIX can be
deleted; VS Code keeps its own installed copy.

Alternatively use **Extensions → … → Install from VSIX…**. Reload VS Code after
updating. If the extension was disabled, enable it explicitly when ready to test.
Disabling or uninstalling it does not change the Python package or build files.

If the Python package is installed, these CLI commands are also available:

```bash
python -m mpyricalspace vscode              # build in a temporary directory and install
python -m mpyricalspace vscode build        # build the VSIX in the current directory
python -m mpyricalspace vscode status       # show installation status and log location
```

The first command cleans up its temporary build directory automatically. The
`build` command optionally accepts an existing output directory. These package
entry points remain supported for compatibility, but can trigger an editable
rebuild; use the direct script above to avoid importing the library. The observer
never calls them. The legacy automatic installation attempt runs only in a VS Code
environment with an available editor CLI. Set `MPYRICALSPACE_NO_VSCODE=1` to opt
out of that automatic attempt. Uninstall with:

```bash
code --uninstall-extension mpyricalspace.mpyricalspace-build-status
```

## Development and validation

`modules.json` maps wrapper names to folders, matching `_build.py`'s registry;
update both when adding a wrapper. `probe.py` is standalone and ships inside the
VSIX. `observer.js` handles serialization, cancellation, and stale results;
`extension.js` handles VS Code events and decorations.

Run from the repository root, without importing the package:

```bash
python -I -S -B tests/vscode/test_probe.py
node --test tests/vscode/observer.test.js tests/vscode/extension.test.js
```

Set `MPY_TEST_PYTHON` to a Python executable for the Node process tests; it defaults
to `python3`. The optional acceptance check below loads real local artifacts and
verifies that build-file sizes and modification times remain unchanged. It does
not compile missing wrappers and fails if any wrapper cannot load:

```bash
python -I -S -B tests/vscode/live_check.py
```

When checking editable isolation, do not use `python -m mpyricalspace doctor`:
that is a separate manual package diagnostic and may trigger a build.
