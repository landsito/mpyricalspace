# mpyricalspace build status (VS Code)

An **optional** editor nicety bundled with the `mpyricalspace` package. It runs

```
python -m mpyricalspace doctor --json
```

and, in the Explorer, tints each `src/<model>/` folder:

| | |
|---|---|
| green + `✓` | the model's compiled Fortran extension is built and imports |
| red + `✗`   | it is missing / fails to import (tooltip shows the error) |

A model can live nested more than one level under `src/` (e.g. `src/iri/iri12/`)
-- the exact path just needs to match `doctor --json`'s own `subdir` for that
extension. Any folder that isn't itself a model (e.g. `src/iri/`, which several
IRI versions live under) gets a rolled-up `n/m` summary instead of no color at
all, tallied from every model folder underneath it.

A status-bar item (bottom-left) shows the overall `n/m` count; click it to
refresh. It also refreshes when a `.so` changes or a build task ends.

100% editor-only: it writes nothing, touches neither the package nor the build,
and if `mpyricalspace` is not installed it just shows "not installed" and
decorates nothing.

## Install

It is **attempted automatically, once**, the first time you run
`python -m mpyricalspace ...` or `import mpyricalspace` from a VS Code terminal
(wheels have no post-install hook). Failures are logged to
`~/.cache/mpyricalspace/vscode-ext.log` and never surface.

Do it (or redo it) by hand:

```bash
python -m mpyricalspace vscode            # build the .vsix + code --install-extension
python -m mpyricalspace vscode status     # installed?  was it tried?
python -m mpyricalspace vscode build DIR  # just write the .vsix into DIR
```

**Fully quit VS Code (Cmd/Ctrl+Q) and reopen** after the first install -- a window
reload is not always enough. If the `code` CLI is missing, VS Code ->
Cmd/Ctrl+Shift+P -> "Shell Command: Install 'code' command in PATH".

Opt out of the auto-attempt entirely with `MPYRICALSPACE_NO_VSCODE=1`.
Uninstall: `code --uninstall-extension mpyricalspace.mpyricalspace-build-status`.

## Settings

| setting | default | meaning |
|---|---|---|
| `mpyricalspace.pythonPath` | `""` | interpreter for `doctor`; empty = the Python extension's selected interpreter, else `python3` |
| `mpyricalspace.decorateBadges` | `true` | show the `✓`/`✗` badge (colour is applied either way) |

## Hack on it

Open this folder in VS Code and press **F5** (Extension Development Host).
