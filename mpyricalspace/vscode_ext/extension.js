/*
 * mpyricalspace build status -- an optional VS Code nicety.
 *
 * Runs `python -m mpyricalspace doctor --json` and, in the Explorer, tints each
 * src/<model>/ folder green (extension compiled and loads) or red (missing /
 * import error), with a matching badge and tooltip -- and any ancestor folder that
 * isn't itself a model (e.g. src/iri/, which several IRI versions live under) gets
 * a rolled-up n/m summary instead, so nothing under src/ is ever left uncolored.
 * A status-bar item shows the overall n/m count.  Everything here is editor-only:
 * nothing is written, and if
 * mpyricalspace is not installed the extension simply shows "not installed" and
 * decorates nothing.
 */
const vscode = require('vscode');
const cp = require('child_process');
const fs = require('fs');
const path = require('path');
const util = require('util');

const execFile = util.promisify(cp.execFile);

let provider;
let statusBar;

function pythonInterpreter() {
  const configured = vscode.workspace.getConfiguration('mpyricalspace').get('pythonPath');
  if (configured) return configured;
  try {
    const py = vscode.extensions.getExtension('ms-python.python');
    const details = py && py.isActive && py.exports && py.exports.settings
      && py.exports.settings.getExecutionDetails
      ? py.exports.settings.getExecutionDetails()
      : null;
    if (details && Array.isArray(details.execCommand) && details.execCommand.length) {
      return details.execCommand[0];
    }
  } catch (_e) { /* fall through */ }
  return 'python3';
}

function projectRoot() {
  return (vscode.workspace.workspaceFolders || [])
    .map((f) => f.uri.fsPath)
    .find((p) => fs.existsSync(path.join(p, 'src', 'meson.build'))) || null;
}

class BuildStatusProvider {
  constructor() {
    this._emitter = new vscode.EventEmitter();
    this.onDidChangeFileDecorations = this._emitter.event;
    this.root = null;
    this.bySubdir = {};      // "hltwim" -> { ok, error, label, ext }
    this.byPrefix = {};      // "iri" -> { ok, total } -- rollup for ancestor folders, see _rollup()
    this.summary = null;     // { ok, total } | { note }
    this._timer = null;
  }

  scheduleRefresh() {
    if (this._timer) clearTimeout(this._timer);
    this._timer = setTimeout(() => this.refresh(), 400);   // debounce a build touching many .so
  }

  async refresh() {
    this.root = projectRoot();
    if (!this.root) { this._apply({}, { note: 'no src/meson.build in the workspace' }); return; }

    const py = pythonInterpreter();
    try {
      const { stdout } = await execFile(py, ['-m', 'mpyricalspace', 'doctor', '--json'],
        { cwd: this.root, timeout: 30000, maxBuffer: 4 * 1024 * 1024 });
      const data = JSON.parse(stdout);
      const bySubdir = {};
      for (const [ext, info] of Object.entries(data.extensions || {})) {
        bySubdir[info.subdir] = { ok: info.ok, error: info.error, label: info.label, ext };
      }
      this._apply(bySubdir, { ok: data.ok, total: data.total });
    } catch (e) {
      const msg = /No module named mpyricalspace/i.test(String(e && e.stderr))
        ? 'mpyricalspace is not installed in ' + py
        : String((e && e.message) || e).split('\n')[0];
      this._apply({}, { note: msg });
    }
  }

  // Every ancestor prefix of a leaf subdir (e.g. "iri" for "iri/iri12") gets its own
  // ok/total tally across every leaf under it, so a container folder (nothing itself
  // builds into an extension, e.g. src/iri/) still gets a summary decoration without
  // needing to expand it -- same idea as VS Code's own Git decorations propagating a
  // changed file's color up to its parent folders.
  _rollup(bySubdir) {
    const byPrefix = {};
    for (const [subdir, info] of Object.entries(bySubdir)) {
      const parts = subdir.split('/');
      for (let i = 1; i < parts.length; i++) {
        const prefix = parts.slice(0, i).join('/');
        const agg = byPrefix[prefix] || (byPrefix[prefix] = { ok: 0, total: 0 });
        agg.total += 1;
        if (info.ok) agg.ok += 1;
      }
    }
    return byPrefix;
  }

  _apply(bySubdir, summary) {
    this.bySubdir = bySubdir;
    this.byPrefix = this._rollup(bySubdir);
    this.summary = summary;
    this._updateStatusBar();
    this._emitter.fire(undefined);          // re-decorate everything
  }

  _updateStatusBar() {
    if (!statusBar) return;
    const s = this.summary || {};
    if (s.note) {
      statusBar.text = '$(circle-slash) mpyricalspace';
      statusBar.tooltip = 'mpyricalspace build status: ' + s.note;
    } else {
      const bad = (s.total || 0) - (s.ok || 0);
      statusBar.text = `${bad ? '$(error)' : '$(check)'} mpyricalspace ${s.ok}/${s.total}`;
      statusBar.tooltip = bad
        ? `${bad} compiled extension(s) missing -- run  python -m mpyricalspace doctor`
        : 'every compiled Fortran extension loads';
    }
    statusBar.command = 'mpyricalspace.refreshBuildStatus';
    statusBar.show();
  }

  provideFileDecoration(uri) {
    if (!this.root) return undefined;
    // subdir values from `doctor --json` are src/-relative and not always a direct
    // child (e.g. the IRI versions live under src/iri/<version>/), so match on the
    // relative path as-is rather than rejecting anything with a path separator.
    const rel = path.relative(path.join(this.root, 'src'), uri.fsPath);
    if (!rel || rel.startsWith('..')) return undefined;

    const badges = vscode.workspace.getConfiguration('mpyricalspace').get('decorateBadges');
    const info = this.bySubdir[rel];
    if (info) {
      if (info.ok) {
        return new vscode.FileDecoration(badges ? '✓' : undefined,
          `${info.ext}  --  compiled, loads OK`,
          new vscode.ThemeColor('gitDecoration.addedResourceForeground'));
      }
      return new vscode.FileDecoration(badges ? '✗' : undefined,
        `${info.ext}  --  NOT built: ${info.error || 'missing'}`,
        new vscode.ThemeColor('gitDecoration.deletedResourceForeground'));
    }

    // Not a leaf itself (e.g. src/iri/, nothing builds into "iri") -- roll up whatever
    // leaves live under it instead, so it's not left with no indicator at all.
    const agg = this.byPrefix[rel];
    if (!agg) return undefined;
    const bad = agg.total - agg.ok;
    if (!bad) {
      return new vscode.FileDecoration(badges ? '✓' : undefined,
        `${agg.ok}/${agg.total} compiled, load OK`,
        new vscode.ThemeColor('gitDecoration.addedResourceForeground'));
    }
    return new vscode.FileDecoration(badges ? '✗' : undefined,
      `${bad}/${agg.total} NOT built`,
      new vscode.ThemeColor('gitDecoration.deletedResourceForeground'));
  }
}

function activate(context) {
  provider = new BuildStatusProvider();
  statusBar = vscode.window.createStatusBarItem(vscode.StatusBarAlignment.Left, 5);

  context.subscriptions.push(
    statusBar,
    vscode.window.registerFileDecorationProvider(provider),
    vscode.commands.registerCommand('mpyricalspace.refreshBuildStatus', () => provider.refresh()),
    vscode.tasks.onDidEndTask(() => provider.scheduleRefresh()),
    vscode.workspace.onDidChangeConfiguration((e) => {
      if (e.affectsConfiguration('mpyricalspace')) provider.refresh();
    })
  );

  const watcher = vscode.workspace.createFileSystemWatcher('**/*.{so,pyd}');
  watcher.onDidCreate(() => provider.scheduleRefresh());
  watcher.onDidDelete(() => provider.scheduleRefresh());
  watcher.onDidChange(() => provider.scheduleRefresh());
  context.subscriptions.push(watcher);

  provider.refresh();
}

function deactivate() {}

module.exports = { activate, deactivate };
