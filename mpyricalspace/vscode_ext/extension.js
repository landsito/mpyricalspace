// Observe existing native artifacts. Never import the package or start a build.
const vscode = require('vscode');
const fs = require('fs');
const path = require('path');
const { ProbeRunner, Observer } = require('./observer');
const registry = require('./modules.json');
let provider;

function projectRoot() {
  return (vscode.workspace.workspaceFolders || []).map(f => f.uri.fsPath)
    .find(p => fs.existsSync(path.join(p, 'mpyricalspace', '_build.py')) && fs.existsSync(path.join(p, 'src', 'meson.build')));
}

async function pythonInterpreter(root, subscribe) {
  const configured = vscode.workspace.getConfiguration('mpyricalspace', vscode.Uri.file(root)).get('pythonPath');
  if (configured) return configured;
  const extension = vscode.extensions.getExtension('ms-python.python');
  if (extension) {
    const api = extension.isActive ? extension.exports : await extension.activate();
    if (api.environments && api.environments.onDidChangeActiveEnvironmentPath) subscribe(api.environments);
    const resource = vscode.Uri.file(root);
    const selected = api.environments && api.environments.getActiveEnvironmentPath(resource);
    if (selected) {
      const environment = await api.environments.resolveEnvironment(selected);
      if (environment && environment.executable && environment.executable.uri) return environment.executable.uri.fsPath;
    }
    const details = api.settings && api.settings.getExecutionDetails(resource);
    if (details && details.execCommand && details.execCommand.length) return details.execCommand[0];
  }
  return 'python3';
}

class BuildStatusProvider {
  constructor(context) {
    this.context = context; this.root = projectRoot(); this.modules = {}; this.watchers = [];
    this.observer = null; this.disposed = false; this.restartId = 0; this.restarts = Promise.resolve();
    this.emitter = new vscode.EventEmitter(); this.onDidChangeFileDecorations = this.emitter.event;
    this.status = vscode.window.createStatusBarItem(vscode.StatusBarAlignment.Left, 5);
    this.status.command = 'mpyricalspace.refreshBuildStatus';
  }
  restart() {
    const id = ++this.restartId;
    // Serialize replacement as well as probes: the previous process must close first.
    this.restarts = this.restarts.then(async () => {
      if (this.disposed || id !== this.restartId) return;
      if (this.observer) {
        const child = this.observer.runner.child;
        const closed = child ? new Promise(resolve => child.once('close', resolve)) : Promise.resolve();
        this.observer.dispose(); this.observer = null;
        await closed;
      }
      this.clearWatchers(); this.modules = {}; this.root = projectRoot();
      this.update({}, 'Selecting Python interpreter.');
      if (!this.root) { this.update({}, 'No mpyricalspace source project in this workspace.'); return; }
      try {
        const python = await pythonInterpreter(this.root, environments => {
          if (!this.pythonEvents && !this.disposed) this.pythonEvents = environments.onDidChangeActiveEnvironmentPath(() => this.restart());
        });
        if (this.disposed || id !== this.restartId) return;
        const config = vscode.workspace.getConfiguration('mpyricalspace', vscode.Uri.file(this.root));
        const args = ['--root', this.root];
        for (const [key, flag] of [['buildDirectory', '--build-dir'], ['packageDirectory', '--package-dir']]) {
          const value = config.get(key);
          if (value) args.push(flag, path.resolve(this.root, value));
        }
        const runner = new ProbeRunner(python, path.join(__dirname, 'probe.py'), config.get('probeTimeout') * 1000);
        this.observer = new Observer(runner, args, (modules, note) => this.update(modules, note), modules => this.watch(modules));
        this.watch({}); this.observer.refresh();
      } catch (error) { this.update({}, error.message); }
    });
    return this.restarts;
  }
  refresh() { if (this.observer) this.observer.refresh(); else this.restart(); }
  clearWatchers() { this.watchers.forEach(w => w.dispose()); this.watchers = []; }
  watch(modules) {
    this.clearWatchers();
    const directories = new Set(Object.values(modules).filter(m => m.path).map(m => path.dirname(m.path)));
    // Watch the project for new artifacts/install plans, including previously missing modules.
    const patterns = [new vscode.RelativePattern(this.root, '**/*.{so,pyd}'),
      new vscode.RelativePattern(this.root, '**/meson-info/intro-install_plan.json')];
    const config = vscode.workspace.getConfiguration('mpyricalspace', vscode.Uri.file(this.root));
    for (const key of ['buildDirectory', 'packageDirectory']) {
      const value = config.get(key);
      if (value) {
        const base = path.resolve(this.root, value);
        patterns.push(new vscode.RelativePattern(base, '**/*.{so,pyd}'));
        if (key === 'buildDirectory') patterns.push(new vscode.RelativePattern(base, 'meson-info/intro-install_plan.json'));
      }
    }
    for (const dir of directories) {
      const relative = path.relative(this.root, dir);
      if (relative.startsWith('..') || path.isAbsolute(relative)) patterns.push(new vscode.RelativePattern(dir, '*.{so,pyd}'));
    }
    for (const pattern of patterns) {
      const watcher = vscode.workspace.createFileSystemWatcher(pattern);
      const changed = uri => {
        if (!this.observer) return;
        const filename = uri.fsPath;
        if (path.basename(filename) === 'intro-install_plan.json') { this.observer.refresh(); return; }
        if (filename.endsWith('.json')) return;
        if (Object.values(this.modules).some(m => m.path === filename)) this.observer.changed(filename);
        else {
          const name = path.basename(filename).split('.')[0];
          if (registry[name]) this.observer.refresh();
        }
      };
      watcher.onDidCreate(changed); watcher.onDidChange(changed); watcher.onDidDelete(changed);
      this.watchers.push(watcher);
    }
  }
  update(modules, note) {
    if (this.disposed) return;
    this.modules = modules;
    const values = Object.values(modules), ok = values.filter(m => m.state === 'ok').length;
    this.status.text = `${note ? '$(circle-slash)' : '$(beaker)'} mpyricalspace ${ok}/${values.length}`;
    this.status.tooltip = note || `${ok}/${values.length} existing wrappers loaded. Click to verify again; no build is started.`;
    this.status.show(); this.emitter.fire(undefined);
  }
  provideFileDecoration(uri) {
    if (!this.root) return undefined;
    const rel = path.relative(path.join(this.root, 'src'), uri.fsPath).split(path.sep).join('/');
    if (!rel || rel.startsWith('..')) return undefined;
    const entries = Object.entries(this.modules).filter(([, m]) => m.subdir === rel || m.subdir.startsWith(rel + '/'));
    if (!entries.length) return undefined;
    const badges = vscode.workspace.getConfiguration('mpyricalspace', uri).get('decorateBadges');
    const failed = entries.filter(([, m]) => m.state === 'failed');
    const ok = entries.filter(([, m]) => m.state === 'ok').length;
    const allOK = ok === entries.length;
    const tooltip = entries.length === 1 ? `${entries[0][0]}: ${allOK ? 'Python loaded this compiled file.' : entries[0][1].error || 'Not verified.'}`
      : `${ok}/${entries.length} wrappers loaded; ${failed.length} failed; ${entries.length - ok - failed.length} not verified.`;
    return new vscode.FileDecoration(badges ? (failed.length ? '✗' : allOK ? '✓' : undefined) : undefined,
      tooltip, failed.length ? new vscode.ThemeColor('gitDecoration.deletedResourceForeground')
        : allOK ? new vscode.ThemeColor('gitDecoration.addedResourceForeground') : undefined);
  }
  dispose() {
    this.disposed = true; this.restartId++;
    if (this.observer) this.observer.dispose();
    this.clearWatchers();
    if (this.pythonEvents) this.pythonEvents.dispose();
    this.status.dispose(); this.emitter.dispose();
  }
}

function activate(context) {
  provider = new BuildStatusProvider(context);
  context.subscriptions.push(provider, vscode.window.registerFileDecorationProvider(provider),
    vscode.commands.registerCommand('mpyricalspace.refreshBuildStatus', () => provider.refresh()),
    vscode.workspace.onDidChangeConfiguration(e => {
      if (e.affectsConfiguration('mpyricalspace') || e.affectsConfiguration('python')) provider.restart();
    }), vscode.workspace.onDidChangeWorkspaceFolders(() => provider.restart()));
  provider.restart();
}
function deactivate() { if (provider) provider.dispose(); }
module.exports = { activate, deactivate };
