const test = require('node:test');
const assert = require('node:assert/strict');
const vm = require('vm');
const fs = require('fs');
const path = require('path');
const { Observer } = require('../../mpyricalspace/vscode_ext/observer');
const directory = path.resolve(__dirname, '../../mpyricalspace/vscode_ext');
const sleep = ms => new Promise(resolve => setTimeout(resolve, ms));
async function until(predicate) {
  for (let i = 0; i < 250; i++) { if (predicate()) return; await sleep(10); }
  throw new Error('Integration condition timed out');
}
function fixture() {
  const watchers = [], calls = [], instances = [];
  const disposable = { dispose() {} };
  const settings = { pythonPath: '/selected/python', probeTimeout: 15, decorateBadges: true };
  class Runner {
    constructor(python, helper) { this.python = python; this.helper = helper; instances.push(this); }
    async run(args) {
      calls.push(args); await sleep(5);
      if (args[0] === 'discover') return {
        rocsatf: { path: '/project/build/src/rocsatf.so', subdir: 'rocsat' },
        hwm14f: { path: '/project/build/src/hwm14f.so', subdir: 'hwm/hwm14' },
      };
      return { state: args[1] === 'hwm14f' ? 'failed' : 'ok', error: args[1] === 'hwm14f' ? 'load failed' : null };
    }
    dispose() { this.disposed = true; }
  }
  const vscode = {
    Uri: { file: fsPath => ({ fsPath }) },
    EventEmitter: class { constructor() { this.event = () => disposable; } fire() {} dispose() {} },
    FileDecoration: class { constructor(badge, tooltip, color) { Object.assign(this, { badge, tooltip, color }); } },
    ThemeColor: class { constructor(id) { this.id = id; } },
    RelativePattern: class { constructor(base, pattern) { Object.assign(this, { base, pattern }); } },
    StatusBarAlignment: { Left: 1 },
    workspace: {
      workspaceFolders: [{ uri: { fsPath: '/project' } }],
      getConfiguration: () => ({ get: key => settings[key] }),
      onDidChangeConfiguration: () => disposable, onDidChangeWorkspaceFolders: () => disposable,
      createFileSystemWatcher(pattern) {
        const watcher = { pattern, handlers: {}, disposed: false,
          onDidCreate(fn) { this.handlers.create = fn; }, onDidChange(fn) { this.handlers.change = fn; },
          onDidDelete(fn) { this.handlers.delete = fn; }, dispose() { this.disposed = true; } };
        watchers.push(watcher); return watcher;
      },
    },
    window: { createStatusBarItem: () => ({ show() {}, dispose() {} }), registerFileDecorationProvider: () => disposable },
    commands: { registerCommand: () => disposable },
  };
  const context = { subscriptions: [] };
  const sandbox = { __dirname: directory, module: { exports: {} },
    require: name => name === 'vscode' ? vscode : name === 'fs' ? { existsSync: () => true }
      : name === './observer' ? { Observer, ProbeRunner: Runner }
        : name === './modules.json' ? require(directory + '/modules.json') : require(name) };
  vm.runInNewContext(fs.readFileSync(directory + '/extension.js', 'utf8'), sandbox);
  sandbox.module.exports.activate(context);
  return { context, watchers, calls, instances, provider: context.subscriptions[0], api: sandbox.module.exports };
}
test('VS Code wiring decorates per wrapper, invalidates checks, ignores unrelated files and disposes', async () => {
  const f = fixture();
  try {
    await until(() => f.provider.modules.hwm14f?.state === 'failed');
    assert.equal(f.instances[0].python, '/selected/python');
    assert.equal(f.provider.provideFileDecoration({ fsPath: '/project/src/rocsat' }).badge, '✓');
    assert.equal(f.provider.provideFileDecoration({ fsPath: '/project/src/hwm' }).badge, '✗');
    const watcher = f.watchers.findLast(w => !w.disposed && w.pattern.pattern === '**/*.{so,pyd}');
    const before = f.calls.length;
    watcher.handlers.change({ fsPath: '/project/build/unrelated.so' });
    await sleep(20); assert.equal(f.calls.length, before);
    watcher.handlers.change({ fsPath: '/project/build/src/rocsatf.so' });
    assert.equal(f.provider.provideFileDecoration({ fsPath: '/project/src/rocsat' }).badge, undefined);
    await until(() => f.provider.modules.rocsatf.state === 'ok');
    assert.equal(f.calls.length, before + 1);
    assert.equal(f.calls.at(-1)[1], 'rocsatf');
    assert.ok(f.calls.every(args => ['discover', 'check'].includes(args[0])));
  } finally { f.api.deactivate(); }
  assert.ok(f.watchers.every(w => w.disposed)); assert.ok(f.instances.every(i => i.disposed));
});
