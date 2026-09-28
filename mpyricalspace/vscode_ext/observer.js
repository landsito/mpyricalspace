// No VS Code dependency: process isolation and scheduling can be tested independently.
const cp = require('child_process');
const path = require('path');

class ProbeRunner {
  constructor(python, helper, timeout = 15000) {
    this.python = python;
    this.helper = helper;
    this.timeout = timeout;
    this.child = null;
    this.disposed = false;
  }
  run(args) {
    if (this.disposed) return Promise.reject(new Error('Observer disposed.'));
    if (this.child) return Promise.reject(new Error('A probe is already running.'));
    return new Promise((resolve, reject) => {
      const env = { ...process.env, MPYRICALSPACE_NO_VSCODE: '1', PYTHONDONTWRITEBYTECODE: '1' };
      // -S prevents site/.pth startup, including meson-python editable loaders.
      const child = cp.spawn(this.python, ['-I', '-S', '-B', this.helper, ...args], {
        env, cwd: path.dirname(this.helper), windowsHide: true,
        detached: process.platform !== 'win32', stdio: ['ignore', 'pipe', 'pipe'],
      });
      this.child = child;
      let stdout = '', stderr = '', failure;
      const timer = setTimeout(() => { failure = new Error('Verification timed out.'); this.cancel(); }, this.timeout);
      const collect = (kind, chunk) => {
        if (kind === 'out') stdout += chunk.toString(); else stderr += chunk.toString();
        if (stdout.length + stderr.length > 1024 * 1024) {
          failure = new Error('Probe output exceeded 1 MiB.'); this.cancel();
        }
      };
      child.stdout.on('data', chunk => collect('out', chunk));
      child.stderr.on('data', chunk => collect('err', chunk));
      child.on('error', error => { failure = error; });
      child.on('close', (code, signal) => {
        clearTimeout(timer);
        this.child = null;
        if (this.disposed) return reject(new Error('Observer disposed.'));
        if (failure) return reject(failure);
        if (code !== 0) return reject(new Error(`Probe exited ${signal || code}: ${stderr.trim().slice(-4000)}`));
        try {
          const frame = stdout.lastIndexOf('MPYRICALSPACE_PROBE=');
          if (frame < 0) throw new Error('Probe returned no result.');
          resolve(JSON.parse(stdout.slice(frame + 'MPYRICALSPACE_PROBE='.length).split('\n')[0]));
        } catch (error) { reject(error); }
      });
    });
  }
  cancel() {
    if (!this.child) return;
    try {
      if (process.platform !== 'win32' && this.child.pid) process.kill(-this.child.pid, 'SIGKILL');
      else this.child.kill('SIGKILL');
    } catch (_) { /* process already exited */ }
  }
  dispose() { this.disposed = true; this.cancel(); }
}

class Observer {
  constructor(runner, discoverArgs, onChange, onDiscover = () => {}, delay = 400) {
    this.runner = runner; this.discoverArgs = discoverArgs;
    this.onChange = onChange; this.onDiscover = onDiscover; this.delay = delay;
    this.modules = {}; this.pending = new Set(); this.versions = new Map();
    this.full = false; this.busy = false; this.disposed = false; this.timer = null;
    this.note = null; this.generation = 0;
  }
  emit() { if (!this.disposed) this.onChange(this.modules, this.note); }
  refresh() {
    if (this.disposed) return;
    this.full = true; this.generation++;
    for (const name of Object.keys(this.modules)) this.invalidate(name);
    this.emit(); this.schedule();
  }
  invalidate(name) {
    this.versions.set(name, (this.versions.get(name) || 0) + 1);
    Object.assign(this.modules[name], { state: 'unknown', error: 'Pending verification.' });
    this.pending.add(name);
  }
  changed(filename) {
    if (this.disposed) return;
    const names = Object.keys(this.modules).filter(name => this.modules[name].path === filename);
    if (!names.length) return; // Unrelated artifacts do not trigger probes.
    names.forEach(name => this.invalidate(name));
    this.emit(); this.schedule();
  }
  schedule() {
    clearTimeout(this.timer);
    this.timer = setTimeout(() => { this.timer = null; this.drain(); }, this.delay);
  }
  async drain() {
    if (this.busy || this.disposed) return;
    this.busy = true;
    try {
      while (!this.disposed && (this.full || this.pending.size)) {
        if (this.timer) break;
        if (this.full) {
          this.full = false;
          const generation = this.generation;
          try {
            const data = await this.runner.run(['discover', ...this.discoverArgs]);
            if (this.disposed) break;
            if (generation !== this.generation) continue;
            if (data.error) throw new Error(data.error);
            this.note = null;
            this.modules = Object.fromEntries(Object.entries(data).map(([name, info]) =>
              [name, { ...info, state: 'unknown', error: info.error || 'Pending verification.' }]));
            this.pending.clear();
            for (const [name, info] of Object.entries(this.modules)) if (info.path) this.pending.add(name);
            this.onDiscover(this.modules); this.emit();
          } catch (error) {
            if (generation !== this.generation) continue;
            this.note = error.message; this.pending.clear(); this.emit();
          }
        } else {
          const name = this.pending.values().next().value;
          this.pending.delete(name);
          const info = this.modules[name];
          if (!info || !info.path) continue;
          const version = this.versions.get(name);
          const generation = this.generation;
          let result;
          try { result = await this.runner.run(['check', name, info.path]); }
          catch (error) { result = { state: 'failed', error: error.message }; }
          if (this.disposed) break;
          if (generation !== this.generation || version !== this.versions.get(name)) continue;
          if (!['ok', 'failed', 'unknown'].includes(result.state)) result = { state: 'unknown', error: result.error || 'Invalid probe result.' };
          Object.assign(info, result); this.emit();
        }
        // Wait for a burst of events to settle before starting another probe.
        if (this.timer) break;
      }
    } finally { this.busy = false; }
  }
  dispose() {
    this.disposed = true; clearTimeout(this.timer); this.pending.clear(); this.full = false;
    this.runner.dispose();
  }
}
module.exports = { ProbeRunner, Observer };
