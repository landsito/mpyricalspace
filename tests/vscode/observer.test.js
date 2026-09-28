const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('fs');
const os = require('os');
const path = require('path');
const { Observer, ProbeRunner } = require('../../mpyricalspace/vscode_ext/observer');
const delay = ms => new Promise(resolve => setTimeout(resolve, ms));
async function until(predicate) {
  for (let i = 0; i < 300; i++) { if (predicate()) return; await delay(10); }
  throw new Error('Timed out waiting for test condition');
}
function fixture() {
  const calls = []; let active = 0, maximum = 0;
  const runner = {
    disposed: false,
    async run(args) {
      calls.push(args); active++; maximum = Math.max(maximum, active);
      await delay(40); active--;
      if (args[0] === 'discover') return {
        rocsatf: { path: '/build/rocsatf.so', subdir: 'rocsat' },
        hwm14f: { path: '/build/hwm14f.so', subdir: 'hwm/hwm14' },
      };
      return { state: 'ok', error: null };
    },
    dispose() { this.disposed = true; },
  };
  return { runner, calls, maximum: () => maximum };
}
test('bursts coalesce, only changed module is rechecked, no overlap', async () => {
  const f = fixture(); const observer = new Observer(f.runner, [], () => {}, () => {}, 5);
  observer.refresh();
  await until(() => observer.modules.hwm14f?.state === 'ok');
  const count = f.calls.length;
  observer.changed('/unrelated/test.so'); await delay(20); assert.equal(f.calls.length, count);
  for (let i = 0; i < 20; i++) observer.changed('/build/rocsatf.so');
  assert.equal(observer.modules.rocsatf.state, 'unknown');
  await until(() => observer.modules.rocsatf.state === 'ok');
  assert.equal(f.calls.length, count + 1); assert.equal(f.maximum(), 1);
  observer.dispose();
});
test('a file change during a check discards the old result', async () => {
  const f = fixture(); const observer = new Observer(f.runner, [], () => {}, () => {}, 5);
  observer.refresh(); await until(() => observer.modules.hwm14f?.state === 'ok');
  observer.changed('/build/rocsatf.so');
  await until(() => observer.busy); observer.changed('/build/rocsatf.so');
  await until(() => !observer.busy && observer.modules.rocsatf.state === 'ok');
  assert.equal(f.calls.filter(a => a[0] === 'check' && a[1] === 'rocsatf').length, 3);
  assert.equal(f.maximum(), 1); observer.dispose();
});
test('manual refreshes coalesce while busy, disposal cancels pending work', async () => {
  const f = fixture(); let updates = 0;
  const observer = new Observer(f.runner, [], () => updates++, () => {}, 5);
  observer.refresh(); await until(() => observer.busy);
  for (let i = 0; i < 10; i++) observer.refresh();
  await until(() => observer.modules.hwm14f?.state === 'ok');
  assert.equal(f.maximum(), 1); assert.equal(f.calls.filter(a => a[0] === 'discover').length, 2);
  observer.changed('/build/rocsatf.so'); observer.dispose();
  const before = updates, calls = f.calls.length; await delay(70);
  assert.equal(updates, before); assert.equal(f.calls.length, calls); assert.equal(f.runner.disposed, true);
});
const python = process.env.MPY_TEST_PYTHON || 'python3';
test('runner timeout kills the isolated process and permits a later probe', async () => {
  const dir = fs.mkdtempSync(path.join(os.tmpdir(), 'mpy-probe-test-'));
  const helper = path.join(dir, 'helper.py');
  fs.writeFileSync(helper, 'import time\ntime.sleep(30)\n');
  const runner = new ProbeRunner(python, helper, 200);
  try {
    await assert.rejects(runner.run([]), /timed out/); assert.equal(runner.child, null);
    fs.writeFileSync(helper, 'print(\'MPYRICALSPACE_PROBE={"state":"ok"}\')\n');
    runner.timeout = 3000; assert.equal((await runner.run([])).state, 'ok');
  } finally { runner.dispose(); fs.rmSync(dir, { recursive: true, force: true }); }
});
test('runner disposal stops an active probe and rejects future requests', async () => {
  const dir = fs.mkdtempSync(path.join(os.tmpdir(), 'mpy-probe-test-'));
  const helper = path.join(dir, 'helper.py'); fs.writeFileSync(helper, 'import time\ntime.sleep(30)\n');
  const runner = new ProbeRunner(python, helper, 3000);
  try {
    const promise = runner.run([]); runner.dispose();
    await assert.rejects(promise, /disposed/); assert.equal(runner.child, null);
    await assert.rejects(runner.run([]), /disposed/);
  } finally { runner.dispose(); fs.rmSync(dir, { recursive: true, force: true }); }
});
test('native-process failure and malformed output stay inside the runner', async () => {
  const dir = fs.mkdtempSync(path.join(os.tmpdir(), 'mpy-probe-test-'));
  const helper = path.join(dir, 'helper.py'); const runner = new ProbeRunner(python, helper, 3000);
  try {
    fs.writeFileSync(helper, 'import os\nos._exit(7)\n');
    await assert.rejects(runner.run([]), /exited 7/);
    fs.writeFileSync(helper, 'print("not JSON")\n');
    await assert.rejects(runner.run([]), /no result/);
    fs.writeFileSync(helper, 'print(\'native output\\nMPYRICALSPACE_PROBE={"state":"ok"}\')\n');
    assert.equal((await runner.run([])).state, 'ok');
  } finally { runner.dispose(); fs.rmSync(dir, { recursive: true, force: true }); }
});
test('a failed module does not stop checks of other modules', async () => {
  const f = fixture(), original = f.runner.run;
  f.runner.run = async args => {
    const result = await original(args);
    if (args[0] === 'check' && args[1] === 'rocsatf') throw new Error('could not load');
    return result;
  };
  const observer = new Observer(f.runner, [], () => {}, () => {}, 5);
  observer.refresh(); await until(() => observer.modules.hwm14f?.state === 'ok');
  assert.equal(observer.modules.rocsatf.state, 'failed'); observer.dispose();
});
test('POSIX timeout also stops descendants in the probe process group', { skip: process.platform === 'win32', timeout: 4000 }, async () => {
  const dir = fs.mkdtempSync(path.join(os.tmpdir(), 'mpy-probe-tree-'));
  const heartbeat = path.join(dir, 'heartbeat');
  const helper = path.join(dir, 'helper.py');
  const childCode = `import pathlib,time\np=pathlib.Path(${JSON.stringify(heartbeat)})\nwhile True:\n p.write_text(str(time.time()))\n time.sleep(0.02)`;
  fs.writeFileSync(helper, `import subprocess,sys,time\nsubprocess.Popen([sys.executable,'-I','-S','-c',${JSON.stringify(childCode)}])\ntime.sleep(30)\n`);
  const runner = new ProbeRunner(python, helper, 600);
  try {
    const ended = assert.rejects(runner.run([]), /timed out/);
    await until(() => fs.existsSync(heartbeat)); await ended;
    const last = fs.readFileSync(heartbeat, 'utf8'); await delay(100);
    assert.equal(fs.readFileSync(heartbeat, 'utf8'), last);
    assert.equal(runner.child, null);
  } finally { runner.dispose(); fs.rmSync(dir, { recursive: true, force: true }); }
});
