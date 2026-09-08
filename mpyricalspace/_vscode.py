'''
Companion VS Code extension -- build & install helper.

The extension (mpyricalspace/vscode_ext/) tints each ``src/<model>/`` folder in
the Explorer green / red by whether its compiled Fortran extension loads.  It is
entirely optional and editor-only.

Python wheels have no post-install hook, so installation is *attempted*, once,
best-effort:

* on the first ``python -m mpyricalspace ...`` command, and
* on the first ``import mpyricalspace`` from inside a VS Code terminal.

Every failure is swallowed and appended to ``~/.cache/mpyricalspace/vscode-ext.log``;
a sentinel (``vscode-ext.tried``) stops it retrying.  Run it by hand any time:

    python -m mpyricalspace vscode              # build the .vsix + code --install-extension
    python -m mpyricalspace vscode build [DIR]  # just build the .vsix
    python -m mpyricalspace vscode status       # installed? was it tried?

Opt out completely with  MPYRICALSPACE_NO_VSCODE=1.
'''
import os
import sys

EXT_ID = "mpyricalspace.mpyricalspace-build-status"

_SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "vscode_ext")
_STATE = os.path.expanduser(os.environ.get("MPYRICALSPACE_DATA_DIR", "~/.cache/mpyricalspace"))
_LOG = os.path.join(_STATE, "vscode-ext.log")
_SENTINEL = os.path.join(_STATE, "vscode-ext.tried")


def _now():
    import datetime
    return datetime.datetime.now().isoformat(timespec="seconds")


def _log(msg):
    try:
        os.makedirs(_STATE, exist_ok=True)
        with open(_LOG, "a") as fh:
            fh.write("%s  %s\n" % (_now(), msg))
    except OSError:
        pass


def find_code():
    '''Path to the VS Code (or VSCodium) CLI, or None.'''
    import shutil
    for name in ("code", "code-insiders", "codium"):
        exe = shutil.which(name)
        if exe:
            return exe
    guesses = [
        "/Applications/Visual Studio Code.app/Contents/Resources/app/bin/code",
        "/Applications/VSCodium.app/Contents/Resources/app/bin/codium",
        "/usr/share/code/bin/code", "/usr/bin/code", "/snap/bin/code",
        os.path.expanduser("~/AppData/Local/Programs/Microsoft VS Code/bin/code.cmd"),
    ]
    return next((g for g in guesses if os.path.isfile(g)), None)


def _pkg():
    import json
    with open(os.path.join(_SRC, "package.json")) as fh:
        return json.load(fh)


def build_vsix(dest_dir):
    '''Write ``<name>-<version>.vsix`` into *dest_dir* (standard library only).  Returns its path.'''
    import zipfile
    import xml.sax.saxutils as sx

    pkg = _pkg()
    out = os.path.join(dest_dir, "%s-%s.vsix" % (pkg["name"], pkg["version"]))
    esc = lambda s: sx.escape(str(s))
    manifest = (
        '<?xml version="1.0" encoding="utf-8"?>\n'
        '<PackageManifest Version="2.0.0" xmlns="http://schemas.microsoft.com/developer/vsx-schema/2011">\n'
        '  <Metadata>\n'
        '    <Identity Language="en-US" Id="%s" Version="%s" Publisher="%s" />\n'
        '    <DisplayName>%s</DisplayName>\n'
        '    <Description xml:space="preserve">%s</Description>\n'
        '    <Categories>Other</Categories>\n'
        '    <Properties><Property Id="Microsoft.VisualStudio.Code.Engine" Value="%s" /></Properties>\n'
        '  </Metadata>\n'
        '  <Installation><InstallationTarget Id="Microsoft.VisualStudio.Code" /></Installation>\n'
        '  <Dependencies/>\n'
        '  <Assets><Asset Type="Microsoft.VisualStudio.Code.Manifest" Path="extension/package.json" Addressable="true" /></Assets>\n'
        '</PackageManifest>\n'
        % (esc(pkg["name"]), esc(pkg["version"]), esc(pkg["publisher"]),
           esc(pkg["displayName"]), esc(pkg["description"]), esc(pkg["engines"]["vscode"]))
    )
    content_types = (
        '<?xml version="1.0" encoding="utf-8"?>\n'
        '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">\n'
        '  <Default Extension=".json" ContentType="application/json" />\n'
        '  <Default Extension=".js" ContentType="application/javascript" />\n'
        '  <Default Extension=".md" ContentType="text/markdown" />\n'
        '  <Default Extension=".vsixmanifest" ContentType="text/xml" />\n'
        '</Types>\n'
    )
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("extension.vsixmanifest", manifest)
        z.writestr("[Content_Types].xml", content_types)
        for f in ("package.json", "extension.js", "README.md"):
            z.write(os.path.join(_SRC, f), "extension/" + f)
    return out


def is_installed():
    '''True if VS Code reports the extension installed.'''
    import subprocess
    code = find_code()
    if not code:
        return False
    try:
        r = subprocess.run([code, "--list-extensions"], capture_output=True, text=True, timeout=20)
        return EXT_ID in r.stdout
    except (OSError, subprocess.SubprocessError):
        return False


def install(verbose=True):
    '''Build the .vsix and ``code --install-extension`` it.  Returns True on success.'''
    import tempfile
    import subprocess

    code = find_code()
    if not code:
        _log("install: 'code' CLI not found")
        if verbose:
            print("VS Code 'code' CLI not found.  In VS Code run  Cmd/Ctrl+Shift+P -> "
                  "'Shell Command: Install code command in PATH'  and try again.")
        return False
    try:
        with tempfile.TemporaryDirectory() as td:
            vsix = build_vsix(td)
            r = subprocess.run([code, "--install-extension", vsix, "--force"],
                               capture_output=True, text=True, timeout=90)
        ok = r.returncode == 0
        tail = (r.stderr or r.stdout or "").strip().replace("\n", " ")[:300]
        _log("install: %s via %s | %s" % ("ok" if ok else "FAILED", code, tail))
        if verbose:
            msg = (r.stdout or r.stderr or "").strip()
            if msg:
                print(msg)
            if ok:
                print("installed.  Fully quit VS Code (Cmd/Ctrl+Q) and reopen the first time.")
        return ok
    except (OSError, subprocess.SubprocessError) as e:
        _log("install: exception %r" % (e,))
        if verbose:
            print("install failed: %s" % e)
        return False


def _in_vscode():
    return (os.environ.get("TERM_PROGRAM") == "vscode"
            or "VSCODE_PID" in os.environ
            or "VSCODE_GIT_IPC_HANDLE" in os.environ)


def maybe_autoinstall(background=False):
    '''First-run, best-effort, silent.  Does nothing outside a VS Code context, and
    writes a sentinel so it is attempted only once.'''
    if os.environ.get("MPYRICALSPACE_NO_VSCODE") or os.path.exists(_SENTINEL):
        return
    if not _in_vscode() or not find_code():
        return                                          # not VS Code -> touch nothing, don't burn the sentinel
    try:
        os.makedirs(_STATE, exist_ok=True)
        with open(_SENTINEL, "w") as fh:
            fh.write(_now())
    except OSError:
        return
    _log("auto: first-run attempt (%s)" % ("background" if background else "inline"))
    if background:
        import subprocess
        try:
            subprocess.Popen(
                [sys.executable, "-c",
                 "from mpyricalspace._vscode import install; install(verbose=False)"],
                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, start_new_session=True)
        except OSError:
            pass
    else:
        install(verbose=False)


def status_text():
    tried = _now_from(_SENTINEL)
    lines = ["VS Code extension  (%s)" % EXT_ID,
             "  source     : %s" % _SRC,
             "  code CLI   : %s" % (find_code() or "not found"),
             "  installed  : %s" % ("yes" if is_installed() else "no"),
             "  auto-tried : %s" % (tried or "no"),
             "  log        : %s" % (_LOG if os.path.exists(_LOG) else "(none)")]
    return "\n".join(lines)


def _now_from(path):
    try:
        with open(path) as fh:
            return fh.read().strip()
    except OSError:
        return None
