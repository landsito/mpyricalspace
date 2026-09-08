"""
build_pyf.py -- regenerate an IGRF version's f2py .pyf signature file from
its own igrf<N>.f. Run by hand (not part of the meson build) whenever a
version's source changes enough that its .pyf needs regenerating.

igrf<N>.f declares igrf<N>syn's 10 arguments with no intent(in)/intent(out)
annotation, so a plain `f2py -h` scan produces the wrong calling convention
(all 10 args in/out, instead of
`x, y, z, f = igrf<N>syn(isv, date, itype, alt, colat, elong)`). This script
runs the scan, then adds the correct intent split directly to the generated
.pyf -- the vendored igrf<N>.f is left untouched.

Usage:  python3 build_pyf.py igrf09 9
    (from src/igrf/, or anywhere -- writes igrf09/igrf09f.pyf; the second
    argument is the generation number used in the .f's own subroutine name,
    igrf<N>syn -- not necessarily the same as the zero-padded directory name)
"""
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

IN_ARGS = {"isv", "date", "itype", "alt", "colat", "elong"}
OUT_ARGS = {"x", "y", "z", "f"}


def build_pyf(version, gen):
    d = HERE / version
    modname = f"{version}f"
    fortran_file = f"igrf{gen}.f"
    subname = f"igrf{gen}syn"
    pyf = d / f"{modname}.pyf"
    subprocess.run(
        [sys.executable, "-m", "numpy.f2py", "-h", str(pyf), "-m", modname,
         fortran_file, "only:", subname, "--overwrite-signature"],
        cwd=d, check=True,
    )

    def _annotate(m):
        typ, name = m.group(1), m.group(2)
        if name in IN_ARGS:
            return f"{typ} intent(in) :: {name}"
        if name in OUT_ARGS:
            return f"{typ} intent(out) :: {name}"
        return m.group(0)

    text = pyf.read_text()
    text = re.sub(r"^(\s*(?:integer|double precision))\s*::\s*(\w+)\s*$",
                  _annotate, text, flags=re.M)
    pyf.write_text(text)
    print(f"wrote {pyf}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(f"usage: python3 {sys.argv[0]} <version dir, e.g. igrf09> <generation number, e.g. 9>")
    build_pyf(sys.argv[1], sys.argv[2])
