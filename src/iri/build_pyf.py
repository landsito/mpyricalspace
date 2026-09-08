"""
build_pyf.py -- regenerate an IRI version's f2py .pyf signature file from
its own irisub.for/irifun.for/iri_batch.f90. Run by hand (not part of the
meson build) whenever a version's .for set changes enough that its .pyf
needs regenerating (e.g. irifun.for gaining read_ig_rz/readapf107 as
separate subroutines).

irisub.for has a multi-block COMMON statement (`common /const/ umr
/const1/ humr,dumr /argexp/... /iounit/konsol`) that a plain f2py -h scan
mis-parses into broken C (`error: use of undeclared identifier`) -- garbage
declarations like `real :: ut0    /block1/hmf2` plus a malformed `common
/const/ ...` echo. None of that affects the wrapper (only iri_sub/
iri_batch's actual arguments matter for the calling convention), so this
strips those lines from the generated .pyf instead of fixing the parse.

Usage:  python3 build_pyf.py iri12
    (from src/iri/, or anywhere -- writes iri<version>/iri<version>f.pyf)
"""
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def build_pyf(version):
    d = HERE / version
    modname = f"{version}f"
    pyf = d / f"{modname}.pyf"
    subprocess.run(
        [sys.executable, "-m", "numpy.f2py", "-h", str(pyf), "-m", modname,
         "irisub.for", "irifun.for", "iri_batch.f90",
         "only:", "iri_sub", "read_ig_rz", "readapf107", "iri_batch",
         "--overwrite-signature"],
        cwd=d, check=True,
    )
    lines = pyf.read_text().splitlines()
    out = []
    for ln in lines:
        if re.match(r"\s*common\s*/", ln):
            continue  # drop the malformed multi-block COMMON echo entirely
        m = re.match(r"^(\s*(?:real|integer|logical|character)\b[^:]*::\s*\w+(?:\([^)]*\))?)(\s*/.*)?\Z", ln)
        if m and m.group(2):
            ln = m.group(1)  # strip the erroneous '/block/nextvar' suffix
        out.append(ln)
    pyf.write_text("\n".join(out) + "\n")
    print(f"wrote {pyf}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(f"usage: python3 {sys.argv[0]} <version, e.g. iri12>")
    build_pyf(sys.argv[1])
