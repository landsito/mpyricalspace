import os

# never let the test run trigger the companion VS Code extension install
os.environ.setdefault("MPYRICALSPACE_NO_VSCODE", "1")
