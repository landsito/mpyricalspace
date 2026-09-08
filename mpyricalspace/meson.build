py3.install_sources([
  '__init__.py',
  '__main__.py',
  'Predictor.py',
  'DataManager.py',
  'models.py',
  'survey.py',
  '_grid.py',
  '_build.py',
  '_vscode.py',
],
  subdir: 'mpyricalspace'
)

# companion VS Code extension source (built into a .vsix on demand by _vscode.py)
install_data(
  ['vscode_ext/package.json', 'vscode_ext/extension.js', 'vscode_ext/README.md'],
  install_dir: py3.get_install_dir() / 'mpyricalspace' / 'vscode_ext',
)
