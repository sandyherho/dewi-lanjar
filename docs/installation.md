# Installation

dewi-lanjar needs Python 3.10 or newer with NumPy, SciPy, Matplotlib and
Pillow (installed automatically).

```bash
pip install dewi-lanjar
```

From source:

```bash
git clone https://github.com/sandyherho/dewi-lanjar
cd dewi-lanjar
pip install .
```

For development, tests and documentation:

```bash
pip install -e ".[test,docs]"
pytest                      # full test suite
pytest -m "not slow"        # skip the long nonlinear test
sphinx-build -b html docs docs/_build/html
```

Check the installation:

```bash
dewi-lanjar --version
dewi-lanjar validate --quick
```

Rendering uses Matplotlib's non-interactive Agg backend, so no display is
needed; dewi-lanjar runs on servers and in batch jobs.
