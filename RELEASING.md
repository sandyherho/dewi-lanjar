# Publishing dewi-lanjar 0.0.1: GitHub, Read the Docs and PyPI

## 1. Check locally

```bash
cd dewi-lanjar
python -m venv .venv && source .venv/bin/activate      # Windows: .venv\Scripts\activate
pip install -e ".[test,docs]"
pytest
sphinx-build -b html -W --keep-going docs docs/_build/html   # open docs/_build/html/index.html
```

## 2. Create the GitHub repository

1. On github.com click **New repository**.
2. Owner `sandyherho`, name `dewi-lanjar`, Public.
3. Leave "Add a README", ".gitignore" and "license" **unchecked**; they are already in the folder.
4. Click **Create repository**.

## 3. Push the code

```bash
cd dewi-lanjar
git init -b main
git add .
git commit -m "dewi-lanjar 0.0.1"
git remote add origin https://github.com/sandyherho/dewi-lanjar.git
git push -u origin main
```

Optionally tag the release:

```bash
git tag -a v0.0.1 -m "dewi-lanjar 0.0.1"
git push origin v0.0.1
```

## 4. Connect Read the Docs

1. Go to <https://readthedocs.org> and choose **Sign up** (or **Log in**) **with GitHub**.
   Authorize Read the Docs to read your repositories.
2. On the dashboard click **Add project**.
3. Type `dewi-lanjar`, select `sandyherho/dewi-lanjar`, and click **Continue**.
   If Read the Docs asks you to install its GitHub App, install it for the
   `dewi-lanjar` repository and come back.
4. Check the pre-filled fields: name `dewi-lanjar` (this gives
   `dewi-lanjar.readthedocs.io`), default branch `main`, language English.
   Click **Next**.
5. Read the Docs looks for `.readthedocs.yaml`. It is already in the repository
   root, so click **This file exists**.
6. The first build starts on its own. Open **Builds** to follow the log; it
   installs `docs/requirements.txt`, installs the library with pip, and runs
   Sphinx with `docs/conf.py`. It takes about two minutes.
7. When the build is green, click **View docs**.

After this, every push to `main` rebuilds the documentation automatically.

### If the build fails

* **Project name taken**: choose another slug (for example `dewi-lanjar-gfd`)
  and update the documentation URL in `README.md`, `pyproject.toml` and
  `CITATION.cff`.
* **No build after a push**: in the project's **Settings → Integrations**,
  check that the GitHub integration exists. Re-sync it, or create one if it is
  missing.
* **Import or Sphinx errors**: the log shows the failing command. Reproduce it
  locally with the `sphinx-build` command from step 1.

### Optional settings

* **Settings → Pull request builds**: on, so each pull request gets a preview.
* **Versions**: activate the `v0.0.1` tag so a frozen copy of the docs stays
  online next to `latest`.

## 5. Publish on PyPI (optional)

```bash
pip install build twine
python -m build                  # creates dist/dewi_lanjar-0.0.1.tar.gz and .whl
twine check dist/*
twine upload dist/*              # username __token__, password: your PyPI API token
```

`pip install dewi-lanjar` then works anywhere. A version number on PyPI can
never be reused, so upload 0.0.1 only after the GitHub and docs checks pass.
