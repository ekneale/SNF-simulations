# Development guide

This guide will outline the basic process for developing the project locally, and the recommended tools and workflows to use when making changes to the codebase.

## Development setup

:::{important}
If you only want to use the package, the simplest option is to install it from PyPI ([https://pypi.org/project/snf-simulations/](https://pypi.org/project/snf-simulations/)).

The setup steps below are intended for anyone wanting to make new additions or changes to the codebase.
:::

### Recommended: uv

The recommended development setup uses [uv](https://docs.astral.sh/uv/). You'll need to install uv first, then you can use it to create a local development environment with all the necessary dependencies.

If you're unfamiliar with uv, the `uv sync` command creates or updates a local virtual environment and installs exactly what dependencies are declared for the project, so it effectively combines `venv` and `pip`. If you need to add a new dependency, use `uv add <package>` for runtime dependencies, or `uv add --group dev <package>` for development-only tools. The generated `uv.lock` file records the fully resolved dependency versions so all contributors and CI use the same pinned set of packages for reproducible environments.

To set up the development environment with `uv`:

```bash
git clone https://github.com/ekneale/SNF-simulations.git
cd SNF-simulations
uv sync --group dev
```

Using the `dev` group ensures you have all the tools needed for development, including linters, type checkers, testing frameworks, and documentation tools. You can see the list of packages included in the `dev` group in the `pyproject.toml` file.

### Alternative: venv + pip

You can still install the package for development using the traditional tools. Creating a virtual environment is highly recommended to isolate the development environment from the system Python. Note the `pip` command installs the package in editable mode, again with the `dev` dependency group to include development-only dependencies.

```bash
git clone https://github.com/ekneale/SNF-simulations.git
cd SNF-simulations
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -e .[dev]
```

## Repository layout

This project has the following structure:

- `src/snf_simulations/`: installable package code.
- `src/snf_simulations/data/`: packaged data files used at runtime.
- `src/snf_simulations/scripts/`: command-line scripts (accessed through entry points defined in `pyproject.toml`).
- `tests/`: unit and integration tests.
- `docs/`: Sphinx + MyST documentation source.
- `.github/workflows/`: GitHub Actions CI workflows.

## Development tools

### Ruff (linting and formatting)

[Ruff](https://docs.astral.sh/ruff/) is configured in `pyproject.toml` and used for style and static lint checks. Ruff will be run automatically on every commit (using `prek`) and pull request (using GitHub Actions), but you can also run it locally to check for any issues before pushing your changes.

If you have the development environment set up with uv, you can run:
```bash
uv run ruff check
```
to check for any linting issues in the source code, and
```bash
uv run ruff format
```
to automatically fix any formatting issues. Using the `uv run` prefix ensures you're using the version of Ruff installed in the development environment, which is important for consistency across contributors and CI.

If not, you can run `ruff` directly from the command line (the same is true for the other tools mentioned below), for example:

```bash
ruff check
```

### ty (type checking)

Type checking is handled with [ty](https://github.com/getdecoded/ty).

With uv you can run the type checker with:
```bash
uv run ty check
```

:::{important}
Ty is still in early development, which is why it is not included in the automatic pre-commit or CI checks, but it is still recommended to run it locally when developing new classes or functions.
:::


### pytest (testing)

Tests are run with [pytest](https://pytest.org), with optional coverage reporting using [pytest-cov](https://github.com/pytest-dev/pytest-cov).

With uv, you can run the tests with:
```bash
uv run pytest tests/
```
This will run all the tests in the `tests/` directory, checking that they all complete correctly and flagging any failures. You can also specify a particular test file or test function to run, for example:
```bash
uv run pytest tests/test_example.py::test_function
```
This can be useful for running a smaller subset of tests while developing new code or debugging specific issues.

In order to include code coverage information, you can run:
```bash
uv run pytest tests/ --cov=snf_simulations
```
This will generate a coverage report showing which lines of code are covered by the tests, which can help identify any gaps in test coverage that may need to be addressed.


### Prek (pre-commit checks)

[prek](https://prek.j178.dev) automatically runs various checks before changes are committed to Git, which which means only compliant code will be committed.

The included hooks are defined in `prek.toml`. Prek comes with several built-in checks, for things like ensuring correct file endings and preventing large files from being committed (see [https://prek.j178.dev/builtin/](https://prek.j178.dev/builtin/)). The current config also includes a hook for `ruff check` and `ruff format` for linting and formatting (see [https://github.com/astral-sh/ruff-pre-commit](https://github.com/astral-sh/ruff-pre-commit)). Currently type checking with `ty` is not yet included, as there is no pre-commit hook available (see [https://github.com/astral-sh/ty/issues/269](https://github.com/astral-sh/ty/issues/269)).

As the test functions using `pytest` can be time-consuming to run, following standard practice they are also not included in the pre-commit hooks. It is recommended you run pytest manually before pushing any changes, although they will be run automatically by the CI pipeline on GitHub.

To run the pre-commit checks manually, you can use the `prek` command directly or through `uv`. By default `prek run` will run the checks on all staged files, but you can also include the `--all-files` flag to run the checks on all files in the repository:

```bash
uv run prek run
```

In practice, prek is best installed as a pre-commit hook into the local Git repository, so that it runs automatically whenever you commit changes. To set up the pre-commit hook, you can run:

```bash
uv run prek install
```

Once this is installed, the checks will run automatically whenever you commit changes with Git. For any of the formatting hooks, attempted fixes will be automatically staged for commit, so you can just review the changes and commit as normal. If any checks fail, the commit will be blocked until the issues are resolved.


### GitHub Actions (continuous integration)

Source control, issues, pull requests, and releases are managed on GitHub, as described in [the contributing guide](contributing.md).

The project currently uses [GitHub Actions](https://docs.github.com/en/actions) to automate workflows for linting, testing, documentation, and publishing. The workflow files are located in `.github/workflows/`. The main continuous integration (CI) workflows run on every pull request or change to the codebase, and merging PRs will be blocked if any checks fail.

The outcomes of all the triggered actions can be viewed in the "Actions" tab of the GitHub repository [https://github.com/ekneale/SNF-simulations/actions](https://github.com/ekneale/SNF-simulations/actions), where you can see the logs for each workflow run. When a pull request is opened, you can also see the status of the checks directly on the PR page, and click through to view the logs for any failed checks.

A notable action that is not defined in the repository is `pages-build-deployment`, which handles the automatic documentation build and deployment to GitHub Pages (see below). Whenever changes are merged into the main branch, the documentation will be automatically rebuilt and published to [https://ekneale.github.io/SNF-simulations/](https://ekneale.github.io/SNF-simulations/), so any changes to the codebase that affect the documentation should also be reflected in the live documentation site.


### Sphinx (documentation)

Documentation is built with [Sphinx](https://www.sphinx-doc.org/) and written in Markdown using [MyST](https://myst-parser.readthedocs.io/).

To build the documentation locally, you can run:

```bash
cd docs
make html
```

The generated HTML files will be located in `docs/_build/html/`, and you can open `index.html` in a web browser to view the documentation.

[The API documentation](../apidocs/index.rst) is automatically generated from docstrings using Sphinx's [autodoc2](https://sphinx-autodoc2.readthedocs.io/) extension, so any new functions or classes you add should be properly documented with docstrings to ensure they are included in the API reference pages.

As described above, when changes are made to the source code on GitHub the documentation will be automatically rebuilt and published to a site hosted on GitHub Pages: [https://ekneale.github.io/SNF-simulations/](https://ekneale.github.io/SNF-simulations/). This ensures that the live documentation site is always up to date with the latest code changes.

The documentation pages use the [PyData Sphinx Theme](https://pydata-sphinx-theme.readthedocs.io/).


## Releasing new versions

:::{important}
This is only relevant for maintainers of the project.
:::

Package releases are created through the [GitHub Releases](https://docs.github.com/en/repositories/releasing-projects-on-github/about-releases) system. New releases are published to PyPI automatically by the [GitHub Actions](https://docs.github.com/en/actions) `publish.yml` workflow.

Versions are marked as tags created through the GitHub Release flow (rather than creating tags separately), so release metadata, tags, and automated publishing stay aligned. The local Python package versioning is managed through [setuptools-scm](https://github.com/pypa/setuptools_scm).
