# Contributing to MICOM

Thank you for your interest in contributing to **MICOM**! We welcome contributions ranging from bug fixes and documentation improvements to new simulation algorithms and visualization features.

To maintain a high standard of code quality, numerical stability, and maintainability, please review the guidelines below before submitting issues or pull requests.

For an overview of the system architecture, mathematical principles, and component design, refer to [DESIGN.md](DESIGN.md).

---

## 1. Code of Conduct

We are committed to providing a friendly, safe, and welcoming environment for all contributors. Please treat others with respect and constructive feedback throughout all discussions, issues, and code reviews.

---

## 2. Getting Started & Development Setup

### 2.1 Prerequisites & Supported Python Versions

- **Supported Python Versions**: MICOM supports CPython **3.11 through 3.14**. CI tests all four versions across Linux, macOS, and Windows.
- **Linear/Quadratic Programming Solvers**:
  - MICOM relies on LP and QP solvers via `optlang`.
  - CPLEX is not a runtime dependency. It is installed by the test dependency group; its PyPI community edition is sufficient for the test suite. Users working with larger models can configure a separately licensed CPLEX or Gurobi installation.

### 2.2 Setting Up a Development Environment

1. **Fork and clone the repository**:
   ```bash
   git clone https://github.com/<your-username>/micom.git
   cd micom
   git checkout -b feature/my-new-feature
   ```

2. **Install uv** using the [official installation instructions](https://docs.astral.sh/uv/getting-started/installation/).

3. **Create the project environment and install development dependencies**:
   ```bash
  uv sync --all-groups
   ```

    This installs MICOM in editable mode along with CPLEX and the test, lint, and documentation tools.

---

## 3. Code Standards and Conventions

To preserve codebase consistency across all modules, all contributions must adhere to the following standards:

### 3.1 Code Style and Formatting

- **Formatting Tool**: Ruff formats code with a maximum line length of **88 characters**:
  ```bash
    uv run ruff format micom tests
  ```
- **Linter**: Ruff checks for Python errors and undefined names:
  ```bash
    uv run ruff check .
  ```

### 3.2 Docstring Format (NumPy Standard)

MICOM follows the **NumPy docstring convention**. Every public function, class, and method should have a comprehensive docstring containing:

- A one-line summary followed by an optional detailed description.
- `Parameters`: Each parameter typed, followed by a clear description.
- `Returns`: The return type and description (or `Nothing` / `None`).
- `Raises`: Any specific exceptions thrown.
- `Notes` or `References`: Mathematical details, publications, or algorithmic nuances where appropriate.

#### Example:
```python
def cooperative_tradeoff(
    self,
    min_growth: float = 0.0,
    fraction: float = 1.0,
    fluxes: bool = False,
    pfba: bool = False,
    atol: float = None,
    rtol: float = None,
) -> Union[CommunitySolution, pd.Series]:
    """Find the best tradeoff between community and individual growth.

    Finds the set of growth rates which maintain a particular community
    growth and spread growth across all taxa as evenly as possible.
    This is achieved by minimizing the L2 norm of individual growth rates
    subject to a minimum community growth constraint.

    Parameters
    ----------
    min_growth : float or array-like, optional
        The minimal growth rate required for each individual taxon.
    fraction : float or list of floats in [0, 1]
        The minimum percentage of the community growth rate that must be
        maintained. For instance, 0.8 preserves 80% of maximal community growth.
    fluxes : bool, optional
        Whether to calculate and return all reaction fluxes.
    pfba : bool, optional
        Whether to obtain fluxes by parsimonious FBA rather than classical FBA.
    atol : float, optional
        Absolute tolerance for the growth rates. If None, uses solver feasibility.
    rtol : float, optional
        Relative tolerance for the growth rates. If None, uses solver feasibility.

    Returns
    -------
    micom.CommunitySolution or pandas.Series of CommunitySolution
        The solution object after optimization, or a series of solutions if
        multiple fraction values were passed.

    Raises
    ------
    OptimizationError
        If no feasible solution can be found within the given tolerances.
    """
```

### 3.3 Type Hints & Path Handling

- Use standard Python type hints (`typing.Union`, `typing.Optional`, `typing.Self`, `typing.Dict`, `pathlib.Path`).
- **Path Coercion**: Whenever a function or method accepts file system paths, decorate it with [`@pathify`](micom/types.py#L17-L63). This decorator automatically converts strings or `os.PathLike` arguments into `pathlib.Path` instances:
  ```python
  from micom.types import pathify
  from pathlib import Path
  from typing import Union

  @pathify
  def load_custom_medium(file_path: Union[str, Path]) -> pd.DataFrame:
      # file_path is guaranteed to be a Path object
      return pd.read_csv(file_path)
  ```
- **Configuration Models**: All settings in [`micom.batch.configuration`](micom/batch/configuration.py) must be validated via **Pydantic v2** (`pydantic.BaseModel`) with explicit type annotations, default factories, and docstrings.

---

## 4. Architectural Rules for Contributions

Before contributing code, please ensure you respect the separation of concerns defined in [DESIGN.md](DESIGN.md):

1. **Maintain the Clean Separation between APIs**:
   - [`Community`](micom/community.py#L38-L1272) is a single-sample in-memory COBRA model subclassing `cobra.Model`. Do not add multi-sample batch orchestration logic directly into `Community`.
   - [`Batch`](micom/batch/batch.py#L26-L723) orchestrates multi-sample cohorts, delegating single-sample builds and simulations to worker tasks. Do not manipulate reaction equations or optlang expressions directly in `Batch`.
2. **Prevent Solver Memory Leaks in Parallel Workflows**:
   - All multi-sample parallel routines must use [`workflow`](micom/batch/core.py#L14-L67), which enforces `multiprocessing.get_context("spawn").Pool(maxtasksperchild=1)`. Never use shared-memory multiprocessing or unconstrained threads that pool `optlang` solver pointers across tasks.
3. **Use Standardized Result Containers**:
   - Batch simulation routines must return [`GrowthResults`](micom/batch/results.py#L20-L116) (`growth_rates`, `exchanges`, `annotations`).
4. **Preserve Validation Invariants**:
   - Input taxonomy tables must pass [`check_taxonomy(df)`](micom/types.py#L103-L157).
   - Input media tables must pass [`check_medium(df)`](micom/types.py#L65-L100).

---

## 5. Testing Guidelines

MICOM has extensive test coverage across Linux, macOS, and Windows. Every new feature or bug fix must include corresponding tests.

### 5.1 Running the Test Suite

Run tests locally using `uv run pytest`:
```bash
# Run all tests
uv run pytest tests/

# Run tests with coverage report
uv run pytest --cov=micom --cov-report=term-missing tests/

# Run a specific test module
uv run pytest tests/test_batch.py
```

### 5.2 Test Placement

- **Single Model & Optimization Tests**: Add to [`tests/test_community.py`](tests/test_community.py), [`tests/test_tradeoff.py`](tests/test_tradeoff.py), [`tests/test_optimizations.py`](tests/test_optimizations.py).
- **Batch API & Configuration Tests**: Add to [`tests/test_batch.py`](tests/test_batch.py), [`tests/test_config.py`](tests/test_config.py), [`tests/test_results.py`](tests/test_results.py).
- **Media & Database Tests**: Add to [`tests/test_media.py`](tests/test_media.py), [`tests/test_db_media.py`](tests/test_db_media.py), [`tests/test_db.py`](tests/test_db.py).
- **Host & Interaction Tests**: Add to [`tests/test_host.py`](tests/test_host.py), [`tests/test_interaction.py`](tests/test_interaction.py).

### 5.3 Numerical Tolerances in Tests

When writing assertions for flux rates or growth rates, always use `pytest.approx` or compare against `atol`/`rtol` tolerances (e.g. `assert growth_rate == pytest.approx(0.5, abs=1e-5)`), as floating-point results may vary slightly across solvers and architectures.

---

## 6. Documentation

Documentation is generated using **Sphinx**, **Furo**, and **nbsphinx** from markdown, reStructuredText, and Jupyter notebooks in [`docs/source/`](docs/source/).

### Building Documentation Locally

```bash
uv sync --all-groups
uv run make -C docs html
```
The compiled HTML will be placed in `docs/_build/html/index.html` (or `docs/index.html`).

When adding a new feature or tutorial:
- Provide clear explanatory prose and code examples in a corresponding Jupyter notebook under `docs/source/`.
- Ensure notebooks are executed and can run end-to-end without unhandled errors.

---

## 7. Versioning

The canonical package version is `project.version` in `pyproject.toml`. Bump it with uv before preparing a release:

```bash
uv version --bump patch
# Use --bump minor or --bump major instead when appropriate.
```

`micom.__version__` reads the installed distribution metadata, so no second version string needs to be updated.

## 8. Submitting a Pull Request (PR)

1. **Sync with Main**:
   Ensure your branch is up to date with the latest `main` branch:
   ```bash
   git fetch origin
   git rebase origin/main
   ```
2. **Pre-PR Quality Checklist**:
    - [ ] All tests pass locally (`uv run pytest tests/`).
    - [ ] Code is formatted (`uv run ruff format --check micom tests`).
    - [ ] Linter passes with no errors (`uv run ruff check .`).
   - [ ] New public functions and classes have NumPy-formatted docstrings and type hints.
   - [ ] If user-facing changes were made, add a brief bullet to [`NEWS.md`](NEWS.md).
3. **Open the Pull Request**:
   - Provide a clear, descriptive PR title and summary explaining *why* the change is needed and *what* was altered.
   - Link any related issues (e.g., `Fixes #123`).
   - CI workflows will automatically run the test suite across operating systems and Python versions. Verify that all status checks turn green.

---

## Questions and Discussions

If you have questions about whether an idea fits into MICOM, or need guidance on solver integration:
- Open a topic in [GitHub Discussions](https://github.com/micom-dev/micom/discussions).
- Join the cobrapy community on [Gitter](https://gitter.im/opencobra/cobrapy).
- For QIIME 2 plugin questions, visit the [QIIME 2 Community Forum](https://forum.qiime2.org/c/community-plugin-support/).
