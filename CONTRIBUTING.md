# Contributing to MICOM

Thank you for your interest in contributing to **MICOM**! We welcome contributions ranging from bug fixes and documentation improvements to new simulation algorithms and visualization features.

To maintain a high standard of code quality, numerical stability, and maintainability, please review the guidelines below before submitting issues or pull requests.

For an overview of the system architecture, mathematical principles, and component design, refer to [DESIGN.md](file:///home/cdiener/code/micom/DESIGN.md).

---

## 1. Code of Conduct

We are committed to providing a friendly, safe, and welcoming environment for all contributors. Please treat others with respect and constructive feedback throughout all discussions, issues, and code reviews.

---

## 2. Getting Started & Development Setup

### 2.1 Prerequisites & Supported Python Versions

- **Supported Python Versions**: MICOM supports all **"alive" Python versions**—specifically, all CPython versions that are currently receiving active bug fixes or security maintenance from the Python Software Foundation (following the official CPython lifecycle / NEP 29). As Python versions reach their official End of Life (EOL), support is dropped in subsequent MICOM releases. CI continuously tests across all supported alive releases (currently testing through Python **3.14** across Linux, macOS, and Windows).
- **Linear/Quadratic Programming Solvers**:
  - MICOM relies on LP and QP solvers via `optlang`.
  - While basic linear operations work with GLPK or HiGHS (`highspy`), quadratic problems (e.g., `cooperative_tradeoff` and `knockout_taxa`) benefit significantly from QP-capable solvers with interior-point methods, such as **CPLEX** or **Gurobi**.
  - We recommend installing CPLEX (`pip install cplex`) or configuring Gurobi for full functionality.

### 2.2 Setting Up a Development Environment

1. **Fork and clone the repository**:
   ```bash
   git clone https://github.com/<your-username>/micom.git
   cd micom
   git checkout -b feature/my-new-feature
   ```

2. **Create and activate a virtual environment**:
   ```bash
   python3 -m venv .venv
   source .venv/bin/activate  # On Windows: .venv\Scripts\activate
   ```

3. **Install MICOM in editable mode with development dependencies**:
   ```bash
   pip install --upgrade pip wheel
   pip install -e .
   pip install pytest pytest-cov flake8 black pydocstyle
   ```

---

## 3. Code Standards and Conventions

To preserve codebase consistency across all modules, all contributions must adhere to the following standards:

### 3.1 Code Style and Formatting

- **Formatting Tool**: [`black`](file:///home/cdiener/code/micom/pyproject.toml#L8-L28) is used for code formatting with a maximum line length of **88 characters**:
  ```bash
  black --line-length 88 micom tests
  ```
- **Linter**: [`flake8`](file:///home/cdiener/code/micom/setup.cfg#L99-L106) is enforced in continuous integration:
  ```bash
  # Check for syntax errors and undefined names
  flake8 . --count --select=E9,F63,F7,F82 --show-source --statistics
  # Check general style compliance
  flake8 . --count --max-line-length=88 --statistics
  ```

### 3.2 Docstring Format (NumPy Standard)

MICOM strictly follows the **NumPy docstring convention** (`convention = numpy` in [`setup.cfg`](file:///home/cdiener/code/micom/setup.cfg#L107-L109)). Every public function, class, and method must have a comprehensive docstring containing:

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
- **Path Coercion**: Whenever a function or method accepts file system paths, decorate it with [`@pathify`](file:///home/cdiener/code/micom/micom/types.py#L17-L63). This decorator automatically converts strings or `os.PathLike` arguments into `pathlib.Path` instances:
  ```python
  from micom.types import pathify
  from pathlib import Path
  from typing import Union

  @pathify
  def load_custom_medium(file_path: Union[str, Path]) -> pd.DataFrame:
      # file_path is guaranteed to be a Path object
      return pd.read_csv(file_path)
  ```
- **Configuration Models**: All settings in [`micom.batch.configuration`](file:///home/cdiener/code/micom/micom/batch/configuration.py) must be validated via **Pydantic v2** (`pydantic.BaseModel`) with explicit type annotations, default factories, and docstrings.

---

## 4. Architectural Rules for Contributions

Before contributing code, please ensure you respect the separation of concerns defined in [DESIGN.md](file:///home/cdiener/code/micom/DESIGN.md):

1. **Maintain the Clean Separation between APIs**:
   - [`Community`](file:///home/cdiener/code/micom/micom/community.py#L38-L1272) is a single-sample in-memory COBRA model subclassing `cobra.Model`. Do not add multi-sample batch orchestration logic directly into `Community`.
   - [`Batch`](file:///home/cdiener/code/micom/micom/batch/batch.py#L26-L723) orchestrates multi-sample cohorts, delegating single-sample builds and simulations to worker tasks. Do not manipulate reaction equations or optlang expressions directly in `Batch`.
2. **Prevent Solver Memory Leaks in Parallel Workflows**:
   - All multi-sample parallel routines must use [`workflow`](file:///home/cdiener/code/micom/micom/batch/core.py#L14-L67), which enforces `multiprocessing.get_context("spawn").Pool(maxtasksperchild=1)`. Never use shared-memory multiprocessing or unconstrained threads that pool `optlang` solver pointers across tasks.
3. **Use Standardized Result Containers**:
   - Batch simulation routines must return [`GrowthResults`](file:///home/cdiener/code/micom/micom/batch/results.py#L20-L116) (`growth_rates`, `exchanges`, `annotations`).
4. **Preserve Validation Invariants**:
   - Input taxonomy tables must pass [`check_taxonomy(df)`](file:///home/cdiener/code/micom/micom/types.py#L103-L157).
   - Input media tables must pass [`check_medium(df)`](file:///home/cdiener/code/micom/micom/types.py#L65-L100).

---

## 5. Testing Guidelines

MICOM has extensive test coverage across Linux, macOS, and Windows. Every new feature or bug fix must include corresponding tests.

### 5.1 Running the Test Suite

Run tests locally using `pytest`:
```bash
# Run all tests
pytest tests/

# Run tests with coverage report
pytest --cov=micom --cov-report=term-missing tests/

# Run a specific test module
pytest tests/test_batch.py
```

### 5.2 Test Placement

- **Single Model & Optimization Tests**: Add to [`tests/test_community.py`](file:///home/cdiener/code/micom/tests/test_community.py), [`tests/test_tradeoff.py`](file:///home/cdiener/code/micom/tests/test_tradeoff.py), [`tests/test_optimizations.py`](file:///home/cdiener/code/micom/tests/test_optimizations.py).
- **Batch API & Configuration Tests**: Add to [`tests/test_batch.py`](file:///home/cdiener/code/micom/tests/test_batch.py), [`tests/test_config.py`](file:///home/cdiener/code/micom/tests/test_config.py), [`tests/test_results.py`](file:///home/cdiener/code/micom/tests/test_results.py).
- **Media & Database Tests**: Add to [`tests/test_media.py`](file:///home/cdiener/code/micom/tests/test_media.py), [`tests/test_db_media.py`](file:///home/cdiener/code/micom/tests/test_db_media.py), [`tests/test_db.py`](file:///home/cdiener/code/micom/tests/test_db.py).
- **Host & Interaction Tests**: Add to [`tests/test_host.py`](file:///home/cdiener/code/micom/tests/test_host.py), [`tests/test_interaction.py`](file:///home/cdiener/code/micom/tests/test_interaction.py).

### 5.3 Numerical Tolerances in Tests

When writing assertions for flux rates or growth rates, always use `pytest.approx` or compare against `atol`/`rtol` tolerances (e.g. `assert growth_rate == pytest.approx(0.5, abs=1e-5)`), as floating-point results may vary slightly across solvers and architectures.

---

## 6. Documentation

Documentation is generated using **Sphinx**, **Furo**, and **nbsphinx** from markdown, reStructuredText, and Jupyter notebooks in [`docs/source/`](file:///home/cdiener/code/micom/docs/source/).

### Building Documentation Locally

```bash
pip install "sphinx>=6.0" "nbsphinx>=0.9.0" furo sphinx-autoapi recommonmark
cd docs
make html
```
The compiled HTML will be placed in `docs/_build/html/index.html` (or `docs/index.html`).

When adding a new feature or tutorial:
- Provide clear explanatory prose and code examples in a corresponding Jupyter notebook under `docs/source/`.
- Ensure notebooks are executed and can run end-to-end without unhandled errors.

---

## 7. Submitting a Pull Request (PR)

1. **Sync with Main**:
   Ensure your branch is up to date with the latest `main` branch:
   ```bash
   git fetch origin
   git rebase origin/main
   ```
2. **Pre-PR Quality Checklist**:
   - [ ] All tests pass locally (`pytest tests/`).
   - [ ] Code is formatted with `black --line-length 88`.
   - [ ] Linter passes with no errors (`flake8 .`).
   - [ ] New public functions and classes have NumPy-formatted docstrings and type hints.
   - [ ] If user-facing changes were made, add a brief bullet to [`NEWS.md`](file:///home/cdiener/code/micom/NEWS.md).
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
