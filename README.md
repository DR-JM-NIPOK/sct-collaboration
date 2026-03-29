# tests/ — Unit and Integration Tests

## Running Tests

```bash
# Full test suite (no external dependencies required)
pip install pytest pytest-cov
pytest tests/ -v

# With coverage report
pytest tests/ -v --cov=sct_core --cov=likelihoods --cov-report=term-missing

# Single test class
pytest tests/test_sct_core.py::TestSoundHorizon -v
```

## Test Modules

| File | Tests | Coverage |
|---|---|---|
| `test_sct_core.py` | 28 tests | `sct_core.py` — all functions |
| `test_likelihoods.py` | 18 tests | `likelihoods/` — all modules |

## Expected Output

```
========================= test session starts ==========================
tests/test_sct_core.py::TestRb::test_R_b0_value              PASSED
tests/test_sct_core.py::TestRb::test_R_b_z_evolution         PASSED
tests/test_sct_core.py::TestSoundSpeed::test_cs2_z0           PASSED
tests/test_sct_core.py::TestSoundSpeed::test_cs2_larger_than_lcdm PASSED
...
tests/test_sct_core.py::TestModelComparison::test_BIC_lcdm_params PASSED
tests/test_sct_core.py::TestModelComparison::test_BIC_wrong_k_caught PASSED
...
========================= 46 passed in 4.2s ===========================
```

## Notes

- All tests use mock data — no real survey data files required
- Tests marked `[SKIP]` require optional packages (CAMB, CLASS, PolyChord)
- The `test_BIC_wrong_k_caught` test deliberately demonstrates the original
  paper's BIC error (k=6 instead of k=48) for transparency
