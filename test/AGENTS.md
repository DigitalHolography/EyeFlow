# Test policy

Use pytest to protect meaningful scientific, schema, scheduling, and
user-visible contracts.

## Default: add no test

Use the fewest tests that uniquely protect a change. Before adding one:

1. Search for an existing test of the same contract or failure mode.
2. Extend an existing test or fixture when it can cover the risk clearly.
3. Add a test only when it catches a distinct regression that existing tests
   cannot catch.
4. State that unique risk in the test name or change summary.

Crossing another module does not by itself justify another test.

## No redundant coverage

- Test a rule at the lowest layer that owns it.
- Add an integration test only when the wiring or boundary is itself at risk.
- Do not repeat the same happy path at unit, pipeline, UI, and end-to-end levels.
- Do not split tests for different constants. Parameterize only when each case
  represents a distinct branch or contract.
- Do not test private steps, library behavior, trivial forwarding, or every
  field of an otherwise identical mapping.
- Do not generate production-scale artifact families in a test. Mock artifact
  export unless the artifact itself is the contract; an artifact test should
  write only the minimum representative files and must not repeat full PNG/EPS
  output matrices across parameter cases.
- Merge or remove an older test when a new test subsumes it.
- Keep numerical reference cases only when each protects a distinct scientific
  invariant.

## Verification scope

- Documentation or instructions: no test.
- Name, path, layout, PDF, or UI: one focused owner/presentation test.
- Formula, axes, NaN, sign, or units: test the numerical owner; test a consumer
  only if its interpretation changed.
- Dependency or option: one DAG test; test the runner only if runtime behavior
  differs.
- Cross-workflow behavior: integrate only when the boundary itself changed.
- Run the smallest relevant test or node first. Expand only after a failure or
  when a concrete risk remains uncovered.
- No individual test should take longer than one minute. Treat a slower test as
  a scoping or performance defect instead of accepting it as normal coverage.
- Do not run the full suite unless requested or focused coverage cannot verify
  the behavior.

Common commands:

```powershell
python -m pytest test/test_file.py::test_name
python -m pytest test/test_file.py
python -m pytest --collect-only -q
```
