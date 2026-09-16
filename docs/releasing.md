# Releasing MobsPy

The release workflow follows the tag-triggered workflow used by
[BioDisCo/pintrs](https://github.com/BioDisCo/pintrs/blob/main/.github/workflows/publish.yml).
MobsPy builds one portable Python wheel and a source archive.

## One-time setup

Create a GitHub Actions environment named `pypi`. In the existing MobsPy project's
PyPI **Publishing** settings, register this Trusted Publisher:

| Field | Value |
| --- | --- |
| Repository owner | `ROBACON` |
| Repository | `mobspy` |
| Workflow filename | `release.yml` |
| Environment | `pypi` |

No PyPI API token is needed. See the
[PyPI Trusted Publisher instructions](https://docs.pypi.org/trusted-publishers/adding-a-publisher/).
This configuration lives in the service accounts and is not created by committing
the workflow file.

## Publish a version

1. Update `project.version` in `pyproject.toml` and the changelog. Use a new version;
   PyPI releases cannot be replaced. For the rewrite, start with `3.0.0rc1`.
2. Merge the version commit to `main` after CI passes.
3. Create and push the matching tag:

   ```bash
   git tag -a v3.0.0rc1 -m "MobsPy 3.0.0 release candidate 1"
   git push origin v3.0.0rc1
   ```

The workflow checks that the tag targets a commit on `main` and matches package
metadata. Reusable CI checks lint, types, documentation, and the built wheel on
Python 3.11–3.14, including every tutorial notebook and MobsPy example script.
External reference models run in a separate Python 3.12 job; see
[verification](verification.md). The tested artifacts are published to PyPI, then attached to a
GitHub Release. Release candidates are marked as prereleases.

For a manual run, select the Release workflow, provide an existing tag, and enter
`RELEASE` in `confirm_publish`. The same validation and CI apply. If publishing
succeeded but GitHub Release creation failed, rerun only the failed job.

After release-candidate feedback, update the version to `3.0.0`, merge, and push
`v3.0.0`. Users can try a candidate with `pip install --pre mobspy`.

## Integrating the rewrite

At the reviewed tips (`167543d` on the rewrite and `a6f4310` on `main`),
the rewrite is 31 commits ahead and 5 behind. Main's intervening changes add
Codespell and fix spelling; they do not change simulation behavior. Recheck the
remote tips before integration.

1. Commit the rewrite fixes on the rewrite branch. Preserve a `maintenance/2.x`
   branch from the old `main` if that series needs continued fixes.
2. Merge `origin/main` into the rewrite branch. Keep the new module layout and
   rewritten notebooks when resolving edits to deleted legacy files. Carry
   spelling corrections into their replacement files. Preserve `.codespellrc`
   and `.github/workflows/codespell.yml` from main.
3. Run the full CI suite and the incoming Codespell check on the merged result.
   Review the migration guide and representative scientific model outputs in
   the integration PR. The regression suite includes analytic decay, unit
   conversions, events, sweeps, compositions, and journal examples.
4. Squash-merge the integration PR into `main`, making the rewrite one reviewable
   change that can be reverted together. Do not cherry-pick individual file moves
   into main: the compiler, DSL, units, and runtime changes depend on each other.
5. Create `v3.0.0rc1` on the merged main commit only after CI is green. Collect
   feedback before publishing `3.0.0`; update the version and changelog for every
   subsequent candidate.

Keep code integration and tag creation as separate steps so publication is explicit.
The release workflow will reject a tag whose commit is outside main's history.
