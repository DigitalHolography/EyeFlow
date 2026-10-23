"""Exercise release merge handling with real Git conflicts in isolated repos."""

from __future__ import annotations

import os
from pathlib import Path
import re
import shutil
import subprocess
import textwrap

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
POWERSHELL = shutil.which("pwsh") or shutil.which("powershell")
pytestmark = pytest.mark.skipif(POWERSHELL is None, reason="PowerShell is unavailable")


def _git(repo, *args, check=True):
    return subprocess.run(
        ["git", "-C", str(repo), *args],
        capture_output=True, text=True, check=check,
    )


@pytest.mark.parametrize("scenario", [
    "clean", "pyproject", "both", "unexpected", "bad_project_version", "bad_lock_version",
])
def test_release_merge_handles_only_expected_version_conflicts(tmp_path, scenario):
    repo = tmp_path / "repo"
    _git(tmp_path, "init", "--initial-branch=main", str(repo))
    for key, value in [
        ("user.name", "Release test"), ("user.email", "release-test@example.invalid"),
        ("commit.gpgsign", "false"), ("tag.gpgsign", "false"), ("core.autocrlf", "false"),
    ]:
        _git(repo, "config", key, value)
    project = (
        '[project]\nname = "EyeFlow"\nversion = "1.23.0"\n\n'
        'description = "Coh5rt-analysis engine"\n'
    )
    lock = (
        'version = 1\n\n[[package]]\nname = "eyeflow"\nversion = "1.23.0"\n'
        'source = { editable = "." }\ndependencies = ["old"]\n'
    )
    (repo / "pyproject.toml").write_text(project, encoding="utf-8")
    (repo / "uv.lock").write_text(lock, encoding="utf-8")
    (repo / "source.py").write_text('VALUE = "base"\n', encoding="utf-8")
    _git(repo, "add", ".")
    _git(repo, "commit", "-m", "Previous release")
    _git(repo, "checkout", "-b", "dev")

    if scenario != "clean":
        project = project.replace("Coh5rt", "Cohort")
    if scenario == "both":
        lock = lock.replace("1.23.0", "1.23.1").replace('["old"]', '["new"]')
    (repo / "pyproject.toml").write_text(project, encoding="utf-8")
    (repo / "uv.lock").write_text(lock, encoding="utf-8")
    (repo / "source.py").write_text('VALUE = "dev"\n', encoding="utf-8")
    _git(repo, "add", ".")
    _git(repo, "commit", "-m", "Development changes")
    dev_sha = _git(repo, "rev-parse", "HEAD").stdout.strip()
    _git(repo, "update-ref", "refs/remotes/origin/dev", dev_sha)
    _git(repo, "tag", "v1.24.0")

    # The workflow squashes tagged dev onto main, then changes the version.
    _git(repo, "checkout", "main")
    _git(repo, "restore", f"--source={dev_sha}", "--staged", "--worktree", "--", ".")
    project_version = "1.25.0" if scenario == "bad_project_version" else "1.24.0"
    # Match the release preparation step, including its whitespace normalization.
    release_project = re.sub(
        r'(?m)^version\s*=\s*"[^"]+"\s*$', f'version = "{project_version}"', project, count=1,
    )
    release_lock = re.sub(
        r'(?m)^version = "[^"]+"$',
        'version = "1.25.0"' if scenario == "bad_lock_version" else 'version = "1.24.0"', lock,
    )
    (repo / "pyproject.toml").write_text(release_project, encoding="utf-8")
    (repo / "uv.lock").write_text(release_lock, encoding="utf-8")
    if scenario == "unexpected":
        (repo / "source.py").write_text('VALUE = "conflicting release"\n', encoding="utf-8")
    _git(repo, "add", ".")
    _git(repo, "commit", "-m", "Version v1.24.0")
    release_sha = _git(repo, "rev-parse", "HEAD").stdout.strip()

    preview = _git(repo, "merge-tree", "--write-tree", dev_sha, release_sha, check=False)
    if scenario in {"pyproject", "both"}:
        assert "CONFLICT (content): Merge conflict in pyproject.toml" in preview.stdout
    if scenario == "both":
        assert "CONFLICT (content): Merge conflict in uv.lock" in preview.stdout

    workflow = (REPO_ROOT / ".github/workflows/release.yml").read_text(encoding="utf-8")
    step = workflow.split("      - name: Atomically publish Git refs\n", 1)[1]
    step = step.split("\n      - name:", 1)[0]
    script = textwrap.dedent(step.split("        run: |\n", 1)[1])
    # Execute the real preparation/merge/validation/commit code, with no push.
    script = script.split("\ngit push --atomic origin", 1)[0]
    script_path = tmp_path / "release-merge.ps1"
    script_path.write_text(script, encoding="utf-8")
    env = {**os.environ, "TAG_NAME": "v1.24.0", "TAG_IS_ANNOTATED": "false",
           "RELEASE_COMMIT": release_sha}
    result = subprocess.run(
        [POWERSHELL, "-NoProfile", "-NonInteractive", "-ExecutionPolicy", "Bypass",
         "-File", str(script_path)],
        cwd=repo, env=env, capture_output=True, text=True,
    )
    if scenario in {"unexpected", "bad_project_version", "bad_lock_version"}:
        assert result.returncode != 0
        expected = {
            "unexpected": "Release merge has unexpected conflicts:",
            "bad_project_version": "Merged pyproject.toml does not contain EyeFlow version 1.24.0",
            "bad_lock_version": "Merged uv.lock does not contain EyeFlow version 1.24.0",
        }[scenario]
        assert expected in result.stderr
        if scenario == "unexpected":
            assert "source.py" in result.stderr
        assert _git(repo, "rev-parse", "HEAD").stdout.strip() == dev_sha
        assert not (repo / ".git/MERGE_HEAD").exists()
    else:
        assert result.returncode == 0, result.stdout + result.stderr
        assert _git(repo, "show", "-s", "--format=%P").stdout.split() == [dev_sha, release_sha]
        assert (repo / "pyproject.toml").read_text(encoding="utf-8") == release_project
        assert (repo / "uv.lock").read_text(encoding="utf-8") == release_lock
        assert (repo / "source.py").read_text(encoding="utf-8") == 'VALUE = "dev"\n'
    assert not _git(repo, "status", "--porcelain").stdout.strip()
