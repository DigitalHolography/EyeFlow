"""Installer configuration and safe cleanup contracts."""
from pathlib import Path
from subprocess import CompletedProcess
from unittest.mock import patch

import pytest

from installer.build_installer import clean_directory, project_metadata, verify_cli, write_inno_script

REPO = Path(__file__).resolve().parents[1]


def test_every_top_level_source_module_is_in_the_setuptools_manifest():
    pyproject = (REPO / "pyproject.toml").read_text(encoding="utf-8")
    setuptools_section = pyproject.split("[tool.setuptools]", 1)[1].split(
        "[tool.setuptools.packages.find]", 1
    )[0]
    source_modules = {
        path.stem
        for path in (REPO / "src").glob("*.py")
        if path.name != "__init__.py"
    }

    for module in source_modules:
        assert f'"{module}"' in setuptools_section


def test_metadata_uses_project_section_and_accepts_windows_bom(tmp_path):
    (tmp_path / "pyproject.toml").write_text(
        '[unrelated]\nname = "wrong"\nversion = "0"\n'
        '[project]\nname = "EyeFlow"\nversion = "2.3.4"\n', encoding="utf-8-sig")
    assert project_metadata(tmp_path) == ("EyeFlow", "2.3.4")


def test_clean_refuses_repository_root_and_outside_directory(tmp_path):
    repo = tmp_path / "repo"
    repo.mkdir()
    outside = tmp_path / "outside"
    outside.mkdir()
    for target in (repo, outside, repo / ".." / "outside"):
        with pytest.raises(ValueError, match="Refusing to clean"):
            clean_directory(repo, target)
        assert repo.is_dir() and outside.is_dir()


def test_clean_removes_generated_output_but_preserves_build_environment(tmp_path):
    repo = tmp_path
    stage = repo / "build/installer"
    stage.mkdir(parents=True)
    (stage / "old.exe").write_bytes(b"old")
    environment = repo / "build/installer-env"
    environment.mkdir()
    clean_directory(repo, stage)
    assert not stage.exists()
    assert environment.is_dir()


def test_inno_configuration_preserves_install_location_and_artifact_name(tmp_path):
    script = write_inno_script(REPO, tmp_path, tmp_path / "bundle", tmp_path / "dist", "EyeFlow", "2.3.4")
    text = script.read_text(encoding="utf-8")
    assert '#define MyAppExeName "EyeFlow 2.3.4.exe"' in text
    assert '#define MyAppDisplayName "EyeFlow 2.3.4"' in text
    assert '#define MyVersionDirName "v2.3.4"' in text
    assert '#define MySetupBaseName "EyeFlow-Setup-2.3.4"' in text
    assert 'VersionInfoVersion=2.3.4.0' in text
    assert 'DefaultDirName={localappdata}\\Programs\\{#MyInstallRootName}\\{#MyVersionDirName}' in text
    assert 'PrivilegesRequired=lowest' in text
    assert '@BUNDLE_DIR@' not in text
    assert 'venv' not in text


def _cli_probe_result(command, **kwargs):
    if command[1:] == ["--help"]:
        return CompletedProcess(command, 0, "usage: EyeFlow --data DATA", "")
    if command[1:] == ["--data"]:
        return CompletedProcess(command, 2, "", "--data: expected one argument")
    assert not Path(command[-1]).exists()
    return CompletedProcess(command, 1, "", "Error: missing-pipelines.txt")


def test_cli_build_checks_capture_streams_and_validate_all_exit_codes(tmp_path):
    executable = tmp_path / "EyeFlow.exe"
    environment = {"PATH": "system-path"}
    with patch("installer.build_installer.subprocess.run", side_effect=_cli_probe_result) as run:
        verify_cli(executable, cwd=tmp_path, env=environment)
    assert run.call_count == 3
    for call in run.call_args_list:
        assert call.args[0][0] == str(executable)
        assert call.kwargs["cwd"] == tmp_path
        assert call.kwargs["env"] == environment
        assert call.kwargs["capture_output"] and call.kwargs["text"]
        assert call.kwargs["timeout"] == 120


@pytest.mark.parametrize("regression", ["missing_stdout", "missing_stderr", "discarded_exit_status"])
def test_cli_build_checks_reject_installer_regressions(tmp_path, regression):
    def run_probe(command, **kwargs):
        result = _cli_probe_result(command, **kwargs)
        if regression == "missing_stdout":
            result.stdout = ""
        elif regression == "missing_stderr":
            result.stderr = ""
        elif result.returncode == 1:
            result.returncode = 0
        return result

    with patch("installer.build_installer.subprocess.run", side_effect=run_probe):
        with pytest.raises(RuntimeError, match="Frozen CLI check failed"):
            verify_cli(tmp_path / "EyeFlow.exe", cwd=tmp_path, env={})
