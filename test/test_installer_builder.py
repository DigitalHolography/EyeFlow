"""Installer configuration and safe cleanup contracts."""
from pathlib import Path

import pytest

from installer.build_installer import clean_directory, project_metadata, write_inno_script

REPO = Path(__file__).resolve().parents[1]


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
