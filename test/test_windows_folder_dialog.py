"""Exercise real Windows COM objects without displaying an interactive dialog."""

import ctypes
import sys
from unittest.mock import Mock

import pytest

from ui import windows_folder_dialog as native

pytestmark = pytest.mark.skipif(sys.platform != "win32", reason="Windows shell COM interfaces")


@pytest.fixture
def libraries(monkeypatch):
    ole32, shell32 = native._load_libraries()
    ole32.CoUninitialize = Mock(wraps=ole32.CoUninitialize)
    ole32.CoTaskMemFree = Mock(wraps=ole32.CoTaskMemFree)
    monkeypatch.setattr(native, "_load_libraries", lambda: (ole32, shell32))
    return ole32, shell32


@pytest.mark.parametrize("status", [native._CANCELLED, 0x80004005])
def test_real_dialog_flags_initial_folder_and_cleanup_without_showing(
    tmp_path, monkeypatch, libraries, status,
):
    ole32, _shell32 = libraries
    real_call = native._ComPointer.call
    calls = []

    def intercept(self, slot, argtypes, *args, check=True):
        calls.append(slot)
        if slot == 3:  # Only intercept Show; setup and release call the actual Windows API.
            flags = native._DWORD()
            real_call(self, 10, (ctypes.POINTER(native._DWORD),), ctypes.byref(flags))
            required = native._FOLDER_OPTIONS | native._EXISTING_PATH_OPTIONS
            assert flags.value & required == required
            return native._HRESULT(status).value
        return real_call(self, slot, argtypes, *args, check=check)

    monkeypatch.setattr(native._ComPointer, "call", intercept)
    if status == native._CANCELLED:
        assert native.askdirectories(initialdir=str(tmp_path), title="Choose folders") == ()
    else:
        with pytest.raises(OSError):
            native.askdirectories(initialdir=str(tmp_path))
    assert 12 in calls  # The remembered initial folder reached IFileDialog::SetFolder.
    assert 27 not in calls  # No GetResults after cancellation or error.
    assert calls.count(2) == 2  # Initial folder and dialog both released.
    ole32.CoUninitialize.assert_called_once_with()


def test_com_initialization_is_balanced_if_dialog_creation_fails(monkeypatch, libraries):
    ole32, _shell32 = libraries
    monkeypatch.setattr(ole32, "CoCreateInstance", lambda *_args: native._HRESULT(0x80004005).value)
    with pytest.raises(OSError):
        native.askdirectories()
    ole32.CoUninitialize.assert_called_once_with()


def _shell_array(paths, ole32, shell32):
    """Create a real multi-item array from Unicode filesystem paths for GetResults."""
    shell32.SHParseDisplayName.argtypes = [
        ctypes.c_wchar_p, native._POINTER, native._OUT_POINTER,
        native._DWORD, ctypes.POINTER(native._DWORD),
    ]
    shell32.SHParseDisplayName.restype = native._HRESULT
    shell32.SHCreateShellItemArrayFromIDLists.argtypes = [
        native._DWORD, ctypes.POINTER(native._POINTER), native._OUT_POINTER,
    ]
    shell32.SHCreateShellItemArrayFromIDLists.restype = native._HRESULT
    pidls = (native._POINTER * len(paths))()
    try:
        for index, path in enumerate(paths):
            pidl = native._POINTER()
            native._check_hresult(shell32.SHParseDisplayName(
                str(path), None, ctypes.byref(pidl), 0, None,
            ))
            pidls[index] = pidl.value
        array = native._ComPointer()
        native._check_hresult(shell32.SHCreateShellItemArrayFromIDLists(
            len(paths), pidls, ctypes.byref(array.pointer),
        ))
        return array
    finally:
        for pidl in pidls:
            if pidl:
                ole32.CoTaskMemFree(pidl)


@pytest.mark.parametrize("fail_on_second_item", [False, True])
def test_native_result_paths_and_resources_are_handled_as_a_batch(
    tmp_path, monkeypatch, libraries, fail_on_second_item,
):
    ole32, shell32 = libraries
    paths = [tmp_path / "first folder", tmp_path / "second é folder"]
    for path in paths:
        path.mkdir()
    native._check_hresult(ole32.CoInitializeEx(None, 2 | 4))
    try:
        with _shell_array(paths, ole32, shell32) as array:
            ole32.CoTaskMemFree.reset_mock()
            real_call = native._ComPointer.call
            releases = []

            def intercept(self, slot, argtypes, *args, check=True):
                if slot == 3:
                    return 0  # Pretend the user accepted the native picker.
                if slot == 27:
                    real_call(array, 1, (), check=False)  # GetResults transfers a new reference.
                    ctypes.cast(args[0], native._OUT_POINTER)[0] = array.pointer.value
                    return 0
                if slot == 8 and args[0] == 1 and fail_on_second_item:
                    raise OSError("synthetic enumeration failure")
                if slot == 2:
                    releases.append(self.pointer.value)
                return real_call(self, slot, argtypes, *args, check=check)

            with monkeypatch.context() as context:
                context.setattr(native._ComPointer, "call", intercept)
                if fail_on_second_item:
                    with pytest.raises(OSError, match="enumeration failure"):
                        native.askdirectories(initialdir=str(tmp_path))
                else:
                    assert native.askdirectories(initialdir=str(tmp_path)) == tuple(map(str, paths))
            expected_items = 1 if fail_on_second_item else 2
            assert len(releases) == expected_items + 3  # folder, results, dialog, each item
            assert ole32.CoTaskMemFree.call_count == expected_items
            ole32.CoUninitialize.assert_called_once_with()
    finally:
        ole32.CoUninitialize()
