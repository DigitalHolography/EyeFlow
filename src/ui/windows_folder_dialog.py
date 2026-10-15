"""Windows folder multi-selection using the system IFileOpenDialog.

Tk's directory dialog has no multi-select option. Keep the COM boundary here
so controllers still receive ordinary filesystem paths, without a new dependency.
Interface slots follow Microsoft's ShObjIdl_core.h declarations.
"""

from __future__ import annotations

import ctypes
from uuid import UUID

_DWORD = ctypes.c_uint32
_HRESULT = ctypes.c_int32
_POINTER = ctypes.c_void_p
_OUT_POINTER = ctypes.POINTER(_POINTER)
_CANCELLED = 0x800704C7
_CHANGED_APARTMENT = 0x80010106
_FOLDER_OPTIONS = 0x20 | 0x40 | 0x200 | 0x8  # folders, filesystem, multi-select, keep cwd
_EXISTING_PATH_OPTIONS = 0x800 | 0x1000
_SIGDN_FILESYSPATH = 0x80058000


class _GUID(ctypes.Structure):
    _fields_ = [
        ("data1", ctypes.c_uint32),
        ("data2", ctypes.c_uint16),
        ("data3", ctypes.c_uint16),
        ("data4", ctypes.c_uint8 * 8),
    ]

    @classmethod
    def parse(cls, value: str):
        return cls.from_buffer_copy(UUID(value).bytes_le)


_CLSID_FILE_OPEN_DIALOG = _GUID.parse("dc1c5a9c-e88a-4dde-a5a1-60f82a20aef7")
_IID_FILE_OPEN_DIALOG = _GUID.parse("d57c7288-d4ad-4768-be02-9d969532d960")
_IID_SHELL_ITEM = _GUID.parse("43826d1e-e718-42ee-bc55-a1e261c37bfe")


def _check_hresult(result: int) -> None:
    if result < 0:
        raise ctypes.WinError(result)


class _ComPointer:
    def __init__(self) -> None:
        self.pointer = _POINTER()

    def __enter__(self):
        return self

    def __exit__(self, *_exc) -> None:
        if self.pointer.value:
            self.call(2, (), check=False)  # IUnknown::Release
            self.pointer.value = None

    def call(self, slot: int, argtypes: tuple, *args, check: bool = True) -> int:
        vtable = ctypes.cast(self.pointer, ctypes.POINTER(ctypes.POINTER(_POINTER))).contents
        method = ctypes.WINFUNCTYPE(_HRESULT, _POINTER, *argtypes)(vtable[slot])
        result = method(self.pointer, *args)
        if check:
            _check_hresult(result)
        return result


def _load_libraries():
    ole32 = ctypes.WinDLL("ole32")
    shell32 = ctypes.WinDLL("shell32")
    ole32.CoInitializeEx.argtypes = [_POINTER, _DWORD]
    ole32.CoInitializeEx.restype = _HRESULT
    ole32.CoUninitialize.argtypes = []
    ole32.CoUninitialize.restype = None
    ole32.CoCreateInstance.argtypes = [
        ctypes.POINTER(_GUID), _POINTER, _DWORD, ctypes.POINTER(_GUID), _OUT_POINTER,
    ]
    ole32.CoCreateInstance.restype = _HRESULT
    ole32.CoTaskMemFree.argtypes = [_POINTER]
    ole32.CoTaskMemFree.restype = None
    shell32.SHCreateItemFromParsingName.argtypes = [
        ctypes.c_wchar_p, _POINTER, ctypes.POINTER(_GUID), _OUT_POINTER,
    ]
    shell32.SHCreateItemFromParsingName.restype = _HRESULT
    return ole32, shell32


def _owner_handle(parent) -> int | None:
    if parent is None:
        return None
    parent.update_idletasks()
    user32 = ctypes.WinDLL("user32")
    user32.GetAncestor.argtypes = [_POINTER, _DWORD]
    user32.GetAncestor.restype = _POINTER
    handle = parent.winfo_id()
    return user32.GetAncestor(handle, 2) or handle  # GA_ROOT: own the native top-level window.


def askdirectories(
    *, initialdir: str = "", title: str = "Select folders", parent=None, mustexist: bool = True,
) -> tuple[str, ...]:
    """Show the native multi-folder picker; cancellation returns an empty tuple."""
    ole32, shell32 = _load_libraries()
    initialized = ole32.CoInitializeEx(None, 2 | 4)  # STA, disable legacy OLE DDE
    if (initialized & 0xFFFFFFFF) != _CHANGED_APARTMENT:
        _check_hresult(initialized)
    try:
        with _ComPointer() as dialog:
            _check_hresult(ole32.CoCreateInstance(
                ctypes.byref(_CLSID_FILE_OPEN_DIALOG), None, 1,
                ctypes.byref(_IID_FILE_OPEN_DIALOG), ctypes.byref(dialog.pointer),
            ))
            options = _DWORD()
            dialog.call(10, (ctypes.POINTER(_DWORD),), ctypes.byref(options))  # GetOptions
            flags = (options.value & ~_EXISTING_PATH_OPTIONS) | _FOLDER_OPTIONS
            if mustexist:
                flags |= _EXISTING_PATH_OPTIONS
            dialog.call(9, (_DWORD,), flags)  # SetOptions
            dialog.call(17, (ctypes.c_wchar_p,), title)  # SetTitle
            if initialdir:
                with _ComPointer() as folder:
                    status = shell32.SHCreateItemFromParsingName(
                        str(initialdir), None, ctypes.byref(_IID_SHELL_ITEM),
                        ctypes.byref(folder.pointer),
                    )
                    # A removed/offline starting directory must not prevent browsing.
                    if status >= 0:
                        dialog.call(12, (_POINTER,), folder.pointer)  # SetFolder
            status = dialog.call(3, (_POINTER,), _owner_handle(parent), check=False)  # Show
            if (status & 0xFFFFFFFF) == _CANCELLED:
                return ()
            _check_hresult(status)
            return _selected_paths(dialog, ole32)
    finally:
        # Balance successful S_OK/S_FALSE initializations, not a changed apartment.
        if initialized >= 0:
            ole32.CoUninitialize()


def _selected_paths(dialog: _ComPointer, ole32) -> tuple[str, ...]:
    paths: list[str] = []
    with _ComPointer() as items:
        dialog.call(27, (_OUT_POINTER,), ctypes.byref(items.pointer))  # GetResults
        count = _DWORD()
        items.call(7, (ctypes.POINTER(_DWORD),), ctypes.byref(count))  # GetCount
        for index in range(count.value):
            with _ComPointer() as item:
                items.call(8, (_DWORD, _OUT_POINTER), index, ctypes.byref(item.pointer))
                name = _POINTER()
                try:
                    item.call(5, (_DWORD, _OUT_POINTER), _SIGDN_FILESYSPATH, ctypes.byref(name))
                    paths.append(ctypes.wstring_at(name))
                finally:
                    if name.value:
                        ole32.CoTaskMemFree(name)
    return tuple(paths)
