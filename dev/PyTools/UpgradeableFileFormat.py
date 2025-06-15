import struct
from pathlib import Path
from typing import Dict, List, Union

class UpgradeableFileParser:
    _DTYPE_MAP = {
        'uint8': '<B', 'int8': '<b',
        'uint16': '<H', 'int16': '<h',
        'uint32': '<I', 'int32': '<i',
        'uint64': '<Q', 'int64': '<q',
        'float32': '<f', 'float64': '<d',
    }

    def __init__(self, path: Path):
        self._sections: Dict[str, bytes] = {}
        with open(path, 'rb') as f:
            num_sections, = struct.unpack('<I', f.read(4))
            for _ in range(num_sections):
                name_len, = struct.unpack('<I', f.read(4))
                name = f.read(name_len).decode('utf-8')
                data_size, = struct.unpack('<Q', f.read(8))
                self._sections[name] = f.read(data_size)

    def list_section_names(self) -> List[str]:
        return list(self._sections.keys())

    def get_section(self, name: str, dtype: str) -> List[Union[int, float]]:
        """
        Retrieve section data as a list of values of the given dtype.
        dtype must be one of:
          'uint8', 'int8', 'uint16', 'int16',
          'uint32', 'int32', 'uint64', 'int64',
          'float32', 'float64'
        """
        fmt = self._DTYPE_MAP.get(dtype)
        if fmt is None:
            raise ValueError(f"Unsupported dtype '{dtype}'")
        data = self._sections.get(name)
        if data is None:
            raise KeyError(f"Section '{name}' not found")
        elem_size = struct.calcsize(fmt)
        if len(data) % elem_size:
            raise ValueError(
                f"Section '{name}' size ({len(data)}) is not a multiple of element size {elem_size}"
            )
        return [v[0] for v in struct.iter_unpack(fmt, data)]
