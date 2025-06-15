from typing import Tuple
import torch

from torch.utils.data import TensorDataset, DataLoader, random_split
from PyTools.UpgradeableFileFormat import UpgradeableFileParser
import numpy as np

# -----------------------------------------------------------------------------
# Data loading
# -----------------------------------------------------------------------------

class WaterDataset(TensorDataset):
    def __init__(self, file_path: str, smallSet : bool = False):
        file = UpgradeableFileParser(file_path)

        numAtoms = file.get_section("numAtoms", 'int32')[0]
        numFrames = file.get_section("numFrames", 'int32')[0]
        print("Num frames ", numFrames, " num atoms ", numAtoms)

        data = np.array(file.get_section("trajectory", 'float32'), dtype=np.float32)

        # Sanity check
        expected_size = numFrames * numAtoms * 3
        if data.size != expected_size:
            raise ValueError(f"Expected {expected_size} floats, got {data.size}")

        # Reshape: (numFrames *  numMolecules, 9)
        data_np = data.reshape(numFrames * (numAtoms // 3), 9)

        self.data = torch.from_numpy(data_np.astype(np.float32))
        if (smallSet):
            self.data = self.data[:100000, :]

        device = torch.device("cuda")
        self.data = self.data.to(device)
        print("Data loaded, shape: {} size: {:6.4f} GB".format(
            self.data.shape, self.data.nelement() * self.data.element_size() *1e-9))


    def __len__(self) -> int:
        return len(self.data)

    def __getitem__(self, idx: int) -> torch.Tensor:
        # returns a 9-element float Tensor
        return self.data[idx]

    def Split(self):
        val_ratio = 0.2
        total = len(self)
        val_len = int(total * val_ratio)
        train_len = total - val_len

        train_ds, val_ds = random_split(self, [train_len, val_len])
        return train_ds, val_ds