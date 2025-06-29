import numpy as np
import json
from pathlib import Path
from src._trajectory import _Trajectory


class Trajectory_Memmap(_Trajectory):
    def __init__(self, sim_params: dict=None, directory='.trajectories', mode='r',
                 clear_first=False, trajectory_dir=None,
                 shape=None, dtype='float32'):
        """
        shape = (steps, positions.shape)
        """

        self.shape = shape
        self.dtype = np.dtype(dtype)
        self._mmap = None
        self._metadata_file = None

        super().__init__(sim_params, directory, mode, clear_first, trajectory_dir)

        self._trajectory_file = self.trajectory_dir / "trajectory.dat"
        self._metadata_file = self.trajectory_dir / "metadata.json"

        # if self.mode == 'w':
        #     if shape is None:
        #         raise ValueError("Shape must be provided in write mode.")
        #     if clear_first: 
        #         self.steps = []
        #         self.mode_memmap = 'w+'
        #         self._mmap = np.memmap(self._trajectory_file, dtype=self.dtype, 
        #                                mode=self.mode_memmap, shape=self.shape, order='F')
        #         self._save_metadata()
        #     else: 
        #         self._load_metadata()
        #         self.mode_memmap = 'r+'
        #         self._mmap = np.memmap(self._trajectory_file, dtype=self.dtype, 
        #                                mode=self.mode_memmap, shape=self.shape, order='F')
        # elif self.mode == 'r':
        #     self.mode_memmap = 'r'
        #     self._init_reader()
        # else:
        #     raise ValueError("Mode must be 'r' or 'w'.")

        # Clear folder only if file exists and clear_first=True in writing mode,
        # otherwise open in writing and reading mode
        if mode == 'w':
            if self._metadata_file.exists() and not clear_first:
                    self.mode_memmap = 'r+'
                    self._load_metadata()
            else:
                self.mode_memmap = 'w+'
                self.steps = []
                self._save_metadata()

        # Check that file exists and then load the metadata in reading mode
        elif mode == 'r':
            if not self._metadata_file.exists():
                raise FileNotFoundError(f"No metadata found for parameters: {sim_params}")
            self.mode_memmap = 'r'
            self._load_metadata()

        else:
            raise ValueError("Mode must be 'r' or 'w'.")

        # Finally initialize memmap
        self._mmap = np.memmap(self._trajectory_file, dtype=self.dtype, mode=self.mode_memmap, shape=self.shape, order='F')



    def _save_metadata(self):
        meta = {
            "shape": self.shape,
            "dtype": str(self.dtype),
            "steps": self.steps
        }
        with open(self._metadata_file, 'w') as f:
            json.dump(meta, f)

    def _load_metadata(self):
        with open(self._metadata_file, 'r') as f:
            meta = json.load(f)
        self.shape = tuple(meta["shape"])
        self.dtype = np.dtype(meta["dtype"])
        self.steps = meta["steps"]



    def _write(self, positions, step):
        frame_idx = len(self.steps)                         # mmmmmm
        if frame_idx >= self.shape[0]:
        # if step >= self.shape[0]:
            raise IndexError("Trajectory is full.")
        # self._mmap[step] = positions
        self._mmap[frame_idx] = positions
        # self._save_metadata()                               # mmmmmmmmmm

    def _read(self, frame):
        return np.array(self._mmap[frame])  # Return a copy to avoid modifying mmap directly
        # return self._mmap[frame].copy()

    def flush(self):
        if self._mmap is not None:
            self._mmap.flush()
        self._save_metadata()

    def _close(self):
        if self._mmap is not None:
            self.flush()
            self._mmap._mmap.close()  # Force closing the internal mmap object
            self._mmap = None
