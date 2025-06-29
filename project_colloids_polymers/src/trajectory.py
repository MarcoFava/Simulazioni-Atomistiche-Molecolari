import os
import json
import hashlib
import shutil
from pathlib import Path
import numpy as np

class Trajectory:
    def __init__(self, sim_params: dict=None, directory='.trajectories', mode='r', clear_first=False, trajectory_dir=None):
        self.sim_params = sim_params
        self.directory = Path(directory)
        self.mode = mode
        self.filepath = self._generate_filename(sim_params)
        
        if trajectory_dir is None:
            self.trajectory_dir = self.directory / self.filepath
        else:
            self.trajectory_dir = Path(trajectory_dir)

        if mode == 'w':
            if self.trajectory_dir.exists() and clear_first:
                shutil.rmtree(self.trajectory_dir)
            self.trajectory_dir.mkdir(parents=True, exist_ok=True)
            self.steps = []
        elif mode == 'r':
            if not self.trajectory_dir.exists():
                raise FileNotFoundError(f"No trajectory found for parameters: {sim_params}")
            self._load_steps()
        else:
            raise ValueError("Mode must be 'r' or 'w'")

    def _generate_filename(self, params):
        """
        Generate a reproducible folder name from simulation parameters.
        """
        stringified = json.dumps(params, sort_keys=True)
        return hashlib.md5(stringified.encode()).hexdigest()

    def _load_steps(self):
        """Load list of available steps from directory"""
        self.steps = sorted(
            int(p.stem.split('_')[-1])
            for p in self.trajectory_dir.glob('frame_*.npy')
        )

    def write(self, positions: np.ndarray, step: int):
        assert self.mode == 'w', "Trajectory not opened in write mode."
        frame_path = self.trajectory_dir / f'frame_{step}.npy'
        np.save(frame_path, positions)
        self.steps.append(step)

    def read(self, frame_index: int):
        assert self.mode == 'r', "Trajectory not opened in read mode."
        step = self.steps[frame_index]
        frame_path = self.trajectory_dir / f'frame_{step}.npy'
        # return np.load(frame_path), step
        return np.load(frame_path)

    def __len__(self):
        return len(self.steps)

    def __iter__(self):
        for i in range(len(self)):
            yield self.read(i)

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        pass
