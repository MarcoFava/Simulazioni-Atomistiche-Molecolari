import json
import hashlib
import shutil
from pathlib import Path

class _Trajectory:

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
        else:
            raise ValueError("Mode must be 'r' or 'w'")

    def _generate_filename(self, params):
        """
        Generate a reproducible folder name from simulation parameters.
        """
        stringified = json.dumps(params, sort_keys=True)
        return hashlib.md5(stringified.encode()).hexdigest()

    def read(self, frame):
        """Return the `System` at the given `frame`"""
        assert self.mode == 'r', "Trajectory not opened in read mode."
        system = self._read(frame)
        return system

    def write(self, system, step):
        """Write the `system` at the given `step`"""
        assert self.mode == 'w', "Trajectory not opened in write mode."
        self._write(system, step)
        self.steps.append(step)


    def __len__(self):
        return len(self.steps)

    def __iter__(self):
        for i in range(len(self)):
            yield self.read(i)

    def __enter__(self):
        return self

    def __exit__(self, type=None, value=None, traceback=None):
        self._close()
