import numpy as np
from src._trajectory import _Trajectory

class Trajectory_Disk(_Trajectory):
    def __init__(self, sim_params: dict=None, directory='.trajectories', mode='r', 
                 clear_first=False, trajectory_dir=None):
        
        super().__init__(sim_params=sim_params, directory=directory, mode=mode, clear_first=clear_first, trajectory_dir=trajectory_dir)
        self._load_steps()



    def _load_steps(self):
        """Load list of available steps from directory"""
        self.steps = sorted(
            int(p.stem.split('_')[-1])
            for p in self.trajectory_dir.glob('frame_*.npy')
        )

    def _read(self, frame_index: int):
        step = self.steps[frame_index]
        frame_path = self.trajectory_dir / f'frame_{step}.npy'
        return np.load(frame_path)

    def _write(self, positions: np.ndarray, step: int):
        frame_path = self.trajectory_dir / f'frame_{step}.npy'
        np.save(frame_path, positions)

    
    def _close(self):
        pass
    
    # def __enter__(self):
    #     return super().__enter__()
    
    # def __exit__(self, type, value, traceback):
    #     return super().__exit__(type, value, traceback)