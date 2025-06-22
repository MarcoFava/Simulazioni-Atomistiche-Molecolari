import sys
import time

class Timer:

    """Timer class inspired by John Paulett's stopwatch class."""

    def __init__(self, iterations=None, output=sys.stdout):
        self.__start_cpu = None
        self.__start_wall = None
        self.output = output
        self.cpu_time = 0.0
        self.wall_time = 0.0
        self.iterations = iterations
        try:
            self._wall_time_func = MPI.Wtime
        except:
            self._wall_time_func = time.time

    def __str__(self):
        return 'timer wall time [s]: {:.2f}, cpu time [s]: {:.2f}'.format(self.wall_time, self.cpu_time)

    def __repr__(self):
        return 'timer wall time [s]: {:.2f}, cpu time [s]: {:.2f}'.format(self.wall_time, self.cpu_time)

    def start(self):
        self.__start_cpu = self.__now_cpu()
        self.__start_wall = self.__now_wall()

    def stop(self):
        if self.__start_cpu is None:
            raise ValueError("Timer not started")
        self.cpu_time += self.__now_cpu() - self.__start_cpu
        self.wall_time += self.__now_wall() - self.__start_wall
        # self.iterations += 1

    def __now_cpu(self):
        try:
            return time.clock()
        except AttributeError:
            return time.process_time()

    def __now_wall(self):
        return self._wall_time_func()

    def __enter__(self):
        self.start()
        return self

    def __exit__(self, *args):
        self.stop()
        if self.output is not None:
            print(self, file=self.output)