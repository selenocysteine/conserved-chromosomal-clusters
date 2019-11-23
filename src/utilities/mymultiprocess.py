"""@authors: Chiara E. Cotroneo, Benjamin Roques"""

import traceback
from multiprocessing import cpu_count, Pool, Manager
import sys
import logging
logger = logging.getLogger(__name__)
manager = Manager()

def wrapper_worker(file, mymultiprocess):
    try:
        mymultiprocess.target(file, *mymultiprocess.args)
    except Exception as e:
        mymultiprocess.errors.append([file, mymultiprocess.target.__name__, e])
        traceback.print_exception(type(e), e, e.__traceback__)
        logger.warning("Thread crashed on file {}: {}\n".format(file, e))
        raise ValueError

class MyMultiProcess:
    def __init__(self, threads, target, input, args=list(), maxtasksperchild=1, destroy=False):
        """ Initialises the class."""
        if threads > cpu_count():
            logger.error("Error: max {} threads allowed.".format(cpu_count()))
            sys.exit("Error: max {} threads allowed.".format(cpu_count()))

        # Number of threads to initialise
        self.num_worker_threads = threads
        # Initialises and fills the queue
        self.q = list(input)
        # Target is the function to be run
        self.target = target
        # Optional arguments of the function target
        self.args = args
        # List in shared memory to memorise errors in case one+ threads crash
        self.errors = manager.list()
        # Maximum tasks to allow per each process before respawning
        self.maxtasksperchild = maxtasksperchild
        # Kill the whole pool if one process fails
        self.destroy = destroy

    def run_all_working(self):
        pool = Pool(processes=self.num_worker_threads, maxtasksperchild=self.maxtasksperchild)
        try:
            for file in self.q:
                pool.apply_async(wrapper_worker, args=(file, self))
        except Exception:
            logger.error("A file wasn't processed correctly, we need to stop here")
            pool.terminate()
        else:
            pool.close()
            pool.join()

    def run_one_may_crash(self):
        pool = Pool(processes=self.num_worker_threads, maxtasksperchild=self.maxtasksperchild)
        try:
            for file in self.q:
                pool.apply_async(wrapper_worker, args=(file, self))
        finally:
            pool.close()
            pool.join()


    def run(self):
        """Method to run the function func in multithread."""
        if self.destroy is True:
            self.run_all_working()
        else:
            self.run_one_may_crash()

        if len(self.errors) > 0:
            logger.warning("Error: the following files were not processed.")
            for err in self.errors:
                logger.warning("{}, function {} failed with error message: {}".format(err[0], err[1], err[2]))
            sys.exit()