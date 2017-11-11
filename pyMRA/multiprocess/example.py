import multiprocessing
import sys
import re

class ProcessWorker(multiprocessing.Process):
    """
    This class runs as a separate process to execute worker's commands in parallel
    Once launched, it remains running, monitoring the task queue, until "None" is sent
    """

    def __init__(self, task_q, result_q):
        multiprocessing.Process.__init__(self)
        self.task_q = task_q
        self.result_q = result_q
        return

    def run(self):
        """
        Overloaded function provided by multiprocessing.Process.  Called upon start() signal
        """
        proc_name = self.name
        print( '%s: Launched' % (proc_name))
        while True:
            next_task_list = self.task_q.get()
            if next_task is None:
                # Poison pill means shutdown
                print('%s: Exiting' % (proc_name))
                self.task_q.task_done()
                break
            next_task = next_task_list[0]
            print( '%s: %s' % (proc_name, next_task))
            args = next_task_list[1]
            kwargs = next_task_list[2]
            answer = next_task(*args, **kwargs)
            self.task_q.task_done()
            self.result_q.put(answer)
        return
# End of ProcessWorker class

class Worker(object):
    """
    Launches a child process to run commands from derived classes in separate processes,
    which sit and listen for something to do
    This base class is called by each derived worker
    """
    def __init__(self, config, index=None):
        self.config = config
        self.index = index

        # Launce the ProcessWorker for anything that has an index value
        if self.index is not None:
            self.task_q = multiprocessing.JoinableQueue()
            self.result_q = multiprocessing.Queue()

            self.process_worker = ProcessWorker(self.task_q, self.result_q)
            self.process_worker.start()
            print( "Got here")
            # Process should be running and listening for functions to execute
        return

    def enqueue_process(target):  # No self, since it is a decorator
        """
        Used to place an command target from this class object into the task_q
        NOTE: Any function decorated with this must use fetch_results() to get the
        target task's result value
        """
        def wrapper(self, *args, **kwargs):
            self.task_q.put([target, args, kwargs]) # FAIL: target is a class instance method and can't be pickled!
        return wrapper

    def fetch_results(self):
        """
        After all processes have been spawned by multiple modules, this command
        is called on each one to retreive the results of the call.
        This blocks until the execution of the item in the queue is complete
        """
        self.task_q.join()                          # Wait for it to to finish
        return self.result_q.get()                  # Return the result

    @enqueue_process
    def run_long_command(self, command):
        print( "I am running number % as process "%number, self.name )

        # In here, I will launch a subprocess to run a  long-running system command
        # p = Popen(command), etc
        # p.wait(), etc
        return 

    def close(self):
        self.task_q.put(None)
        self.task_q.join()

if __name__ == '__main__':
    config = ["some value", "something else"]
    index = 7
    workers = []
    for i in range(5):
        worker = Worker(config, index)
        worker.run_long_command("ls /")
        workers.append(worker)
    for worker in workers:
        worker.fetch_results()

    # Do more work... (this would actually be done in a distributor in another class)

    for worker in workers:
        worker.close() 
