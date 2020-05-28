# internal lib
from lib.settings import Commands

# Note this is also used (but not referenced) by the swagger API as an enum
task_list = Commands.all_commands


class TaskReg(object):
    """
    Singleton, keeps track of what task is ahead
    """
    @staticmethod 
    def get_instance():
        if TaskReg.__instance is None:
            TaskReg()
        return TaskReg.__instance

    def __init__(self):
        if TaskReg.__instance is not None:
            return
        else:
            self.tasks = {}
            TaskReg.__instance = self

    def get_up_task(self):
        return self.tasks

    def set_up_task(self, task, subtask, other={}):
        other["task"] = task
        other["subtask"] = subtask
        self.tasks = other
        return self.tasks




