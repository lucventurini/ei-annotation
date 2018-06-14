from eiannot.workflow.abstract import AtomicOperation, EIWorfkflow
import snakemake
import tempfile
import os
from itertools import chain


class Writer(AtomicOperation):

    __name__ = "writer"

    def __init__(self, inputs=dict(), outputs=dict(), message=None,
                 touch=True, cmd=None, index=0):

        super().__init__()
        self.input = inputs
        self.output = outputs
        self.message = message
        self.touch = touch
        self.cmd = cmd
        self.__index = index

    @property
    def rulename(self):
        return "{}_{}".format(self.__name__, self.__index)

    @property
    def loader(self):
        return []


class Writer2(Writer):

    __name__ = "writer2"


class Outer(AtomicOperation):

    __name__ = "outer"

    def __init__(self, inputs=dict(), outputs=dict(), message=None,
                 touch=True, cmd=None, index=None):
        super().__init__()
        self.input = inputs
        self.output = outputs
        self.message = message
        self.touch = touch
        self.cmd = cmd
        self.__index = index

    @property
    def rulename(self):
        return "{}_{}".format(self.__name__, self.__index)

    @property
    def loader(self):
        return []


if __name__ == '__main__':

    wf = EIWorfkflow()
    tasks = []
    print("TempDir: {}".format(tempfile.gettempdir()))
    for el in range(10):
        tasks.append(Writer(inputs={"tmpdir": tempfile.gettempdir()},
                            outputs={"written": os.path.join(tempfile.gettempdir(), "{}.txt".format(el))},
                            index=el))
    targ = Outer(inputs={"written": list(chain([_.output["written"] for _ in tasks]))}, outputs={"flag": "outer.done"},
                 touch=True)
    wf.add_edges_from([(_, targ) for _ in tasks])
    newtasks = []
    for el in range(5):
        newtasks.append(Writer2(inputs=targ.output,
                                outputs={"written": os.path.join(tempfile.gettempdir(), "{}_1.txt".format(el))},
                                index=el))
    wf.add_edges_from([(targ, _) for _ in newtasks])
    targ = Outer(inputs={"written": list(chain([_.output["written"] for _ in newtasks]))}, outputs={"flag": "finished.done"},
                 index=1, touch=True)
    wf.add_edges_from([(_, targ) for _ in newtasks])

    view = wf.graph
    # import networkx as nx
    # print([_.rulename for _ in nx.algorithms.dfs_postorder_nodes(view)])

    snakename = os.path.join(tempfile.gettempdir(), "snake.smk")
    with open(snakename, "wt") as snake:
        print(wf.__str__(source=None), file=snake)

    snakemake.snakemake(snakefile=snakename)
