from eiannot.workflow import AtomicOperation, EIWorfkflow
import snakemake
import tempfile
import os
from itertools import chain


class Writer(AtomicOperation):

    __name__ = "writer"

    def __init__(self, inputs=(), outputs=(), message=None,
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


class Writer2(Writer):

    __name__ = "writer2"


class Outer(AtomicOperation):

    __name__ = "outer"

    def __init__(self, inputs=(), outputs=(), message=None,
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


if __name__ == '__main__':

    wf = EIWorfkflow()
    tasks = []
    for el in range(10):
        tasks.append(Writer(inputs=tempfile.gettempdir(),
                            outputs=os.path.join(tempfile.gettempdir(), "{}.txt".format(el)),
                            index=el))
    targ = Outer(inputs=list(chain(*[_.output for _ in tasks])), outputs="outer.done", touch=True)
    wf.add_edges_from([(_, targ) for _ in tasks])
    newtasks = []
    for el in range(5):
        newtasks.append(Writer2(inputs=targ.output,
                                outputs=os.path.join(tempfile.gettempdir(), "{}_1.txt".format(el)),
                                index=el))
    wf.add_edges_from([(targ, _) for _ in newtasks])
    targ = Outer(inputs=list(chain(*[_.output for _ in newtasks])), outputs="finished.done", index=1, touch=True)
    wf.add_edges_from([(_, targ) for _ in newtasks])

    view = wf.graph
    import networkx as nx
    print([_.rulename for _ in nx.algorithms.dfs_postorder_nodes(view)])

    snakename = os.path.join(tempfile.gettempdir(), "snake.smk")
    with open(snakename, "wt") as snake:
        print(wf.__str__(source=None), file=snake)

    # snakemake.snakemake(snakefile=snakename)
