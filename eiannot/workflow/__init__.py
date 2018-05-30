import abc  # The single instances of the
import networkx as nx
from eicore.external_process.snakemake_helper import loadPreCmd


class AtomicOperation(metaclass=abc.ABCMeta):

    __name__ = "atomic"
    __loader = []

    @abc.abstractmethod
    def __init__(self):
        """Each task must have:
        - input
        - output
        - shell command / touch
        Optional:
        - message"""

        self.__outputs = dict()
        self.__inputs = dict()
        self.__message = None
        self.__touch = False
        self.__cmd = None
        self.__log = None
        self.__params = dict()
        self.__configuration = dict()

    @property
    @abc.abstractmethod
    def rulename(self):
        return self.__name__

    @property
    def message(self):
        return self.__message

    @message.setter
    def message(self, message):
        if isinstance(message, bytes):
            message = message.decode()
        elif message is not None and not isinstance(message, str):
            raise TypeError("Invalid type of message: {}".format(type(message)))
        self.__message = message

    @property
    def log(self):
        return self.__log

    @log.setter
    def log(self, log):
        if isinstance(log, bytes):
            log = log.decode()
        elif log is not None and not isinstance(log, str):
            raise TypeError("Invalid type of log: {}".format(type(log)))
        self.__log = log

    @property
    def touch(self):
        return self.__touch

    @touch.setter
    def touch(self, touch):
        if touch not in (True, False):
            raise ValueError("Invalid touch: {}".format(touch))
        self.__touch = touch

    @property
    def cmd(self):
        return self.__cmd

    @cmd.setter
    def cmd(self, cmd):
        if isinstance(cmd, bytes):
            cmd = cmd.decode()
        elif cmd is not None and not isinstance(cmd, str):
            raise TypeError("Invalid type of log: {}".format(type(cmd)))
        self.__cmd = cmd

    @property
    def input(self):
        return self.__inputs

    @input.setter
    def input(self, inputs):

        """In the EI workflow, inputs *must* be dictionaries."""

        if not isinstance(inputs, dict):
            raise TypeError("Inputs must be dictionaries")
        elif not all(isinstance(_[0], (str, bytes)) and isinstance(_[1], (str, bytes)) for _ in
                       inputs.items()):
            raise ValueError("Input dictionaries must contain only strings")
        self.__inputs = inputs

    @property
    def output(self):
        return self.__outputs

    @output.setter
    def output(self, outputs):

        if not isinstance(outputs, dict):
            raise TypeError("Inputs must be dictionaries")
        elif not all(isinstance(_[0], (str, bytes)) and isinstance(_[1], (str, bytes)) for _ in
                       outputs.items()):
            raise ValueError("Input dictionaries must contain only strings")
        self.__outputs = outputs

    def add_to_params(self, key, value):
        """This method is used to add a key/value to params.
        We do not want to expose directly the dictionary to the outside.
        Params is useful for containing information such as memory, job name, etc."""

        if not isinstance(key, str):
            if isinstance(key, bytes):
                key = key.decode()
            else:
                raise TypeError("Invalid type: {}".format(type(key)))
        if not isinstance(value, str):
            if isinstance(value, bytes):
                value = value.decode()
            else:
                raise TypeError("Invalid type: {}".format(type(key)))
        self.__params[key] = value

    @property
    def load(self):
        """This property defines how to "load" a module for use in cluster environments."""

        # TODO: this will have to be revised once we have the configuration
        runs = [self.__configuration["programs"][_]["runs"] for _ in self.__loader]

        return loadPreCmd(*runs)

    def delete_from_params(self, key):
        if key in self.__params:
            del self.__params[key]

    def __str__(self):
        """This will create the SnakeMake-like block of code."""

        d = dict()
        key = "rule {}".format(self.rulename)
        d[key] = dict()
        if self.message:
            d[key]["message"] = self.message
        if self.log:
            d[key]["log"] = self.log
        d[key]["input"] = self.input
        d[key]["output"] = self.output
        if self.cmd:
            d[key]["shell"] = self.cmd
        string = ["rule {}:".format(self.rulename)]
        if isinstance(self.input, str):
            string.append("  input: \"{}\"".format(self.input))
        else:
            string.append("  input:")
            for obj in self.input[:-1]:
                string.append("    \"{}\",".format(obj))
            string.extend(["    \"{}\"".format(self.input[-1])])
        if self.touch is True:
            if len(self.output) == 1:
                string.append("  output: touch(\"{}\")".format(self.output[0]))
            else:
                string.append("  output:")
                for obj in self.output[-1]:
                    string.append("    touch(\"{}\"),".format(obj))
                string.append("  output: touch(\"{}\")".format(self.output[-1]))
        else:
            if isinstance(self.output, str):
                string.append("  output: \"{}\"".format(self.output))
            else:
                string.append("  output:")
                for obj in self.output[:-1]:
                    string.append("    \"{}\",".format(obj))
                string.append("    \"{}\"".format(self.output[-1]))
        if self.message:
            string.append("  message: \"{}\" """.format(self.message))
        if self.log:
            string.append("  log: \"{}\"".format(self.log))
        if self.cmd:
            string.append("  shell: \"{}\"".format(self.cmd))

        return "\n".join(string) + "\n"


class EIWorfkflow:

    """The workflow is, at its core, a nx.DiGraph.
    We use the dfs_postorder_nodes to get the correct order of printing.
    The goal is to print out a completely correct SnakeMake workflow, where all the
    computation necessary to determine input/output files and how to create them has been already pre-calculated."""

    def __init__(self):

        self.__graph = nx.DiGraph()

    def check_graph(self):
        """This method will check that the graph is valid, before printing."""

        pass

    def add_node(self, node: AtomicOperation):

        if not isinstance(node, AtomicOperation):

            raise ValueError("Only AtomicOperations are valid nodes")
        self.__graph.add_node(node)

    def __add__(self, node):

        self.add_node(node)

    def add_edge(self, in_edge, out_edge):
        """Here we are going to exploit the fact that adding an edge will automatically add the nodes
        in the graph."""
        if not isinstance(in_edge, AtomicOperation) or not isinstance(out_edge, AtomicOperation):
            raise TypeError("Invalid types for adding edges - only AtomicOperations are allowed")

        self.__graph.add_edge(in_edge, out_edge)

    def add_edges_from(self, tasks):

        # First check that everything is actually an AtomicOperations
        for el in tasks:
            if len(el) != 2:
                raise TypeError("Invalid task edge")
            self.add_edge(*el)
        return

    def __str__(self, source=None):

        self.check_graph()

        snake = ""

        for task in nx.algorithms.dfs_postorder_nodes(self.__graph, source=source):
            snake += str(task)

        return snake

    def connect(self, inbound: AtomicOperation, outbound: AtomicOperation, name: [str|bytes]):

        self.add_edge(inbound, outbound)
        outbound.input[name] = inbound.output[name]

    @property
    def graph(self):
        return self.__graph.copy(as_view=True)

    def merge(self, others):

        self.__graph = nx.compose_all([self.__graph] + others)

    def __iter__(self):
        return iter(self.__graph)

    @property
    def adj(self):
        return self.__graph.adj


class EIWrapper(EIWorfkflow):

    @property
    def output(self):

        """EI wrappers must have *ONE* output"""

        nodes = [_ for _ in self if not self.adj[_]]
        if not len(nodes) == 1 and len(nodes[0].output) == 1:
            raise ValueError("EI wrappers must have at most one output!")
        return nodes[0].output
