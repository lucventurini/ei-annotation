import abc  # The single instances of the
import networkx as nx
from eicore.external_process.snakemake_helper import loadPreCmd
import os
from frozendict import frozendict
import re
import copy


class Sample(metaclass=abc.ABCMeta):

    def __init__(self, label, read_dir):
        self.__label = label
        self.__read_dir = read_dir
        self.__strandedness = None
        if not os.path.exists(self.read_dir):
            os.makedirs(self.read_dir)

    @property
    def read_dir(self):
        return self.__read_dir

    def __str__(self):
        return self.label

    @property
    def label(self):
        return self.__label

    @property
    def stranded(self):
        return self.strandedness is not None and self.strandedness != "fr-unstranded"

    @property
    def strandedness(self):
        """Returns the exact strandedness of the sample."""
        # We will check the exact type of the
        return self.__strandedness


class LongSample(Sample):

    def __init__(self, readfile, label, read_dir, strandedness):

        super().__init__(label=label, read_dir=read_dir)
        suffix = readfile.split(".")[-1]
        if suffix in ("fa", "fna", "fasta"):
            suffix = ".fa"
        elif suffix in ("fq", "fastq"):
            suffix = ".fq"
        else:
            suffix = ".{}".format(suffix)
        rout = os.path.join(self.read_dir,
                            "{label}.long{suffix}".format(**locals()))

        if not os.path.islink(rout):
            os.symlink(os.path.abspath(readfile), rout)
        self.__readfile = rout
        self.__strandedness = strandedness

    @property
    def readfile(self):
        return self.__readfile

    @property
    def read1(self):
        return self.readfile

    @property
    def read2(self):
        """This will always return None; long-read samples only have one set of reads."""
        return None

    @property
    def strandedness(self):
        return self.__strandedness


class ShortSample(Sample):
    """This simple class defines the input reads for a given sample."""

    def __init__(self, read1, read2, label, read_dir, strandedness=None):

        super(ShortSample, self).__init__(label=label, read_dir=read_dir)

        self.__strandedness = strandedness
        suffix = read1.split(".")[-1]
        if suffix not in ("gz", "bz2"):
            suffix = ""
        else:
            suffix = ".{}".format(suffix)
        self.__suffix = suffix

        r1out = os.path.join(read_dir, "{label}.R1.fq{suffix}".format(**locals()))

        if not os.path.islink(r1out):
            os.symlink(os.path.abspath(read1), r1out)
        self.__read1 = r1out

        if read2 is not None:
            r2out = os.path.join(read_dir, "{label}.R2.fq{suffix}".format(**locals()))
            if not os.path.islink(r2out):
                os.symlink(os.path.abspath(read2), r2out)
            self.__read2 = r2out

    @property
    def strandedness(self):
        return self.__strandedness

    @property
    def read1(self):
        return self.__read1

    @property
    def read2(self):
        return self.__read2

    @property
    def paired(self):
        return self.read2 is not None

    @property
    def suffix(self):
        return self.__suffix

    @property
    def reader(self):
        if not self.suffix:
            return ""
        elif self.suffix == ".gz":
            return "zcat"
        elif self.suffix == ".bz2":
            return "bzcat"


class AtomicOperation(metaclass=abc.ABCMeta):

    __name__ = "atomic"

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
        self.__threads = None
        self.__configuration = dict()
        self.__temps = []
        self.__touchers = []

    def __hash__(self):
        """The hash will only consider the rule name. In a SnakeMake graph,
        a rule should be present only once anyway."""
        return hash(self.rulename)

    def __eq__(self, other):

        """In order to make sure that we add a rule only once,
        we are going to compare the hashes, ie, the rule names."""

        return hash(self) == hash(other)

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

    @property
    def temps(self):

        return self.__temps

    @temps.setter
    def temps(self, temps):
        if not isinstance(temps, (tuple, set, list)):
            raise TypeError
        if not all(_ in self.output for _ in temps):
            raise KeyError
        self.__temps = temps

    @property
    def touchers(self):

        return self.__touchers

    @touchers.setter
    def touchers(self, touchers):
        if not isinstance(touchers, (tuple, set, list)):
            raise TypeError
        if not all(_ in self.output for _ in touchers):
            raise KeyError
        self.__touchers = touchers

    @input.setter
    def input(self, inputs):

        """In the EI workflow, inputs *must* be dictionaries."""

        if not isinstance(inputs, dict):
            raise TypeError("Inputs must be dictionaries")
        elif not all(isinstance(_[0], (str, bytes)) and isinstance(_[1], (str, bytes, list)) for _ in
                       inputs.items()):
            raise ValueError("Input dictionaries must contain only strings")
        self.__inputs = copy.deepcopy(inputs)

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
        self.__outputs = copy.deepcopy(outputs)

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
    @abc.abstractmethod
    def loader(self):
        """Each atomic operation has to specify which are the modules to load for its functioning."""
        return []

    @property
    def load(self):
        """This property defines how to "load" a module for use in cluster environments."""

        # TODO: this will have to be revised once we have the configuration
        try:
            to_load = [self.configuration["programs"].get(_, dict()).get("load", '') for _ in self.loader]
        except KeyError:
            raise KeyError("{}:\n{}".format(self.rulename, self.configuration))
        cmd = loadPreCmd(*to_load)
        if cmd:
            return cmd
        else:
            return ''

    def delete_from_params(self, key):
        if key in self.__params:
            del self.__params[key]

    def __str__(self):
        """This will create the SnakeMake-like block of code."""

        d = dict()
        rulename = re.sub("\.", "_", re.sub("-", "_", self.rulename))
        string = ["rule {}:".format(rulename)]
        # Inputs now are always dictionaries
        if self.input:
            string.append("  input:")
            for key, value in self.input.items():
                if isinstance(value, bytes):
                    value = value.decode()
                if isinstance(value, list):
                    string.append("    {key}={value},".format(**locals()))
                else:
                    string.append("    {key}=\"{value}\",".format(**locals()))
            string[-1] = string[-1].rstrip(",")  # Remove trailing comma
        string.append("  output:")
        for key, value in self.output.items():
            if self.touch is True or key in self.touchers:
                assert len(self.output) == 1, (self.rulename, self.input, self.output)
                string.append("    {key}=touch(\"{value}\")".format(**locals()))
            elif isinstance(value, list):
                string.append("    {key}={value},".format(**locals()))
            elif key in self.temps:
                string.append("    {key}=temp({value})".format(**locals()))
            else:
                string.append("    {key}=\"{value}\",".format(**locals()))
        string[-1] = string[-1].rstrip(",")
        if self.message:
            string.append("  message: \"{}\" """.format(self.message))
        if self.log:
            string.append("  log: \"{}\"".format(self.log))
        if self.threads > 1:
            string.append("  threads: {}".format(self.threads))

        # Add resources
        if self.resources:
            string.append("  resources:")
            for resource in self.resources:
                string.append('    {key}={value}'.format(key=resource, value=self.resources[resource]))

        if self.cmd:
            string.append("  shell: \"{}\"".format(self.cmd))

        return "\n".join(string) + "\n"

    @property
    def configuration(self):
        return frozendict(self.__configuration)

    @configuration.setter
    def configuration(self, d):
        if not isinstance(d, (dict, frozendict)):
            raise TypeError("Invalid configuration type: {}".format(type(d)))
        self.__configuration = d

    @property
    def threads(self):
        # Single location. We can change this whenever we desire.
        if self.__threads is None:
            return self.configuration.get("threads", 1)
        else:
            return self.__threads

    @threads.setter
    def threads(self, threads):
        if not isinstance(threads, int) or threads < 1:
            raise ValueError
        self.__threads = threads

    @property
    def species(self):
        """This property returns the name of the organism to be used."""
        # TODO: check consistency of name
        return self.configuration.get("name", "eiannot")

    @property
    def genome(self):
        return os.path.join(self.configuration["outdir"], "inputs", "reference", "genome.fa")

    @property
    def masked_genome(self):
        return os.path.join(self.configuration["outdir"], "repeats", "masker", "genome.masked.fa")

    @property
    def min_intron(self):
        return max(self.configuration["reference"]["min_intron"], 20)

    @property
    def max_intron(self):
        return self.configuration["reference"]["max_intron"]

    @property
    def resources(self):

        """This property will define the resources to be used for each rule."""

        default = self.configuration["resources"]["default"]

        if self.step and self.step in self.configuration["resources"]:
            default.update(self.configuration["resources"][self.__step])

        return default

    @property
    def step(self):
        """This property is used to retrieve the resource definition from the dictionary"""
        return None


class Linker(AtomicOperation):

    """Basic class to link an input into an output, when no operation is needed
    but still the pipeline needs two different files."""

    def __init__(self, input_file, output_file, key_in, key_out, rulename):
        super().__init__()
        self.__key_in = None
        self.key_in = key_in
        self.__key_out = None
        self.key_out = key_out
        self.input[self.key_in] = input_file
        self.output[self.key_out] = output_file
        self.__rulename = rulename

    @property
    def key_in(self):
        return self.__key_in

    @key_in.setter
    def key_in(self, key):
        assert isinstance(key, str)
        self.__key_in = key

    @property
    def key_out(self):
        return self.__key_out

    @key_out.setter
    def key_out(self, key):
        assert isinstance(key, str)
        self.__key_out = key

    @property
    def cmd(self):
        outdir = os.path.dirname(self.output[self.key_out])
        input, output = self.input[self.key_in], self.output[self.key_out]
        link_src = os.path.relpath(input, start=outdir)
        link_name = os.path.basename(output)

        return "mkdir -p {outdir} && cd {outdir} && ln -s {link_src} {link_name}".format(**locals())

    @property
    def rulename(self):
        return self.__rulename

    @property
    def loader(self):
        return []


class EIWorfkflow:

    """The workflow is, at its core, a nx.DiGraph.
    We use the dfs_postorder_nodes to get the correct order of printing.
    The goal is to print out a completely correct SnakeMake workflow, where all the
    computation necessary to determine input/output files and how to create them has been already pre-calculated."""

    def __init__(self, configuration=None):

        self.__graph = nx.DiGraph()
        self.__configuration = dict()
        self.configuration = configuration

    def check_graph(self):
        """This method will check that the graph is valid, before printing."""

        pass

    def add_node(self, node: AtomicOperation):

        # TODO: we have to verify whether we should add the nodes as *strings*
        # TODO: another possibility is to have their hash function to return the string, to avoid ambiguities

        if not isinstance(node, AtomicOperation):
            raise ValueError("Only AtomicOperations are valid nodes")
        self.__graph.add_node(node)

    def __add__(self, node):

        self.add_node(node)

    def __contains__(self, item):

        if isinstance(item, AtomicOperation):
            return item in self.graph
        elif isinstance(item, (str, bytes)):
            return any(rule.rulename == item for rule in self)
        else:
            return False

    def add_edge(self, in_edge, out_edge):
        """Here we are going to exploit the fact that adding an edge will automatically add the nodes
        in the graph."""
        if isinstance(in_edge, EIWrapper):
            self.merge([in_edge])
            in_edge = in_edge.exit
        if in_edge == out_edge:  # Silently avoid linking with oneself
            return
        elif not isinstance(in_edge, AtomicOperation):
            raise TypeError("Input edge must be an Atomic operation, not {}!".format(type(in_edge)))
        if isinstance(out_edge, EIWrapper):
            self.merge([out_edge])
            self.add_edges_from([(in_edge, entry) for entry in out_edge.entries])
            return
        elif not isinstance(out_edge, AtomicOperation):
            raise TypeError("Outgoing node must be an Atomic operation, not {}".format(type(out_edge)))
        self.validate_operation(in_edge)
        self.validate_operation(out_edge)

        self.__graph.add_edge(in_edge, out_edge)

    @staticmethod
    def validate_operation(operation: AtomicOperation):

        assert isinstance(operation, AtomicOperation)
        inputs = []
        for value in operation.input.values():
            if isinstance(value, list):
                inputs.extend(value)
            else:
                inputs.append(value)
        outputs = []
        for value in operation.output.values():
            if isinstance(value, list):
                outputs.extend(value)
            else:
                outputs.append(value)
        from collections import Counter
        count_inps = Counter(inputs)
        count_outs = Counter(outputs)
        if count_inps.most_common()[0][1] > 1:
            raise ValueError("Rule {operation.rulename} has repeated inputs: {count_inps}".format(**locals()))
        elif count_outs.most_common()[0][1] > 1:
            raise ValueError("Rule {operation.rulename} has repeated outputs: {count_outs}".format(**locals()))

        common = set.intersection(set(inputs), set(outputs))
        if common:
            raise ValueError(
                "Rule {operation.rulename} has inputs and outputs in common, \
                which causes a cyclic dependency: {common}.\n{operation.input}\n{operation.output}".format(
                    **locals()
                ))
        return True

    def add_edges_from(self, tasks):

        # First check that everything is actually an AtomicOperations
        for el in tasks:
            if len(el) != 2:
                raise TypeError("Invalid task edge: {}".format(el))
            self.add_edge(*el)
        return

    def __str__(self, source=None):

        self.check_graph()

        snake = ""

        for task in nx.algorithms.dfs_postorder_nodes(self.__graph, source=source):
            snake += str(task)

        return snake

    @property
    def resources(self):

        resources = []
        for rule in self.nodes:
            if rule.resources:
                resources.append(rule.resources)

        return "\n".join(resources)

    def __len__(self):
        return len(self.graph)

    def connect(self, inbound: AtomicOperation, outbound: AtomicOperation, name):

        self.add_edge(inbound, outbound)
        outbound.input[name] = inbound.output[name]

    @property
    def graph(self):
        return self.__graph.copy(as_view=True)

    def merge(self, others):

        self.__graph = nx.compose_all([self.__graph] + [other.graph for other in others])

    def __iter__(self):
        return iter(self.__graph.nodes)

    def remove_nodes(self, nodes):
        self.__graph.remove_nodes_from(nodes)

    def remove_node(self, node):
        self.__graph.remove_node(node)

    @property
    def adj(self):
        return self.__graph.adj

    @property
    def configuration(self):
        return self.__configuration

    @configuration.setter
    def configuration(self, conf):
        # TODO: this should be powered by JSON schemas
        if conf is None or isinstance(conf, (frozendict, dict)):
            self.__configuration = conf
        else:
            raise TypeError

    @property
    def nodes(self):
        return self.graph.nodes

    @property
    def edges(self):
        return self.graph.edges

    def __find_exits(self):
        """A workflow can have more than one final point. This method will identify and collect them all."""

        nodes = [_ for _ in self if not self.adj[_]]
        if len(nodes) == 0 and len(self) > 0:
            raise ValueError("Unexpected: no end nodes in a non-empty graph!")
        return nodes

    def add_final_flag(self, flag):
        """This method will find all the end nodes of the pipeline and connect them to the final flag."""
        finals = {"inputs": []}

        nodes = self.__find_exits()
        for node in nodes:
            for key, item in node.output.items():
                if isinstance(item, list):
                    finals["inputs"].extend(item)
                else:
                    finals["inputs"].append(item)

        # Remove redundancies
        finals["inputs"] = list(set(finals["inputs"]))
        flag = FinalFlag(finals, flag)
        self.add_edges_from([(node, flag) for node in nodes])


class FinalFlag(AtomicOperation):

    def __init__(self, inputs, flag):

        super().__init__()
        self.input = inputs
        self.output["flag"] = flag
        self.touch = True

    @property
    def rulename(self):
        return "all"

    @property
    def loader(self):
        return []


class EIWrapper(EIWorfkflow):

    @property
    def output(self):

        """EI wrappers must have *ONE* output"""

        return self.exit.output

    @property
    def exit(self):
        nodes = [_ for _ in self if not self.adj[_]]
        if len(nodes) == 0:
            raise ValueError("Something has gone wrong. Total nodes:\n{}".format(
                '\n'.join(self)
            ))

        if not len(nodes) == 1 and len(nodes[0].output) == 1:
            raise ValueError("EI wrappers must have at most one output! Output nodes:\n{}".format(
                "\n".join(_.rulename for _ in nodes)
            ))
        return nodes[0]

    @property
    def entries(self):
        """This property returns the roots of the workflow. Contrary to exits, there can be multiple
        entry points into a wrapper."""

        # Code snippet example:
        # g.add_edges_from([(0,1), (0,2), (2,3), (1,3)])
        # >>> list(zip(*g.edges))
        # [(0, 0, 1, 2), (1, 2, 3, 3)]
        # >>> set.difference(*[set(_) for _ in list(zip(*g.in_edges))])
        # {0}  # CVD

        in_edges, out_edges = zip(*self.graph.edges)
        return set.difference(set(in_edges), set(out_edges))

    def add_flag_to_inputs(self, flag, flag_name, key):
        assert isinstance(flag, (AtomicOperation, EIWrapper))

        self.add_edges_from([(flag, entry) for entry in self.entries])
        if isinstance(flag, EIWrapper):
            flag = flag.exit
        assert key in flag.output
        preds = nx.ancestors(self.graph, flag)
        for edge in self.edges:
            if edge[0] == flag and edge[1] not in preds:
                edge[1].input[flag_name] = flag.output[key]
