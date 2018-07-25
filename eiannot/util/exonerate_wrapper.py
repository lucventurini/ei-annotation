#!/usr/bin/env python3

import argparse
import subprocess as sp
import multiprocessing as mp
import socket
import os
import time
import sys
from threading import Thread
from time import sleep
import tempfile


__doc__ = """This wrapper will execute Exonerate queries by first pre-loading the exonerate database into memory.
As this allows to use multi-threading AND a precomputed index, it should cut down computing times massively."""


class ExoStreamThread(Thread):

    def __init__(self, buffer, out_handle):
        Thread.__init__(self)
        self.buffer = buffer
        self.out = out_handle
        self.__listening = False
        self.__exit = False

    def run(self):

        while 1:
            line = self.buffer.readline().decode()
            print(line, end="", file=self.out)
            if "listening on port" in line:
                self.__listening = True
            if self.__exit is True:
                break
        # self.join()

    @property
    def listening(self):
        return self.__listening

    def send_exit(self):
        self.__exit = True


def get_open_port():
    # StackOverflow source: https://goo.gl/SrPDqh
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.bind(("", 0))
    s.listen(1)
    port = s.getsockname()[1]
    # s.close()
    return s, port


def percent(string: str):
    if not string.isdigit():
        raise ValueError("Invalid percent type")
    string = float(string)
    if not 0 < string <= 100:
        raise ValueError("Invalid percent, it should be between 0 and 100")
    return int(round(string, 0))


def pos(string: str):
    if not string.isdigit():
        raise ValueError
    string = int(string)
    if not string > 0:
        raise ValueError
    return string


def main():
    # Keep all data on local disks.
    # Apply the highest acceptable score thresholds using a combination of --score, --percent and --bestn.
    # Repeat mask and dust the genomic (target) sequence. (Softmask these sequences and use --softmasktarget).
    # Increase the --fsmmemory option to allow more query multiplexing.
    # Increase the value for --seedrepeat
    # When using an alignment model containing introns, set --geneseed as high as possible.
    # If you are compiling exonerate yourself, see the README file supplied with the source code for details of compile-time optimisations.

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("-t", "--threads", default=mp.cpu_count(), type=int)
    parser.add_argument("-M", "--memory", default=2000, type=float)
    parser.add_argument("-g", "--geneseed", default=250, type=pos)
    parser.add_argument("--bestn", default=10, type=pos)
    parser.add_argument("--hspfilter", default=100, type=pos)
    # parser.add_argument("--identity", required=True, type=percent)
    parser.add_argument("--serverlog", required=True)
    parser.add_argument("--log", required=True)
    parser.add_argument("-ir", "--intron-range", dest="intron_range", nargs=2, type=pos, required=True)
    parser.add_argument("genome")
    parser.add_argument("query")
    parser.add_argument("outfile")
    args = parser.parse_args()

    args.memory = abs(int(round(args.memory, 0)))

    if not args.genome.endswith(".esi"):
        raise ValueError("This wrapper requires a compiled genome!")

    args.intron_range = sorted(args.intron_range)
    assert args.intron_range[0] != args.intron_range[1]

    sock, port = get_open_port()
    log = args.serverlog
    server_cmd = "exonerate-server --preload --maxconnections {args.threads} "
    server_cmd += " --port {port} --input {args.genome}"
    server_cmd = server_cmd.format(**locals())
    sock.close()
    print("Starting the server on port {port}".format(**locals()), file=sys.stderr)
    server = sp.Popen(server_cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    log_handle = open(log, "wt")
    stderr = ExoStreamThread(server.stderr, log_handle)
    stdout = ExoStreamThread(server.stdout, log_handle)
    stderr.start(), stdout.start()

    while not any(_.listening is True for _ in [stderr, stdout]):
        sleep(1)

    # print(server_cmd)
    children = []
    # Use one LESS thread for aligning
    threads = max(args.threads - 1, 1)
    logs = []
    outputs = []
    for thread in range(1, threads + 1):

        output = tempfile.NamedTemporaryFile(mode="wt", suffix=".txt", delete=True)
        log = tempfile.NamedTemporaryFile(mode="wt", suffix=".log", delete=True)
        print(output.name, file=sys.stderr)
        logs.append(log)
        outputs.append(output)
        memory = args.memory / args.threads
        cmd = "exonerate --model protein2genome --showtargetgff yes --showvulgar yes "
        # Specify the chunk number
        cmd += " --querychunkid {thread} --querychunktotal {threads} "
        min_intron, max_intron = args.intron_range
        # INCREASE THE MEMORY!
        cmd += " -M {memory} -D {memory}"
        # Experimental
        cmd += " --hspfilter {args.hspfilter} "
        cmd += " --softmaskquery yes --softmasktarget yes --bestn {args.bestn} --minintron {min_intron} "
        cmd += " --maxintron {max_intron} --showalignment no "
        # These options have to be TAILORED, otherwise it will take forever!
        cmd += " --geneseed {args.geneseed} "
        cmd += " --query {args.query} --target localhost:{port} "
        cmd += " --ryo \">%qi\\tlength=%ql\\talnlen=%qal\\tscore=%s\\tpercentage=%pi\\nTarget>%ti\\tlength=%tl\\talnlen=%tal\\n\" "
        cmd = cmd.format(**locals())
        print(cmd, file=sys.stderr)
        child = sp.Popen(cmd, shell=True, stdout=output, stderr=log)
        children.append(child)

    statuses = [None for _ in range(threads)]
    while any(_ is None for _ in statuses):
        statuses = [_.poll() for _ in children]
        time.sleep(1)

    server.terminate()
    stderr.send_exit()
    stdout.send_exit()

    # Merge log files
    [log.flush() for log in logs]
    with open(args.log, "at") as final_log:
        sp.call(["cat"] + [log.name for log in logs], stdout=final_log)
        final_log.flush()
    [log.close() for log in logs]

    if any(_ != 0 for _ in statuses):
        print("Error in executing exonerate. Please check the logs.", file=sys.stderr)
        sys.exit(1)

    [output.flush() for output in outputs]
    with open(args.outfile, "wt") as out:
        sp.call(["cat"] + [output.name for output in outputs], stdout=out)
    [output.close() for output in outputs]

    # Merge output files

    sys.exit(0)


if __name__ == "__main__":
    main()
