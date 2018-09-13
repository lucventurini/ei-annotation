import magic
import sqlite3
from Mikado.exceptions import CorruptIndex
import collections
from Mikado.utilities.intervaltree import IntervalTree
from Mikado.loci import Gene
import json


def build_pos_index(index_name):

    positions = dict()

    wizard = magic.Magic(mime=True)

    if wizard.from_file("{0}".format(index_name)) == b"application/gzip":
        raise CorruptIndex("Invalid index file")
    try:
        conn = sqlite3.connect("{0}".format(index_name))
        cursor = conn.cursor()
        tables = cursor.execute("SELECT name FROM sqlite_master WHERE type='table';").fetchall()
        if sorted(tables) != sorted([("positions",), ("genes",)]):
            raise CorruptIndex("Invalid database file")

    except sqlite3.DatabaseError:
        raise CorruptIndex("Invalid database file")

    try:
        gene_positions = dict()
        for obj in cursor.execute("SELECT * from positions"):
            chrom, start, end, gid = obj
            if chrom not in positions:
                positions[chrom] = collections.defaultdict(list)
            positions[chrom][(start, end)].append(gid)
            gene_positions[gid] = (chrom, start, end)
    except sqlite3.DatabaseError:
        raise CorruptIndex("Invalid index file. Rebuilding.")

    indexer = collections.defaultdict(list).fromkeys(positions)
    for chrom in positions:
        indexer[chrom] = IntervalTree.from_tuples(
            positions[chrom].keys()
        )

    genes = dict()
    try:
        for obj in cursor.execute("SELECT * from genes"):
            gid, blob = obj
            genes[gid] = Gene(None)
            genes[gid].load_dict(json.loads(blob))
    except sqlite3.DatabaseError:
        raise CorruptIndex("Invalid index file. Rebuilding.")

    return indexer, positions, gene_positions, genes