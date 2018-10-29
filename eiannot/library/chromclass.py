from sqlalchemy import Column
from sqlalchemy import String
from sqlalchemy import Integer
from sqlalchemy import ForeignKey
from sqlalchemy import CHAR
from sqlalchemy import Index
from sqlalchemy import Float
from sqlalchemy import select
from sqlalchemy.orm import relationship, column_property
from sqlalchemy.engine.reflection import Inspector
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.engine import create_engine, Engine
from sqlalchemy import event
from sqlalchemy_utils import database_exists, create_database
import numbers
import sqlite3

Inspector = Inspector
AugBase = declarative_base()


class Chrom(AugBase):
    """
    Simple serialization class for chromosomes.
    """

    __tablename__ = "chrom"
    __table_args__ = {"extend_existing": True}

    chrom_id = Column(Integer, primary_key=True)
    name = Column(String(200), unique=True)
    length = Column(Integer, nullable=True)

    def __init__(self, name, length=None):
        self.name = name
        if length is not None:
            assert isinstance(length, int)
        self.length = length


class Chunk(AugBase):
    """
    Class to specify the span of each chunk.
    """

    __tablename__ = "chunk"

    prog_id = Column(Integer, primary_key=True)
    chunk_id = Column(Integer, primary_key=False, nullable=False, unique=False)
    chrom_id = Column(Integer, ForeignKey(Chrom.chrom_id), primary_key=False, nullable=False, unique=False)
    start = Column(Integer, primary_key=False, nullable=False, unique=False)
    end = Column(Integer, primary_key=False, nullable=False, unique=False)

    chrom_object = relationship(Chrom, uselist=False)
    chrom = column_property(select([Chrom.name]).where(chrom_id == Chrom.chrom_id))
    __table_args__ = (Index("chunk_id", "chrom_id", "start", "end"), {"extend_existing": True})

    def __init__(self, chunk_id, chrom_id, start, end):

        for key, obj in (("chunk_id", chunk_id),
                         ("chrom_id", chrom_id),
                         ("start", start),
                         ("end", end)):
            if not isinstance(obj, numbers.Number):
                raise ValueError("Invalid value for {}: {}".format(key, obj))
            if not isinstance(obj, numbers.Integral):
                raise ValueError("Invalid non-integer number for {}: {}".format(key, obj))
            obj = int(obj)
            if obj < 0:
                raise ValueError("Invalid number")
            self.__dict__[key] = obj
