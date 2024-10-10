import shutil
from os import mkdir
from os.path import dirname, abspath, join, exists
from unittest import TestCase

from seqteleporter.config import CODON_TABLE
from seqteleporter.partitioner.compute_best_partitions import compute_best_partitions



