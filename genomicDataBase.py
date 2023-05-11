from typing import Dict, Set
import copy
import random

# Function that creates a sequence and returns the sequence itself
Sequence = str
Kmer = str
Label = str

def generateSequence(length) -> Sequence:
    return ''.join(random.choices(['A', 'C', 'G', 'T'], k=length))

def getKmersFromSequence(sequence, kmer_size):
    if len(sequence) < kmer_size:
        return []
    return [sequence[i: i+kmer_size] for i in range(len(sequence) - kmer_size + 1)]

class GenomicDB(object):
    def __init__(
        self,
        kmer_size
    ) -> None:
        self.kmer_size = kmer_size
        self.orig_db: Dict[Kmer, Set[Label]] = {} # Set for efficient search
        self.db = copy.deepcopy(self.orig_db)
        self.threshold = -1

    def setThreshold(
        self,
        thresh
    ) -> None:
        # Sanity check needed
        self.threshold = thresh
        self.db = copy.deepcopy(self.orig_db)
        if self.threshold == -1:
            return
        kmers_to_trim = []
        for kmer in self.db:
            if len(self.db[kmer]) > self.threshold:
                kmers_to_trim.append(kmer)
        for kmer in kmers_to_trim:
            self.db.pop(kmer)

    def addToDB(
        self,
        sequence: Sequence,
        name: str
    ) -> None:
        kmers = getKmersFromSequence(sequence, self.kmer_size)
        for kmer in kmers:
            if kmer in self.orig_db:
                self.orig_db[kmer].add(name)
            else:
                self.orig_db[kmer] = set([name])
        self.setThreshold(self.threshold)


    def querySequence(
        self,
        sequence) -> None:
        results = {}
        true_label = None
        true_label_count = 0
        for kmer in getKmersFromSequence(sequence, self.kmer_size):
            if kmer in self.db:
                for label in self.db[kmer]:
                    if not label in results:
                        results[label] = 0
                    results[label] += 1
                    if results[label] > true_label_count:
                        true_label = label
                        true_label_count = results[label]
        return true_label
