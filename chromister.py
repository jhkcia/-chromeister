import numpy as np

class KmerHashFunction:

    def __init__(self, keyLength, z):
        self.keyLength = keyLength
        self.z = z
        self.CHAR_TO_POW4 = {
            'A': [0 for i in range(33)],
            'C': [1 * (4 ** i) for i in range(33)],
            'G': [2 * (4 ** i) for i in range(33)],
            'T': [3 * (4 ** i) for i in range(33)],
        }

    def hash(self, kmer):
        return self.hashKey(kmer), self.hashValue(kmer)

    def hashKey(self, kmer):
        return kmer[0: self.keyLength]

    def hashValue(self, kmer):
        value = 0
        kmer_len = len(kmer)

        for i in range(0, kmer_len, self.z):
            value += self.CHAR_TO_POW4[kmer[i]][kmer_len - (i + 1)]
        return value
# todo add chromister configuration file


class HitMatrixElement:
    def __init__(self, hashValue, dbIndex):
        self.hashValue = hashValue

        self.dbIndex = dbIndex
        self.hitCount = 0

        self.queryIndex = None
        self.repetition = False

    def __repr__(self):
        return "[Hash={} |DB: Index= {}, count= {}|Query: Index= {}, repitition={}]".format(self.hashValue, self.dbIndex, self.hitCount, self.queryIndex, self.repetition)


class Chromeister:
    def __init__(self, dbSequenceDir, querySequenceDir, kmerSize=32, hashKeyLength=12, Z=4):
        self.kmerSize = kmerSize
        self.hashFunction = KmerHashFunction(hashKeyLength, Z)
        self.hitMatrix = dict()
        self.Z=Z
        self.COMPLEMENT = {
            'A': 'T',
            'C': 'G',
            'G': 'C',
            'T': 'A',
        }
        self.dbSequenceDir = dbSequenceDir
        self.querySequenceDir = querySequenceDir

    def index(self):
        print("Loading database")
        self.__generateDictionary()
        print("Database loaded and of length "+str(self.dbSequenceSize))
        self.__matchInExactWords()
        print("Query length "+str(self.querySequenceSize))

    def __generateDictionary(self):
        file = open(self.dbSequenceDir)
        chars = ['' for _ in range(self.kmerSize)]
        kmers = dict()
        current_kmer_size = 0
        index = 0
        for line in file.readlines():
            if not line.startswith('>'):
                for c in line:
                    C = c.upper()
                    if C == "A" or C == "T" or C == "C" or C == "G":
                        index += 1
                        chars[current_kmer_size] = C
                        current_kmer_size += 1

                        if current_kmer_size == self.kmerSize:
                            kmer = ''.join(chars)
                            current_kmer_size = 0
                            key = self.hashFunction.hashKey(kmer)
                            if(key in self.hitMatrix):
                                self.hitMatrix[key].repetition = True
                            else:
                                value = self.hashFunction.hashValue(kmer)
                                self.hitMatrix[key] = HitMatrixElement(
                                    value, index)
                    else:
                        if C != '\n' and C != '>':
                            index += 1
                            current_kmer_size = 0
            else:
                current_kmer_size = 0
        self.dbSequenceSize = index

    def __addKmerToHitMatrix(self, kmer, index):
        key = self.hashFunction.hashKey(kmer)
        if key in self.hitMatrix and self.hitMatrix[key].repetition == False:
            value = self.hashFunction.hashValue(kmer)
            if self.hitMatrix[key].hashValue == value:
                self.hitMatrix[key].hitCount += 1
                # which index? first or last?
                self.hitMatrix[key].queryIndex = index

    def __matchInExactWords(self):
        file = open(self.querySequenceDir)
        chars = ['' for _ in range(self.kmerSize)]
        complement_chars = ['' for _ in range(self.kmerSize)]
        kmers = dict()
        current_kmer_size = 0
        index = 0
        kmer = ""
        kmer_complement = ""
        for line in file.readlines():
            if not line.startswith('>'):
                for c in line:
                    C = c.upper()
                    if C == "A" or C == "T" or C == "C" or C == "G":
                        index += 1
                        kmer += C
                        kmer_complement = self.COMPLEMENT[C]+kmer_complement
                        current_kmer_size += 1
                        if current_kmer_size == self.kmerSize:
                            self.__addKmerToHitMatrix(kmer, index)
                            self.__addKmerToHitMatrix(kmer_complement, index)
                            current_kmer_size -= 1
                            kmer = kmer[1:]
                            kmer_complement = kmer_complement[:-1]
                    else:
                        if C != '\n' and C != '>':
                            index += 1
                            current_kmer_size = 0
                            kmer = ""
                            kmer_complement = ""
            else:
                current_kmer_size = 0
                kmer = ""
                kmer_complement = ""
        self.querySequenceSize = index

    def downSample(self, dimension):
        self.downSampleDimension=dimension
        print("Generating a {}x{} Matrix.".format(dimension, dimension))
        representation = [[0 for _ in range(dimension)]
                          for _ in range(dimension)]
        queryPixelSize = dimension / self.querySequenceSize
        queryRatio = self.querySequenceSize / dimension

        dbPixelSize = dimension / self.dbSequenceSize
        dbRatio = self.dbSequenceSize / dimension

        i_r_fix = max(1.0, self.kmerSize * queryPixelSize)
        j_r_fix = max(1.0, self.kmerSize * dbPixelSize)

        print("Ratios: Q [{}] D [{}]. Lenghts: Q [{}] D [{}]".format(
            queryRatio, dbRatio, self.querySequenceSize, self.dbSequenceSize))
        print("Pixel size: Q [{}] D [{}].".format(queryPixelSize, dbPixelSize))
        print("Computing absolute hit numbers")
        print("Scanning hits table.")
        ## or iterate over kmers
        for key in self.hitMatrix.keys():
            if self.hitMatrix[key].hitCount == 1:
                dbRedir = int(self.hitMatrix[key].dbIndex / dbRatio)
                queryRedir = int(self.hitMatrix[key].queryIndex / queryRatio)
                i_r, j_r = i_r_fix, j_r_fix

                while (int(i_r) >= 1) and (int(j_r) >= 1):
                    if queryRedir > int(i_r) and dbRedir > int(j_r):
                        representation[queryRedir -
                                       int(i_r)][dbRedir - int(j_r)] += 1
                    else:
                        representation[queryRedir][dbRedir] += 1
                        break

                    i_r -= min(1.0, queryPixelSize)
                    j_r -= min(1.0, dbPixelSize)
        return representation


    