__author__ = 'Joseph Korpela'

import time
import os
import cPickle as pickle

class CNV:
    """
    Locates Copy Number Variations in a donor genome using a reference genome and short reads from a
    high throughput sequencer as input.
    """
    def __init__(self, allowable_error, key_size, cnv_key_size, ref_file, read_file, ans_file):
        self._allowable_error = allowable_error #max number of alleles that can differ from reference
        self._genome_id = '' #used when  matching answers to answer key during grading
        self._cnv_key_size = cnv_key_size #key size to use in hash maps
        self._read_file = read_file #file object with the reference sequence
        self._ans_file = ans_file #file object which will hold the answers
        self._ref_genome = self.read_ref_file(ref_file) #reference genome used when mapping short reads
        self._count_genome = [0] * len(self._ref_genome) #used to count depth of coverage
        self._count_reads = [0] * len(self._ref_genome)

        self._cnv_index = {}
        self._cnv_list = []
        self._cnv_dict = {}
        self.build_cnv_index(store_index=False)
        self.find_cnv()

        self._index = {} #maps k-mers to positions in the reference genome
        self._key_size = key_size #k-mer size to use for mapping sections of dna to the index
        self.build_index(store_index=False)

    def build_index(self, store_index):
        """
        Builds the index used when mapping reads to the reference genome
        :param store_index: Flags whether the index should be stored to file for later reuse.
        """
        file_name = self._genome_id[1:].rstrip() + '.p' #file name to use when storing the index
        if store_index and os.path.isfile(file_name):
            self._index = pickle.load(open(file_name, 'rb'))
        else:
            for i in range(len(self._ref_genome)):
                if i + self._key_size <=  len(self._ref_genome):
                    next_key = self._ref_genome[i:i+self._key_size]
                    self.add_posn(next_key, i)
            if store_index:
                pickle.dump(self._index, open(file_name, 'wb'))

    def add_posn(self, key, posn):
        """
        Adds k-mers and their locations in the reference genome to the index
        :param key: k-mer that is being mapped
        :param posn: position of that k-mer in the reference genome
        """
        if self._index.has_key(key):
            self._index[key].add(posn)
        else:
            self._index[key] = set()
            self._index[key].add(posn)

    def process_reads(self):
        """
        Processes all high throughput sequencer short reads
        """
        self._read_file.seek(0)
        for line in self._read_file:
            line = line.rstrip()
            for one_read in line.split(','): #short reads are given as paired end reads
                self.process_one_read(one_read)

    def process_one_read(self, one_read):
        """
        Processes a single high throughput sequencer short read
        :param one_read:
        :return:
        """
        match_sets = []
        max_count = 0
        first_start = -1
        number_splits = len(one_read) / self._key_size
        cnv_posn_set = set() #holds all CNV positions that were matched in this read
        cnv_list = [] #list of all matched CNVs for this read
        for cnv, posn_list in self._cnv_list: #find any portion of this read that matches a CNV
            for posn in posn_list:
                for i in range(len(cnv)-self._key_size):
                    cnv_posn_set.add(posn+i)
                cnv_list.append((cnv, posn, posn+i)) #cnv, cnv stt posn in ref genome, last posn in ref genome that matches
        for i in range(number_splits): #match all portions of this read to possible locations in ref genome
            stt = i * self._key_size
            stop = stt + self._key_size
            this_set = self._index.get(one_read[stt:stop])
            if this_set:
                match_sets.append(this_set)
            else:
                match_sets.append(set())
        #put all the returned positions into one set and then check if any of the positions matched with the CNV
        combined_set = set()
        for next_set in match_sets:
            combined_set = combined_set.union(next_set)
        #if the intersection is > 1, then consider this as containing part of the CNV, so use this read to map a
        #position for the given CNV
        this_cnv_set = combined_set.intersection(cnv_posn_set)
        combined_set = combined_set.difference(cnv_posn_set)
        if len(this_cnv_set) > 1:
            cnv_stt = 0
            cnv_end = 0
            cnv_seq = ''
            for seq, stt, end in cnv_list:
                if min(this_cnv_set) < end and max(this_cnv_set) > stt:
                    cnv_seq = seq
                    cnv_stt = stt
                    cnv_end = end
            #use the remaining matches from this read to map out possible boundary positions for the CNV
            last_i = -1
            for i in range(number_splits): #using each split of the read as the starting posn
                this_list = []
                if len(match_sets[i]) > 0:
                    for j in match_sets[i]:
                        if j not in this_cnv_set:
                            this_list.append(j)
                if len(this_list) > 0:
                    for j in this_list:
                        if cnv_seq in self._cnv_dict:
                            self._cnv_dict[cnv_seq].append(j)
                        else:
                            self._cnv_dict[cnv_seq] = []
                            self._cnv_dict[cnv_seq].append(j)

    def clean_cnv(self):
        """
        Takes the raw list of CNVs and possible boundary positions and removes probable false positives to
        leave the likely list of CNVs along with their boundary positions
        :return:
        """
        for key, value in list(self._cnv_dict.items()):
            raw_list = sorted(value)

            grouped_list = []
            grouped_list.append([])
            prev_posn = raw_list[0]
            sublist = 0
            grouped_list[sublist].append(prev_posn)

            for i in range(1, len(raw_list)):
                #all CNVs in this dataset are smaller than 50, so only boundary positions within a range
                #of 100 of each other are considered as being in a single set
                if raw_list[i] - prev_posn < 100:
                    grouped_list[sublist].append(raw_list[i]) #each sublist represents positions around 1 CNV loc
                else:
                    sublist += 1
                    grouped_list.append([])
                    grouped_list[sublist].append(raw_list[i])
                prev_posn = raw_list[i]

            #remove any list that is too small
            i = 0
            while i < len(grouped_list):
                #too few items to determine cnv
                if len(grouped_list[i]) < 3:
                    grouped_list.pop(i)
                #too small a range to determine cnv
                elif grouped_list[i][len(grouped_list[i])-1] - grouped_list[i][0] < 30:
                    grouped_list.pop(i)
                else:
                    i += 1

            cnv_stt_posns = []
            cnv_len = len(key)

            #for each possible CNV location, check for most likely gap position in the list to use
            #as the CNV boundary positions
            for j in range(len(grouped_list)):
                prev_posn = grouped_list[j][0]
                max_gap = 0
                gap_stt = 0
                this_posn = -1
                for i in range(1, len(grouped_list[j])):
                    posn_diff = grouped_list[j][i] - prev_posn
                    if posn_diff > cnv_len + 10:
                        this_posn = prev_posn + 10
                    elif posn_diff > cnv_len:
                        this_posn = grouped_list[j][i] - cnv_len
                    else:
                        if posn_diff > max_gap:
                            max_gap = posn_diff
                            gap_stt = prev_posn
                    #if no large gap was found, then use the largest seen as the cnv start posn
                    if i == len(grouped_list[j]) - 1 and this_posn == -1:
                        this_posn = gap_stt
                    prev_posn = grouped_list[j][i]
                if this_posn != 0:
                    cnv_stt_posns.append(this_posn)
            #if enough valid copies were found, then treat as a cnv, else delete this entry
            if len(cnv_stt_posns) > 4:
                self._cnv_dict[key] = cnv_stt_posns
            else:
                del self._cnv_dict[key]

    def read_ref_file(self, ref_file):
        """
        Reads the reference genome in from file and saves it as a string
        :param ref_file:
        :return:
        """
        genome = ''
        for line in ref_file:
            if line.startswith('>'):
                if self._genome_id == '':
                    self._genome_id = line
            else:
                line = line.rstrip()
                genome += line
        return genome

    def build_cnv_index(self, store_index):
        """
        Creates an index mapping k-mers to positions in the reference genome, similar to the normal index
        but allows a different sized key to be used when determining the depth of coverage over alleles
        in the reference genome
        :param store_index: whether to store a copy of the index for later reuse
        """
        file_name = self._genome_id[1:].rstrip() + '.p'
        if store_index and os.path.isfile(file_name):
            self._cnv_index = pickle.load(open(file_name, 'rb'))
        else:
            for i in range(len(self._ref_genome)):
                if i + self._cnv_key_size <=  len(self._ref_genome):
                    next_key = self._ref_genome[i:i+self._cnv_key_size]
                    self.add_cnv_idx_posn(next_key, i)
            if store_index:
                pickle.dump(self._cnv_index, open(file_name, 'wb'))

    def add_cnv_idx_posn(self, key, posn):
        """
        Adds a single k-mer position to the CNV index
        :param key: k-mer from the reference genome
        :param posn: corresponding position of that k-mer
        """
        if self._cnv_index.has_key(key):
            self._cnv_index[key].add(posn)
        else:
            posn_set = set()
            posn_set.add(posn)
            self._cnv_index[key] = posn_set

    def count_reads(self):
        """
        Processes reads to count their depth of coverage over alleles in the reference genome
        """
        for line in self._read_file:
            if line.startswith('>'):
                if self._genome_id == '':
                    self._genome_id = line
            else:
                line = line.rstrip()
                for one_read in line.split(','):
                    self.count_one_read(one_read)

    def count_one_read(self, one_read):
        """
        Processes a single read for counting depth of coverage
        :param one_read:
        """
        posn_set = set()
        #count = 0
        for i in range(len(one_read) - self._cnv_key_size):
            key = one_read[i:i+self._cnv_key_size]
            if key in self._cnv_index:
                posn_set = posn_set.union(self._cnv_index[key])
        for posn in posn_set:
            self._count_reads[posn] += 1

    def find_cnv(self):
        """
        Determines possible CNV locations based on depth of coverage of reads over the reference genome
        :return:
        """
        self.count_reads()
        this_cnv = ''
        this_stt = 0
        for i in range(len(self._count_reads)):
            if self._count_reads[i] > 19:
                if this_cnv == '':
                    this_cnv += self._ref_genome[i:i+self._cnv_key_size]
                    this_stt = i
                else:
                    this_cnv += self._ref_genome[i+self._cnv_key_size]
            else:
                if len(this_cnv) > 19:
                    if not self.str_check(this_cnv,3,5):
                        self._cnv_list.append((this_cnv, [this_stt]))
                    this_cnv = ''
                    this_stt = 0
                else:
                    this_cnv = ''
                    this_stt = 0

    def str_check(self, sequence, min_str, max_str):
        """
        The repetitive nature of STRs results in many false positives, so these are being checked for and removed
        :param sequence: dna sequence of a possible CNV
        :param min_str: min length allowed for STR
        :param max_str: max length allowed for STR
        :return:
        """
        while len(sequence) / max_str < 4:
            if max_str > min_str:
                max_str -= 1
            else:
                return None #since this sequence is too short to id an str
        #check for repeating patterns of several alleles (STRs) starting from all points in this sequence
        for i in range(len(sequence)):
            for str_size in reversed(range(min_str, max_str+1)):
                if i + (4 * str_size) < len(sequence):
                    if sequence[i:i+str_size] == sequence[i+str_size:i+2*str_size] and \
                        sequence[i+str_size:i+2*str_size] == sequence[i+2*str_size:i+3*str_size] and \
                        sequence[i+2*str_size:i+3*str_size] == sequence[i+3*str_size:i+4*str_size]:
                        return (sequence[i:i+str_size], i)
        return None

    def create_answer_file(self):
        """
        Used to submit possible CNVs and their boundary positions for evaluation
        :return:
        """
        self._ans_file.write(self._genome_id)
        self._ans_file.write('>COPY')
        for key in self._cnv_dict:
            if len(self._cnv_dict[key]) > 0:
                self._ans_file.write('\n1,' + key)
                for posn in self._cnv_dict[key]:
                    self._ans_file.write(',' + str(posn))

def main():
    ref_file_path = 'data/ref_genomeE1.txt'
    reads_file_path = 'data/reads_genomeE1.txt'
    ans_file_path = 'output/ans_genomeE1.txt'
    with open(ref_file_path, 'r') as ref_file:
        with open(reads_file_path, 'r') as read_file:
            with open(ans_file_path, 'w') as answer_file:
                start = time.time()
                cnv = CNV(2, 10, 25, ref_file, read_file, answer_file)
                cnv.process_reads()
                cnv.clean_cnv()
                cnv.create_answer_file()
                print 'Sequencing time: ' + str(time.time() - start)

if __name__ == '__main__':
    main()
