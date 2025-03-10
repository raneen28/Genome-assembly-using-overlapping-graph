import itertools
import sys
import random
import argparse
import time

def read_fasta(filename):
    """
    Read the first sequence from a FASTA file
    Returns the sequence as a string (uppercase)
    """
    with open(filename, 'r') as f:
        lines = f.readlines()

    seq_lines =[line.strip() for line in lines if not line.startswith('>')]
    genome_seq = "".join(seq_lines).upper()
    return genome_seq

def Generate_error_free_reads(genome_seq, N, l):
    """
    in this function, given genome sequence we are generating N error-free reads of length l.
    the reads are randomly sampled across the genome to achieve uniform coverage (on average)
    :param genome_seq: string, the full PhiX genome
    :param N: int, the number of reads
    :param l: int, the length of each read
    :return: list of read strings
    """
    Genome_len = len(genome_seq)
    circular_seq = genome_seq + genome_seq
    reads =[]

    for i in range(N):
        start_position = random.randint(0, Genome_len - 1)
        curr_read = circular_seq[start_position : start_position + l]
        reads.append(curr_read)
    return reads

def mismatch_base(base):
    """
    given a base return a different base
    """
    bases =['A', 'C', 'G', 'T']
    bases.remove(base)
    return random.choice(bases)
def generate_error_prone_reads(genome_seq, N, l, p):
    """
    in this function, given genome sequence we are generating N error-prone reads of length l.
    where each base is mutated with probability p.
    :param genome_seq: string, the full PhiX genome
    :param N: int, the number of reads
    :param l: int, the length of each read
    :param p: float, the probability of mutating each base
    :return: list of read strings
    """
    Genome_len = len(genome_seq)
    circular_seq = genome_seq + genome_seq
    reads = []
    for i in range(N):
        start_position = random.randint(0, Genome_len - 1)
        curr_read_list = list(circular_seq[start_position: start_position + l])

        for j in range(l):
            if random.random() < p:
                curr_read_list[j] = mismatch_base(curr_read_list[j])
        reads.append("".join(curr_read_list))
    return reads
def generate_reads(genome_seq, N, l, p=0.0):
    """
    If p=0, return error-free reads.
    Otherwise, return error-prone reads.
    """
    if p == 0.0:
        return Generate_error_free_reads(genome_seq, N, l)
    else:
        return generate_error_prone_reads(genome_seq, N, l, p)


def overlap(read1, read2, min_overlap=1):
    """
    given two reads, Return the length of the maximum overlap between
    the suffix of read1 and the prefix of read2. If no overlap of at least min_overlap return 0
    """
    maximum_len = min(len(read1), len(read2))
    for length in range(maximum_len, min_overlap -1, -1):
        if read1.endswith(read2[:length]):
            return length
    return 0


def build_overlap_graph(reads, min_overlap=1):
    """
    Build the overlap graph from the list of reads.
    Returns a dictionary adjacency list: graph[i] = list of (j, overlap_length)
    where there's an edge from read i to read j if they overlap by >= min_overlap.
    """
    graph = {i: [] for i in range(len(reads))}
    for i in range(len(reads)):
        for j in range(len(reads)):
            if i != j:
                overlap_len = overlap(reads[i], reads[j], min_overlap=min_overlap)
                if overlap_len > 0:
                    graph[i].append((j, overlap_len))
    return graph


def greedy_assemble(reads, min_overlap=1):
    """
    in this function we build the overlap graph, collect all the edges (i,j,overlap_length)
    pick the pair with the largest overlap, Merge, then rebuild the read list and rebuild the graph
    repeat until no overlaps are bigger than min_overlap or only one read is left
    """
    reads = list(reads)

    while True:
        # Build overlap graph for the current reads
        graph = build_overlap_graph(reads, min_overlap)

        edges =[]
        for i,adj in graph.items():
            for (j, overlap_len) in adj:
                edges.append((i, j, overlap_len))

        if not edges:
            break

        max_overlap = 0
        best_edge =(None,None)
        for (i, j, overlap_len) in edges:
            if overlap_len > max_overlap:
                max_overlap = overlap_len
                best_edge = (i, j)

        if max_overlap < min_overlap:
            break

        i,j = best_edge
        merged_read = reads[i] + reads[j][max_overlap:]

        new_reads=[]
        for index in range(len(reads)):
            if index != i and index != j:
                new_reads.append(reads[index])
        new_reads.append(merged_read)

        reads = new_reads

    return reads


def compute_coverage(N, l, genome_len):
    """
    Compute average coverage: (N * l) / genome_length
    """
    return (N * l) / float(genome_len)


def basic_performance_metrics(contigs, reference):
    num_contigs = len(contigs)
    lengths = [len(c) for c in contigs]
    longest = max(lengths) if lengths else 0
    total_assembled_length = sum(lengths)
    return {
        "num_contigs": num_contigs,
        "longest_contig": longest,
        "total_assembled_length": total_assembled_length,
        "reference_length": len(reference),
    }



def main():
    fasta_file = sys.argv[1]
    reference_genome = read_fasta(fasta_file)
    G = len(reference_genome)
    print(f"Loaded reference genome of length {G} from {fasta_file}")

    # 1. Define parameter ranges for experimentation

    N_values = [250,500,1000]
    l_values =[50,100,150]
    p_values = [0.0, 0.01]
    min_overlap = 10

  
    print("\n RESULTS Naive Overlap Assembly")
    print("N\tl\tp\tCoverage\t#Contigs\tLongest\tTotalAsm\tTimeBuild(s)\tTimeAssemble(s)")

    
    for (N, l, p) in itertools.product(N_values, l_values, p_values):
        coverage_est = compute_coverage(N, l, G)

        # Generate reads
        start_time = time.time()
        reads = generate_reads(reference_genome, N, l, p)
        gen_time = time.time() - start_time

        # Build overlap graph
        start_time = time.time()
        graph = build_overlap_graph(reads, min_overlap=min_overlap)
        build_time = time.time() - start_time

        # Greedy Assembly
        start_time = time.time()
        contigs = greedy_assemble(reads, min_overlap=min_overlap)
        assemble_time = time.time() - start_time

        #Metrics
        metrics = basic_performance_metrics(contigs, reference_genome)
        num_contigs = metrics["num_contigs"]
        longest_contig = metrics["longest_contig"]
        total_asm_len = metrics["total_assembled_length"]

        # 4. Print or store results in a table-like format
        print(f"{N}\t{l}\t{p}\t{coverage_est:.2f}\t"
              f"{num_contigs}\t{longest_contig}\t{total_asm_len}\t"
              f"{build_time:.2f}\t{assemble_time:.2f}")



if __name__ == "__main__":
    main()






