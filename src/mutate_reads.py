import gzip
import random
import argparse
from Bio import SeqIO
from Bio.Seq import Seq

def count_records(file_name):
    """Count the number of records in a gzipped FASTQ file."""
    with gzip.open(file_name, "rt") as file:
        parser = SeqIO.parse(file, "fastq")
        count = sum(1 for _ in parser)
    return count

def modify_base_randomly(sequence):
    """Modify a random base in the sequence to a random alternative base."""
    position = random.randint(0, len(sequence) - 1)
    current_base = sequence[position]
    possible_bases = {'A', 'C', 'G', 'T'}
    if current_base in possible_bases:
        possible_bases.remove(current_base)
    new_base = random.choice(list(possible_bases))
    return sequence[:position] + new_base + sequence[position+1:]

def process_fastq(input_file, output_file, num_reads_to_modify, seed):
    num_records = count_records(input_file)
    random.seed(seed)
    indices_to_modify = set(random.sample(range(num_records), min(num_reads_to_modify, num_records)))

    with gzip.open(input_file, "rt") as handle_in, gzip.open(output_file, "wt") as handle_out:
        records = SeqIO.parse(handle_in, "fastq")
        for i, record in enumerate(records):
            if i in indices_to_modify:
                record.seq = Seq(modify_base_randomly(str(record.seq)))
            SeqIO.write([record], handle_out, "fastq")

def main():
    parser = argparse.ArgumentParser(description="Modify random bases in specified number of reads in a FASTQ file.")
    parser.add_argument("-f1", required=True, help="File/Pair 1 forward read")
    parser.add_argument("-r1", help="File/Pair 1 reverse read")
    parser.add_argument("-f2", required=True, help="File/Pair 2 forward read")
    parser.add_argument("-r2", help="File/Pair 2 reverse read")
    parser.add_argument("--num_reads_to_modify", type=int, default=1, help="Number of reads to randomly modify")
    parser.add_argument("--seed", type=int, default=None, help="Seed for random number generator for reproducibility")
    args = parser.parse_args()

    # Process each pair
    n_forward = random.randint(0, args.num_reads_to_modify)
    n_reverse = args.num_reads_to_modify - n_forward
    process_fastq(args.f1, args.f2, n_forward, args.seed)
    if args.r1 and args.r2:
        process_fastq(args.r1, args.r2, n_reverse, args.seed)

if __name__ == "__main__":
    main()
