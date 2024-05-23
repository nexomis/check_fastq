#!/usr/bin/env python

import gzip
import argparse
import hashlib
import time
import psutil
import random

def parse_fastq(handle, handle2, crypt_algo, num_bytes):
  """Yields sequences from a FASTQ file with optional cryptographic hashing."""
  while True:
    header = handle.readline().strip()
    if not header:
      break
    seq = handle.readline().strip()
    handle.readline()  # Skip line
    handle.readline()  # Skip line

    if handle2 is not None:
      handle2.readline()  # Skip line
      seq += "+" + handle2.readline().strip()
      handle2.readline()  # Skip line
      handle2.readline()  # Skip line

    seq = encode_seq(seq, crypt_algo, num_bytes)

    yield seq

def modify_base_randomly(sequence):
  """Modify a random base in the sequence to a random alternative base."""
  position = random.randint(0, len(sequence) - 1)
  current_base = sequence[position]
  possible_bases = {'A', 'C', 'G', 'T'}
  if current_base in possible_bases:
    possible_bases.remove(current_base)
  new_base = random.choice(list(possible_bases))
  return sequence[:position] + new_base + sequence[position+1:]

def count_records(file_name):
    """Count the number of records in a gzipped FASTQ file."""
    with gzip.open(file_name, "rt") as file:
      line_count = sum(1 for _ in file)
    return int(line_count / 4)

def encode_seq(seq, crypt_algo, num_bytes):
  seq = seq.encode('utf-8')
  if crypt_algo:
    hash_object = hashlib.new(crypt_algo)
    hash_object.update(seq)
    seq = hash_object.digest()
    if num_bytes != 0:
      seq = seq[:num_bytes]
  return(seq)

def parse_fastq_test(handle, handle2, crypt_algo, num_bytes, target_reads, n_mut_per_test):
  """Yields sequences from a FASTQ file with optional cryptographic hashing."""
  count = 0
  while True:
    n_diff = 0
    header = handle.readline().strip()
    if not header:
      break
    seq = handle.readline().strip()
    handle.readline()  # Skip line
    handle.readline()  # Skip line

    if handle2 is not None:
      handle2.readline()  # Skip line
      seq2 = handle2.readline().strip()
      handle2.readline()  # Skip line
      handle2.readline()  # Skip line
    else:
      seq2 = None

    if seq2:
      rseq = seq + "+" + seq2
    else:
      rseq = seq

    rseq = encode_seq(rseq, crypt_algo, num_bytes)
    
    if count in target_reads:
      for _ in range(n_mut_per_test):
        _seq = seq
        if seq2:
          _seq2 = seq2
          if random.sample([True, False]):
            _seq2 = modify_base_randomly(_seq2)
          else:
            _seq = modify_base_randomly(_seq)
        else:
          _seq = modify_base_randomly(_seq)
        _seq = encode_seq(_seq, crypt_algo, num_bytes)
        if _seq != rseq:
          n_diff += 1
    
    yield n_diff

def open_file(filename):
  """Opens a file, gzips it if necessary."""
  if filename.endswith(".gz"):
    return gzip.open(filename, "rt")
  else:
    return open(filename, "r")

def main(args):
  seq_count = {}
  same_seqs = True
  break_reason = "N/A"
  if args.test:
    total_test = args.n_test * args.mut_per_test
    if not ( (args.f1 == args.f2) and (args.r1 == args.r2) ):
      print("Test mode only works with identical files/pairs")

  start_time1 = time.time()  # Start timer for batch 1
  if not args.test:
    with open_file(args.f1) as handle1:
      handle1_r = open_file(args.r1) if args.r1 else None
      # Read first batch
      for seq in parse_fastq(handle1, handle1_r, args.hash, args.n_bytes):
        if seq in seq_count:
          seq_count[seq][0] += 1
        else:
          seq_count[seq] = [1, 0]
      if handle1_r:
        handle1_r.close()
  batch1_duration = time.time() - start_time1  # End timer for batch 1

  start_time2 = time.time()  # Start timer for batch 2

  if args.test:
    n_records = count_records(args.f1)
    target_reads = random.sample(range(n_records),args.n_test)

  with open_file(args.f2) as handle2:
    handle2_r = open_file(args.r2) if args.r2 else None
    # Read second batch
    if args.test:
      seq_diff = 0
      for n_diff in parse_fastq_test(handle2, handle2_r, args.hash, args.n_bytes, target_reads, args.mut_per_test):
        seq_diff += n_diff
    else:
      for seq in parse_fastq(handle2, handle2_r, args.hash, args.n_bytes):
        if seq in seq_count:
          seq_count[seq][1] += 1
        else:
          same_seqs = False
          print(seq_count)
          break_reason = "new key in second File/Pair"
          break
    if handle2_r:
      handle2_r.close()

  batch2_duration = time.time() - start_time2  # End timer for batch 2

  # Final check on counts
  if not args.test:
    if same_seqs:
      for counts in seq_count.values():
        if counts[0] != counts[1]:
          same_seqs = False
          break_reason = "mismatched occurences"
          break

  memory_usage = psutil.Process().memory_info().rss / (1024**2)  # Memory in MB
  if args.test:
    output_result(f'{seq_diff}/{total_test}', args.hash, args.n_bytes, break_reason, batch1_duration, batch2_duration, int(memory_usage))
  else:
    output_result(same_seqs, args.hash, args.n_bytes, break_reason, batch1_duration, batch2_duration, int(memory_usage))

def output_result(same_seqs, hash_algo, num_bytes, reason, time1, time2, memory):
  """Prints results in YAML format."""
  print(f"sequences_equal: {same_seqs}")
  print(f"hash_algorithm: {hash_algo}")
  print(f"max_bytes: {num_bytes}")
  print(f"reason: {reason}")
  print(f"batch1_duration_seconds: {time1:.2f}")
  print(f"batch2_duration_seconds: {time2:.2f}")
  print(f"resident_memory_MB: {memory}")

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description="Compare set of sequences based on their occurrences in two FASTQ files (or pair of files for paired-end).")
  parser.add_argument("-f1", required=True, help="File/Pair 1 forward read")
  parser.add_argument("-r1", help="File/Pair 1 reverse read")
  parser.add_argument("-f2", required=True, help="File/Pair 2 forward read")
  parser.add_argument("-r2", help="File/Pair 2 reverse read")
  parser.add_argument("--test", action="store_true", help="Test mode: only compare a random subset of reads")
  parser.add_argument("--n_test", type=int, default=10000, help="Number of reads to test")
  parser.add_argument("--mut_per_test", type=int, default=15, help="Number of mutations per test read")
  parser.add_argument("--n_bytes", type=int, default=0, help="Number of bytes of the hash to store (0 for full hash)")
  parser.add_argument("--hash", type=str, default="", help="Hashing algorithm to use, to reduce memory footprint (must be available in hashlib, e.g md5)")
  args = parser.parse_args()
  main(args)

