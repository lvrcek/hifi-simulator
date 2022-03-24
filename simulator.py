import argparse
import os
import re
from random import random

from Bio import SeqIO
import numpy as np


def get_step(reference_length, length_mean, depth):
    num_reads = reference_length // length_mean * depth
    step_mean = (reference_length - length_mean) / (num_reads - 1)
    step_mean = round(step_mean / 100) * 100
    step_std = step_mean * 0.1
    return step_mean, step_std


def introduce_errors(read, subs, indel):
    i = 0
    while True:
        if i >= len(read):
            break

        if random() < subs:
            r = random()
            if r < 0.25:
                c = 'A'
            elif r < 0.5:
                c = 'C'
            elif r < 0.75:
                c = 'G'
            else:
                c = 'T'
            read = read[:i] + c + read[i+1:]
            i += 1
            continue

        if random() < indel:
            r = random()
            if r < indel / 2:
                read = read[:i] + read[i+1:]
                continue
            if r > 1 - indel / 2:
                read = read[:i] + 2 * read[i] + read[i+1:]  # Only homopolymer insertion, not a random one
                i += 1
                continue
        i += 1

    return read


def sample_strand(reference, reads_list, length_mean, length_std, step_mean, step_std, subs, indel, strand):
    idx = len(reads_list)
    position = 0
    stop = len(reference)

    while position < stop:
        length = int(np.random.normal(length_mean, length_std))
        if position + length < len(reference):
            read = reference[position:position+length]
        else:
            break
        read.id = str(idx)
        read.seq = introduce_errors(read.seq, subs, indel)
        if strand == '+':
            read.description = f'idx={idx}, strand=+, start={position}, end={position+length}'
        else:
            read.description = f'idx={idx}, strand=-, start={len(reference)-position}, end={len(reference)-position-length}'

        # read.letter_annotations = {'phred_quality': [50] * len(read)}
        reads_list.append(read)
        step = int(np.random.normal(step_mean, step_std))
        position += step
        idx += 1

    # return reads_list


def main(args):
    reference_path = os.path.abspath(args.ref)
    out_path = os.path.abspath(args.out)
    assert re.findall('.*\.([A-Za-z]*)', out_path)[-1] in ('fastq', 'fq', 'fasta', 'fa'), \
            "Output path should be in the FASTA/Q format"
    depth = args.depth  # That's actually depth per strand, it has to be divided by 2 (see below for step)
    length_mean = args.length_mean
    length_std = args.length_std if args.length_std is not None else length_mean * 0.075

    subs = args.subs
    indel = args.indel

    reference = next(SeqIO.parse(reference_path, 'fasta'))
    reference_rc = reference.reverse_complement()

    step_mean, step_std = get_step(len(reference), length_mean, depth // 2)
    reads_list = []

    # Sample positive and negative strand
    sample_strand(reference, reads_list, length_mean, length_std, step_mean, step_std, subs, indel, strand='+')
    sample_strand(reference_rc, reads_list, length_mean, length_std, step_mean, step_std, subs, indel, strand='-')
    SeqIO.write(reads_list, out_path, 'fasta')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Takes reference path and path where to store the generated reads.')
    parser.add_argument('--ref', type=str, default='debug/reference.fasta', help='reference path')
    parser.add_argument('--out', type=str, default='debug/reads.fasta', help='path where to store the reads')
    parser.add_argument('--length-mean', type=int, default=20000, help='mean length of the simulated reads')
    parser.add_argument('--length-std', type=int, help='standard deviation in length of the simulated reads')
    parser.add_argument('--depth', type=int, default=20, help='sequencing depth to be simulated')
    parser.add_argument('--subs', type=float, default=0.0, help='substitution percentage')
    parser.add_argument('--indel', type=float, default=0.0, help='indel percentage')
    args = parser.parse_args()
    main(args)

