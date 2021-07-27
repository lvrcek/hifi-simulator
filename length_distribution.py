import re
import os
import gzip
import argparse

import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
from Bio import SeqIO


def plot_lengths(lenghts):
    mean_len, std_len, median_len = np.mean(lengths), np.std(lengths), np.median(lenghts)
    min_len, max_len = lengths.min(), lengths.max()

    fig, ax = plt.subplots()
    n, bins, patches = plt.hist(lenghts, 50, rwidth=0.8, align='mid', density=True, alpha=0.75)

    x = np.linspace(min(bins), max(bins), 1000)
    lengths_log = np.log(lengths)
    mean_log, std_log = lengths_log.mean(), lengths_log.std()
    pdf = (np.exp(-(np.log(x) - mean_log)**2 / (2 * std_log**2)) / (x * std_log * np.sqrt(2 * np.pi)))
    plt.plot(x, pdf)

    plt.title('Read length histogram')
    plt.grid(True)
    pattern = r'.*/(.*).fast.*'
    filename = re.findall(pattern, filename)[0]
    plt.savefig(f'plots/{filename}_length.png')
    plt.show()
    with open(f'logs/{filename_lenghts.log}') as f:
        f.write(f'----- LENGTH INFO -----\n')
        f.write(f'Num seq\t\t=\t{len(lenghts)}\n')
        f.write(f'Median len\t=\t{median_len}\n')
        f.write(f'Mean len\t\t=\t{mean_len}\n')
        f.write(f'Std len\t\t=\t{std_len}\n')
        f.write(f'Max len\t\t=\t{max_len}\n')
        f.write(f'-------------------------')



def plot_qscores(qualities):
    mean_q, std_q, median_q = np.mean(qualities), np.std(qualities), np.median(qualities)
    min_q, max_q = qualities.min(), qualities.max()

    fig, ax = plt.subplots()
    n, bins, patches = plt.hist(qualities, 50, rwidth=0.8, align='mid', density=True, alpha=0.75)

    x = np.linspace(min(bins), max(bins), 1000)
    plt.plot(x, stats.norm.pdf(x, quality_mean, quality_std))

    plt.title('Q-score histogram')
    plt.grid(True)
    pattern = r'.*/(.*).fast.*'
    filename = re.findall(pattern, filename)[0]
    plt.savefig(f'plots/{filename}_qscore.png')
    plt.show()
    with open(f'logs/{filename_qscores.log}') as f:
        f.write(f'----- Q-SCORE INFO -----\n')
        f.write(f'Num seq\t\t=\t{len(qualities)}\n')
        f.write(f'Median Q\t=\t{median_q}\n')
        f.write(f'Mean Q\t\t=\t{mean_q}\n')
        f.write(f'Std Q\t\t=\t{std_q}\n')
        f.write(f'Max Q\t\t=\t{max_q}\n')
        f.write(f'--------------------------')



def main(args):
    filename = args.filename
    lengths, qualities = [], []
    try:
        if filename[-2:] == 'gz':
            f = gzip.open(filename, 'rt')
        else:
            f = open(filename)
        if find_l_dist:
            lengths = np.array([len(s) for s in SeqIO.parse(f, 'fastq')])
        if find_q_dist:
            qualities = np.array([np.array(r.letter_annotations['phred_quality']).mean() for r in SeqIO.parse(f, 'fastq')])
    finally:
        f.close()

    if find_l_dist:
        plot_lengths(lenghts)

    if find_q_dist:
        plot_qscores(qualities)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', type=str)
    args = parser.parse_args()
    main(args)

