import re
import os
import gzip
import argparse

import numpy as np
import scipy.stats as stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from Bio import SeqIO


def plot_lengths(lengths, name):
    mean_len, std_len, median_len = np.mean(lengths), np.std(lengths), np.median(lengths)
    min_len, max_len = lengths.min(), lengths.max()

    fig, ax = plt.subplots()
    n, bins, patches = plt.hist(lengths, 50, rwidth=0.8, align='mid', density=True, alpha=0.75)

    x = np.linspace(min(bins), max(bins), 1000)
    lengths_log = np.log(lengths)
    mean_log, std_log = lengths_log.mean(), lengths_log.std()
    pdf = (np.exp(-(np.log(x) - mean_log)**2 / (2 * std_log**2)) / (x * std_log * np.sqrt(2 * np.pi)))
    plt.plot(x, pdf)

    plt.title('Read length histogram')
    plt.grid(True)
    plt.savefig(f'plots/{name}_length.png')
    with open(f'logs/{name}_length.log', 'w') as f:
        f.write(f'----- LENGTH INFO -----\n')
        f.write(f'Num seq\t\t=\t{len(lengths)}\n')
        f.write(f'Median len\t=\t{median_len}\n')
        f.write(f'Mean len\t\t=\t{mean_len}\n')
        f.write(f'Std len\t\t=\t{std_len}\n')
        f.write(f'Max len\t\t=\t{max_len}\n')
        f.write(f'-------------------------')



def plot_qscores(qualities, name):
    mean_q, std_q, median_q = np.mean(qualities), np.std(qualities), np.median(qualities)
    min_q, max_q = qualities.min(), qualities.max()

    fig, ax = plt.subplots()
    n, bins, patches = plt.hist(qualities, 50, rwidth=0.8, align='mid', density=True, alpha=0.75)

    x = np.linspace(min(bins), max(bins), 1000)
    plt.plot(x, stats.norm.pdf(x, mean_q, std_q))

    plt.title('Q-score histogram')
    plt.grid(True)
    plt.savefig(f'plots/{name}_qscore.png')
    with open(f'logs/{name}_qscores.log', 'w') as f:
        f.write(f'----- Q-SCORE INFO -----\n')
        f.write(f'Num seq\t\t=\t{len(qualities)}\n')
        f.write(f'Median Q\t=\t{median_q}\n')
        f.write(f'Mean Q\t\t=\t{mean_q}\n')
        f.write(f'Std Q\t\t=\t{std_q}\n')
        f.write(f'Max Q\t\t=\t{max_q}\n')
        f.write(f'--------------------------')


def plot_accuracies(accuracies, name):
    mean_acc, std_acc, median_acc = np.mean(accuracies), np.std(accuracies), np.median(accuracies)
    min_acc, max_acc = accuracies.min(), accuracies.max()

    fig, ax = plt.subplots()
    n, bins, patches = plt.hist(accuracies, 50, rwidth=0.8, align='mid', density=True, alpha=0.75)

    x = np.linspace(min(bins), max(bins), 1000)
    plt.plot(x, stats.norm.pdf(x, mean_acc, std_acc))

    plt.title('Accuracy histogram')
    plt.grid(True)
    plt.savefig(f'plots/{name}_acc.png')
    with open(f'logs/{name}_acc.log', 'w') as f:
        f.write(f'----- Q-SCORE INFO -----\n')
        f.write(f'Num seq\t\t=\t{len(accuracies)}\n')
        f.write(f'Median acc\t=\t{median_acc}\n')
        f.write(f'Mean acc\t\t=\t{mean_acc}\n')
        f.write(f'Std acc\t\t=\t{std_acc}\n')
        f.write(f'Max acc\t\t=\t{max_acc}\n')
        f.write(f'--------------------------')


def get_acc(q):
    return 1.0 - 10 ** -(q / 10)


def main(args):
    filename = args.filename
    find_l_dist = args.length
    find_q_dist = args.qscore
    find_acc_dist = args.accuracy

    try:
        if filename[-2:] == 'gz':
            f = gzip.open(filename, 'rt')
        else:
            f = open(filename)
        if find_l_dist:
            lengths = np.array([len(s) for s in SeqIO.parse(f, 'fastq')])
            f.seek(0)
        if find_q_dist:
            qualities = np.array([np.array(r.letter_annotations['phred_quality']).mean() \
                                  for r in SeqIO.parse(f, 'fastq')])
            f.seek(0)
        if find_acc_dist:
            accuracies = np.array([np.array(list(map(get_acc, r.letter_annotations['phred_quality']))).mean() \
                                   for r in SeqIO.parse(f, 'fastq')])
            f.seek(0)
    finally:
        f.close()

    pattern = r'.*/(.*).fast.*'
    name = re.findall(pattern, filename)[0]

    if find_l_dist:
        print('Processing lengths...')
        plot_lengths(lengths, name)

    if find_q_dist:
        print('Processing q-scores...')
        plot_qscores(qualities, name)

    if find_acc_dist:
        print('Processing accuracies...')
        plot_accuracies(accuracies, name)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot distributions')
    parser.add_argument('filename', type=str)
    parser.add_argument('--length', action='store_true', default=False)
    parser.add_argument('--qscore', action='store_true', default=False)
    parser.add_argument('--accuracy', action='store_true', default=False)
    args = parser.parse_args()
    main(args)

