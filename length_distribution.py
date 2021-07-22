import os
import gzip
import argparse

import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
from Bio import SeqIO


def get_distribution(filename):
    try:
        if filename[-2:] == 'gz':
            f = gzip.open(filename, 'rt')
        else:
            f = open(filename)
        length_list = [len(s) for s in SeqIO.parse(f, 'fastq')]
        length_list = np.array(length_list)
        mean, std = length_list.mean(), length_list.std()
        min_len, max_len = length_list.min(), length_list.max()
    finally:
        f.close()

    fig, ax = plt.subplots()
    n, bins, patches = plt.hist(length_list, 50, rwidth=0.8, align='mid', density=True, alpha=0.75)

    x = np.linspace(min_len, max_len, 1000)
    plt.plot(x, stats.norm.pdf(x, mean, std))

    plt.title('Read length histogram')
    plt.grid(True)
    plt.show()
    plt.savefig('plots/length_hist.png')
    print('Length distribution info:')
    print('-' * 80)
    print(f'{mean = }\n{std = }\nmin = {min_len}\nmax = {max_len}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', type=str)
    args = parser.parse_args()
    filename = args.filename
    get_distribution(filename)

