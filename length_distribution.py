import re
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
        # length_list = [len(s) for s in SeqIO.parse(f, 'fastq')]
        # length_list = np.array(length_list)
        # mean, std = length_list.mean(), length_list.std()
        # min_len, max_len = length_list.min(), length_list.max()

        quality = np.array([np.array(r.letter_annotations['phred_quality']).mean() for r in SeqIO.parse(f, 'fastq')])
        quality_mean, quality_std = quality.mean(), quality.std()
    finally:
        f.close()

    # fig, ax = plt.subplots()
    # n, bins, patches = plt.hist(length_list, 50, rwidth=0.8, align='mid', density=True, alpha=0.75)

    # x = np.linspace(min_len, max_len, 1000)
    # lengths_log = np.log(length_list)
    # mean_log, std_log = lengths_log.mean(), lengths_log.std()
    # pdf = (np.exp(-(np.log(x) - mean_log)**2 / (2 * std_log**2)) / (x * std_log * np.sqrt(2 * np.pi)))
    # plt.plot(x, pdf)

    # plt.title('Read length histogram')
    # plt.grid(True)
    # pattern = r'.*/(.*).fast.*'
    # figname = re.findall(pattern, filename)[0]
    # plt.savefig(f'plots/{figname}.png')
    # plt.show()
    # print('Length distribution info:')
    # print('-' * 80)
    # print(f'{mean = }\n{std = }\nmin = {min_len}\nmax = {max_len}')

    fig, ax = plt.subplots()
    n, bins, patches = plt.hist(length_list, 50, rwidth=0.8, align='mid', density=True, alpha=0.75)
    x = np.linspace(min(bins), max(bins), 1000)
    plt.plot(x, stats.norm.pdf(x, quality_mean, quality_std))
    plt.title('Q-score histogram')
    plt.grid(True)
    pattern = r'.*/(.*).fast.*'
    figname = re.findall(pattern, filename)[0]
    plt.savefig(f'plots/{figname}_quality.png')
    plt.show()



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', type=str)
    args = parser.parse_args()
    filename = args.filename
    get_distribution(filename)

