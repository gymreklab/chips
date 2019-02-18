#!/usr/bin/env python3
import pysam
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy.stats as spst
import scipy.special as sps
import math
import argparse
parser = argparse.ArgumentParser(
    description="Predict the fragment length distribution.")
parser.add_argument('bamfile')
parser.add_argument('bedfile')
parser.add_argument('threshold', type=float)
parser.add_argument('outdir', type=str)

args = parser.parse_args()
bamfile = args.bamfile
bedfile = args.bedfile
OUTDIR = args.outdir
threshold = args.threshold

peaks = []
frag_lens = []

with open(bedfile) as f:
    for line in f:
        arr = line.split("\t")
        if float(arr[6]) < threshold:
            break
        peaks.append([arr[0], int(arr[1]), int(arr[2])])

#bam file
with pysam.AlignmentFile(bamfile, "rb") as samfile:
    agg_starts = []
    agg_ends = []
    for (rn, start, stop) in peaks: #consider a peak
        starts = []
        ends = []
        for frag in samfile.fetch(rn, start, stop):
            if frag.is_paired and frag.is_proper_pair and\
               (not frag.is_unmapped) and (not frag.mate_is_unmapped)\
               and (not frag.is_unmapped):
                if (not frag.is_read1) and (not frag.is_read2):
                    print (frag)
                    exit(1)
                if frag.infer_query_length() != frag.reference_length:
                    continue
                #discard those fragments that touch the exterior of the region
                if frag.is_reverse:
                    frag_start = frag.next_reference_start
                    frag_end = frag.reference_start + frag.infer_query_length()
                else:
                    frag_start = frag.reference_start
                    frag_end = frag.next_reference_start + \
                               frag.infer_query_length()
                if frag_start < start or frag_end > stop:
                    continue
                frag_lens.append(abs(frag.template_length))
                if np.random.randint(2) == 0: #discard half at random
                    continue
                if frag.is_reverse:
                    ends.append(frag_end)
                else:
                    starts.append(frag_start)
        if len(starts) == 0 or len(ends) == 0:
            continue
        avg_starts = np.mean(starts)
        avg_ends = np.mean(ends)
        pos = (avg_ends + avg_starts) / 2
        agg_ends.extend([x - pos for x in ends])
        agg_starts.extend([x - pos for x in starts])

def find_mv_act(frag_lens):
    mu_act = np.mean(frag_lens)
    var_act = spst.moment(frag_lens, 2)
    return (mu_act, var_act)

def find_mv(lefts, rights):
    avg_starts = np.mean(lefts)
    avg_ends = np.mean(rights)
    lmin, lmax = math.floor(min(lefts)), math.ceil(max(lefts))
    rmin, rmax = math.floor(min(rights)), math.ceil(max(rights))
    r = range(800)
    mu = avg_ends - avg_starts #avg_len

    len_lefts = len(lefts)
#    lmax = int(lmax)
#    lmin = int(lmin)
    lefts_pdf = np.zeros(lmax-lmin+10)
    for x in lefts:
        lefts_pdf[math.ceil(x)] += 1
    lefts_cdf = np.zeros(lmax-lmin+10)
    for i in range(lmin-1,lmax):
        lefts_cdf[i+1] = lefts_cdf[i] + lefts_pdf[i+1]
    normed_lcdf = lefts_cdf / len_lefts

    ledf_r = np.zeros(lmax-lmin+10)
    for x in lefts:
        ledf_r[math.floor(x)] += x - math.floor(x)
    for i in range(lmax,lmin,-1):
        ledf_r[i-1] += ledf_r[i] + len_lefts - lefts_cdf[i]
    ledf_l = np.zeros(lmax-lmin+10)
    for x in lefts:
        ledf_l[math.ceil(x)] += math.ceil(x) - x
    for i in range(lmin-1,lmax):
        ledf_l[i+1] += ledf_l[i] + lefts_cdf[i]
    normed_ledf = (ledf_l + ledf_r) / len_lefts

    len_rights = len(rights)
    rights_pdf = np.zeros(rmax-rmin+10)
    for x in rights:
        rights_pdf[math.ceil(x)] += 1
    rights_cdf = np.zeros(rmax-rmin+10)
    for i in range(rmin-1,rmax):
        rights_cdf[i+1] = rights_cdf[i] + rights_pdf[i+1]
    normed_rcdf = rights_cdf / len_rights

    redf_r = np.zeros(rmax-rmin+10)
    for x in rights:
        redf_r[math.floor(x)] += x - math.floor(x)
    for i in range(rmax,rmin,-1):
        redf_r[i-1] += redf_r[i] + len_rights - rights_cdf[i]
    redf_l = np.zeros(rmax-rmin+10)
    for x in rights:
        redf_l[math.ceil(x)] += math.ceil(x) - x
    for i in range(rmin-1,rmax):
        redf_l[i+1] += redf_l[i] + rights_cdf[i]
    normed_redf = (redf_l + redf_r) / len_rights

    len_rights = len(rights)

    lefts_pdf /= len_lefts
    rights_pdf /= len_rights
    del lefts_cdf
    del rights_cdf

    def lefts_cdf(x):
        if x >= lmax:
            return 1
        elif x <= lmin:
            return 0
        elif math.floor(x) == x:
            return normed_lcdf[int(x)]
        else:
            return (math.ceil(x) - x) * normed_lcdf[math.floor(x)] \
                   + (x - math.floor(x)) * normed_lcdf[math.ceil(x)]

    def rights_cdf(x):
        if x >= rmax:
            return 1
        elif x <= rmin:
            return 0
        elif math.floor(x) == x:
            return normed_rcdf[int(x)]
        else:
            return (math.ceil(x) - x) * normed_rcdf[math.floor(x)] \
                   + (x - math.floor(x)) * normed_rcdf[math.ceil(x)]

    def lefts_edf(x):
        if x <= lmin:
            return -mu/2 - x
        elif x >= lmax:
            return x + mu/2
        else:
            xf = math.floor(x)
            return normed_ledf[xf] + (x - xf) * (2 * normed_lcdf[xf] - 1)

    def rights_edf(x):
        if x <= rmin:
            return mu/2 - x
        elif x >= rmax:
            return x - mu/2
        else:
            xf = math.floor(x)
            return normed_redf[xf] + (x - xf) * (2 * normed_rcdf[xf] - 1)

    gs1 = sum([lefts_pdf[i] * rights_edf((i+mu/2)/4+mu/2)
                   for i in range(lmin,lmax)])
    gs2 = sum([rights_pdf[i] * lefts_edf((i-mu/2)/4-mu/2)
                   for i in range(rmin,rmax)])
    gs = max(gs1, gs2)
    def score(v):
        shape = mu * mu / v
        scale = v / mu
        probs = spst.gamma.pdf(r, shape, loc=0, scale=scale)
        eexpl1 = sum([probs[i] * lefts_edf(-i/2) for i in r])
        eexpl2 = sum([probs[i] * rights_edf(i/2) for i in r])
        return min(eexpl1, eexpl2) - gs

    def search(low, high):
        res = score((low+high)/2)
        if res > 0:
            return (low, (low+high)/2)
        else:
            return ((low+high)/2, high)

    low = 200
    high = 8000
    for i in range(35):
        low, high = search(low, high)
    return (mu, (low+high)/2)

def plot_gamma(mu, v, mu_act, var_act):
    shape = mu * mu / v
    scale = v / mu
    shape_act = mu_act * mu_act / var_act
    scale_act = var_act / mu_act
    x, bins, z = plt.hist(frag_lens, 100, normed=True)
    plt.plot(bins, spst.gamma(shape, 0, scale).pdf(bins), color='r')
    plt.plot(bins, spst.gamma(shape_act, 0, scale_act).pdf(bins), color='g')
    plt.savefig(OUTDIR + "/" + bamfile.split("/")[-1] + "_" + str(int(threshold)) + "_" + str(np.random.randint(1000)) + ".png")

def print_frags():
    f = open(OUTDIR+ "/" + bamfile.split("/")[-1] + "_"+ str(int(threshold)) + "_" + str(np.random.randint(1000)) + ".frags.txt", "w")
    for fl in frag_lens: f.write("%s\n"%fl)
    f.close()

mu_act, var_act = find_mv_act(frag_lens)
mu, v = find_mv(agg_starts, agg_ends)
print(mu_act, var_act)
print(mu, v)
plot_gamma(mu, v, mu_act, var_act)
print_frags()
