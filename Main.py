import sys
import dotPlotDrawer
from chromister import Chromeister
import time
import IOUtils
import synteniesIdentification
import numpy as np
import argparse
from os import path

def main(args):
    chromister = Chromeister(args.db, args.query,
                             args.k, args.hash_key_len, args.z)
    chromister.index()
    representationMatrix = chromister.downSample(args.dim)

    IOUtils.storeRepresentationMatrix(
        representationMatrix, args.o+"/hit.mat")
    IOUtils.storeXYMatrix(representationMatrix, args.o+"/hit.xy.mat")

    dotPlotDrawer.drawDotPlot(representationMatrix,
                              args.o+"/hit_raw.png", "HIT_MATRIX_RAW")
    cleanData = synteniesIdentification.cleanData(
        np.array(representationMatrix), args.diag_len)
    score = synteniesIdentification.calculateScore(cleanData)
    dotPlotDrawer.drawDotPlot(
        cleanData, args.o+"/hit_clean.png", "HIT_MATRIX [SCORE= "+str(score)+" ]")
    dsData = synteniesIdentification.downsample(cleanData, args.diag_len)
    dotPlotDrawer.drawDotPlot(
        dsData, args.o+"/hit_downsampled.png", "HIT_MATRIX_DOWNSAMPLED")

    hsps = synteniesIdentification.growing_regions(dsData, args.reward, args.penalty, args.sidePenalty, args.max_hsps,args.gr_th, args.w_size )
    events = synteniesIdentification.detect_events(hsps, args.event_sampling_size)
    IOUtils.storeEvents(events, args.o+"/events.txt")


if __name__ == "__main__":
    start = time.time()
    parser = argparse.ArgumentParser(description='CHROMEISTER ALGORITHME')
    parser.add_argument('-db', type=str, help='database sequence fasta file path',required=True)
    parser.add_argument('-query', required=True,type=str, help='database sequence fasta file path')
    parser.add_argument('-k',type=int, default=32, help='k-mer length')
    parser.add_argument('-hash-key-len', type=int, default=12, help='hash function key length')
    parser.add_argument('-z', type=int, default=4, help='hash function Z parameter')
    parser.add_argument('-dim', type=int, default=1000, help='hit matrix dimension')
    parser.add_argument('-o',type=str, default='./output', help='output directory')
    parser.add_argument('-diag-len', type=int, default=4, help='diagonal length')
    parser.add_argument('-event-sampling-size', type=int, default=4, help='event sampling size')
    parser.add_argument('-reward', type=int, default=6, help='reward for growing regions')
    parser.add_argument('-penalty', type=int, default=15, help='penalty for growing regions')
    parser.add_argument('-sidePenalty', type=int, default=3, help='sidePenalty for growing regions')
    parser.add_argument('-max-hsps', type=int, default=500, help='max HSPS size for growing regions')
    parser.add_argument('-gr-th', type=int, default=5, help='threshold for growing regions')
    parser.add_argument('-w-size', type=int, default=7, help='window size for growing regions')

    args = parser.parse_args()

    main(args)
    print("Total time is "+str(time.time()-start)+"ms")

