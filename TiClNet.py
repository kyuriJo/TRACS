import argparse
import sys
import string
import os
import numpy as np
import GPRegression as GPR
import gapGP
import genNetwork as genN
import matplotlib.pyplot as plt
import Clustering

# General assumption: 
# numS = 2 ?
# x for each condition is same

class inputC :
        args = {}
        organism = ""
        expFile = ""
        outDir = ""
        numTP = 0
        numRe = 0
        DEGCut = 0
        clusterCut = 0
        timeLag = 0

def arg_parse():
        inputs = inputC()
        parser = argparse.ArgumentParser()
        parser.add_argument('-f', '--file', help="text file of expression value", required=True)
        parser.add_argument('-o', '--outdir', help="output directory", required=True)
        parser.add_argument('-s', '--organism', help="organism code for KEGG", required=True)
        parser.add_argument('-nt', '--numTP', type=int, help="number of time points", required=True)
        parser.add_argument('-nr', '--numRe', type=int, help="number of replicates", required=True)
        parser.add_argument('-c2', '--clustercut', type=float, default=1.0, help="cutoff for clustering", required=True)
        parser.add_argument('-c1', '--degcut', type=float, default=1.0, help="cutoff for finding DEGs", required=True)
        parser.add_argument('-l', '--timelag', type=float, default=1.0, help="time-lag factor", required=True)

        args = vars(parser.parse_args())
        inputs.args = args
        inputs.expFile = args['file']
        inputs.outDir = args['outdir']
        inputs.organism = args['organism']
        inputs.numTP = args['numTP']
        inputs.numRe = args['numRe']
        inputs.clusterCut = args['clustercut']
        inputs.DEGCut = args['degcut']
        inputs.timeLag = args['timelag']
        return inputs
 
def addNoise(series, reps, s):
  start = 0
  for i in range(len(reps)):
    stds = np.std(series[start:start+reps[i]], axis=0)
    for j in range(reps[i]) :
      series[start+j]+=np.random.normal(np.zeros(stds.shape), s*stds)
    start += reps[i]
  return series

def loadData(inF, reps, delim, cols, norm=False, noise=0.0) :
  with open(inF, "r") as f:
    lines = [l.split(delim) for l in filter(lambda x: len(x)>0, f.readlines())]
    labels = np.array([int(float(Y[cols[2]])) for Y in lines])
    names = np.array([Y[cols[0]] for Y in lines])
    probes = np.array([Y[cols[1]] for Y in lines])
    series = np.array([[float(_.strip()) for _ in Y[cols[3][0]:cols[3][1]]] for Y in lines])
    tmp = np.hsplit(series, sum(reps))
    tmp = addNoise(tmp, reps, noise)

    # labels should be [0,k-1]
    if np.min(labels)==1 : labels-=1
    # Normalization
    if (norm):
      for i in range(sum(reps)) :
        tmp[i] = np.apply_along_axis(lambda x: (x-np.average(x))/np.std(x), 1, tmp[i])
    series_avg = []
    start = 0
    for i in range(len(reps)) :
      series_avg.append(np.average(tmp[start:start+reps[i]], axis=0))
      start += reps[i]
    series_avg = np.concatenate(np.array(series_avg), axis=1)
  return labels, names, probes, series, series_avg

def main(inputs) :
  # Configuration - change according to the dataset
  ks = conf[0]
  reps = conf[1]
  X = np.array(conf[2])
  expFile = conf[3]
  k = conf[4]
  norm = conf[5]
  delim = conf[6]
  cols = conf[7:] # cols = (name, probe, lab, series(0:4) )

  numS = len(reps)
  interval = 1
  X_pred = np.array([X[0]])
  while (X_pred[len(X_pred)-1]<max(X)) :
    X_pred = np.append(X_pred, X_pred[len(X_pred)-1]+interval)
  if not os.path.exists(inputs.pjname) :
    os.mkdir(inputs.pjname)  

  # Load expression data
  labels, names, probes, series, series_avg = loadData(expFile, reps, delim, cols, norm)
  print 'series', series.shape
  print 'series_avg', series_avg.shape
  print np.min(series_avg), np.max(series_avg)

  # assume there is one condition
  c, std, c_long, std_long, p, k, labels = gapGP.gap(X, X_pred, series_avg, met, inputs.pjname, numS, ks)
  print 'optK', k
  genN.genNetwork(inputs) 

if __name__ == "__main__" :
  inputs = arg_parse()
  main(inputs)
