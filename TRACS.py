import argparse
import sys
import string
import os
import numpy as np
import pandas as pd
import GPRegression as GPR
import gapGP
import genNetwork as genN
import matplotlib.pyplot as plt

def arg_parse():
  parser = argparse.ArgumentParser()
  parser.add_argument('-f', '--file', help="gene expression data (tab-delimited)", required=True)
  parser.add_argument('-tn', '--timenums', type=int, help="the number of time points in gene expresson data", required=True)
  parser.add_argument('-tp', '--timepoints', help="a list of time points in gene expresson data", required=True)
  parser.add_argument('-rn', '--repnums', type=int, help="the number of replicates in gene expression data", required=True)
  parser.add_argument('-og', '--organism', default="NA", help="organism of interest")
  parser.add_argument('-o', '--outdir', help="output directory", required=True)
  parser.add_argument('-m', '--method', default="KM", help="method for clustering ")
  parser.add_argument('-ks', '--kstart', type=int, default=1, help="Starting k for testing the number of clusters")
  parser.add_argument('-ke', '--kend', type=int, default=10, help="Ending k for testing the number of clusters")

  args = vars(parser.parse_args())
  return args

def addNoise(series, reps, s):
  start = 0
  for i in range(len(reps)):
    stds = np.std(series[start:start+reps[i]], axis=0)
    for j in range(reps[i]) :
      series[start+j]+=np.random.normal(np.zeros(stds.shape), s*stds)
    start += reps[i]
  return series

def loadData(inF, reps, delim, norm=False, noise=0.0) :
  with open(inF, "r") as f:
    lines = [l.split(delim) for l in filter(lambda x: len(x)>0, f.readlines())]
    lines = lines[1:]
    names = np.array([Y[0] for Y in lines])
    series = np.array([[float(_.strip()) for _ in Y[1:]] for Y in lines])
    tmp = np.hsplit(series, sum(reps))
    tmp = addNoise(tmp, reps, noise)
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
  return names, series, series_avg

def main(inputs) :
  # Assume that there is one phenotype
  ks = list(range(inputs['kstart'], inputs['kend']+1))
  reps = [inputs['repnums']]
  X = np.array(map(int, inputs['timepoints'].split(',')))
  expFile = inputs['file']
  norm = True
  delim = '\t'

  numS = len(reps)
  interval = 1
  X_pred = np.array([X[0]])
  while (X_pred[len(X_pred)-1]<max(X)) :
    X_pred = np.append(X_pred, X_pred[len(X_pred)-1]+interval)
  if not os.path.exists(inputs['outdir']) :
    os.mkdir(inputs['outdir'])  

  # Load expression data
  names, series, series_avg = loadData(expFile, reps, delim, norm)
  print '# of Genes:', series_avg.shape[0], '\t# of time points:', series_avg.shape[1]

  # Calculate gap statistics
  c, std, c_long, std_long, p, k, labels = gapGP.gap(X, X_pred, series_avg, inputs['method'], inputs['outdir'], numS, ks)
  k = 12
  print 'Optimal K predicted:', k
  print 'Generating a network ...'
  genN.genNetwork(inputs, k) 
  print 'Done'

if __name__ == "__main__" :
  inputs = arg_parse()
  main(inputs)
