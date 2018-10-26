import sys
import string
import os
import numpy as np
import GPRegression as GPR
import gapGP
import matplotlib.pyplot as plt
import Clustering

# General assumption: 
# numS = 2 ?
# x for each condition is same

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

def run(pjname, conf, met) :
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
  if not os.path.exists(pjname) :
    os.mkdir(pjname)  

  # Load expression data
  labels, names, probes, series, series_avg = loadData(expFile, reps, delim, cols, norm)
  print 'series', series.shape
  print 'series_avg', series_avg.shape
  print np.min(series_avg), np.max(series_avg)

  # assume there is one condition
  c, std, c_long, std_long, p, k, labels = gapGP.gap(X, X_pred, series_avg, met, pjname, numS, ks)
  print 'optK', k 

def main(pjname, met) :
  # ks, reps, X, expFile, k, norm, col of name, probe, lab, series(0:4)
  conf22875 = [range(1,21), [3,1], [0,8,16,24,48,96], 'GSE22875/GSE22875.txt', ]
  conf69667 = [range(1,11), [2], [0,6,12,24,36,48,72,96], 'GSE69667/GSE69667.txt', 3, True, '\t', 0, 19, 18, (2, 18)]
  confCho = [range(2, 21), [1], range(0, 170, 10), 'Cho/Cho_expr.txt', 5, True, '\t', 0, 0, 1, (2, 19)]
  confCho2 = [range(1, 11), [1], [0, 10, 20, 40, 80, 160], 'Cho/Cho_expr2.txt', 5, True, '\t', 0, 0, 1, (2, 8)]
  confLeaf = [range(1, 11), [1], range(427), 'UCR/OSULeaf/OSULeaf', 6, False, ',', 0, 0, 0, (1, 428)]
  confZNF =  [range(1, 11), [3], [0,1,2,3,4,5,6,8,12,24], 'ZNF217/GSE78169_expr_DEG.txt', 3, True, '\t', 0, 1, 2, (4, 34)]
  run('GSE69667_AC_C_1', conf69667, met)
  #run(pjname, confCho, met)
  #run('Cho_T6', confCho2)
  #run('UCR_Leaf', confLeaf)
  #run('Cho_T6-51', confCho2)
  #run(pjname, confZNF, met)

if __name__ == "__main__" :
  main(sys.argv[1], sys.argv[2])
