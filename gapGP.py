# (c) 2013 Mikael Vejdemo-Johansson
# BSD License
#
# SciPy function to compute the gap statistic for evaluating k-means clustering.
# Gap statistic defined in
# Tibshirani, Walther, Hastie:
#  Estimating the number of clusters in a data set via the gap statistic
#  J. R. Statist. Soc. B (2001) 63, Part 2, pp 411-423

import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.cluster.vq
import scipy.stats
import scipy.spatial.distance
import pickle
import GPRegression as GPR
import Clustering
import sklearn.metrics
from scipy.interpolate import interp1d
from sklearn.cluster import KMeans

def dst(t1, t2) :
  return sklearn.metrics.mean_squared_error(t1, t2)*t1.shape[0]

def calDist(data, m, std):
  return sum([scipy.stats.norm.pdf(p, m[i], std[i]) for (i, p) in enumerate(data)])

def calDist_l(data, x, xl, ml, stdl): # When used in GPR (already split data)
  f = interp1d(x, data)
  data = f(xl)
  return sum([scipy.stats.norm.pdf(p, ml[i], stdl[i]) for (i, p) in enumerate(data)])

def calDist_long(data, x, xl, ml, stdl, numS):
  datal = []
  datas = np.hsplit(data, numS)
  for s in range(numS) :
    f = interp1d(x, datas[s])
    datal.append(f(xl))
  datal = np.concatenate(datal) # Assume len(xl)*2 = len(ml)
  return sum([scipy.stats.norm.pdf(p, ml[i], stdl[i]) for (i, p) in enumerate(datal)])

def optK(ks, gaps, stds) :
  opt_i = len(ks)-1 
  for (i, k) in enumerate(ks[:len(ks)-1]) :
    if (gaps[i]>=(gaps[i+1]-stds[i+1])) :
      opt_i = i
      break
  return opt_i

def optK2(ks, gap, std) :
  max_i = np.argmax(gap)
  opt_i = max_i
  cur = max_i-1
  while (cur>=0) :
    if (gap[max_i]-std[max_i]<=gap[cur]) and (gap[max_i]+std[max_i]>=gap[cur]) :
      opt_i = cur
    cur -=1
  return opt_i

def calCent(X, X_pred, data, pjname, subname, labels, numS, k, bound=(0.1,0.1)) :
  minmax = (np.floor(np.min(data)), np.ceil(np.max(data)))
  data_s = np.hsplit(data, numS)
  kmc = np.zeros((numS, k, X.shape[0]))
  kmstd = np.zeros((numS, k, X.shape[0]))
  kmc_long = np.zeros((numS, k, X_pred.shape[0]))
  kmstd_long = np.zeros((numS, k, X_pred.shape[0]))
  ps = np.zeros((numS, k, 3)) 
  for i in range(k) : # labels starts from 0
    for cond in range(numS) :
      temp_s = data_s[cond][labels==i]
      if temp_s.shape[0]==0 :
        m = np.zeros(X.shape[0])+5.
        v = np.zeros(X.shape[0])+1.0
        ml = np.zeros(X_pred.shape[0])+5.
        vl = np.zeros(X_pred.shape[0])+1.0
        p = [0,0,0] 
      else :
        m, v, ml, vl, p = GPR.multipleReg(X, X_pred, temp_s, pjname, subname, cond+1, i+1, minmax, bound)
      kmc[cond][i]=m 
      kmstd[cond][i]=v
      kmc_long[cond][i]=ml
      kmstd_long[cond][i]=vl
      ps[cond][i]=p
  kmc = np.concatenate(kmc, axis=1)
  kmstd = np.concatenate(kmstd, axis=1)
  kmc_long = np.concatenate(kmc_long, axis=1)
  kmstd_long = np.concatenate(kmstd_long, axis=1)
  ps = np.concatenate(ps, axis=1)
  return kmc, kmstd, kmc_long, kmstd_long, ps

def gap(X, X_pred, data, met, pjname, numS, ks, bound=(0.1,0.1)):
  # Generate reference dist of the origianl dataset
  nrefs = 10
  shape = data.shape
  refs = None
  if refs==None:
    tops = data.max(axis=0)
    bots = data.min(axis=0)
    dists = scipy.matrix(scipy.diag(tops-bots))
    rands = scipy.random.random_sample(size=(shape[0],shape[1],nrefs))
    for i in range(nrefs):
      rands[:,:,i] = rands[:,:,i]*dists+bots
  else:
    rands = refs
 
  # Calculate gap for each k
  res = []
  gaps = scipy.zeros((len(ks),))
  stds = np.zeros((len(ks),))
  gaps_d = scipy.zeros((len(ks),))
  stds_d = np.zeros((len(ks),))
  for (i,k) in enumerate(ks):
    kmc, kml = Clustering.clustering(k, X, data, met)
    # Added procedure : calculate GP mean and var rather than cluster centers
    # Calculate distance by likelihood of Gaussian rather than Euclidean
    disp_d = sum([dst(data[m,:],kmc[kml[m],:]) for m in range(shape[0])])
    subname = 'K'+'{:02d}'.format(k)
    np.savetxt(pjname+'/'+subname+'_labels.txt', kml, fmt='%s', delimiter='\n')
    kmc, kmstd, kmc_long, kmstd_long, p = calCent(X, X_pred, data, pjname, subname, kml, numS, k, bound)
    res.append((kmc, kmstd, kmc_long, kmstd_long, p, kml))
    disp = sum([calDist_long(data[m,:], X, X_pred, kmc_long[kml[m],:], kmstd_long[kml[m],:], numS) for m in range(shape[0])])
    # Calculate reference gap
    refdisps = scipy.zeros((rands.shape[2],))
    refdisps_d = scipy.zeros((rands.shape[2],))
    for j in range(rands.shape[2]):
      kmc, kml = Clustering.clustering(k, X, rands[:,:,j], met)
      refdisps_d[j] = sum([dst(rands[m,:,j],kmc[kml[m],:]) for m in range(shape[0])])
      subname = 'R'+'{:02d}'.format(k)
      kmc, kmstd, kmc_long, kmstd_long, p = calCent(X, X_pred, rands[:,:,j], pjname, subname, kml, numS, k, bound)
      refdisps[j] = sum([calDist_long(rands[m,:,j], X, X_pred, kmc_long[kml[m],:], kmstd_long[kml[m],:], numS) for m in range(shape[0])])
    gaps[i] = scipy.log(disp)-scipy.mean(scipy.log(refdisps))
    stds[i] = np.std(scipy.log(refdisps))*np.sqrt(1+1/float(nrefs))
    gaps_d[i] = scipy.mean(scipy.log(refdisps_d))-scipy.log(disp_d)
    stds_d[i] = np.std(scipy.log(refdisps_d))*np.sqrt(1+1/float(nrefs))
    print 'Gap(GP) for', k, 'is', scipy.log(disp), '-', scipy.mean(scipy.log(refdisps)), '=', gaps[i]
    print 'Gap(ED) for', k, 'is', scipy.log(disp_d), '-', scipy.mean(scipy.log(refdisps_d)), '=', gaps_d[i]

  # Find the optimal k by std of log(refdisps)
  opt_i = optK2(ks, gaps, stds)
  opt_id = optK2(ks, gaps_d, stds_d)
  # Visualize gap statistics
  plt.errorbar(ks, gaps, yerr=stds)
  plt.errorbar(ks, gaps_d, yerr=stds_d)
  plt.savefig(pjname+'/GapStatistics.png')
  plt.close()
  pickle.dump((gaps, stds, gaps_d, stds_d), open(pjname+'/gaps.dump', 'w'))
  c, std, c_long, std_long, p, labels = res[opt_i]

  return c, std, c_long, std_long, p, ks[opt_i], labels
