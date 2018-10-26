import numpy as np
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri as numpy2ri
numpy2ri.activate()
dtw = importr("dtwclust")

def SBD(a, b):
  a = a+np.abs(np.min(a))
  b = b+np.abs(np.min(b))
  l = a.shape[0]
  dist = dtw.SBD(a, b, znorm=False)
  cc = dtw.NCCc(a, b)
  cc = 1-np.array(cc)
  delay = np.argmin(cc)-(l-1)
  return dist[0][0], delay

def KShape(k, data):
  res = dtw.tsclust(data, type='partitional', k=k, distance='sbd', centroid='shape')
  label = np.array(res.slots['cluster'])
  label -= 1 
  return label
