import numpy as np
import GPC, DTW
from sklearn.cluster import KMeans
from sklearn.cluster import AgglomerativeClustering as AC

def clustering(k, x, series_avg, met):
  if met == "KM" :
    res = KMeans(k).fit(series_avg)
    return res.cluster_centers_, res.labels_
  elif met == "AC" :
    res = AC(n_clusters=k, linkage="complete").fit(series_avg)
    cent = np.array([np.mean(series_avg[res.labels_==i], axis=0) for i in range(k)])
    return cent, res.labels_
  elif met == "GP" :
    label = GPC.run(k, x, 1, series_avg)
    cent = np.array([np.mean(series_avg[label==i], axis=0) for i in range(k)])
    #print ','.join(map(str, label))
    return cent, label
  elif met == "KS" :
    label = DTW.KShape(k, series_avg)
    cent = np.array([np.mean(series_avg[label==i], axis=0) for i in range(k)])
    return cent, label

