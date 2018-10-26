import sys, os
import scipy.stats
import numpy as np
import TiClNet
import rpvote
import GPRegression as GPR
import DTW

def calSBD(m):  # order of the cluster / 
  inputs = []
  k = m.shape[0]
  inputs.append('* '+' '.join(map(str, range(1, k+1))))
  for i in range(k) :
    for j in range(i+1, k):
      dist, delay = DTW.SBD(m[i], m[j])
      print i+1, j+1, dist, delay
      if delay<0 : inputs.append(str(i+1)+' '+str(j+1))
      elif delay>0 : inputs.append(str(j+1)+' '+str(i+1))     
  order = rpvote.vote(inputs)
  print order
  return order

def genNames(org, names):
  gedict = {}
  gidict = {}
  with open('kegg_'+org+'_geneID.txt','r') as f:
    for l in f.readlines() :
      tp = l.split('\t')
      e = tp[0][4:]
      gs = tp[1].split(';')[0].split(', ')
      for i, g in enumerate(gs) : 
        if g not in gedict : 
          gedict[g]=[]
          gidict[g]=np.array([])
        gedict[g].append(e)
        gidict[g] = np.append(gidict[g], i)
  new_names = []
  for n in names:
    new_names.append(gedict[n][np.argmin(gidict[n])] if n in gedict else '00000')
  return np.array(new_names)

def genPGdict(org):
  pgdict = {}
  genes = set()
  with open('kegg_'+org+'_pathway_gene.txt', 'r') as f :
    for l in f.readlines() :
      tp = l.split()
      p = tp[0][5:]
      g = tp[1][4:]
      if p not in pgdict :
        pgdict[p]=set()
      pgdict[p].add(g)
      genes.add(g)
  return pgdict, genes

def calEnrich(pjname, org, k, names, lab):
  pgDict, genes = genPGdict(org)
  print len(genes)
  # If gene names are official symbols, generate a corresponding array of Entrez ID
  if len(genes & set(names))==0 : 
    names = genNames(org, names)
  print names
  # gen matrix of pvals for cluster x pathway
  pathways = np.array(pgDict.keys())
  cn = np.array([len(np.where(lab==i)[0]) for i in range(k)])
  pn = np.array([len(pgDict[p]) for p in pathways])
  cp = np.array([[scipy.stats.hypergeom.sf(len(pgDict[pathways[p]]&set(names[lab==c])), N, pn[p], cn[c]) for p in range(len(pn))] for c in range(k)])
  #print cn
  #print [[len(pgDict[pathways[p]]&set(names[lab==c])) for p in range(len(pn))] for c in range(k)]
  for i in range(k):
    for j in range(i+1, k):
      pcc, pval = scipy.stats.pearsonr(cp[i], cp[j])
      cp_int = np.array([scipy.stats.hypergeom.sf(len(pgDict[pathways[p]]&set(names[lab==i]))+len(pgDict[pathways[p]]&set(names[lab==j])), N, pn[p], cn[i]+cn[j]) for p in range(len(pn))])
      print i+1, j+1, pcc, pval, pathways[np.logical_and(np.logical_and(np.logical_and(np.logical_and(cp[i]<0.05, cp[j]<0.05), pn>10), cp_int<cp[i]), cp_int<cp[j])]
  return

def genNetwork(pjname, conf, k, cond, org) :
  # Configuration - change according to the dataset
  ks = conf[0]
  reps = conf[1]
  X = np.array(conf[2])
  expFile = conf[3]
  k_orig = conf[4]
  norm = conf[5]
  delim = conf[6]
  cols = conf[7:] # cols = (name, probe, lab, series(0:4) )

  numS = len(reps)
  interval = 1
  X_pred = np.array([X[0]])
  while (X_pred[len(X_pred)-1]<max(X)) :
    X_pred = np.append(X_pred, X_pred[len(X_pred)-1]+interval)

  labels, names, probes, series, series_avg = TiClNet.loadData(expFile, reps, delim, cols, norm)
  subname = 'K'+'{:02d}'.format(k)
  labels = np.genfromtxt(pjname+'/'+subname+'_labels.txt', delimiter='\t', dtype=np.str, usecols=[0])
  labels = labels.astype(np.int)
 
  minmax = (np.floor(np.min(series_avg)), np.ceil(np.max(series_avg)))
  temp = np.hsplit(series_avg, len(reps))
  means = []
  for i in range(k) :
    temp_i = temp[cond-1][labels==i]
    pred_m_X, pred_v_X, pred_m, pred_v, params = GPR.multipleReg(X, X_pred, temp_i, pjname, 'test', cond, i+1, minmax)
    means.append(pred_m)
  means = np.array(means)
  order = calSBD(means)
  calEnrich(pjname, org, k, names, labels)
  return
