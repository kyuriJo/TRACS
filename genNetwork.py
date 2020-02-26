import sys, os
import scipy.stats
import numpy as np
import TRACS
import rpvote
import GPRegression as GPR
import DTW
import mkjson

def calSBD(m):  # order of the cluster 
  inputs = []
  k = m.shape[0]
  inputs.append('* '+' '.join(map(str, range(1, k+1))))
  pairs = {}
  for i in range(k) :
    for j in range(i+1, k):
      dist, delay = DTW.SBD(m[i], m[j])
      pairs[(str(i+1), str(j+1))]=(dist, delay)
      #print i+1, j+1, dist, delay
      if delay<0 : inputs.append(str(i+1)+' '+str(j+1))
      elif delay>0 : inputs.append(str(j+1)+' '+str(i+1))
  order = rpvote.vote(inputs)
  return order, pairs

def genNames(org, names):
  gedict = {}
  gidict = {}
  with open('kegg/kegg_'+org+'_geneID.txt','r') as f:
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
  with open('kegg/kegg_'+org+'_pathway_gene.txt', 'r') as f :
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
  N = len(genes)
  # If gene names are official symbols, generate a corresponding array of Entrez ID
  if len(genes & set(names))==0 : 
    names = genNames(org, names)
  # gen matrix of pvals for cluster x pathway
  pathways = np.array(pgDict.keys())
  cn = np.array([len(np.where(lab==i)[0]) for i in range(k)])
  pn = np.array([len(pgDict[p]) for p in pathways])
  cp = np.array([[scipy.stats.hypergeom.sf(len(pgDict[pathways[p]]&set(names[lab==c])), N, pn[p], cn[c]) for p in range(len(pn))] for c in range(k)])

  enRes = {}
  for i in range(k):
    for j in range(i+1, k):
      pcc, pval = scipy.stats.pearsonr(cp[i], cp[j])
      cp_int = np.array([scipy.stats.hypergeom.sf(len(pgDict[pathways[p]]&set(names[lab==i]))+len(pgDict[pathways[p]]&set(names[lab==j])), N, pn[p], cn[i]+cn[j]) for p in range(len(pn))])
      enRes[(str(i+1), str(j+1))]=[pcc, pval, pathways[np.logical_and(np.logical_and(np.logical_and(np.logical_and(cp[i]<0.05, cp[j]<0.05), pn>10), cp_int<cp[i]), cp_int<cp[j])]]
  return enRes

def makeNetFiles(pjname, subname, order, pairs, enRes):
  nodeF = open(pjname+'/'+subname+'_nodes.txt','w')
  edgeF = open(pjname+'/'+subname+'_edges.txt','w')
  for i in range(len(order)-1) :
    if (order[i],order[i+1]) in pairs :
      dist, delay = pairs[(order[i], order[i+1])]
      ens = enRes[(order[i], order[i+1])]
    else :
      dist, delay = pairs[(order[i+1], order[i])]
      ens = enRes[(order[i+1], order[i])]
    edgeF.write('C'+order[i]+'\tC'+order[i+1]+'\t'+str(abs(delay))+'\t'+','.join(ens[2])+'\tNeighbor\n')
  for i in range(len(order)) :
    nodeF.write('C'+str(i+1)+'\t'+subname+'P1cluster'+str(i+1)+'.png'+'\n')
  index = open(pjname+'/index.html', 'w')
  index.write('<!DOCTYPE>\n<html>\n  <head>\n    <title>TRACS</title>\n    <meta name="viewport" content="width=device-width, user-scalable=no, initial-scale=1, maximum-scale=1">\n    <link href="style.css" rel="stylesheet" />\n    <script src="https://cdn.polyfill.io/v2/polyfill.min.js?features=Promise,fetch"></script>\n    <script src="../cytoscape.min.js"></script>\n  </head>\n  <body>\n    <h1>TRACS</h1>\n    <div id="cy"></div>\n    <script src="code.js"></script>\n  </body>\n</html>\n')
  js = open(pjname+'/code.js', 'w')
  jscode = """fetch('jsonfile', {mode: 'no-cors'})
.then(function(res) {
  return res.json()
})
  .then(function(data) {
    var cy = window.cy = cytoscape({
      container: document.getElementById('cy'),

      boxSelectionEnabled: false,
      autounselectify: true,

      layout: {
        name: 'preset'
      },

      style: [
        {
          selector: 'node',
          css: {
            'height': 40,
            'width': 40,
            'label': 'data(id)',
            'border-width': 2,
            'border-color': 'black',
            'font-size': '10px'
          }
        },
        {
          selector: 'edge',
          css: {
            'width': 1,
            'font-size':'6px',
            'label': 'data(label)',
            'curve-style': 'bezier',
            'target-arrow-shape': 'triangle',
            'text-wrap':'wrap',
            'text-max-width': 30
          }
        }
      ],

      elements: data
    });
  });"""
  jscode = jscode.replace('jsonfile', subname+'_data.json') 
  js.write(jscode)
  return

def genNetwork(inputs, k) :
  reps = [inputs['repnums']]
  X = np.array(map(int, inputs['timepoints'].split(',')))
  expFile = inputs['file']
  pjname = inputs['outdir']
  org = inputs['organism']
  norm = True
  delim = '\t'
  cond = 1

  numS = len(reps)
  interval = 1
  X_pred = np.array([X[0]])
  while (X_pred[len(X_pred)-1]<max(X)) :
    X_pred = np.append(X_pred, X_pred[len(X_pred)-1]+interval)

  names, series, series_avg = TRACS.loadData(expFile, reps, delim, norm)
  subname = 'K'+'{:02d}'.format(k)
  labels = np.genfromtxt(pjname+'/'+subname+'_labels.txt', delimiter='\t', dtype=np.str, usecols=[0])
  labels = labels.astype(np.int)
 
  minmax = (np.floor(np.min(series_avg)), np.ceil(np.max(series_avg)))
  temp = np.hsplit(series_avg, len(reps))
  means = []
  for i in range(k) :
    temp_i = temp[cond-1][labels==i]
    pred_m_x, pred_v_x, pred_m, pred_v, params = GPR.multipleReg(X, X_pred, temp_i, pjname, 'rand', cond, i+1, minmax)
    means.append(pred_m)
  means = np.array(means)
  order, pairs = calSBD(means)
  enRes = calEnrich(pjname, org, k, names, labels)
  makeNetFiles(pjname, subname, order, pairs, enRes)
  mkjson.visNetwork(pjname, subname)
  return
