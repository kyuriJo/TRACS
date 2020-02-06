import numpy as np
import json

def visNetwork(pjname, subname) :
  #Read node.txt
  Ndata = []
  Nlines = open(pjname+'/'+subname+'_nodes.txt','r').readlines()
  for line in Nlines:
      Ndata.append(line.replace('\n','').split('\t'))
  #Read edge.txt
  Edata = []
  Elines = open(pjname+'/'+subname+'_edges.txt','r').readlines()
  for line in Elines:
      Edata.append(line.replace('\n','').split('\t'))
  Delay = np.array(Edata)[:,2]

  #define position of node
  Pos = {}
  Pos[Edata[0][0]] = np.array([0,0], dtype=np.int32)
  if Delay[0] == '0':
      for i, data in enumerate(Edata):
          if data[4]=='1' : continue
          start, end = data[:2]
          if i==0:
              Pos[end] = Pos[start] + np.array([0,100])
          else:
              Pos[end] = Pos[start] + np.array([100,0])
              Pos[end][1] = 50
  elif '0' in Delay and Delay[0] != '0':
      Pos[Edata[0][0]] = np.array([0,50], dtype=np.int32)
      for i, data in enumerate(Edata):
          if data[4]=='1' : continue
          start, end = data[:2]
          if Delay[i] == '0':
              Pos[start][1] = 0
              Pos[end] = Pos[start] + np.array([0,100])
          else:
              Pos[end] = Pos[start] + np.array([100,0])
              Pos[end][1] = 50
  else:
      for data in Edata:
          if data[4]=='1' : continue
          start, end = data[:2]
          Pos[end] = Pos[start] + np.array([100,0])

  #combine position data with node data
  for key in Pos.keys():
      for i, data in enumerate(Ndata):
          if key == data[0]:
              Ndata[i] = Ndata[i] + [val for val in Pos[key]]
  datalist = []

  for i in range(len(Ndata)):
      datalist.append({})
      datalist[i]['data'] = {'id':Ndata[i][0]}
      datalist[i]['position'] = {'x':int(Ndata[i][2]), 'y':int(Ndata[i][3])}
      datalist[i]['style'] = {'background-image':Ndata[i][1],'background-fit':'cover'}
      datalist[i]['group'] = 'nodes'

  for i in range(len(Edata)):
      datalist.append({})
      datalist[len(Ndata)+i]['data']={'id':Edata[i][0]+Edata[i][1], 'source':Edata[i][0], 'target':Edata[i][1],'label': Edata[i][3], 'pathways':Edata[i][3]}
      datalist[len(Ndata)+i]['group'] = 'edges'

  #make data.json
  with open(pjname+'/'+subname+'_data.json','w') as outfile:
      json.dump(datalist, outfile, indent=2)
