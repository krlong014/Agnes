import pprint

files = (
  (3, 'strang4'),
  (4, 'strang5'),
  (5, 'strang7'),
  (6, 'strang9')
)

data = {}

for f in files:
  order = f[0]
  name = f[1]

  wFile = open('{}_w.txt'.format(name), 'r')

  w = []
  x = []
  
  with open('{}_w.txt'.format(name), 'r') as wFile:
    for line in wFile:
      L = line.strip(' \r\n')
      if len(L)==0:
        continue
      w.append(float(L))

  with open('{}_x.txt'.format(name), 'r') as xFile:
    for line in xFile:
      L = line.strip(' \r\n')
      if len(L)==0:
        continue
      xy = L.split()
      x.append((float(xy[0]), float(xy[1])))

  data[order] = {'x' : x, 'w' : w}

pprint.pp(data) 