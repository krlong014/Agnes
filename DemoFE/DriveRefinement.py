import os

meshName = 'triExample'
maxLevel = 25

for i in range(maxLevel):
    a = 1.0/2.0**(i+1)
    cmd = 'triangle -rqe -a%-20.15f %s.%d' % (a, meshName, i)
    print(cmd)
    os.system(cmd)
