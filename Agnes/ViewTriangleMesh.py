from LoadableMesh import *
from MPLMeshViewer import *
from TriangleMeshReader import *
import sys

import argparse
parser = argparse.ArgumentParser(description='Triangle mesh viewer')


parser.add_argument('--f', action='store', default=None)
parser.add_argument('--r', action='store', default=0.1)
parser.add_argument('--fs', action='store', default=14)

args = parser.parse_args()


if args.f==None:
    print('Usage: python3 ./ViewTriangleMesh.py --f=[filename]')
    sys.exit(0)

reader = TriangleMeshReader(args.f)
mesh = reader.getMesh()

viewer = MPLMeshViewer(fontSize=int(args.fs), vertRad=float(args.r))
viewer.show(mesh)
