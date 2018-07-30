# This program is to identify cis Proline residue in protein
# You need a PDB file, and the program out put the cis Proline
# residue with residue number and chain id
# Usage: python cis-pro-in-pdb.py 3ui4     (where the pdb file is 3ui4.pdb)

import os
import string
import sys
import math
from sys import argv

try: 
   input = sys.argv[1]
except:
   print "no Input"

pdb = open(input+'.pdb', "r").readlines()

chains = []
resname = []
resnum = []
atoms = []
coor_x = []
coor_y = []
coor_z = []

total_atom = 0
count = 0
rpi = 180/3.14159

for line in range(len(pdb)):
    N = pdb[line]
    id1 = N[:4]
    if id1 == 'ATOM':
       chains.append(N[21])
       resname.append(N[17:20])
       resnum.append(N[23:26])
       atoms.append(N[12:16])
       coor_x.append(N[31:38])
       coor_y.append(N[39:46])
       coor_z.append(N[47:54])
       total_atom = total_atom + 1

for i in range(1,total_atom):
    if resname[i] == 'PRO' and atoms[i] == ' N  ':
       pro_resnum = resnum[i]
       pro_chain = chains[i]
       pro_nx = coor_x[i]
       pro_ny = coor_y[i]
       pro_nz = coor_z[i]

    if resname[i] == 'PRO' and atoms[i] == ' CA ' and chains[i] == pro_chain:
       pro_cax = coor_x[i]
       pro_cay = coor_y[i]
       pro_caz = coor_z[i]

       pro_resnum_a = int(pro_resnum) - 1
       for j in range(1,total_atom):
           if pro_resnum_a == int(resnum[j]) and atoms[j] == ' CA ' and chains[j] == pro_chain:
              xpro_cax = coor_x[j]
              xpro_cay = coor_y[j]
              xpro_caz = coor_z[j]
           if pro_resnum_a == int(resnum[j]) and atoms[j] == ' C  ' and chains[j] == pro_chain:
              xpro_cx = coor_x[j]
              xpro_cy = coor_y[j]
              xpro_cz = coor_z[j]

       xij = float(xpro_cax) - float(xpro_cx)
       yij = float(xpro_cay) - float(xpro_cy)
       zij = float(xpro_caz) - float(xpro_cz)
       xkj = float(xpro_cx) - float(pro_nx)
       ykj = float(xpro_cy) - float(pro_ny)
       zkj = float(xpro_cz) - float(pro_nz)
       xkl = float(pro_nx) - float(pro_cax)
       ykl = float(pro_ny) - float(pro_cay)
       zkl = float(pro_nz) - float(pro_caz)

       dx = -yij*zkj+zij*ykj
       dy = -zij*xkj+xij*zkj
       dz = -xij*ykj+yij*xkj
       gx = zkj*ykl-ykj*zkl
       gy = xkj*zkl-zkj*xkl
       gz = ykj*xkl-xkj*ykl

       norm1 = math.sqrt(dx*dx+dy*dy+dz*dz)
       norm2 = math.sqrt(gx*gx+gy*gy+gz*gz)

       dangle = (dx*gx+dy*gy+dz*gz)/(norm1*norm2)

       if dangle > 1.0:
          dangle = 1.0
       
       if dangle < -1.0:
          dangle = -1.0

       dangle = math.acos(dangle)*rpi
       ssum = xkj*(-dz*gy+dy*gz)+ykj*(-dx*gz+dz*gx)+zkj*(-dy*gx+dx*gy)
       if ssum < 0.0:
          dangle = -dangle

       if dangle >= -50.0 and dangle <= 50.0:
          print "cis-PRO", pro_resnum, pro_chain




