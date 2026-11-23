#!/usr/bin/python
###!/local/python-2.5.1-login-gnu/bin/python
# Some topography data has zero values in land locations
# This small script puts "huge" values in those places -
# - need for interaction with standard hycom partit
import struct 
import sys
import shutil

if len(sys.argv)!= 3 :
   print "addhuge requires 2 arguments - infile outfile"
   print "infile and outfile will have .a and .b appended" 
   sys.exit(-1)
(tmp,infile, outfile)=sys.argv;
huge=2**100

# Get grid size from regional.grid.b
gridinb = open('regional.grid.b',"r")
line=gridinb.readline()
idm=int((line.split())[0])
line=gridinb.readline()
jdm=int((line.split())[0])
gridinb.close()

n2drec = ((idm*jdm+4095)/4096)*4096
print 'idm, jdm and n2drec: ',idm , jdm, n2drec
amax=-1e10
amin= 1e10
depthina  = open(infile + '.a',"rb")
depthouta = open(outfile + '.a',"wb")
counthuge=0
countzero=0
countocn=0
for i in range(n2drec) :
   line=depthina.read(4)          # four byte real
   [ b ]=struct.unpack(">f",line)
   a=b
   if (b > 1e20 or b<0.001 ) :
      a=huge
      if (b > 1e20 ) :
         counthuge=counthuge+1
      if (b < 0.001 ) :
         countzero=countzero+1
   else :
      countocn=countocn+1
      amax=max(a,amax)
      amin=min(a,amin)
   lineout=struct.pack(">f",a)
   depthouta.write(lineout)
depthina.close()
depthouta.close()

# Copy .b file 
depthinb = open(infile + '.b',"r")
line=depthinb.readline()
depthinb.close()

depthoutb = open(outfile + '.b',"w")
depthoutb.writelines(line.strip() + "(modified by addhuge.py)\n");
depthoutb.writelines("\n")
depthoutb.writelines("\n")
depthoutb.writelines("\n")
depthoutb.writelines("\n")
depthoutb.writelines("min,max depth = %11.3f%11.3f\n" % (amin , amax))
depthoutb.close()

print "Number of huge  values in original file:" , counthuge
print "Number of zero  values in original file:" , countzero
print "Number of ocean values in original file:" , countocn
print "Zero values reset to huge in original file"
print "New file name is :",outfile+'.[ab]'

