#!/local/python-2.5.1-login-gnu/bin/python
# Script converts "nersc" bathy file to standard hycom bathy file
import struct 
import sys
import shutil

if len(sys.argv)!= 5 :
   print "depthtoregional requires 4 arguments - infile outfile idm  jdm"
   print "outfile will have .a and .b appended" 
   sys.exit(-1)
(tmp,infile, outfile, idm,jdm)=sys.argv;
idm=int(idm)
jdm=int(jdm)
huge=2**100
n2drec = ((idm*jdm+4095)/4096)*4096

# Read binary data from infile + .a , put huge values
# where land is present, and dump to outfile + .a
amax=-1e10
amin= 1e10
myfileina  = open(infile ,"rb")
myfileouta = open(outfile + '.b',"wb")
line=myfileina.read(4) # 1st 4 bytes are fortran header
counthuge=0
countzero=0
countocn=0
for i in range(1,n2drec+1) :
   if (i<=idm*jdm) :
      line=myfileina.read(8)
      [ a ]=struct.unpack(">d",line)
      if (a<0.001 or a > huge/2. ) :
         if (a > huge/2  ) :
            counthuge=counthuge+1
         if (a < 0.001 ) :
            countzero=countzero+1
         a=huge
      else :
         countocn=countocn+1
         amax=max(a,amax)
         amin=min(a,amin)
   else :
      a=huge
      myfileouta.write(struct.pack(">f",a))
myfileina.close()
myfileouta.close()

# Make .b file 
myfileoutb = open(outfile + '.b',"w")
myfileoutb.writelines("Depth Grid converted from " + infile + " by depthtoregional.py\n")
myfileoutb.writelines("\n")
myfileoutb.writelines("\n")
myfileoutb.writelines("\n")
myfileoutb.writelines("\n")
myfileoutb.writelines("min,max depth = %11.3f%11.3f \n" % (amin , amax))
myfileoutb.close()


print "Number of huge  values in original file:" , counthuge
print "Number of zero  values in original file:" , countzero
print "Number of ocean values in original file:" , countocn
print "Zero values reset to huge in original file"
print "New file name is :",outfile+'.[ab]'
