import re
import sys

InFa = sys.argv[1]
temp = sys.argv[2]
OutDir=sys.argv[3]
print(sys.argv)

InFile = open(InFa,"r")
dat = InFile.read()
InFile.close()

fa2d = dict(re.findall(">(\S+)\s\S+\n([A-Z]+)\n",dat))

TempFile = open(temp,"r")
Template = TempFile.read()
TempFile.close()

for k,v in fa2d.items():
    Primer = open(OutDir+k+".txt","w+")
    Out = Template.replace("Seq_ID",k)
    Out = Out.replace("Seq_InFo",v)
    Primer.write(Out)
    Primer.close()