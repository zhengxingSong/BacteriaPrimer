import pandas as pd
import re
import sys
import glob

InFileD = sys.argv[1]
OutFile = sys.argv[2]
left = sys.argv[3]
right = sys.argv[4]

print(sys.argv)
Files = glob.glob(InFileD)
df = pd.DataFrame()
for i in Files:
    InFile=open(i,"r")
    dat=InFile.read()
    InFile.close()
    if dat:
        gene_id = re.findall("SEQUENCE_ID=(\S+)\n",dat)[0]
        left_d = dict(re.findall("PRIMER_LEFT_(\d+)_SEQUENCE=(\S+)\n",dat))
        right_d = dict(re.findall("PRIMER_RIGHT_(\d+)_SEQUENCE=(\S+)\n",dat))
        left_TM = dict(re.findall("PRIMER_LEFT_(\d+)_TM=(\S+)\n",dat))
        right_TM = dict(re.findall("PRIMER_RIGHT_(\d+)_TM=(\S+)\n",dat))
        left_GC = dict(re.findall("PRIMER_LEFT_(\d+)_GC_PERCENT=(\S+)\n",dat))
        right_GC = dict(re.findall("PRIMER_RIGHT_(\d+)_GC_PERCENT=(\S+)\n",dat))

        out_list = []
        for i in left_d.keys():
            out_list.append([gene_id,i,left_d[i],right_d[i],left_TM[i],right_TM[i],left_GC[i],right_GC[i]])
        
        df_tmp = pd.DataFrame(out_list,columns=["gene_id","Index","LP","RP","TMl","TMr","GCl","GCr"])
        df = df.append(df_tmp)
    else:
        print(i+"is Failed!")

df["primer"] = df.apply(lambda x:"%s_%s"%(x["gene_id"],x["Index"]),axis=1)
df = df.reset_index().drop(["index"],axis=1)
df.to_csv(OutFile,sep="\t",index=False)

Left_primer = open(left,"w+")
Right_primer = open(right,"w+")
for n in df.index:
    nam = df.loc[n,"primer"]
    Left_primer.writelines(">%s_L\n%s\n"%(nam,df.loc[n,"LP"]))
    Right_primer.writelines(">%s_R\n%s\n"%(nam,df.loc[n,"RP"]))
Left_primer.close()
Right_primer.close()