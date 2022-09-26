import pandas as pd
import re
import sys

def USAGE():
    print("%s InFile RefFile OutFile"%sys.argv[0])


def Fa2dic(data):
    dic = {}
    anno_dic = {}
    dat_pre = data.split(">lcl|")
    dat2 = [x.split("\n") for x in dat_pre if x]
    for i in dat2:
        k = i[0].split(" ")[0]
        anno = re.findall("(?<=\[protein=)([A-Za-z\s]*)]",i[0])
        dic[k] = "".join(i[1:])
        anno_dic[k] = anno[0] if anno else "None"
    return(dic,anno_dic)


def main():   
    InFile =sys.argv[1]
    RefFile = sys.argv[2]
    OutFile=sys.argv[3]

    df = pd.read_table(InFile,header=None,names=["ID","start","end"])
    df = df.sort_values("ID").reset_index().drop("index",axis=1)


    Inf = open(RefFile,"r")
    dat = Inf.read()
    dic_seq,dic_anno = Fa2dic(dat)

    tmp = ""
    with open(OutFile,"w+") as out:
        for i in df.index:
            qualifier = df.loc[i,"ID"]
            start = int(df.loc[i,"start"])
            end = int(df.loc[i,"end"])
            anno = dic_anno[qualifier]
            if qualifier !=tmp:
                n =1
                tmp = qualifier
            else:
                n+=1
            out.writelines(">%s_%d\t%s\n%s\n"%(qualifier,n,anno,dic_seq[qualifier][(start-1):end]))
            print("%s is done!"%qualifier)
    Inf.close()

try:
    main()
except:
    USAGE()