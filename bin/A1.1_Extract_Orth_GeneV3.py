import pandas as pd
import sys

def USAGE():
    print("%s InFile OutFile RefFile Num Identity MInLen MaxLen"%sys.argv[0])

def Construct_Region(dat,id,st,end):
    if dat[st] < dat[end]:
        qst,qend = dat[st],dat[end]  
    else:
        qst,qend =  dat[end],dat[st]
    return("%s:%d-%d"%(dat[id],qst,qend))

if len(sys.argv) <3:
    USAGE()
    sys.exit(1)

def Construct_Region(dat,id,st,end):
    if dat[st] < dat[end]:
        qst,qend = dat[st],dat[end]  
    else:
        qst,qend =  dat[end],dat[st]
    return("%s:%d-%d"%(dat[id],qst,qend))

def Cds2Group(x,search_d):
    #cds_nam = "_".join(x.split("_")[0:2])
    group = search_d.get(x,"Unknown")
    return(group)



print(sys.argv)

InFile = sys.argv[1]
OutFile = sys.argv[2]
RefFile = sys.argv[3]
Num = int(sys.argv[4])
Identity = 95 if len(sys.argv) <6 else int(sys.argv[5])
min_len  = 100 if len(sys.argv) <7  else int(sys.argv[6])
max_len  = 1450 if len(sys.argv) <8 else int(sys.argv[7])

Ref = open(RefFile,"r")
dat = Ref.readlines()
Ref.close()
dat_list = [x.strip().split("\t")  for x in dat]
group_dict = dict(zip([x[1].split("|")[1] for x in dat_list],[x[0] for x in dat_list]))

colname = ["qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"]
df = pd.read_table(InFile,header=None,names=colname)
df.iloc[:,2:12] = df.iloc[:,2:12].astype(int)

df_con = df[(df["pident"] >= Identity) &(df["length"] >min_len)&(df["length"] < max_len)]
df_con["qsp"] = df_con["qseqid"].apply(lambda x:group_dict.get(x.split("|")[1],"Unknown"))

all_sp = df_con["qsp"].tolist()
df_con_new = df_con

df_con_new = df_con_new.reset_index().drop(["index"],axis=1)
# region define
#df_con["q-reg"] = df_con.apply(lambda x:Construct_Region(x,"qid","qstart","qend"),axis=1)
df_con_new["s-reg"] = df_con_new.apply(lambda x:Construct_Region(x,"sseqid","sstart","send"),axis=1)
df_con_re = (df_con_new.loc[:,"s-reg"]).to_frame(name="Uniq_reg")
df_con_re.drop_duplicates(inplace=True)  ## duplicate reg,left uniq region

## reconstruct region datafram
df_con_re["ID"] = df_con_re["Uniq_reg"].apply(lambda x:x.split(":")[0])
df_con_re["start"] = df_con_re["Uniq_reg"].apply(lambda x:x.split(":")[1].split("-")[0])
df_con_re["end"] = df_con_re["Uniq_reg"].apply(lambda x:x.split(":")[1].split("-")[1])
df_con_re["length"] = df_con_re.apply(lambda x:int(x["end"]) - int(x["start"]),axis=1)
df_con_re = df_con_re.reset_index().drop(["index"],axis=1)
df_con_re.iloc[:,2:] = df_con_re.iloc[:,2:].astype(int)

ref_dic = {}
for i in df_con_new.index:
    ref_key = df_con_new.loc[i,"s-reg"]
    if ref_key in ref_dic.keys():
        ref_dic[ref_key] = ref_dic[ref_key]+","+df_con_new.loc[i,"qsp"]
    else:
        ref_dic[ref_key] = df_con_new.loc[i,"qsp"]


new_ref= {}
for k,v in ref_dic.items():
    ref_dic[k] = v.split(",")
    new_ref[k] = set(v.split(","))

df_stat = pd.DataFrame(ref_dic.keys(),columns=["Ref_Gene"])
All_sp_len = len(set(df_con_new["qsp"]))
for i in df_stat.index:
    ref = df_stat.loc[i,"Ref_Gene"]
    df_stat.loc[i,"Exist Num"] = len(ref_dic[ref])
    df_stat.loc[i,"Uni Num"] = len(new_ref[ref])
df_stat["dup rate"] = df_stat.apply(lambda x:round(int(x["Exist Num"])/All_sp_len,3),axis=1)
df_stat["perc"] = df_stat.apply(lambda x:round(int(x["Uni Num"])/All_sp_len,3),axis=1)
df_stat = df_stat[df_stat["perc"] > 0.1]

SelectNum = Num if len(df_stat) >Num else len(df_stat)
df_final = df_stat.sort_values(["perc"],ascending=False)[:SelectNum]

df_final = pd.merge(df_final,df_con_re,left_on="Ref_Gene",right_on="Uniq_reg").drop(["Uniq_reg"],axis=1)
df_final.to_csv(OutFile,sep="\t",index=False)