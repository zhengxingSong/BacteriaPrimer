import pandas as pd
import re
import sys


def USAGE():
    print("%s InFile(psl) TaxData PrimerGenus OutDir"%sys.argv[0])

if len(sys.argv) <5:
    USAGE()
    sys.exit(1)

def PuteRate(inf,pat):
    pat_list = inf.split(",")
    need_list = [x for x in pat_list if re.findall(pat,x)]
    rate = round(len(need_list)/len(pat_list),2)
    return(rate)

def ComRate(inf,pat,Spdict):
    pat_list = inf.split(",")
    need_list = [x for x in pat_list if re.findall(pat,x)]
    SpciesNum = int(Spdict[pat])
    rate = round(len(need_list)/SpciesNum,2)
    return(rate)

InFile = sys.argv[1]
DataTax = sys.argv[2]
PrimerGenus = sys.argv[3]
OutDir = sys.argv[4]
Pute = 0.3 if len(sys.argv) <6 else float(sys.argv[5])
Com = 0.1 if len(sys.argv) < 7 else float(sys.argv[6])
Num = 200

df_tax = pd.read_table(DataTax)
id2sp = dict(zip(df_tax["BACT"],df_tax["Species"]))
GenuSpecies = df_tax[["Genus","Species"]].drop_duplicates()
Sp2Ge = dict(zip(GenuSpecies["Species"],GenuSpecies["Genus"]))

Species_counts = pd.DataFrame(df_tax.drop_duplicates(subset=["Species"])["Genus"].value_counts()).reset_index()
Strain_count = pd.DataFrame(df_tax["Species"].value_counts()).reset_index()
Species_counts.columns = ["Genus","SpCount"]
Strain_count.columns = ["Species","StCount"]

count = pd.merge(GenuSpecies,Species_counts)
count = pd.merge(count,Strain_count)
SpNum = dict(zip(count["Genus"],count["SpCount"]))


df_psl = pd.read_table(InFile,header=None)
df_primer = df_psl.loc[:,[9,13]]
df_primer["BACT"] = df_primer[9].apply(lambda x:x.split("|")[0])
df_primer["ref"] = df_primer["BACT"].apply(lambda x:id2sp.get(x,"None"))
# df_primer["status"] = df_primer["ref"].apply(lambda x: 1 if re.findall(PrimerGenus.replace("_"," "),x) else 0)
# df_primer["ori"] = df_primer[13].apply(lambda x:" ".join(x.split("_")[0:2]))
# df_primer1 = df_primer[df_primer["status"] == 1]
# df_not_uni = df_primer1[df_primer1["ref"] != df_primer1["ori"]]

ref_dic = {}
for i in df_primer[13].drop_duplicates().tolist():
    df_tmp = df_primer[df_primer[13] ==i]
    ori = ",".join(set(df_tmp["ref"].tolist()))
    ref_dic[i] = ori
df_filter = pd.DataFrame(ref_dic.items(),columns=["Primer","Ref_Sp"])

pat = PrimerGenus.split("_")
df_filter["Target"] = df_filter["Ref_Sp"].apply(lambda x:len([y for y in x.split(",") if re.findall("\[*%s\]* %s"%(pat[0],pat[1]),y)]))
df_filter["UnTarget"] = df_filter["Ref_Sp"].apply(lambda x:len([y for y in x.split(",") if not re.findall("\[*%s\]* %s"%(pat[0],pat[1]),y)]))

df_not_uni = df_filter[df_filter["UnTarget"] != 0]
df_not_uni["Pute"] = df_not_uni["Ref_Sp"].apply(lambda x:PuteRate(x,Sp2Ge[PrimerGenus.replace("_"," ")]))
df_not_uni["Complete"] = df_not_uni["Ref_Sp"].apply(lambda x:ComRate(x,Sp2Ge[PrimerGenus.replace("_"," ")],SpNum))

df_not_uni.to_csv("%s/Total_%s_Not_Uni_Gene.txt"%(OutDir,PrimerGenus),sep="\t",index=False)
df_Genus_primer = df_not_uni[(df_not_uni["Pute"] > Pute)&(df_not_uni["Complete"] > Com)]
SelectNum = Num if len(df_Genus_primer) > 200 else len(df_Genus_primer)
df_Genus_primer = df_Genus_primer.sort_values(["Pute","Complete"],ascending=False)[:SelectNum]
df_Genus_primer.to_csv("%s/Total_%s_Genus_primer.txt"%(OutDir,PrimerGenus),sep="\t",index=False)

OutFile=open("%s/Genus_Primer_ID.txt"%OutDir,"w+")
OutFile.writelines("\n".join(df_Genus_primer["Primer"].tolist()))
OutFile.close()