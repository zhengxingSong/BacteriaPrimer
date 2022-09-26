import pandas as pd
import sys
import re
import numpy as np

def sam2df(Inpath):
    with open(Inpath,"r") as f:
        dat = f.readlines()

    dat = [x.strip().split("\t") for x in dat]
    dat = pd.DataFrame(dat)
    dat = dat[dat[8]!=0]
    dat = dat.loc[:,[0,2,8]]
    dat.columns = ["Primer","ref","length"]
    dat["length"] = abs(dat["length"].astype(int))
    dat["length"] = dat["length"].apply(lambda x:str(x))
    dat["Orient"] = dat["Primer"].apply(lambda x:x.split("_")[-1])
    dat["primer"] = dat["Primer"].apply(lambda x:"_".join(x.split("_")[0:-1]))
    return(dat)

def Sp2GeNum(df):
    GenuSpecies = df[["Genus","Species"]].drop_duplicates()
    Species_counts = pd.DataFrame(df.drop_duplicates(subset=["Species"])["Genus"].value_counts()).reset_index()
    Strain_count = pd.DataFrame(df["Species"].value_counts()).reset_index()
    Species_counts.columns = ["Genus","SpCount"]

    Sp2Ge = dict(zip(GenuSpecies["Species"],GenuSpecies["Genus"]))
    SpNum = dict(zip(Species_counts["Genus"],Species_counts["SpCount"]))
    return(Sp2Ge,SpNum)

def PuteRate(inf,pat):
    pat_list = inf.split(",")
    need_list = [x for x in pat_list if re.findall(pat,x)]
    rate = round(len(need_list)/len(pat_list),2)
    return(rate)

def ComRate(inf,pat,Spdict):
    pat_list = inf.split(",")
    need_list = [x for x in pat_list if re.findall(pat,x)]
    SpciesNum = int(Spdict.get(pat,10000))
    rate = round(len(need_list)/SpciesNum,2)
    return(rate)

def LenStd(num):
    if num == "0":
        std = 0
    elif "," in num:
        num_list = num.split(",")
        std = np.std([int(x) for x in num_list])
    else:
        std = 0
    return(std)

def main():
    PrimerFile =sys.argv[1]
    TaxHmp= sys.argv[2]
    TaxMgyg=sys.argv[3]
    HmpSam = sys.argv[4]
    MgygSam =sys.argv[5]
    OutDir = sys.argv[6]
    PrimerGenus=sys.argv[7]

    df_primer = pd.read_table(PrimerFile)
    df_tax_hmp = pd.read_table(TaxHmp)
    df_tax_mgyg = pd.read_table(TaxMgyg)

    HSp2Ge,HSp2Num = Sp2GeNum(df_tax_hmp)
    MSp2Ge,MSp2Num = Sp2GeNum(df_tax_mgyg)

    df_hmp = sam2df(HmpSam)
    df_hmp["BACT"] = df_hmp["ref"].apply(lambda x:x.split("|")[0])
    df_merge_hmp = pd.merge(df_hmp,df_tax_hmp,how="left")
    df_final_hmp = df_final_mgyg = df_merge_hmp.groupby(["primer"])["Orient","length","Species"].agg(lambda x:x.drop_duplicates().str.cat(sep=",")).reset_index()
    df_final_hmp = df_final_hmp[df_final_hmp["length"]!="0"]

    df_mgyg = sam2df(MgygSam)
    df_mgyg["Species_rep"] = df_mgyg["ref"].apply(lambda x:"MGYG-HGUT-%s"%x.split("_")[0][-5:])
    df_merge_mgyg = pd.merge(df_mgyg,df_tax_mgyg,how="left")
    df_final_mgyg = df_merge_mgyg.groupby(["primer"])["Orient","length","Species"].agg(lambda x:x.drop_duplicates().str.cat(sep=",")).reset_index()
    df_final_mgyg = df_final_mgyg[df_final_mgyg["length"]!="0"]

    df_F = pd.merge(df_final_hmp,df_final_mgyg,how="outer",on=("primer","Orient"))
    df_F["length_y"].fillna("0",inplace=True)
    df_F["length_x"].fillna("0",inplace=True)
    df_F["Species_y"].fillna("None",inplace=True)
    df_F["Species_x"].fillna("None",inplace=True)
    df_F.loc[df_F["Species_y"] == "","Species_y"] = PrimerGenus.replace("_"," ")
    #df_F = df_F.dropna(subset=["length_y","length_x"])
    df_F["sd_x"] = df_F["length_x"].apply(lambda x:LenStd(x))
    df_F["sd_y"] = df_F["length_y"].apply(lambda x:LenStd(x))
    df_F["sd"] = df_F.apply(lambda x:(x["sd_x"]+x["sd_y"])/2,axis=1)

    df_F["Pure_H"] = df_F["Species_x"].apply(lambda x:PuteRate(x,HSp2Ge.get(PrimerGenus.replace("_"," "),"Unknown")))
    df_F["Complete_H"] = df_F["Species_x"].apply(lambda x:ComRate(x,HSp2Ge.get(PrimerGenus.replace("_"," "),"Unknown"),HSp2Num))
    df_F["Pure_M"] = df_F["Species_y"].apply(lambda x:PuteRate(x,MSp2Ge.get(PrimerGenus.replace("_"," "),"Unknown")))
    df_F["Complete_M"] = df_F["Species_y"].apply(lambda x:ComRate(x,MSp2Ge.get(PrimerGenus.replace("_"," "),"Unknow"),MSp2Num))

    df_F = df_F.drop(columns=["sd_x","sd_y"])
    df_all = pd.merge(df_primer,df_F,how="inner")
    df_all = df_all.reset_index().drop("index",axis=1)
    df_all = df_all.drop(["gene_id","Index"],axis=1)
    df_all = df_all.drop_duplicates(["primer"])

    if len(df_all) != 0:
        df_all.to_csv("%s/%s_Genus_Primer.txt"%(OutDir,PrimerGenus),index=False,sep="\t")
    else:
        print("Sorry,This Species %s don't have primer"%PrimerGenus)

if __name__ == "__main__":
    main()
