import numpy as np
import pandas as pd
import re
import random
import os
import argparse

def QKtrans(tax_list,data,col_name):
    QK_dict = {}
    true_exist = []
    for k in tax_list:
        k_list = k.strip().split(" ")
        pat = "\[*%s\]*"%k_list[0] if len(k_list)==1 else "\[*%s\]* %s"%(k_list[0],"*".join(k_list[1:]))
        #print(pat)
        trans_list = [x for x in set(data[col_name].dropna()) if re.findall(pat,x)]
        true_exist.extend(trans_list)
        if trans_list != 0:
            for v in trans_list:
                QK_dict[v] = k
    data = data[data[col_name].isin(true_exist)]
    data["key"] = data[col_name].apply(lambda x:QK_dict[x])
    return(data)

def SpeStat(DList,Trans,HMP,MGYG,RefSeq,GeneB):
    MGYG_Key = Trans(DList,MGYG,"Species")
    HMP_Key = Trans(DList,HMP,"Species")
    assemble_genebank_Key = Trans(DList,GeneB,"organism_name")
    assemble_ref_Key = Trans(DList,RefSeq,"organism_name")
    ## stats data cal
    assemble_ref_need = pd.DataFrame(assemble_ref_Key["key"].value_counts())
    assemble_genebak_need = pd.DataFrame(assemble_genebank_Key["key"].value_counts())
    MGYG_need = pd.DataFrame(MGYG_Key["key"].value_counts())
    HMP_need = pd.DataFrame(HMP_Key["key"].value_counts())

    df_combine = [assemble_ref_need,assemble_genebak_need,HMP_need,MGYG_need]
    merge_df = pd.concat(df_combine,axis=1,sort=False,join="outer").fillna(0)
    merge_df = merge_df.reset_index()
    merge_df.columns = ["Species","Refseq Num","Genebank Num","HMP Num","MGYG Num"]
    return(merge_df)

def Sample(RequireL,Num):
    OutL = []
    for i in range(Num):
        tmp = random.choice(RequireL)
        RequireL.remove(tmp)
        OutL.append(tmp)
    return(OutL)

def DownPath(DownL,refseq,OutPath,SeqNum):
    for i in DownL:
        nam = i.replace(" ","_")
        tmp = refseq[refseq["key"] == i]
        DPath_df = tmp[tmp["refseq_category"] == "representative genome"]
        if len(tmp) > max(SeqNum):
            ProteinL = tmp[tmp["refseq_category"] != "representative genome"].index.tolist()
            IndexL = Sample(ProteinL,max(SeqNum))
            DPath_df = DPath_df.append(tmp.loc[IndexL,])
        else:
            DPath_df = tmp
        DPath_df = DPath_df[~DPath_df["ftp_path"].isin([x for x in tmp["ftp_path"] if x.endswith("fasta")])]
        if len(DPath_df) != 0:
            if not os.path.exists(os.path.join(OutPath,nam)):
                os.mkdir(os.path.join(OutPath,nam))
            refer_status = [x for x in DPath_df["refseq_category"].tolist() if re.findall("representative genome",x)]
            if not refer_status:
                DPath_df["seq_rel_date"] = pd.to_datetime(DPath_df["seq_rel_date"])
                DPath_df["seq_rel_date"] = DPath_df["seq_rel_date"].apply(lambda x:100*x.year+x.month)
                refseq_id = DPath_df.loc[DPath_df["seq_rel_date"] == DPath_df["seq_rel_date"].max()].index.tolist()
                if len(refseq_id) != 1:
                    refseq_id = random.choice(refseq_id)
                    DPath_df.loc[refseq_id,"refseq_category"] = "representative genome"
                else:
                    DPath_df.loc[refseq_id,"refseq_category"]="representative genome"
            else:
                pass
            DPath_df["cds"] = DPath_df["ftp_path"].apply(lambda x:x.split("/")[-1])
            DPath_df[["refseq_category","cds","assembly_accession","ftp_path","organism_name"]].to_csv("%s/%s.txt"%(os.path.join(OutPath,nam),nam),index=False,sep="\t")

            

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        default=None,
        type=str,
        required=True,
        help="The input data. Should contain the .tsv files (or other data files) for the task.",
    )
    parser.add_argument(
        "--ref",
        default=None,
        type=str,
        required=True,
        help="The reference data of HMP/MGYG taxanomy",
    )
    parser.add_argument(
        "--assemble",
        default=None,
        type=str,
        required=True,
        help="The assemble down list.",
    )
    parser.add_argument(
        "--out",
        default=None,
        type=str,
        required=True,
        help="The output dir ,should include the stats and each species ftp_path file.",
    )
    parser.add_argument(
        "--s",
        dest="spe",
        action="store_true",
        help="Species level genome will download.",
    )
    parser.add_argument(
        "--g",
        dest="gen",
        action="store_true",
        help="Genus level genome will download.",
    )
    args = parser.parse_args()

    ## constant params
    assemble_columns = ['assembly_accession','bioproject','biosample','wgs_master','refseq_category','taxid','species_taxid','organism_name','infraspecific_name','isolate','version_status','assembly_level','release_type','genome_rep','seq_rel_date','asm_name','submitter','gbrs_paired_asm','paired_asm_comp','ftp_path','excluded_from_refseq','relation_to_type_material','asm_not_live_date']
    ref_columns = ["BACT_ID","Phylum","Class","Order","Family","Genus","Species","Strain"]
    SeqNum = [1,1000]

    ## load database
    ref = np.load(args.ref,allow_pickle=True)
    assemble= np.load(args.assemble,allow_pickle=True)
    MGYG = pd.DataFrame(ref["MGYG"],columns=ref_columns[:-1])
    HMP=pd.DataFrame(ref["HMP"],columns=ref_columns)
    assemble_ref = pd.DataFrame(assemble["ASR"],columns=assemble_columns)
    assemble_genebank = pd.DataFrame(assemble["ASG"],columns=assemble_columns)

    species = args.spe
    genus = args.gen

    down = open(args.input,"r")
    down_data = down.readlines()
    down.close()
    down_load = [x.strip() for x in down_data]
    print(len(down_load))
    if species:
        assemble_ref_Key = QKtrans(down_load,assemble_ref,"organism_name")
        MergeDF = SpeStat(down_load,QKtrans,HMP,MGYG,assemble_ref,assemble_genebank)
        MergeDF = MergeDF[MergeDF["Refseq Num"] > min(SeqNum)]
        MergeDF.insert(0,"Genus",pd.Series([i]*len(MergeDF)))
        MergeDF.to_csv(os.path.join(args.out,"Species_stat.txt"),index=False,sep="\t")
        DownL = MergeDF["Species"].tolist()
        ## save all species download list
        DownPath(DownL,assemble_ref_Key,args.out,SeqNum)
    
    if genus:
        hmp_Ge2Sp = pd.DataFrame(HMP.groupby(["Genus"])["Species"].agg(lambda x:x.drop_duplicates().str.cat(sep=",")).reset_index())
        Ge2Sp = dict(zip(hmp_Ge2Sp["Genus"],hmp_Ge2Sp["Species"]))

        HMP_Key = QKtrans(down_load,HMP,"Genus")
        print(len(HMP_Key["Genus"].drop_duplicates()))
        for i in HMP_Key["Genus"].unique().tolist():
            GenusPath = os.path.join(args.out,i.replace(" ","_"))
            if not os.path.exists(GenusPath):
                os.mkdir(GenusPath)
            else:
                pass
            DList = Ge2Sp.get(i,"Unaccessable").split(",")
            MergeDF = SpeStat(DList,QKtrans,HMP,MGYG,assemble_ref,assemble_genebank)
            MergeDF= MergeDF[MergeDF["Refseq Num"] > min(SeqNum)]
            MergeDF.insert(0,"Genus",pd.Series([i]*len(MergeDF)))
            MergeDF.to_csv(os.path.join(args.out,"%s_Stat.txt"%i.replace(" ","_")),sep="\t",index=False)
            DownL = MergeDF["Species"].tolist()
            assemble_ref_Key = QKtrans(DownL,assemble_ref,"organism_name")
            DownPath(DownL,assemble_ref_Key,GenusPath,SeqNum)

if __name__ == "__main__":
    try:
        main()
    except:
        print()
