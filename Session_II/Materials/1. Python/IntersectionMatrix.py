"""
@author: Tatiana Lenskaia
Updated: May 7, 2021
"""


import core_methods_mlmb as cm;
import time as tm;

# Data preprocessing

# (1) Extract phage genomes from multifasta file
print("Extracting individual genome from multifasta file...")
fInName = "172phages.fasta"
fColName = cm.MFastaCutter(fInName)
t_phs = cm.GetListFromFile(fColName)


#fRowName = "101_E.coli.txt"
#t_bac = cm.GetListFromFile(fRowName)
#fColName = "172phages.txt"



# (2) Specify bacterial genome
t_bac = ["U00096.3", "BA000007.3"]






# (3) Creating metadictionary for phage genomes and checking data quality
sep = ","

d = {}
m = 40


fOutName = "Info_"+fInName.rsplit(".",1)[0]+".csv"
fOut = open(fOutName, "w")
fOut.write("NCBI_ID"+sep+"Genome_size_bp"+sep+"GC"+sep+"Unrecognized_bases"+"\n")

print("Reading-in data and checking data quality...")
t1 = tm.time()
for it in t_phs:
    fVName = it+".fasta"
    text_v = cm.GetText(fVName)
    res = cm.CountGC(text_v)
    
    fOut.write(it+sep+str(len(text_v))+sep+str(res[0])+sep+str(res[1])+"\n")

    d_v = cm.CreateDictLoc(text_v, m)
    if it not in d:
        d[it] = d_v
    else:
        print("Duplicate id!")
        
fOut.close()

t2 = tm.time()
print("Metadictionary for", str(len(d)),"phages is created in", round(t2-t1,4),"sec.")


# (4) Compute phage fingerprints for bacterial strains
print("Computing phage fingerprints for bacterial strains...")
ct = 0
M = []
for bac in t_bac:
    t3 = tm.time()
    fBacName = bac+".fasta"
    text_b = cm.GetText(fBacName)
    d_b = cm.CreateDictLoc(text_b, m)
    tt = []
    ct = ct+1
    for v in d:
        t_ints = cm.FindIntersection(d[v], d_b)
        tt.append(len(t_ints)/len(d[v]))
    M.append(tt)
    t4 = tm.time()
    print(ct)
    print(bac, tt, t4-t3)
    


cm.PrintMatrix(m,"_2E.coli_pathogen", t_bac, t_phs, M, sep = ",")

