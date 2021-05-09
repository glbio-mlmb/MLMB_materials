
"""
@author: Tatiana Lenskaia

Updated: May 7, 2021
"""

import random;



def MFastaCutter(fInName):
    #Multi fasta cutter
    
    t_names = []
    text = ""
    name = ""
    
    fIn = open(fInName,"r")
    for line in fIn:
        #line = line.strip()
        if ">" in line:
            if text != "" and name != "":
                fOut = open(name+".fasta","w")
                fOut.write(text)
                fOut.close()
            text = line
            name = line[1:].split(" ")[0]
            t_names.append(name)
        else:
            text = text+line
    #print(len(t_names)) 
    
    if text != "" and name != "":
        fOut = open(name+".fasta","w")
        fOut.write(text)
        fOut.close()
    
    fname = fInName.rsplit(".",1)[0]+".txt"
    fOut = open(fname, "w")
    s = ""
    for it in t_names:
        s = s+it+"\n"
    fOut.write(s[:(-1)])
    fOut.close()
    return fname


def CountGC(seq):
    seq = seq.lower();
    n_seq = len(seq)
    
    n_a = seq.count("a");
    n_c = seq.count("c");
    n_g = seq.count("g");
    n_t = seq.count("t");
    
    n_bases = n_a + n_c + n_g + n_t
    
    n_other = n_seq - n_bases;
    
    if n_bases != 0:
        gc = round(1.0*(n_c+n_g)/n_bases*100,2)
    else:
        gc = "N/A"
    
    return [gc, n_other]



def GetListFromFile(fInName):
    t = []
    fIn = open(fInName, "r")
    for line in fIn:
        line = line.strip()
        if line != "":
            if line not in t:
                t.append(line)
    #print(len(t))
    fIn.close()
    return t


def GetDictFromFile(fInName, sep, header, unique_col = 0):
    fIn = open(fInName,"r")
    lines = fIn.readlines()
    if header == 1:
        header_line = lines[0]
        lines = lines[1:]
    
    d = {}
    t = []
    
    n = len(lines[0].split(sep))
    
    for line in lines:
        line = line.strip()
        if line != "":
            t_line = line.split(sep)
            if len(t_line) != n:
                print("Check format!", line, t_line)
                if t_line[0] not in d:
                    t.append(t_line[0])
                    d[t_line[0]] = t_line[1:]
                
            else:
                
                col_id = t_line[unique_col]
                tt = t_line[0:unique_col]+t_line[unique_col+1:]
                
                if col_id not in d:
                    t.append(col_id)
                    d[col_id] = tt
                else:
                    print("First colum has non-unique values!")
    fIn.close()
    return [t,d]








#Updated: May 4, 2019
# Needs update: base on the file extention (fasta or gb)    
def GetText(finName):
    '''Extracts text from a single fasta file'''
    fin = open(finName, 'r')
    text = ''
    for line in fin:
        line = line.strip()
        if (line != "") and (line[0] != '>'):
            text = text + line
    fin.close()
    return text


def CreateDictLocD_upper_count(text, m, gtp = "l"):
    if (m <= 0) or (m > len(text)):
        print("(n = "+str(m)+") is not a valid window size for this genome!!!");
        return {}
	
    d_g = dict()
    nn = len(text);
    gtype = gtp.lower();
    gtype = gtype[0];
	
    if gtype == "c":
        text = text + text[0:(m-1)];
        lastpos = nn;
    elif gtype == "l":
        lastpos = nn-m+1;
    else:
        print("Is this genome linear or circular?");
        return d_g;
		
    for ii in range (lastpos):
        bl = text[ii:(ii+m)]; 
        bl = bl.upper();
        if bl in d_g :
            d_g[bl] = d_g[bl] + 1
        else:
            d_g[bl] = 1
    return d_g;	



# Create a dictionary of unique strings of length mm with frequencies and locations for texttext
#Linear genome or circular genome!!!
def CreateDictLoc(text, mm, gtp = "l"):
	if (mm< 0) or (mm > len(text)):
		#print "N is bigger than genome size!!!";
		return {};
	
	
	d_g = dict()
	nn = len(text);
	gtype = gtp.lower();
	gtype = gtype[0];
	
	if gtype == "c":
		text = text + text[0:(mm-1)];
		lastpos = nn;
	elif gtype == "l":
			lastpos = nn-mm+1;
			#print lastpos;
	else:
		#print "Is this genome linear or circular?";
		return d_g;

		
	for ii in range (lastpos):
		bl = text[ii:(ii+mm)]; 
		#bl = bl.lower();
		#print ii, bl;
		if bl in d_g :
			tt = d_g[bl];
			#tt[0] = tt[0] + 1;
			tt.append(ii);
			d_g[bl] = tt;
		else:
			tt = list();
			#tt.append(1);
			tt.append(ii);
			d_g[bl] = tt;
			
	return d_g;	




def FindIntersection(d_g11, d_g22):
	
	t_ints = list()
	nn1 = len(d_g11);
	nn2 = len(d_g22);
	
	if nn1 <= nn2:
		d_first = d_g11;
		d_last = d_g22;
	else:
		d_first = d_g22;
		d_last = d_g11;
	
	for s in d_first:
		if s in d_last:
			t_ints.append(s);
	return t_ints;


def CheckIntersection(fInName1, fInName2, m = 40):
    text1 = GetText(fInName1)
    d1 =  CreateDictLocD_upper_count(text1, m, gtp = "l")
    text2 = GetText(fInName2)
    d2 =  CreateDictLocD_upper_count(text2, m, gtp = "l")
    
    t_ints = FindIntersection(d1,d2)
    
    return[len(t_ints), len(d1), len(d2)]
    
    
    
    
    
    
def PrintMatrix(k, pref, row_index, col_index, matrix, sep = ","):
    fOut = open(str(k)+pref+".csv","w");
    fOut.write(str(k)+sep+sep.join(str(x) for x in col_index)+"\n")
    for i in range(len(matrix)):
        s = str(row_index[i])
        for el in matrix[i]:
            s = s+sep+str(el)
        fOut.write(s+"\n")
        
    fOut.close()
    return

def PrintSubMatrix(k, pref, row_names, row_index, col_names, col_index, matrix, sep = ","):
    if row_index != [] and col_index != []:
        fOut = open(str(k)+"_"+pref+".csv","w");
        #fOut.write(str(k)+sep+sep.join(str(x) for x in col_index)+"\n")
        s = str(k)
        for jj in range(len(col_index)):
            s = s+sep+col_names[col_index[jj]]
        fOut.write(s+"\n")

        
        for ii in range(len(row_index)):
            i = row_index[ii]
            s = str(row_names[row_index[ii]])
            
            for jj in range(len(col_index)):
                j = col_index[jj]
                
                if i == j:
                    sign = "-"
                else:
                    sign = ""
                    
                s = s+sep+sign+str(matrix[i][j])
            
            fOut.write(s+"\n")
            
        fOut.close()
    return