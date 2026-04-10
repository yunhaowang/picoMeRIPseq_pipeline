import sys

if len(sys.argv) != 3:
	print ("Usage: python py_assign_strand2peak.py input(PeakAnno/Liver_IP_rep1.txt) output(PeakAnno/Liver_IP_rep1.bed)")
	sys.exit()

in_fl = open(sys.argv[1])
out_fl = open(sys.argv[2],"w")

dic = {}
peak_lst = []
for L in in_fl:	
	if L.strip().endswith("."): continue
	chrom,s,e,peak,qv,strand = L.strip().split(" ")[-1].split("\t")
	chr_s_e = "\t".join([chrom,s,e,peak,qv])
	dic[chr_s_e] = []
	peak_lst.append(chr_s_e)
in_fl.close()

in_fl = open(sys.argv[1])
for L in in_fl:
	if L.strip().endswith("."): continue
	chrom,s,e,peak,qv,strand = L.strip().split(" ")[-1].split("\t")
	chr_s_e = "\t".join([chrom,s,e,peak,qv])
	dic[chr_s_e].append(strand)
in_fl.close()

for peak in peak_lst:
	if len(dic[peak]) == 1:
		out_fl.write(peak + "\t" + dic[peak][0] + "\n")
out_fl.close()
