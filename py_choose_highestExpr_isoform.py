import sys

if len(sys.argv) != 4:
	print ("Usage: python py_choose_highestExpr_isoform.py input_1(StringTie/Liver_Input.iso.txt) input_2(MetaGene/Liver_IP_rep1.dist.txt) output(MetaGene/Liver_IP_rep1.dist.highIso.txt)")
	sys.exit()

in_exp = open(sys.argv[1])
in_fl = open(sys.argv[2])
out_fl = open(sys.argv[3],"w")

dic_iso_exp = {}
head = 1
for L in in_fl:
	if head:
		head -= 1
		continue
	chrom,coord,gene_name,refseqID,rel_location,utr5_st,utr5_end,cds_st,cds_end,utr3_st,utr3_end,utr5_size,cds_size,utr3_size = L.strip().split("\t")
	dic_iso_exp[refseqID] = 0
in_fl.close()

for L in in_exp:
	iso,gene,gene_n,cov,fpkm,tpm = L.strip().split("\t")
	tpm = float(tpm)
	dic_iso_exp[iso] = tpm
in_exp.close()

in_fl = open(sys.argv[2])
head = 1
head1 = 1
for L in in_fl:
	if head:
		out_fl.write(L.strip()+"\n")
		head -= 1
		continue
	chrom,coord,gene_name,refseqID,rel_location,utr5_st,utr5_end,cds_st,cds_end,utr3_st,utr3_end,utr5_size,cds_size,utr3_size = L.strip().split("\t")
	exp = dic_iso_exp[refseqID]
	chr_d = chrom + "\t" + coord
	if head1:
		info = L.strip()
		exp_1 = exp
		chr_d_1 = chr_d
		head1 -= 1
		continue
	if chr_d != chr_d_1:
		out_fl.write(info + "\n")
		info = L.strip()
		exp_1 = exp
		chr_d_1 = chr_d
	else:
		if exp > exp_1:
			info = L.strip()
			exp_1 = exp 
			chr_d_1 = chr_d
		else:
			continue
in_fl.close()
out_fl.write(info)
