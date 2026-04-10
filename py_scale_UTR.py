import sys,numpy

# take CDS as 1-2
# 5UTR will be X-1
# 3UTR will be 2+Y
# See https://github.com/olarerin/metaPlotR for more details about scale factor

if len(sys.argv) != 5:
	print ("Usage: python py_scale_UTR.py 5UTR_start(0.8) 3UTR_end(2.5) input(MetaGene/Liver_IP_rep1.dist.highIso.txt) output(MetaGene/Liver_IP_rep1.dist.highIso.sc.txt)")
	sys.exit()

s5 = float(sys.argv[1])
e3 = float(sys.argv[2])
in_fl = open(sys.argv[3])
out_fl = open(sys.argv[4],"w")

def scale_number(unscaled, to_min, to_max, from_min, from_max):
	return (to_max-to_min)*(unscaled-from_min)/(from_max-from_min)+to_min

head = 1
for L in in_fl:
	if head:
		out_fl.write(L.strip() + "\n")
		head -= 1
		continue
	chrom,coord,gene_name,refseqID,rel_location,utr5_st,utr5_end,cds_st,cds_end,utr3_st,utr3_end,utr5_size,cds_size,utr3_size = L.strip().split("\t")	
	rel_location = float(rel_location)

	if rel_location < 1:
		rel_location_n = str(scale_number(rel_location, s5, 1, 0, 1))
	elif rel_location >= 1 and rel_location <= 2:
		rel_location_n = str(rel_location)
	elif rel_location >2 and rel_location <=3:
		rel_location_n = str(scale_number(rel_location, 2, e3, 2, 3))
	else:
		print ("Error:\n" + L.strip())
		sys.exit()
	rel_location = rel_location_n
	out_fl.write("\t".join([chrom,coord,gene_name,refseqID,rel_location,utr5_st,utr5_end,cds_st,cds_end,utr3_st,utr3_end,utr5_size,cds_size,utr3_size]) + "\n")
in_fl.close()
out_fl.close()
