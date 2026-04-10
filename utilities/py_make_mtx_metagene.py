import sys,re

if len(sys.argv) != 3:
	print ("Usage: python py_make_mtx_metagene.py input(MetaGene/metagene.fofn) output(MetaGene/metagene.txt)")
	sys.exit()

in_fofn = open(sys.argv[1])
out_fl = open(sys.argv[2],"w")

out_fl.write("Sample\tRel_Location\n")

for L in in_fofn:
	p = L.strip().split("/")[-1].split(".")[0]
	in_fl = open(L.strip())
	head = 1
	for ln in in_fl:
		if head:
			head -= 1
			continue
		rel_loc = ln.strip().split("\t")[4]
		out_fl.write("\t".join([p,rel_loc]) + "\n")
	in_fl.close()
in_fofn.close()
out_fl.close()
