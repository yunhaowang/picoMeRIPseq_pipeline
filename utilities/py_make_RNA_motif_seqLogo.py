import sys
import pandas
import seqlogo

#Note: seqlogo depends on weblogo package. The new version of weblogo may cause error. The version 3.6.0 of weblogo works in our test.

if len(sys.argv) != 4:
	print ("Usage: python py_make_RNA_motif_seqLogo.py input(homer2_motif_file) output_1(new_motif.mtx) output_2(new_motif.pdf)")
	sys.exit()

in_fl = open(sys.argv[1])
out_fl = open(sys.argv[2],"w")

head = 1
for L in in_fl:
	if head:
		print ("A\tC\tG\tU", file=out_fl)
		head -= 1
		continue
	a,c,g,u = [float(i) for i in L.strip().split("\t")] # sum = 1
	u = 1.0-(a+c+g)
	print ("\t".join([str(i) for i in [a,c,g,u]]), file=out_fl)
in_fl.close()
out_fl.close()

df = pandas.read_table(sys.argv[2])

ppm = seqlogo.Ppm(df, alphabet_type = "RNA")

seqlogo.seqlogo(ppm, ic_scale = True, format = 'pdf', size = 'medium',filename=sys.argv[3])
