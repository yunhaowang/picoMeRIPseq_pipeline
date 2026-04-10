import sys

if len(sys.argv) != 5:
	print ("Usage: python py_extract_stringtie_expr.py input_1(StringTie/Liver_Input.gtf) input_2(StringTie/Liver_Input.gene.tab) output_1(StringTie/Liver_Input.gene.txt) output_2(StringTie/Liver_Input.iso.txt)")
	sys.exit()

in_gtf = open(sys.argv[1])
in_gene = open(sys.argv[2])
out_gene = open(sys.argv[3],"w")
out_iso = open(sys.argv[4],"w")

dic_iso = {}
dic_gene = {}
dic = {}
gene_list = []
for L in in_gtf:
	if L.startswith("#"): continue
	fte = L.strip().split("\t")[2]
	if fte != "transcript": continue

	cov = L.strip().split("cov \"")[1].split("\"")[0]
	fpkm = L.strip().split("FPKM \"")[1].split("\"")[0]
	tpm = L.strip().split("TPM \"")[1].split("\"")[0]

	gene_id = L.strip().split("gene_id \"")[1].split("\"")[0]
	iso_id = L.strip().split("transcript_id \"")[1].split("\"")[0]
	gene_name = L.strip().split("gene_name \"")[1].split("\"")[0]

	dic_iso[iso_id] = [gene_id,gene_name,cov,fpkm,tpm]
	
	if gene_id not in gene_list:
		dic_gene[gene_id] = [iso_id]
		dic[gene_id] = gene_name
		gene_list.append(gene_id)
	else:
		dic_gene[gene_id].append(iso_id)
in_gtf.close()

head = 1
dic_cov = {}
for L in in_gene:
	if head:
		head -= 1
		continue
	gene_id,gene_name,chrom,strand,s,e,cov,fpkm,tpm = L.strip().split("\t")
	dic_cov[gene_id] = cov
in_gene.close()

gene_list.sort()
for gene in gene_list:
	fpkm_list = []
	tpm_list = []
	iso_list = dic_gene[gene]
	iso_list.sort()
	iso_list_new = []

	for iso in iso_list:
		gene_id,gene_name,cov,fm,tm = dic_iso[iso]
		iso_list_new.append(":".join([iso,cov,fm,tm]))
		out_iso.write("\t".join([iso,gene_id,gene_name,cov,fm,tm]) + "\n")
		fpkm_list.append(float(fm))
		tpm_list.append(float(tm))

	fpkm = sum(fpkm_list)
	tpm = sum(tpm_list)
	out_gene.write("\t".join([gene,dic[gene],dic_cov[gene_id],str(fpkm),str(tpm),",".join(iso_list_new)]) + "\n")
out_iso.close()
out_gene.close()
