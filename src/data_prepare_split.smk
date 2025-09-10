import os
SAMPLES = ["43_E"]
WORKDIR = "/s1/project/SpatialOmics/"
os.chdir(WORKDIR)
Spotnum = 1886
READLIST = [10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000]

rule all:
	input:
		expand("data/rawdata_split/{sample}_{readlist}_R1.fastq.gz", sample=SAMPLES, readlist=READLIST),
		expand("data/rawdata_split/{sample}_{readlist}_R2.fastq.gz", sample=SAMPLES, readlist=READLIST),
		expand("data/rawdata_split/{sample}_{readlist}/{sample}_{readlist}_S1_L001_R1_001.fastq.gz", sample=SAMPLES, readlist=READLIST),
		expand("data/rawdata_split/{sample}_{readlist}/{sample}_{readlist}_S1_L001_R2_001.fastq.gz", sample=SAMPLES, readlist=READLIST)

rule fetch_R2id_form_R1:
	input:
		fq_R1 = "data/{sample}_R1.fastq.gz"
	output:
		new_fq_R1 = temp("data/rawdata_split/{sample}_{readlist}_R1.fastq"),
		new_fq_R1_gz = "data/rawdata_split/{sample}_{readlist}_R1.fastq.gz",
		R2_id = "data/rawdata_split/{sample}_{readlist}_R2id.list"
	params:
		num = lambda wildcards: Spotnum * int(wildcards.readlist),
		python_script = "src/python/fetch_n_read_from_R1.py"
	shell:"""
		python {params.python_script} -i {input.fq_R1} -n {params.num} -o1 {output.R2_id} -o2 {output.new_fq_R1}
		pigz -p 100 -k -f {output.new_fq_R1}
	"""

rule fetch_R2_rawfq_byR2id:
	input:
		fq_R2 = "data/{sample}_R2.fastq.gz",
		R2_name = rules.fetch_R2id_form_R1.output.R2_id
	output:
		new_fq_R2 = temp("data/rawdata_split/{sample}_{readlist}_R2.fastq"),
		new_fq_R2_gz = "data/rawdata_split/{sample}_{readlist}_R2.fastq.gz"
	shell:"""
		seqtk subseq {input.fq_R2} {input.R2_name} > {output.new_fq_R2}
		pigz -p 100 -k -f {output.new_fq_R2}
	"""

rule split_barcode:
	input:
		new_fq_R1 = "data/rawdata_split/{sample}_{readlist}_R1.fastq.gz",
		new_fq_R2 = "data/rawdata_split/{sample}_{readlist}_R2.fastq.gz"
	output:
		fq_R1_un = temp("data/rawdata_split/{sample}_{readlist}/{sample}_{readlist}_S1_L001_R1_001.fastq"),
		fq_R1_gz = "data/rawdata_split/{sample}_{readlist}/{sample}_{readlist}_S1_L001_R1_001.fastq.gz",
		fq_R2_gz = "data/rawdata_split/{sample}_{readlist}/{sample}_{readlist}_S1_L001_R2_001.fastq.gz",
	params:
		split_script = "src/python/Split_BC.py"
	shell:"""
		python {params.split_script} -i {input.new_fq_R2} -o {output.fq_R1_un}
		pigz -p 100 -k -f {output.fq_R1_un}
		cp {input.new_fq_R1} {output.fq_R2_gz}
	"""
