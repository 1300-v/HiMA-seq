import os
SAMPLES = ["2024seq2988_1010_22_E1"]
WORKDIR = "/s1/project/SpatialOmics/"
os.chdir(WORKDIR)

rule all:
	input:
		expand("data/rawdata_cellranger/{sample}/{sample}_S1_L001_R1_001.fastq.gz", sample=SAMPLES),
		expand("data/rawdata_cellranger/{sample}/{sample}_S1_L001_R2_001.fastq.gz", sample=SAMPLES)

rule split_barcode:
	input:
		fq_R1 = "data/rawdata/RNA/{sample}_R1.fastq.gz",
		fq_R2 = "data/rawdata/RNA/{sample}_R2.fastq.gz"
	output:
		fq_R1_un = temp("data/rawdata_cellranger/{sample}/{sample}_S1_L001_R1_001.fastq"),
		fq_R1_gz = "data/rawdata_cellranger/{sample}/{sample}_S1_L001_R1_001.fastq.gz",
		fq_R2_gz = "data/rawdata_cellranger/{sample}/{sample}_S1_L001_R2_001.fastq.gz",
	params:
		split_script = "src/python/Split_BC.py"
	shell:"""
		python {params.split_script} -i {input.fq_R2} -o {output.fq_R1_un}
		pigz -p 100 -k -f {output.fq_R1_un}
		cp {input.fq_R1} {output.fq_R2_gz}
	"""