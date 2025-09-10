include: "data_prepare_split.smk"
REFERENCE = "/home/project/software/CellRanger/refdata-gex-GRCm39-2024-A/"
RESULTDIR  = "/s1/project/SpatialOmics/result_split/"

if os.path.exists(RESULTDIR):
	pass
else:
	os.mkdir(RESULTDIR)

rule CellRangerprocess:
	input:
		expand("result_split/{sample}_{readlist}/{sample}_{readlist}.mri.tgz", sample=SAMPLES, readlist=READLIST)

rule cellRanger_mapping:
	input:
	output:
		mri = "result_split/{sample}_{readlist}/{sample}_{readlist}.mri.tgz"
	params:
		fastqpath="data/rawdata_split/{sample}_{readlist}",
		outdir="result_split/{sample}_{readlist}/",
		id_sample="{sample}_{readlist}"
	shell:"""
		cellranger count --id={params.id_sample} \
		--output-dir {params.outdir} \
		--transcriptome={REFERENCE} \
		--fastqs={params.fastqpath} \
		--sample={params.id_sample} \
		--create-bam=true --localcores=50 --localmem=64
	"""