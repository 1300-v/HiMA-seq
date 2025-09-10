include: "data_prepare.smk"
REFERENCE = "/home/project/software/CellRanger/refdata-gex-GRCm39-2024-A/"
RESULTDIR  = "/s1/project/SpatialOmics/result_cellranger/"

if os.path.exists(RESULTDIR):
	pass
else:
	os.mkdir(RESULTDIR)

rule barcodepro:
	input:
		expand("result_cellranger/{sample}/{sample}.mri.tgz", sample=SAMPLES)

rule cellRanger_mapping:
	input:
	output:
		mri = "result_cellranger/{sample}/{sample}.mri.tgz"
	params:
		fastqpath = "data/rawdata_cellranger/{sample}",
		outdir = "result_cellranger/{sample}/"
	shell:"""
		cellranger count --id={wildcards.sample} \
		--output-dir {params.outdir} \
		--transcriptome={REFERENCE} \
		--fastqs={params.fastqpath} \
		--sample={wildcards.sample} \
		--create-bam=true --localcores=50 --localmem=64
	"""