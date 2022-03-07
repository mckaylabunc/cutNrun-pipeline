import pandas as pd
import preProcessSampleConfig as pre

configfile: 'config.json'

file_info_path = config['sampleInfo']
basename_columns = config['baseNameColumns']
pool_basename_columns = config['poolBaseNameColumns']

REFGENOME = config['refGenome']
SPIKEGENOME = config['spikeGenome']
REFGENOMEPATH = config['genome'][REFGENOME]['bowtie']
controlDNAPath  = config['genome'][REFGENOME]['controlDNAPath']
chromSize_Path  = config['genome'][REFGENOME]['chrSize']

genomeSize = config['genome'][REFGENOME]['genomeSize']
readLen = config['readLen']

modules = config['module']
#########
# Validation 

if os.path.exists(file_info_path) == False:
	sys.exit('Error: {name} does not exist. Be sure to set `sampleInfo` in config.json.'.format(name = file_info_path))

if type(REFGENOME) is not str:
	sys.exit('Error: refGenome must be a string. Currently set to: {}. Double check `refGenome` in config.json.'.format(REFGENOME))

if type(SPIKEGENOME) is str:
	# Convert spikegenome to list if not already
	# allows users to pass single genomes as e.g. "spikeGenome" = "sacCer3" in config.json
	# hopefully cutting down on errors for forgetting []
	# Also ensures type safety for SPIKEGENOME
	SPIKEGENOME = [SPIKEGENOME]

if type(SPIKEGENOME) is not list:
	sys.exit('Error: spikeGenome entry in config.json is invalid. Entry must be a string or an array of strings. Currently set to: {}'.format(SPIKEGENOME))

#TODO: check that spikegenomes are inside config['genome'] & has all necessary files
#TODO: check that refgenome is inside config['genome'] & has all necessary files

#########
# Generating sampleSheet outputs
speciesList  = [REFGENOME] + SPIKEGENOME
#speciesList.append(SPIKEGENOME)

combinedGenome = '-'.join(speciesList)
#TODO: remove
#indexDict    = {REFGENOME: REFGENOMEPATH, SPIKEGENOME: SPIKEGENOMEPATH}

# TODO: fix _spikeNorm -> _{species}-spikeNorm
fragTypes    = ['allFrags', '20to120', '150to700']
normTypeList = ['', '_rpgcNorm'] + ["_{}-spikeNorm".format(species) for species in SPIKEGENOME]

sampleInfo, sampleSheet = pre.makeSampleSheets(file_info_path, basename_columns, "-", fileDelimiter = config['sampleInfoDelimiter'])
poolSampleSheet = sampleSheet.copy()

sampleSheet['fastq_trim_r1'] = expand("results/Fastq/{sample}_R{num}_trim.fastq.gz", sample = sampleSheet.baseName, num = ['1'])
sampleSheet['fastq_trim_r2'] = expand("results/Fastq/{sample}_R{num}_trim.fastq.gz", sample = sampleSheet.baseName, num = ['2'])
sampleSheet['bam']           = expand("results/Bam/{sample}_{species}_trim_q5_dupsRemoved.{ftype}", sample = sampleSheet.baseName, species = REFGENOME, ftype = {"bam"})

for frag, norm in zip(fragTypes, normTypeList):
	# Add column per peak call list
	peak_colName = 'peak_{frag}'.format(frag = frag)
	sampleSheet[peak_colName] = expand("results/Peaks/{sample}_{species}_trim_q5_dupsRemoved_{fragType}_peaks.narrowPeak", sample = sampleSheet.baseName, species = REFGENOME, fragType = frag)
	
	bed_colName = 'bed_{frag}'.format(frag = frag)
	sampleSheet[bed_colName] = expand('results/Bed/{sample}_{species}_trim_q5_dupsRemoved_{fragType}.bed', sample = sampleSheet.baseName, species = REFGENOME, fragType = frag)
	

# Add columns for all combinations of fragType & normType
## Create column names:
fragTypeCombn = fragTypes * len(normTypeList)
fragTypeCombn.sort()

normTypeCombn = normTypeList * len(fragTypes)

fragNormCombn = [''.join(frag + norm) for frag, norm in zip(fragTypeCombn, normTypeCombn)]

## Add columns to sampleSheet
for fragNorm in fragNormCombn:

	# Add column per bigwig
	bw_colName = 'bigwig_{fragNorm}'.format(fragNorm = fragNorm.replace("-", "_"))
	sampleSheet[bw_colName] = expand("results/BigWig/{sample}_{species}_trim_q5_dupsRemoved_{fragNorm}.bw", sample = sampleSheet.baseName, species = REFGENOME, fragNorm = fragNorm)

	# Add column per zNorm bigwig
	bw_colName = 'zNorm_bigwig_{fragNorm}'.format(fragNorm = fragNorm.replace("-", "_"))
	sampleSheet[bw_colName] = expand("results/BigWig/{sample}_{species}_trim_q5_dupsRemoved_{fragNorm}_zNorm.bw", sample = sampleSheet.baseName, species = REFGENOME, fragNorm = fragNorm)

	# Threshold peakcalls:
	thresh_colName = 'threshold_peaks_{fragNorm}'.format(fragNorm = fragNorm.replace("-", "_"))
	sampleSheet[thresh_colName] = expand('results/Threshold_PeakCalls/{sample}_{species}_trim_q5_dupsRemoved_{fragNorm}_thresholdPeaks.bed', sample = sampleSheet.baseName, species = REFGENOME, fragNorm = fragNorm)

sampleSheet.to_csv('sampleSheet.tsv', sep = "\t", index = False)


####
# BEGIN PIPELINE:
####

# TODO: remove
localrules: all, collect_genome_align_stats

rule all:
	input:
		expand("results/Fastq/{sample}_R{num}_trim.fastq.gz", sample = sampleSheet.baseName, num = ['1','2']),
		expand("results/Sam/{sample}_{species}_trim.sam", sample = sampleSheet.baseName, species = combinedGenome),
		expand("results/Bam/{sample}_{species}_trim_q5_dupsRemoved.{ftype}", sample = sampleSheet.baseName, species = speciesList, ftype = {"bam", "bam.bai"}),
		expand("results/Logs/{sample}_{species}_trim_q5_dupsRemoved_genomeStats.tsv", sample = sampleSheet.baseName, species = combinedGenome),
		expand("results/BigWig/{sample}_{species}_trim_q5_dupsRemoved_{fragType}{normType}.{ftype}", sample = sampleSheet.baseName, species = REFGENOME, fragType = fragTypes, normType = normTypeList, ftype = {"bw", "bg"}),
		expand("results/Peaks/{sample}_{species}_trim_q5_dupsRemoved_{fragType}_peaks.narrowPeak", sample = sampleSheet.baseName, species = REFGENOME, fragType = fragTypes),
		expand('results/Threshold_PeakCalls/{sample}_{species}_trim_q5_dupsRemoved_{fragType}{normType}_thresholdPeaks.bed', sample = sampleSheet.baseName, species = REFGENOME, fragType = fragTypes, normType = normTypeList),
		expand('results/FastQC/{sample}_R1_fastqc.html', sample = sampleSheet.baseName),
		expand('results/FastQC/{sample}_R1_trim_fastqc.html', sample = sampleSheet.baseName),
		expand('results/FQscreen/{sample}_R1_trim_screen.txt', sample = sampleSheet.baseName),
		expand('results/FQscreen/{sample}_R1_trim_screen.html', sample = sampleSheet.baseName),
		"multiqc_report.html",
		expand('results/Plots/FragDistInPeaks/{sample}_{REFGENOME}_trim_q5_allFrags_fragDistPlot.png', sample = sampleSheet.baseName, REFGENOME = REFGENOME),
		expand('results/BigWig/{sample}_{REFGENOME}_trim_q5_dupsRemoved_{fragType}_rpgcNorm_zNorm.bw', sample = sampleSheet.baseName, REFGENOME = REFGENOME, fragType = fragTypes),
		expand("results/AlignmentStats/{sample}_{species}_trim.tsv", sample = sampleSheet.baseName, species = combinedGenome),
		expand("results/AlignmentStats/{sample}_{species}_trim_q5.tsv", sample = sampleSheet.baseName, species = combinedGenome),
		expand("results/AlignmentStats/{sample}_{species}_trim_q5_dupsRemoved.tsv", sample = sampleSheet.baseName, species = combinedGenome)


rule combine_technical_reps:
	input:
		r1 = lambda wildcards : sampleInfo[sampleInfo.baseName == wildcards.sample].fastq_r1,
		r2 = lambda wildcards : sampleInfo[sampleInfo.baseName == wildcards.sample].fastq_r2
	output:
		r1 = 'results/Fastq/{sample}_R1.fastq.gz',
		r2 = 'results/Fastq/{sample}_R2.fastq.gz'
	shell:
		"""
		cat {input.r1} > {output.r1} &&
		cat {input.r2} > {output.r2}
		"""

rule fastQC:
	input:
		rules.combine_technical_reps.output.r1
	output:
		'FastQC/{sample}_R1_fastqc.html'
	envmodules:
		modules['fastqcVer']
	shell:
		"""
		fastqc -o ./FastQC/ -f fastq {input}
		"""

rule trim_adapter:
	input:
		r1 = rules.combine_technical_reps.output.r1,
		r2 = rules.combine_technical_reps.output.r2
	output:
		r1 = "results/Fastq/{sample}_R1_trim.fastq.gz",
		r2 = "results/Fastq/{sample}_R2_trim.fastq.gz"
	log:
		adapterStats = 'Logs/{sample}_adapterStats',
		trimStats = 'Logs/{sample}_trimStats'
	envmodules:
		modules['bbmapVer']
	shell:
		"""
		bbduk.sh in1={input.r1} in2={input.r2} out1={output.r1} out2={output.r2} ktrim=r ref=adapters rcomp=t tpe=t tbo=t hdist=1 mink=11 stats={log.adapterStats} > {log.trimStats}
		"""
		#bbduk.sh in1={input.r1} in2={input.r2} out1={output.r1} out2={output.r2} ktrim=r ref=adapters rcomp=t tpe=t tbo=t hdist=1 mink=11

rule fastQC_trim:
	input:
		rules.trim_adapter.output.r1
	output:
		'results/FastQC/{sample}_R1_trim_fastqc.html'
	envmodules:
		modules['fastqcVer']
	shell:
		"""
		fastqc -o ./FastQC/ -f fastq {input}
		"""

rule fastqScreen:
	input:
		rules.trim_adapter.output.r1
	output:
		txt = 'FQscreen/{sample}_R1_trim_screen.txt',
		html = 'FQscreen/{sample}_R1_trim_screen.html'
	params:
		fqscreenPath = modules['fqscreenPath'],
		fqscreenConf = modules['fqscreenConf']
	threads: 4
	shell:
		"""
		{params.fqscreenPath} --threads {threads} --force --aligner bowtie2 -conf {params.fqscreenConf} {input} --outdir ./FQscreen/
		"""

def get_genome_fastas(config, speciesList):
	"""
	Takes config info & list of species as input, returns dict where keys = species name, values = genome fasta paths
	"""
	get_fasta = lambda species : (species, config['genome'][species]['fasta'])
	fastas = {genome:fasta for (genome, fasta) in [get_fasta(s) for s in speciesList]}
	return(fastas)

genome_fastas = get_genome_fastas(config, speciesList)

rule bowtie2index:
    	input:
	    	genome_fastas.values(),
	output:
	    	"results/Bowtie2Index/" + combinedGenome +  ".fa",
	    	expand("results/Bowtie2Index/{genome}.{num}.bt2", genome = combinedGenome, num = ["1", "2", "3", "4", "rev.1", "rev.2"])
	params:
		module = modules['bowtie2Ver'],
		prefix = "Bowtie2Index/" + combinedGenome,
	    	combined_fa = "Bowtie2Index/" + combinedGenome +  ".fa"
	run:
	    init = True
	    for genome, fasta in genome_fastas.items():
		    if init:
			    shell("sed -e 's/>/>{genome}_/' {fasta} > {fa}".format(genome = genome, fasta = fasta, fa = params.combined_fa))
			    init = False
		    else:
			    shell("sed -e 's/>/>{genome}_/' {fasta} >> {fa}".format(genome = genome, fasta = fasta, fa = params.combined_fa))
	    shell("module purge && module load {module} && bowtie2-build {fasta} {prefix}".format(module = params.module, fasta = params.combined_fa, prefix = params.prefix))



rule align:
	input:
	    	expand("results/Bowtie2Index/{genome}.{num}.bt2", genome = combinedGenome, num = ["1", "2", "3", "4", "rev.1", "rev.2"]),
		r1 = rules.trim_adapter.output.r1,
		r2 = rules.trim_adapter.output.r2
	output:
		sam = "results/Sam/{sample}_{species}_trim.sam",
		logInfo = "Logs/{sample}_{species}_bowtie2.txt"
	threads: 8
	params:
		#refgenome = lambda wildcards: indexDict[wildcards.species],
		refgenome = "Bowtie2Index/" + combinedGenome
	envmodules:
		modules['bowtie2Ver']
	resources:
	    time="1:30:00"
	shell:
		"""
		bowtie2 --seed 123 -p {threads} -q --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 -x {params.refgenome} -1 {input.r1} -2 {input.r2} -S {output.sam} 2> {output.logInfo}
		"""
# TODO: merge w/ align rule
rule convertToBam:
	input:
		rules.align.output.sam
	output:
		'results/Bam/{sample}_{species}_trim.bam'
	envmodules:
		modules['samtoolsVer']
	resources:
	    time="1:00:00"
	shell:
		"""
		samtools view -b {input} > {output}
		"""

rule qFilter:
	input:
		'results/Bam/{sample}_' + combinedGenome + '_trim.bam'
	output:
		'results/Bam/{sample}_' + combinedGenome + '_trim_q5.bam'
	envmodules:
		modules['samtoolsVer']
	shell:
		"""
		samtools view -@ 4 -bq 5 {input} > {output}
		"""

rule markDups:
	input:
		'results/Bam/{sample}_{species}_trim_q5.bam'
	output:
		sorted = 'Bam/{sample}_{species}_trim_q5_sorted.bam',
		markedDups = 'Bam/{sample}_{species}_trim_q5_dupsMarked.bam',
		PCRdups = "PCRdups/{sample}_{species}_trim_PCR_duplicates"
	params:
		picardPath = modules['picardPath']
	envmodules:
		modules['picardVer']
	shell:
		"""
		java -Xmx8g -jar {params.picardPath} SortSam INPUT= {input} OUTPUT= {output.sorted} SORT_ORDER=coordinate &&
		java -Xmx8g -jar {params.picardPath} MarkDuplicates INPUT= {output.sorted} OUTPUT= {output.markedDups} METRICS_FILE= {output.PCRdups} REMOVE_DUPLICATES= false ASSUME_SORTED= true
		"""


# TODO: remove dups in markDups rule
rule removeDups:
	input:
		'results/Bam/{sample}_' + combinedGenome + '_trim_q5_dupsMarked.bam'
	output:
		bam = 'results/Bam/{sample}_' +  combinedGenome + '_trim_q5_dupsRemoved.bam',
		index = 'results/Bam/{sample}_' + combinedGenome + '_trim_q5_dupsRemoved.bam.bai'
	envmodules:
		modules['samtoolsVer']
	shell:
		"""
		samtools view -@ 4 -bF 0x400 {input} > {output.bam} &&
		samtools index {output.bam} {output.index}
		"""

rule collect_genome_align_stats:
	input:
		bam = 'results/Bam/{sample}_' +  combinedGenome + '_trim_q5_dupsRemoved.bam',
		index = 'results/Bam/{sample}_' + combinedGenome + '_trim_q5_dupsRemoved.bam.bai'
	output:
		stats = 'Logs/{sample}_' + combinedGenome + '_trim_q5_dupsRemoved_genomeStats.tsv'
	params:
		mlr = modules['mlrPath']
	envmodules:
		modules['samtoolsVer']
	shell:
		"""
		samtools idxstats {input.bam} | \
		 sed '$d' | \
		 sed 's/_/,/;s/\t/,/g' | \
		 {params.mlr} --csv --ofs tab --implicit-csv-header cat then group-by 1 then stats1 -g 1 -a sum -f 4 then rename 1,genome,4_sum,n then filter '$genome != "*"' | \
		 sed -e 's|$|\t{input.bam}|' > {output.stats}
		"""
		# sed '$d' removes last line of idxstats, which is the * line. Since we use --nounal in bowtie2, no unmapped reads will propagate, so we ignore this line.
		# NOTE: I use | instead of / in the last sed command since the input file uses a / in the name which '' expands into a sed command, causing an error
		# it's typical to use @ in this context instead of /, but since
		# I'm matching $, $@ is expanded by '' into ARGV, so can't use
		# that either. I guess I could use "" since snakemake will eval
		# first and no bash variable expansion will occurr, but I'll keep it this way in case this gets moved to a script somewhere.



rule splitSpecies:
    	input:
	    	bam = 'results/Bam/{sample}_' + combinedGenome + '_trim_q5_dupsRemoved.bam',
	    	index = 'results/Bam/{sample}_' + combinedGenome + '_trim_q5_dupsRemoved.bam.bai',
	output:
		bam = expand('results/Bam/{{sample}}_{species}_trim_q5_dupsRemoved.bam', species = speciesList),
		index = expand('results/Bam/{{sample}}_{species}_trim_q5_dupsRemoved.bam.bai', species = speciesList)
	params:
		#prefix = lambda wildcards : ["{}_".format(wildcards.species)]
		prefix = ["{}_".format(species) for species in speciesList],
		module = modules['samtoolsVer']
	run:
	    for prefix, output_bam in zip(params.prefix, output.bam):
		    #print("prefix: {}, bam: {}".format(prefix, output_bam))
		    #shell("module purge && module load {} && ".format(params.module))
		    shell("module purge && module load {} && sh scripts/get_bam_reads_prefix.sh {} {} {}".format(params.module, input.bam, prefix, output_bam))
		    shell("module purge && module load {} && samtools index {bam} {index}".format(params.module, bam = output_bam, index = output_bam + ".bai"))
		#TODO: IS SORT ORDER OF prefix and output preserved??? can I zip or do I need to do something more complex?
		# - to test this, I included a print statement above
		# - CONFIRMED: order is preserved based on speciesList order
		#TODO: consider merging sortBam rule with this on

rule sortBam:
	input:
		'results/Bam/{sample}_' + REFGENOME + '_trim_q5_dupsRemoved.bam'
	output:
		bam = 'results/Bam/{sample}_' + REFGENOME + '_trim_q5_dupsRemoved_sorted.bam',
		idx = 'results/Bam/{sample}_' + REFGENOME + '_trim_q5_dupsRemoved_sorted.bam.bai'
	envmodules:
		modules['samtoolsVer']
	threads: 4
	shell:
		"""
		samtools sort -@ {threads} -o {output.bam} {input}
		samtools index {output.bam}	
		"""

rule nameSortBam:
	input:
		'results/Bam/{sample}_{species}_trim_q5_dupsRemoved_sorted.bam'
	output:
		'results/Bam/{sample}_{species}_trim_q5_dupsRemoved_nameSorted.bam'
	envmodules:
		modules['samtoolsVer']
	shell:
		"""
		samtools sort -n {input} -o {output}
		"""

rule convertBamToBed:
	input:
		bam = 'results/Bam/{sample}_{species}_trim_q5_dupsRemoved_nameSorted.bam',
	output:
		'results/Bed/{sample}_{species}_trim_q5_dupsRemoved.bed'
	envmodules:
		modules['bedtoolsVer']
	shell:
		"""
		bedtools bamtobed -bedpe -i {input.bam} | sort -k 1,1 -k 2,2n > {output}
		"""

rule splitFragments:
	input:
	    	rules.convertBamToBed.output
	output:
		allFrags = 'results/Bed/{sample}_{REFGENOME}_trim_q5_dupsRemoved_allFrags.bed',
		smallFrags = 'results/Bed/{sample}_{REFGENOME}_trim_q5_dupsRemoved_20to120.bed',
		bigFrags = 'results/Bed/{sample}_{REFGENOME}_trim_q5_dupsRemoved_150to700.bed'
	group:
		"fragment"
	resources:
	    mem_mb=16000
	shell:
		"""
		cut -f 1,2,6,7 {input} | awk -F '\t' '{{print $0, ($3-$2)}}' - > {output.allFrags}
		awk -v OFS='\t' '($5>20) && ($5<120) {{print $0}}' {output.allFrags} > {output.smallFrags}
		awk -v OFS='\t' '($5>150) && ($5<700) {{print $0}}' {output.allFrags} > {output.bigFrags}
		"""

# TODO: get readLen from sampleInfo (reference column like sampleSheet.readLen[sampleSheet.sample == wildcards.sample])
rule makeFragmentBedGraphs:
	input:
		ref   = lambda wildcards : 'results/Bed/' + wildcards.sample + '_' + REFGENOME + '_trim_q5_dupsRemoved_' + wildcards.fragType + '.bed'
	output:
		unNorm    = 'results/BigWig/{sample}_{REFGENOME}_trim_q5_dupsRemoved_{fragType}.bg',
		rpgcNorm  = 'results/BigWig/{sample}_{REFGENOME}_trim_q5_dupsRemoved_{fragType}_rpgcNorm.bg'
	params:
		genomeSize = genomeSize,
		chromSize_Path = chromSize_Path,
		readLen = readLen
	group:
		"fragment"
	envmodules:
		modules['bedtoolsVer']
	resources:
	    mem_mb=16000
	shell:
		"""
		# Count reads in spike-in & inputs for normalization
		readCount=$(wc -l {input.ref} | sed -e 's/^  *//' -e 's/  */,/g' | cut -d , -f 1)
		rpgcScale=$(echo "scale=5; {params.genomeSize}/(${{readCount}} * {params.readLen})" | bc)

		bedtools genomecov -i {input.ref} -bga -g {params.chromSize_Path} > {output.unNorm}
		bedtools genomecov -i {input.ref} -bga -g {params.chromSize_Path} -scale ${{rpgcScale}} > {output.rpgcNorm}
		"""

rule makeSpikeNormFragmentBedGraphs:
	input:
		ref   = lambda wildcards : 'results/Bed/' + wildcards.sample + '_' + REFGENOME + '_trim_q5_dupsRemoved_' + wildcards.fragType + '.bed',
		#spike = lambda wildcards : 'Bam/' + wildcards.sample + '_' + SPIKEGENOME + '_trim_q5_dupsRemoved.bam'
		spike = 'results/Bam/{sample}_{spikeGenome}_trim_q5_dupsRemoved.bam'
	output:
		spikeNorm = 'results/BigWig/{sample}_{REFGENOME}_trim_q5_dupsRemoved_{fragType}_{spikeGenome}-spikeNorm.bg'
	params:
		genomeSize = genomeSize,
		chromSize_Path = chromSize_Path
	group:
		"fragment"
	envmodules:
		modules['bedtoolsVer']
	resources:
	    mem_mb=16000
	shell:
		"""
		# Count reads in spike-in & inputs for normalization
		spikeCount=$(samtools view -c {input.spike})
		spikeScale=$(echo "scale=5; 10000/${{spikeCount}}/" | bc)

		bedtools genomecov -i {input.ref} -bga -g {params.chromSize_Path} -scale ${{spikeScale}} > {output.spikeNorm}
		"""

rule convertToBigWig:
	input:
		'results/BigWig/{sample}_{REFGENOME}_trim_q5_dupsRemoved_{fragType}{normType}.bg'
	output:
		'results/BigWig/{sample}_{REFGENOME}_trim_q5_dupsRemoved_{fragType}{normType}.bw'
	params:
		chromSize_Path = chromSize_Path
	group:
		"bigwig"
	envmodules:
		modules['ucscVer']
	shell:
		"""
		wigToBigWig {input} {params.chromSize_Path} {output}
		"""

rule zNormBigWig:
	input:
		'results/BigWig/{sample}_{REFGENOME}_trim_q5_dupsRemoved_{fragType}_rpgcNorm.bw'
	output:
		zNorm = 'results/BigWig/{sample}_{REFGENOME}_trim_q5_dupsRemoved_{fragType}_rpgcNorm_zNorm.bw',
		zStats = 'Logs/{sample}_{REFGENOME}_trim_q5_dupsRemoved_{fragType}.zNorm'
	group:
		"bigwig"
	envmodules:
		modules['rVer']
	shell:
		"""
		Rscript --vanilla scripts/zNorm.r {input} {output.zNorm} > {output.zStats}
		"""

rule callThresholdPeaks:
	input:
		'results/BigWig/{sample}_{REFGENOME}_trim_q5_dupsRemoved_{fragType}{normType}.bw'
	output:
		'results/Threshold_PeakCalls/{sample}_{REFGENOME}_trim_q5_dupsRemoved_{fragType}{normType}_thresholdPeaks.bed'
	group:
		"bigwig"
	envmodules:
		modules['rVer']
	shell:
		"""
		Rscript --vanilla scripts/callThresholdPeaks.R {input} {output}
		"""
	
rule callPeaks:
	input:
		'results/Bed/{sample}_{REFGENOME}_trim_q5_dupsRemoved_{fragType}.bed'
	output:
		'results/Peaks/{sample}_{REFGENOME}_trim_q5_dupsRemoved_{fragType}_peaks.narrowPeak'
	params:
		control = controlDNAPath,
		prefix = 'Peaks/{sample}_{REFGENOME}_trim_q5_dupsRemoved_{fragType}'
	envmodules:
		modules['macsVer']
	shell:
		"""
		macs2 callpeak -f BEDPE -c {params.control} -n {params.prefix} -g 121400000 -t {input}  --nomodel --seed 123
		"""

rule qcReport:
	input:
		expand("results/Bam/{sample}_{species}_trim_q5_dupsRemoved.{ftype}", sample = sampleSheet.baseName, species = speciesList, ftype = ['bam', 'bam.bai']),
		expand('results/FQscreen/{sample}_R1_trim_screen.txt', sample = sampleSheet.baseName),
		expand('results/FastQC/{sample}_R1_trim_fastqc.html', sample = sampleSheet.baseName),
		expand("results/AlignmentStats/{sample}_{species}_trim.tsv", sample = sampleSheet.baseName, species = combinedGenome),
		expand("results/AlignmentStats/{sample}_{species}_trim_q5.tsv", sample = sampleSheet.baseName, species = combinedGenome),
		expand("results/AlignmentStats/{sample}_{species}_trim_q5_dupsRemoved.tsv", sample = sampleSheet.baseName, species = combinedGenome)
	output:
		"multiqc_report.html"
	envmodules:
	    modules['multiqcVer']
	shell:
		"""
		multiqc . -f -x *.out -x *.err
		"""

rule makeFragmentSizePlots_inPeaks:
	input:
		bed = 'results/Bed/{sample}_{REFGENOME}_trim_q5_dupsRemoved_allFrags.bed',
		peaks = 'results/Peaks/{sample}_{REFGENOME}_trim_q5_dupsRemoved_allFrags_peaks.narrowPeak'
	output:
		'results/Plots/FragDistInPeaks/{sample}_{REFGENOME}_trim_q5_allFrags_fragDistPlot.png'
	envmodules:
		modules['rVer']
	shell:
		"""
		Rscript --vanilla scripts/makeFragsizePlot.R {input.bed} {input.peaks} {output}
		"""

rule alignmentStats:
    	input:
	    	trim = "results/Bam/{sample}_" + combinedGenome + "_trim.bam",
		trim_q5 = "results/Bam/{sample}_" + combinedGenome + "_trim_q5.bam",
		q5_dupsRemoved = "results/Bam/{sample}_" + combinedGenome + "_trim_q5_dupsRemoved.bam"
    	output:
	    	trim = "results/AlignmentStats/{sample}_" + combinedGenome + "_trim.tsv",
		trim_q5 = "results/AlignmentStats/{sample}_" + combinedGenome + "_trim_q5.tsv",
		q5_dupsRemoved = "results/AlignmentStats/{sample}_" + combinedGenome + "_trim_q5_dupsRemoved.tsv"
    	envmodules:
	    	modules['samtoolsVer']
    	shell:
    		"""
		samtools flagstat {input.trim} > {output.trim} &&
		samtools flagstat {input.trim_q5} > {output.trim_q5} &&
		samtools flagstat {input.q5_dupsRemoved} > {output.q5_dupsRemoved}
		"""
	# TODO: use newer samtools and use json output
	#samtools flagstat -O json {input.q5_dupsRemoved} > {output.q5_dupsRemoved}

#rule makeFragmentSizePlots:
#	input:
#		bed = 'Bed/{sample}_{REFGENOME}_trim_q5_dupsRemoved_allFrags.bed'
#	output:
#		'Plots/FragDist/{sample}_{REFGENOME}_trim_q5_allFrags_cumulativeDistPlot.png'
#	params:
#		srcDirectory = srcDirectory
#	envmodules:
#		modules['rVer']
#	shell:
#		"""
#		Rscript --vanilla scripts/plotFragsizeDist.R {input.bed}
#		"""
