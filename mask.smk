import os 
import sys
import re
import re
import pysam 
import pandas as pd
from datetime import date


# delete if not debug
DEBUG=True
def tempd(fname):
	if(DEBUG):
		return(fname)
	return(temp(fname))
SDIR=os.path.dirname(workflow.snakefile)
#print (SDIR)
#sys.exit()
CWD=os.getcwd()
shell.prefix(f"source {SDIR}/env.cfg ; set -eo pipefail; ")
shell.prefix(f"source ~/.bashrc_qsub")
shell.prefix("source activate annotation")

FASTA = os.path.abspath( config["fasta"] )
FAI = FASTA + ".fai"
assert os.path.exists(FAI), f"Index must exist. Try: samtools faidx {FASTA}"

# WILDCARDS
NIDS = min(100, len(open(FAI).readlines()) )
IDS = [ "{:08}".format(ID+1) for ID in range(NIDS) ]

SM = "asm"
if("sample" in config): SM = config["sample"]
SPECIES = "human"
if("species" in config): SPECIES = config["species"]
THREADS = 8
if("threads" in config): THREADS = config["threads"]
#RM_LIB = config["rm_lib"]

SPECIES = '"'+SPECIES+'"'
print (SPECIES)
#sys.exit()
SMS = [SM]

wildcard_constraints:
	SM="|".join(SMS),
	ID="\d+",

# temp output
FASTA_FMT = f"Masked/temp/{SM}_{{ID}}.fasta"

#
# OUTPUTS
#
RM = os.path.abspath(f"Masked/{SM}_repeatmasker.out")
RMD_lib = os.path.abspath(f"{SM}-families.fa")
RMD = os.path.abspath(f"Masked/{SM}_repeatmodeler.out")
RMBED = os.path.abspath(f"Masked/{SM}_repeatmasker.out.bed")
RMDBED = os.path.abspath(f"Masked/{SM}_repeatmodeler.out.bed")
TRFBED = os.path.abspath(f"Masked/{SM}_trf.bed")
ULTRABED  = os.path.abspath(f"Masked/{SM}_ultra.bed")

rule all:
	input:
		RM,
		RMBED,
		#RMD,
		#RMDBED,
		#TRFBED,
		#ULTRABED,

rule spliit_fasta:
	input:
		fasta = FASTA,
	output:
		fastas = tempd(expand(FASTA_FMT, ID=IDS)),
	threads: 1
	resources:
		mem=8,
		hrs=7,
	run:
		fasta = pysam.FastaFile(input["fasta"])
		outs = [open(f,"w+") for f in output.fastas]
		outidx = 0
		for name in fasta.references:
			seq = fasta.fetch(name)
			outs[outidx].write( ">{}\n{}\n".format(name, seq) )
			outidx += 1
			if(outidx == NIDS): outidx = 0 

		for out in outs: 
			out.close()


####################################################################
#################### REPEAT MASKER #################################
####################################################################

rule RunRepeatMasker:
	input:
		fasta = FASTA_FMT,
	output:
		out = tempd(FASTA_FMT + ".out"),
		msk = tempd(FASTA_FMT + ".masked"),
	resources:
		mem=8,
		hrs=72,
	#conda:
	#	"envs_conda.yaml"
	threads: THREADS/2 
	shell:"""
echo "RM on {input.fasta}"
RepeatMasker \
	-s \
	-xsmall \
	-e ncbi \
	-species {SPECIES} \
	-dir $(dirname {input.fasta}) \
	-pa {threads} \
	{input.fasta} 
	
if [ -f "{output.msk}" ]; then
    echo "masked fasta exists"
else 
    echo "No repeats found, copying unmasked fasta to masked fasta"
	cp {input.fasta} {output.msk}
fi
"""

#
# RepeatMasker merge output
#
rule mergeRM:
	input:
		outs = expand( rules.RunRepeatMasker.output.out, ID=IDS, SM=SM),
	output:
		out=RM,
	resources:
		mem=4,
		hrs=2,
	shell:"""
head -n 3 {input.outs[0]} > {output.out}  && \
	tail -q -n +4 {input.outs} >> {output.out}
"""

rule mergeRMbed:
	input:
		out = rules.mergeRM.output.out,
	output:
		bed=RMBED,
	resources:
		mem=4,
		hrs=2,
	threads: 1
	shell:"""
python {SDIR}/script/out2bed.py {input.out} {output.bed}
"""

rule RepeatMasker:
    input:
        out = rules.mergeRM.output.out,
        bed = rules.mergeRMbed.output.bed,


####################################################################
################# REPEAT MODELER ###################################
####################################################################



rule RunRepeatModeler:
	input:
		fasta = FASTA,
	output:
		lib = RMD_lib,
	resources:
		mem=4,
		hrs=72,
	#conda:
	#	"envs_conda.yaml"
	threads: THREADS*3 
	shell:"""
	source activate annotation
	BuildDatabase -name {SM} -engine ncbi {input.fasta}
	RepeatModeler -database {SM} -engine ncbi -pa {threads}
	"""

rule RunRepeatMasker2:
	input:
		fasta = FASTA_FMT + ".masked",
		rm_db = rules.RunRepeatModeler.output.lib,
	output:
		out = tempd(FASTA_FMT + ".masked.out"),
		msk = tempd(FASTA_FMT + ".masked.masked"),
	resources:
		mem=8,
		hrs=72,
	#conda:
	#	"envs_conda.yaml"
	threads: THREADS/2 
	shell:"""
echo "RM on {input.fasta}"
RepeatMasker \
	-s \
	-xsmall \
	-e ncbi \
	-lib {input.rm_db} \
	-dir $(dirname {input.fasta}) \
	-pa {threads} \
	{input.fasta} 
	
if [ -f "{output.msk}" ]; then
    echo "masked fasta exists"
else 
    echo "No repeats found, copying unmasked fasta to masked fasta"
	cp {input.fasta} {output.msk}
fi
"""

	
#
# RepeatModeler merge output
#
rule mergeRM2:
	input:
		outs = expand( rules.RunRepeatMasker2.output.out, ID=IDS, SM=SM),
	output:
		out=RMD,
	resources:
		mem=4,
		hrs=2,
	shell:"""
head -n 3 {input.outs[0]} > {output.out}  && \
	tail -q -n +4 {input.outs} >> {output.out}
"""

rule mergeRM2bed:
	input:
		outs = rules.mergeRM2.output.out,
	output:
		bed=RMDBED,
	resources:
		mem=8,
		hrs=2,
	threads: 1
	shell:"""
python {SDIR}/script/out2bed.py {input.outs} {output.bed}
"""


rule RepeatMasker2:
    input:
        out = rules.mergeRM2.output.out,
        bed = rules.mergeRM2bed.output.bed,


####################################################################
####################### TRF MASKER #################################
####################################################################


rule run_ultra:
	input:
		fasta = FASTA_FMT,
	output:
		ultra = tempd(FASTA_FMT + ".ultra"),
		ultra_out = tempd(FASTA_FMT + ".ultra.bed")
	benchmark:
		FASTA_FMT + ".bench2"
	resources:
		mem=16,
		hrs=120,
	#conda:
	#	"envs_conda.yaml"
	threads: 8
	shell:"""
/net/eichler/vol26/home/hsjeong/nobackups/programs/ULTRA/ultra -ml 2000 -mu 3 -n {threads} -p 4000 -ws 16384 -f {output.ultra} {input.fasta}
python {SDIR}/json2bed.py {output.ultra} {output.ultra_out}
"""

rule ultra_merge:
	input:
		dats = expand(rules.run_ultra.output.ultra_out, ID=IDS, SM=SM),
	output:
		bed = ULTRABED,
	resources:
		mem=12,
		hrs=24,
	threads: 4
	shell:"""
cat {input.dats} | sort -k1,1 -k2,2n > {output.bed}
"""

rule run_trf:
	input:
		fasta = FASTA_FMT,
	output:
		dat = tempd(FASTA_FMT + ".dat")
	benchmark:
		FASTA_FMT + ".bench"
	resources:
		mem=24,
		hrs=72,
	#conda:
	#	"envs_conda.yaml"
	threads: 1
	shell:"""
trf {input.fasta} 2 7 7 80 10 50 500 -l 20 -h -ngs > {output.dat}
"""

rule trf_bed:
	input:
		dats = expand(rules.run_trf.output.dat, ID=IDS, SM=SM),
	output:
		bed = TRFBED,
	resources:
		mem=8,
		hrs=72,
	threads: 1
	run:
		trf = []
		header = '#chr start end PeriodSize CopyNumber ConsensusSize PercentMatches PercentIndels Score A C G T Entropy Motif Sequence'.split()
		for datf in input.dats:
			chrom = None
			sys.stderr.write( "\r" + datf )			
			with open(datf, 'r') as dat:
				for line in dat:
					splitline = line.split()
					if( line.startswith("Sequence:") ):
						chrom = int(line.split()[1].strip())
						#sys.stderr.write(chrom + "\n")
					elif( line.startswith("@") ):
						chrom = splitline[0][1:].strip() # grab everything after the @ in the first word
					else:
						# Catch index errors when line is blank
						try:
							# Check if in header sequence (all non-header lines start with an int: start pos)
							try:
								int(splitline[0])
							except ValueError:
								continue
							trf.append([chrom] + splitline[ 0: (len(header)-1) ] )
						except IndexError:
							pass
		trf = pd.DataFrame(trf, columns=header)
		print(trf.shape)
		
		trf["start"] = trf["start"].astype(int)
		trf.sort_values(by=["#chr", "start"], inplace=True)
		print("done sorting trf")

		trf.to_csv(output.bed, sep="\t", index=False)

rule trf:
    input:
        bed = rules.trf_bed.output.bed


####################################################################
####################### MASKED FASTA ###############################
####################################################################
# TODO


