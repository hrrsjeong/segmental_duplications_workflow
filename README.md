# *Segmental duplications workflow
a snakemake workflow for calling / postprocessing segmental duplications from genome assembly. 

# Running pipelines
```
mkdir -p log; snakemake --ri --jobname "{rulename}.{jobid}" --drmaa " -l centos=7 -V -cwd -j y -o ./log -e ./log -l h_rt={resources.hrs}:00:00 -l mfree={resources.mem}G -pe serial {threads} -w n -S /bin/bash" -w 60 -s mask.smk -k --restart-times 1 --use-conda --use-envmodules -j 50 -p --config sample=SAMPLE fasta=$(realpath path/to/fasta.fa) species='Simiiformes'
```
