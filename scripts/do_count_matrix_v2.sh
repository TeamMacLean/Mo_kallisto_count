#!/bin/bash


function usage {
  cat <<EOM
usage: $(basename "$0") [OPTION]...
  -h        Display help
  dryrun    Do a dry run of the pipeline and get a report on what steps need to be done.
  unlock    Remove the lockfile preventing unauthorised run after failure of process.
  dag       Generate a graphviz dot file of the process DAG

Requires a file 'config.yaml' with 4 entries, e.g

 scratch: "/some/place/in/your/scratch"
 projroot: "/some/place/in/your/home/count"
 reference_genome: "/tsl/data/sequences/fungi/magnaporthe/ensemble_genome/Magnaporthe_oryzae.MG8.50.cdna.all.fa"
EOM
  exit 2
}

if [ -z "$1" ]
then
    sbatch -J count \
    -o count.log \
    --wrap="source snakemake-5.5.3; snakemake -s scripts/counts_v2.snakefile counts --cluster 'sbatch --partition={params.queue} -c {threads} --mem={params.mem}' -j 20 --latency-wait 60"
elif [ $1 = 'unlock' ]
then
    sbatch -J unlock \
        -o count.log \
        --wrap="source snakemake-5.5.3; snakemake -s scripts/counts_v2.snakefile --unlock" 
elif [ $1 = "dryrun" ]
then
    sbatch -J count \
    -o count.log \
    --wrap="source snakemake-5.5.3; snakemake -s scripts/counts_v2.snakefile -n"
elif [ $1 = "dag" ]
then
    sbatch -J dag \
    -o count.log \
    --wrap="source snakemake-5.5.3; snakemake --dag -s scripts/counts_v2.snakefile  > count.dot" \
    --partition="tsl-short" \
    --mem="16G"
elif [ $1 = '-h' ]
then
  usage
fi
