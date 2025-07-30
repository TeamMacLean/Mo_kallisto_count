# Mo_kallisto_count - Perform Kallisto counts on HPC for Magnaporthe oryzae

This repo contains a snakemake workflow that will run the Magnaporthe oryzae RNA-sequencing reads through `kallisto count` generating count files for each sample.

## Why snakemake?

`Snakemake` is a workflow manager. Crucially, it can restart where it left off, so if one step of the pipeline fails, not all previous steps need to be re-run. It manages the steps that need to be run for you. This makes it very useful for bioinformatics pipelines especially this one, where we want to divide some steps into many similar small ones. Managing a `snakemake` process is easier than managing the hundreds (literally) of files that come out from the splitting and classifying and counting steps in a pipeline like this, especially when one goes wrong.


## The pipeline:
    
    1. count extracted reads per Mo transcript per sample (kallisto)
    2. process count files into count matrix


### `kallisto` step

The M. oryzae reads are put through the `kallisto count` pipeline to generate estimates of transcript abundance in each read file. A single output file is generated for each sample. A count matrix is compiled.


## Running the pipeline

The pipeline is just a collection of scripts - see the folder `scripts`. Do not run these directly, they are run from a master script! To install the pipeline, pull this repo into a folder in your home area and `cd` into that folder. This will be your working directory.

### Source Files

#### `config.yaml`

`snakemake` pipelines are very tightly linked to the filesystem, so you need to tell it where everything it needs is. The `config.yaml` file does this.

```
 scratch: "/some/place/in/your/scratch"
 projroot: "/some/place/in/your/home/split_and_count"
 reference_genome: "/tsl/data/sequences/fungi/magnaporthe/ensemble_genome/Magnaporthe_oryzae.MG8.50.cdna.all.fa"
```

 * `scratch` is the place in your scratch that you want the intermediate files to go to
 * `projroot` is the place in your home area that you are working from, it should be the place you pulled the repo to
* `reference_genome` is the file of reference sequences that `kallisto` will use, this is currently correct for Mo genes, change it for Os genes (you will need a cultivar specific file).
 
 ### `lib/samples_to_files.csv`

 The pipeline needs to know how samples map to the input files. The file `lib/samples_to_files.csv` does this. The structure of the file is inherited from an earlier project. It looks like this

| name      | timepoint | biorep | samplename    | alt_samplename | misc_info | fq1          | fq2          |
|-----------|-----------|--------|---------------|----------------|-----------|--------------|--------------|
| co39_leaf | 0         | 1      | co39_leaf_0_1 | B1             | NA        | /path/to/fq1 | /path/to/fq2 |

You must modify this for your info, but the column positions and names must stay the same.

The key columns are

 * `name` which must contain the name
 * `samplename` which must contain a unique sample identifier - usually the concatenation of `name` `timepoint` and `biorep`.
 * `fq1` and `fq2` which contain the full path to the fastq files in the hpc.

Note these _must_ be in the 1st, 4th, 7th and 8th columns of a .csv file with a header, otherwise the data structure cannot be made.


### The submission script

Submitting the snakemake pipeline to the cluster is done from the script `bash scripts/do_count_matrix.sh`.

    1. See the help - `bash scripts/do_count_matrix_v2.sh -h`
    2. Do a dryrun (check what remains of the pipeline to run and check everything is in place) - `bash scripts/do_count_matrix_v2.sh dryrun`
    3. Check the log to see that all is ok - `less count.log`
    4. Submit the pipeline - `bash scripts/do_count_matrix_v2.sh`
    
## Checking the pipeline status once submitted

The pipeline creates a master job called `count` and that creates subjobs, this will persist for as long as the pipeline runs. Output from the main process goes
into a file `count.log`. Other job output goes into job id specific `slurm**.out` files.  


## Output folder

In this pipeline all output goes into the home working directory in a folder `results`.  Expect the following files

 * "results/all_samples_tpm_matrix.txt" - a tpm matrix
 * "results/run_metadata.txt" - the metadata file needed for `sleuth`
 * "results/kallisto_abundances.gz"  - a compressed file of abundance files for `sleuth`
 
## Restarting the pipeline

If the pipeline fails (usually because some file wasn't created somewhere), once a fix is made, the pipeline can be restarted from the point it failed. Just redo the submission step `bash scripts/do_count_matrix_v2.sh` . The pipeline will go from as far as it got previously, so if you need to redo earlier steps, you need to delete their output files. If the process crashes a lock will sometimes be generated on the folder, remove it with `bash scripts/do_count_matrix_v2.sh unlock`
