# Circle-Map-Cpp_Pipeline
This is a snakemake pipeline used for detecting eccDNA from circle-seq data based on Circle-Map-Cpp.
# overview
![workflow](https://github.com/panxiaoguang/Circle-Map-Cpp_Pipeline/blob/main/myWorkFlow.png)
# usage
If submission , this pipeline can submit to different nodes according to your required resources. If total memory is more than 50 GB, it will use the supermem node automatically.

## Download software
Please note that all the software, including `bwa`,`samtools`,`fastqc`,`fastp` and `Circle-Map-Cpp` has been pre-built into singularity images. You don't need to be puzzled by any environment or dependence!

So please download all the softwares from this [link](https://bgitech-my.sharepoint.com/:f:/g/personal/panxiaoguang_genomics_cn/EgR6UVyGBPpBqZ4Rw8pkeYkBMno6UiRg1UTmSxDc4om6jg?e=mAbytU) before using this pipeline, and then 
put all of them into `workflow/softwares`
## Configure 

#### Prepare sample tables 

Please prepare tab-delimited files that contain all the information about your samples. The file should be like this without a header:

|  unknown |  rawData/unknown_circle_reads_1.fastq |  rawData/unknown_circle_reads_2.fastq |
| ------------ | ------------ | ------------ |
|  ... | ...  |  ... |

#### Modify the config.yml

Change `samples`,`reference` and resources according to your requirements.

It should be known that `reference` should be a file path in the singularity container using "--bind."

#### Modify the run script

You will find a script called "run_workflow.sh" in the main folder; it's a script used for running the whole pipeline. Change `REF_PATH` to your reference directory; `CONTAINER_REF_PATH` should be consistent with
`reference` in config.yml. For example, if your `reference` is `/GRCh38/hg38full.fa`, `CONTAINER_REF_PATH` must be its parent directory, which is `/GRCh38/` .

`JOB_NUM` means how many jobs will be submitted to the SGE system.
`chmod a+x workflow/qsub.py` means giving an executable ability to the` qsub.py`.

```bash
snakemake \
  -s workflow/run_circle_map.smk \
  -p \
  --use-singularity \
  --singularity-args "--bind $REF_PATH:$CONTAINER_REF_PATH" \
  --cluster "workflow/qsub.py" \
  -j $JOB_NUM
  ```
  
This is the command used to run the pipeline.

#### option: modify qsub.py

Please modify qsub command to which you need.

```bash
if resources < 50:
    os.system("qsub -cwd -l vf={resources},num_proc={threads} -P P21Z25400N0107 -q st.q {script}".format(resources = resources, threads=threads, script=jobscript))
else:
    os.system("qsub -cwd -l vf={resources},num_proc={threads} -P st_supermem -q st_supermem.q {script}".format(resources = resources, threads=threads, script=jobscript))
```

# test

### Downloading the raw data

We can download the left and right Illumina reads using the following commands:

```bash
wget https://raw.githubusercontent.com/iprada/Circle-Map/master/tutorial/unknown_circle_reads_1.fastq
wget https://raw.githubusercontent.com/iprada/Circle-Map/master/tutorial/unknown_circle_reads_2.fastq
```
### Install snakemake

```
conda create -c conda-forge -c bioconda -n snakemake snakemake
conda activate snakemake
snakemake --help
```

### Run pipeline

If you want to submit your work using qsub, you can use 

```
sh run_workflow.sh
```
otherwise, you can use

```
## run
REF_PATH="?"
CONTAINER_REF_PATH="?"
snakemake -s workflow/run_circle_map.smk --cores 1 --use-singularity --singularity-args "--bind $REF_PATH:$CONTAINER_REF_PATH"
## dry-run
snakemake -s workflow/run_circle_map.smk --cores 1 -np --use-singularity --singularity-args "--bind $REF_PATH:$CONTAINER_REF_PATH"
## run and printing shell cmd
snakemake -s workflow/run_circle_map.smk --cores 1 -p --use-singularity --singularity-args "--bind $REF_PATH:$CONTAINER_REF_PATH"
```
