# Squirrel_ diet

Mapping of shotgun sequencing data against the 2025 version of the UNITE database.

# Usage

## Clone repository

```sh
git clone https://github.com/alberdilab/squirrel_diet.git
```

## Add reads and UNITE

Place all reads in: **resources/reads**

```sh
cd squirrel_diet/resources/reads
sh getreads.sh
cd ../..
```

Place the decompressed fasta file of the UNITE database (10.15156/BIO/3301232) in: **resources/database**

```sh
cd squirrel_diet/resources/database
wget https://s3.hpc.ut.ee/plutof-public/original/b02db549-5f04-43fc-afb6-02888b594d10.tgz
tar zxvf b02db549-5f04-43fc-afb6-02888b594d10.tgz
cd ../..
```

## Run the pipeline

```sh
screen -S squirrel_diet
cd squirrel_diet
snakemake --workflow-profile profile/slurm
```
