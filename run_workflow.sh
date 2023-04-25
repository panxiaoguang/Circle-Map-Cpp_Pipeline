REF_PATH="/dellfsqd2/ST_LBI/USER/panxiaoguang/app/data_repo/GRCh38/"
CONTAINER_REF_PATH="/GRCh38/"
JOB_NUM=1
chmod a+x workflow/qsub.py

snakemake \
  -s workflow/run_circle_map.smk \
  -p \
  --use-singularity \
  --singularity-args "--bind $REF_PATH:$CONTAINER_REF_PATH" \
  --cluster "workflow/qsub.py" \
  -j $JOB_NUM
