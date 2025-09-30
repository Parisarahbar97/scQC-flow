export NXF_WORK="/rds/general/user/pr422/ephemeral/work/scQC-flow"
mkdir -p "$NXF_WORK" 
rm -r /rds/general/user/pr422/ephemeral/work/scQC-flow/*

nextflow drop Parisarahbar97/scQC-flow



nextflow run Parisarahbar97/scQC-flow -r parisa/auto-expected-cells \
  --mapping_dirs /rds/general/user/pr422/projects/puklandmarkproject/live/Users/Parisa/EPILEP/diseased/qc/input/mapping_dirs.csv \
  --outputDir /rds/general/user/pr422/projects/puklandmarkproject/live/Users/Parisa/EPILEP/diseased/qc/output_cb_auto \
  -w "$NXF_WORK" \
  --cellbender true --gpu true \
  -profile imperial \
  -with-report -with-trace -with-timeline
