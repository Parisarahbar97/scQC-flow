
export NXF_WORK="/rds/general/user/pr422/ephemeral/work/scQC-flow"
mkdir -p "$NXF_WORK" 
rm -r /rds/general/user/pr422/ephemeral/work/scQC-flow/*

nextflow drop Parisarahbar97/scQC-flow



# BONN epileptic samples:
nextflow run Parisarahbar97/scQC-flow -r main \
  --mapping_dirs /rds/general/user/pr422/projects/puklandmarkproject/live/Users/Parisa/EPILEP/diseased/qc/input/mapping_dirs.csv \
  --outputDir /rds/general/user/pr422/projects/puklandmarkproject/live/Users/Parisa/EPILEP/diseased/qc/output_latest \
  -w "$NXF_WORK" \
  --cellbender true --gpu true \
  -profile imperial \
  -with-report -with-trace -with-timeline

# MRC control samples:
nextflow run Parisarahbar97/scQC-flow -r main \
  --mapping_dirs /rds/general/user/pr422/projects/puklandmarkproject/live/Users/Parisa/EPILEP/healthy/qc/input/mapping_dirs.csv \
  --outputDir /rds/general/user/pr422/projects/puklandmarkproject/live/Users/Parisa/EPILEP/healthy/qc/output_latest \
  -w "$NXF_WORK" \
  --cellbender true --gpu true \
  -profile imperial \
  -with-report -with-trace -with-timeline

# UCL epileptic samples:
nextflow run Parisarahbar97/scQC-flow -r main \
  --mapping_dirs /rds/general/user/pr422/projects/puklandmarkproject/live/Users/Parisa/EPILEP/diseased/UCL/qc/input/mapping_dirs.csv \
  --outputDir   /rds/general/user/pr422/projects/puklandmarkproject/live/Users/Parisa/EPILEP/diseased/UCL/qc/output \
  -w "$NXF_WORK" \
  --cellbender true --gpu true \
  -profile imperial \
  -with-report -with-trace -with-timeline