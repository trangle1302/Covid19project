Spatialprotein covid screen
===============================

### Steps to preprocess plates for annotation-tool:
Step 1. Convert file names and max-projection
`python generate_covid_dataset-max-projection.py`

Step 2. Segmentation with HPACellSegmentator and convert mask to geojson file
`python generate_covid_dataset-max-projection_segmentation.py`
`python mask_to_json.py`

### Steps to analyse intensity hits:
Step 3. Violin plot
`python s1_covid_quantify_well.py`
`python s2_combine_quantification_df.py`
`python s3_define_cellsize_and_virus_threshold.py`
`python s4_violinplot_adjustment.py`

Step 4. Fold change, volcano and circle plot
`python covid19plots.py`