#!/usr/bin/env bash
####################
# Title: Watershed pipeline 
# Author: Sherry Yang
# Date: 2026-02-16
# Description: SLURM job script for generating WS annotation, intermediate connecting steps, WS training and evaluation, and performance and result visualization
####################

####################
# SLURM Configuration
####################
#SBATCH --job-name=watershed_pipeline
#SBATCH --output=slurm_%x_%j.out
#SBATCH --error=slurm_%x_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=750G

#logging
# if [[ -z "$WSH_LOGGING" ]]; then
#     log_time=$(date +"%Y%m%d_%H%M%S")
#     export WSH_LOGGING=1
#     exec "$0" "$@" > "driver_run_logs/try_various_cell_type_${log_time}.log" 2>&1
# fi

## script starts here
eval $(pixi shell-hook)

working_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) 
# first navigate to where this shell script is present, or where all the relevant data and annotations are stored 
base_dir=$working_dir/Watershed-SV/
vcf_path=$base_dir/ADRC.cohort_combined.GRCh38.sniffles_jointcall.nonduplicated_samples.vcf
dataset="adrc_scRNA"
collapse_instruction=$base_dir/collapse_annotation_instructions.tsv

## this is to run at different pcs for result comparison or multiple ome in a loop 
for i in 10; do 

    annotation_dir="$base_dir/generate_annotation/output_bulk"
    annotation_done="$annotation_dir/.generate_annotations_ABC.done"

    # first step is to generate annotations that will be collected in one directory, specify modes in the generate_annotations_ABC executable as well as input files including vcf, output directory, and references. 
    if [[ ! -f "$annotation_done" ]]; then
    generate_annotations_ABC -x $base_dir/input -p population -l None -v $vcf_path -f PASS -k 100000 -r 0.01 -o $working_dir/generate_annotation/output_bulk -b $working_dir/GRCh38_no_alt_analysis_set_GCA_000001405.15.genome_bound_file.genome -g $working_dir/gencode.v32.annotation.gtf.gz -c $base_dir/vep_cache -e False -i True
        touch "$annotation_done"
    else
        echo "Checkpoint: $annotation_done exists, skipping generate_annotations_ABC"
    fi

    # the next step is to combine the annotations generated in step 1 
    dir_name="$base_dir/run_watershed_combinations/output_with_e2g_annotation/adrc_scRNA_cd8t_geno1_expr${i}"
    
    mkdir -p "$dir_name"
    # Define dynamic file names
    output_filename_gene_level="$dir_name/combined_annotation_pre_merge_noimpute.medZ.gene_level.csv"
    output_filename_sv_level="$dir_name/combined_annotation_pre_merge_noimpute.medZ.csv"

    expressions_filename="$working_dir/formatted_data/adrc_scRNA_cd8t_scaled_geno1_expr${i}.tsv"
    combined_annotation="$dir_name/watershed_input_combined_annotations_adrc_scRNA_cd8t_geno1_expr${i}.csv"
    echo "$expressions_filename"
    
    # specify the name of the columns
    expression_field="expression_zscore"
    expression_id_field="subjectid"

    # Gene level annotation checkpoint
    gene_level_done="$output_filename_gene_level.done"

    #modify the modes and parameters as specified in the "combine_all_annotations_ABC_polars" script 
    if [[ ! -s "$output_filename_gene_level" || ! -f "$gene_level_done" ]]; then
        combine_all_annotations_ABC_polars --vcf $vcf_path --genotypes "$annotation_dir/intermediates/pipeline_input_genotypes.tsv" --genes "$annotation_dir/intermediates/genes.bed" --gene-sv "$annotation_dir/intermediates/gene_sv_slop.100000.bed" --annotation-dir "$annotation_dir" --outfile "$output_filename_gene_level" --expressions "$expressions_filename" --expression-field "$expression_field" --expression-id-field "$expression_id_field" --maf-mode upload --maf-file "$annotation_dir/intermediates/pipeline_maf.tsv" --length-mode extract --length-field SVLEN --CN-mode extract --collapse-mode gene --remove-control-genes --filter-rare --minimum-support-tissue-count 1 --flank 100000 --collapse-annotation-instruction $collapse_instruction 
        touch "$gene_level_done"
    else
        echo "Checkpoint: $output_filename_gene_level exists, skipping gene level annotation"
    fi

    # SV level annotation checkpoint
    sv_level_done="$output_filename_sv_level.done"

    # same process as above, with the difference that the output is based on sv level instead of collapsed to gene_level
    if [[ ! -s "$output_filename_sv_level" || ! -f "$sv_level_done" ]]; then
        combine_all_annotations_ABC_polars --vcf $vcf_path --genotypes "$annotation_dir/intermediates/pipeline_input_genotypes.tsv" --genes "$annotation_dir/intermediates/genes.bed" --gene-sv "$annotation_dir/intermediates/gene_sv_slop.100000.bed" --annotation-dir "$annotation_dir" --outfile "$output_filename_sv_level" --expressions "$expressions_filename" --expression-field "$expression_field" --expression-id-field "$expression_id_field" --maf-mode upload --maf-file "$annotation_dir/intermediates/pipeline_maf.tsv" --length-mode extract --length-field SVLEN --CN-mode extract --collapse-mode gene-sv --remove-control-genes --filter-rare --minimum-support-tissue-count 1 --flank 100000 --collapse-annotation-instruction $collapse_instruction
        touch "$sv_level_done"
    else
        echo "Checkpoint: $output_filename_sv_level exists, skipping SV level annotation"
    fi

    # Watershed prep checkpoint, this step is to merge the annotations from step 1 with omic level input and assigning N2 pairs 
    combined_annotation_done="$combined_annotation.done"
    if [[ ! -s "$combined_annotation" || ! -f "$combined_annotation_done" ]]; then
        eval_watershed_prep --gene-sv-annotation "$output_filename_sv_level" --gene-annotation "$output_filename_gene_level" --collapse-instructions $collapse_instruction --output "$combined_annotation"
        touch "$combined_annotation_done"
    else
        echo "Checkpoint: $combined_annotation exists, skipping eval_watershed_prep"
    fi
        
    echo "finished prepping, now run Watershed evaluation"

    model="Watershed_exact"  # Can take on "Watershed_exact", "Watershed_approximate", "RIVER"
    number_of_dimensions="1" # Can take on any real number greater than or equal to one, represent # of ome as input
    input_file="$combined_annotation" # Input file
    output_prefix_evaluate="${dir_name}/evaluate_model_${model}_number_of_dimensions_${number_of_dimensions}_on_rna_pc_${i}"
    output_prefix_prediction="${dir_name}/prediction_model_${model}_number_of_dimensions_${number_of_dimensions}_on_rna_pc_${i}"
    dirichlet_prior=10.0
    l2_prior=NA
    pvalue_threshold=0.01
    # default fraction of outliers is 0.05
    pvalue_fraction=0.05

    # Watershed evaluation checkpoint
    evaluate_done="${output_prefix_evaluate}_watershed_evaluate.done"
    if [[ ! -f "${output_prefix_evaluate}_watershed_evaluate.done" ]]; then
        Rscript evaluate_watershed.R --input $input_file --number_dimensions $number_of_dimensions --output_prefix $output_prefix_evaluate --model_name $model --dirichlet_prior_parameter $dirichlet_prior --l2_prior_parameter $l2_prior --binary_pvalue_threshold $pvalue_threshold --n2_pair_pvalue_fraction $pvalue_fraction
        touch "$evaluate_done"
    else
        echo "Checkpoint: $evaluate_done exists, skipping evaluate_watershed.R"
    fi
    
    echo "done evaluting, now prediction" 

    # Prediction checkpoint
    prediction_done="${output_prefix_prediction}_watershed_predict.done"
    if [[ ! -f "$prediction_done" ]]; then
        Rscript predict_watershed.R --training_input $input_file --prediction_input $input_file --number_dimensions $number_of_dimensions --output_prefix $output_prefix_prediction --model_name $model
        touch "$prediction_done"
    else
        echo "Checkpoint: $prediction_done exists, skipping predict_watershed.R"
    fi

    echo "prediction done"

    # Visualization checkpoint
    annotation_weight_png="$base_dir/watershed_evaluate_and_visualization/${dataset}_rna_pc_${i}_annotation_weights_model_${model}_dimensions_${number_of_dimensions}.png"
    annotation_weight_done="${annotation_weight_png}.done"
    if [[ ! -f "$annotation_weight_done" ]]; then
        Rscript $working_dir/watershed_model_weights_rna.R --model_object_rds ${output_prefix_prediction}_prediction_object.rds --combined_annotation_input $input_file --output_file_png "$annotation_weight_png"
        touch "$annotation_weight_done"
    else
        echo "Checkpoint: $visualization_done exists, skipping visualization"
    fi

    model_performance_comparison_png="$base_dir/watershed_evaluate_and_visualization/${dataset}_rna_pc_${i}_annotation_weights_model_${model}_dimensions_${number_of_dimensions}.png"
    model_performance_comparison_done="${model_performance_comparison}.done"
    if [[ ! -f "$model_performance_comparison_done" ]]; then
        Rscript $working_dir/watershed_model_weights_GAM_WS_performance_comparison.R --eval_object_rds ${output_prefix_evaluate}_evalute_object.rds  --output_file_png "$model_performance_comparison_png" 
    
done