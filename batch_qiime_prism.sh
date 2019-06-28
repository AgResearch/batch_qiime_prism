#!/bin/bash

declare -a files_array

function get_opts() {

   DRY_RUN=no
   DEBUG=no
   HPC_TYPE=slurm
   OUT_DIR=
   MAX_TASKS=1
   FORCE=no
   taxonomy_blast_database=
   taxonomy_lookup_file=
   seqlength_min=200
   similarity=0.97
   help_text="
\n
./batch_qiime_prism.sh  [-h] [-n] [-d] [-b blast_data -t tax_data -m min_length -s similarity  -O outdir [-C local|slurm ] input_file_names\n
\n
\n
example:\n
./batch_qiime_prism.sh -n -b /dataset/datacache/archive/metagenomics/build/combined_092016_silva.fasta -t /dataset/datacache/archive/metagenomics/build/combined_092016_silva.taxonomy -m 200 -s 0.99  -O /dataset/PJ_MGS00116/ztmp/Output_BWA01_R1_99   /dataset/PJ_MGS00116/scratch/MGS00116/Output_BWA01/*_R1_*.fastq.trimmed.gz\n
\n
"

   # defaults:
   while getopts ":nhdfO:C:b:t:m:s:" opt; do
   case $opt in
       n)
         DRY_RUN=yes
         ;;
       d)
         DEBUG=yes
         ;;
       h)
         echo -e $help_text
         exit 0
         ;;
       f)
         FORCE=yes
         ;;
       O)
         OUT_DIR=$OPTARG
         ;;
       C)
         HPC_TYPE=$OPTARG
         ;;
       b)
         taxonomy_blast_database=$OPTARG
         ;;
       t)
         taxonomy_lookup_file=$OPTARG
         ;;
       m)
         seqlength_min=$OPTARG
         ;;
       s)
         similarity=$OPTARG
         ;;
       \?)
         echo "Invalid option: -$OPTARG" >&2
         exit 1
         ;;
       :)
         echo "Option -$OPTARG requires an argument." >&2
         exit 1
         ;;
     esac
   done

   shift $((OPTIND-1))

   FILE_STRING=$@

   # this is needed because of the way we process args a "$@" - which 
   # is needed in order to parse parameter sets to be passed to the 
   # aligner (which are space-separated)
   declare -a files="(${FILE_STRING})";
   NUM_FILES=${#files[*]}
   for ((i=0;$i<$NUM_FILES;i=$i+1)) do
      files_array[$i]=${files[$i]}     
   done
}


function check_opts() {
   if [  -z "$OUT_DIR" ]; then
      echo "must specify OUT_DIR ( -O )"
      exit 1
   fi
   if [ ! -d $OUT_DIR ]; then
      echo "OUT_DIR $OUT_DIR not found"
      exit 1
   fi
   if [[ $HPC_TYPE != "local" && $HPC_TYPE != "slurm" ]]; then
      echo "HPC_TYPE must be one of local, slurm"
      exit 1
   fi
   if [ ! -f ${taxonomy_blast_database}.nin  ]; then
      echo "bad blast database (cant see ${taxonomy_blast_database}.nin ) (you might need to supply the full path ?)"
      exit 1
   fi
   if [ ! -f $taxonomy_lookup_file ]; then
      echo "tax lookup file ($taxonomy_lookup_file ) not found (you might need to supply the full path  ?) "
      exit 1
   fi
   python -c "print float('$similarity')" >/dev/null 2>&1
   if [ $? != 0 ]; then
      echo "looks like similarity requested ( $similarity ) is not a number"
      exit 1
   fi
   python -c "print float('$seqlength_min')" >/dev/null 2>&1
   if [ $? != 0 ]; then
      echo "looks like min seq length requested ( $seqlength_min ) is not a number"
      exit 1
   fi
}

function echo_opts() {
  echo OUT_DIR=$OUT_DIR
  echo DRY_RUN=$DRY_RUN
  echo DEBUG=$DEBUG
  echo HPC_TYPE=$HPC_TYPE
  echo seqlength_min=$seqlength_min
  echo taxonomy_blast_database=$taxonomy_blast_database
  echo taxonomy_lookup_file=$taxonomy_lookup_file
  echo similarity=$similarity
}

#
# edit this method to set required environment (or set up
# before running this script)
#
function configure_env() {
   export CONDA_ENVS_PATH=$CONDA_ENVS_PATH:/dataset/bioinformatics_dev/active/conda-env

   cd $BATCH_QIIME_PRISM_BIN
   cp ./batch_qiime_prism.sh $OUT_DIR
   cp ./batch_qiime_prism.mk $OUT_DIR
   cp ./add_sample_name.py $OUT_DIR
   cp ./tax_summary_heatmap.r $OUT_DIR
   cat >$OUT_DIR/tardis.toml <<EOF
seqlength_min=$seqlength_min
EOF
   echo "
export CONDA_ENVS_PATH=$CONDA_ENVS_PATH
conda activate biopython
PATH="$OUT_DIR:\$PATH"
PYTHONPATH="$OUT_DIR:\$PYTHONPATH"
" > $OUT_DIR/configure_biopython_env.src
   echo "
export CONDA_ENVS_PATH=$CONDA_ENVS_PATH
conda activate bioconductor
PATH="$OUT_DIR:\$PATH"
PYTHONPATH="$OUT_DIR:\$PYTHONPATH"
" > $OUT_DIR/configure_bioconductor_env.src
   echo "
export CONDA_ENVS_PATH=$CONDA_ENVS_PATH
conda activate qiime_1
PATH="$OUT_DIR:\$PATH"
PYTHONPATH="$OUT_DIR:\$PYTHONPATH"
export QIIME_CONFIG_FP=$OUT_DIR/qiime_config
" > $OUT_DIR/configure_qiime_env.src

   # set up a custom qiime config , so as to use a different tmpdir (default can be too small)
   cp $BATCH_QIIME_PRISM_BIN/etc/qiime_config_template $OUT_DIR/qiime_config
   echo "temp_dir $OUT_DIR" | awk '{printf("%s\t%s\n",$1,$2);}' - >> $OUT_DIR/qiime_config
   echo "run will use qiime config file $OUT_DIR/qiime_config"

   cd $OUT_DIR
}


function check_env() {
   if [ -z "$SEQ_PRISMS_BIN" ]; then
      echo "SEQ_PRISMS_BIN not set - exiting"
      exit 1
   fi
   if [ -z "$BATCH_QIIME_PRISM_BIN" ]; then
      echo "BATCH_QIIME_PRISM_BIN not set - exiting"
      exit 1
   fi
}

function get_targets() {

   rm -f $OUT_DIR/batch_qiime_targets.txt
   rm -f $OUT_DIR/processed_file_targets.txt
   rm -f $OUT_DIR/input_file_list.txt
   touch $OUT_DIR/processed_file_targets.txt

   for ((j=0;$j<$NUM_FILES;j=$j+1)) do
      file=${files_array[$j]}
      echo $file >> $OUT_DIR/input_file_list.txt
      file_base=`basename $file`
      parameters_moniker=`basename $taxonomy_blast_database`
      parameters_moniker=${parameters_moniker}.`basename taxonomy_lookup_file`
      parameters_moniker=${parameters_moniker}.len${seqlength_min}
      parameters_moniker=${parameters_moniker}.sim${similarity}

      SUMMARY_TARGETS="$SUMMARY_TARGETS $OUT_DIR/${file_base}.${parameters_moniker}.1"
      moniker=${file_base}.${parameters_moniker}
      echo $OUT_DIR/${moniker}.batch_qiime_prism >> $OUT_DIR/batch_qiime_targets.txt

      # generate wrapper
      script_filename=$OUT_DIR/${moniker}.sh

      if [ -f $script_filename ]; then
         if [ ! $FORCE == yes ]; then
            echo "found existing script $script_filename - will re-use (use -f to force rebuild ) "
            continue
         fi
      fi

      echo "#!/bin/bash
#prepare length-filtered fasta files with sequence naming syntax suitable for qiime
tardis -d $OUT_DIR -c 999999999  cat _condition_fastq2fasta_input_$file \| $OUT_DIR/add_sample_name.py $file \> _condition_uncompressedtext_output_$OUT_DIR/${file_base}.combined.fasta
" > $script_filename
      chmod +x $script_filename 


      # append output to processed targets 
      echo $OUT_DIR/${file_base}.combined.fasta >> $OUT_DIR/processed_file_targets.txt
   done 
}


function fake_prism() {
   echo "dry run ! "
   make -n -f batch_qiime_prism.mk -d -k  --no-builtin-rules -j 16 `cat $OUT_DIR/batch_qiime_targets.txt` > $OUT_DIR/batch_qiime_prism.log 2>&1
   echo "dry run : summary commands are 
   "
   exit 0
}

function run_prism() {
   # this prepares each file
   make -f batch_qiime_prism.mk -d -k  --no-builtin-rules -j 16 `cat $OUT_DIR/batch_qiime_targets.txt` > $OUT_DIR/batch_qiime_prism.log 2>&1

   # if necessary (re)create the single output file 
   target_count=`wc $OUT_DIR/processed_file_targets.txt | awk '{print $1}' -`
   if [[ ( $target_count > 0 ) || ( ! -f $OUT_DIR/combined.fa ) ]]; then
      echo "(re)building combined.fa"
      rm -f $OUT_DIR/combined.fa
      for filename in `cat $OUT_DIR/processed_file_targets.txt`; do
         cat $filename >> $OUT_DIR/combined.fa
      done
   else
      echo "(not rebuilding $OUT_DIR/combined.fa )"
   fi

   # pick OTU
   if [ -f $OUT_DIR/qiime_uclust/combined_clusters.uc ]; then
      echo "skipping otu picking as $OUT_DIR/qiime_uclust/combined_clusters.uc already exists"
   else
      tardis --shell-include-file $OUT_DIR/configure_qiime_env.src pick_otus.py -m uclust -s $similarity -i $OUT_DIR/combined.fa -o $OUT_DIR/qiime_uclust > $OUT_DIR/pick_otu.log 2>&1
      if [ $? != 0 ]; then
         echo "** error code returned from pick_otus.py job, giving up **"
         exit 1
      fi
   fi


   # pick rep_set
   if [ -f $OUT_DIR/qiime_uclust/combined_rep_set.txt  ]; then
      echo "skipping rep_set picking as $OUT_DIR/qiime_uclust/combined_rep_set.txt already exists"
   else
      tardis --shell-include-file $OUT_DIR/configure_qiime_env.src pick_rep_set.py -i $OUT_DIR/qiime_uclust/combined_otus.txt -f $OUT_DIR/combined.fa -o $OUT_DIR/qiime_uclust/combined_rep_set.txt > $OUT_DIR/pick_rep_set.log 2>&1
      if [ $? != 0 ]; then
         echo "** error code returned from pick_rep_set.py job, giving up **"
         exit 1
      fi
   fi

   # assign taxonomy
   if [ -f $OUT_DIR/qiime_uclust/combined_rep_set_tax_assignments.txt  ]; then
      echo "skipping assign_taxonomy as $OUT_DIR/qiime_uclust/combined_rep_set_tax_assignments.txt already exists"
   else
      tardis --shell-include-file $OUT_DIR/configure_qiime_env.src assign_taxonomy.py -i _condition_fasta_input_$OUT_DIR/qiime_uclust/combined_rep_set.txt -m blast -t $taxonomy_lookup_file  -b $taxonomy_blast_database -o . '_condition_uncompressedtext_product_combined_rep_set\.\d{5}_tax_assignments\.txt,'$OUT_DIR'/qiime_uclust/combined_rep_set_tax_assignments.txt' 1\>_condition_text_output_$OUT_DIR/qiime_uclust/combined_rep_set_tax_assignments.stdout 2\>_condition_text_output_$OUT_DIR/qiime_uclust/combined_rep_set_tax_assignments.stderr
      if [ $? != 0 ]; then
         echo "** error code returned from assign_taxonomy, giving up **"
         exit 1
      fi
   fi

   # make otu table 
   if [ -f $OUT_DIR/qiime_uclust/combined_rep_set_otu_table.txt  ]; then
      echo "skipping otu table as $OUT_DIR/qiime_uclust/combined_rep_set_otu_table.txt already exists"
   else
      tardis --shell-include-file $OUT_DIR/configure_qiime_env.src make_otu_table.py -i $OUT_DIR/qiime_uclust/combined_otus.txt -t $OUT_DIR/qiime_uclust/combined_rep_set_tax_assignments.txt -o $OUT_DIR/qiime_uclust/combined_rep_set_otu_table.txt  1\>$OUT_DIR/qiime_uclust/make_otu_table.stdout 2\>$OUT_DIR/qiime_uclust/make_otu_table.stderr
      if [ $? != 0 ]; then
         echo "** error code returned from make_otu_table, giving up **"
         exit 1
      fi
   fi

   # summarise taxa
   if [ -f $OUT_DIR/qiime_uclust/combined_rep_set_otu_table_L10.txt  ]; then
      echo "skipping summarise taxa as landmark $OUT_DIR/qiime_uclust/combined_rep_set_otu_table_L10.txt already exists"
   else
      tardis --shell-include-file $OUT_DIR/configure_qiime_env.src summarize_taxa.py -i $OUT_DIR/qiime_uclust/combined_rep_set_otu_table.txt -o $OUT_DIR/qiime_uclust -a -L 4,6,7,10 1\>$OUT_DIR/qiime_uclust/summarize_taxa.stdout 2\>$OUT_DIR/qiime_uclust/summarize_taxa.stderr
      tardis --shell-include-file $OUT_DIR/configure_qiime_env.src biom convert -i $OUT_DIR/qiime_uclust/combined_rep_set_otu_table.txt -o $OUT_DIR/qiime_uclust/table.from_biom_w_taxonomy.txt --to-tsv --header-key taxonomy
      if [ $? != 0 ]; then
         echo "** error code returned from summarize_taxa, giving up **"
         exit 1
      fi
   fi


   # do plots
   if [ -f $OUT_DIR/qiime_uclust/combined_rep_set_otu_table_L10.jpg  ]; then
      echo "skipping plots as landmark $OUT_DIR/qiime_uclust/combined_rep_set_otu_table_L10.jpg already exists"
   else
      tardis --shell-include-file $OUT_DIR/configure_bioconductor_env.src Rscript --vanilla ./tax_summary_heatmap.r num_profiles=30 moniker=combined_rep_set_otu_table_L10 datafolder=$OUT_DIR/qiime_uclust 1\>$OUT_DIR/qiime_uclust/tax_summary_heatmap.stdout 2\>$OUT_DIR/qiime_uclust/tax_summary_heatmap.stderr
      tardis --shell-include-file $OUT_DIR/configure_bioconductor_env.src Rscript --vanilla ./tax_summary_heatmap.r num_profiles=30 moniker=combined_rep_set_otu_table_L4 datafolder=$OUT_DIR/qiime_uclust 1\>\>$OUT_DIR/qiime_uclust/tax_summary_heatmap.stdout 2\>\>$OUT_DIR/qiime_uclust/tax_summary_heatmap.stderr
      tardis --shell-include-file $OUT_DIR/configure_bioconductor_env.src Rscript --vanilla ./tax_summary_heatmap.r num_profiles=30 moniker=combined_rep_set_otu_table_L6 datafolder=$OUT_DIR/qiime_uclust 1\>\>$OUT_DIR/qiime_uclust/tax_summary_heatmap.stdout 2\>\>$OUT_DIR/qiime_uclust/tax_summary_heatmap.stderr
      tardis --shell-include-file $OUT_DIR/configure_bioconductor_env.src Rscript --vanilla ./tax_summary_heatmap.r num_profiles=30 moniker=combined_rep_set_otu_table_L7 datafolder=$OUT_DIR/qiime_uclust 1\>\>$OUT_DIR/qiime_uclust/tax_summary_heatmap.stdout 2\>\>$OUT_DIR/qiime_uclust/tax_summary_heatmap.stderr
   fi


}

function clean() {
   rm -rf $OUT_DIR/tardis_*
}


function html_prism() {
   echo "tba" > $OUT_DIR/batch_qiime_prism.html 2>&1
}


function main() {
   get_opts "$@"
   check_opts
   echo_opts
   check_env
   configure_env
   get_targets
   if [ $DRY_RUN != "no" ]; then
      fake_prism
   else
      run_prism
      if [ $? == 0 ] ; then
         clean
         html_prism
      else
         echo "error state from batch_qiime run - skipping clean and html page generation"
         exit 1
      fi
   fi
}


set -x
main "$@"
set +x

