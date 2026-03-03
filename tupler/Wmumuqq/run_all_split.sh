#!/bin/bash

#rm -rf local_run*

nargs=$#
if (( nargs < 2 )); then
    echo ERROR: Please provide at least two arguments!
    echo Example: ./run_all_split.sh nFilesPerJob [data_or_mc=data/MC] [LFNs_file] [options]
    exit
fi

nFilesPerJob=$1
data_or_mc=$2
LFNs_file_full=sample_list/LFNs_42912016_20GeVCtau0ps.txt
options=options.py

if (( nargs > 2 )); then
    LFNs_file_full=$3
fi

if (( nargs > 3 )); then
    options=$4
fi

LFNs_file=$(basename "$LFNs_file_full")

dir_suffix=$(date +"%Y%m%d_%H%M%S")
dir_prefix="${options%.py}"

work_dir=local_run_"${LFNs_file%.*}"_${dir_prefix}_$dir_suffix
hadd_filename="${LFNs_file%.*}"_$dir_suffix

echo INFO: running with arguments: $nargs, $options, $nFilesPerJob, $data_or_mc, $LFNs_file, $work_dir

mkdir $work_dir
cp $LFNs_file_full $work_dir/
cd $work_dir
cp ../*.py .

##### MC part
if [[ "$data_or_mc" == "MC" ]] || [[ "$data_or_mc" == "mc" ]]; then
    echo "INFO: running MC davinci"
   
    line_counter=0
    counter=0

    while IFS= read -r line; do
        if (( line_counter % nFilesPerJob == 0 )); then
            ((counter++))
            cp ../info_MC_template.txt info_MC_${counter}.yaml
            echo "ntuple_file: dv_ntuple_mc_${counter}.root" >> info_MC_${counter}.yaml
            echo "input_files:" >> info_MC_${counter}.yaml
        fi
        echo "  - $line" >> info_MC_${counter}.yaml
        ((line_counter++))
    done < $LFNs_file
    
    echo "Split into $counter jobs"
    for ((iJob=1; iJob<=counter; iJob++)); do
        echo "Running MC job $iJob"
        echo "nohup lb-run DaVinci/v66r6 lbexec $dir_prefix:main info_MC_${iJob}.yaml > log_nohup_MC_${iJob}.log 2>&1 &"
        nohup lb-run DaVinci/v66r6 lbexec $dir_prefix:main info_MC_${iJob}.yaml > log_nohup_MC_${iJob}.log 2>&1 &
    done
    #prepare hadd script
    echo "hadd -k -f $hadd_filename.root dv_ntuple_mc_*.root" >> cmd_hadd.sh
    chmod 744 cmd_hadd.sh
fi

##### data part
if [[ "$data_or_mc" == "data" ]]; then
    echo "INFO: running data davinci"
   
    line_counter=0
    counter=0

    while IFS= read -r line; do
        if (( line_counter % nFilesPerJob == 0 )); then
            ((counter++))
            cp ../info_data_template.txt info_data_${counter}.yaml
            echo "ntuple_file: dv_ntuple_data_${counter}.root" >> info_data_${counter}.yaml
            echo "input_files:" >> info_data_${counter}.yaml
        fi
        echo "  - $line" >> info_data_${counter}.yaml
        ((line_counter++))
    done < $LFNs_file
    
    echo "Split into $counter jobs"
    for ((iJob=1; iJob<=counter; iJob++)); do
        echo "Running data job $iJob"
        echo "nohup lb-run DaVinci/v66r6 lbexec $dir_prefix:main info_data_${iJob}.yaml > log_nohup_data_${iJob}.log 2>&1 &"
        nohup lb-run DaVinci/v66r6 lbexec $dir_prefix:main info_data_${iJob}.yaml > log_nohup_data_${iJob}.log 2>&1 &
    done
    #prepare hadd script
    echo "hadd -k -f $hadd_filename.root dv_ntuple_data_*.root" >> cmd_hadd.sh
    chmod 744 cmd_hadd.sh
fi


cd ..

