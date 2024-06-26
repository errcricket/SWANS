set -e

#ANNOTATIONS-------------------
#unpack pangeloDB and cellmarker
reference_data="reference_data/"
cd $reference_data
pang="panglaoDB/"
cell="CellMarker/"

if [ ! -d "$pang" ]; then
	tar -xvf "panglaoDB.tar"
fi

if [ ! -d "$cell" ]; then
	tar -xvf "CellMarker.tar"
fi

cd -
#------------------------------

sample_list_file="samples.sample_list"

#confirm sample list file exists
if [ ! -e "$sample_list_file" ]; then
	echo "------------------------------"
	echo "AUCHTUNG AUCHTUNG AUCHTUNG"
	echo "------------------------------"
	echo "$sample_list_file does not exist."
	echo ""
	echo "Create $sample_list_file according to the instructions & run this script again."
	exit 1;
fi

config_file="configs/local_configs.yaml"

#confirm config file exists
if [ ! -e "$config_file" ]; then
	echo "------------------------------"
	echo "AUCHTUNG AUCHTUNG AUCHTUNG"
	echo "------------------------------"
	echo "$config_file does not exist."
	echo ""
	echo "Create $config_file according to the instructions & run this script again."
	exit 1;
fi

project=`python cache.py PROJECT:`

echo "Creating project directory (if it does not exist)"
path="data/endpoints/"$project

mkdir -p $path
python setup.py $project #this will make project_name/sample_name/10X for each sample
python make_symbolic_links.py $project #this will make a bash script
sh make_symbolic_links.sh #this will make the symbolic links in each 10X directory
echo "Setting up project directory"

#CHECK FOR MISSING DATA
missing=`python check_setup.py $project`  #this check to see if you have 10X before running Snakefile

for i in ${missing[*]}; do
	if [[ $i == 0 ]]; then 
		echo 'You have not placed any data in your 10X folders'
		echo 'Add data to your 10X folders and run this again'
		exit 1;
fi
done

date=$(date '+%a_%d_%Y')

filename=$date"_Rsession_info.txt"

#sessionInfo
#Rscript src/scripts/session.R $rpath > $filename

image="workflow_targets.png"
#snakemake --forceall -r -s Snakefile --dag | dot -Tpng > $image
#snakemake -j300 --snakefile Snakefile --printshellcmds --dryrun 
#snakemake -j300 --snakefile Snakefile --printshellcmds --touch --dryrun --reason --forceall

# run 3 times: 1) sig 2) sc_objects 3)
for  i in 1 2 3 
do
	#snakemake -j300 --snakefile Snakefile --printshellcmds --dryrun 
	snakemake -j300 --snakefile Snakefile --printshellcmds

	# create html file of all output
	if [[ i -eq 2 ]]
	then
		sh create_htmlReport.sh
	fi

	#echo $i
	#echo $?
	if [[ ! $? = 0 ]]
	then
		exit
	fi
done
