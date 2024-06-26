set -e
#!/bin/bash

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

project=`python cache.py PROJECT:`

#Creating output directory for html files
path_html="html_report_"$project
mkdir -p $path_html

#Creating html file with links to all output files
#output_name=$path_html"/"$project"_file_list.html"
#output_name=$path_html"/"$project"_file_list.html"

output_name=$project"_file_list.html"
summarized_output_name=$project"_summarized_file_list.html"
project_path="../data/endpoints/"$project"/"

cd $path_html
tree $project_path -H $project_path -T $project -o $output_name
tree $project_path -d -H $project_path -T $project -o $summarized_output_name
cd -

#Creating file continaing recovered cells for each sample.
#python web_scraper.py $project

#Creating .Rmd files for each sample
python html_sample_script.py $project

#getting list of samples as a string (for iteration)
#samples=`python sample_string.py $project`

echo "\n\nYou will need to create (e.g., Knit) html files from the Rmd files for each sample ($samples) using Rstudio.\n"
echo "After the html files for each sample have been created, use the \"Run Document\" command (in Rstudio) on the Interactive_report.Rmd file to create the interactive Shiny app.\n"
echo "The html output is stored here: $path_html"

#echo "Creating .html files from .Rmd scripts (for each sample)"
#for s in $samples
#do
#	html_file="html_report/"$s".Rmd"
#	echo $s
#	echo $html_file
#	Rscript -e 'library(rmarkdown); rmarkdown::render("html_report/"$s".Rmd", "html_document")'
#	#Rscript -e 'library(rmarkdown); rmarkdown::render(paste0("html_report/", s, ".Rmd"), "html_document")'
#	#Rscript -e 'library(rmarkdown); rmarkdown::render("html_report/"$s".Rmd", "html_document")'
#done
#rmarkdown::render("htmltest.Rmd", "html_document")

#echo $project_path
#project_upper="${project^}"
#project_upper="${project~}"
