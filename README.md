# MULTICOV-AB_Publication
R code related to the MULTICOV-AB publication in Nature Communications:
"Exploring beyond clinical routine SARS-CoV-2 serology using MultiCoV-Ab to evaluate endemic coronavirus cross-reactivity"

The code in the file "Analysis_submission.R" used to generate the graphs used in the paper.
In general graphs are exported to an output folder in .pdf format, 
from where they were further processed into the finalized figures using a vector graphics editing software such as InkScape.

R version used was version 3.6.1 (2019-07-05) -- "Action of the Toes" within RStudio 1.2.5001

Required are an input folder at "/input" with the file containing the raw measurement data 
"Raw_Data_for_analysis.csv" with annotations, the raw data for the parallelism and Quality Control performance 
"Parallelism.csv.txt", "QC_IgG.csv.txt" and "QC_IgA.csv.txt"
