
How to use the epigenomics model:

Download the required additional data files at: https://www.dropbox.com/sh/yhvkps4g1uucryn/AAB7Kp8UengB9F_YfNl284D7a?dl=0
Put epigenomic_regluation.R in the same directory as these data files.  Maintain the structure of all subdirectories.

Open epigenomic_regulation.R in Rstudio or any text editor.
This R script requires two files from user input:
	File 1) default: "user_input/lung_cancer_meth450.Rdata".  This is an Rdata file containing a data frame of differential methylation results as processed by limma.  Load the default file into R to see an example.  Name this variable "cpgs_limma"
	File 2) default: "user_input/transcript_regulation.Rdata".  This is an Rdata file containing a matrix of mRNA transcripts and whether they are upregulated or downregulated.  Column 1 is the transcript name, column 2 is the upregulated/downregulated example.  Load the default file into R to see an example.  Name this variable "transcript_regulation".

Once these two file names have been set, go ahead and run the R script.
A resultant csv file will be produced, named finalModel_forWEKA.csv.
This csv file is immediately ready for feature selection and cross-fold model evaluation in weka.
Optionally, the user can split the data table in pieces, witholding a portion as a "testing set".

USING WEKA:

WEKA can be downloaded at: http://www.cs.waikato.ac.nz/ml/weka/downloading.html

Launch WEKA.
Click "Explorer".
Load the training set (finalModel_forWEKA.csv, by default) into WEKA (click Open file...".
REMOVE the first entry titled "transcript_id".
Click the "Select attributes" tab at the top.
	Here you can choose your feature selector and any parameters.  For our model, we use the "ReliefFAttributeEval" attribute evaluator and "Ranker" search method.  We set the Ranker numToSelect parameter to 67, leaving everything else as default.
Click start.
When feature selection finishes, right click the result in "Result list" and click "Save reduced data...".  Save this file with the .arff extension.
Go back to the "Preprocess" tab and load the new reduced data set into WEKA.
Go to the "Classify" tab.
Here you can choose your classifier and set any parameters.  For our model, we use "RandomForest" and set numTrees to 100.
IF YOU ARE RUNNING CROSS VALIDATION:
	The model is ready to run, click start.

IF YOU CREATED A HOLDOUT TESTING SET:
	In order for weka to run, the column names of your training and testing set must be identical.
	Unfortunately, reliefF/ranker feature selection reorders column names.  You must go back into your testing set and change the order of columns to match that of your training set.  This can be done in R.  You must also remove the "transcript_id" column in R, as there will be no option to do it in WEKA.
	Once your testing set is ordered correctly, click "Supplied test set" and choose your testing set.
	The model is ready to run, click start.
	
SAMPLE DATA:
We provide 4 sample files in the sample_results folder:
The two files that have not gone through feature selection are full models - feel free to apply different feature selection methods to compare results.
The two files that have gone through reliefF feature selection are ready for cross validation or evaluation of a testing set.  The testing set is already ordered in the manner described above, so feel free to load both training and testing into R to evaluate results.

We also provide 3 sample files in the user_input folder:
lung_cancer_meth450 is our differential methylation data, as processed by limma
transcript_regulation is the set of transcripts and their regulation information we use to build the training set
sample_testingset is the set of transcripts and their regulation information we use to build the testing set