Thanks for downloading SynthMLM. These repository contains a few files that demonstrate the use cases of SynthMLM.

Before use, download all files and folders and put them all in the same parent folder. 

Two additional files are needed: 'nocodazole_1' and the model file 'synthmlm_analysis.jld2' which can be obtained from this link:

DATASET:
https://bham-my.sharepoint.com/personal/l_gall_bham_ac_uk/_layouts/15/guestaccess.aspx?share=IQCYxNtLeDAQS72oKgQc4K-CAapvFvztziluSq1UcCLLt8E&e=ncLM9o

And

MODEL:
https://bham-my.sharepoint.com/personal/l_gall_bham_ac_uk/_layouts/15/guestaccess.aspx?share=IQDfzZ0WgJl1TquMVPdj3Hx3AT_b4Y-7QrOfMexzT-rHilI&e=D3bYmm

To use these files, you will need a working version of Julia (tested on 1.10 and 1.12).
Packages can be obtained by activating the environment by typing 
"] activate synthmlmenv"
followed by "] instantiate".

'synthmlm base.jl' should not be run. It contains functions needed by the other two files.

'synthmlm data generation.jl' demonstrates the synthetic data generation functionality of SynthMLM and has the following functions:

generateall(j,savepp) generates j synthetic images containing all elements
generatefibres(j,savepp) generate j synthetic images containing fibres
generatebackground(j,savepp), generaterings(j,savepp), generateclusters(j,savepp), generateclustersfibres(j,savepp) are also available.

j defaults to 25 and savepp defaults to false. If savepp = true, the raw localisation coordinate files are saved too.
j and savepp can be omitted to use the default values.



'synthmlm analysis and cloning' demonstrates the dataset analysis and cloning component of SynthMLM

It has two functions:
+ analysedataset(location) to analyse a dataset stored in the folder 'location'.
+ clonedataset(j,location) to clone a dataset stored at location, creating j synthetic images.

Default parameters have the dataset location as nocodazole_1 (also in this Github repo), which when placed in the same folder as this file will allow the code to analyse and/or clone this data.

You can analyse/clone your own data by replacing this folder. It will need to be formated in the way nano-org formats its datasets.
To convert an SMLM dataset into this format, upload the dataset to nano-org (nano-org.bham.ac.uk) and click 'Download Gridded Images' after processing.
Place these gridded images into a folder and replace 'location' in the function call above with this directory.

Analysing and producing descriptor-matched synthetic clones of datasets is fully integrated into nano-org.bham.ac.uk!
