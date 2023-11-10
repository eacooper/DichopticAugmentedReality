# DichopticAugmentedReality
Data and code associated with manuscript "The effect of interocular contrast differences on the appearance of augmented reality imagery"
ACM Transactions on Applied Perception, in press

Folder: data --> for both experiments are in the data folder (Expt 1 and Expt2_AR subfolders)
- data structure fields
	- reference = whether the reference stim was the top one or the bottom one
	- adjust = the adjustable one was top or bottom
	- stim = columns of stimulus info 
	      --col1 is stimulus type (also listed in the .pattern field)
	      --col2 is left eye contrast
	      --col3 is right eye contrast
              --col4 (Expt 1 only) is surround contrast (average of left and right eye)
	- pattern = what stimulus pattern was coded by stim field col1
	      -- example: for Expt 1, stim col1=1 means grating, col1=5 was bandpassed noise
	- resp = responses
	      --col1 is matching result
	      --col2 is exact match response (1 = yes, 2 = no)
	      --col3 to col7 are perceptual effect responses: brightness, contrast, luster, rivalry, depth
		  ---1 = top selected, 2 = bottom selected, 4 = same/neither, 5 = unsure
*********************************************************************************************
Folder: Expt2 matching analysis --> matching task analysis code for Expt2 only, see other repo for Expt1
- code files for making Figure 8: 
	- plot_ARmatch_byICR.m --> determine high contrast stim weight for each interocular contrast ratio condition 
	- plot_ARmatch_byPattern.m --> determine high contrast stim weight for each stimulus pattern condition
	- genBino.m --> helper function for the simple weighted combination model
- subfolder: eye dominance --> eye dominance analysis based on matching
	- parseByEye.m --> determine eye dominance status for each subject, and save result as 'ARexpt_eyedom.mat'
	- ARmatch_weight.m --> determine the weight for dominant eye high contrast trials vs. nondominant eye, and make plots
- subfolder: R stats --> processed data files and R code to run ANOVA/t-test for the matching task, Table 4

*********************************************************************************************
Folder: Perceptual question analysis --> scripts for analyzing perceptual question responses for both Expt 1 and 2
- main plot and analysis code
	- plot_exactMatch_Expt*.m --> plot Figure 6 and Figure 8 (proportion of exact match)
	- doGLME_Match_Expt*.m --> run mixed effect logistic regression model, Table 1, 2, 5, 6 
	- plot_propEffect_Expt*.m --> plot Figure 7 and Figure 9 (proportion of each perceptual effect)
	- doGLME_Effect_Expt*.m --> run logistic regression for different perceptual effects, Table 3, Table 7
- misc. analysis
	- check_luster_rivalry.m --> check if dichoptic stim was selected to be lustrous or rivalrous
	- check_unsure.m --> check the number of unsure responses
	- check_numEffectPerTrial.m --> check effect co-occurence
	- check_expt1_vs_2_ttest.m --> compare expt 1 and 2 for the 5 effects (uses expt*_effect.mat saved from plot_propEffect_Expt*.m)