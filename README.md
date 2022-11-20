# eMPT

 The eMPT Software and User Guide is contained in file eMPT_user_guide_release.pdf. For the full introduction to and in-depth description of its algorithms and design, refer to the submitted draft paper contained in eMPT_paper_AA_submit.pdf.
  
QUICK START GUIDE
 
Run the eMPT using the provided test data set by executing the shell script 'batch_00.sh' on the command line:

% more batch_00.sh:

#! /bin/bash
./ipa 00.conf
./k_make 00.conf
./k_clean 00.conf
./m_make 00.conf
./m_sort 00.conf
./m_pick 00.conf
./m_check 00.conf
./m_check_regions 00.conf

% ./batch_00.sh

This script runs each of the core eMPT Fortran modules in the required strict sequence shown, that both read and write back to the controlling configuration file 00.conf that contains all parameter settings, both default and user-provided. When all modules have finished running, the configuration file serves as both a record of the trial run of the eMPT software, and also the place to continue to make manual updates to input parameter values -- both at the top of the file and in several places further down in the file where specific module output can be examined and changed, in which case all downstream modules should then be re-run.


**The required user input parameters include:**

 - RA, Dec & PA_AP of nominal pointing in decimal degrees [RA 0:360 - Dec -90:90 - PA_AP 0:360] 
 - name of the eMPT-formatted target input catalog, which must contain Target ID/No., RA, Dec, and Priority Class columns, in that order; all additional    (rightward) columns will be ignored. Header line(s) and all non-catalog lines must start with '#'. Target  priority class assignments must include '1' 
   in order for the software to be able to identify the highest priority class.
 - the number of the trial, from 00:100
 - the number of dithers, if other than the default value of 1


To become acquainted with the various required and default parameter settings of the eMPT software, as well as the plentiful module output printed to the screen, printed to the controlling configuration file, and that is deposited in ASCII and PS files to the subdirectories it creates matching each investigative trial run, it is recommended to continue manually running the software in this way, in single trials.   

For experienced users who are comfortable with Python, a Python template wrapper script is provided, 'empt_explore.py', that was designed to be opened and directly edited by the user to enable a quick and automatic execution of a batch of many exploratory trials to be compared and contrasted. The Python wrapper script may also be utilized for single runs of the software, as described in the introduction section of the script. 






