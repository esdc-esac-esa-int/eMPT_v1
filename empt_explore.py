import numpy as np 
import sys,time,os,subprocess

def execute(name, confname):
        
    shell_command =  './'+name+' '+confname
    c = subprocess.call(shell_command, shell=True)

    if int(c) > 0:
      raise IOError('The Fortran '+name+' module failed to execute.')


def parse_config_file(confname):

       inpars={} 
       pars=[]
       pars_ind=[]
       contents_arr = []
       
       f = open(confname, "r")
       contents = f.readlines()

       for i,line in enumerate(contents):
         
         if not line.startswith("#"):
            pars.append(line)
            pars_ind.append(i)
           
       inpars["n_trial"]      = pars[0] 
       inpars["disperser"]    = pars[1]
       inpars["raccx"]        = pars[2].split(" ")[0]
       inpars["raccy"]        = pars[2].split(" ")[1]
       inpars["sthresh"]      = pars[3]
       inpars["catfile"]      = pars[4]
       inpars["fitsref"]      = pars[5]
       inpars["segmap"]       = pars[6]
       inpars["cra"]          = pars[7].split(" ")[0]
       inpars["cdec"]         = pars[7].split(" ")[1]
       inpars["cpa_ap"]       = pars[7].split(" ")[2]
       inpars["szone_x"]      = pars[8].split(" ")[0]
       inpars["szone_y"]      = pars[8].split(" ")[1]
       inpars["true_ang"]     = pars[9]
       inpars["max_c_score1"] = pars[10]
       inpars["max_c_score2"] = pars[11]
       inpars["n_dither"]     = pars[12]

       return np.array(contents),inpars,np.array(pars_ind)
    

def edit_config_file(confname,n_trial,contents_arr,inpars,pars_ind,user_args=None):

         if user_args:
           for k,v in user_args.items():
             print("Reading user input parameter values: ")
             print(k,v)   
             inpars[k] = v

         contents_arr[pars_ind] = [inpars["n_trial"],inpars["disperser"],str(inpars["raccx"])+" "+str(inpars["raccy"]),inpars["sthresh"],inpars["catfile"],inpars["fitsref"],inpars["segmap"],str(inpars["cra"])+" "+str(inpars["cdec"])+" "+str(inpars["cpa_ap"]),str(inpars["szone_x"])+" "+str(inpars["szone_y"]),inpars["true_ang"],inpars["max_c_score1"],inpars["max_c_score2"],inpars["n_dither"]]

         if int(str(n_trial))<9:
           n_trial_str = "0"+str(int(n_trial))
         else:
           n_trial_str = str(n_trial)
           
         updated_config_filename = n_trial_str+".conf"
         updated_config_file = open(updated_config_filename, 'w')

         for line in contents_arr:
            if not str(line).endswith("\n"):
              updated_config_file.write(str(line)+"\n")
            else:
              updated_config_file.write(line)

         updated_config_file.close()  
         return updated_config_filename

# -----------------------------------------------------------------------------------------------------------------------------------------

if __name__ == "__main__":
  
 # The eMPT is written in Fortran due to its superior ability to handle computationally intensive tasks, very quickly.
 #
 # For users who are comfortable with Python and who wish to run the software automatically in many repeated
 # trial runs to explore the eMPT parameters space, we provide this template script with an example 'for' loop
 # that can easily be adapted to the user's specific needs.
 #
 # If you choose to run this script for a single trial --  as opposed to editing the 'for' loop starting at line 225
 # to run through many trials, with parameter settings changing with each pass through the loop -- then EITHER enter the
 # required parameter values plus any additional preferred changes to template.conf, provided with the eMPT package; OR
 # specify them on the command line, as shown below.
 #
 # The script will look to both template.conf and the command line for user input for a single run of the software, and for
 # the starting parameters for a batch of user-specified runs; if both are provided, the information read from the command
 # line will take precedence.

 # =========================================================================================================================================
 #
 # Examples of usage:
 #
 # OPTION 1: Batch mode
 #
 #  Uncomment and edit the example 'for' loop starting on line 225 of this script to define the batch of trials you would like to run.
 #  Define the starting set of parameters on the command line or in template.conf, both of which are searched by the script for user input
 #  (with the command line taking precendence where entries have been updated in both places). To take advantage of this setup, keep default
 #  settings, or the set of parameters that won't change from trial to trial, in template.conf, and only update the parameter values
 #  that are being explored on the command line or inside the script itself.
 #
 #  % python ipa_empt_template_script.py mode=batch
 #
 # OPTION 2 : Single run interactive mode (command line parameter entry OR in template.conf)
 #
 #  Execute the script on the command line with required parameters (plus any others to change from default settings), ensuring
 #  parameter names are entered exactly as shown in the "in_pars" dictionary below.
 #
 #  % python ipa_empt_template_script.py mode=interactive cra=53.159 cdec=-27.80 cpa_ap=194.0 n_trial=1 n_dither=3 catfile=my_msa_targets.cat
 #
 #
 #  Execute the script on the command line without arguments, setting parameters preferences instead in template.conf read by the script.
 #
 #   % python ipa_empt_template_script.py mode=interactive
 # 
 # ==========================================================================================================================================

 confname="template.conf"                                     
 in_pars = parse_config_file(confname)[1]       # read default parameters from template configuration file

 #in_pars={}
 in_pars["mode"]         = 'batch'              # 'interactive' or 'batch' mode
 #in_pars["confname"]     =  "00.conf"          # trial log (Fortran-formatted configuration file) 
 #in_pars["n_trial"]      =  "1"                # Unique ID number of run nn [00:99]  
 #in_pars["disperser"]    = "PRISM/CLEAR"       # NIRSpec disperser employed [PRISM/CLEAR,G140M/F070LP,G140M/F110LP,G235M/F170LP,G395MF290LP,G140H/F070LP,G140H/F110LP, G235H/F170LP, G395HF290LP]
 #in_pars["raccx"]        = "0.343"             # Acceptance Zone half extents in MSA X and Y in units of MSA facets
 #in_pars["raccy"]        = "0.420"             # (shaved full open area) 
 #in_pars["sthresh"]      = "3.5"               # Minimum vertical spectral separation threshold in MSA shutter Y facets
 #in_pars["catfile"]      = "test_trial_00.cat" # Name of eMPT-formatted target input catalog
 #inpars["fitsref"]       = "none"              # N/A (unused in release v1)  
 #inpars["segmap"]        = "none"              # N/A (unused in release v1)
 #in_pars["cra"]          = "53.159"            # RA, Dec & PA_AP of nominal pointing in decimal degrees [RA 0:360 - Dec -90:90 - PA_AP 0:360]  
 #in_pars["cdec"]         = "-27.80"
 #in_pars["cpa_ap"]       = "194"
 #in_pars["szone_x"]      = "25.0"              # Allowable absolute deviation from nominal pointing in arcsec in 
 #in_pars["szone_y"]      = "25.0"              # X and Y on MSA [Default +/-25.0 x +/-25.0]   
 #in_pars["szone_y"]      = "80.7829"           # Angle between JWST spacecraft velocity vector and the nominal NIRSpec pointing in
                                                # decimal degrees to be needed to calculate the change in MSA plate scale caused by
                                                # Differential Velocity Aberration.
 #in_pars["max_c_score1"] = "100 100 100"       # Maximum acceptable contamination score for Priority Class 1 targets 
 #in_pars["max_c_score2"] = "1 2 2"             # and Priority Class 2 and larger targets                 
 #in_pars["n_dither"]     = "1"                 # Number of dithered pointings sought in above search box [1, 2 or 3]



 global_var = ("mode","confname","n_trial","disperser","raccx","raccy","sthresh","catfile","fitsref","segmap","cra","cdec","cpa_ap","szone_x","szone_y","max_c_score1","max_c_score2","n_dither")
 required_args = ('mode','cra','cdec','cpa_ap','catfile','n_trial','n_dither')

 arglist  = sys.argv[1:]

 if len(arglist) > 0:

   cl_args = [a.split("=")[0] for a in arglist if "=" in a]
   cl_arg_vals = [a.split("=")[1] for a in arglist if "=" in a]
  
   for a in required_args:
     if not a.lower() in cl_args:
         Warning("\nRequired input parameter "+a+" not found (or entered incorrectly) on command line; taking the value you may have entered directly to this script, or directly to template.conf (if you failed to enter it anywhere, the example value from the original template.conf file will be used.)")
         if not a.lower() in in_pars: 
           raise IOError("\nPlease enter a starting value for the required input parameter "+a+". \nRequired input parameters: nominal pointing coordinates and roll angle (cra, cdec, cpa_ap); \nthe name of the file containing the input target catalog (catfile); \nthe ID number for the current trial (n_trial, 0:99); \nand the number of dithers (n_dither) for the configuration (1 for a single pointing, 2 or 3).\n All remaining parameter values will remain at their default values, listed in the par_entry.py module.\n")

   for i,val in enumerate(cl_args):
     if not val.lower() in global_var:
       raise IOError('\nParameter name '+val+' not recognized; choose from among the provided list of recognized parameters: \n\n'+str(global_var)+'\n')
     else:
       in_pars[val] = cl_arg_vals[i]
       
 else:
      print("\nNo command line arguments found. Taking starting input parameter values from template.conf.\n")
      

 updated_config = edit_config_file(confname,in_pars["n_trial"],parse_config_file(confname)[0],in_pars,parse_config_file(confname)[2])

 print("\n\n\nPreparing to run the eMPT with the following input parameter values and configuration file "+updated_config+":\n")

 for k,v in in_pars.items(): 
   print(k+": "+v)

 time.sleep(10)
 

 if not in_pars["mode"]=="batch":

     mod_message={}
     mod_message["ipa"]     = "You are free to edit the candidate pointing list automatically selected by the ipa, contained in the newly created 00.conf file -- that is actively being updated during this run of the modules -- to see if there are any you would like to remove (by commenting out) to avoid having them processed by the downstream modules. When this interactive run finishes, you can go back and edit the ipa pointing list in 00.conf and re-run only the subsequent modules accordingly."
     mod_message["k_clean"] = "For dithers (i.e. n_dither > 1) , the default placement sequence shown in the k_clean output assumes that you wish to maximize target commonality between all pairs of dithered pointings (and thereby the average exposure time of the observed targets) by first placing all the common targets in sequence of increasing Priority Class, followed by the unique targets in the same sequence. Should you instead wish to maximize the total number of different targets covered by a pair or triplet of pointings, then the two/three rows should be swapped in the placement sequence."
     mod_message["m_make"]  = "You are free to edit the weights in the m_make output to better reflect the scientific objectives of the program in question, keeping in mind the placement sequence specified above."
     mod_message["m_sort"]  = "You are free to select a different preferred optimal solution than the one output by m_sort by 'commenting out' the proposed solution with a hashtag at the start of the line and 'uncommenting' a different line in the configuration file -  or selecting an altogether different preferred solution by copying and pasting the relevant entry in the complete ranked FOM listings produced by the m_sort module."
             
     for mod in ["ipa","k_make","k_clean","m_make","m_sort","m_pick","m_check","m_check_regions"]:
         
         print("\nRunning the "+mod+" module...\n")
         time.sleep(4)
         
         if mod in mod_message:
           print(mod_message[mod])
           if  mod in ["ipa","k_clean","m_sort"]:
             time.sleep(20)
           else:
             time.sleep(10)  

         execute(mod, confname=updated_config)
         

 else:

     print("\nSince you have specified (or left the default) 'batch' run mode setting, please uncomment and customize the portion of this script starting at line 225 to define your batch runs of the eMPT software; an example is provided to get you started. Then, simply re-run the script in batch mode.")
     time.sleep(10)
     
   
     # --------------------------------------------------------------------------------------------------------------------------------------------
     # 
     # Example 'for' loop that explores two different nominal (i.e. starting) (Ra,Dec) locations around which to search for optimal pointings, at 5
     # different roll angles. Edit the parameter names accordingly, add or subtract a nested loop, etc., to loop through your eMPT parameters in
     # which you are interested in a similar fashion, to explore the eMPT parameter space and produce the most optimal MSA masks for your science program.
     # 
     # --------------------------------------------------------------------------------------------------------------------------------------------
     #
     #user_pars = {}
     #trial_map={}  


     #test_nom_pos = [(53.14187071428572,-27.809516857142857), (53.141582,-27.81504)]
     #test_pa_v3 = range(25,225,5)
     #test_pa_ap = np.array([float(p) for p in test_pa_v3]) + 138.492


     #ntrial=1        #00:99 allowed
     #ntrial_str="01"
    
     
     #for pos in test_nom_pos:
     
     #  cra,cdec= pos[0],pos[1]
     
     #  for i,pa in enumerate(test_pa_ap):

     #    if ntrial<10:
     #      ntrial_str = "0"+str(ntrial)
     
     #    else:
     #      ntrial_str=str(ntrial)
       
     #    trial_map[ntrial] = (cra,cdec,pa)

     #    user_pars['cra']    = str(cra)
     #    user_pars['cdec']   = str(cdec)
     #    user_pars['cpa_ap']  = str(pa)
     #    user_pars['n_trial'] = ntrial_str
         
     #    updated_config = edit_config_file(confname,ntrial_str,parse_config_file(confname)[0], parse_config_file(confname)[1], parse_config_file(confname)[2], user_pars)
  
     #    print("\nStarting trial "+ntrial_str+" with RA, Dec, PA_AP: "+str(trial_map[ntrial]))

     #    for mod in  ["ipa","k_make","k_clean","m_make","m_sort","m_pick","m_check","m_check_regions"]:
     #      execute(mod, confname=updated_config) 

     #    if ntrial<99:
     #      Warning("Batch runs may not exceed 100 trials. Quitting after next trial finishes.")

     #    else:
     #       raise ValueError("Quitting...batch runs may not exceed 100 trials.")

     #    ntrial+=1
 
 
     #for k in trial_map:
     #   print("trial "+str(k)+": "+str(trial_map[k])+" PA_AP,RA,Dec")


            

quit()

