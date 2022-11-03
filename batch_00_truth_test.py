import filecmp, sys, subprocess

if sys.argv[1].lower() in ["y","yes"]:
  open_plots=True

else:
  open_plots=False
 
modules={}
modules['ipa_output'] = ['grouped_optimal_pointings.txt']
modules['k_clean_output'] = ['k_list_mod.txt']
modules['k_make_output'] = ['k_list_raw.txt']
modules['m_make_output'] = ['single_m_list.txt']
modules['m_sort_output'] = ['single_list_fom2.txt']
modules['m_pick_output/pointing_100/'] = ['observed_targets.cat','pointing_summary.txt','shutter_mask.csv','slitlet_usage_stats.txt','target_list.txt']

truth_dir = './trial_00_ref/'
test_dir = './trial_00/'

for m in modules:
  for f in modules[m]:
    print("Comparing "+truth_dir+m+"/"+f+" and "+test_dir+m+"/"+f+" for T or F match:")
    print(filecmp.cmp(truth_dir+m+"/"+f,test_dir+m+"/"+f))
    filecmp.clear_cache()
    
     

if open_plots:
  for m in ['ipa_output','m_pick_output/pointing_100', 'm_check_output/pointing_100']:
    subprocess.call("open "+test_dir+m+"/*.ps", shell=True)



quit()

