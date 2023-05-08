from astropy.io import fits as pyfits
import numpy as np

# Access the most current MSA shutter status FITS tables from the Instrument Model directory:

with open('model.conf', 'r') as model_config_file:
  model_config_file_lines = [line for line in model_config_file.readlines() if not str(line).startswith("#")]
model_config_file.close()

msa_chk_fitsobj = pyfits.open(model_config_file_lines[1].strip())       # NRS_MSOP_CHK_2.0.1_20220312T0000.msl
msa_vig_fitsobj = pyfits.open(model_config_file_lines[2].strip())       # NRS_MSOP_VIG_1.0.0_20160106T0000.msl
msa_shortmask_fitsobj = pyfits.open(model_config_file_lines[3].strip()) # NRS_MSOP_SM_3.0.0_20220326T0000.msl


esa_map = {}

for k in range(1,5):
 
  # Initiate quadrant shutter values as all operable (value 2): 
  esa_map['Q'+str(k)] = np.array([2]*171*365)
  
  # Read the lists of inoperable shutters to now assign to the shutter map: 
  # NO = i+(j−1)×365, where (i,j) respectively span the intervals [1, 365] and [1, 171]

  msa_chk_shutter_no = (msa_chk_fitsobj['Q'+str(k)]).data.field('NO')  
  
  if not len(esa_map['Q'+str(k)]) == len(msa_chk_shutter_no):
      raiseValueError("Number of elements in output esa_map and 'NO' field of NRS_MSOP_CHK MSL files does no match.")
  else:  
          
    # Grab the shutter values:  
    msa_chk_shutter_val = np.array((msa_chk_fitsobj['Q'+str(k)]).data.field('STATUS'))
    
    # Locate the indexes in the array of status values where the following conditions apply, and later
    # use to apply to the esa_map with matching indexes:
    msa_chk_failed_closed_perm = np.where(msa_chk_shutter_val == 2)[0] 
    msa_chk_failed_closed_inter = np.where(msa_chk_shutter_val == 4)[0]
    msa_chk_failed_closed_vign = np.where(msa_chk_shutter_val == 8)[0]
    msa_chk_failed_open_perm = list(np.where(msa_chk_shutter_val == 256)[0]) + list(np.where(msa_chk_shutter_val == 258)[0])
    print(msa_chk_failed_open_perm)
    msa_chk_failed_open_inter = np.where(msa_chk_shutter_val == 512)[0]

    #integer check(4,365,171),short(4,365,171),fstop(4,365,171),esa_map(4,365,171)
    #  call readtable(filename_chk,check)
    #		if (check(k,i,j).eq.0)   esa_map(k,i,j)=2 !(functional)
    #		if (check(k,i,j).eq.2)   esa_map(k,i,j)=0 !(failed closed - permanent)
    #		if (check(k,i,j).eq.4)   esa_map(k,i,j)=0 !(failed closed - intermittent)
    #		if (check(k,i,j).eq.8)   esa_map(k,i,j)=0 !(failed closed - vignetted)
    #		if (check(k,i,j).eq.256) esa_map(k,i,j)=1 !(failed open - permanent)
    #		if (check(k,i,j).eq.258) esa_map(k,i,j)=1 !(failed open - permanent)
    #		if (check(k,i,j).eq.512) esa_map(k,i,j)=1 !(failed open - intermittent)

  msa_shortmask_shutter_no = (msa_shortmask_fitsobj['Q'+str(k)]).data.field('NO')
 
  if not len(esa_map['Q'+str(k)]) == len(msa_shortmask_shutter_no):
      raiseValueError("Number of elements in output esa_map and 'NO' field of NRS_MSOP_SM MSL files does no match.")
  else: 
    
    msa_shortmask_shutter_val = np.array((msa_shortmask_fitsobj['Q'+str(k)]).data.field('STATUS'))
    msa_shortmask_failed_closed = np.where((msa_shortmask_shutter_val >= 16) & (msa_shortmask_shutter_val <= 255))[0] 
                                 #np.where(msa_shortmask_shutter_val > 0)[0]
    #call readtable(filename_sht,short)
    #		if ((short(k,i,j).ge.16).and.(short(k,i,j).le.255))   esa_map(k,i,j)=0 !(failed closed - masked off)
  
  msa_vig_shutter_no = (msa_vig_fitsobj['Q'+str(k)]).data.field('NO')
  
  if not len(esa_map['Q'+str(k)]) == len(msa_vig_shutter_no):
      raiseValueError("Number of elements in output esa_map and 'NO' field of NRS_MSOP_VIG MSL files does no match.")
  else:
    
    msa_vig_shutter_val = (msa_vig_fitsobj['Q'+str(k)]).data.field('STATUS')
    msa_vig_failed_closed = np.where(msa_vig_shutter_val == 8)[0]
 
 
  # Check
  print('Q'+str(k)+': Found '+str(len(msa_chk_failed_closed_perm))+' shutters with status: NRS_MSOP_CHK failed closed permanent')
  print('Q'+str(k)+': Found '+str(len(msa_chk_failed_closed_inter))+' shutters with status: NRS_MSOP_CHK failed closed intermittent')
  print('Q'+str(k)+': Found '+str(len(msa_chk_failed_closed_vign))+' shutters with status: NRS_MSOP_CHK failed closed vignetted')
  print('Q'+str(k)+': Found '+str(len(msa_chk_failed_open_perm))+' shutters with status: NRS_MSOP_CHK failed open permanent')
  print('Q'+str(k)+': Found '+str(len(msa_chk_failed_open_inter))+' shutters with status: NRS_MSOP_CHK failed open intermittent')
  print('Q'+str(k)+': Found '+str(len(msa_shortmask_failed_closed))+' shutters with status: NRS_MSOP_SM failed closed masked off')
  print('Q'+str(k)+': Found '+str(len(msa_vig_failed_closed))+' shutters with status: NRS_MSOP_VIG failed closed permanent\n\n')
  
  # Reassign shutter values as 'failed closed/open' where applicable, in the following order:
  
  msa_chk_failed_closed = list(msa_chk_failed_closed_perm) + list(msa_chk_failed_closed_inter) + list(msa_chk_failed_closed_vign)
  msa_chk_failed_open = list(msa_chk_failed_open_perm) + list(msa_chk_failed_open_inter)
  
  esa_map['Q'+str(k)][np.array(msa_chk_failed_closed)] = 0
  esa_map['Q'+str(k)][msa_shortmask_failed_closed] = 0
  esa_map['Q'+str(k)][np.array(msa_chk_failed_open)] = 1
  esa_map['Q'+str(k)][msa_vig_failed_closed] = 0
  


with open('esa_msa_map.dat', 'w') as file_writer:

  for k in range(1,5):
     for s in esa_map['Q'+str(k)]:
        file_writer.write('           '+str(s))
     file_writer.write('\n')         

print("Wrote shutter map values to esa_msa_map.dat ") 
file_writer.close()



quit()
