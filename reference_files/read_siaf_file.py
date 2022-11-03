from astropy.io import fits as pyfits
import numpy as np

# Access the most current SIAF FITS tables from the Instrument Model directory:

with open('model.conf', 'r') as model_config_file:
  model_config_file_lines = [line for line in model_config_file.readlines() if not str(line).startswith("#")]
model_config_file.close()

msa_siaf_fitsobj = pyfits.open(model_config_file_lines[4].strip()) # NRS_SIAF_1.0.0_20220609T1004.fits

data = msa_siaf_fitsobj[1].data  
 
RefXPOSKY_NRS_FULL_MSA = data["RefXPOSKY"][36]
RefYPOSKY_NRS_FULL_MSA = data["RefYPOSKY"][36]
AngleV3_NRS_FULL_MSA = data["AngleV3"][36]

with open('msa_ref_v2v3.ascii', 'w') as file_writer:
  file_writer.write('# NRS_MSA_FULL V2,V3 reference point\n')
  file_writer.write(' \n')
  file_writer.write('*V2Ref\n')
  file_writer.write(str(RefXPOSKY_NRS_FULL_MSA)+'\n')
  file_writer.write(' \n')
  file_writer.write('*V3Ref\n')
  file_writer.write(str(RefYPOSKY_NRS_FULL_MSA)+'\n')
  file_writer.write(' \n')
  file_writer.write('*AngV3\n')
  file_writer.write(str(AngleV3_NRS_FULL_MSA))


print("Wrote V2Ref,V3Ref MSA fiducial reference point coordinates to msa_ref_v2v3.ascii") 
file_writer.close()


quit()
