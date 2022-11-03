from astropy.io import fits as pyfits
import numpy as np

# Access the most current MSA.msa FITS table containing MSA and SLIT metrology parameters from the Instrument Model directory:

with open('model.conf', 'r') as model_config_file:
  model_config_file_lines = [line for line in model_config_file.readlines() if not str(line).startswith("#")]
  msa_fits_table = model_config_file_lines[0].strip()+'/Description/MSA.msa'
  
model_config_file.close()

# Read the MSA.msa FITS table into a Python FITS object:
fitsobj = pyfits.open(msa_fits_table)
#primaryheader = (fitsobj[0]).header

# Initiate Python dictionary variables that will store relevant geometrical parameters of the NIRSpec micro-shutter arrays.

Qxref={} # position of the reference shutter (1,1) along the x-axis in MSA-plane coordinates (in m)
Qyref={} # position of the reference shutter (1,1) along the y-axis in MSA-plane coordinates (in m)
Qrot={} # rotation angle around the reference shutter (1,1) in rad (positive anticlockwise)

Qno={} # list of 4 elements, each containing the list of micro-shutter indexes (no = i + (j-1) * nx) for each quadrant.
Qxc={} # list of 4 elements, each containing an array of the x-axis position of the center of each micro-shutter (in m)
Qyc={} # list of 4 elements, each containing an array of the y-axis position of the center of each micro-shutter (in m)

Qnx={} # number of micro-shutters along the x-axis for each quadrant
Qny={} # number of micro-shutters along the y-axis for each quadrant

Qpitchx = {} # list of 4 elements, each containing an array of the x-axis pitch size of the micro-shutters (single value,   
             # common to all the shutters of a quadrant; spectral direction; in m)
Qpitchy = {} # list of 4 elements, each of them containing an array of the x-axis pitch size of the micro-shutters (single               
             # value, common to all the shutters of a quadrant; spectral direction; in m)

Qxopen = {} # x-axis size of each micro-shutter (spectral direction; in m)
Qyopen = {} # x-axis size of each micro-shutter (spectral direction; in m)

# Read and write the geometrical parameter values to an ASCII file for the eMPT reference files to access:
with open('MSA.msa.ascii', 'w') as file_writer:

  for q in range(1,5):

    header = (fitsobj['Q'+str(q)]).header

    Qxref[q] = float(header['QUADXREF'])
    Qyref[q] = float(header['QUADYREF'])
    Qrot[q] = float(header['QUADROT'])

    data = (fitsobj['Q'+str(q)]).data

    Qno[q] = data.field('NO')
    Qxc[q] = np.copy(data.field('XC'))
    Qyc[q] = np.copy(data.field('YC')).astype(float)

    Qnx[q] = 365
    Qny[q] = 171
    Qxopen[q] = '0.078 !mm'
    Qyopen[q] = '0.178'  

    Qpitchx[q] = Qxc[q][1] - Qxc[q][0]
    Qpitchy[q] = Qyc[q][Qnx[q]] - Qyc[q][0]
  
    file_writer.write('*Q'+str(q)+'xref\n')
    file_writer.write(str(Qxref[q])+' !mm\n\n')
    file_writer.write('*Q'+str(q)+'yref\n')
    file_writer.write(str(Qyref[q])+'\n\n')
    file_writer.write('*Q'+str(q)+'rot\n')
    file_writer.write(str(Qrot[q])+' !rad\n\n')
    file_writer.write('*Q'+str(q)+'xpitch\n')
    file_writer.write(str(Qpitchx[q]*1000.)+' !mm\n\n')
    file_writer.write('*Q'+str(q)+'ypitch\n')
    file_writer.write(str(Qpitchy[q]*1000.)+'\n\n')
    file_writer.write('*Q'+str(q)+'xopen\n')
    file_writer.write(Qxopen[q]+'\n\n')
    file_writer.write('*Q'+str(q)+'yopen\n')
    file_writer.write(Qyopen[q]+'\n\n')

print("Wrote MSA.msa metrology parameters to MSA.msa.ascii") 

file_writer.close()
fitsobj.close()


quit()
