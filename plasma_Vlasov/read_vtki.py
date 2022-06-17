
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import glob
import sys
from matplotlib.ticker import ScalarFormatter
from matplotlib import rc
from vtk import vtkXMLPImageDataReader
from vtk.util import numpy_support

def Vtki_PhaseSpaceData(n):

	'''
	
	Return phase space data from vtki-files
	
	n : filenumber (n \in [0, 100])
	
	'''
	
	field_name = 'f_e'
	loadpath = f'two_stream_64x64_pfc/phase_space/output{n}.pvti'
	
	reader = vtkXMLPImageDataReader()
	reader.SetFileName(loadpath)
	reader.Update()
	img = reader.GetOutput()
	
	dims = [img.GetDimensions()[0]-1,img.GetDimensions()[1]-1,img.GetDimensions()[2]-1]
	dx = 4.*np.pi/dims[0]
	dv = 10./dims[1]
	
	x = np.array([i*dx for i in range(dims[0])])
	y = np.array([-5. + i*dv for i in range(dims[1])])
	z = numpy_support.vtk_to_numpy(img.GetCellData().GetArray(field_name))
	z = np.reshape(z,(dims[1],dims[0]))
	
	return x,y,z

def Vtki_PhaseSpacePlot(n):
	'''
	
	Plot phase space data in vtki
	
	n : filenumber (n \in [0, 100])
	
	'''
	
	field_name = 'f_e'
	loadpath = f'two_stream_64x64_pfc/phase_space/output{n}.pvti'
	
	#fig = plt.figure()
	
	reader = vtkXMLPImageDataReader()
	reader.SetFileName(loadpath)
	reader.Update()
	img = reader.GetOutput()
	
	dims = [img.GetDimensions()[0]-1,img.GetDimensions()[1]-1,img.GetDimensions()[2]-1]
	dx = 4.*np.pi/dims[0]
	dv = 10./dims[1]
	
	x = np.array([i*dx for i in range(dims[0])])
	y = np.array([-5. + i*dv for i in range(dims[1])])
	z = numpy_support.vtk_to_numpy(img.GetCellData().GetArray(field_name))
	z = np.reshape(z,(dims[1],dims[0]))
	
	z_min = z.min()
	z_max = z.max()
	
	im = plt.pcolor(x, y, z, cmap='magma', vmin=z_min, vmax=z_max, rasterized=True, shading='auto')
	plt.colorbar(im, pad = 0.06, aspect=25)#, shrink = .7)
	plt.title(n)

	plt.show()
	savepath = loadpath.replace('.pvti', '.png')
	#plt.savefig(savepath)
	
	return 0
	
#Vtki_PhaseSpacePlot(30)
x, y, fe = Vtki_PhaseSpaceData(30)

print('xshape : ', np.shape(x))
print('yshape : ', np.shape(y))
print('feshape : ', np.shape(fe))
