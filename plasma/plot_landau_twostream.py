from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import glob
import sys
from matplotlib.ticker import ScalarFormatter
from matplotlib import rc
from vtk import vtkXMLPImageDataReader
from vtk.util import numpy_support

rc('text', usetex=True)
rc('font', family='serif')

writepath = "landau_twostream.eps"

############
output_directory_path = "lin_landau_32x32_pfc"

path = output_directory_path + "/csv/*.csv"
data_sum = np.loadtxt(path[:-5]+"output_0_0_0.csv", delimiter=' ')
num_rows = data_sum.shape[0]
data_sum[:,1:] = 0.

files = glob.glob(path)
num_files = len(files)

i = 0
for fname in files:
  data = np.loadtxt(fname, delimiter=' ')
  num_rows = min(num_rows,data.shape[0])
  data_sum[:num_rows,1:] += data[:num_rows,1:]
  i += 1

time = data_sum[:num_rows,0]
electric_energy_pfc = data_sum[:num_rows,1]
############
output_directory_path = "lin_landau_32x32_pfc_deltaf"

path = output_directory_path + "/csv/*.csv"
data_sum = np.loadtxt(path[:-5]+"output_0_0_0.csv", delimiter=' ')
num_rows = data_sum.shape[0]
data_sum[:,1:] = 0.

files = glob.glob(path)
num_files = len(files)

i = 0
for fname in files:
  data = np.loadtxt(fname, delimiter=' ')
  num_rows = min(num_rows,data.shape[0])
  data_sum[:num_rows,1:] += data[:num_rows,1:]
  i += 1

electric_energy_deltaf = data_sum[:num_rows,1]
############
output_directory_path = "nonlin_landau_32x64_pfc"

path = output_directory_path + "/csv/*.csv"
data_sum = np.loadtxt(path[:-5]+"output_0_0_0.csv", delimiter=' ')
num_rows = data_sum.shape[0]
data_sum[:,1:] = 0.

files = glob.glob(path)
num_files = len(files)

i = 0
for fname in files:
  data = np.loadtxt(fname, delimiter=' ')
  num_rows = min(num_rows,data.shape[0])
  data_sum[:num_rows,1:] += data[:num_rows,1:]
  i += 1

time_nonlin = data_sum[:num_rows,0]
electric_energy_nonlin_pfc = data_sum[:num_rows,1]
############
output_directory_path = "nonlin_landau_32x64_pfc_deltaf"

path = output_directory_path + "/csv/*.csv"
data_sum = np.loadtxt(path[:-5]+"output_0_0_0.csv", delimiter=' ')
num_rows = data_sum.shape[0]
data_sum[:,1:] = 0.

files = glob.glob(path)
num_files = len(files)

i = 0
for fname in files:
  data = np.loadtxt(fname, delimiter=' ')
  num_rows = min(num_rows,data.shape[0])
  data_sum[:num_rows,1:] += data[:num_rows,1:]
  i += 1

electric_energy_nonlin_deltaf = data_sum[:num_rows,1]




fig, axs = plt.subplots(2, 2, figsize=(7.,4.5))

axs[1,0].set_xlabel(r'$t\ /\ \omega_{p}^{-1}$')
axs[1,1].set_xlabel(r'$x / \lambda_{D,0}$')

axs[0,0].set_ylabel(r'$\epsilon_0 \mathbf{E}^2 / 2\ /\ (k_B T_{0})$')
axs[1,0].set_ylabel(r'$\epsilon_0 \mathbf{E}^2 / 2\ /\ (k_B T_{0})$')

axs[0,1].set_ylabel(r'$v / v_{t,0}$')
axs[1,1].set_ylabel(r'$v / v_{t,0}$')


# Log plot electric energy
axs[0,0].semilogy(time[:], electric_energy_pfc[:], color="black")
axs[0,0].semilogy(time[:], electric_energy_deltaf[:], color="0.6", linestyle="--")
axs[0,0].semilogy(time[:], 0.0009*np.exp(-0.3066*time[:]), color="orange")
 
# Log plot electric energy
axs[1,0].semilogy(time_nonlin[:], electric_energy_nonlin_pfc[:], color="black")
axs[1,0].semilogy(time_nonlin[:], electric_energy_nonlin_deltaf[:], color="0.6", linestyle="--")






field_name = 'f_e'
# PFC
loadpath = 'two_stream_64x64_pfc/phase_space/output30.pvti'

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

im = axs[0,1].pcolor(x, y, z, cmap='magma', vmin=z_min, vmax=z_max, rasterized=True, shading='auto')
fig.colorbar(im, ax=axs[0,1], pad = 0.06, aspect=25)#, shrink = .7)


# deltaf
loadpath = 'two_stream_64x64_pfc_deltaf/phase_space/output30.pvti'

reader.SetFileName(loadpath)
reader.Update()
img = reader.GetOutput()

z = numpy_support.vtk_to_numpy(img.GetCellData().GetArray(field_name))
z = np.reshape(z,(dims[1],dims[0]))
z_min = z.min()
z_max = z.max()

im = axs[1,1].pcolor(x, y, z, cmap='magma', vmin=z_min, vmax=z_max, rasterized=True, shading='auto')
fig.colorbar(im, ax=axs[1,1], pad = 0.06, aspect=25)#, shrink = .7)



axs[0,0].legend([r'Standard PFC', r'Moment Fitting'], loc='upper center', prop={'size':7.5})
axs[1,0].legend([r'Standard PFC', r'Moment Fitting'], loc='upper center', prop={'size':7.5})
axs[0,0].set_title(r'(a) Linear Landau Damping')
axs[1,0].set_title(r'(b) Non-linear Landau Damping')
axs[0,1].set_title(r'(c) Two-Stream Standard PFC')
axs[1,1].set_title(r'(d) Two-Stream Moment Fitting')

plt.tight_layout()
plt.savefig(writepath, bbox_inches='tight', dpi=300)

