#%%
import sys
import os
from pathlib import Path
os.environ['MPLCONFIGDIR'] = "./tmp"
import numpy as np
import matplotlib
import matplotlib as mpl
from matplotlib import rc
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator
import pandas as pd
import cantera as ct

#%%
# 0. Initialize parameter ---------------------------------------
Re_j   = 4913
njet   = 2
P      = 10.0 * ct.one_atm
D_j    = 5.0E-4; A_j = 0.25 * np.pi * D_j * D_j
intv   = 3.0;
Lx     = 26 * D_j;
Ly     = (2*intv+2) * D_j;
Lz     = 12 * D_j;
A_c    = Ly * Lz
T_j    = 300.;
T_c    = 750.;
X_j    = {}; X_j["H2"] = 1.0; X_j["N2"] = 1 - X_j["H2"]
X_c    = {}; X_c["O2"] = 0.21; X_c["N2"] = 0.79
mech   = "Wen30/chem.yaml";
freq   = 1000

#%%
if (True):
  matplotlib.rcParams['mathtext.fontset'] = 'custom'
  matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
  matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
  matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
  matplotlib.rcParams['mathtext.fontset'] = 'stix'
  matplotlib.rcParams['font.family'] = 'STIXGeneral'
  Tc = 1750.
  dT = 5.
  Ts = np.linspace(300, 3000, 2700)
  thd = 0.5 * (1 + np.tanh((Ts - Tc) / dT))
  fig, ax = plt.subplots(figsize = (4.5, 4))
  ax.plot(Ts, thd, "b", linewidth = 2.5)
  ax.set_xlabel(r"$T~\mathrm{[K]}$", fontsize = 24, color = "black")
  ax.set_xlim([1400, 2100])
  ax.set_ylim([-0.0, 1.02])
  ax.set_ylabel(r"$\dot{\omega}_{\alpha}(T)$", fontsize = 24, color = "blue")
  ax.tick_params(axis='both', which='major', labelsize=22, )
  ax.tick_params(axis='both', which='minor', labelsize=22, )
  ax.tick_params(axis='y', colors = "blue")
  ax.grid()

  ax2 = ax.twinx()
  xder = (Ts - Tc) / dT
  deriv = 0.5 * (1 - np.tanh(xder)**2) / dT
  ax2.plot(Ts, deriv, "r--", linewidth = 2.5)
  ax2.set_ylabel(r"$\partial \dot{\omega}_{\alpha} / \partial T~\mathrm{[1/K]}$", fontsize = 24, color = "red")
  ax2.tick_params(axis='both', which='major', labelsize=22)
  ax2.tick_params(axis='both', which='minor', labelsize=22)
  ax2.set_ylim([0, 0.2])
  ax2.set_yticks([0, 0.1, 0.2])
  ax2.tick_params(axis='y', colors = "red")
  plt.savefig("./fT.png", bbox_inches="tight", dpi=400)

#%%
if (True):
  nb = 16
  ny = 32; d0 = Ly / ny;
  nz = (int(Lz / d0 / nb) + 0) * nb;
  Lz = nz * d0;
  nx = (int(Lx / d0 / nb) + 1) * nb;
  Lx = nx * d0;
  U_c   = 17.5
  A_c    = Ly * Lz
  print("nx=",str(nx).rjust(3), " , Lx=", "%10.6E"%Lx, " , Lx/Dj=", Lx/D_j, " , dx=", Lx/nx)
  print("ny=",str(ny).rjust(3), " , Ly=", "%10.6E"%Ly, " , Ly/Dj=", Ly/D_j, " , dy=", Ly/ny)
  print("nz=",str(nz).rjust(3), " , Lz=", "%10.6E"%Lz, " , Lz/Dj=", Lz/D_j, " , dz=", Lz/nz)

#%%
gas_j = ct.Solution(mech)
gas_c = ct.Solution(mech)
spn_out = gas_j.species_names

id_H2 = gas_j.species_index("H2")
id_O2 = gas_j.species_index("O2")

gas_j.TPX = T_j, P, X_j
gas_c.TPX = T_c, P, X_c

rho_j = gas_j.density
rho_c = gas_c.density

nu_j = gas_j.viscosity / gas_j.density
nu_c = gas_c.viscosity / gas_c.density

U_j = Re_j * nu_j / D_j

m_j = njet * A_j * U_j * rho_j
m_c = A_c * U_c * rho_c
J = (rho_j * U_j**2) / (rho_c * U_c**2)

gas_mix = ct.Solution(mech)
h_mix   = m_j*gas_j.enthalpy_mass + m_c*gas_c.enthalpy_mass
Y_mix   = (m_j/(m_j+m_c))*gas_j.Y + (m_c/(m_j+m_c))*gas_c.Y
gas_mix.HPY = h_mix, P, Y_mix
equiv = gas_mix.X[id_H2] / (2*gas_mix.X[id_O2])

# Calculate zst
Yst = 4. / (4 + 32)
YO_O2 = gas_c.Y[id_O2]
YF_O2 = gas_j.Y[id_O2]
YF_H2 = gas_j.Y[id_H2]
zst = Yst*YO_O2 / (YF_H2 - Yst*(YF_O2-YO_O2))
print("Stoic z:", zst)

print("Density:     {:6.3E}, {:6.3E}".format(rho_j, rho_c))
print("Viscosity:   {:6.3E}, {:6.3E}".format(nu_j, nu_c))
print("Velocity:    {:6.3E}, {:6.3E}".format(U_j, U_c))
print("Sound speed: {:6.3E}, {:6.3E}".format(gas_j.sound_speed, gas_c.sound_speed))
print("Mach number: {:6.3E}, {:6.3E}".format(U_j/gas_j.sound_speed, U_c/gas_c.sound_speed))
print("mdot:        {:6.3E}, {:6.3E}".format(m_j, m_c))
print("Rv:          {:6.3E}".format(U_j/U_c))
print("J:           {:6.3E}".format(J))
print("Gloabl eqv:  {:6.3E}".format(equiv))

print("Domain x:    {:6.3E} - {:6.3E}".format(0-1.5*D_j, Lx-1.5*D_j))
print("Domain y:    {:6.3E} - {:6.3E}".format(-Ly/2., Ly/2.))
print("Domain z:    {:6.3E} - {:6.3E}".format(0.0, Lz))
print("Jet y coord: {:6.3E} ".format(0.5*(intv+1)*D_j))
print("Flow through time: {:6.3E}".format(Lx/U_c))
print("Oscillation:       {:6.3E}".format(1/freq))

#%%
nzbin = 100; zmin = 0.0; zmax = 1.0
zbins = np.linspace(zmin, zmax, nzbin)
gas_s = ct.Solution(mech)
res0D = []
ncol = 4 + len(spn_out)
nrow = nzbin
res_pmf = np.zeros((nrow, ncol))
for iz in range(0, nzbin):
  zs = zbins[iz]
  Ys = zs*gas_j.Y + (1-zs)*gas_c.Y
  Hs = zs*gas_j.enthalpy_mass + (1-zs)*gas_c.enthalpy_mass
  gas_s.HPY = Hs, P, Ys
  gas_s.equilibrate("HP")
  res0D.append(gas_s.T)

  res_pmf[iz, 0] = zbins[iz]
  res_pmf[iz, 1] = gas_s.T
  res_pmf[iz, 2] = 0.0
  res_pmf[iz, 3] = gas_s.density
  for isp, spn in enumerate(spn_out):
    igas = gas_s.species_index(spn)
    res_pmf[iz, 4+isp] = gas_s.Y[igas]
  if True:
    print(iz, np.amin(res_pmf[iz,4:]), np.sum(res_pmf[iz,4:]))

s = 'VARIABLES = "X" "temp" "u" "rho"'
for isp, spn in enumerate(spn_out):
        s = s + ' "' + spn + '"'
s = s + "\n"
s = s + "ZONE I=" + str(nrow) + " FORMAT=POINT SPECFORMAT=MASS"
np.savetxt("./initeq.dat", res_pmf, fmt='%.18e', header=s, comments="", delimiter="     ")

Tz = np.array(res0D)
#Tz[6:-4] = Tz[6]
fig, ax4 = plt.subplots()
ax4.plot(zbins, res0D)
ax4.set_ylabel("Tad")

#%%
Nx     = 128
Ny     = 56
Nz     = 64
jet_xs = [0.0, 0.0]
jet_ys = [-1.4E-3, 1.4E3]
xmin = -1.5*D_j; xmax = xmin + Lx
ymin = -Ly / 2.; ymax = ymin + Ly
zmin = 0.;       zmax = zmin + Lz
H = Lz
x_ = np.linspace(xmin, xmax, Nx)
y_ = np.linspace(ymin, ymax, Ny)
z_ = np.linspace(zmin, zmax, Nz)
xx, yy, zz = np.meshgrid(x_, y_, z_, indexing='ij')
mf = np.zeros((Nx,Ny,Nz))
temp = np.zeros((Nx,Ny,Nz)) + T_c

ix0 = np.argmax(x_>jet_xs[0])
iy = np.argmax(y_>jet_ys[0])
phi2D = mf[:,iy,:]
figunit_x = 10
fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(2*figunit_x*(Lz/Lx),figunit_x))
ax = axs[0]
ax2 = axs[1]
# Jet centerplane trajectory
xc = x_[(ix0-1):]; xc[0] = jet_xs[0]
b = 0.002 * np.power(xc/H, 0.1)
zc = H * 0.35 * np.power(xc/H, 0.45) * np.exp(-b)
xw = np.argmax(x_>jet_xs[0]-D_j/0.5)
wp = H * 0.48 * np.power(xc/H, 0.38)
wn = H * 0.12 * np.power(xc/H, 0.25)

ax2.plot(xc, zc)
ax2.plot(xc, wp, "g--")
ax2.plot(xc, wn, "g-.")

theta_mix = m_j / (m_j + m_c)
theta_c1 = 1 - 0.2 * np.power(xc/H, 2)
theta_c2 = theta_mix + (1-theta_mix) * np.power(0.6 / (xc/H), 0.5)
theta_c = theta_c1
indx = (xc/H > 1.0);
theta_c[indx] = (theta_c2[indx] + theta_c1[indx])/2.

cn = 0.2 * np.power(xc/H, -0.6) * np.exp(4.0 * np.power(xc/H, 2))
rcn = 1 - np.exp(-cn)

cp = 0.5 * np.power(xc/H, 1.0) #4.
rcp = 1 - np.exp(-cp)

mf_j = 1.0
mf_c = 0.0
iy = 0
theta_norm = np.zeros_like(mf)
for ix in range(0, Nx):
  for iz in range(0, Nz):
    ic = ix - ix0 +1
    if ix < ix0:
      theta_norm[ix,iy,iz] = 0.0
      continue
    zdis = z_[iz] - zc[ic]
    if zdis > 0:
      bot = -np.log(2)*zdis**2
      bot = bot / (wp[ic]+1E-30)**2
      theta_norm[ix,iy,iz] = np.exp(bot)
      theta_p = theta_c[ic] * rcp[ic]
      mf[ix,iy,iz] = theta_p + theta_norm[ix,iy,iz] * (theta_c[ic]-theta_p)
    else:
      bot = -np.log(2)*zdis**2
      bot = bot / (wn[ic]+1E-30)**2
      theta_norm[ix,iy,iz] = np.exp(bot)
      theta_n = theta_c[ic] * rcn[ic]
      mf[ix,iy,iz] = theta_n + theta_norm[ix,iy,iz] * (theta_c[ic]-theta_n)

    mf[ix,iy,iz] = mf[ix,iy,iz] * 0.08

    nzbin = len(res0D)
    dZ = 1.0 / nzbin
    ibin = int(mf[ix,iy,iz]/dZ)
    ibin = np.amax((0,ibin))
    ibin = np.amin((nzbin-2,ibin))
    z0 = ibin*dZ; r0 = mf[ix,iy,iz] - z0
    temp[ix,iy,iz] = res0D[ibin]*(1-r0) + res0D[ibin+1]*r0
for iy in range(1, Ny):
  mf[:,iy,:] = mf[:,0,:]
  temp[:,iy,:] = temp[:,0,:]

phi2D = temp[:,Ny-1,:]
vmin=0.; vmax=0.25
vmin=400; vmax=2700
#im = ax.imshow(phi2D.transpose(), extent=[xmin, xmax, zmin, zmax], origin="lower",
#         vmin = vmin, vmax = vmax,
#          cmap=cm.viridis)
#ctrz = ax.contourf(mf[:,5,:].transpose(), extent=[xmin, xmax, zmin, zmax], origin="lower",
#          vmin=vmin, vmax=vmax, levels=np.linspace(vmin,vmax,20),
#          cmap=cm.viridis)
ctr = ax.imshow(phi2D.transpose(), extent=[xmin, xmax, zmin, zmax], origin="lower",
          vmin=vmin, vmax=vmax, #levels=np.linspace(vmin,vmax,20),
          cmap=cm.viridis, aspect="equal")
#ax.contour(ctrz, levels=[zst], colors='r')
ax3 = ax2.twinx()
ax3.plot(xc, theta_c, "orange")
ax3.plot(xc, rcn, color="m", linestyle="--")
ax3.plot(xc, rcp, color="m", linestyle="-.")
ax3.set_xlim([xmin,xmax])
ax3.set_ylim([0, 1])
#ax.set_ylim([zmin,zmax])

#%% Calculate 1D flame
if True:
  Zs = [0.0252]
  lfs = []
  for iz, z in enumerate(Zs):
    Y_in        = {}
    Y_in["H2"]  = gas_j.Y[id_H2] * z + gas_c.Y[id_H2] * (1-z)
    Y_in["O2"]  = gas_j.Y[id_O2] * z + gas_c.Y[id_O2] * (1-z)
    Y_in["N2"]  = 1.0 - Y_in["H2"] - Y_in["O2"]
    h_mix = gas_j.enthalpy_mass * z + gas_c.enthalpy_mass * (1-z)
    gas_zst = ct.Solution(mech)
    gas_zst.HPY = h_mix, P, Y_in
    print(gas_zst.T)

    f = ct.FreeFlame(gas_zst)
    loglevel = 0
    f.set_max_grid_points(domain=f.domains[1], npmax=10000)
    f.set_refine_criteria(ratio=3.0, slope=0.1, curve=0.1)
    f.transport_model = 'mixture-averaged'
    f.solve(loglevel, auto=True)
    f.transport_model = 'multicomponent'
    f.set_refine_criteria(ratio=2.0, slope=0.05, curve=0.05)
    f.solve(loglevel)  # don't use 'auto' on subsequent solves

    gradf = np.zeros_like(f.T)
    gradf[1:-1] = (f.T[2:]-f.T[0:-2]) / (f.flame.grid[2:]-f.flame.grid[0:-2])
    max_gradf = np.amax(gradf)
    lf = (f.T[-1]-f.T[0]) / max_gradf
    lfs.append(lf)
    print("Flame thickness: ", z, lf)
    print("Flame speed: ", z, f.velocity[0])
    print("nu:", nu_j)

  fig, ax = plt.subplots()
  ax.plot(Zs, lfs)
# %%
patm = [6, 7, 8, 9, 10, 11, 12]
nu = [1.831E-5, 1.569E-5, 1.373E-5, 1.221E-5, 1.099E-5, 9.987E-6, 9.155E-6]
lf = [5.405E-5, 4.443E-5, 3.757E-5, 3.245E-5, 2.850E-5, 2.535E-5, 2.284E-5]
sl = [6.985E+0, 6.810E+0, 6.640E+0, 6.481E+0, 6.329E+0, 6.172E+0, 6.031E+0]

labelsize=20
fig, ax = plt.subplots()
ax.plot(patm, nu, color="b")
ax.set_ylabel(r"$\nu_\mathrm{st} \; [m^2/s]$", fontsize=20)
ax.set_xlabel(r"$P\; [atm]$", fontsize=20)

# %%
