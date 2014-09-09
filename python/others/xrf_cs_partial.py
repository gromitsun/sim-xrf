import sys, os
sys.path.append(os.path.realpath('./core/'))
from sim_xraylib import *
import matplotlib.pyplot as plt

# # # 2D plot of partial XRF cross sections as a function of E and Z

Z = np.arange(1,99)
ev0 = np.arange(1e3,41e3,1e3) # in ev

im_k =  np.array([mac_pe_line_kissel_cascade(ev0.reshape(1,-1), Z.reshape(-1,1), line) for line in K_LINES]).sum(axis=0)
im_l =  np.array([mac_pe_line_kissel_cascade(ev0.reshape(1,-1), Z.reshape(-1,1), line) for line in L_LINES]).sum(axis=0)
im_m =  np.array([mac_pe_line_kissel_cascade(ev0.reshape(1,-1), Z.reshape(-1,1), line) for line in M_LINES]).sum(axis=0)
im_n =  np.array([mac_pe_line_kissel_cascade(ev0.reshape(1,-1), Z.reshape(-1,1), line) for line in N_LINES]).sum(axis=0)

# K lines
plt.figure('2d_k')
plt.title('X-ray fluorescence production cross section (with full cascade) \n K-lines (cm$^2$/g)')
plt.imshow(im_k, extent=[1,40,1,99],aspect = 'auto')
plt.colorbar()
plt.xlabel('Energy (KeV)')
plt.ylabel('Atomic number (Z)')

# L lines
plt.figure('2d_l')
plt.title('X-ray fluorescence production cross section (with full cascade) \n L-lines (cm$^2$/g)')
plt.imshow(im_l, extent=[1,40,1,99],aspect = 'auto')
plt.colorbar()
plt.xlabel('Energy (KeV)')
plt.ylabel('Atomic number (Z)')

# M lines
plt.figure('2d_m')
plt.title('X-ray fluorescence production cross section (with full cascade) \n M-lines (cm$^2$/g)')
plt.imshow(im_m, extent=[1,40,1,99],aspect = 'auto')
plt.colorbar()
plt.xlabel('Energy (KeV)')
plt.ylabel('Atomic number (Z)')

# N lines
plt.figure('2d_n')
plt.title('X-ray fluorescence production cross section (with full cascade) \n N-lines (cm$^2$/g)')
plt.imshow(im_n, extent=[1,40,1,99],aspect = 'auto')
plt.colorbar()
plt.xlabel('Energy (KeV)')
plt.ylabel('Atomic number (Z)')


# # # 1D plot of partial XRF cross sections as a function of E
Z = symbol2number('Zn')
ev0 = np.arange(1e3,41e3,1e3) # in ev


cs_k = np.array([mac_pe_line_kissel_cascade(ev0, Z, line) for line in K_LINES]).sum(axis=0)
cs_l = np.array([mac_pe_line_kissel_cascade(ev0, Z, line) for line in L_LINES]).sum(axis=0)
cs_m = np.array([mac_pe_line_kissel_cascade(ev0, Z, line) for line in M_LINES]).sum(axis=0)
cs_n = np.array([mac_pe_line_kissel_cascade(ev0, Z, line) for line in N_LINES]).sum(axis=0)

plt.figure('1d_e')
plt.title('X-ray fluorescence production cross section (with full cascade)')
plt.plot(ev0/1.e3,cs_k,label = 'K lines')
plt.plot(ev0/1.e3,cs_l,label = 'L lines')
plt.plot(ev0/1.e3,cs_m,label = 'M lines')
plt.plot(ev0/1.e3,cs_n,label = 'N lines')
plt.xlabel('Energy (KeV)')
plt.ylabel('Cross section (cm$^2$/g)')

plt.legend()



# # # 1D plot of partial XRF cross sections as a function of Z
Z = np.arange(1,99)
ev0 = 40e3 # in ev


cs_k = np.array([mac_pe_line_kissel_cascade(ev0, Z, line) for line in K_LINES]).sum(axis=0)
cs_l = np.array([mac_pe_line_kissel_cascade(ev0, Z, line) for line in L_LINES]).sum(axis=0)
cs_m = np.array([mac_pe_line_kissel_cascade(ev0, Z, line) for line in M_LINES]).sum(axis=0)
cs_n = np.array([mac_pe_line_kissel_cascade(ev0, Z, line) for line in N_LINES]).sum(axis=0)

plt.figure('1d_z')
plt.title('X-ray fluorescence production cross section (with full cascade)')
plt.plot(Z,cs_k,label = 'K lines')
plt.plot(Z,cs_l,label = 'L lines')
plt.plot(Z,cs_m,label = 'M lines')
plt.plot(Z,cs_n,label = 'N lines')
plt.xlabel('Atomic number (Z)')
plt.ylabel('Cross section (cm$^2$/g)')

plt.legend()

plt.show()