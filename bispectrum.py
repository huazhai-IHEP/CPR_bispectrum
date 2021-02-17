#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 16 20:54:08 2020

@author: zhaihua
"""

import numpy as np
#from camb_cl import CL
import calculator_bispectrum

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

import matplotlib as mpl
from matplotlib.ticker import ScalarFormatter

# usage: utils.Wigner3j( l2, l3, m1, m2, m3)

# bispectrum order: 0:TTB, 1:TTE, 2:TBB, 3:TEB,4:TEE,5:BBB, 6:EBB,7:EEB,8:EEE
import time

time_start=time.time()


# LCDM parameters
ombh2 = 0.02236
omch2=0.1186
tau = 0.034
ns =0.9674
As = 2.9943
H0 = 67.32
Alens = 1.0

r=0.01
nt=0
pivot_tensor=0.05


paramdic = {'ombh2':ombh2, 'omch2':omch2, 'tau':tau, 'ns':ns, 'As':As, 'H0':H0, 'Alens':Alens, \
            'r':r, 'nt':nt, 'pivot_tensor':pivot_tensor}

lmax_teb = 2500
lmax_start=5
lmax_alpha=1000
lmax_a = 2000
lmax_bis = 2000

lmax_computed_cl = 4000
alpha_bar =  1.

'''
# index 0 correspond to real zero multipole, 
CTE = np.zeros(lmax_teb)
CEE = np.zeros(lmax_teb)
cmb_cls = np.transpose(CL( paramdic, 'unlensed_total', lmax=lmax_teb, raw=0 ))
CTE = cmb_cls[3][:lmax_teb]
CEE = cmb_cls[1][:lmax_teb]
'''

loc = './input'
# no lensing
cmb_cls = np.transpose(np.loadtxt(loc+"/planck_2018_totCls.dat"))
DTE = cmb_cls[4][:lmax_teb]
DEE = cmb_cls[2][:lmax_teb]
DBB = cmb_cls[3][:lmax_teb]

CTE = np.zeros(lmax_teb)
CEE = np.zeros(lmax_teb)
CBB = np.zeros(lmax_teb)
for l in range(1, lmax_teb):
    fac1 = l*(l+1.)/2./np.pi
    CTE[l] = DTE[l-1]/fac1
    CEE[l] = DEE[l-1]/fac1
    CBB[l] = DBB[l-1]/fac1

ell = np.arange(2, lmax_a+1)
llfac = 2*np.pi/(ell*(ell+1.))
## caa, cta, cea calculated by camb, in units of mk rad, or rad^2
loc = './input'
cmb_cls_1 = np.transpose(np.loadtxt(loc+"/cpt_power_1_scalCls.dat"))
Daa_1 = cmb_cls_1[4][:lmax_teb]
DTa_1 = cmb_cls_1[5][:lmax_teb]
DEa_1 = cmb_cls_1[6][:lmax_teb]
Caa_1=np.zeros(lmax_a)
CTa_1=np.zeros(lmax_a)
CEa_1=np.zeros(lmax_a)
Caa_1[1:lmax_a] = Daa_1[:lmax_a-1]*llfac
CTa_1[1:lmax_a] = DTa_1[:lmax_a-1]*llfac
CEa_1[1:lmax_a] = DEa_1[:lmax_a-1]*llfac

cmb_cls_01 = np.transpose(np.loadtxt(loc+"/cpt_power_01_scalCls.dat"))
Daa_01 = cmb_cls_01[4][:lmax_teb]
DTa_01 = cmb_cls_01[5][:lmax_teb]
DEa_01 = cmb_cls_01[6][:lmax_teb]
Caa_01=np.zeros(lmax_a)
CTa_01=np.zeros(lmax_a)
CEa_01=np.zeros(lmax_a)
Caa_01[1:lmax_a] = Daa_01[:lmax_a-1]*llfac
CTa_01[1:lmax_a] = DTa_01[:lmax_a-1]*llfac
CEa_01[1:lmax_a] = DEa_01[:lmax_a-1]*llfac

cmb_cls_001 = np.transpose(np.loadtxt(loc+"/cpt_power_001_scalCls.dat"))
Daa_001 = cmb_cls_001[4][:lmax_teb]
DTa_001 = cmb_cls_001[5][:lmax_teb]
DEa_001 = cmb_cls_001[6][:lmax_teb]
Caa_001=np.zeros(lmax_a)
CTa_001=np.zeros(lmax_a)
CEa_001=np.zeros(lmax_a)
Caa_001[1:lmax_a] = Daa_001[:lmax_a-1]*llfac
CTa_001[1:lmax_a] = DTa_001[:lmax_a-1]*llfac
CEa_001[1:lmax_a] = DEa_001[:lmax_a-1]*llfac
    


# set for temp bispectrum integration.

# CPR angle: isotropic and anisotropic parameters

#bispectrum in form of l-1,l,l+1 ,  l=even for parity even, l=odd  odd parity.
# order, TTB, TTE, TBB, TEB, TEE, BBB, EBB, EEB, EEE

#bispectra_cpt_even = calculator_bispectrum.calculator_bispectrum('T ', 'even', alpha_bar, Caa_01, CEa_01, CTa_01, CEE, CTE,lmax_alpha, lmax_start, lmax_bis, lmax_computed_cl)
bispectra_cpt_odd = calculator_bispectrum.calculator_bispectrum('T ', 'odd', alpha_bar, Caa_1, CEa_1, CTa_1, CEE, CTE, lmax_alpha, lmax_start, lmax_bis, lmax_computed_cl)


#------------------------------------------------------------------------------------------------------------------------------


loc = './input'
fnl = 30.
bispectra_fnl_even = np.transpose(np.loadtxt(loc+"/planck2018_bispectrum_fnl_base_4_delta_4.dat"))
fnllabels = ['TTT', 'TTE', 'TET', 'TEE', 'ETT', 'ETE', 'EET', 'EEE']

bispectra_len_odd = np.transpose(np.loadtxt(loc+"/planck2018_bispectrum_lensing_base_4_delta_3.dat"))
lenlabels= ['TTB', 'TEB', 'TBT', 'TBE', 'ETB', 'EEB', 'EBT', 'EBE', 'BTT', 'BTE', 'BET', 'BEE']
loc_store = './output'

'''
bis_arr = np.transpose(np.loadtxt("./cpt_bispectrum_2000even.dat"))
bispectra_cpt_even = bis_arr[1:]

bis_arr = np.transpose(np.loadtxt("./cpt_bispectrum_2000odd.dat"))
bispectra_cpt_odd = bis_arr[1:]

'''

# plot l1(l1+1)l3(l3+1)bispectrum(l1,l2,l3)/(2pi)^2, or in here case (l2-1)l2(l2+1)(l2+2)..

ell = np.arange(1,lmax_bis+1)

new_fontsize = 22
ols, ils, als  = ['solid', 'solid', 'dashed']
olw, ilw, alw = [2.5, 2, 2]
sub_ratio = [1, 2, 2]
sub1, sub2 =[1, 0] 
zs, ts=50, 50

blabels = ['TTB', 'TTE', 'TBB', 'TEB','TEE','BBB','EBB','EEB','EEE']
lcs = ['blue', 'red', 'green', 'yellow',  'cyan', '#BF3EFF','#CDAD00', '#1E90FF', '#8B2323']
#lcs=['b']*9
lts = [(0, (5, 1, 5, 1)), (0, (7,1,7,1)), (0, (3, 3)), (0, (5, 5)), (0, (7, 7)), \
       (0, (1,3,1,3)), (0, (2, 5, 2, 5)),(0, (3,7,3,7)), (0, (4,9,4,9))]
       #(0, (5, 1, 5, 1)),(0, (7,1,7,1)),(0, (1,3,1,3)), (0, (3,7,3,7))]
#    (0, (1, 2, 2, 5)), (0, (1, 5, 3, 1)), (0, (1,3, 1, 3)) , (0, (1,2, 1, 5)) ]
 

plt.figure(figsize=(20, 10))
plt.ylabel('$\ell_1(\ell_1+1)\ell_3(\ell_3+1)b_{\ell_1\ell_2\ell_3}/(2\pi)^2 [\mu \mathrm{ K}^3]$')
plt.title( 'odd')
plt.axhline(y=0, c= 'grey', ls= 'dashed')
##plt.plot(ell_in, DBB1, label='orig-len, r=0.1')
plt.xlim(6, lmax_bis)
plt.yscale('log')
plt.xscale('log')
if (alpha_bar!=0):
    for i in range(5):
        plt.plot( ell, np.abs(ell**4*bispectra_cpt_odd[i]),  c=lcs[i], label=blabels[i], lw=2, ls = lts[i])
plt.legend() 
plt.savefig('./output/test.png') 
'''       
plt.figure(figsize=(20, 10))
plt.ylabel('$\ell_1(\ell_1+1)\ell_3(\ell_3+1)b_{\ell_1\ell_2\ell_3}/(2\pi)^2 [\mu \mathrm{ K}^3]$')
plt.title( 'even')
plt.axhline(y=0, c= 'grey', ls= 'dashed')
##plt.plot(ell_in, DBB1, label='orig-len, r=0.1')
plt.xlim(6, lmax_bis)
if (alpha_bar!=0):
    for i in range(9):
        plt.plot( ell, np.abs(ell**4*bispectra_cpt_even[i]),  c=lcs[i], label=blabels[i], lw=2, ls = lts[i])
        # primordial compare: fnl=, TTT, TTE, TEE
        if(i==0 or i==7):
            plt.plot( bispectra_fnl_even[0], fnl*np.abs(bispectra_fnl_even[0]**4*bispectra_fnl_even[i+1]),  c='grey', label=fnllabels[i], lw=1)
else:
    plt.plot( bispectra_fnl_even[0], fnl*np.abs(bispectra_fnl_even[0]**4*bispectra_fnl_even[1]),  c='pink', label='TTT', lw=1)
    plt.plot( bispectra_fnl_even[0], fnl*np.abs(bispectra_fnl_even[0]**4*bispectra_fnl_even[8]),  c='gold', label='EEE', lw=1)
    for i in range(9):
        #only TTB, TEB AND POLARIZATION
        if (i!=1 and i!=2 and i!=4):
            plt.plot( ell, np.abs(ell**4*bispectra_cpt_even[i]),  c=lcs[i], label=blabels[i], lw=3, ls = lts[i])

#plt.plot(ell[1:lmax_bis:2], np.abs(out_bispectrum[k, 1:lmax_bis:2]),  c='green')
plt.tick_params(which='both', direction='in', width = 1)
plt.tick_params(labelsize = new_fontsize)
#plt.set_xlabel(r'$\ell$', fontsize= new_fontsize)
plt.xscale('log')
plt.yscale('log')
plt.legend()
#plt.ylim(1e-25, 1e4)
plt.xlim(4, lmax_bis)  
#plt.set_ylabel('$l_1(l_1+1)l_3(l_3+1)b_{l_1l_2l_3}/(2\pi)^2 $')
plt.legend(loc=[0.03, 0.2])

#plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
#mpl.ticker.ScalarFormatter(useMathText = True)
#plt.yaxis.offsetText.set_fontsize(new_fontsize)

plt.tick_params(labelsize = new_fontsize)
plt.savefig(loc_store+'/bispectrum_even_fnl_%dx%dx%d.png'%(lmax_teb, lmax_a, lmax_bis), dpi=400)
#plt.show()



plt.figure(figsize=(20, 10))

plt.ylabel('$\ell_1(\ell_1+1)\ell_3(\ell_3+1)b_{\ell_1\ell_2\ell_3}/(2\pi)^2 [\mu \mathrm{ K}^3]$')

##plt.plot(ell_in, DBB1, label='orig-len, r=0.1')
plt.title( 'odd')
plt.axhline(y=0, c= 'grey', ls= 'dashed')

plt.xlim(6, lmax_bis)
if (alpha_bar!=0):
    for i in range(9):
        #plt.scatter(ell[2:lmax_bis-1:2], np.abs(out_bispectrum[i, 2:lmax_bis-1:2]), \
         #           c=lcs[i], marker='s', s=ts, zorder=zs, label=blabels[i] )
        plt.plot(ell, np.abs(ell**4*bispectra_cpt_odd[i]),  c=lcs[i], label=blabels[i], lw=3, ls = lts[i])
        # lensing odd :TTB, TEB, EEB
        if(i<2):
            plt.plot(bispectra_len_odd[0], np.abs(bispectra_len_odd[0]**4*bispectra_len_odd[i+1]),  c='grey', label=lenlabels[i], lw=1)
    #plt.plot(ell[0:lmax_bis:2], np.abs(out_bispectrum[k, 0:lmax_bis:2]),  c='r')
else:
    plt.plot(bispectra_len_odd[0], np.abs(bispectra_len_odd[0]**4*bispectra_len_odd[1]),  c='pink', label='TTB', lw=1)
    plt.plot(bispectra_len_odd[0], np.abs(bispectra_len_odd[0]**4*bispectra_len_odd[2]),  c='gold', label='TEB', lw=1)
    for i in range(9):
        # only TTE, TBB, TEE and polarization
        if(i!=0 and i!=3):
            plt.plot(ell, np.abs(ell**4*bispectra_cpt_odd[i]),  c=lcs[i], label=blabels[i], lw=3, ls = lts[i])
plt.tick_params(which='both', direction='in', width = 1)
plt.tick_params(labelsize = new_fontsize)
#plt.set_xlabel(r'$\ell$', fontsize= new_fontsize)
#plt.set_ylabel('$l_1(l_1+1)l_3(l_3+1)b_{l_1l_2l_3}/(2\pi)^2 $')
plt.xscale('log')
plt.yscale('log')
#plt.ylim(1e-23, 1e5)
plt.xlim(3, lmax_bis)

#ax = plt.subplot(111)
#box = ax.get_position()
#ax.set_position([box.x0, box.y0, box.width, box.height])
#legend_x = 1
#legend_y = 1.
#plt.legend(loc='best', bbox_to_anchor=(legend_x, legend_y))
plt.legend(loc=[0.03, 0.2])
#plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
#mpl.ticker.ScalarFormatter(useMathText = True)
#plt.yaxis.offsetText.set_fontsize(new_fontsize)

plt.tick_params(labelsize = new_fontsize)
plt.savefig(loc_store+'/bispectrum_odd_len_%dx%dx%d.png'%(lmax_teb, lmax_a, lmax_bis), dpi=400)
#plt.show()
'''
