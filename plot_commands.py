#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  4 17:05:48 2021

@author: joris
"""
import matplotlib.pyplot as plt
import numpy as np
#from settings_suds import *

def plot_contact_freqs(matrix_to_plot,genomic_length):
    plt.imshow(matrix_to_plot[::-1], cmap=plt.cm.BuPu, interpolation='none', extent=[0,genomic_length,0,genomic_length],vmin=0,vmax=0.02)
    plt.xlabel('Genome position [Mb]', fontsize = 15)
    plt.ylabel('Genome position [Mb]', fontsize = 15)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    #for x in expt_locations_domains:
        #plt.axvline(x/1000,0,0.03,color='black')
        #plt.axvline(x/1000,color='black',linestyle='dashed',linewidth=0.5,alpha=0.5)
    cbar = plt.colorbar() 
    cbar.ax.get_yaxis().labelpad = 15
    cbar.set_ticks(np.arange(0,0.025,0.005))
    cbar.ax.tick_params(labelsize=12)
    cbar.ax.set_ylabel('Contact frequency', rotation=270,fontsize=15)
    #    plt.savefig('/media/joris/raid_data/Chromosome_Maxent/Figures_paper/Appendix/Contacts_lin_oritether.png', bbox_inches='tight',dpi=300)
    #plt.savefig('/media/joris/raid_data/Chromosome_Maxent/Contacts_compare_pres_all.png', bbox_inches='tight',dpi=600)
    plt.show()

def plot_correlations(matrix_to_plot,genomic_length):
    fig = plt.figure(figsize=(8, 6))
    plt.imshow(matrix_to_plot[::-1], cmap='bwr', interpolation='none', extent=[0,genomic_length,0,genomic_length],vmin=-1.5,vmax=1.5)
    plt.xlabel('Genomic position [Mb]', fontsize = 15)
    plt.ylabel('Genomic position [Mb]', fontsize = 15)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    cbar = plt.colorbar() 
    cbar.ax.get_yaxis().labelpad = 15
    cbar.ax.tick_params(labelsize=15)
    cbar.ax.set_ylabel('Correlation of long axis position', rotation=270,fontsize=15)
    #plt.savefig('/media/joris/raid_data/Chromosome_Maxent/Figures_paper/Corr_longax_wt_freep.png', bbox_inches='tight',dpi=1200)
    plt.show()
    plt.close(fig)



def plot_mean_sudsizes(sizes_old,sizes_new,genomic_length,n_bins):
    plt.plot([x*genomic_length/n_bins for x in range(n_bins)],sizes_old,color='C0')
    plt.plot([x*genomic_length/n_bins for x in range(n_bins)],sizes_new,color='C1')
    plt.xlabel('Genomic position [Mb]', fontsize = 15)
    plt.ylabel('Mean SuD size [Mb]', fontsize = 15)
    plt.xlim(0,genomic_length)
    #plt.ylim(2.2*lattice_spacing,3.65*lattice_spacing)
    plt.ylim(0,1.1*max([max(sizes_old),max(sizes_new)]))
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.show()

def plot_extensions(extensions_old,extensions_stds_lower_old,extensions_stds_upper_old,extensions_new,extensions_stds_lower_new,extensions_stds_upper_new,reference_extensions,genomic_length,n_bins,lin_length,lattice_spacing):
    plt.plot([x*genomic_length/n_bins for x in range(n_bins)][:lin_length+1],[x*lattice_spacing for x in extensions_old][:lin_length+1],color='C0')
    plt.fill_between([x*genomic_length/n_bins for x in range(n_bins)][:lin_length+1],[(x-y)*lattice_spacing for x,y in zip(extensions_old,extensions_stds_lower_old)][n_bins-lin_length-1:],[(x+y)*lattice_spacing for x,y in zip(extensions_old,extensions_stds_upper_old)][n_bins-lin_length-1:],color='C0',alpha=0.3)
    
    plt.plot([x*genomic_length/n_bins for x in range(n_bins)][n_bins-lin_length-1:],[x*lattice_spacing for x in extensions_old][n_bins-lin_length-1:],color='C0')
    plt.fill_between([x*genomic_length/n_bins for x in range(n_bins)][n_bins-lin_length-1:],[(x-y)*lattice_spacing for x,y in zip(extensions_old,extensions_stds_lower_old)][n_bins-lin_length-1:],[(x+y)*lattice_spacing for x,y in zip(extensions_old,extensions_stds_upper_old)][n_bins-lin_length-1:],color='C0',alpha=0.3)
    
    plt.plot([x*genomic_length/n_bins for x in range(n_bins)][lin_length:n_bins-lin_length],[x*lattice_spacing for x in extensions_old][lin_length:n_bins-lin_length],color='k')
    plt.fill_between([x*genomic_length/n_bins for x in range(n_bins)][lin_length:n_bins-lin_length],[(x-y)*lattice_spacing for x,y in zip(extensions_old,extensions_stds_lower_old)][lin_length:n_bins-lin_length],[(x+y)*lattice_spacing for x,y in zip(extensions_old,extensions_stds_upper_old)][lin_length:n_bins-lin_length],color='k',alpha=0.3)
    
    plt.plot([x*genomic_length/n_bins for x in range(n_bins)][:lin_length+2],[x*lattice_spacing for x in extensions_new][:lin_length+2],color='C1')
    plt.fill_between([x*genomic_length/n_bins for x in range(n_bins)][:lin_length+1],[(x-y)*lattice_spacing for x,y in zip(extensions_new,extensions_stds_lower_new)][:lin_length+1],[(x+y)*lattice_spacing for x,y in zip(extensions_new,extensions_stds_upper_new)][:lin_length+1],color='C1',alpha=0.3)
    
    plt.plot([x*genomic_length/n_bins for x in range(n_bins)][n_bins-lin_length-2:],[x*lattice_spacing for x in extensions_new][n_bins-lin_length-2:],color='C1')
    plt.fill_between([x*genomic_length/n_bins for x in range(n_bins)][n_bins-lin_length-1:],[(x-y)*lattice_spacing for x,y in zip(extensions_new,extensions_stds_lower_new)][n_bins-lin_length-1:],[(x+y)*lattice_spacing for x,y in zip(extensions_new,extensions_stds_upper_new)][n_bins-lin_length-1:],color='C1',alpha=0.3)
    
    plt.plot([x*genomic_length/n_bins for x in range(n_bins)],[x*lattice_spacing for x in reference_extensions],linestyle=':',color='k')
    
    plt.xlabel('Genomic position [Mb]', fontsize = 15)
    plt.ylabel('Local extension [$\mu$m]', fontsize = 15)
    plt.xlim(0,genomic_length)
    #plt.ylim(2.2*lattice_spacing,3.65*lattice_spacing)
    plt.ylim(1*lattice_spacing,5.7*lattice_spacing)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.show()

def plot_local_extension_ratio(extensions_old,extensions_new,reference_extensions,genomic_length,n_bins):
    plt.plot([x*genomic_length/n_bins for x in range(n_bins)],[x/y for x,y in zip(extensions_old,reference_extensions)])
    plt.plot([x*genomic_length/n_bins for x in range(n_bins)],[x/y for x,y in zip(extensions_new,reference_extensions)])
    plt.axhline(1,color='k',linestyle='--',alpha=0.5)
    plt.xlabel('Genomic position [Mb]', fontsize = 15)
    plt.ylabel('Local extension ratio', fontsize = 15)
    plt.xlim(0,genomic_length)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.show()

def calc_smooth(input_vector,smoothwidth):
    len_input = len(input_vector)
    smooth_vector = [0]*len_input
    for i in range(len_input):
        smooth_vector[i] = np.average([input_vector[j%len_input] for j in range(i-smoothwidth,i+smoothwidth+1)])
    return smooth_vector

def plot_smooth_extension_ratio(extensions_old,extensions_new,reference_extensions,genomic_length,n_bins):
    smooth_old_frac = calc_smooth([x/y for x,y in zip(extensions_old,reference_extensions)],15)
    smooth_new_frac = calc_smooth([x/y for x,y in zip(extensions_new,reference_extensions)],15)
    plt.plot([x*genomic_length/n_bins for x in range(n_bins)],smooth_old_frac)
    plt.plot([x*genomic_length/n_bins for x in range(n_bins)],smooth_new_frac)
    plt.axhline(1,color='k',linestyle='--',alpha=0.5)
    plt.xlabel('Genomic position [Mb]', fontsize = 15)
    plt.ylabel('Local extension ratio', fontsize = 15)
    plt.xlim(0,genomic_length)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.show()

def plot_distr_oldpol(distr_oldpol,oriented_oldpol,min_l,max_l,genomic_length,n_bins):
    plt.imshow(np.flip(np.transpose(distr_oldpol),0),aspect='auto', cmap ='Blues',extent=[0,genomic_length,0,1])
    plt.plot([x*genomic_length/n_bins for x in range(n_bins)],[(x+min_l)/(min_l+max_l) for x in oriented_oldpol],'C0')
    plt.xlabel('Genomic position [Mb]', fontsize = 15)
    plt.ylabel('Scaled long-axis position', fontsize = 15)
    plt.ylim(0,1)
    plt.xlim(0,genomic_length)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    #plt.savefig('/media/joris/raid_data/Chromosome_Maxent/4dchrom_crescentus/Figures/localization_separations_old_'+ str(stage)+'.png', bbox_inches='tight',dpi=300)
    plt.show()

def plot_distr_newpol(distr_newpol,oriented_newpol,min_l,max_l,genomic_length,n_bins):
    plt.imshow(np.flip(np.transpose(distr_newpol),0),aspect='auto', cmap ='Oranges',extent=[0,genomic_length,0,1])
    plt.plot([x*genomic_length/n_bins for x in range(n_bins)],[(x+min_l)/(min_l+max_l) for x in oriented_newpol],'C1')
    plt.xlabel('Genomic position [Mb]', fontsize = 15)
    plt.ylabel('Scaled long-axis position', fontsize = 15)
    plt.ylim(0,1)
    plt.xlim(0,genomic_length)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    #plt.savefig('/media/joris/raid_data/Chromosome_Maxent/4dchrom_crescentus/Figures/localization_separations_new_'+ str(stage)+'.png', bbox_inches='tight',dpi=300)
    plt.show()
    
def plot_means_oriented(oriented_oldpol,oriented_newpol,lin_length,min_l,max_l,genomic_length,n_bins,add_separations=False,loc_separations_newpol=[],loc_separations_oldpol=[]):
    plt.plot([x*genomic_length/n_bins for x in range(n_bins)][:lin_length],[(x+min_l)/(min_l+max_l)for x in oriented_oldpol][:lin_length],linestyle='--',color='C0')
    plt.plot([x*genomic_length/n_bins for x in range(n_bins)][n_bins-lin_length:],[(x+min_l)/(min_l+max_l)for x in oriented_oldpol][n_bins-lin_length:],linestyle='--',color='C0')
    
    plt.plot([x*genomic_length/n_bins for x in range(n_bins)][:lin_length],[(x+min_l)/(min_l+max_l) for x in oriented_newpol][:lin_length],linestyle='--',color='C1')
    plt.plot([x*genomic_length/n_bins for x in range(n_bins)][n_bins-lin_length:],[(x+min_l)/(min_l+max_l) for x in oriented_newpol][n_bins-lin_length:],linestyle='--',color='C1')
    
    plt.plot([x*genomic_length/n_bins for x in range(n_bins)][lin_length:n_bins-lin_length],[(x+min_l)/(min_l+max_l) for x in oriented_newpol][lin_length:n_bins-lin_length],linestyle='--',color='k')
    
    if add_separations==True:
        plt.plot([x*genomic_length/n_bins for x in range(n_bins)],loc_separations_newpol,color='C1')
        plt.plot([x*genomic_length/n_bins for x in range(n_bins)],loc_separations_oldpol,color='C0')
        plt.plot([x*genomic_length/n_bins for x in range(n_bins)][lin_length:n_bins-lin_length],loc_separations_oldpol[lin_length:n_bins-lin_length],color='k')
    plt.xlabel('Genomic position [Mb]', fontsize = 15)
    plt.ylabel('Scaled long-axis position', fontsize = 15)
    plt.ylim(0,1)
    plt.xlim(0,genomic_length)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    #plt.savefig('/media/joris/raid_data/Chromosome_Maxent/4dchrom_crescentus/Figures/localization_comparison_old_new_'+ str(stage)+'.png', bbox_inches='tight',dpi=300)
    plt.show()

def plot_means_nearfar(means_nearpol,means_farpol,lin_length,min_l,max_l,genomic_length,n_bins,add_separations=False,loc_separations_near=[],loc_separations_far=[]):
    plt.plot([x*genomic_length/n_bins for x in range(n_bins)][:lin_length],[(x+min_l)/(min_l+max_l) for x in means_nearpol][:lin_length],linestyle='--',color='C0')
    plt.plot([x*genomic_length/n_bins for x in range(n_bins)][n_bins-lin_length:],[(x+min_l)/(min_l+max_l) for x in means_nearpol][n_bins-lin_length:],linestyle='--',color='C0')
    
    plt.plot([x*genomic_length/n_bins for x in range(n_bins)][:lin_length],[(x+min_l)/(min_l+max_l) for x in means_farpol][:lin_length],linestyle='--',color='C1')
    plt.plot([x*genomic_length/n_bins for x in range(n_bins)][n_bins-lin_length:],[(x+min_l)/(min_l+max_l) for x in means_farpol][n_bins-lin_length:],linestyle='--',color='C1')
    
    plt.plot([x*genomic_length/n_bins for x in range(n_bins)][lin_length:n_bins-lin_length],[(x+min_l)/(min_l+max_l) for x in means_farpol][lin_length:n_bins-lin_length],linestyle='--',color='k')
    
    if add_separations==True:
        plt.plot([x*genomic_length/n_bins for x in range(n_bins)],loc_separations_near,color='C0')
        plt.plot([x*genomic_length/n_bins for x in range(n_bins)],loc_separations_far,color='C1')
        plt.plot([x*genomic_length/n_bins for x in range(n_bins)][lin_length:n_bins-lin_length],loc_separations_far[lin_length:n_bins-lin_length],color='k')
    
    plt.xlabel('Genomic position [Mb]', fontsize = 15)
    plt.ylabel('Scaled long-axis position', fontsize = 15)
    plt.ylim(0,1)
    plt.xlim(0,genomic_length)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    #plt.savefig('/media/joris/raid_data/Chromosome_Maxent/4dchrom_crescentus/Figures/localization_comparison_near_far_'+ str(stage)+'.png', bbox_inches='tight',dpi=300)
    plt.show()

def plot_distr_near(distr_nearpol,means_nearpol,min_l,max_l,genomic_length,n_bins):
    plt.imshow(np.flip(np.transpose(distr_nearpol),0),aspect='auto', cmap ='Blues',extent=[0,genomic_length,0,1])
    plt.plot([x*genomic_length/n_bins for x in range(n_bins)],[(x+min_l)/(min_l+max_l) for x in means_nearpol],'C0')
    plt.xlabel('Genomic position [Mb]', fontsize = 15)
    plt.ylabel('Scaled long-axis position', fontsize = 15)
    plt.ylim(0,1)
    plt.xlim(0,genomic_length)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    #plt.savefig('/media/joris/raid_data/Chromosome_Maxent/4dchrom_crescentus/Figures/localization_separations_near_'+ str(stage)+'.png', bbox_inches='tight',dpi=300)
    plt.show()

def plot_distr_far(distr_farpol,means_farpol,min_l,max_l,genomic_length,n_bins):
    plt.imshow(np.flip(np.transpose(distr_farpol),0),aspect='auto', cmap ='Oranges',extent=[0,genomic_length,0,1])
    plt.plot([x*genomic_length/n_bins for x in range(n_bins)],[(x+min_l)/(min_l+max_l) for x in means_farpol],'C1')
    plt.xlabel('Genomic position [Mb]', fontsize = 15)
    plt.ylabel('Scaled long-axis position', fontsize = 15)
    plt.ylim(0,1)
    plt.xlim(0,genomic_length)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    #plt.savefig('/media/joris/raid_data/Chromosome_Maxent/4dchrom_crescentus/Figures/localization_separations_far_'+ str(stage)+'.png', bbox_inches='tight',dpi=300)
    plt.show()

