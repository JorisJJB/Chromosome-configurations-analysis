
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 23 01:25:39 2021

@author: joris
"""

import matplotlib.pyplot as plt
import numpy as np
import math
import plot_commands

compute_twopoint = False
add_separations = True
calculate_suds = True

save_twopoint = False
save_localizations = False
save_label_positions = False

#Replication stage to be calculated
stage = 5

#System properties
pol_length = 1620
reduction_factor = 4

n_bins = int(pol_length/reduction_factor)
ori_sim = int(1220/reduction_factor)
ori_new = 0
ter = int(n_bins/2)
#ori=0
ter_sim = int(ori_sim-n_bins/2)
genomic_length = 4.017
lattice_spacing = 0.088
n_configs = 10

averaging_domain_suds = 3
cutoff_suds = 3.0
extension_width = 2

#lin_lengths = [0,98,300,454,606,758]
#lin_lengths = [0,100,312,472,632,792]

#fork_positions_new
lin_lengths = [0,72,292,456,624,788]
#cell_lengths = [29,31,38,42,47,52]
cell_lengths = [23,25,29,32,36,40]

#cell_lengths_new
#cell_lengths = [25,26,29,31,35,38]


labelled_positions = [0,0.43,0.96,2.03,2.67,3.03]
label_binnrs = [int(round(n_bins*x/genomic_length)) for x in labelled_positions]
n_labels = len(label_binnrs)

#sim_files = ["","2021-08-23_013854_10","2021-08-23_013828_30","2021-08-23_013759_45","2021-08-23_013729_60","2021-08-23_013653_75"]
#sim_files = ["2021-09-12_235852","2021-09-13_000249","2021-09-13_000520","2021-09-13_000608","2021-09-13_000641"]
#sim_files = ["2021-10-17_233139_0","2021-10-17_233202_10","2021-10-17_233224_30","2021-10-17_233244_45","2021-10-17_233304_60","2021-10-17_233319_75"]
#sim_files = ["2021-10-18_124333_0","2021-10-18_124348_10","2021-10-18_124407_30","2021-10-18_124421_45","2021-10-18_124436_60","2021-10-18_124447_75"]
file_path = "/media/joris/raid_data/Chromosome_Maxent/4dchrom_crescentus/forward/"
sim_files = ["2021-10-20_091804_0","2021-10-20_091849_10","2021-10-20_091929_30","2021-10-20_092003_45","2021-10-20_092042_60","2021-10-23_144008_75"]

#grzes_path = "forward_runs/ori_only/"
#sim_files = ["2021-09-10_160439","2021-09-10_160522","2021-09-10_210520","2021-09-10_210606","2021-09-10_210705"]

#grzes_path = "forward_runs/HiC_only/"
#sim_files = ["2021-09-11_210557","2021-09-11_210639","2021-09-14_202127","2021-09-14_202204","2021-09-14_202316"]

#file_path = "/media/joris/raid_data/Chromosome_Maxent/4dchrom_crescentus/Grzes_sim/hic_only/"
#sim_files = ["2021-10-21_093011_0","2021-10-21_093107_10","2021-10-21_093207_30","2021-10-21_093307_45","2021-10-21_093405_60","2021-10-21_093503_75"]

#file_path = "/media/joris/raid_data/Chromosome_Maxent/4dchrom_crescentus/Grzes_sim/ori_only/"
#sim_files = ["2021-10-20_171249_0","2021-10-20_171321_10","2021-10-20_171404_30","2021-10-20_171506_45","2021-10-20_171603_60","2021-10-20_171706_75"]

#Factor 2 data
#sim_files = ["","2021-10-29_183659_10","2021-10-29_183724_30","2021-10-29_183748_45","2021-10-29_183808_60","2021-10-29_183830_75"]

lin_length = int(lin_lengths[stage]/reduction_factor)
cell_length = cell_lengths[stage]
min_l = int(math.floor((cell_length+cell_length%2-1)/2))
max_l = int(math.floor((cell_length)/2))
sim_file = sim_files[stage]

reference_extensions = np.loadtxt("/media/joris/raid_data/Chromosome_Maxent/Dataprocess_check_review/Raw_r0/Forward/stiffness_distributions.txt")
reference_extensions = reference_extensions[n_bins:2*n_bins]

reference_positions_file = '/media/joris/raid_data/Chromosome_Maxent/4dchrom_crescentus/Separations/Localization_data/'

def preprocess_data(data_long,ori):
    data = [[[] for x in range(2)] for y in range(n_configs)]
    for configindex, config in enumerate(data_long):
        for strandindex,strand in enumerate(config):
            for monindex, mon in enumerate(strand):
                if monindex%4==0:
                    data[configindex][strandindex].append(mon)
            data[configindex][strandindex] = np.roll(data[configindex][strandindex],-ori,axis=0)
          
    for configindex, config in enumerate(data):
        strand = config[1]
        for monindex, mon in enumerate(strand):
            if lin_length > n_bins:
                if((lin_length)%n_bins < monindex<(lin_length)%n_bins):
                    data[configindex][1][monindex]  = data[configindex][0][monindex] 
            else:
                if(((lin_length) < monindex) and (monindex<(lin_length)%n_bins)):
                    data[configindex][1][monindex] = data[configindex][0][monindex]
    return data

def dist(vector1,vector2):
    return(math.sqrt(sum([(x-y)**2 for x,y in zip(vector1,vector2)])))

def overlap(vector1,vector2):
    if dist(vector1,vector2)==0:
        return 1
    else:
        return 0

def compute_longax_corr(matrix,avs,sq_avs,avs2=0,sq_avs2=0):
    n_monomers = len(matrix)
    corr_matrix = [[0 for x in range(n_monomers)] for y in range(n_monomers)]
    if avs2==0:
        for mon1 in range(n_monomers):
            for mon2 in range(n_monomers):
                corr_matrix[mon1][mon2] = (matrix[mon1][mon2]-avs[mon1]*avs[mon2])/(math.sqrt(sq_avs[mon1]-avs[mon1]**2)*math.sqrt(sq_avs[mon2]-avs[mon2]**2))
    else:
        for mon1 in range(n_monomers):
            for mon2 in range(n_monomers):
                corr_matrix[mon1][mon2] = (matrix[mon1][mon2]-avs[mon1]*avs2[mon2])/(math.sqrt(sq_avs[mon1]-avs[mon1]**2)*math.sqrt(sq_avs2[mon2]-avs2[mon2]**2))
    return corr_matrix


def save_localization_data(data,file_name):
    file = open("/media/joris/raid_data/Chromosome_Maxent/4dchrom_crescentus/"+file_name+".txt","w") 
    for x in range(n_bins):
        file.write(str(data[x]))
        if x !=n_bins-1:
            file.write(' ')
    file.close()

def save_twopoint_data(data_matrix,file_name):
    file = open("/media/joris/raid_data/Chromosome_Maxent/4dchrom_crescentus/"+file_name+".txt","w") 
    for x in range(n_bins):
        for y in range(n_bins):
            file.write(str(data_matrix[x][y]))
    #        print(shifted_matrix[i][j])
            if y !=n_bins-1:
                file.write(' ')
        if x!=n_bins-1:
            file.write('\n')
    file.close()
    
class ReferenceLocalizations:
    def __init__(self):
        self.ref_means_near = [0]*n_bins
        self.ref_means_far = [0]*n_bins
        self.ref_means_oldpol = [0]*n_bins
        self.ref_means_newpol = [0]*n_bins
    def read_ref_loc(self,file):
        self.loc_separations_near = np.loadtxt(file+"oriented_near_"+str(stage)+'.txt', dtype=float)
        self.loc_separations_far = np.loadtxt(file+"oriented_far_"+str(stage)+'.txt', dtype=float)
        self.loc_separations_oldpol = np.loadtxt(file+"oriented_oldpol_"+str(stage)+'.txt', dtype=float)
        self.loc_separations_newpol = np.loadtxt(file+"oriented_newpol_"+str(stage)+'.txt', dtype=float)

class SudStats:
    global n_configs
    global averaging_domain_suds
    global cutoff_suds
    global n_bins
    def __init__(self):
        self.mean_domainsizes = [0]*n_bins
#        self.all_domainsizes = []
    
#    def dist(vector1,vector2):
#        return(math.sqrt(sum([(x-y)**2 for x,y in zip(vector1,vector2)])))
    
    def calc_sud_properties(self, strand):
        strand_len = len(strand)
        domainsizes = [0]*strand_len
        matrixofdists = [[0 for i in range(strand_len)] for j in range(strand_len)]
        for i in range(strand_len):
            for j in range(strand_len):
              matrixofdists[i][j] = dist(strand[i],strand[j]) 
        for i in range(strand_len):
            distset = averaging_domain_suds
            stop=0

            while stop==0:
                currav =0
                ndata = 0
                for loc_dist in range(distset-averaging_domain_suds,distset):
                    for r in [-loc_dist,loc_dist]:
                        currav += matrixofdists[i][(i+r)%strand_len]
                        ndata +=1
                currav  = currav/ndata

                if currav < cutoff_suds:
                    distset +=1
                else:
                    stop=1
                    domainsizes[i] = loc_dist-averaging_domain_suds
            distset = averaging_domain_suds
            stop=0
            while stop==0:
                currav =0
                ndata = 0
                for loc_dist in range(distset-averaging_domain_suds,distset):
                    for r in [-loc_dist,loc_dist]:
                        currav += matrixofdists[i][(i+r)%strand_len]
                        ndata +=1
                currav  = currav/ndata
                if currav < cutoff_suds:
                    distset +=1
                else:
                    stop=1
                    domainsizes[i] = loc_dist-averaging_domain_suds
                if loc_dist >=strand_len:
                    stop=1
        
        for i in range(strand_len):
            for j in range(strand_len):
                if domainsizes[i]<domainsizes[j]-min(abs(j-i),strand_len-max(i,j)+min(i,j)):
                    domainsizes[i] = domainsizes[j]-min(abs(j-i),strand_len-max(i,j)+min(i,j))
       
        areasize = 20
        peaks = []
        values = []
        for i in range(strand_len):
            countmax = 0
            for l in range(-areasize,areasize+1):
                if domainsizes[(i-l)%strand_len] <= domainsizes[i]:
                    countmax+=1
                if countmax==2*areasize+1:
                    peaks.append(i)
                    values.append(domainsizes[i])
        peaks_reduced = []
        x=0
        
        while x< len(values):
            compare = values[x]
            countn = 0
            while compare ==values[x]: #runs until last matching element is found. Last matching element is at x+countn-1
                if x +countn == len(values)-1:
                    countn += 1
                    break
                if peaks[x+countn+1] >peaks[x+countn]+areasize:
                    countn += 1
                    break
                countn +=1
                compare = values[x+countn]
                
            peaks_reduced.append(peaks[int(x+countn-1)])
            x += countn
         
#        for y in peaks_reduced:
#            self.all_domainsizes.append(domainsizes[y])
        
#        centers_longax = []
#        dists_longax = []
#        occ_left = [[0 for x in range(width)] for y in range(length)]
#        occ_left_allsuds = [[[0 for x in range(width)] for y in range(length)] for z in range(len(peaks_reduced))]
#        occ_left_allsuds_hires = [[[0 for x in range(width*resolution_factor)] for y in range(length*resolution_factor)] for z in range(len(peaks_reduced))]
#        occ_right = [[0 for x in range(width)] for y in range(length)]
#        occ_right_allsuds = [[[0 for x in range(width)] for y in range(length)] for z in range(len(peaks_reduced))]
#        occ_right_allsuds_hires = [[[0 for x in range(width*resolution_factor)] for y in range(length*resolution_factor)] for z in range(len(peaks_reduced))]
#        centers_left = []
#        centers_right = []
#        for y in range(len(peaks_reduced)):
#            #count long-axis positions for all sites that are part of cluster
#            longax_occ = [0]*length
#            for x in range(peaks_reduced[y]-domainsizes[peaks_reduced[y]],peaks_reduced[y]+domainsizes[peaks_reduced[y]]):
#                longax_occ[int(configs[m][reduction_factor*x%polymer_length][2])] +=1
#                if peaks_reduced[y]<polymer_length/(2*reduction_factor):
#                    occ_left[int(configs[m][reduction_factor*x%polymer_length][2])][int(configs[m][reduction_factor*x%polymer_length][1])] +=1
#                    occ_left_allsuds[y][int(configs[m][reduction_factor*x%polymer_length][2])][int(configs[m][reduction_factor*x%polymer_length][1])] +=1
#                else:
#                    occ_right[int(configs[m][reduction_factor*x%polymer_length][2])][int(configs[m][reduction_factor*x%polymer_length][1])] +=1
#                    occ_right_allsuds[y][int(configs[m][reduction_factor*x%polymer_length][2])][int(configs[m][reduction_factor*x%polymer_length][1])] +=1
#            if peaks_reduced[y]<polymer_length/(2*reduction_factor):
#                centers_left.append([int(configs[m][reduction_factor*peaks_reduced[y]%polymer_length][2]),int(configs[m][reduction_factor*peaks_reduced[y]%polymer_length][1])])
#            else:
#                centers_right.append([int(configs[m][reduction_factor*peaks_reduced[y]%polymer_length][2]),int(configs[m][reduction_factor*peaks_reduced[y]%polymer_length][1])])
#            dists_longax.append(longax_occ)
#            centers_longax.append(np.argmax(longax_occ))
#            
#            higher_resolution(occ_left_allsuds[y],occ_left_allsuds_hires[y])
#            occ_left_allsuds_hires[y] = gaussian_filter(occ_left_allsuds_hires[y], sigma,output='float')
#            
#            higher_resolution(occ_right_allsuds[y],occ_right_allsuds_hires[y])
#            occ_right_allsuds_hires[y] = gaussian_filter(occ_right_allsuds_hires[y], sigma,output='float')
#        
#        sums_dists_all.append(np.sum(domainsizes)/n_regions)
#        max_sizes.append(np.max(domainsizes))
#        
#        n_leftarm = 0
#        n_rightarm = 0
#        for y in peaks_reduced:
#            for i in range(y-5,y+6):
#                peak_positions[i%n_regions]+=1/n_snaps
#            if y<n_regions/2:
#                n_leftarm += 1
#            else:
#                n_rightarm +=1
#        hist_leftarm[n_leftarm] +=1/n_snaps
#        hist_rightarm[n_rightarm] += 1/n_snaps
#        
        peaksizes_config = [0]*strand_len
        in_cluster = [0]*strand_len
        for y in peaks_reduced:
            for xcoord in range((y-domainsizes[y]),(y+domainsizes[y])):
                x = xcoord%strand_len
                in_cluster[x] = 1
                if peaksizes_config[x]==0:
                    peaksizes_config[x] = domainsizes[y]
#                    ndata_peaksizes[x] += 1
                else:
                    if domainsizes[y]>peaksizes_config[x]:
                        peaksizes_config[x] = domainsizes[y]
                        
#        frac_in_domain += np.sum(in_cluster)/(strand_len*n_snaps) 
                  
        for i in range(strand_len):
            self.mean_domainsizes[i] += peaksizes_config[i]/n_configs
#            
#        for y in peaks_reduced:
#            if domainsizes[y]>min_csize:
#                for x in range(y-domainsizes[y],(y+domainsizes[y])%n_regions):
#                    dist_center_c[min(abs(y-x),abs(n_regions+max(y,x)-min(y,x)))] += matrixofdists[m][y][x]
#                    ndata_center_c[min(abs(y-x),abs(n_regions+max(y,x)-min(y,x)))] +=1
#                    for z in range(x,(y+domainsizes[y])%n_regions):
#                        dists_all_c[min(abs(z-x),abs(n_regions+max(z,x)-min(z,x)))] += matrixofdists[m][z][x]
#                        ndata_dist_all_c[min(abs(z-x),abs(n_regions+max(z,x)-min(z,x)))] += 1
#                
#                rand_site = np.random.randint(0,n_regions-1)
#                for x in range(rand_site-domainsizes[y],(rand_site+domainsizes[y])%n_regions):
#                    rand_dists[min(abs(rand_site-x),abs(n_regions+max(rand_site,x)-min(rand_site,x)))] += matrixofdists[m][rand_site][x]
#                    for z in range(x,(rand_site+domainsizes[y])%n_regions):
#                        dists_all_random[min(abs(z-x),abs(n_regions+max(z,x)-min(z,x)))] += matrixofdists[m][z][x]
#                for i in range(2*min_csize):
#                    for j in range(2*min_csize):
#                        r_corr_rand[i][j] += matrixofdists[m][(rand_site+i-min_csize)%n_regions][rand_site]*matrixofdists[m][(rand_site+j-min_csize)%n_regions][rand_site]
#                        if theta_defined(configs[m][reduction_factor*((rand_site+i-min_csize)%n_regions)],configs[m][reduction_factor*((rand_site+j-min_csize)%n_regions)], configs[m][reduction_factor*rand_site]) ==1:
#                            theta_corr_rand[i][j] += calc_theta(configs[m][reduction_factor*((rand_site+i-min_csize)%n_regions)],configs[m][reduction_factor*((rand_site+j-min_csize)%n_regions)], configs[m][reduction_factor*rand_site])
#                            ndata_theta_rand[i][j] +=1
#                        
#                        if phi_defined(configs[m][reduction_factor*((rand_site+i-min_csize)%n_regions)],configs[m][reduction_factor*((rand_site+j-min_csize)%n_regions)], configs[m][reduction_factor*rand_site]) ==1:
#                            phi_corr_rand[i][j] += calc_phi(configs[m][reduction_factor*((rand_site+i-min_csize)%n_regions)],configs[m][reduction_factor*((rand_site+j-min_csize)%n_regions)], configs[m][reduction_factor*rand_site])
#                            ndata_phi_rand[i][j] +=1
#                    r_av_rand[i] += matrixofdists[m][(rand_site+i-min_csize)%n_regions][rand_site]
#                    r_sq_av_rand[i] += matrixofdists[m][(rand_site+i-min_csize)%n_regions][rand_site]**2        
#                        
#        peaks_included_config = 0
#        for y in peaks_reduced:
#            if domainsizes[y]>min_csize:
#                peaks_included_config += 1
#                for i in range(2*min_csize):
#                    for j in range(2*min_csize):
#                        distances_cluster[i][j] += matrixofdists[m][(y+i-min_csize)%n_regions][(y+j-min_csize)%n_regions]
#                ncluster_distances += 1
#                
#                if (min_csize <y<(n_regions/2-min_csize)) or ((n_regions/2+min_csize)<y<(n_regions-min_csize)):
#                    if y<n_regions/2:
#                        for i in range(2*min_csize):
#                            longax_positions[i] += configs[m][reduction_factor*(y-min_csize+i)][2]-configs[m][reduction_factor*y][2]
#                    if y>n_regions/2:
#                        for i in range(2*min_csize):
#                            longax_positions[i] += configs[m][reduction_factor*(y+min_csize-i)][2]-configs[m][reduction_factor*y][2]
#                    n_longax_pos +=1
#                for i in range(2*min_csize):
#                    for j in range(2*min_csize):
#                        r_corr[i][j] += matrixofdists[m][(y+i-min_csize)%n_regions][y]*matrixofdists[m][(y+j-min_csize)%n_regions][y]
#                        if theta_defined(configs[m][reduction_factor*((y+i-min_csize)%n_regions)],configs[m][reduction_factor*((y+j-min_csize)%n_regions)], configs[m][reduction_factor*y]) ==1:
#                            theta_corr[i][j] += calc_theta(configs[m][reduction_factor*((y+i-min_csize)%n_regions)],configs[m][reduction_factor*((y+j-min_csize)%n_regions)], configs[m][reduction_factor*y])
#                            ndata_theta[i][j] +=1
#                        
#                        if phi_defined(configs[m][reduction_factor*((y+i-min_csize)%n_regions)],configs[m][reduction_factor*((y+j-min_csize)%n_regions)], configs[m][reduction_factor*y]) ==1:
#                            phi_corr[i][j] += calc_phi(configs[m][reduction_factor*((y+i-min_csize)%n_regions)],configs[m][reduction_factor*((y+j-min_csize)%n_regions)], configs[m][reduction_factor*y])
#                            ndata_phi[i][j] +=1
#    #                    print(calc_theta(configs[m][reduction_factor*(y+i-min_csize)%n_regions],configs[m][reduction_factor*(y+j-min_csize)%n_regions], configs[m][reduction_factor*y]))
#    #                    print(calc_phi(configs[m][reduction_factor*(y+i-min_csize)%n_regions],configs[m][reduction_factor*(y+j-min_csize)%n_regions]))
#                    r_av[i] += matrixofdists[m][(y+i-min_csize)%n_regions][y]
#                    r_sq_av[i] += matrixofdists[m][(y+i-min_csize)%n_regions][y]**2
#                n_domains_corr += 1
#    
#        npeaks_all.append(len(peaks_reduced))



class StatsNearFar:
    global n_labels
    def __init__(self):
        self.avs_labels = [[[] for x in range(2)] for y in range(n_labels)]
        self.stds_labels_low = [[[] for x in range(2)] for y in range(n_labels)]
        self.stds_labels_high = [[[] for x in range(2)] for y in range(n_labels)]
        self.means_nearpol = [0]*n_bins
        self.means_farpol = [0]*n_bins
        self.distr_nearpol = [[0 for y in range(cell_length)] for x in range(n_bins)]
        self.distr_farpol = [[0 for y in range(cell_length)] for x in range(n_bins)]
        self.distr_all = [[0 for y in range(cell_length)] for x in range(n_bins)]
    
    def compute_stats(self,data):
        dists = [[[] for x in range(2)] for y in range(n_labels)]
        for configindex, config in enumerate(data):
            for strandindex,strand in enumerate(config):
                if strandindex==0:
                    for labelindex, label in enumerate(label_binnrs):
                        if label>lin_length and label<n_bins-lin_length:
                            dists1 = [abs(-min_l-strand[label][2]),abs(max_l-strand[label][2])]
                            dists[labelindex][0].append(min(dists1))
                        else:
                           dists12 = np.array([[abs(-min_l-strand[label][2]),abs(max_l-strand[label][2])],[abs(-min_l-config[1][label][2]),abs(max_l-config[1][label][2])]])
                           mincoord = np.unravel_index(np.argmin(dists12, axis=None), dists12.shape)
        #                   print(label,dists12,mincoord)
                           dists[labelindex][0].append(dists12[mincoord])
                           dists[labelindex][1].append(dists12[(mincoord[0]+1)%2][mincoord[1]])
        #                   print(dists[labelindex][0][-1],dists[labelindex][1][-1])
                           
                        
                for monindex, mon in enumerate(strand):
                    if strandindex==0:
                        self.means_nearpol[monindex] += mon[2]/n_configs
                        self.distr_nearpol[monindex][mon[2]+min_l] += 1/n_configs
                    else:
                        self.means_farpol[monindex] += mon[2]/n_configs
                        self.distr_farpol[monindex][mon[2]+min_l] += 1/n_configs
        for xind, x in enumerate(self.means_farpol):
            if (xind>lin_length and xind<n_bins-lin_length):
                self.means_farpol[xind] = self.means_nearpol[xind]
                self.distr_farpol[(xind)%n_bins] = self.distr_nearpol[(xind)%n_bins]
        self.avs_labels = [[np.average(x[0])/cell_length,np.average(x[1])/cell_length] for x in dists]
        self.stds_labels_low =  [[math.sqrt(np.average([(i-self.avs_labels[xindex][0]*cell_length)**2 for i in x[0] if i <self.avs_labels[xindex][0]*cell_length]))/cell_length,math.sqrt(np.average([(i-self.avs_labels[xindex][1]*cell_length)**2 for i in x[1] if i <self.avs_labels[xindex][1]*cell_length]))/cell_length] for xindex, x in enumerate(dists)]
        self.stds_labels_high = [[math.sqrt(np.average([(i-self.avs_labels[xindex][0]*cell_length)**2 for i in x[0] if i >self.avs_labels[xindex][0]*cell_length]))/cell_length,math.sqrt(np.average([(i-self.avs_labels[xindex][1]*cell_length)**2 for i in x[1] if i >self.avs_labels[xindex][1]*cell_length]))/cell_length] for xindex, x in enumerate(dists)]

class StatsOldNew:
    global stats_old
    global ter
    global compute_twopoints
    global n_bins
    def __init__(self):
        self.extensions_old_all = [[] for x in range(n_bins)]
        self.extensions_new_all = [[] for x in range(n_bins)]
        self.extensions_old = [0]*n_bins
        self.extensions_stds_upper_old = [0]*n_bins
        self.extensions_stds_lower_old = [0]*n_bins
        self.extensions_new = [0]*n_bins
        self.extensions_stds_upper_new = [0]*n_bins
        self.extensions_stds_lower_new = [0]*n_bins
        self.cfs_intra_old = [[0 for x in range(n_bins)] for y in range(n_bins)]
        self.cfs_intra_new = [[0 for x in range(n_bins)] for y in range(n_bins)]
        self.cfs_inter = [[0 for x in range(n_bins)] for y in range(n_bins)]
        self.means_oldpol = [0]*n_bins
        self.means_newpol = [0]*n_bins
        self.longax_corr_old = [[0 for i in range(n_bins)] for j in range(n_bins)]
        self.longax_corr_new = [[0 for i in range(n_bins)] for j in range(n_bins)]
        self.longax_corr_inter = [[0 for i in range(n_bins)] for j in range(n_bins)]
        self.longax_av_old = [0 for i in range(n_bins)]
        self.longax_sq_av_old = [0 for i in range(n_bins)]
        self.longax_av_new = [0 for i in range(n_bins)]
        self.longax_sq_av_new = [0 for i in range(n_bins)]
        self.distr_oldpol = [[0 for y in range(cell_length)] for x in range(n_bins)]
        self.distr_newpol = [[0 for y in range(cell_length)] for x in range(n_bins)]
        self.distr_all = [[0 for y in range(cell_length)] for x in range(n_bins)]
    
    def ComputeStats(self,data,extension_width,sud_stats_old,sud_stats_new):
        flip_factor = 1
        for configindex, config in enumerate(data):
            if (config[0][0][2]<0 and config[1][0][2]>0 and config[0][ter][2]<0):
                flip_factor = -1
            else:
                flip_factor = 1
            for strandindex,strand in enumerate(config):
                if (strandindex==0 and flip_factor==1) or (strandindex==1 and flip_factor==-1):
                    if calculate_suds == True:
                        sud_stats_old.calc_sud_properties(strand)
                    for monindex, mon in enumerate(strand):
                        self.means_oldpol[monindex] += mon[2]*flip_factor/n_configs
                        self.longax_av_old[monindex] += mon[2]*flip_factor/n_configs
                        self.longax_sq_av_old[monindex] += mon[2]**2/n_configs
                        self.distr_oldpol[monindex][mon[2]*flip_factor+[max_l,min_l][int((flip_factor+1)/2)]] += 1/n_configs
                        self.extensions_old_all[monindex].append(dist(strand[monindex-extension_width],strand[(monindex+extension_width)%n_bins]))
                        
                        #calculate Hi-C maps in this representation -> from this obtain p(s) curves etc
                        if compute_twopoint==True:
                            for mon2index in range(monindex,n_bins):
                                self.cfs_intra_old[monindex][mon2index] += overlap(mon,config[strandindex][mon2index])/n_configs
                                self.cfs_intra_old[mon2index][monindex] = self.cfs_intra_old[monindex][mon2index]
                                
                                self.cfs_intra_new[monindex][mon2index] += overlap(config[(strandindex+1)%2][monindex],config[(strandindex+1)%2][mon2index])/n_configs
                                self.cfs_intra_new[mon2index][monindex] = self.cfs_intra_new[monindex][mon2index]
                                
                                self.cfs_inter[monindex][mon2index] += overlap(mon,config[(strandindex+1)%2][mon2index])/n_configs
                                self.cfs_inter[mon2index][monindex] += overlap(config[strandindex][mon2index],config[(strandindex+1)%2][monindex])/n_configs
        
                                self.longax_corr_old[monindex][mon2index] += mon[2]*config[strandindex][mon2index][2]/n_configs
                                self.longax_corr_old[mon2index][monindex]=self.longax_corr_old[monindex][mon2index]
                                
                                self.longax_corr_new[monindex][mon2index] += config[(strandindex+1)%2][monindex][2]*config[(strandindex+1)%2][mon2index][2]/n_configs
                                self.longax_corr_new[mon2index][monindex] = self.longax_corr_new[monindex][mon2index]
                                
                                self.longax_corr_inter[monindex][mon2index] += mon[2]*config[(strandindex+1)%2][mon2index][2]/n_configs
                                self.longax_corr_inter[mon2index][monindex] += config[(strandindex+1)%2][monindex][2]*config[strandindex][mon2index][2]/n_configs
                else:
                    if calculate_suds == True:
                        sud_stats_new.calc_sud_properties(strand)
                    for monindex, mon in enumerate(strand):
                        self.means_newpol[monindex] += mon[2]*flip_factor/n_configs 
                        self.longax_av_new[monindex] += mon[2]*flip_factor/n_configs
                        self.longax_sq_av_new[monindex] += mon[2]**2/n_configs
                        self.distr_newpol[monindex][mon[2]*flip_factor+[max_l,min_l][int((flip_factor+1)/2)]] += 1/n_configs
                        self.extensions_new_all[monindex].append(dist(strand[monindex-extension_width],strand[(monindex+extension_width)%n_bins]))
        for xind, x in enumerate(self.means_newpol):
            if (xind>lin_length and xind<n_bins-lin_length):
                self.means_newpol[xind] = self.means_oldpol[xind]
                self.distr_newpol[(xind+ori_new)%n_bins] = self.distr_oldpol[(xind+ori_new)%n_bins]
        self.extensions_old = [np.average(x) for x in self.extensions_old_all]
        self.extensions_stds_upper_old = [math.sqrt(np.average([(i-self.extensions_old[xindex])**2 for i in x if i>self.extensions_old[xindex]])) for xindex, x in enumerate(self.extensions_old_all)]
        self.extensions_stds_lower_old = [math.sqrt(np.average([(i-self.extensions_old[xindex])**2 for i in x if i<self.extensions_old[xindex]])) for xindex, x in enumerate(self.extensions_old_all)]
        
        self.extensions_new = [np.average(x) for x in self.extensions_new_all]
        self.extensions_stds_upper_new = [math.sqrt(np.average([(i-self.extensions_new[xindex])**2 for i in x if i>self.extensions_new[xindex]])) for xindex, x in enumerate(self.extensions_new_all)]
        self.extensions_stds_lower_new = [math.sqrt(np.average([(i-self.extensions_new[xindex])**2 for i in x if i<self.extensions_new[xindex]])) for xindex, x in enumerate(self.extensions_new_all)]
        if compute_twopoint==True:
            self.longax_corr_old = compute_longax_corr(self.longax_corr_old,self.longax_av_old,self.longax_sq_av_old)
            self.longax_corr_new = compute_longax_corr(self.longax_corr_new,self.longax_av_new,self.longax_sq_av_new)
            self.longax_corr_inter = compute_longax_corr(self.longax_corr_inter,self.longax_av_old,self.longax_sq_av_old,self.longax_av_new,self.longax_sq_av_new)


data_long = np.loadtxt(file_path +sim_file+'/snaps_configs_0.txt', dtype=int,usecols=range(3*pol_length*2*n_configs))
data_long = data_long.reshape(n_configs,2,pol_length,3)
data = preprocess_data(data_long,ori_sim)
reference_localizations = ReferenceLocalizations()
reference_localizations.read_ref_loc(reference_positions_file)

#Calculation of statistics in the near ori/far ori reference frame
stats_near_far = StatsNearFar()
stats_near_far.compute_stats(data)

#The two commands below need to be run before oriented_stats are run
sud_stats_old = SudStats()
sud_stats_new = SudStats()

#If calculate_suds is true, sud properties will be calculated via oriented_stats. This construction allows the sud properties to be calculated efficiently
oriented_stats = StatsOldNew()
oriented_stats.ComputeStats(data,extension_width,sud_stats_old,sud_stats_new)


#Commands for plotting results
plot_commands.plot_means_nearfar(stats_near_far.means_nearpol,stats_near_far.means_farpol,lin_length,min_l,max_l,genomic_length,n_bins,add_separations=True,loc_separations_near=reference_localizations.loc_separations_near,loc_separations_far=reference_localizations.loc_separations_far)
plot_commands.plot_distr_near(stats_near_far.distr_nearpol,stats_near_far.means_nearpol,min_l,max_l,genomic_length,n_bins)
plot_commands.plot_distr_far(stats_near_far.distr_farpol,stats_near_far.means_farpol,min_l,max_l,genomic_length,n_bins)

if calculate_suds == True:
    plot_commands.plot_mean_sudsizes(sud_stats_old.mean_domainsizes,sud_stats_new.mean_domainsizes,genomic_length,n_bins)

if compute_twopoint==True:
    plot_commands.plot_contact_freqs(oriented_stats.cfs_intra_old,genomic_length)
    plot_commands.plot_contact_freqs(oriented_stats.cfs_intra_new,genomic_length)
    plot_commands.plot_contact_freqs(oriented_stats.cfs_inter,genomic_length)

    plot_commands.plot_correlations(oriented_stats.longax_corr_old,genomic_length)
    plot_commands.plot_correlations(oriented_stats.longax_corr_new,genomic_length)
    plot_commands.plot_correlations(oriented_stats.longax_corr_inter,genomic_length)

plot_commands.plot_extensions(oriented_stats.extensions_old,oriented_stats.extensions_stds_lower_old,oriented_stats.extensions_stds_upper_old,oriented_stats.extensions_new,oriented_stats.extensions_stds_lower_new,oriented_stats.extensions_stds_upper_new,reference_extensions,genomic_length,n_bins,lin_length,lattice_spacing)
plot_commands.plot_local_extension_ratio(oriented_stats.extensions_old,oriented_stats.extensions_new,reference_extensions,genomic_length,n_bins)
plot_commands.plot_smooth_extension_ratio(oriented_stats.extensions_old,oriented_stats.extensions_new,reference_extensions,genomic_length,n_bins)

plot_commands.plot_means_oriented(oriented_stats.means_oldpol,oriented_stats.means_newpol,lin_length,min_l,max_l,genomic_length,n_bins,add_separations=add_separations,loc_separations_newpol=reference_localizations.loc_separations_newpol,loc_separations_oldpol=reference_localizations.loc_separations_oldpol)
plot_commands.plot_distr_oldpol(oriented_stats.distr_oldpol,oriented_stats.means_oldpol,min_l,max_l,genomic_length,n_bins)
plot_commands.plot_distr_newpol(oriented_stats.distr_newpol,oriented_stats.means_newpol,min_l,max_l,genomic_length,n_bins)


#Saving of calculation results
if save_localizations==True:
    save_localization_data([(x+min_l)/(min_l+max_l) for x in oriented_stats.oriented_near],"Localization_data/oriented_near_"+str(stage))
    save_localization_data([(x+min_l)/(min_l+max_l) for x in oriented_stats.oriented_far],"Localization_data/oriented_far_"+str(stage))
    save_localization_data([(x+min_l)/(min_l+max_l) for x in oriented_stats.means_oldpol],"Localization_data/oriented_oldpol_"+str(stage))
    save_localization_data([(x+min_l)/(min_l+max_l) for x in oriented_stats.means_newpol],"Localization_data/oriented_newpol_"+str(stage))
##Save Hi-C maps since they take long to calculate

if save_twopoint==True and compute_twopoint==True:
    save_twopoint_data(oriented_stats.cfs_intra_old,"Hi-C_maps_maxent_inter_intra/intra_old_"+str(stage))
    save_twopoint_data(oriented_stats.cfs_intra_new,"Hi-C_maps_maxent_inter_intra/intra_new_"+str(stage))
    save_twopoint_data(oriented_stats.cfs_inter,"Hi-C_maps_maxent_inter_intra/inter_"+str(stage))
    #
    save_twopoint_data(oriented_stats.longax_corr_old,"Correlations_inter_intra/intra_old_"+str(stage))
    save_twopoint_data(oriented_stats.longax_corr_new,"Correlations_inter_intra/intra_new_"+str(stage))
    save_twopoint_data(oriented_stats.longax_corr_inter,"Correlations_inter_intra/inter_"+str(stage))

if save_label_positions==True:
    with open('/media/joris/raid_data/Chromosome_Maxent/4dchrom_crescentus/distr_labels_sim/new_averages_'+str(stage)+'.txt','wb') as f:
        np.savetxt(f, stats_near_far.avs_labels, fmt='%.2f')
    with open('/media/joris/raid_data/Chromosome_Maxent/4dchrom_crescentus/distr_labels_sim/new_stds_low_'+str(stage)+'.txt','wb') as f:
        np.savetxt(f, stats_near_far.stds_labels_low, fmt='%.2f')
    with open('/media/joris/raid_data/Chromosome_Maxent/4dchrom_crescentus/distr_labels_sim/new_stds_high_'+str(stage)+'.txt','wb') as f:
        np.savetxt(f, stats_near_far.stds_labels_high, fmt='%.2f')