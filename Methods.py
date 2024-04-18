import os
import platform
import csv
import math
import shutil 
import psutil
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.libqsturng import qsturng, psturng
import numpy as np
import numpy_indexed as npi
import pandas as pd
from pyopenms import *
from datetime import *
import pandas as pd
import pingouin as pg
from hypothetical.descriptive import var
from scipy import stats, signal
from scipy.stats import zscore, t
from scipy.optimize import curve_fit
from scipy.integrate import quad, cumulative_trapezoid
from scipy.signal import find_peaks, peak_widths, welch, savgol_filter
from scipy.ndimage import gaussian_filter1d
import scikit_posthocs as sp
from itertools import combinations, combinations_with_replacement, permutations, product
import warnings
import xlsxwriter
import peakutils
import time
from pathlib import Path

def sys_checks(parent_dir,prots,centroided,sel_cutoff,p_score_cutoff,peak_RT,peps_sec,reps,nr_EICs,LOD,mzML):
    '''
    Checks system and inputs for compatibility - primarily ensures that relevant metrics are not negative 
    and that there are files to analyze
    '''
    openms_path = os.path.join(os.getcwd(), 'venv', 'Lib', 'site-packages', 'pyopenms', 'share', 'OpenMS')
    if platform.system() == 'Windows':     # Windows (either 32-bit or 64-bit)
        os.system(f'set OPENMS_DATA_PATH={openms_path};%OPENMS_DATA_PATH%')
        os.environ["OPENMS_DATA_PATH"] = openms_path # try 2 ways so maybe it doesn't get mad
    elif platform.system() == "Linux" or "Darwin":      # linux or Mac OS (Darwin)
        print(f'Warning: MSConvert not compatible with {platform.system()}')
    else:
        raise Exception("Error identifying operating system")

    try:
        ### Check inputs and raise warnings/exceptions
        mem = psutil.virtual_memory().total
        if mem < 16E9:
            warnings.warn(f'16GB of RAM is recommended, current system has {mem}')
        check_files = []
        for file in sorted(os.listdir(parent_dir)):
            if file.endswith('.raw') and any(prot in file for prot in prots):
                check_files.append(file)
        print(f'RAW files found to analyze: {check_files}')
        for prot in prots:
            if not any(prot in file for file in check_files):
                warnings.warn(f'RAW files for {prot} not found in data directory')
        if not check_files and not centroided:
            mess = 'Data directory does not contain any RAW files'
            warnings.warn('Input error:',mess)
            raise Exception(mess)
        if sel_cutoff > 1 or sel_cutoff < 0:
            mess = 'Selectivity score cutoff out of range (0 to 1)'
            warnings.warn('Input error:',mess)
            raise Exception(mess)
        if p_score_cutoff < 0 or p_score_cutoff > 1:
            mess = 'P score cutoff out of range (0 to 1)'
            warnings.warn('Input error:',mess)
            raise Exception(mess)
        if peak_RT < 0:
            mess = 'Peak retention time is less than 0'
            warnings.warn('Input error:',mess)
            raise Exception(mess)
        if peps_sec < 1:
            mess = "Max peptides sequenced by Orbi less than 1"
            warnings.warn('Input error:',mess)
            raise Exception(mess)
        if reps < 1:
            mess = 'Number of replicates is less than 1'
            warnings.warn('Input error:',mess)
            raise Exception(mess)
        if nr_EICs < 1:
            mess = 'Number of EICs specified is less than 1'
            warnings.warn('Input error:',mess)
            raise Exception(mess)
        if LOD < 0:
            mess = 'Minimum intensity for peak detection below 0'
            warnings.warn('Input error:',mess)
            raise Exception(mess)
    except:
        print('Error in inputs! Please check your entries')
        
def prime_OpenMS(primer_dir):
    '''
    Use this to initialize OpenMS, as the kernel will instantly die without doing it for an unknown reason
    Requires file named "Primer.mzML" in the primer directory (primer_dir)
    '''
    exp_raw = MSExperiment()
    MzMLFile().load(os.path.join(primer_dir,'Primer.mzML'), exp_raw)
    exp_centroid = MSExperiment()
    pickme = PeakPickerHiRes()
    try:
        pickme.pickExperiment(exp_raw, exp_centroid,True)
        MzMLFile().store(os.path.join(primer_dir,'Primer Centroid.mzML'), exp_centroid)  
    except:
        pickme.pickExperiment(exp_raw, exp_centroid)
        MzMLFile().store(os.path.join(primer_dir,'Primer Centroid.mzML'), exp_centroid)          
        
class Setup:
    def __init__(self):
        pass
    
    def directory_setup(self,parent_dir,folder,save_dir,prots,centroided,PRTC_check):
        '''
        Uses the input directories to set-up the location of the RAW data, as well as the proper directories
        to save all outputs 
        '''
        if not os.path.isdir(save_dir):         # Checks if data is already centroided, and if so it will create save directory
            os.mkdir(save_dir)                  # (which is same as data directory if centroided = False)

        png_dir = os.path.join(parent_dir,folder,'Data Checks and Visualizations')
        if not os.path.isdir(png_dir):
            os.mkdir(png_dir)
            
        misc_dir = os.path.join(parent_dir,folder,'Supporting Info')
        if not os.path.isdir(misc_dir):
            os.mkdir(misc_dir)

        results_dir = os.path.join(parent_dir,folder,'Results')
        if not os.path.isdir(results_dir):
            os.mkdir(results_dir)
        
        if PRTC_check:
            PRTC_dir = os.path.join(parent_dir,folder,'PRTC Data')
            if not os.path.isdir(PRTC_dir):
                os.mkdir(PRTC_dir)
        else:
            PRTC_dir = None

        return png_dir,PRTC_dir,results_dir, misc_dir
    
    def export_inputs(self,save_dir,prots,selection,sel_cutoff,p_score_cutoff,peak_RT,peak_RT_search,peps_sec,n_reps,
                     lib_size,conc_lib_start,default_thresh,area_limit,sel_minqual,gradient_time,peps_min,LOD,centroided,
                     centroid_dir,check_data,full_out,nr_EICs,sub_EICs,RT_start,elute_window):
        '''
        Takes the relevant inputs and packages them into a text file for reference
        '''
        inputs_file_path = os.path.join(save_dir, "Inputs.txt")
        with open(inputs_file_path, 'w') as f:
            f.write("Proteins Analyzed: {}\n".format(prots))
            f.write("Name of selection: '{}'\n".format(selection))
            f.write("Minimum enrichment score cutoff: {}\n".format(sel_cutoff))
            f.write("Maximum p value cutoff: {}\n".format(p_score_cutoff))
            f.write("Estimated max peak retention time width (seconds): {}\n".format(peak_RT))
            f.write("Time window used to find feature peak (seconds): {}\n".format(peak_RT_search))
            f.write("Estimated largest number of peptides that can coelute and be identified: {}\n".format(peps_sec))
            f.write("Number of replicates per protein: {}\n".format(n_reps))
            f.write("Length of peptides in library: {}\n".format(lib_size))
            f.write("Starting concentration of library peptides (pM/member): {}\n".format(conc_lib_start))
            f.write("Peak detection threshold (counts): {}\n".format(default_thresh))
            f.write("Peak area limit: {}\n".format(area_limit))
            f.write("Minimum enrichment score: {}\n".format(sel_minqual))
            f.write("Total time of gradient (min): {}\n".format(gradient_time))
            f.write("Time of MS on (min): {}\n".format(RT_start))
            f.write("Time window used for inclusion list retention time and EICs: {}\n".format(elute_window))
            f.write("Estimated number of peptides that can be analyzed by spectrometer per minute: {}\n".format(peps_min))
            f.write("Limit of detection of spectrometer (counts): {}\n".format(LOD))
            f.write("Is the data centroided (necessary for usage on Unix): {}\n".format(centroided))
            f.write("Location of centroided data: '{}'\n".format(centroid_dir))
            f.write("Do data visualizations: {}\n".format(check_data))
            f.write("Export full results (including features below score thresholds): {}\n".format(full_out))
            f.write("Number of EICs performed per protein: {}\n".format(nr_EICs))
            f.write("Perform additional EICs based only on enrichment score: {}\n".format(sub_EICs))
        print("Inputs exported to: {}".format(inputs_file_path))

    def data_centroiding(self,centroided,parent_dir,folder,PRTC_dir,data_dir,prots,centroid_dir = None,refname = 'Beads'):
        '''
        Performs centroiding on profile mzML data
        '''
        if centroided:
            if not centroid_dir:
                print('Error: requires input for centroid_dir')
            data_dir = centroid_dir
            PRTC_dir = parent_dir
            return data_dir
        else:
            if not os.path.isdir(data_dir):         # Checks that the data directory exists
                os.mkdir(data_dir)       

            if os.path.exists(os.path.join(parent_dir,folder,'centroided')):   # If the same data directory exists, delete and remake it 
                shutil.rmtree(os.path.join(parent_dir,folder,'centroided'))    #   (though it shouldn't thanks to the time stamp)
            os.mkdir(os.path.join(parent_dir,folder, 'centroided'))

            count = 1
            for file in sorted(os.listdir(parent_dir)):               # Iterate through each folder and centroid the data
                if file.endswith('.mzML') and any(prot in file for prot in prots + [refname]):
                    exp_raw = MSExperiment()
                    MzMLFile().load(os.path.join(parent_dir,file), exp_raw)

                    print(f'{file} loaded')
                    exp_centroid = MSExperiment()
                    pickme = PeakPickerHiRes()
                    params_pick = pickme.getParameters()
                    params_pick.setValue('signal_to_noise',0.0)  # 0 (disabled) is default
                    params_pick.setValue('spacing_difference_gap',4.0) # 4 is default
                    params_pick.setValue('spacing_difference',1.5) # 1.5 is default
                    params_pick.setValue('missing',1) # 1 is default
                    params_pick.setValue('report_FWHM','true') # false is default
                    params_pick.setValue('ms_levels',[1])
                    pickme.setParameters(params_pick)
                    try:
                        pickme.pickExperiment(exp_raw, exp_centroid)
                    except:
                        pickme.pickExperiment(exp_raw, exp_centroid, True)
                    print('peaks picked')
                    if 'PRTC_i' in file:
                        MzMLFile().store(os.path.join(PRTC_dir, file), exp_centroid)       
                    else:
                        MzMLFile().store(os.path.join(data_dir, 'centroided', file), exp_centroid)       
                    del exp_raw 
                    print('another one done')
                    print(count)
                    count -= -1
            data_dir = os.path.join(parent_dir,folder,'centroided')
            return data_dir

    def data_visualization(self,check_data,data_dir,png_dir,parent_dir):
        '''
        Plots the TIC for each run of both the profile and centroided versions of the data
        '''
        rt_cent = []
        int_cent = []
        if check_data:
            ## plotting TIC 
            for file in sorted(os.listdir(data_dir)):
                if file.endswith('.mzML'):
                    exp_test = MSExperiment()
                    MzMLFile().load(os.path.join(data_dir,file), exp_test)
                    print('loaded')

                    # choose one of the following three methods to access the TIC data
                    # 1) recalculate TIC data with the calculateTIC() function
                    tic = exp_test.calculateTIC()
                    retention_times, intensities = tic.get_peaks()
                    rt_cent.append(retention_times)
                    int_cent.append(intensities)

                    # plot retention times and intensities and add labels
                    fig, ax = plt.subplots(1,1,figsize=(36,12))
                    ax.plot(retention_times, intensities)
                    print(f'Fifth percentile = {np.percentile(intensities,5,axis=0)}')
                    ax.set_title(f'TIC of {file} (Centroided)',fontsize=36)
                    ax.set_xlabel('time (s)',fontsize=30)
                    ax.set_ylabel('intensity (cps)',fontsize=30)
                    ax.tick_params(axis='both',labelsize=24)
                    plt.show()
                    fig.savefig(os.path.join(png_dir,file[:-5] + ' TIC (centroided).png'))


    ## plotting TIC 
        rt_pro = []
        int_pro = []
        names = []
        if check_data:
            for file in sorted(os.listdir(parent_dir)):
                if file.endswith('.mzML'):
                    exp_test = MSExperiment()
                    MzMLFile().load(os.path.join(data_dir,file), exp_test)
                    print('loaded')

                    # choose one of the following three methods to access the TIC data
                    # 1) recalculate TIC data with the calculateTIC() function
                    tic = exp_test.calculateTIC()
                    retention_times, intensities = tic.get_peaks()
                    rt_pro.append(retention_times)
                    int_pro.append(intensities)
                    names.append(file[:-5])

                    # plot retention times and intensities and add labels
                    fig, ax = plt.subplots(1,1,figsize=(36,12))
                    ax.plot(retention_times, intensities)
                    print(f'Fifth percentile = {np.percentile(intensities,5,axis=0)}')
                    ax.set_title(f'TIC of {file} (Profile)',fontsize=26)
                    ax.set_xlabel('time (s)',fontsize=30)
                    ax.set_ylabel('intensity (cps)',fontsize=30)
                    ax.tick_params(axis='both',labelsize=24)
                    fig.savefig(os.path.join(png_dir,file[:-5] + ' TIC (profile).png'))
                    plt.show()
            for i in range(len(names)):
                x = rt_pro[i]
                y = [a_i - b_i for a_i, b_i in zip(int_pro[i], int_cent[i])]
                plt.plot(x,y)
                plt.xlabel('Retention time (sec)')
                plt.ylabel('Residual')
                plt.title('TIC Residuals')
                plt.show()
                plt.savefig(os.path.join(png_dir,names[i] + ' Residuals.png'), bbox_inches="tight")

        '''centroided'''
        if check_data:
            for file in sorted(os.listdir(data_dir)):
                if file.endswith('mzML'):
                    exp_2D = MSExperiment()
                    MzMLFile().load(os.path.join(data_dir,file),exp_2D)
                    self.plot_spectra_2D_overview(exp_2D,file[:-5],'centroided',png_dir)

        '''profile'''
        if check_data:
            for file in sorted(os.listdir(parent_dir)):
                if file.endswith('mzML'):
                    exp_2D = MSExperiment()
                    MzMLFile().load(os.path.join(data_dir,file),exp_2D)
                    self.plot_spectra_2D_overview(exp_2D,file[:-5],'profile',png_dir)


    def plot_spectra_2D_overview(self,exp,title,other,png_dir):
        '''
        Utilizes bilinear interpolation to graph a 2D RT vs m/z plot.
        '''                                          
        rows = 200.0                            
        cols = 200.0 
        exp.updateRanges()

        bilip = BilinearInterpolation()
        tmp = bilip.getData()
        tmp.resize(int(rows), int(cols), float())
        bilip.setData(tmp)
        bilip.setMapping_0(0.0, exp.getMinRT()/60, rows-1, exp.getMaxRT()/60)
        bilip.setMapping_1(0.0, exp.getMinMZ(), cols-1, exp.getMaxMZ())
        print('collecting peak data...')
        for spec in exp:
            if spec.getMSLevel() == 1:
                mzs, ints = spec.get_peaks()
                rt = spec.getRT()/60
                for i in range(0, len(mzs)):
                    bilip.addValue(rt, mzs[i], ints[i])

        data = np.ndarray(shape=(int(cols), int(rows)), dtype=np.float64)
        for i in range(int(rows)):
            for j in range(int(cols)):
                data[i][j] = bilip.getData().getValue(i,j)
        mean = np.mean(data)
        std_dev = np.std(data)
        z_score = 2
        cutoff = mean + z_score*std_dev
        plt.rcParams["figure.figsize"] = (6,6)
        plt.imshow(data, cmap='Blues',aspect=len(data[0])/len(data),vmin=0,vmax=cutoff)
        plt.ylabel('retention time (min)',fontsize=14)
        plt.xlabel('m/z',fontsize=14)
        plt.title(title,fontsize=18)
        plt.yticks(np.linspace(0,int(rows),7, dtype=int),
                np.linspace(0,120,7, dtype=int),fontsize=12)
        plt.xticks(np.linspace(0,int(cols),7, dtype=int),
                np.linspace(200,1400,7, dtype=int),fontsize=12)
        filename = title + ' 2D map ' + other + '.png'
        plt.savefig(os.path.join(png_dir,filename), bbox_inches="tight")
        plt.show()
        print('showing plot...')

    def alignment_check(self,consensus_map,feature_maps,png_dir,new_map_RTs,files,original_RTs_dict):
        '''
        Checks the movement of identified features in a consensus map to see alignment across runs
        '''
        feature_map_RTs_aligned = []
        consensus_RTs_aligned = []
        consensus_mzs_aligned = [cf.getMZ() for cf in consensus_map]
        for i,fm in enumerate(feature_maps):
            consensus_RTs_aligned.append([feat[i] for feat in new_map_RTs])
            feature_map_RTs_aligned.append([f.getRT() for f in fm])
            plt.plot([f.getRT() for f in fm],[f.getMZ() for f in fm],'k.',markersize=10,label='feature map RTs')
            plt.plot([feat[i] for feat in new_map_RTs],consensus_mzs_aligned,'y.',label='aligned feature RTs')
            plt.xlabel('Retention Time (min)')
            plt.ylabel('m/z')
            plt.savefig(os.path.join(png_dir,'Alignment Check Should Overlap ' + files[i] + '.png'), bbox_inches="tight")
            plt.show()

        fixed_RT_maps = []
        for i,RT_map in enumerate(consensus_RTs_aligned):
            fixed_RT_map = RT_map
            feature_map_RTs_original = original_RTs_dict[i]
            feature_map_RTs_aligned = [f.getRT() for f in feature_maps[i]]
            mzs_feature_map = [f.getMZ() for f in feature_maps[i]]
            for j,RT in enumerate(RT_map):
                if RT in feature_map_RTs_aligned:
                    idx = feature_map_RTs_aligned.index(RT)
                    fixed_RT_map[j] = feature_map_RTs_original[idx]
            fixed_RT_maps.append(fixed_RT_map)

        for i,RT_map in enumerate(fixed_RT_maps):
            plt.plot(RT_map,consensus_mzs_aligned,'k.',label='aligned feature RTs')
            plt.plot([f.getRT() for f in feature_maps[i]],[f.getMZ() for f in feature_maps[i]],'y.',label='original feature RTs')
            plt.xlim([0,max(RT_map)])
            plt.ylim([200,900])
            plt.xlabel('Retention Time (min)')
            plt.ylabel('m/z')
            plt.savefig(os.path.join(png_dir,'Alignment Check Should Not Overlap ' + files[i] + '.png'), bbox_inches="tight")
            plt.show()

        fixed_RT_maps_formatted = []
        for i in range(len(new_map_RTs)):
            feat_RTs = [RT[i] for RT in fixed_RT_maps]
            fixed_RT_maps_formatted.append(feat_RTs)

        return fixed_RT_maps_formatted

    def PRTC_stats(self,areas_ref_full,files,ref_feature_mass,png_dir):
        '''
        Performs analysis of doped peptide retention time calibration mix to be used for signal normalization
        '''
        for i,areas in enumerate(areas_ref_full):
            z_score = t.ppf(1-0.05,len(areas))
            std = np.std(areas)
            upper = [area+z_score*std/len(areas) for area in areas]
            lower = [area-z_score*std/len(areas) for area in areas]
            plt.plot(files,areas, label=f'PRTC mass {np.round(ref_feature_mass[i],2)}')
            plt.plot(files,lower,'k--',label=f'95% CI Lower bound')
            plt.plot(files,upper,'k--',label=f'95% CI Upper bound')
            plt.xlabel('Replicate Name')
            plt.ylabel('Area')
            plt.title(f'PRTC Normalization Areas, RSD = {std/np.average(areas)*100:.1f}%'.format())
            plt.legend(bbox_to_anchor=(1.04,0.5), loc="center left", borderaxespad=0)
            plt.xticks(rotation = 90)
            plt.savefig(os.path.join(png_dir,f'PRTC Normalization {np.round(ref_feature_mass[i],3)}.png'),bbox_inches='tight')
            plt.show()
            
    def feature_finding(self,data_dir,save_dir,prots):
        '''
        Performs feature finding on all identified centroided mzML files. Individual parameters were optimized to give the 
        highest number of peptide feature identifications
        '''
        column_names = []
        feature_maps = []
        files = []
        files_full = []
        # paper on feature finder: https://doi.org/10.1021/pr300992u

        for file in sorted(os.listdir(data_dir)):
            if file.endswith('.mzML') and any(prot in file for prot in prots):
                print(file)
                if not 'PRTC' in file:
                    files.append(file[:-5])
                    files_full.append(file)
                exp = MSExperiment()
                out = FeatureXMLFile()
                MzMLFile().load(os.path.join(data_dir,file),exp)
                exp.updateRanges()

                feature_finder = FeatureFinder()       
                # Get and set parameters 
                params = feature_finder.getParameters('centroided')
                params.setValue('mass_trace:mz_tolerance',0.004) # default 0.004, old opt 0.004
                params.setValue('isotopic_pattern:mz_tolerance',0.010) # default 0.005, old opt 0.010
                params.setValue('isotopic_pattern:charge_low',2)
                params.setValue('isotopic_pattern:charge_high',5)
                params.setValue('feature:max_rt_span',3.0) # default 2.5 - try 1.5,2,3 - done, old opt 3.0
                params.setValue('mass_trace:min_spectra',9) # default 10, try 9,8,7,6,5 - done, old opt 9 - 7 did nothing
                params.setValue('feature:rt_shape', 'asymmetric')  # default symmetric - done, old opt asymmetric
                params.setValue('seed:min_score',0.5)  # default 0.8, try 0.5,0.6,0.7 - done, old opt 0.5
                params.setValue('feature:min_score',0.5) # default 0.7, try 0.6, 0.5 - done, old opt 0.5
                params.setValue('mass_trace:max_missing',4)   # 1 is default, try 2,3,4,5 gave no features - done, old opt 4
                print(params.items())

                # Run feature finder and store as featureXML and in array
                if 'PRTC' in file:
                    PRTC_feature_map = FeatureMap()
                    feature_finder.run('centroided', exp, PRTC_feature_map, params, FeatureMap())
                    PRTC_feature_map.setPrimaryMSRunPath([str.encode(file[:-5])])
                    fXML = FeatureXMLFile()
                    file_name = file[:-5] + '.featureXML'
                    fXML.store(os.path.join(PRTC_dir,file_name),PRTC_feature_map)
                else:
                    feature_map = FeatureMap()
                    feature_finder.run('centroided', exp, feature_map, params, FeatureMap())
                    feature_map.setPrimaryMSRunPath([str.encode(file[:-5])])
                    feature_maps.append(feature_map)
                    fXML = FeatureXMLFile()
                    file_name = file[:-5] + '.featureXML'
                    fXML.store(os.path.join(save_dir,file_name),feature_map)
                print('feature map done')
                print(' ')
        return column_names, feature_maps, files, files_full
    
    def feature_alignment(self, feature_maps, files_all, inj_level):
        '''
        Aligns the individual feature maps in retention time based on the experiment with the greatest number of 
        features. 
        '''
        # Use replicate with most features as reference
        ref_index = [i[0] for i in sorted(enumerate([fm.size() for fm in feature_maps]), key=lambda x: x[1])][-1]
        ref_RTs = [f.getRT() for f in feature_maps[ref_index]]
        aligner = MapAlignmentAlgorithmPoseClustering()
        aligner.setReference(feature_maps[ref_index])
        align_params = aligner.getParameters()
        align_params.setValue('superimposer:mz_pair_max_distance', 0.5) 
        align_params.setValue('pairfinder:distance_RT:max_difference',300.00)
        align_params.setValue('superimposer:max_shift',2000.0)
        aligner.setParameters(align_params)

        for feature_map in feature_maps[:ref_index] + feature_maps[ref_index + 1:]:
            trafo = TransformationDescription()
            aligner.align(feature_map, trafo)
            transformer = MapAlignmentTransformer()
            transformer.transformRetentionTimes(feature_map, trafo, True)  # store original RT as meta value

        original_RTs_dict = {}
        original_RTs_dict_all = {}
        j = 0
        for i,fm in enumerate(feature_maps):
            if i == ref_index:
                original_RTs_dict_all[i] = ref_RTs
                if inj_level in files_all[i]:
                    original_RTs_dict[j] = ref_RTs
                    j += 1
            else:
                original_RT = [f.getMetaValue('original_RT') for f in fm]
                original_RTs_dict_all[i] = original_RT
                if inj_level in files_all[i]:
                    original_RTs_dict[j] = original_RT
                    j += 1

        alignment_dict = {}  # contains matrices of form [original_RT, aligned_RT]
        for i,dic in enumerate(original_RTs_dict_all.values()):
            feature_map = feature_maps[i]
            RTs_array = []
            for j,RT_orig in enumerate(dic):
                pair = [RT_orig,feature_map[j].getRT()]
                RTs_array.append(pair)
            alignment_dict[files_all[i]] = RTs_array
            
        return original_RTs_dict, alignment_dict
    
    def feature_grouping(self,feature_maps,save_dir):
        '''
        Groups features into a consensus map. Options used based on recommendations from DOI: 10.1021/pr300992u 
        '''
        # 
        feature_grouper = FeatureGroupingAlgorithmQT()  # Uses a quality threshold feature grouper
        grouper_params = feature_grouper.getParameters()
        grouper_params.setValue('distance_MZ:max_difference',0.01)
        grouper_params.setValue('distance_RT:max_difference',150.0)
        feature_grouper.setParameters(grouper_params)

        consensus_map = ConsensusMap()
        file_descriptions = consensus_map.getColumnHeaders()

        for i, feature_map in enumerate(feature_maps):
            file_description = file_descriptions.get(i, ColumnHeader())
            file_description.filename = feature_map.getMetaValue('spectra_data')[0].decode()
            file_description.size = feature_map.size()
            file_description.unique_id = feature_map.getUniqueId()
            file_descriptions[i] = file_description

        consensus_map.setColumnHeaders(file_descriptions)
        feature_grouper.group(feature_maps, consensus_map)

        ConsensusXMLFile().store(os.path.join(save_dir,'consensusmap.consensusXML'),consensus_map)
        return consensus_map, file_descriptions
    
    def consensus_alignment_check(self, consensus_map, feature_maps):
        '''
        Checks alignment of consensus map features based on the feature alignment
        '''
        con_feat = []               
        for cf in consensus_map:
            RT_cf = cf.getRT()
            Qual_cf = cf.getQuality()
            Charge_cf = cf.getCharge()
            MZ_cf = cf.getMZ()
            con_feat.append([cf.getRT(),cf.getCharge(),cf.getMZ(),cf.getQuality()])

        map_idxs = []
        map_RTs = []

        for cf in consensus_map:
            map_idx = []
            map_RT = []
            for f in cf.getFeatureList():
                map_idx.append(f.getMapIndex())
                map_RT.append(f.getRT())
            map_idxs.append(map_idx)
            map_RTs.append(map_RT)

        new_map_idxs = []
        new_map_RTs = []   # these are the RTs from each feature map post-alignment

        for i in range(len(map_RTs)):
            feat_idxs = map_idxs[i]
            feat_RTs = map_RTs[i]
            new_map_idx = []
            new_map_RT = []
            for i in range(len(feature_maps)):
                new_map_idx.append(i)
                if i in feat_idxs:
                    RT = feat_RTs[feat_idxs.index(i)]
                    new_map_RT.append(RT)
                if i not in feat_idxs:
                    new_map_RT.append(0)   # this is where the problem is - if feature not in replicate, puts 0 for RT
            new_map_idxs.append(new_map_idx)
            new_map_RTs.append(new_map_RT)
            
        return con_feat, map_idxs, map_RTs, new_map_idxs, new_map_RTs
    
    def consensus_feature_details(self,feature_maps,con_feat,fixed_RT_maps_formatted,new_map_RTs,parent_dir,folder):
        '''
        Exports the consensus features to a csv file for review
        '''
        col_names = []
        for fm in feature_maps:
            col_names.append(fm.getMetaValue('spectra_data')[0].decode())

        cf_df_1 = pd.DataFrame(con_feat,columns=['RT consensus','charge','m/z','quality'])
        cf_df_3 = pd.DataFrame(fixed_RT_maps_formatted,columns=col_names)
        cf_df_2 = pd.DataFrame(new_map_RTs,columns=col_names)

        cf_df_final = pd.concat([cf_df_1,cf_df_2,cf_df_3],axis=1)
        cf_df_final_np = cf_df_final.to_numpy()
        cf_df_export = pd.concat([cf_df_1,cf_df_2],axis=1)

        cf_df_export.to_csv(os.path.join(parent_dir,folder,'Consensus Features.csv'),index=False)

        feat_RT_cf_combined = []
        feat_mz_combined = []
        feat_z_combined = []
        feat_quality_combined = []
        feat_obs_mass_combined = []
        RTs = []

        for entry in cf_df_final_np:
            feat_RT_cf_combined.append(entry[0])
            feat_mz_combined.append(entry[2])
            feat_z_combined.append(entry[1])
            feat_obs_mass_combined.append(entry[2]*entry[1]-entry[1]*1.0074)
            feat_quality_combined.append(entry[3])
            RTs.append(entry[4:])

        
        return cf_df_final, cf_df_final_np, col_names, feat_RT_cf_combined,feat_mz_combined,feat_z_combined,feat_quality_combined,feat_obs_mass_combined,RTs
    
    def consensus_feature_details_massfilt(self,feat_obs_mass_combined,PRTC_mass_round,cf_df_final_np,fixed_RT_maps_formatted,lib_size, PRTC_check = False):
        '''
        Exports list of consensus features after using a crude mass filter based on a minimum mass of a peptide consisting of 
        only Gly and a maximum mass of a peptide consisting of only Trp
        '''
        feat_RT_cf_combined_massfiltered = []
        feat_RT_orig_massfiltered = []
        feat_mz_combined_massfiltered = []
        feat_z_combined_massfiltered = []
        feat_quality_combined_massfiltered = []
        feat_obs_mass_combined_massfiltered = []
        RTs_massfiltered = []

        ref_feature_mass = []
        ref_feature_RT = []
        ref_feature_RT_orig = []
        refs_found = []
        removed = 0

        for i,obs_mass in enumerate(feat_obs_mass_combined):
            entry = cf_df_final_np[i]
            if PRTC_check:
                if (np.round(obs_mass,2) in PRTC_mass_round) and (entry[1] == 2):
                    ref_feature_mass.append(entry[2])
                    ref_feature_RT.append(entry[0])
                    ref_feature_RT_orig.append(fixed_RT_maps_formatted[i])
                    refs_found.append(np.round(obs_mass,2))
            if obs_mass > 57.021*lib_size + 17.98 and obs_mass < 185.94*lib_size + 18.29:  # mass filter, all gly and all trp
                feat_RT_cf_combined_massfiltered.append(entry[0])
                feat_RT_orig_massfiltered.append(fixed_RT_maps_formatted[i])
                feat_mz_combined_massfiltered.append(entry[2])
                feat_z_combined_massfiltered.append(entry[1])
                feat_obs_mass_combined_massfiltered.append(entry[2]*entry[1]-entry[1]*1.0074)
                feat_quality_combined_massfiltered.append(entry[3])
                RTs_massfiltered.append(entry[4:])
            else:
                removed += 1
        print(f'{removed} features removed due to mass out of {len(feat_obs_mass_combined)} features')
        return feat_RT_cf_combined_massfiltered,feat_RT_orig_massfiltered,feat_mz_combined_massfiltered,feat_z_combined_massfiltered,feat_z_combined_massfiltered,feat_quality_combined_massfiltered,feat_obs_mass_combined_massfiltered,RTs_massfiltered,ref_feature_mass,ref_feature_RT,ref_feature_RT_orig,refs_found

class Calculations:
    def __init__(self):
        pass
    
    def games_howell_scipy(self,names,groups,data,reps=3,alpha=0.05):  
        '''
        Games-howell post-hoc testing, does not assume equal variance or sample size between populations.
        Implementation from https://aaronschlegel.me/games-howell-post-hoc-multiple-comparisons-test-python.html
        '''
        k = len(names)
        group_ints = {name:data[i*reps:i*reps+reps] for i,name in enumerate(names)}
        group_means = dict(npi.group_by(groups, data, np.mean))
        group_obs = dict(npi.group_by(groups, data, len))
        group_variance = dict(npi.group_by(groups, data, var))
        combs = list(combinations(names, 2))

        group_comps = []
        mean_differences = []
        degrees_freedom = []
        t_values = []
        p_values = []
        std_err = []
        up_conf = []
        low_conf = []
        A = []
        B = []

        for comb in combs:
            # Mean differences of each group combination
            diff = group_means[comb[1]] - group_means[comb[0]]

            t_val, p_val = stats.ttest_ind(group_ints[comb[0]],group_ints[comb[1]],equal_var=False)

            # Standard error of each group combination
            se = np.sqrt(0.5 * (group_variance[comb[0]] / group_obs[comb[0]] + 
                                group_variance[comb[1]] / group_obs[comb[1]]))

            # Append the computed values to their respective lists.
            mean_differences.append(diff)
            t_values.append(t_val)
            p_values.append(p_val)
            std_err.append(se)
            group_comps.append(str(comb[0]) + ' : ' + str(comb[1]))
            A.append(comb[0])
            B.append(comb[1])

        return pd.DataFrame({'A': A,
                          'B': B,
                          'mean_difference': mean_differences,
                          'std_error': std_err,
                          't_value': t_values,
                          'p_value': p_values})

    def checkList(self,lst):
        '''
        Checks labels of a certain peptide feature and assigns specificity 
        '''
        if not lst: return ('Nonspecific','Nonspecific')
        else:
            ele = lst[0]
            chk = True

            # Comparing each element with first item 
            for item in lst:
                if ele != item:
                    chk = False
                    break

            if (chk == True): return (item,'Specific')
            else: return ('Nonspecific','Nonspecific')

    def find_nearest_tol(self,array, value, tol):
        '''
        Finds the element in an array nearest the specified value, within a certain tolerance. 
        If no value is close enough, returns -1
        '''
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        if (np.abs(array[idx] - value)).min() > tol:
            return -1,-1
        else:
            return array[idx],idx

    def find_nearest(self,array, value):
        '''
        Finds the element in an array nearest the specified value
        '''
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return array[idx],idx
    
    def smoother(self,ints,kernal_size):  
        '''
        Smooths a given curve using a given kernal size
        '''
        kernal = np.ones(kernal_size)/kernal_size
        return np.convolve(ints,kernal,mode='same')

    def signaltonoise(self,a, axis=0, ddof=0):
        '''
        Calculates the signal-to-noise ratio of a given series along a specified axis
        '''
        a = np.asanyarray(a)
        m = a.mean(axis)
        sd = a.std(axis=axis, ddof=ddof)
        return np.where(sd == 0, 0, m/sd)

    def get_data(self,directory,prots,inj_levels = None,refname = 'Beads'):
        '''
        Extracts retention time, intensity and m/z data from all mzML files in the given directory
        '''
        RTs = []
        ints = []
        mzs = []
        if inj_levels:
            files = []
            for file in sorted(os.listdir(directory)):
                if file.endswith('.mzML') and any(prot in file for prot in prots) and (any(level in file for level in inj_levels) or refname in file):
                    print(file[:-5])
                    files.append(file[:-5])
                    exp = MSExperiment()
                    MzMLFile().load(os.path.join(directory,file),exp)
                    RT_list_rep = []
                    int_list_rep = []
                    RTs_full = []
                    mzs_full = []
                    ints_full = []
                    for spec in exp:
                        RTs_full.append(spec.getRT())
                        mzs_scan, ints_scan = spec.get_peaks()  
                        mzs_full.append(mzs_scan)
                        ints_full.append(ints_scan)
                    RTs.append(RTs_full)
                    ints.append(ints_full)
                    mzs.append(mzs_full)   
            return mzs,RTs,ints,files
        else:
            for file in sorted(os.listdir(directory)):
                if file.endswith('.mzML') and any(prot in file for prot in prots):
                    print(file[:-5])
                    exp = MSExperiment()
                    MzMLFile().load(os.path.join(directory,file),exp)
                    RT_list_rep = []
                    int_list_rep = []
                    RTs_full = []
                    mzs_full = []
                    ints_full = []
                    for spec in exp:
                        RTs_full.append(spec.getRT())
                        mzs_scan, ints_scan = spec.get_peaks()  
                        mzs_full.append(mzs_scan)
                        ints_full.append(ints_scan)
                    RTs.append(RTs_full)
                    ints.append(ints_full)
                    mzs.append(mzs_full) 
            return mzs,RTs,ints,files
    
    def get_data_volcano(self,directory,inj_level,prots):
        '''
        Extracts retention time, intensity and m/z data from all mzML files in the given directory
        '''
        RTs = []
        ints = []
        mzs = []
        for file in sorted(os.listdir(directory)):
            if file.endswith('.mzML') and inj_level in file and any(prot in file for prot in prots):
                print(file[:-5])
                exp = MSExperiment()
                MzMLFile().load(os.path.join(directory,file),exp)
                RT_list_rep = []
                int_list_rep = []
                RTs_full = []
                mzs_full = []
                ints_full = []
                for spec in exp:
                    RTs_full.append(spec.getRT())
                    mzs_scan, ints_scan = spec.get_peaks()  
                    mzs_full.append(mzs_scan)
                    ints_full.append(ints_scan)
                RTs.append(RTs_full)
                ints.append(ints_full)
                mzs.append(mzs_full)         
        return mzs,RTs,ints

    def feature_int_extractor(self,m_z_feature_list,RT_feature_list,RTs_orig_list,mzs,RTs,ints,LOD=1.5E4,
                              noise_level=150,peak_range=120,decimal=2,percent=1,subtract=False):
        '''
        Uses a list of feature m/z and retention times to take an extracted ion chromatogram and report the highest intensity
        of the specified ion
        '''
        ints_windows = []
        rt_windows = []
        backs = []
        max_int_features = []
        for j,mz in enumerate(m_z_feature_list):       # loop through each feature
            if j % 1000 == 0:
                print(j)
            RT = RT_feature_list[j]               # each feature has mz, RT to reference
            RTs_orig = RTs_orig_list[j]

            int_list = []
            max_int_list = []
            back_feature = []
            rt_feature_window = []
            for k in range(len(RTs)):          # loop through each replicate
                RT_rep_feature = RTs[k]        # get list of retention times for given feature
                Int_rep_feature = ints[k]      # get list of ints arrays for whole replicate
                Mz_rep_feature = mzs[k]
                if RTs_orig[k] != 0:
                    RT = RTs_orig[k]
                RT_idx_low = self.find_nearest(RT_rep_feature,RT-peak_range/2)[1]
                RT_idx_high = self.find_nearest(RT_rep_feature,RT+peak_range/2)[1]  
                RT_window = RT_rep_feature[RT_idx_low:RT_idx_high+1]
                ints_slice = Int_rep_feature[RT_idx_low:RT_idx_high+1]    # should take a slice of arrays for both of these
                mzs_slice = Mz_rep_feature[RT_idx_low:RT_idx_high+1]
                mz_idx = [self.find_nearest_tol(entry,mz,10**(-decimal)/2)[1] for entry in mzs_slice] # find idx where the feature mz is in each scan
                max_int_candidates = [array[mz_idx[o]] if mz_idx[o] > 0 else 1 + np.random.randint(1,100000)/100000 for o,array in enumerate(ints_slice)] # use 15000 for baseline value
                int_feature = np.amax(max_int_candidates)
                idx_max = np.where(max_int_candidates == int_feature)[0][0]
                background_val = np.percentile(ints_slice[idx_max],percent,axis=0)    # fetches background value for the specific scan of the max
                if subtract:
                    int_feature = int_feature - background_val   
                if int_feature < LOD:
                    int_feature = LOD
                back_feature.append(background_val)
                int_list.append(max_int_candidates)
                max_int_list.append(int_feature)
                rt_feature_window.append(RT_window)
            rt_windows.append(rt_feature_window)
            ints_windows.append(int_list)
            max_int_features.append(max_int_list)
            backs.append(back_feature)

        return max_int_features,rt_windows,ints_windows,backs 
    
    def feature_int_extractor_average(self,m_z_feature_list,RT_feature_list,mzs,RTs,ints,LOD=1.5E4,
                              noise_level=150,peak_range=120,decimal=2,percent=1,subtract=False):
        '''
        Uses a list of feature m/z and retention times to take an extracted ion chromatogram and report the highest intensity
        of the specified ion
        '''
        ints_windows = []
        rt_windows = []
        backs = []
        max_int_features = []
        for j,mz in enumerate(m_z_feature_list):       # loop through each feature
            if j % 1000 == 0:
                print(j)
            RT = RT_feature_list[j]               # each feature has mz, RT to reference

            int_list = []
            max_int_list = []
            back_feature = []
            rt_feature_window = []
            for k in range(len(RTs)):          # loop through each replicate
                RT_rep_feature = RTs[k]        # get list of retention times for given feature
                Int_rep_feature = ints[k]      # get list of ints arrays for whole replicate
                Mz_rep_feature = mzs[k]
                RT_idx_low = self.find_nearest(RT_rep_feature,RT-peak_range/2)[1]
                RT_idx_high = self.find_nearest(RT_rep_feature,RT+peak_range/2)[1]  
                RT_window = RT_rep_feature[RT_idx_low:RT_idx_high+1]
                ints_slice = Int_rep_feature[RT_idx_low:RT_idx_high+1]    # should take a slice of arrays for both of these
                mzs_slice = Mz_rep_feature[RT_idx_low:RT_idx_high+1]
                mz_idx = [self.find_nearest_tol(entry,mz,10**(-decimal)/2)[1] for entry in mzs_slice] # find idx where the feature mz is in each scan
                max_int_candidates = [array[mz_idx[o]] if mz_idx[o] > 0 else 1 + np.random.randint(1,100000)/100000 for o,array in enumerate(ints_slice)] # use 15000 for baseline value
                int_feature = np.amax(max_int_candidates)
                idx_max = np.where(max_int_candidates == int_feature)[0][0]
                background_val = np.percentile(ints_slice[idx_max],percent,axis=0)    # fetches background value for the specific scan of the max
                if subtract:
                    int_feature = int_feature - background_val   
                if int_feature < LOD:
                    int_feature = LOD
                back_feature.append(background_val)
                int_list.append(max_int_candidates)
                max_int_list.append(int_feature)
                rt_feature_window.append(RT_window)
            rt_windows.append(rt_feature_window)
            ints_windows.append(int_list)
            max_int_features.append(max_int_list)
            backs.append(back_feature)

        return max_int_features,rt_windows,ints_windows,backs  

    def peak_finder_savgol(self,rt_windows,ints_windows,names,plot=False,reps=3,default_thresh=1.5E5,
                           kernal_size=10,width=5,prominence=2,threshold=2,rel_height=0.5):
        '''
        Utilizes peakutils to identify peaks for a given slice of retention time vs intensity in an EIC.
        The data has a Savitzky-Golay filter applied before peak detection to assist in identification
        '''
        all_peaks = []
        for i,feature in enumerate(ints_windows):  # each entry in ints_windows is a feature
            feature_peaks = []
            rt_feature = rt_windows[i]
            is_peaks = False
            if plot:
                fig,axs = plt.subplots(1,len(names)*reps,figsize=(16,12),sharey=True)
            for j,rep in enumerate(feature):       # each feature should have 6 reps (or however many for your selection)
                try:
                    rep_smoothed = savgol_filter(rep,19,9,mode='interp')
                    height = np.percentile(rep_smoothed,99)
                    if height < default_thresh:
                        height = default_thresh
                    peaks,props = find_peaks(rep_smoothed,width=width,prominence=prominence,threshold=threshold,height=height,rel_height=rel_height) #np.percentile(rep,99)
                    feature_peaks.append(peaks)
                except ValueError:
                    rep_smoothed = savgol_filter(rep,9,5,mode='mirror')
                    height = np.percentile(rep_smoothed,99)
                    if height < default_thresh:
                        height = default_thresh
                    peaks,props = find_peaks(rep_smoothed,width=width,prominence=prominence,threshold=threshold,height=height,rel_height=rel_height) #np.percentile(rep,99)
                    feature_peaks.append(peaks)
                if plot:
                    axs[j].plot(rt_feature[j],rep,'k-')
                    axs[j].plot(rt_feature[j],rep_smoothed,'b-')
                    axs[j].plot(rt_feature[j],np.ones(len(rt_feature[j]))*height,'g--')
                    for peak in peaks:
                        axs[j].scatter(rt_feature[j][peak],rep_smoothed[peak],c='red',s=20)
                        is_peaks = True
            if plot: 
                if is_peaks:
                    plt.show()
                plt.close(fig)
            all_peaks.append(feature_peaks)
            if i % 1000 == 0:
                print(i)
        return all_peaks
    
    def peak_filter_savgol(self,peak_regular,feat_mz_combined_massfiltered,feat_RT_cf_combined_massfiltered,feat_RT_orig_massfiltered,rt_windows,int_windows,feat_z_combined_massfiltered,max_ints):
        '''
        Takes the list of identified peaks from peak_finder_savgol and removes any features that had no discernible 
        peak within the retention time slice
        '''
        peak_filtered = []
        feat_mz_filtered = []
        feat_RT_filtered = []
        feat_RT_orig_filtered = []
        feat_z_filtered = []
        check_rt_windows = []
        check_int_windows = []
        max_ints_filtered = []
        for i,peak_group in enumerate(peak_regular):
            test = np.concatenate(peak_group)
            if np.any(test):
                peak_filtered.append(peak_group)
                feat_mz_filtered.append(feat_mz_combined_massfiltered[i])
                feat_RT_filtered.append(feat_RT_cf_combined_massfiltered[i])
                feat_RT_orig_filtered.append(feat_RT_orig_massfiltered[i])
                check_rt_windows.append(rt_windows[i])
                check_int_windows.append(int_windows[i])
                feat_z_filtered.append(feat_z_combined_massfiltered[i])
                max_ints_filtered.append(max_ints[i])

        print(f'Features removed: {len(feat_mz_combined_massfiltered) - len(feat_mz_filtered)}')
        print(f'Features remaining: {len(feat_mz_filtered)}')
        return peak_filtered,feat_mz_filtered,feat_RT_filtered,feat_RT_orig_filtered,feat_z_filtered,check_rt_windows,check_int_windows,max_ints_filtered   

    def feature_area_extractor_savgol(self,rt_windows_filt,int_windows_filt,saveme = None,mz_filtered=None,RT_filtered=None,
                                      check=False, peak_RT = 30, width_start=5,prominence=2,threshold=2,
                                      rel_height=0.5,area_baseline = 1e4):
        '''
        Takes the identified peaks and numerically calculates a peak area using a sum of cumulative trapezoids
        '''
        all_areas = []
        for i,feature in enumerate(int_windows_filt):  # each entry in ints_windows is a feature
            feature_areas = []
            rt_feature = rt_windows_filt[i]
            if check:
                fig, axs = plt.subplots(1,len(feature),figsize=(16,10),sharey=True)
            for j,rep in enumerate(feature):       # each feature should have 6 reps (or however many for your selection)
                is_peaks = False
                x_window = 19
                polyorder = 9
                if len(rep) < x_window:
                    x_window = len(rep)
                    polyorder = int(x_window/2)
                rep_smoothed = savgol_filter(rep,x_window,polyorder,mode='interp')
                height = np.percentile(rep_smoothed,99)
                width = width_start
                while not is_peaks:
                    if width < 0 or height < 0:                        
                        peak_subwindow = rt_feature[j]
                        ints_subwindow = rep_smoothed
                        area_trap = np.random.randint(area_baseline*0.97,area_baseline*1.03) # placeholder non-zero value for scoring
                        feature_areas.append(area_trap)  
                        break
                    peaks,props = find_peaks(rep_smoothed,width=[width,50],prominence=prominence,threshold=threshold,height=height,rel_height=rel_height)                
                    width = width - 1
                    height = height*0.95
                    if np.any(peaks):
                        peak_ints = rep_smoothed[peaks]
                        peak_max_idx = np.where(peak_ints == max(peak_ints))[0][0]
                        peak_max = np.where(rep_smoothed == max(peak_ints))[0]
                        peak_width = math.ceil(props['widths'][peak_max_idx])
                        time_inc = rt_feature[j][1] - rt_feature[j][0] # average time between scans
                        if peak_width*time_inc > peak_RT:   
                            peak_width = int(peak_width/time_inc)
                        peak = peak_max[0]
                        if peak - peak_width < 0:
                            peak_subwindow = rt_feature[j][0:peak+peak_width] # take points around the identified peak
                            ints_subwindow = rep_smoothed[0:peak+peak_width]
                        elif peak + peak_width > len(rep_smoothed):
                            peak_subwindow = rt_feature[j][peak-peak_width:] # take points around the identified peak
                            ints_subwindow = rep_smoothed[peak-peak_width:]
                        else:
                            peak_subwindow = rt_feature[j][peak-peak_width:peak+peak_width] # take points around the identified peak
                            ints_subwindow = rep_smoothed[peak-peak_width:peak+peak_width]
                        try:
                            area_trap = cumulative_trapezoid(ints_subwindow,x=peak_subwindow)[-1]
                        except:
                            print('Area calculation failed, setting to default')
                            area_trap = np.random.randint(1e4*0.97,1e4*1.03) # placeholder non-zero value for scoring
                        feature_areas.append(area_trap)
                        is_peaks = True
                if check:
                    axs[j].plot(rt_feature[j],rep,'k-',linewidth=2)
                    axs[j].plot(rt_feature[j],rep_smoothed,'c--')
                    axs[j].plot(peak_subwindow,ints_subwindow,'m--',label=f'A={area_trap:.2E}',linewidth=4)
#                     
                    axs[j].legend(loc='best')
            if check:                
                if saveme:
                    plt.savefig(os.path.join(saveme,f'EIC {np.round(mz_filtered[i],4)} at {RT_filtered[i]}.png'), bbox_inches = 'tight')
                plt.show()
                plt.close(fig)
            all_areas.append(feature_areas)
            if i % 1000 == 0:
                print(i)
        return all_areas

    def selectivity_scoring(self,max_int_features,ref_val,prot_names,prots,reps=3,p_score_cutoff=0.05):
        '''
        Performs statistical hypothesis testing between proteins to decide whether or not a feature is 
        enriched to a significant degree for a given protein. Features with significant signal for multiple 
        proteins are then defined as nonspecific, and ones with low p values are classified as unclear
        '''
        int_sums = []
        spec_label = []
        specificity = []
        p_vals = []
        for int_array in max_int_features:   
            int_sums.append([np.sum(int_array[reps*l:reps*l+reps]) for l in range(len(prot_names))])
            # this should condense the enrichments from each rep to an averaged one for each protein

        print('Moving on to stats')
        group = np.repeat(prots,reps)
        for i,feat in enumerate(max_int_features):
            if len(feat) != len(group):
                print(f'Mismatch in length: {feat,group}')
            df = pd.DataFrame({'score': feat,'group': group})
            variance_test = pg.homoscedasticity(data=df,dv='score',group='group')  # checking for unequal variance
            var = variance_test['equal_var'].values[0]                             # if equal, Tukey, else Games-Howell
            if var:
                stats_out = pg.welch_anova(dv='score', between='group', data=df)
                pval = stats_out['p-unc'].values
                if pval < p_score_cutoff:
                    stats_out_GH = pg.pairwise_tukey(dv='score', between='group', data=df)
                    GH_ps = stats_out_GH['p-tukey'].values
                    id_list = []
                    for i,p in enumerate(GH_ps):
                        if p <= p_score_cutoff:
                            diff = stats_out_GH['diff'].values[i]
                            if diff > 0:
                                id_list.append(stats_out_GH['A'].values[i])
                            else:
                                id_list.append(stats_out_GH['B'].values[i])
                    if not id_list:
                        spec_label.append('Nonspecific')
                        specificity.append('Nonspecific')
                    else:
                        result = self.checkList(id_list)
                        spec_label.append(result[0])
                        specificity.append(result[1]) 
                else:
                    spec_label.append('Unclear')
                    specificity.append('Unclear')
                p_vals.append(pval)
            else:
                print('UNEQUAL VARIANCES, USING GAMES-HOWELL')
                stats_out = pg.welch_anova(dv='score', between='group', data=df)
                pval = stats_out['p-unc'].values
                if pval < p_score_cutoff:
                    stats_out_GH = self.games_howell_scipy(prots,group,feat,reps=reps)
                    GH_ps = stats_out_GH['p_value'].values
                    id_list = []
                    for i,p in enumerate(GH_ps):
                        if p <= p_score_cutoff:
                            diff = stats_out_GH['mean_difference'].values[i]
                            if diff > 0:
                                id_list.append(stats_out_GH['A'].values[i])
                            else:
                                id_list.append(stats_out_GH['B'].values[i])
                    if not id_list:
                        spec_label.append('Nonspecific')
                        specificity.append('Nonspecific')
                    else:
                        result = checkList(id_list)
                        spec_label.append(result[0])
                        specificity.append(result[1])            
                else:
                    spec_label.append('Unclear')
                    specificity.append('Unclear')
                p_vals.append(pval)

        print('Moving on to scoring')
        sel_scores = []
        sel_scores_intstd = []
        for m in range(len(prot_names)):                          # now we divide and get enrichment scores
            s_prot = []
            s_prot_intstd = []
            for feature in int_sums:
                int_sum_prot = feature[m]
                s_score = int_sum_prot/np.sum(feature)
                s_prot.append(s_score)
                s_score_intstd = int_sum_prot/ref_val
                s_prot_intstd.append(s_score_intstd)
            sel_scores.append(s_prot)
            sel_scores_intstd.append(s_prot_intstd)
        return np.asarray(sel_scores),np.asarray(sel_scores_intstd),p_vals,specificity,spec_label

    def setup_output(self,prots,spec_label,selectivity,selectivity_intstd,pvals,
                     feat_RT_filtered,feat_RT_orig_filtered,feat_mz_filtered,feat_z_filtered,areas_savgol):
        '''
        Organizes output into three separate categories - protein specific features, protein nonspecific features, 
        and unclear features as defined during enrichment scoring (selectivity_scoring)
        '''
        sels_graphing = [ [] for _ in range(len(prots))]
        sels_graphing_intstd = [ [] for _ in range(len(prots))]
        ps_graphing = [ [] for _ in range(len(prots))]
        RTs_graphing = [ [] for _ in range(len(prots))]
        RTs_graphing_orig = [ [] for _ in range(len(prots))]
        mzs_graphing = [ [] for _ in range(len(prots))]
        z_graphing = [ [] for _ in range(len(prots))]
        areas_graphing = [ [] for _ in range(len(prots))]

        sels_nonspecific = []
        sels_nonspecific_intstd = []
        ps_nonspecific = []
        RTs_nonspecific = []
        RTs_nonspecific_orig = []
        mzs_nonspecific = []
        z_nonspecific = []
        areas_nonspecific = []

        sels_unclear = []
        sels_unclear_intstd = []
        ps_unclear = []
        RTs_unclear = []
        RTs_unclear_orig = []
        mzs_unclear = []
        z_unclear = []
        areas_unclear = []


        for i,label in enumerate(spec_label):
            if label in prots:
                prot_name = np.where(np.array(prots) == label)[0][0]  # actually a number lol
                sels_graphing[prot_name].append(selectivity[prot_name][i])
                sels_graphing_intstd[prot_name].append(selectivity_intstd[prot_name][i])
                ps_graphing[prot_name].append(pvals[i])
                RTs_graphing[prot_name].append(feat_RT_filtered[i])
                RTs_graphing_orig[prot_name].append(feat_RT_orig_filtered[i])
                mzs_graphing[prot_name].append(feat_mz_filtered[i])
                z_graphing[prot_name].append(feat_z_filtered[i])
                areas_graphing[prot_name].append(areas_savgol[i])
            elif label == 'Nonspecific':
                sels_nonspecific.append(selectivity[:,i])
                sels_nonspecific_intstd.append(selectivity_intstd[:,i])
                ps_nonspecific.append(pvals[i])
                RTs_nonspecific.append(feat_RT_filtered[i])
                RTs_nonspecific_orig.append(feat_RT_orig_filtered[i])
                mzs_nonspecific.append(feat_mz_filtered[i])
                z_nonspecific.append(feat_z_filtered[i])
                areas_nonspecific.append(areas_savgol[i])
            else:
                sels_unclear.append(selectivity[:,i])
                sels_unclear_intstd.append(selectivity_intstd[:,i])
                ps_unclear.append(pvals[i])
                RTs_unclear.append(feat_RT_filtered[i])
                RTs_unclear_orig.append(feat_RT_orig_filtered[i])
                mzs_unclear.append(feat_mz_filtered[i])
                z_unclear.append(feat_z_filtered[i])
                areas_unclear.append(areas_savgol[i])
        return sels_graphing,sels_graphing_intstd,ps_graphing,RTs_graphing,RTs_graphing_orig,mzs_graphing,z_graphing,areas_graphing,sels_nonspecific,sels_nonspecific_intstd,ps_nonspecific,RTs_nonspecific,RTs_nonspecific_orig, \
        mzs_nonspecific,z_nonspecific,areas_nonspecific, sels_unclear,sels_unclear_intstd,ps_unclear,RTs_unclear,RTs_unclear_orig,mzs_unclear,z_unclear,areas_unclear

    def selectivity_rankings(self,prots,sels_graphing,sels_graphing_intstd,png_dir):
        '''
        Organizes features based on enrichment score, then plots them based on ranking position.
        Also provides a zoomed in version of only the top 1000 features
        '''
        if len(prots) > 2:
            fig,axs = plt.subplots(math.ceil(len(prots)/2),2,figsize=(14,14),sharey=True)
            for i in range(math.ceil(len(prots)/2)):
                for j in range(2):
                    try:
                        scores = sels_graphing[2*i+j]
                        scores = np.sort(scores)[::-1]
                        scores_intstd = sels_graphing_intstd[2*i+j]
                        scores_intstd = np.sort(scores_intstd)[::-1]
                        axs[i,j].plot(scores,'k-',label=f'{prots[2*i+j]} weighted')
                        axs[i,j].set_yscale('linear')
                        axs[i,j].set_title(f'Selectivity of {prots[2*i+j]}', weight = 'bold', fontsize = 18)
                        axs[i,j].set_xlabel('Selectivity rank', weight = 'bold', fontsize = 14)
                        axs[i,j].set_ylabel('Selectivity score', weight = 'bold', fontsize = 14)
                        axs[i,j].grid(True, which='major', color='dimgray', linestyle='-')
                        axs[i,j].grid(True, which='minor', color='lightgray', linestyle='--')
                        axs[i,j].tick_params(axis='both', which='major', labelsize=14)
                        ax2 = axs[i,j].twinx()
                        ax2.plot(scores_intstd,'r-',label=f'{prots[2*i+j]} Ref to PRTC')
                        ax2.set_yscale('log')
                        handles1, labels1 = axs[i,j].get_legend_handles_labels()
                        handles2, labels2 = ax2.get_legend_handles_labels()
                        ax2.legend(handles1 + handles2, labels1 + labels2, loc='upper right')
                    except IndexError:
                        continue
            if len(prots) % 2 != 0:
                fig.delaxes(axs[-1][-1])
            plt.tight_layout()
            figname = os.path.join(png_dir,'Selectivity Ranking.png')
            plt.savefig(figname, bbox_inches="tight")
            plt.show()

            fig,axs = plt.subplots(math.ceil(len(prots)/2),2,figsize=(14,14),sharey=True)
            for i in range(math.ceil(len(prots)/2)):
                for j in range(2):
                    try:
                        scores = sels_graphing[2*i+j]
                        scores = np.sort(scores)[::-1]
                        scores_intstd = sels_graphing_intstd[2*i+j]
                        scores_intstd = np.sort(scores_intstd)[::-1]
                        axs[i,j].plot(scores[0:1000],'k-',label=f'{prots[2*i+j]} weighted')
                        axs[i,j].set_yscale('linear')
                        axs[i,j].set_title(f'Selectivity of {prots[2*i+j]} Zoomed', weight = 'bold', fontsize = 18)
                        axs[i,j].set_xlabel('Selectivity rank', weight = 'bold', fontsize = 14)
                        axs[i,j].set_ylabel('Selectivity score', weight = 'bold', fontsize = 14)
                        axs[i,j].grid(True, which='major', color='dimgray', linestyle='-')
                        axs[i,j].grid(True, which='minor', color='lightgray', linestyle='--')
                        axs[i,j].tick_params(axis='both', which='major', labelsize=14)
                        ax2 = axs[i,j].twinx()
                        ax2.plot(scores_intstd[0:1000],'r-',label=f'{prots[2*i+j]} Ref to PRTC')
                        ax2.set_yscale('log')
                        handles1, labels1 = axs[i,j].get_legend_handles_labels()
                        handles2, labels2 = ax2.get_legend_handles_labels()
                        ax2.legend(handles1 + handles2, labels1 + labels2, loc='upper right')
                    except IndexError:
                        continue
            if len(prots) % 2 != 0:
                fig.delaxes(axs[-1][-1])
            plt.tight_layout()
            figname = os.path.join(png_dir,'Selectivity Ranking Zoomed.png')
            plt.savefig(figname, bbox_inches="tight")
            plt.show()

        else:
            fig,axs = plt.subplots(1,2,figsize=(14,10),sharey=True)
            for j in range(2):
                scores = sels_graphing[j]
                scores = np.sort(scores)[::-1]
                scores_intstd = sels_graphing_intstd[j]
                scores_intstd = np.sort(scores_intstd)[::-1]
                axs[j].plot(scores,'k-',label=f'{prots[j]} weighted')
                axs[j].set_yscale('linear')
                axs[j].set_title(f'Selectivity of {prots[j]}', weight = 'bold', fontsize = 18)
                axs[j].set_xlabel('Selectivity rank', weight = 'bold', fontsize = 14)
                axs[j].set_ylabel('Selectivity score', weight = 'bold', fontsize = 14)
                axs[j].grid(True, which='major', color='dimgray', linestyle='-')
                axs[j].grid(True, which='minor', color='lightgray', linestyle='--')
                ax2 = axs[j].twinx()
                ax2.plot(scores_intstd,'r-',label=f'{prots[j]} Ref to PRTC')
                ax2.set_yscale('log')
                handles1, labels1 = axs[j].get_legend_handles_labels()
                handles2, labels2 = ax2.get_legend_handles_labels()
                ax2.legend(handles1 + handles2, labels1 + labels2, loc='upper right')
            plt.tight_layout()
            figname = os.path.join(png_dir,'Selectivity Ranking.png')
            plt.savefig(figname, bbox_inches="tight")
            plt.show()

            fig,axs = plt.subplots(1,2,figsize=(14,10),sharey=True)
            for j in range(2):
                scores = sels_graphing[j]
                scores = np.sort(scores)[::-1]
                scores_intstd = sels_graphing_intstd[j]
                scores_intstd = np.sort(scores_intstd)[::-1]
                axs[j].plot(scores[0:1000],'k-',label=f'{prots[j]} weighted')
                axs[j].set_yscale('linear')
                axs[j].set_title(f'Selectivity of {prots[j]} Zoomed', weight = 'bold', fontsize = 18)
                axs[j].set_xlabel('Selectivity rank', weight = 'bold', fontsize = 14)
                axs[j].set_ylabel('Selectivity score', weight = 'bold', fontsize = 14)
                axs[j].grid(True, which='major', color='dimgray', linestyle='-')
                axs[j].grid(True, which='minor', color='lightgray', linestyle='--')
                ax2 = axs[j].twinx()
                ax2.plot(scores_intstd[0:1000],'r-',label=f'{prots[j]} Ref to PRTC')
                ax2.set_yscale('log')
                handles1, labels1 = axs[j].get_legend_handles_labels()
                handles2, labels2 = ax2.get_legend_handles_labels()
                ax2.legend(handles1 + handles2, labels1 + labels2, loc='upper right')
            plt.tight_layout()
            figname = os.path.join(png_dir,'Selectivity Ranking Zoomed.png')
            plt.savefig(figname, bbox_inches="tight")
            plt.show()

    def volcano_plotting(self,prots,sels_graphing,sels_nonspecific,sels_unclear,ps_graphing,ps_nonspecific,ps_unclear,
                         feat_mz_combined,selectivity,p_score_cutoff,sel_cutoff,save_dir):
        '''
        Constructs a volcano plot comparing enrichment score on the x axis to p value on the y axis. 
        Color codes points based on protein specificity 
        '''
        if len(prots) > 2:
            fig, axs = plt.subplots(math.ceil(len(prots)/2),2,figsize=(18,18))
            for i in range(math.ceil(len(prots)/2)):
                for j in range(2):
                    try:
                        s_spec = sels_graphing[2*i+j]
                        s_nonspec = [s_val[2*i+j] for s_val in sels_nonspecific]
                        s_unclear = [s_val[2*i+j] for s_val in sels_unclear]

                        p_spec = [-np.log10(p) for p in ps_graphing[2*i+j]]
                        p_nonspec = [-np.log10(p) for p in ps_nonspecific]
                        p_unclear = [-np.log10(p) for p in ps_unclear]
                        axs[i,j].scatter(s_unclear,p_unclear,c='r',label='Unclear')
                        axs[i,j].scatter(s_nonspec,p_nonspec,c='gray',label='Nonspecific')
                        axs[i,j].scatter(s_spec,p_spec,c='b',label='Specific')   
                        axs[i,j].plot([min(s_spec + s_nonspec + s_unclear),max(s_spec + s_nonspec + s_unclear)],[-np.log10(p_score_cutoff),-np.log10(p_score_cutoff)],'k--')
                        axs[i,j].plot([sel_cutoff,sel_cutoff],[min(p_spec+p_nonspec+p_unclear),max(p_spec+p_nonspec+p_unclear)],'k--')
                        axs[i,j].set_xscale('linear')
                        axs[i,j].set_yscale('linear')
                        axs[i,j].set_title(f'P-Value vs Selectivity Score for {prots[2*i+j]}',weight = 'bold',fontsize=20)
                        axs[i,j].set_xlabel('Selectivity Score',weight = 'bold',fontsize=16)
                        axs[i,j].set_ylabel('-Log10(P value)',weight = 'bold',fontsize=16)
                        axs[i,j].tick_params(axis='both', which='major', labelsize=16)
                    except IndexError:
                        continue
            if len(prots) % 2 != 0:
                fig.delaxes(axs[-1][-1])
            plt.tight_layout(rect=[0, 0.03, 1, 0.95])
            figname = 'P-Value vs Selectivity Score.png'
            plt.savefig(os.path.join(save_dir,figname), bbox_inches="tight")
            plt.show()          
        else:
            fig, axs = plt.subplots(1,2,figsize=(14,10))
            for i in range(len(prots)):
                s_spec = sels_graphing[i]
                s_nonspec = [s_val[i] for s_val in sels_nonspecific]
                s_unclear = [s_val[i] for s_val in sels_unclear]

                p_spec = [-np.log10(p) for p in ps_graphing[i]]
                p_nonspec = [-np.log10(p) for p in ps_nonspecific]
                p_unclear = [-np.log10(p) for p in ps_unclear]
                axs[i].scatter(s_unclear,p_unclear,c='r',label='Unclear')
                axs[i].scatter(s_nonspec,p_nonspec,c='gray',label='Nonspecific')
                axs[i].scatter(s_spec,p_spec,c='b',label='Specific')   
                axs[i].plot([min(s_spec + s_nonspec + s_unclear),max(s_spec + s_nonspec + s_unclear)],[-np.log10(p_score_cutoff),-np.log10(p_score_cutoff)],'k--')
                axs[i].plot([1/len(prots),1/len(prots)],[min(p_spec+p_nonspec+p_unclear),max(p_spec+p_nonspec+p_unclear)],'k--')
                axs[i].set_xscale('linear')
                axs[i].set_xlim([0,1])
                axs[i].set_yscale('linear')
                axs[i].set_title(f'P-Value vs Selectivity Score for {prots[i]}',weight = 'bold',fontsize=20)
                axs[i].set_xlabel('Selectivity Score',weight = 'bold',fontsize=16)
                axs[i].set_ylabel('-Log10(P value)',weight = 'bold',fontsize=16)
                axs[i].legend(loc='best',fontsize=16)
                axs[i].tick_params(axis='both',which='major', labelsize=16)
            plt.tight_layout(rect=[0, 0.03, 1, 0.95])
            figname = 'P-Value vs Selectivity Score.png'
            plt.savefig(os.path.join(save_dir,figname), bbox_inches="tight")
            plt.show()

    def volcano_plots_intstd(self,prots,sels_graphing_intstd,sels_nonspecific_intstd,sels_unclear_intstd,
                                 ps_graphing,ps_nonspecific,ps_unclear,p_score_cutoff,sel_cutoff,save_dir):
        '''
        Constructs a volcano plot comparing enrichment score on the x axis to p value on the y axis. 
        Color codes points based on protein specificity. Normalizes enrichment scores based on PRTC signal
        '''
        if len(prots) > 2:
            fig, axs = plt.subplots(math.ceil(len(prots)/2),2,figsize=(18,18))
            for i in range(math.ceil(len(prots)/2)):
                for j in range(2):
                    try:
                        s_spec = sels_graphing_intstd[2*i+j]
                        s_nonspec = [s_val[2*i+j] for s_val in sels_nonspecific_intstd]
                        s_unclear = [s_val[2*i+j] for s_val in sels_unclear_intstd]
                        p_spec = [-np.log10(p) for p in ps_graphing[2*i+j]]
                        p_nonspec = [-np.log10(p) for p in ps_nonspecific]
                        p_unclear = [-np.log10(p) for p in ps_unclear]
                        axs[i,j].scatter(s_unclear,p_unclear,c='r',label='Unclear')
                        axs[i,j].scatter(s_nonspec,p_nonspec,c='gray',label='Nonspecific')
                        axs[i,j].scatter(s_spec,p_spec,c='b',label='Specific')   
                        axs[i,j].plot([min(s_spec + s_nonspec + s_unclear),max(s_spec + s_nonspec + s_unclear)],[-np.log10(p_score_cutoff),-np.log10(p_score_cutoff)],'k--')
                        axs[i,j].plot([sel_cutoff,sel_cutoff],[min(p_spec+p_nonspec+p_unclear),max(p_spec+p_nonspec+p_unclear)],'k--')
                        axs[i,j].set_xscale('log')
                        axs[i,j].set_yscale('linear')
                        axs[i,j].set_xlim([min(s_spec),max(s_spec)])
                        axs[i,j].set_title(f'P-Value vs Selectivity Score for {prots[2*i+j]}',weight = 'bold',fontsize=20)
                        axs[i,j].set_xlabel('Selectivity Score',weight = 'bold',fontsize=16)
                        axs[i,j].set_ylabel('-Log10(P value)',weight = 'bold',fontsize=16)
                        axs[i,j].tick_params(axis='both', which='major', labelsize=16)
                    except IndexError:
                        continue
            if len(prots) % 2 != 0:
                fig.delaxes(axs[-1][-1])            
            plt.tight_layout(rect=[0, 0.03, 1, 0.95])
            figname = 'P-Value vs Selectivity Score normalized.png'
            plt.savefig(os.path.join(save_dir,figname), bbox_inches="tight")
            plt.show()          
        else:
            fig, axs = plt.subplots(1,2,figsize=(14,10))
            for i in range(len(prots)):
                s_spec = sels_graphing_intstd[i]
                s_nonspec = [s_val[i] for s_val in sels_nonspecific_intstd]
                s_unclear = [s_val[i] for s_val in sels_unclear_intstd]

                p_spec = [-np.log10(p) for p in ps_graphing[i]]
                p_nonspec = [-np.log10(p) for p in ps_nonspecific]
                p_unclear = [-np.log10(p) for p in ps_unclear]
                axs[i].scatter(s_unclear,p_unclear,c='r',label='Unclear')
                axs[i].scatter(s_nonspec,p_nonspec,c='gray',label='Nonspecific')
                axs[i].scatter(s_spec,p_spec,c='b',label='Specific')   
                axs[i].plot([min(s_spec + s_nonspec + s_unclear),max(s_spec + s_nonspec + s_unclear)],[-np.log10(p_score_cutoff),-np.log10(p_score_cutoff)],'k--')
                axs[i].plot([1,1],[min(p_spec+p_nonspec+p_unclear),max(p_spec+p_nonspec+p_unclear)],'k--')
                axs[i].set_xscale('log')
                axs[i].set_yscale('linear')
                axs[i].set_xlim([min(s_spec),max(s_spec)])
                axs[i].set_title(f'P-Value vs Selectivity Score for {prots[i]}',weight = 'bold',fontsize=20)
                axs[i].set_xlabel('Selectivity Score',weight = 'bold',fontsize=16)
                axs[i].set_ylabel('-Log10(P value)',weight = 'bold',fontsize=16)
                axs[i].legend(loc='best',fontsize=16)
                axs[i].tick_params(axis='both',which='major', labelsize=16)
            plt.tight_layout(rect=[0, 0.03, 1, 0.95])
            figname = 'P-Value vs Selectivity Score Internal Standardized.png'
            plt.savefig(os.path.join(save_dir,figname), bbox_inches="tight")
            plt.show()

############ Changed ####################

class Report_results:
    def __init__(self):
        pass
    
    def export_results(self,selectivity,sels_graphing,RTs_graphing,mzs_graphing,ps_graphing,z_graphing,sel_cutoff,
                       p_score_cutoff, prots,results_dir,misc_dir,feat_RT_cf_combined,feat_mz_filtered,
                       feat_z_filtered,pvals,spec_label,areas_savgol,full_out):
        '''
        Outputs calculated values into separate csv files for each protein
        '''
        for i in range(len(selectivity)):
            df = pd.DataFrame(columns=('Compound', 'm/z','z','p value','specificity'))
            score = sels_graphing[i]
            RT = RTs_graphing[i]
            mz = mzs_graphing[i]
            ps = ps_graphing[i]
            zs = z_graphing[i]
            score_sorted = np.sort(score)[::-1]            # sort the scores by descending order
            indices_sorted = np.argsort(score)[::-1]        # get indices of scores to get retention time and m/z
            count = 0
            while score_sorted[count] >= sel_cutoff:# and len(df) < 10000:
                index = indices_sorted[count]
                if ps[index] < p_score_cutoff:
                    A = f'Sel Score: {np.round(score_sorted[count],2)} RTime: {np.round(RT[index],3)} Areas: {[float(entry) for entry in areas_savgol[index]]}'
                    B = np.round(mz[index],4)
                    D = zs[index]
                    C = float(ps[index])
                    E = prots[i]
                    df.loc[count] = [A,B,D,C,E]
                count -= -1
                if count >= len(score):
                    break
            df.to_csv(os.path.join(misc_dir,prots[i] + '.csv'),index=False)

        if full_out:
            for i,score in enumerate(selectivity):
                df = pd.DataFrame(columns=('Compound', 'm/z','z','p value','specificity','areas'))
                score = selectivity[i]
                score_sorted = np.sort(score)[::-1]            # sort the scores by descending order
                indices_sorted = np.argsort(score)[::-1]        # get indices of scores to get retention time and m/z
                for j,score_sort in enumerate(score_sorted):
                    index = indices_sorted[j]
                    A = f'Sel Score: {np.round(score_sort,2)} RTime: {np.round(feat_RT_cf_combined[index],3)}'
                    B = np.round(feat_mz_filtered[index],4)
                    D = feat_z_filtered[index]
                    C = float(pvals[index])
                    E = spec_label[index]
                    F = [float(entry) for entry in areas_savgol[index]]
                    df.loc[j] = [A,B,D,C,E,F]
                df.to_csv(os.path.join(misc_dir,prots[i] + ' Full Output.csv'),index=False)

    def inclusion_lists(self,parent_dir,results_dir,misc_dir,folder,prots,gradient_time,peps_min,
                        RT_start=3,RT_window = 10,variable_loading = True):
        ''' 
        Generates inclusion lists based on the specified mass spectrometer capabilities, taking into account
        the maximum number of peptides that can be analyzed per minute over the given gradient. 
        Prioritizes high enrichment score features first.
        Adapted from Inclusion List Generator script, hence the need for importing files that were 
        just generated from the previous functions - need to rectify
        '''
        for i,prot in enumerate(prots):
            file = os.path.join(misc_dir,prot + '.csv')
            file_full = os.path.join(misc_dir,prot + ' Full Output.csv')

            data = []
            with open(file, 'r') as data_file:
                csv_reader = csv.DictReader(data_file,delimiter=',')
                for row in csv_reader:
                    data.append(row)
            RTs = []
            selscores = []
            mzs = []
            ps = []
            zs = []
            specs = []
            for line in data:            
                idx = line['Compound'].find('RTime:')
                idx2 = line['Compound'].find('Areas:')
                selscores.append(float(line['Compound'][11:idx-1]))
                RTs.append(float(line['Compound'][idx+7:idx2-1])/60)   # convert to mins
                mzs.append(float(line['m/z']))
                zs.append(float(line['z']))
                ps.append(float(line['p value']))
                specs.append(line['specificity'])

            RT_ranges = [(max(3,RT - 10),min(RT+10,gradient_time)) for RT in RTs]   # method is 120 mins, MS on at 3 mins

            dataframe = pd.DataFrame([selscores,RTs,mzs,zs,ps,specs]).transpose()
            dataframe.columns = ['SelScore','RT','mz','z','p','spec']
            bins = 120

            fig,axs = plt.subplots(1,2)
            axs[0].hist(RTs,bins=bins)
            axs[0].set_xlabel('Retention time')
            axs[0].set_ylabel('Frequency')
            axs[0].set_title('Retention Time Histogram')
            axs[1].hist(selscores,bins=100)
            axs[1].set_xlabel('Selectivity Score')
            axs[1].set_ylabel('Frequency')
            axs[1].set_title('Selectivity Score Histogram')
            plt.tight_layout()
            plt.savefig(os.path.join(misc_dir,f'{prot} Inclusion List Sel Score Histograms Start.png'), bbox_inches="tight")
            plt.show()

            n_bins = gradient_time - RT_start + 1
            hist,bins = np.histogram(RTs,bins=np.linspace(RT_start,gradient_time,n_bins),density=False)

            bin_members = np.digitize(RTs,bins=np.linspace(RT_start,gradient_time,n_bins))
            unique_members = np.unique(bin_members)
            idxs_filtered = []
            for member in unique_members:
                repeats = np.where(bin_members == member)[0]
                idxs = repeats.copy()
                scores_sub = [selscores[v] for v in repeats]
                if len(repeats) > peps_min:
                    scores_sub_sorted = np.sort(scores_sub)   # sort from lowest to highest
                    while len(scores_sub) > peps_min:
                        delete = np.where(scores_sub == scores_sub_sorted[0])[0][0]
                        scores_sub = np.delete(scores_sub,delete)
                        idxs = np.delete(idxs,delete)
                        scores_sub_sorted = np.delete(scores_sub_sorted,0)
                idxs_filtered.append(idxs)

            idxs_flattened = [val for sublist in idxs_filtered for val in sublist]
            RTs_filtered = [RTs[i] for i in idxs_flattened]
            selscores_filtered = [selscores[i] for i in idxs_flattened]
            mzs_filtered = [mzs[i] for i in idxs_flattened]
            zs_filtered = [zs[i] for i in idxs_flattened]
            ps_filtered = [ps[i] for i in idxs_flattened]
            t_starts = [RT - RT_window for RT in RTs_filtered]
            t_ends = [RT + RT_window for RT in RTs_filtered]                  

            fig, axs = plt.subplots(1,1)
            RT_ranges_filtered = [(max(3,RT - RT_window),min(RT+RT_window,120)) for RT in RTs_filtered]   # method is 120 mins, MS on at 3 mins
            for i,pair in enumerate(RT_ranges_filtered):
                axs.plot([pair[0],pair[1]],[mzs_filtered[i],mzs_filtered[i]])
                axs.scatter(RTs_filtered[i],mzs_filtered[i], s=12)
                axs.grid()
            axs.set_xlabel('Retention Time')
            axs.set_ylabel('m/z')
            axs.set_title('Peptides in Inclusion List \n Visualized by Estimated Elution Profile')
            plt.savefig(os.path.join(misc_dir,f'{prot} Sel Score Filtered Inclusion List.png'), bbox_inches="tight")
            plt.show()

            bins = 120
            fig,axs = plt.subplots(1,2)
            axs[0].hist(RTs_filtered,bins=bins)
            axs[0].set_xlabel('Retention time')
            axs[0].set_ylabel('Frequency')
            axs[0].set_title('Filtered Peptide List \n Retention Time Histogram')
            axs[1].set_xlabel('Selectivity Score')
            axs[1].set_ylabel('Frequency')
            axs[1].set_title('Filtered Peptide List \n Selectivity Score Histogram')
            axs[1].hist(selscores_filtered,bins=100)
            plt.tight_layout()
            plt.savefig(os.path.join(misc_dir,f'{prot} Sel Score Filtered Inclusion List Histogram.png'), bbox_inches="tight")

            df = pd.DataFrame(columns=('Compound', 'Formula','Adduct','m/z','z','t start (min)','t stop (min)'))
            for j,score in enumerate(selscores_filtered):
                A = f'Sel Score: {np.round(score,2)} p value: {"{:.3e}".format(ps_filtered[j])} RTime: {np.round(RTs_filtered[j],3)}'
                B = None
                C = '(no adduct)'
                D = np.round(mzs_filtered[j],4)
                E = int(zs_filtered[j])
                F = max(3,int(RTs_filtered[j] - RT_window))
                G = min(int(RTs_filtered[j] + RT_window),gradient_time)
                df.loc[j] = [A,B,C,D,E,F,G]
            if variable_loading:
                df.to_csv(os.path.join(misc_dir,f'{prot} Inclusion List Sel Score Filtered.csv'),index=False)
            else:
                df.to_csv(os.path.join(results_dir,f'{prot} Inclusion List Sel Score Filtered.csv'),index=False)

            df = pd.DataFrame(columns=('Compound', 'Formula','Adduct','m/z','z','t start (min)','t stop (min)'))
            for j,score in enumerate(selscores):
                A = 'Sel Score: ' + str(np.round(score,2)) + ' p value: ' + "{:.3e}".format(ps[j]) + ' RTime: ' + str(np.round(RTs[j],3))
                B = None
                C = '(no adduct)'
                D = np.round(mzs[j],4)
                E = int(zs[j])
                F = max(3,int(RTs[j] - 10))
                G = min(int(RTs[j] + 10),gradient_time)
                df.loc[j] = [A,B,C,D,E,F,G]
            df.to_csv(os.path.join(misc_dir,f'{prot} Sel Score Inclusion List Not Filtered.csv'),index=False)

    def inclusion_areascorefiltered(self,sels_graphing,RTs_graphing,RTs_graphing_orig,mzs_graphing,ps_graphing,z_graphing,
                                    areas_graphing,sel_cutoff,p_score_cutoff,prots,misc_dir,
                                    n_reps,area_limit,sel_minqual,gradient_time,variable_loading = True):
        ''' 
        Filters the given inclusion list based on a minimum peak area, then formats it for direct import into a 
        Thermo Xcalibur method.
        '''

        for i,score in enumerate(sels_graphing):
            df = pd.DataFrame(columns=('Compound', 'Formula','Adduct','m/z','z','t start (min)','t stop (min)'))
            score_sorted = np.sort(score)[::-1]            # sort the scores by descending order
            indices_sorted = np.argsort(score)[::-1]        # get indices of scores to get retention time and m/z
            for j,score_sort in enumerate(score_sorted):
                index = indices_sorted[j]
                time = RTs_graphing[i][index]/60
                area_prot = areas_graphing[i][index]
                ave_area = np.average(area_prot[n_reps*i:n_reps*i+n_reps])
                if ave_area > area_limit and score_sort > sel_minqual:
                    A = f'Sel Score: {np.round(score_sort,2)} RTime: {np.round(time,3)} Area: {ave_area} All areas: {area_prot} p-value: {ps_graphing[i][index][0]} original RTs: {RTs_graphing_orig[i][index]}'
                    B = None
                    C = '(no adduct)'
                    D = np.round(mzs_graphing[i][index],4)
                    E = z_graphing[i][index]
                    F = max(3,int(time - 10))
                    G = min(int(time + 10),gradient_time)
                    df.loc[j] = [A,B,C,D,E,F,G]
            df.to_csv(os.path.join(misc_dir,prots[i] + ' Inclusion List Area Filtered.csv'),index=False)

    def inclusion_areascorefiltered_excel(self,parent_dir,results_dir,misc_dir,prots,inj_level,
                                          num = 100,reps = 3,RT_start = 3, elute_window = 10):
        ''' 
        Exports all of the identified protein specific features into an Excel sheet, and additionally
        plots an EIC for the given top number of features. 
        '''
        savedir = os.path.join(misc_dir,"EICs Excel")   # Create placeholder directory for EIC pngs
        smallsavedir = os.path.join(savedir,"Small")
        if not os.path.exists(savedir):
            os.mkdir(savedir)
        if not os.path.exists(smallsavedir):
            os.makedirs(smallsavedir)
        all_files = []
        for file in sorted(os.listdir(parent_dir)):
            if file.endswith('.mzML') and inj_level in file:
                all_files.append(file)
        RTs_full = []
        ints_full = []
        mzs_full = []
        names = []
        for file in all_files:
            if file.endswith('.mzML'):
                print(file)
                names.append(file[:-5])
                exp = MSExperiment()
                MzMLFile().load(os.path.join(parent_dir,file),exp)
                RT_list_rep = []
                int_list_rep = []
                RTs = []
                mzs = []
                ints = []
                for spec in exp:
                    RTs.append(spec.getRT())
                    mzs_scan, ints_scan = spec.get_peaks()  
                    mzs.append(mzs_scan)
                    ints.append(ints_scan)
                RTs_full.append(RTs)
                ints_full.append(ints)
                mzs_full.append(mzs)  
        for i,prot in enumerate(prots):
            print(f'Working on {prot}')
            data_list = os.path.join(misc_dir,prots[i] + ' Inclusion List Sel Score Filtered.csv')
            target_file = pd.read_csv(data_list)
            target_file = target_file.fillna(' ')
            m_z_list = target_file['m/z'].tolist()
            tstart_list = target_file['t start (min)'].tolist()
            tend_list = target_file['t stop (min)'].tolist()
            feature_RT = [round(((start + stop) / 2)*60,2) for start, stop in zip(tstart_list, tend_list)]
            score_list = np.zeros(len(m_z_list))
            p_scores = np.zeros(len(m_z_list))
            time = elute_window*2
            peak_range = 30
            decimal = 2
            baseline = 1

            workbook = xlsxwriter.Workbook(os.path.join(misc_dir,f'EICs of {prot} by Sel Score Filter.xlsx'))
            print(os.path.join(results_dir,f'EICs of {prot}.xlsx'))
            worksheet = workbook.add_worksheet()
            for i,col in enumerate(target_file.columns):
                column_name = xlsxwriter.utility.xl_col_to_name(i)
                listing = target_file[col].tolist()
                for j,element in enumerate(listing):
                    if j == 0:
                        worksheet.write(column_name+str(j+1),target_file.columns[i])
                    worksheet.write(column_name+str(j+2),element)

            for j,mz in enumerate(m_z_list):       # loop through each feature
                RT_feature = feature_RT[j]
                s_score = score_list[j]
                p_val = p_scores[j]
                xs = []
                xs_zoom = []
                ys = []
                max_array = []
                for k in range(len(RTs_full)):             # loop through each replicate
                    RTs_rep = RTs_full[k]
                    RT_orig = RT_feature
                    ints_rep = ints_full[k]
                    mzs_rep = mzs_full[k]

                    RT_idx_low = Calculations().find_nearest(RTs_rep,RT_orig-time/2*60)[1]
                    RT_idx_high = Calculations().find_nearest(RTs_rep,RT_orig+time/2*60)[1] 

                    RTs_slice = RTs_rep[RT_idx_low:RT_idx_high+1]

                    max_idx_low = Calculations().find_nearest(RTs_slice,RT_orig-peak_range/2)[1]
                    max_idx_high = Calculations().find_nearest(RTs_slice,RT_orig+peak_range/2)[1]

                    ints_slice = ints_rep[RT_idx_low:RT_idx_high+1]    # should take a slice of arrays for both of these
                    mzs_slice = mzs_rep[RT_idx_low:RT_idx_high+1]
                    mz_idx = [Calculations().find_nearest_tol(entry,mz,10**(-decimal)/2)[1] for entry in mzs_slice] # find idx where the feature mz is in each scan

                    ints_EIC = [array[mz_idx[o]] if mz_idx[o] > 0 else baseline for o,array in enumerate(ints_slice)]

                    try:
                        max_found = np.amax(ints_EIC[max_idx_low:max_idx_high])
                    except ValueError:
                        max_found = np.amax(ints_EIC)
                        #print(ints_EIC[max_idx_low:max_idx_high],ints_EIC)
                    max_array.append(max_found)
                    xs.append(RTs_slice)
                    xs_zoom.append(RTs_slice[max_idx_low:max_idx_high])
                    ys.append(ints_EIC)

                #create the EIC as a small figure to be inserted into the Excel
                fig, axs = plt.subplots(1,int(len(names)),figsize=(2*len(names),2),sharey=True)
                for l,name in enumerate(names):
                    x = [t/60 for t in xs[l]]
                    y = ys[l]
                    axs[l].plot(x,y)
                    axs[l].set_title(f"{np.round(mz,decimal)} in {name}",fontsize=6)
                    axs[l].set_xlabel('Retention time (min)',fontsize=5)
                    axs[l].set_ylabel('Intensity',fontsize=5)
                    yfmt = mticker.ScalarFormatter(useMathText=True)
                    yfmt.set_powerlimits((3, 4))
                    axs[l].yaxis.set_major_formatter(yfmt)
                    axs[l].tick_params(axis='both', which='major', labelsize=7)
                    t1 = axs[l].yaxis.get_offset_text()
                    t1.set_x(-0.04)
                    t1.set_fontsize(5)
                plt.tight_layout()
                fignamesmall = f'EIC of mz of {np.round(mz,decimal)} at {np.round(RT_feature/60,1)} small.png'
                plt.savefig(os.path.join(smallsavedir,fignamesmall),dpi = 150, bbox_inches="tight")
                plt.close()

                # Insert into Excel, after making the row bigger
                worksheet.set_row(j+1,150)
                worksheet.insert_image('H'+str(j+2),os.path.join(smallsavedir,fignamesmall),{'x_offset': 2, 'y_offset': 2})  
                if j > num:
                    workbook.close()
                    break
            workbook.close() # catch in case masses reported < masses requested
        return savedir
    
    def delete_directory(self, directory):
        '''
        Removes the in situ generated EIC folder used to build the excel sheets. Use pathlib to avoid Mac
        resource forking issues while being compatible with Windows
        '''
        try:
            shutil.rmtree(directory)
            print(f'Successfully removed {directory}')
        except OSError as e:
            print(f'Error: {directory} - {e.strerror}')  
            
    def combined_result_export(self,feat_mz_filtered,feat_RT_filtered,feat_z_filtered,spec_label,pvals,selectivity,
                               selectivity_intstd,all_slope_abs,slopes_combined, 
                               areas_by_prot,yints_by_prot,ave_area_ref,prots,misc_dir):
        '''
        Takes both the fold enrichment and the selectivity score and exports into one file
        '''
        colors = ['#377eb8', '#ff7f00', '#4daf4a',
                          '#f781bf', '#a65628', '#984ea3',
                          '#999999', '#e41a1c', '#dede00']
        
        df_ranking = pd.DataFrame()
        df_ranking['m/z'] = feat_mz_filtered
        df_ranking['RT (min)'] = [a/60 for a in feat_RT_filtered]
        df_ranking['z'] = feat_z_filtered
        df_ranking['Obs Mass'] = [mz*z - z*1.00727647 for mz,z in zip(feat_mz_filtered,feat_z_filtered)]
        df_ranking['Specificity'] = spec_label
        df_ranking['P-value'] = [a[0] for a in pvals]        

        for i,prot in enumerate(prots):
            selscores = selectivity[i]
            selscores_intstd = selectivity_intstd[i]
            slope_abs = [entry[i] for entry in all_slope_abs]

            df_ranking[f'SelScore {prot}'] = selscores
            df_ranking[f'SelScore Intstd {prot}'] = selscores_intstd
            df_ranking[f'Ln(Concentration-Dependent Enrichment) {prot}'] = slope_abs     
      
        for i,prot in enumerate(prots):
            df_ranking[f'Areas_toberemoved_{prot}'] = areas_by_prot[prot]
            df_ranking[f'yint_toberemoved_{prot}'] = yints_by_prot[prot]
            df_ranking[f'Slope_toberemoved_{prot}'] = [a[i] for a in slopes_combined]
            
        df_ranking['Area of Beads Only Reference'] = [a[1] for a in ave_area_ref]
        df_ranking.to_csv(os.path.join(misc_dir,f'Combined Data.csv'),index = False)
        return df_ranking
    
    def combined_result_EICs(self,df_ranking,parent_dir,results_dir,misc_dir,prots,
                             num = 100, reps = 3,RT_start = 3, elute_window = 10, gradient_time = 120, intstd = False,
                            sortby_selscore = True):
        '''
        Takes both the fold enrichment and the selectivity score and exports into one file and also 
        includes plots for the fold enrichment and an EIC of the feature
        '''
        
        colors = ['#377eb8', '#ff7f00', '#4daf4a',
                          '#f781bf', '#a65628', '#984ea3',
                          '#999999', '#e41a1c', '#dede00']
        
        def plot_slopes(xvals,slope,y_int):
            return [x*slope + y_int for x in xvals]
        
        def custom_key(x):
            if x is None:
                return (float('inf'),)
            else:
                return (-x,)
        
        savedir = os.path.join(results_dir,"EICs Excel")   # Create placeholder directory for EIC pngs
        smallsavedir = os.path.join(savedir,"Small")
        if not os.path.exists(savedir):
            os.mkdir(savedir)
        if not os.path.exists(smallsavedir):
            os.makedirs(smallsavedir)
        all_files = []
        for file in sorted(os.listdir(parent_dir)):
            if file.endswith('.mzML') and any(sub in file for sub in prots):
                all_files.append(file)
        RTs_full = []
        ints_full = []
        mzs_full = []
        names = []
        for file in all_files:
            if file.endswith('.mzML'):
                print(file)
                names.append(file[:-5])
                exp = MSExperiment()
                MzMLFile().load(os.path.join(parent_dir,file),exp)
                RT_list_rep = []
                int_list_rep = []
                RTs = []
                mzs = []
                ints = []
                for spec in exp:
                    RTs.append(spec.getRT())
                    mzs_scan, ints_scan = spec.get_peaks()  
                    mzs.append(mzs_scan)
                    ints.append(ints_scan)
                RTs_full.append(RTs)
                ints_full.append(ints)
                mzs_full.append(mzs)  
        dfs_by_prot = []
        for i,prot in enumerate(prots):
            df = df_ranking[df_ranking['Specificity'] == prot]
            df = df[df[f'Ln(Concentration-Dependent Enrichment) {prot}'] > 0]                    
            print(f'Working on {prot}')
            
            selscores = df[f'SelScore {prot}']
            selscores_intstd = df[f'SelScore Intstd {prot}']
            slope_abs = df[f'Ln(Concentration-Dependent Enrichment) {prot}']
            selscore_sort = sorted(enumerate(selscores), key=lambda x: custom_key(x[1]))
            selscore_intstd_sort = sorted(enumerate(selscores_intstd), key=lambda x: custom_key(x[1]))
            slope_sort = sorted(enumerate(slope_abs), key=lambda x: custom_key(x[1]))

            selscore_ranks = sorted(list(enumerate([a[0] for a in selscore_sort])), key = lambda x: x[1])
            selscore_intstd_ranks = sorted(list(enumerate([a[0] for a in selscore_intstd_sort])), key = lambda x: x[1])
            slope_ranks = sorted(list(enumerate([a[0] for a in slope_sort])), key = lambda x: x[1])
            df[f'SelScore Rank {prot}'] = [a[0] for a in selscore_ranks]
            df[f'SelScore Intstd Rank {prot}'] = [a[0] for a in selscore_intstd_ranks]
            df[f'Concentration-Dependent Enrichment Rank {prot}'] = [a[0] for a in slope_ranks]

            rankme = [e + s for e,s in zip([a[0] for a in selscore_ranks],[a[0] for a in slope_ranks])]                  # Use average placement to determine ranking
            rankme_norm = [e + s for e,s in zip([a[0] for a in selscore_intstd_ranks],[a[0] for a in slope_ranks])]
            df[f'Total Rank {prot}'] = rankme
            df[f'Total Intstd Rank {prot}'] = rankme_norm 
            if intstd and sortby_selscore:
                df = df.sort_values(by = f'SelScore Intstd Rank {prot}')
                workbook_name = os.path.join(results_dir,f'Final Results for {prot} with EICs Int Std.xlsx')
            if intstd and not sortby_selscore:
                df = df.sort_values(by = f'Total Intstd Rank {prot}')
                workbook_name = os.path.join(results_dir,f'Final Results for {prot} with EICs Int Std.xlsx')
            if not intstd and sortby_selscore:
                df = df.sort_values(by = f'SelScore Rank {prot}')
                workbook_name = os.path.join(results_dir,f'Final Results for {prot} with EICs.xlsx')
            else:
                df = df.sort_values(by = f'Total Rank {prot}')
                workbook_name = os.path.join(results_dir,f'Final Results for {prot} with EICs.xlsx')
            
            df.to_csv(os.path.join(misc_dir,f'Filtered Results {prot}.csv'), index = False)
            dfs_by_prot.append(df)

            workbook = xlsxwriter.Workbook(workbook_name)
            print(os.path.join(results_dir,f'EICs of {prot}.xlsx'))
            worksheet = workbook.add_worksheet()
            columns_report = [c for c in df.columns.to_list() if '_toberemoved_' not in c]
            for i,col in enumerate(columns_report):    
                column_name = xlsxwriter.utility.xl_col_to_name(i)
                listing = df[col].tolist()
                for j,element in enumerate(listing):
                    if j == 0:
                        worksheet.write(column_name+str(j+1),columns_report[i])
                    worksheet.write(column_name+str(j+2),element)
                    
            m_z_list = df['m/z'].to_list()
            RTs = df['RT (min)'].to_list()
            tstart_list = [max(RT_start, RT-elute_window) for RT in RTs]
            tend_list = [min(gradient_time, RT+elute_window) for RT in RTs]
            feature_RT = [round(((start + stop) / 2)*60,2) for start, stop in zip(tstart_list, tend_list)]
            score_list = np.zeros(len(m_z_list))
            p_scores = np.zeros(len(m_z_list))
            time = elute_window*2
            peak_range = 30
            decimal = 2
            baseline = 1

            for j,mz in enumerate(m_z_list):       # loop through each feature
                ## Insert slope plot
                fig, axs = plt.subplots(1,1,figsize=(4,3),sharey=True)
                for i,prot in enumerate(prots):
                    area_prot = df[f'Areas_toberemoved_{prot}'].to_list()[j]
                    yint_prot = df[f'yint_toberemoved_{prot}'].to_list()[j]
                    slope_prot = df[f'Slope_toberemoved_{prot}'].to_list()[j]
                    x_values = [a[0] for a in area_prot]
                    y_values = [a[1] for a in area_prot]
                    plt.scatter(x_values, y_values, color = colors[i], marker='o', label=prot, alpha = 0.6) 
                    x_values.append(0)
                    y_fit = plot_slopes(x_values,slope_prot,yint_prot)
                    plt.plot(x_values,y_fit, color = colors[i], linestyle = '--', label = f'{prot} fit')

                plt.scatter(0,df['Area of Beads Only Reference'].to_list()[j], 
                            color = 'k', marker = 'x', label = 'Beads Only')
                axs.set_title(f"{np.round(mz,decimal)} Areas vs Loading",fontsize=10, weight = 'bold')
                axs.set_xlabel('Percentage Loading',fontsize=8, weight = 'bold')
                axs.set_ylabel('Area',fontsize=8, weight = 'bold')
                axs.legend(loc = 'best', frameon = True)
                yfmt = mticker.ScalarFormatter(useMathText=True)
                yfmt.set_powerlimits((3, 4))
                axs.yaxis.set_major_formatter(yfmt)
                axs.tick_params(axis='both', which='major', labelsize=7)
                t1 = axs.yaxis.get_offset_text()
                t1.set_x(-0.04)
                t1.set_fontsize(8)
                
                plt.tight_layout()
                fignamesmall = f'Enrichment Plot of {np.round(mz,decimal)} small.png'
                plt.savefig(os.path.join(savedir,fignamesmall),dpi = 150, bbox_inches="tight")
                plt.close()
                column_name = xlsxwriter.utility.xl_col_to_name(len(columns_report))
                worksheet.insert_image(column_name+str(j+2),os.path.join(savedir,fignamesmall),{'x_offset': 2, 'y_offset': 2})                                 
                RT_feature = feature_RT[j]
                s_score = score_list[j]
                p_val = p_scores[j]
                xs = []
                xs_zoom = []
                ys = []
                max_array = []
                for k in range(len(RTs_full)):             # loop through each replicate
                    RTs_rep = RTs_full[k]
                    RT_orig = RT_feature
                    ints_rep = ints_full[k]
                    mzs_rep = mzs_full[k]

                    RT_idx_low = Calculations().find_nearest(RTs_rep,RT_orig-time/2*60)[1]
                    RT_idx_high = Calculations().find_nearest(RTs_rep,RT_orig+time/2*60)[1] 

                    RTs_slice = RTs_rep[RT_idx_low:RT_idx_high+1]

                    max_idx_low = Calculations().find_nearest(RTs_slice,RT_orig-peak_range/2)[1]
                    max_idx_high = Calculations().find_nearest(RTs_slice,RT_orig+peak_range/2)[1]

                    ints_slice = ints_rep[RT_idx_low:RT_idx_high+1]    # should take a slice of arrays for both of these
                    mzs_slice = mzs_rep[RT_idx_low:RT_idx_high+1]
                    mz_idx = [Calculations().find_nearest_tol(entry,mz,10**(-decimal)/2)[1] for entry in mzs_slice] # find idx where the feature mz is in each scan

                    ints_EIC = [array[mz_idx[o]] if mz_idx[o] > 0 else baseline for o,array in enumerate(ints_slice)]

                    try:
                        max_found = np.amax(ints_EIC[max_idx_low:max_idx_high])
                    except ValueError:
                        max_found = np.amax(ints_EIC)
                        #print(ints_EIC[max_idx_low:max_idx_high],ints_EIC)
                    max_array.append(max_found)
                    xs.append(RTs_slice)
                    xs_zoom.append(RTs_slice[max_idx_low:max_idx_high])
                    ys.append(ints_EIC)

                #create the EIC as a small figure to be inserted into the Excel
                fig, axs = plt.subplots(1,int(len(names)),figsize=(2*len(names),2),sharey=True)
                for l,name in enumerate(names):
                    x = [t/60 for t in xs[l]]
                    y = ys[l]
                    axs[l].plot(x,y)
                    axs[l].set_title(f"{np.round(mz,decimal)} in {name}",fontsize=6)
                    axs[l].set_xlabel('Retention time (min)',fontsize=5)
                    axs[l].set_ylabel('Intensity',fontsize=5)
                    yfmt = mticker.ScalarFormatter(useMathText=True)
                    yfmt.set_powerlimits((3, 4))
                    axs[l].yaxis.set_major_formatter(yfmt)
                    axs[l].tick_params(axis='both', which='major', labelsize=7)
                    t1 = axs[l].yaxis.get_offset_text()
                    t1.set_x(-0.04)
                    t1.set_fontsize(5)
                plt.tight_layout()
                fignamesmall = f'EIC of mz of {np.round(mz,decimal)} at {np.round(RT_feature/60,1)} small.png'
                plt.savefig(os.path.join(smallsavedir,fignamesmall),dpi = 150, bbox_inches="tight")
                plt.close()

                # Insert into Excel, after making the row bigger
                worksheet.set_row(j+1,225)
                column_name = xlsxwriter.utility.xl_col_to_name(len(columns_report) + 6)
                worksheet.insert_image(column_name+str(j+2),os.path.join(smallsavedir,fignamesmall),{'x_offset': 2, 'y_offset': 2})  
                if j > num:
                    workbook.close()
                    break
            workbook.close() # catch in case masses reported < masses requested
        return savedir, dfs_by_prot
    
    def final_inclusion_lists(self,df,parent_dir,results_dir,misc_dir,folder,prot,gradient_time,peps_min,
                              RT_start=3,RT_window = 10, intstd = False, sortby_selscore = False):
        ''' 
        Generates inclusion lists based on the specified mass spectrometer capabilities, taking into account
        the maximum number of peptides that can be analyzed per minute over the given gradient. 
        Prioritizes high enrichment score features first.
        Adapted from Inclusion List Generator script, hence the need for importing files that were 
        just generated from the previous functions - need to rectify
        '''
        if intstd:
            if sortby_selscore:
                Ranks = df[f'SelScore Intstd Rank {prot}']
            else:
                Ranks = df[f'Total Intstd Rank {prot}'].to_list() 
            selscores = df[f'SelScore Intstd {prot}'].to_list()
            Ln_CDE = df[f'Ln(Concentration-Dependent Enrichment) {prot}'].to_list()
            histstart = os.path.join(misc_dir,f'{prot} Inclusion List Intstd Histograms Start.png')
            plotfilt = os.path.join(results_dir,f'{prot} Final Filtered Inclusion List Intstd.png')
            histend = os.path.join(results_dir,f'{prot} Final Filtered Inclusion List Intstd Histogram.png')
            finallist = os.path.join(results_dir,f'{prot} Final Inclusion List Intstd Filtered.csv')
            finallist_notfilt = os.path.join(misc_dir,f'{prot} Final Inclusion List Intstd Not Filtered.csv')
        else:
            if sortby_selscore:
                Ranks = df[f'SelScore Rank {prot}']
            else:
                Ranks = df[f'Total Rank {prot}'].to_list()
            selscores = df[f'SelScore {prot}'].to_list()
            Ln_CDE = df[f'Ln(Concentration-Dependent Enrichment) {prot}'].to_list()
            histstart = os.path.join(misc_dir,f'{prot} Inclusion List Histograms Start.png')
            plotfilt = os.path.join(results_dir,f'{prot} Final Filtered Inclusion List.png')
            histend = os.path.join(results_dir,f'{prot} Final Filtered Inclusion List Histogram.png')
            finallist = os.path.join(results_dir,f'{prot} Final Inclusion List Filtered.csv')
            finallist_notfilt = os.path.join(misc_dir,f'{prot} Final Inclusion List Not Filtered.csv')
        
        RTs = df['RT (min)'].to_list()
        mzs = df['m/z'].to_list()
        ps = df['P-value'].to_list()
        zs = df['z'].to_list()
        specs = df['Specificity'].to_list()
        slopes = df[f'Ln(Concentration-Dependent Enrichment) {prot}'].to_list()

        RT_ranges = [(max(3,RT - RT_window),min(RT+RT_window,gradient_time)) for RT in RTs]   # method is 120 mins, MS on at 3 mins

        dataframe = pd.DataFrame([Ranks,RTs,mzs,zs,ps,specs]).transpose()
        dataframe.columns = ['Ranks','RT','mz','z','p','spec']
        bins = 120

        fig,axs = plt.subplots(1,3, figsize = (10,4))
        axs[0].hist(RTs,bins=bins)
        axs[0].set_xlabel('Retention time')
        axs[0].set_ylabel('Frequency')
        axs[0].set_title(f'{prot} Retention Time Histogram')
        axs[1].hist(selscores,bins=100)
        axs[1].set_xlabel('Selectivity Score')
        axs[1].set_ylabel('Frequency')
        axs[1].set_title(f'{prot} Selectivity Score Histogram')
        axs[2].hist(Ln_CDE,bins=100)
        axs[2].set_xlabel('Ln(Enrichment)')
        axs[2].set_ylabel('Frequency')
        axs[2].set_title(f'{prot} Concentration-Dependent \n Enrichment Histogram')
        plt.tight_layout()
        plt.savefig(histstart, bbox_inches="tight")
        plt.show()

        n_bins = gradient_time - RT_start + 1
        hist,bins = np.histogram(RTs,bins=np.linspace(RT_start,gradient_time,n_bins),density=False)

        bin_members = np.digitize(RTs,bins=np.linspace(RT_start,gradient_time,n_bins))
        unique_members = np.unique(bin_members)
        idxs_filtered = []
        for member in unique_members:
            repeats = np.where(bin_members == member)[0]
            idxs = repeats.copy()
            scores_sub = [Ranks[v] for v in repeats]
            if len(repeats) > peps_min:
                scores_sub_sorted = np.sort(scores_sub)   # sort from lowest to highest
                while len(scores_sub) > peps_min:
                    delete = np.where(scores_sub == scores_sub_sorted[0])[0][0]
                    scores_sub = np.delete(scores_sub,delete)
                    idxs = np.delete(idxs,delete)
                    scores_sub_sorted = np.delete(scores_sub_sorted,0)
            idxs_filtered.append(idxs)

        idxs_flattened = [val for sublist in idxs_filtered for val in sublist]
        RTs_filtered = [RTs[i] for i in idxs_flattened]
        Ranks_filtered = [Ranks[i] for i in idxs_flattened]
        selscores_filtered = [selscores[i] for i in idxs_flattened]
        Ln_CDEs_filtered = [Ln_CDE[i] for i in idxs_flattened]
        mzs_filtered = [mzs[i] for i in idxs_flattened]
        zs_filtered = [zs[i] for i in idxs_flattened]
        ps_filtered = [ps[i] for i in idxs_flattened]
        slopes_filtered = [slopes[i] for i in idxs_flattened]
        t_starts = [RT - RT_window for RT in RTs_filtered]
        t_ends = [RT + RT_window for RT in RTs_filtered]

        fig, axs = plt.subplots(1,1)
        RT_ranges_filtered = [(max(3,RT - RT_window),min(RT+RT_window,120)) for RT in RTs_filtered]   # method is 120 mins, MS on at 3 mins
        for i,pair in enumerate(RT_ranges_filtered):
            axs.plot([pair[0],pair[1]],[mzs_filtered[i],mzs_filtered[i]])
            axs.scatter(RTs_filtered[i],mzs_filtered[i], s=12)
            axs.grid()
        axs.set_xlabel('Retention Time')
        axs.set_ylabel('m/z')
        axs.set_title(f'{prot} Peptides in Inclusion List \n Visualized by Estimated Elution Profile')
        plt.savefig(plotfilt, bbox_inches="tight")
        plt.show()

        bins = 120
        fig,axs = plt.subplots(1,3,figsize = (10,4))
        axs[0].hist(RTs_filtered,bins=bins)
        axs[0].set_xlabel('Retention time')
        axs[0].set_ylabel('Frequency')
        axs[0].set_title(f'{prot} Filtered Peptide List \n Retention Time Histogram')
        axs[1].set_xlabel('Selectivity Score')
        axs[1].set_ylabel('Frequency')
        axs[1].set_title(f'{prot} Filtered Peptide List \n Selectivity Score Histogram')
        axs[1].hist(selscores_filtered,bins=100)
        axs[2].set_xlabel('Ln(Enrichment)')
        axs[2].set_ylabel('Frequency')
        axs[2].set_title(f'{prot} Filtered Peptide List \n Concentration-Dependent Enrichment Histogram')
        axs[2].hist(Ln_CDEs_filtered,bins=100)
        plt.tight_layout()
        plt.savefig(histend, bbox_inches="tight")

        df = pd.DataFrame(columns=('Compound', 'Formula','Adduct','m/z','z','t start (min)','t stop (min)'))
        for j,score in enumerate(selscores_filtered):
            A = f'Overall rank: {Ranks_filtered[j]} Sel Score: {np.round(score,2)} CDE: {np.round(slopes_filtered[j],4)} p value: {"{:.3e}".format(ps_filtered[j])} RTime: {np.round(RTs_filtered[j],3)}'
            B = None
            C = '(no adduct)'
            D = np.round(mzs_filtered[j],4)
            E = int(zs_filtered[j])
            F = max(3,int(RTs_filtered[j] - RT_window))
            G = min(int(RTs_filtered[j] + RT_window),gradient_time)
            df.loc[j] = [A,B,C,D,E,F,G]
        df.to_csv(finallist,index=False)

        df = pd.DataFrame(columns=('Compound', 'Formula','Adduct','m/z','z','t start (min)','t stop (min)'))
        for j,score in enumerate(selscores):
            A = f'Overall rank: {Ranks[j]} Sel Score: {np.round(score,2)} CDE: {np.round(slopes[j],4)} p value: {"{:.3e}".format(ps[j])} RTime: {np.round(RTs[j],3)}'
            B = None
            C = '(no adduct)'
            D = np.round(mzs[j],4)
            E = int(zs[j])
            F = max(3,int(RTs[j] - 10))
            G = min(int(RTs[j] + 10),gradient_time)
            df.loc[j] = [A,B,C,D,E,F,G]
        df.to_csv(finallist_notfilt,index=False)

class Conc_Dep_Enrichment:
    def __init__(self):
        pass
    
    def slice_list_by_indices(self,data_list, indices):
        return [data_list[i] for i in indices]

    def calculate_slopes(self,coordinate_sets,ref_points):
        '''
        Fits the points to a line with the y-intercept set by the beads only sample
        '''
        def line(x, slope):
            return slope * x

        slopes = []
        yints = []
        for i,coordinates in enumerate(coordinate_sets):
            # Extract x and y values from the coordinate set
            x_values = [coord[0] for coord in coordinates]
            x_values.insert(-1,ref_points[i][0])
            y_values = [coord[1] for coord in coordinates]
            y_values.insert(-1,ref_points[i][1])
            y_values_fit = [a - ref_points[i][1] for a in y_values]  # set y-intercept to 0 for fitting 
            # Fit a line to the data points
            popt, pcov = curve_fit(line, x_values, y_values_fit)
            slopes.append(popt[0])
            yints.append(ref_points[i][1])
        return slopes, yints
    
    def setup_CDE(self,areas_savgol_CDE,files_all,prots,inj_levels,inj_percs,area_ref,
                 num = 100,decimal = 4, refname='BeadsOnly'):
        '''
        Calculate the slopes for all of the specified features using the areas
        '''
        areas_per_file = {}
        for i,file in enumerate(files_all):
            areas_sub = [entry[i] for entry in areas_savgol_CDE]
            areas_per_file[file] = areas_sub

        combinations = list(product(prots, inj_levels))
        combo_locs = {}

        slope_list = []
        for combo in combinations:
            for i,inj in enumerate(inj_levels):
                if inj in combo:
                    slope_list.append(inj_percs[i])

        for combo in combinations:
            combo_idx = []
            for i,key in enumerate(files_all):
                if all(substr in key for substr in combo):
                    combo_idx.append(i)
            combo_locs[combo] = combo_idx

        area_prot_levels = {}
        for i,(key,value) in enumerate(combo_locs.items()):
            areas_per_combo = [self.slice_list_by_indices(entry,value) for entry in areas_savgol_CDE]
            areas_per_combo_ave = [(slope_list[i],np.mean(a)) for a in areas_per_combo]
            area_prot_levels[key] = areas_per_combo_ave

        areas_by_prot = {}
        for prot in prots:
            areas_per_prot = []
            for key,value in area_prot_levels.items():
                if prot in key:
                    areas_per_prot.append(value)
            areas_by_prot[prot] = list(zip(*areas_per_prot))
            
        areas_ref = {key: value for key, value in areas_per_file.items() if refname in key}
        ave_area_ref = [(0,sum(values)/len(values)) for values in zip(*areas_ref.values())]
            
        slopes_by_prot = {}
        yints_by_prot = {}
        for i,(key,value) in enumerate(areas_by_prot.items()):
            slopes, yints = self.calculate_slopes(value,ave_area_ref)
            slopes_by_prot[key] = slopes
            yints_by_prot[key] = yints

        slopes_combined = list(zip(*list(slopes_by_prot.values())))
        return areas_by_prot, yints_by_prot, slopes_combined, ave_area_ref

    def CDE_report_excel(self,areas_by_prot,yints_by_prot,ave_area_ref,slopes_combined,feat_mz_filtered,feat_RT_filtered,
                        feat_z_filtered,results_dir,misc_dir,prots,num = 100,decimal = 4):
        ''' 
        Exports all of the identified protein specific features into an Excel sheet, and additionally
        plots an EIC for the given top number of features. 
        '''
        def plot_slopes(xvals,slope,y_int):
            return [x*slope + y_int for x in xvals]
        
        colors = ['#377eb8', '#ff7f00', '#4daf4a',
                          '#f781bf', '#a65628', '#984ea3',
                          '#999999', '#e41a1c', '#dede00']
        savedir = os.path.join(misc_dir,"CDE Plots Excel")   # Create placeholder directory for EIC pngs
        if not os.path.exists(savedir):
            os.mkdir(savedir)

        fig_refs = {}
            
        all_slope_abs = []
        for l,slope_pair in enumerate(slopes_combined):
            slope_abs = []
            for p,prot in enumerate(prots):
                slope_spec = slope_pair[p]
                if slope_spec > 0:
                    slope_abs.append(np.log(slope_spec))
                elif slope_spec < 0: 
                    slope_abs.append(-np.log(np.abs(slope_spec)))
                else:
                    slope_abs.append(0)
            all_slope_abs.append(slope_abs)

        for p,prot in enumerate(prots):
            workbook = xlsxwriter.Workbook(os.path.join(misc_dir,f'CDE Visualization of {prot}.xlsx'))
            print(os.path.join(results_dir,f'CDE Visualization of {prot}.xlsx'))
            worksheet = workbook.add_worksheet()
            columns = ['m/z','z','RT','Obs Mass']
            for name in prots:
                columns.append(f'CDE Slope with respect to {name}')
            columns.append('Ln(Slope Diff) Prot of Interest vs All Others')
            columns.append('Slope Plot')
            for i,col in enumerate(columns):
                column_name = xlsxwriter.utility.xl_col_to_name(i)
                worksheet.write(column_name+str(1),columns[i])

            indexed_slopes = [(index, value) for index, value in enumerate(all_slope_abs)]
            sorted_slopes = sorted(indexed_slopes, key=lambda x: x[1][p], reverse = True)
            original_indices = [indexed_tuple[0] for indexed_tuple in sorted_slopes]
            sorted_slopes = [indexed_tuple[1] for indexed_tuple in sorted_slopes]

            ## Write all rows, very redundant but whatever
            for l,idx in enumerate(original_indices):
                mz = feat_mz_filtered[idx]
                RT = feat_RT_filtered[idx]
                z = feat_z_filtered[idx]
                obs_mass = mz*z - z*1.0073
                columns_write = columns[:-1]
                data = [mz,z,RT,obs_mass]
                for slope in slopes_combined[idx]:
                    data.append(slope)
                data.append(all_slope_abs[idx][p])
                for i, col in enumerate(columns_write):
                    column_name = xlsxwriter.utility.xl_col_to_name(i)
                    worksheet.write(column_name+str(l+2),data[i])

            ## Make slope plots for easy visualization
            for l,idx in enumerate(original_indices):
                mz = feat_mz_filtered[idx]
                fig, axs = plt.subplots(1,1,figsize=(4,3),sharey=True)
                slope_spec = slopes_combined[idx]
                slope_other = sum(slope_spec[:idx] + slope_spec[idx + 1:])
                for i,(prot,val) in enumerate(areas_by_prot.items()):
                    coordinates = val[idx]
                    x_values, y_values = zip(*coordinates)                  
                    plt.scatter(x_values, y_values, color = colors[i], marker='o', label=prot, alpha = 0.6) 
                    x_values = x_values + (0,)
                    y_fit = plot_slopes(x_values,slope_spec[i],yints_by_prot[prot][idx])
                    plt.plot(x_values,y_fit, color = colors[i], linestyle = '--', label = f'{prot} fit')
                    
                plt.scatter(ave_area_ref[idx][0],ave_area_ref[idx][1], color = 'k', marker = 'x', label = 'Beads Only')
                axs.set_title(f"{np.round(mz,decimal)} Areas vs Loading",fontsize=10, weight = 'bold')
                axs.set_xlabel('Percentage Loading',fontsize=8, weight = 'bold')
                axs.set_ylabel('Area',fontsize=8, weight = 'bold')
#                 axs.set_yscale('log')
                axs.legend(loc = 'best', frameon = True)
                yfmt = mticker.ScalarFormatter(useMathText=True)
                yfmt.set_powerlimits((3, 4))
                axs.yaxis.set_major_formatter(yfmt)
                axs.tick_params(axis='both', which='major', labelsize=7)
                t1 = axs.yaxis.get_offset_text()
                t1.set_x(-0.04)
                t1.set_fontsize(8)
                
                plt.tight_layout()
                fignamesmall = f'CDE Plot of {np.round(mz,decimal)} small.png'
                plt.savefig(os.path.join(savedir,fignamesmall),dpi = 150, bbox_inches="tight")
                plt.close()
                
                fig_refs[str(mz)] = os.path.join(savedir,fignamesmall)

                # Insert into Excel, after making the row bigger
                worksheet.set_row(l+1,225)
                worksheet.insert_image(xlsxwriter.utility.xl_col_to_name(len(columns)-1)+str(l+2),os.path.join(savedir,fignamesmall),{'x_offset': 2, 'y_offset': 2})  
                if l > num:
                    workbook.close()
                    break
            workbook.close()
        return savedir, all_slope_abs, fig_refs
    
    def plot_cde_ranking(self,prots,all_slope_abs,results_dir):
        '''
        Plots the concentration-dependent enrichment based on the ranking from highest to lowest
        '''
        if len(prots) > 2:
            fig,axs = plt.subplots(math.ceil(len(prots)/2),2,figsize=(14,14),sharey=True)
            for i in range(math.ceil(len(prots)/2)):
                for j in range(2):
                    try:
                        slope_abs = [entry[2*i+j] for entry in all_slope_abs]
                    except:
                        continue
                    sorted_slopes = np.sort(slope_abs)[::-1]
                    axs[i,j].plot(sorted_slopes[:1000],label = f'{prots[2*i+j]}')
                    axs[i,j].set_yscale('linear')
                    axs[i,j].set_title(f'CDE Ranking of {prots[2*i+j]}', weight = 'bold', fontsize = 18)
                    axs[i,j].set_xlabel('Rank', weight = 'bold', fontsize = 14)
                    axs[i,j].set_ylabel('Ln(CDE)', weight = 'bold', fontsize = 14)
                    axs[i,j].grid(True, which='major', color='dimgray', linestyle='-')
                    axs[i,j].grid(True, which='minor', color='lightgray', linestyle='--')
                    axs[i,j].legend(loc='upper right')
            if len(prots) % 2 != 0:
                fig.delaxes(axs[-1][-1])
            plt.tight_layout()
            figname = os.path.join(results_dir,'Concentration-Dependent Enrichment Plots Zoomed.png')
            plt.savefig(figname, bbox_inches="tight")

            fig,axs = plt.subplots(math.ceil(len(prots)/2),2,figsize=(14,10),sharey=True)
            for i in range(math.ceil(len(prots)/2)):
                for j in range(2):
                    try:
                        slope_abs = [entry[2*i+j] for entry in all_slope_abs]
                    except:
                        continue
                    sorted_slopes = np.sort(slope_abs)[::-1]
                    axs[i,j].plot(sorted_slopes,label = f'{prots[2*i+j]}')
                    axs[i,j].set_yscale('linear')
                    axs[i,j].set_title(f'CDE Ranking of {prots[2*i+j]}', weight = 'bold', fontsize = 18)
                    axs[i,j].set_xlabel('Rank', weight = 'bold', fontsize = 14)
                    axs[i,j].set_ylabel('Ln(CDE Enrichment)', weight = 'bold', fontsize = 14)
                    axs[i,j].grid(True, which='major', color='dimgray', linestyle='-')
                    axs[i,j].grid(True, which='minor', color='lightgray', linestyle='--')
                    axs[i,j].legend(loc='upper right')
            if len(prots) % 2 != 0:
                fig.delaxes(axs[-1][-1])
            plt.tight_layout()
            figname = os.path.join(results_dir,'Concentration-Dependent Enrichment Plots.png')
            plt.savefig(figname, bbox_inches="tight") 

        else:
            fig,axs = plt.subplots(math.ceil(len(prots)/2),2,figsize=(14,10),sharey=True)
            for i,prot in enumerate(prots):
                slope_abs = [entry[i] for entry in all_slope_abs]
                sorted_slopes = np.sort(slope_abs)[::-1]
                axs[i].plot(sorted_slopes[:1000],label = f'{prot}')
                axs[i].set_yscale('linear')
                axs[i].set_title(f'CDE Ranking of {prot}', weight = 'bold', fontsize = 18)
                axs[i].set_xlabel('Rank', weight = 'bold', fontsize = 14)
                axs[i].set_ylabel('Ln(CDE)', weight = 'bold', fontsize = 14)
                axs[i].grid(True, which='major', color='dimgray', linestyle='-')
                axs[i].grid(True, which='minor', color='lightgray', linestyle='--')
                axs[i].legend(loc='upper right')
            plt.tight_layout()
            figname = os.path.join(results_dir,'Concentration-Dependent Enrichment Plots Zoomed.png')
            plt.savefig(figname, bbox_inches="tight")

            fig,axs = plt.subplots(math.ceil(len(prots)/2),2,figsize=(14,10),sharey=True)
            for i,prot in enumerate(prots):
                slope_abs = [entry[i] for entry in all_slope_abs]
                sorted_slopes = np.sort(slope_abs)[::-1]
                axs[i].plot(sorted_slopes,color = 'k', label = f'{prot}')
                axs[i].set_yscale('linear')
                axs[i].set_title(f'CDE Ranking of {prot}', weight = 'bold', fontsize = 18)
                axs[i].set_xlabel('Rank', weight = 'bold', fontsize = 14)
                axs[i].set_ylabel('Ln(CDE)', weight = 'bold', fontsize = 14)
                axs[i].grid(True, which='major', color='dimgray', linestyle='-')
                axs[i].grid(True, which='minor', color='lightgray', linestyle='--')
                axs[i].legend(loc='upper right')
            plt.tight_layout()
            figname = os.path.join(results_dir,'Concentration-Dependent Enrichment Plots.png')
            plt.savefig(figname, bbox_inches="tight")
            
    def selscore_vs_cde(self,prots,selectivity,selectivity_intstd,all_slope_abs,results_dir, intstd = False,
                    sel_cutoff = 0.5):
        '''
        Plots the selectivity score versus the concentration-dependent enrichment for all of the features
        '''
        for i,prot in enumerate(prots):
            fig,axs = plt.subplots(1,1,figsize = (6,5))
            print(prot)
            selscores = selectivity[i]
            selscores_intstd = selectivity_intstd[i]
            slope_abs = [entry[i] for entry in all_slope_abs]
            if intstd:
                for e,s in zip(selscores_intstd,slope_abs):
                    if e > sel_cutoff and s > 0:
                        axs.scatter(s,e,color = 'b')
                    else:
                        axs.scatter(s,e,color = 'r')
                figname = os.path.join(results_dir,f'Selectivity vs Concentration-Dependent Enrichment {prot} Int Standardized.png')
                axs.set_title(f'CDE vs Selectivity Score Int Std of {prot}', weight = 'bold', fontsize = 16)
            else:
                for e,s in zip(selscores,slope_abs):
                    if e > sel_cutoff and s > 0:
                        axs.scatter(s,e,color = 'b')
                    else:
                        axs.scatter(s,e,color = 'r')
                figname = os.path.join(results_dir,f'Selectivity vs Concentration-Dependent Enrichment {prot}.png')
                axs.set_title(f'CDE vs Selectivity Score of {prot}', weight = 'bold', fontsize = 16)
            axs.plot([min(slope_abs),max(slope_abs)],[sel_cutoff,sel_cutoff],
                      color = 'gray', linestyle = '--', label = 'Selectivity score cutoff')
            axs.set_xlabel('Ln(CDE)', weight = 'bold', fontsize = 14)
            axs.set_ylabel('Selectivity Score', weight = 'bold', fontsize = 14)
            axs.grid(True, which='major', color='dimgray', linestyle='-')
            axs.grid(True, which='minor', color='lightgray', linestyle='--')
            axs.legend(loc='best')
            plt.tight_layout()           
            plt.savefig(figname, bbox_inches="tight") 
            
    def cde_vs_pval(self,prots,all_slope_abs,pvals,results_dir,p_score_cutoff = 0.05):
        '''
        Plots the concentration-dependent enrichment versus the p-value for all of the features, similar to the volcano plot
        '''
        for i,prot in enumerate(prots):
            fig,axs = plt.subplots(1,1,figsize = (6,5))
            print(prot)
            slope_abs = [entry[i] for entry in all_slope_abs]
            cutoff = -np.log10(p_score_cutoff)
            ln_pval = [-np.log10(p[0]) for p in pvals]
            for s,p in zip(slope_abs,ln_pval):
                if p > cutoff and s > 0: 
                    axs.scatter(s,p,color = 'b', label = f'{prot} Int Standardized',s = 2, alpha = 0.5)
                else:
                    axs.scatter(s,p,color = 'r', label = f'{prot} Int Standardized',s = 2, alpha = 0.5)
            axs.plot([min(slope_abs),max(slope_abs)], [-np.log10(p_score_cutoff),-np.log10(p_score_cutoff)],
                     color = 'grey', linestyle = '--', label = 'P score cutoff')
            figname = os.path.join(results_dir,f'P-value vs Concentration-Dependent Enrichment {prot}.png')    
            axs.set_title(f'CDE vs P-Value of {prot}', weight = 'bold', fontsize = 16)
            axs.set_xlabel('Ln(CDE)', weight = 'bold', fontsize = 14)
            axs.set_ylabel('-Log10(P-value)', weight = 'bold', fontsize = 14)
            axs.grid(True, which='major', color='dimgray', linestyle='-')
            axs.grid(True, which='minor', color='lightgray', linestyle='--')
            plt.tight_layout()     
            plt.savefig(figname, bbox_inches="tight")
            plt.show()
            
    def pval_vs_cde_vs_selscore(self,prots,pvals,all_slope_abs,enrichment,results_dir,intstd = False):
        '''
        Plots a 3D graph of the p-value, concentration-dependent enrichment, and selectivity score
        '''
        for i,prot in enumerate(prots):
            fig = plt.figure()
            axs = fig.add_subplot(projection='3d')
            print(prot)
            selscores = enrichment[i]
            ln_pvals = [-np.log10(p[0]) for p in pvals]
            slope_abs = [entry[i] for entry in all_slope_abs]
            if intstd:
                figname = os.path.join(results_dir,f'P-value vs Concentration-Dependent Enrichment vs Selectivity Intstd {prot}.png')  
                axs.set_title(f'P-value vs CDE vs Intstd Selectivity Score of {prot}', weight = 'bold', fontsize = 14)
            else:
                figname = os.path.join(results_dir,f'P-value vs Concentration-Dependent Enrichment vs Selectivity {prot}.png')  
                axs.set_title(f'P-value vs CDE vs Selectivity Score of {prot}', weight = 'bold', fontsize = 14)
            axs.scatter(selscores,slope_abs,ln_pvals, c = 'k', marker = 'o', s = 2, alpha = 0.5)                        
            axs.set_zlabel('-Log10(P-value)', weight = 'bold', fontsize = 14,labelpad=-1.5)
            axs.set_ylabel('Ln(CDE)', weight = 'bold', fontsize = 14)
            axs.set_xlabel('Selectivity score', weight = 'bold', fontsize = 14)
            plt.tight_layout(pad=5)     
            plt.savefig(figname, bbox_inches="tight")
            plt.show()
            
class Generate_FASTA:
    def __init__(self):
        pass
    
    def export_fasta(self, sequences, filtered_permutations, start_index, end_index, save_dir):
        '''
        Writes filtered permutations of potential peptide sequences to a FASTA file
        '''
        with open(os.path.join(save_dir,f'filtered_permutations_{start_index}_to_{end_index}.fasta'), 'w') as file:
            counter = 0
            for i, seq in enumerate(sequences[start_index:end_index]):
                for j,perm in enumerate(list(set(filtered_permutations[i]))):
                    file.write(f'>Sequence_{counter + j}\n')
                    file.write(f'{perm}\n')
                counter += len(filtered_permutations[i])

    def find_combinations(self,AA_set, input_mass, length, ppm_tol = 5, amidated = True, 
                          fixed_mass = [128.094963050]):
        '''
        Finds the possible sequences (within a given ppm tolerance) for a given input mass using a supplied monomer set.
        Length refers to the length of the variable region, any fixed residues can be lumped into "fixed_resis"
        and will be subtracted from the target mass
        Output strings do not have the fixed residues in, they will be added later to prevent permutations 
        of sequences exported to FASTA from getting too large
        '''

        # Get all keys and values from the dictionary
        keys = list(AA_set.keys())
        values = list(AA_set.values())

        # Generate combinations of length + 1 values
        combinations = combinations_with_replacement(values, length)
        target_sum = input_mass - sum(fixed_mass)
        target_tol = target_sum*(ppm_tol/10**6)
        lower = target_sum - target_tol
        upper = target_sum + target_tol
        seqs = []
        if amidated:
            fixed_mass_sum = sum(fixed_mass) + 17.0265
        else:
            fixed_mass_sum = sum(fixed_mass) + 18.0106
        # Check combinations for the desired sum
        for combo in combinations:
            mass = sum(combo) + fixed_mass_sum      
            if lower < mass < upper:
                # Retrieve corresponding keys for the values in the combination
                keys_for_combo = [f'{keys[values.index(val)]}{i+1}' for i,val in enumerate(combo)]
                seqs.append(dict(zip(keys_for_combo, combo)))
        return seqs if seqs else None

    def batch_fasta_export(self, sequences,save_dir, batch_size = 5000, fixed_resis = ['K'],locs = [-1]):
        '''
        Takes a list of sequences, finds relevant permutations, adds a set of fixed residues at their given
        locations, and exports the list of sequences in batches. Batch size of <10000 recommended due to file 
        size constraints.
        '''
        processed = 0

        # Process sequences and save progress
        while processed < len(sequences):
            print(f'Processed {processed} sequences')
            batch_sequences = sequences[processed:processed + batch_size]
            filtered_permutations = []

            for seq in batch_sequences:
                perms = [''.join(p) for p in permutations(seq)]
                perms_withconst = []
                for original_str in perms:
                    modified_str = original_str
                    for letter, position in zip(fixed_resis, locs):
                        modified_str = modified_str[:position] + letter + modified_str[position:]
                    perms_withconst.append(modified_str)
                filtered_permutations.append(perms_withconst)

            self.export_fasta(sequences, filtered_permutations, processed, processed + len(filtered_permutations),save_dir)
            processed += batch_size

    def make_fastas(self,obs_masses,AA_set, save_dir, name, 
                    fixed_resis = ['K'], locs = [-1],
                    batch_size = 5000, num_seqs = 500, length = 9, ppm_tol = 1,fixed_mass = [128.094963050]):
        '''
        Generates potential peptide sequences for a list of masses with a given monomer set, then generates 
        all possible permutations for the sequences and exports them to FASTA files for each protein
        '''
        if not os.path.isdir(os.path.join(save_dir,f'FASTAs of {name}')):         
            os.mkdir(os.path.join(save_dir,f'FASTAs of {name}'))  
        export_dir = os.path.join(save_dir,f'FASTAs of {name}')
        all_possible = [self.find_combinations(AA_set,mass,length, ppm_tol = ppm_tol,fixed_mass = fixed_mass) for mass in obs_masses[0:num_seqs]]
        none_indices = [index for index, value in enumerate(all_possible) if value is None]
        print(f'Number of masses that resulted in no sequences: {len(none_indices)}')
        try:
            combined_list = [item for sublist in all_possible if sublist is not None for item in sublist if item is not None]
            print(f'{len(combined_list)} total sequences')
            combined_strings_list = []
            # Loop through each dictionary in the list
            for dictionary in combined_list:
                # Extract first letters from keys and combine into a string for each dictionary
                first_letters = {key: key[0] for key in dictionary}
                combined_string = ''.join(value for value in first_letters.values())

                # Append the combined string to the list
                combined_strings_list.append(combined_string)

            self.batch_fasta_export(combined_strings_list,
                               export_dir,
                              batch_size = batch_size, fixed_resis = fixed_resis, locs = locs)
        except TypeError as e:
            print(f'Warning - no sequences found, check library length for mass size requirements. {e}')            
