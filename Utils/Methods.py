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
from itertools import combinations
import warnings
import peakutils

def sys_checks(parent_dir,prots,centroided,enr_score_cutoff,p_score_cutoff,peak_RT,peps_sec,reps,nr_EICs,LOD,mzML):
    openms_path = os.path.join(os.getcwd(), 'venv', 'Lib', 'site-packages', 'pyopenms', 'share', 'OpenMS')
    if platform.system() == 'Windows':     # Windows (either 32-bit or 64-bit)
        os.system(f'set OPENMS_DATA_PATH={openms_path};%OPENMS_DATA_PATH%')
        os.environ["OPENMS_DATA_PATH"] = openms_path # try 2 ways so maybe it doesn't get mad
    elif platform.system() == "Linux" or "Darwin":      # linux or Mac OS (Darwin)
        # os.system(f'export OPENMS_DATA_PATH = {openms_path}')
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
        if enr_score_cutoff > 1 or enr_score_cutoff < 0:
            mess = 'Enrichment score cutoff out of range (0 to 1)'
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

def directory_setup(parent_dir,folder,save_dir,prots,centroided):
    if not os.path.isdir(save_dir):         # Checks if data is already centroided, and if so it will create save directory
        os.mkdir(save_dir)                  # (which is same as data directory if centroided = False)

    eic_dirs = []
    eic_dirs_spec = []
    eic_dirs_nonspec = []
    eic_dirs_other = []
    for name in prots:
        eic_dir = os.path.join(parent_dir,folder,'EICs of ' + name)
        eic_dirs.append(eic_dir)
        eic_dir_spec = os.path.join(eic_dir,'EICs of ' + name + ' Specific')
        eic_dirs_spec.append(eic_dir_spec)
        eic_dir_nonspec = os.path.join(eic_dir,'EICs of ' + name + ' Nonspecific')
        eic_dirs_nonspec.append(eic_dir_nonspec)
        eic_dir_other = os.path.join(eic_dir,'EICs of ' + name + ' Unclear Specificity')
        eic_dirs_other.append(eic_dir_other)
        if not os.path.isdir(eic_dir):
            os.mkdir(eic_dir)
        if not os.path.isdir(eic_dir_spec):
            os.mkdir(eic_dir_spec)
        if not os.path.isdir(eic_dir_nonspec):
            os.mkdir(eic_dir_nonspec)
        if not os.path.isdir(eic_dir_other):
            os.mkdir(eic_dir_other)

    png_dir = os.path.join(parent_dir,folder,'PNGs')
    if not os.path.isdir(png_dir):
        os.mkdir(png_dir)
        
    if not centroided:
        PRTC_dir = os.path.join(parent_dir,folder,'PRTC Data')
        if not os.path.isdir(PRTC_dir):
            os.mkdir(PRTC_dir)

    return eic_dirs,eic_dirs_spec,eic_dirs_nonspec,eic_dirs_other,png_dir,PRTC_dir

def data_centroiding(centroided,parent_dir,folder,PRTC_dir,data_dir):
    if centroided:
        data_dir = os.path.join(parent_dir,'PRTC Time Points 2022_07_14-12_36_18_PM','centroided') # would need to be changed if true
        PRTC_dir = os.path.join(parent_dir)
    else:
        if not os.path.isdir(data_dir):         # Checks that the data directory exists
            os.mkdir(data_dir)       

        if os.path.exists(os.path.join(parent_dir,folder,'centroided')):   # If the same data directory exists, delete and remake it 
            shutil.rmtree(os.path.join(parent_dir,folder,'centroided'))    #   (though it shouldn't thanks to the time stamp)
        os.mkdir(os.path.join(parent_dir,folder, 'centroided'))

        count = 1
        for file in sorted(os.listdir(parent_dir)):               # Iterate through each folder and centroid the data
            if file.endswith('.mzML'):
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
                pickme.pickExperiment(exp_raw, exp_centroid)
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

def data_visualization(check_data,data_dir,png_dir,parent_dir):
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
            plt.savefig(os.path.join(png_dir,names[i] + ' Residuals.png'))

    '''centroided'''
    if check_data:
        for file in sorted(os.listdir(data_dir)):
            if file.endswith('mzML'):
                exp_2D = MSExperiment()
                MzMLFile().load(os.path.join(data_dir,file),exp_2D)
                plot_spectra_2D_overview(exp_2D,file[:-5],'centroided')
                
    '''profile'''
    if check_data:
        for file in sorted(os.listdir(parent_dir)):
            if file.endswith('mzML'):
                exp_2D = MSExperiment()
                MzMLFile().load(os.path.join(data_dir,file),exp_2D)
                plot_spectra_2D_overview(exp_2D,file[:-5],'profile')


def plot_spectra_2D_overview(exp,title,other,png_dir):
    '''Utilizes bilinear interpolation to graph the 2D RT vs m/z plot more quickly.'''
    '''
    Example:
    exp_2D = MSExperiment() # opens a new empty class to accept data
    MzMLFile().load(os.path.join(data_dir,file),exp_2D) # loads mzML file into class
    plot_spectra_2D_overview(exp_2D,file[:-5])
    exp is the MS experiment'''                                          
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
    plt.savefig(os.path.join(png_dir,filename))
    plt.show()
    print('showing plot...')

def alignment_check(consensus_map,feature_maps,png_dir,new_map_RTs,files,original_RTs_dict):
    ### OLD VERSION, CHECK THAT IT EVEN WORKS
    feature_map_RTs_aligned = []
    consensus_RTs_aligned = []
    consensus_mzs_aligned = [cf.getMZ() for cf in consensus_map]
    for i,fm in enumerate(feature_maps):
        consensus_RTs_aligned.append([feat[i] for feat in new_map_RTs])
        feature_map_RTs_aligned.append([f.getRT() for f in fm])
        plt.plot([f.getRT() for f in fm],[f.getMZ() for f in fm],'k.',markersize=10,label='feature map RTs')
        plt.plot([feat[i] for feat in new_map_RTs],consensus_mzs_aligned,'y.',label='aligned feature RTs')
        plt.savefig(os.path.join(png_dir,'Alignment Check Should Overlap ' + files[i] + '.png'))
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
        plt.savefig(os.path.join(png_dir,'Alignment Check Should Not Overlap ' + files[i] + '.png'))
        plt.show()
        
    fixed_RT_maps_formatted = []
    for i in range(len(new_map_RTs)):
        feat_RTs = [RT[i] for RT in fixed_RT_maps]
        fixed_RT_maps_formatted.append(feat_RTs)
        
    return fixed_RT_maps_formatted

def PRTC_stats(areas_ref_full,files,ref_feature_mass,png_dir):
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
        plt.show()
        plt.savefig(os.path.join(png_dir,f'PRTC Normalization {np.round(ref_feature_mass[i],3)}.png'),bbox_inches='tight')

def games_howell_scipy(names,groups,data,reps=3,alpha=0.05):  
    #'''stealing this from https://aaronschlegel.me/games-howell-post-hoc-multiple-comparisons-test-python.html
    #   since pingouin keeps failing'''
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

def checkList(lst):
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

def find_nearest_tol(array, value, tol):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    if (np.abs(array[idx] - value)).min() > tol:
        return -1,-1
    else:
        return array[idx],idx
    
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx],idx

def smoother(ints,kernal_size):  # try this
    kernal = np.ones(kernal_size)/kernal_size
    return np.convolve(ints,kernal,mode='same')

def signaltonoise(a, axis=0, ddof=0):
    a = np.asanyarray(a)
    m = a.mean(axis)
    sd = a.std(axis=axis, ddof=ddof)
    return np.where(sd == 0, 0, m/sd)

def get_data(directory):
    RTs = []
    ints = []
    mzs = []
    for file in sorted(os.listdir(directory)):
        if file.endswith('.mzML'):
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

def feature_int_extractor(m_z_feature_list,RT_feature_list,RTs_orig_list,mzs,RTs,ints,LOD=1.5E4,
                          noise_level=150,peak_range=120,decimal=2,percent=1,subtract=False):
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
        for k in range(len(RTs)):             # loop through each replicate
            RT_rep_feature = RTs[k]        # get list of retention times for given feature
            Int_rep_feature = ints[k]      # get list of ints arrays for whole replicate
            Mz_rep_feature = mzs[k]
            if RTs_orig[k] != 0:
                RT = RTs_orig[k]
            RT_idx_low = find_nearest(RT_rep_feature,RT-peak_range/2)[1]
            RT_idx_high = find_nearest(RT_rep_feature,RT+peak_range/2)[1]  
            RT_window = RT_rep_feature[RT_idx_low:RT_idx_high+1]
            ints_slice = Int_rep_feature[RT_idx_low:RT_idx_high+1]    # should take a slice of arrays for both of these
            mzs_slice = Mz_rep_feature[RT_idx_low:RT_idx_high+1]
#             noise = np.random.randint(noise_level-150,noise_level+150,len(RTs))   noise[np.random.randint(1,len(RTs))-1]
            mz_idx = [find_nearest_tol(entry,mz,10**(-decimal)/2)[1] for entry in mzs_slice] # find idx where the feature mz is in each scan
            max_int_candidates = [array[mz_idx[o]] if mz_idx[o] > 0 else 1 for o,array in enumerate(ints_slice)]
            # use 15000 for baseline value
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

def feature_int_extractor_start(m_z_feature_list,RT_feature_list,mzs,RTs,ints,LOD=1.5E4,
                                noise_level=150,peak_range=120,decimal=2,percent=1,subtract=False):
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
        for k in range(len(RTs)):             # loop through each replicate
            RT_rep_feature = RTs[k]        # get list of retention times for given feature
            Int_rep_feature = ints[k]      # get list of ints arrays for whole replicate
            Mz_rep_feature = mzs[k]
            RT_idx_low = find_nearest(RT_rep_feature,RT-peak_range/2)[1]
            RT_idx_high = find_nearest(RT_rep_feature,RT+peak_range/2)[1]  
            RT_window = RT_rep_feature[RT_idx_low:RT_idx_high+1]
            ints_slice = Int_rep_feature[RT_idx_low:RT_idx_high+1]    # should take a slice of arrays for both of these
            mzs_slice = Mz_rep_feature[RT_idx_low:RT_idx_high+1]
#             noise = np.random.randint(noise_level-150,noise_level+150,len(RTs))   noise[np.random.randint(1,len(RTs))-1]
            mz_idx = [find_nearest_tol(entry,mz,10**(-decimal)/2)[1] for entry in mzs_slice] # find idx where the feature mz is in each scan
            max_int_candidates = [array[mz_idx[o]] if mz_idx[o] > 0 else 1 for o,array in enumerate(ints_slice)]
            # use 15000 for baseline value
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

def peak_finder_savgol(rt_windows,ints_windows,names,plot=False,reps=3,default_thresh=1.5E5,
                       kernal_size=10,width=5,prominence=2,threshold=2,rel_height=0.5):
    all_peaks = []
    for i,feature in enumerate(ints_windows):  # each entry in ints_windows is a feature
        feature_peaks = []
        rt_feature = rt_windows[i]
        is_peaks = False
        if plot:
            fig,axs = plt.subplots(1,len(names)*reps,figsize=(16,12),sharey=True)
        for j,rep in enumerate(feature):       # each feature should have 6 reps (or however many for your selection)
            if len(rep) <= 10:
                continue
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

def feature_area_extractor_savgol(rt_windows_filt,int_windows_filt,check=False,
                                  width_start=5,prominence=2,threshold=2,rel_height=0.5):
    all_areas = []
    for i,feature in enumerate(int_windows_filt):  # each entry in ints_windows is a feature
        feature_areas = []
        rt_feature = rt_windows_filt[i]
        if check:
            fig, axs = plt.subplots(1,6,figsize=(16,10),sharey=True)
        for j,rep in enumerate(feature):       # each feature should have 6 reps (or however many for your selection)
            is_peaks = False
            x_window = 19
            polyorder = 9
            if len(rep) <= 9:
                feature_areas.append(np.random.randint(9.9E4,1.1E5))  # baseline area will be around 1E5 
                continue
            if len(rep) < x_window:
                x_window = 9
                polyorder = 5
            rep_smoothed = savgol_filter(rep,x_window,polyorder,mode='interp')
            height = np.percentile(rep_smoothed,99)
            width = width_start
            while not is_peaks:
                if width < 0 or height < 0:
                    feature_areas.append(1)  # placeholder non-zero value for scoring
                    break
                peaks,props = find_peaks(rep_smoothed,width=[width,50],prominence=prominence,threshold=threshold,height=height,rel_height=rel_height)                
                width = width - 1
                height = height*0.95
                if np.any(peaks):
                    peak_ints = rep_smoothed[peaks]
                    peak_max_idx = np.where(peak_ints == max(peak_ints))[0][0]
                    peak_max = np.where(rep_smoothed == max(peak_ints))[0]
                    peak_width = math.ceil(props['widths'][peak_max_idx])
                    if peak_width < 10:
                        peak_width = 10
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
                    area_trap = cumulative_trapezoid(ints_subwindow,x=peak_subwindow)[-1]
                    feature_areas.append(area_trap)
                    is_peaks = True
            if check:
                axs[j].plot(rt_feature[j],rep,'k-',linewidth=2)
                axs[j].plot(rt_feature[j],rep_smoothed,'c--')
                axs[j].plot(peak_subwindow,ints_subwindow,'m--',label=f'A={area_trap:.2E}',linewidth=4)
                axs[j].fill_between(peak_subwindow,ints_subwindow)
                axs[j].legend(loc='best')
        if check:
            plt.show()
            plt.close(fig)
        all_areas.append(feature_areas)
        if i % 1000 == 0:
            print(i)
    return all_areas

def enr_scoring(max_int_features,ref_val,prot_names,prots,reps=3,p_score_cutoff=0.05):
    # redundancy in prot_name and prots?
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
            print(feat,group)
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
                    if p <= 0.05:
                        diff = stats_out_GH['diff'].values[i]
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
        else:
            print('UNEQUAL VARIANCES, USING GAMES-HOWELL')
            stats_out = pg.welch_anova(dv='score', between='group', data=df)
            pval = stats_out['p-unc'].values
            if pval < p_score_cutoff:
                stats_out_GH = games_howell_scipy(prots,group,feat,reps=reps)
                GH_ps = stats_out_GH['p_value'].values
                id_list = []
                for i,p in enumerate(GH_ps):
                    if p <= 0.05:
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
    enr_scores = []
    enr_scores_normalized = []
    for m in range(len(prot_names)):                          # now we divide and get enrichment scores
        e_prot = []
        e_prot_normalized = []
        for feature in int_sums:
            int_sum_prot = feature[m]
            e_score = int_sum_prot/np.sum(feature)
            e_prot.append(e_score)
            e_score_normalized = int_sum_prot/ref_val
            e_prot_normalized.append(e_score_normalized)
        enr_scores.append(e_prot)
        enr_scores_normalized.append(e_prot_normalized)
    return np.asarray(enr_scores),np.asarray(enr_scores_normalized),p_vals,specificity,spec_label

def setup_output(prots,spec_label,enrichmentscores,enrichment_normalized,pvals,
                 feat_RT_filtered,feat_RT_orig_filtered,feat_mz_filtered,feat_z_filtered,areas_savgol):
    es_graphing = [ [] for _ in range(len(prots))]
    es_normalized_graphing = [ [] for _ in range(len(prots))]
    ps_graphing = [ [] for _ in range(len(prots))]
    RTs_graphing = [ [] for _ in range(len(prots))]
    RTs_graphing_orig = [ [] for _ in range(len(prots))]
    mzs_graphing = [ [] for _ in range(len(prots))]
    z_graphing = [ [] for _ in range(len(prots))]
    areas_graphing = [ [] for _ in range(len(prots))]

    es_nonspecific = []
    es_normalized_nonspecific = []
    ps_nonspecific = []
    RTs_nonspecific = []
    RTs_nonspecific_orig = []
    mzs_nonspecific = []
    z_nonspecific = []
    areas_nonspecific = []

    es_unclear = []
    es_normalized_unclear = []
    ps_unclear = []
    RTs_unclear = []
    RTs_unclear_orig = []
    mzs_unclear = []
    z_unclear = []
    areas_unclear = []


    for i,label in enumerate(spec_label):
        if label in prots:
            prot_name = np.where(np.array(prots) == label)[0][0]  # actually a number lol
            es_graphing[prot_name].append(enrichmentscores[prot_name][i])
            es_normalized_graphing[prot_name].append(enrichment_normalized[prot_name][i])
            ps_graphing[prot_name].append(pvals[i])
            RTs_graphing[prot_name].append(feat_RT_filtered[i])
            RTs_graphing_orig[prot_name].append(feat_RT_orig_filtered[i])
            mzs_graphing[prot_name].append(feat_mz_filtered[i])
            z_graphing[prot_name].append(feat_z_filtered[i])
            areas_graphing[prot_name].append(areas_savgol[i])
        elif label == 'Nonspecific':
            es_nonspecific.append(enrichmentscores[:,i])
            es_normalized_nonspecific.append(enrichment_normalized[:,i])
            ps_nonspecific.append(pvals[i])
            RTs_nonspecific.append(feat_RT_filtered[i])
            RTs_nonspecific_orig.append(feat_RT_orig_filtered[i])
            mzs_nonspecific.append(feat_mz_filtered[i])
            z_nonspecific.append(feat_z_filtered[i])
            areas_nonspecific.append(areas_savgol[i])
        else:
            es_unclear.append(enrichmentscores[:,i])
            es_normalized_unclear.append(enrichment_normalized[:,i])
            ps_unclear.append(pvals[i])
            RTs_unclear.append(feat_RT_filtered[i])
            RTs_unclear_orig.append(feat_RT_orig_filtered[i])
            mzs_unclear.append(feat_mz_filtered[i])
            z_unclear.append(feat_z_filtered[i])
            areas_unclear.append(areas_savgol[i])
    return es_graphing,es_normalized_graphing,ps_graphing,RTs_graphing,RTs_graphing_orig,mzs_graphing,z_graphing,areas_graphing,es_nonspecific,es_normalized_nonspecific,ps_nonspecific,RTs_nonspecific,RTs_nonspecific_orig, \
    mzs_nonspecific,z_nonspecific,areas_nonspecific, es_unclear,es_normalized_unclear,ps_unclear,RTs_unclear,RTs_unclear_orig,mzs_unclear,z_unclear,areas_unclear

def enrichment_rankings(prots,es_graphing,es_normalized_graphing,png_dir):
    if len(prots) > 2:
        fig,axs = plt.subplots(math.ceil(len(prots)/2),2,figsize=(14,14),sharey=True)
        for i in range(math.ceil(len(prots)/2)):
            for j in range(2):
                try:
                    scores = es_graphing[2*i+j]
                    scores = np.sort(scores)[::-1]
                    scores_normalized = es_normalized_graphing[2*i+j]
                    scores_normalized = np.sort(scores_normalized)[::-1]
                    axs[i,j].plot(scores,'k-',label=f'{prots[2*i+j]} weighted')
                    axs[i,j].set_yscale('linear')
                    axs[i,j].set_title(f'Enrichment of {prots[2*i+j]}',fontsize=20)
                    axs[i,j].set_xlabel('Enrichment rank',fontsize=16)
                    axs[i,j].set_ylabel('Enrichment score',fontsize=16)
                    axs[i,j].grid(b=True, which='major', color='dimgray', linestyle='-')
                    axs[i,j].grid(b=True, which='minor', color='lightgray', linestyle='--')
                    axs[i,j].tick_params(axis='both', which='major', labelsize=14)
                    ax2 = axs[i,j].twinx()
                    ax2.plot(scores_normalized,'r-',label=f'{prots[2*i+j]} normalized')
                    ax2.set_yscale('log')
                    ax2.legend(loc='best')
                except IndexError:
                    continue
        plt.tight_layout()
        figname = os.path.join(png_dir,'Final Enrichments.png')
        plt.savefig(figname)
        plt.show()

        fig,axs = plt.subplots(math.ceil(len(prots)/2),2,figsize=(14,14),sharey=True)
        for i in range(math.ceil(len(prots)/2)):
            for j in range(2):
                try:
                    scores = es_graphing[2*i+j]
                    scores = np.sort(scores)[::-1]
                    scores_normalized = es_normalized_graphing[2*i+j]
                    scores_normalized = np.sort(scores_normalized)[::-1]
                    axs[i,j].plot(scores[0:1000],'k-',label=f'{prots[2*i+j]} weighted')
                    axs[i,j].set_yscale('linear')
                    axs[i,j].set_title(f'Enrichment of {prots[2*i+j]}',fontsize=20)
                    axs[i,j].set_xlabel('Enrichment rank',fontsize=16)
                    axs[i,j].set_ylabel('Enrichment score',fontsize=16)
                    axs[i,j].grid(b=True, which='major', color='dimgray', linestyle='-')
                    axs[i,j].grid(b=True, which='minor', color='lightgray', linestyle='--')
                    axs[i,j].tick_params(axis='both', which='major', labelsize=14)
                    ax2 = axs[i,j].twinx()
                    ax2.plot(scores_normalized[0:1000],'r-',label=f'{prots[2*i+j]} normalized')
                    ax2.set_yscale('log')
                    ax2.legend(loc='best')
                except IndexError:
                    continue
        plt.tight_layout()
        figname = os.path.join(png_dir,'Final Enrichments Zoomed.png')
        plt.savefig(figname)
        plt.show()

    else:
        fig,axs = plt.subplots(1,2,figsize=(14,10),sharey=True)
        for j in range(2):
            scores = es_graphing[j]
            scores = np.sort(scores)[::-1]
            scores_normalized = es_normalized_graphing[j]
            scores_normalized = np.sort(scores_normalized)[::-1]
            axs[j].plot(scores,'k-',label=f'{prots[j]} weighted')
            axs[j].set_yscale('linear')
            axs[j].set_title(f'Enrichment of {prots[j]}')
            axs[j].set_xlabel('Enrichment rank')
            axs[j].set_ylabel('Enrichment score')
            axs[j].grid(b=True, which='major', color='dimgray', linestyle='-')
            axs[j].grid(b=True, which='minor', color='lightgray', linestyle='--')
            ax2 = axs[j].twinx()
            ax2.plot(scores_normalized,'r-',label=f'{prots[j]} normalized')
            ax2.set_yscale('log')
            ax2.legend(loc='best')
        plt.tight_layout()
        figname = os.path.join(png_dir,'Final Enrichments.png')
        plt.savefig(figname)
        plt.show()

        fig,axs = plt.subplots(1,2,figsize=(14,10),sharey=True)
        for j in range(2):
            scores = es_graphing[j]
            scores = np.sort(scores)[::-1]
            scores_normalized = es_normalized_graphing[j]
            scores_normalized = np.sort(scores_normalized)[::-1]
            axs[j].plot(scores[0:1000],'k-',label=f'{prots[j]} weighted')
            axs[j].set_yscale('linear')
            axs[j].set_title(f'Enrichment of {prots[j]}')
            axs[j].set_xlabel('Enrichment rank')
            axs[j].set_ylabel('Enrichment score')
            axs[j].grid(b=True, which='major', color='dimgray', linestyle='-')
            axs[j].grid(b=True, which='minor', color='lightgray', linestyle='--')
            ax2 = axs[j].twinx()
            ax2.plot(scores_normalized[0:1000],'r-',label=f'{prots[j]} normalized')
            ax2.set_yscale('log')
            ax2.legend(loc='best')
        plt.tight_layout()
        figname = os.path.join(png_dir,'Final Enrichments Zoomed.png')
        plt.savefig(figname)
        plt.show()

def volcano_plotting(prots,es_graphing,es_nonspecific,es_unclear,ps_graphing,ps_nonspecific,ps_unclear,
                     feat_mz_combined,enrichmentscores,p_score_cutoff,enr_score_cutoff,save_dir):
    if len(prots) > 2:
        fig, axs = plt.subplots(math.ceil(len(prots)/2),2,figsize=(18,18))
        for i in range(math.ceil(len(prots)/2)):
            for j in range(2):
                try:
                    e_spec = es_graphing[2*i+j]
                    e_nonspec = [e_val[2*i+j] for e_val in es_nonspecific]
                    e_unclear = [e_val[2*i+j] for e_val in es_unclear]

                    p_spec = [-np.log10(p) for p in ps_graphing[2*i+j]]
                    p_nonspec = [-np.log10(p) for p in ps_nonspecific]
                    p_unclear = [-np.log10(p) for p in ps_unclear]
                    axs[i,j].scatter(e_unclear,p_unclear,c='r',label='Unclear')
                    axs[i,j].scatter(e_nonspec,p_nonspec,c='gray',label='Nonspecific')
                    axs[i,j].scatter(e_spec,p_spec,c='b',label='Specific')   
                    axs[i,j].plot([min(e_spec + e_nonspec + e_unclear),max(e_spec + e_nonspec + e_unclear)],[-np.log10(p_score_cutoff),-np.log10(p_score_cutoff)],'k--')
                    axs[i,j].plot([enr_score_cutoff,enr_score_cutoff],[min(p_spec+p_nonspec+p_unclear),max(p_spec+p_nonspec+p_unclear)],'k--')
                    axs[i,j].set_xscale('linear')
                    axs[i,j].set_yscale('linear')
                    axs[i,j].set_title(f'Volcano plot for {prots[2*i+j]}',fontsize=24)
                    axs[i,j].set_xlabel('Enrichment Score',fontsize=20)
                    axs[i,j].set_ylabel('-Log10(P value)',fontsize=20)
                    axs[i,j].tick_params(axis='both', which='major', labelsize=16)
                except IndexError:
                    continue
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        figname = 'Volcano Plots Specific.png'
        plt.savefig(os.path.join(save_dir,figname))
        plt.show()          
    else:
        fig, axs = plt.subplots(1,2,figsize=(14,10))
        for i in range(len(prots)):
            e_spec = es_graphing[i]
            e_nonspec = [e_val[i] for e_val in es_nonspecific]
            e_unclear = [e_val[i] for e_val in es_unclear]
            
            p_spec = [-np.log10(p) for p in ps_graphing[i]]
            p_nonspec = [-np.log10(p) for p in ps_nonspecific]
            p_unclear = [-np.log10(p) for p in ps_unclear]
            axs[i].scatter(e_unclear,p_unclear,c='r',label='Unclear')
            axs[i].scatter(e_nonspec,p_nonspec,c='gray',label='Nonspecific')
            axs[i].scatter(e_spec,p_spec,c='b',label='Specific')   
            axs[i].plot([min(e_spec + e_nonspec + e_unclear),max(e_spec + e_nonspec + e_unclear)],[-np.log10(p_score_cutoff),-np.log10(p_score_cutoff)],'k--')
            axs[i].plot([1/len(prots),1/len(prots)],[min(p_spec+p_nonspec+p_unclear),max(p_spec+p_nonspec+p_unclear)],'k--')
            axs[i].set_xscale('linear')
            axs[i].set_xlim([0,1])
            axs[i].set_yscale('linear')
            axs[i].set_title(f'Volcano plot for {prots[i]}',fontsize=24)
            axs[i].set_xlabel('Enrichment Score',fontsize=20)
            axs[i].set_ylabel('-Log10(P value)',fontsize=20)
            axs[i].legend(loc='best',fontsize=16)
            axs[i].tick_params(axis='both',which='major', labelsize=16)
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        figname = 'Volcano Plots Specific.png'
        plt.savefig(os.path.join(save_dir,figname))
        plt.show()

def volcano_plots_normalized(prots,es_normalized_graphing,es_normalized_nonspecific,es_normalized_unclear,
                             ps_graphing,ps_nonspecific,ps_unclear,p_score_cutoff,enr_score_cutoff,save_dir):
    if len(prots) > 2:
        fig, axs = plt.subplots(math.ceil(len(prots)/2),2,figsize=(18,18))
        for i in range(math.ceil(len(prots)/2)):
            for j in range(2):
                try:
                    e_spec = es_normalized_graphing[2*i+j]
                    e_nonspec = [e_val[2*i+j] for e_val in es_normalized_nonspecific]
                    e_unclear = [e_val[2*i+j] for e_val in es_normalized_unclear]
                    p_spec = [-np.log10(p) for p in ps_graphing[2*i+j]]
                    p_nonspec = [-np.log10(p) for p in ps_nonspecific]
                    p_unclear = [-np.log10(p) for p in ps_unclear]
                    axs[i,j].scatter(e_unclear,p_unclear,c='r',label='Unclear')
                    axs[i,j].scatter(e_nonspec,p_nonspec,c='gray',label='Nonspecific')
                    axs[i,j].scatter(e_spec,p_spec,c='b',label='Specific')   
                    axs[i,j].plot([min(e_spec + e_nonspec + e_unclear),max(e_spec + e_nonspec + e_unclear)],[-np.log10(p_score_cutoff),-np.log10(p_score_cutoff)],'k--')
                    axs[i,j].plot([enr_score_cutoff,enr_score_cutoff],[min(p_spec+p_nonspec+p_unclear),max(p_spec+p_nonspec+p_unclear)],'k--')
                    axs[i,j].set_xscale('log')
                    axs[i,j].set_yscale('linear')
                    axs[i,j].set_xlim([min(e_spec),max(e_spec)])
                    axs[i,j].set_title(f'Volcano plot for {prots[2*i+j]}',fontsize=24)
                    axs[i,j].set_xlabel('Enrichment Score',fontsize=20)
                    axs[i,j].set_ylabel('-Log10(P value)',fontsize=20)
                    axs[i,j].tick_params(axis='both', which='major', labelsize=16)
                except IndexError:
                    continue
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        figname = 'Volcano Plots Specific normalized.png'
        plt.savefig(os.path.join(save_dir,figname))
        plt.show()          
    else:
        fig, axs = plt.subplots(1,2,figsize=(14,10))
        for i in range(len(prots)):
            e_spec = es_normalized_graphing[i]
            e_nonspec = [e_val[i] for e_val in es_normalized_nonspecific]
            e_unclear = [e_val[i] for e_val in es_normalized_unclear]
            
            p_spec = [-np.log10(p) for p in ps_graphing[i]]
            p_nonspec = [-np.log10(p) for p in ps_nonspecific]
            p_unclear = [-np.log10(p) for p in ps_unclear]
            axs[i].scatter(e_unclear,p_unclear,c='r',label='Unclear')
            axs[i].scatter(e_nonspec,p_nonspec,c='gray',label='Nonspecific')
            axs[i].scatter(e_spec,p_spec,c='b',label='Specific')   
            axs[i].plot([min(e_spec + e_nonspec + e_unclear),max(e_spec + e_nonspec + e_unclear)],[-np.log10(p_score_cutoff),-np.log10(p_score_cutoff)],'k--')
            axs[i].plot([1,1],[min(p_spec+p_nonspec+p_unclear),max(p_spec+p_nonspec+p_unclear)],'k--')
            axs[i].set_xscale('log')
            axs[i].set_yscale('linear')
            axs[i].set_xlim([1E-4,1E4])
            axs[i].set_title(f'Volcano plot for {prots[i]}',fontsize=24)
            axs[i].set_xlabel('Enrichment Score',fontsize=20)
            axs[i].set_ylabel('-Log10(P value)',fontsize=20)
            axs[i].legend(loc='best',fontsize=16)
            axs[i].tick_params(axis='both',which='major', labelsize=16)
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        figname = 'Volcano Plots Specific normalized.png'
        plt.savefig(os.path.join(save_dir,figname))
        plt.show()

def export_results(enrichmentscores,es_graphing,RTs_graphing,mzs_graphing,ps_graphing,z_graphing,enr_score_cutoff,p_score_cutoff,
                   prots,parent_dir,folder,feat_RT_cf_combined,feat_mz_filtered,feat_z_filtered,pvals,spec_label,areas_savgol,full_out):
    for i in range(len(enrichmentscores)):
        df = pd.DataFrame(columns=('Compound', 'm/z','z','p value','specificity'))
        score = es_graphing[i]
        RT = RTs_graphing[i]
        mz = mzs_graphing[i]
        ps = ps_graphing[i]
        zs = z_graphing[i]
        score_sorted = np.sort(score)[::-1]            # sort the scores by descending order
        indices_sorted = np.argsort(score)[::-1]        # get indices of scores to get retention time and m/z
        count = 0
        while score_sorted[count] >= enr_score_cutoff:# and len(df) < 10000:
            index = indices_sorted[count]
            if ps[index] < p_score_cutoff:
                A = 'EScore: ' + str(np.round(score_sorted[count],2)) + ' RTime: ' + str(np.round(RT[index],3))
                B = np.round(mz[index],4)
                D = zs[index]
                C = float(ps[index])
                E = prots[i]
                df.loc[count] = [A,B,D,C,E]
            count -= -1
            if count >= len(score):
                break
        df.to_csv(os.path.join(parent_dir,folder,prots[i] + '.csv'),index=False)
            
    if full_out:
        for i in range(len(enrichmentscores)):
            df = pd.DataFrame(columns=('Compound', 'm/z','z','p value','specificity','areas'))
            score = enrichmentscores[i]
            score_sorted = np.sort(score)[::-1]            # sort the scores by descending order
            indices_sorted = np.argsort(score)[::-1]        # get indices of scores to get retention time and m/z
            for j in range(len(score_sorted)):
                index = indices_sorted[j]
                A = 'EScore: ' + str(np.round(score_sorted[j],2)) + ' RTime: ' + str(np.round(feat_RT_cf_combined[index],3))
                B = np.round(feat_mz_filtered[index],4)
                D = feat_z_filtered[index]
                C = float(pvals[index])
                E = spec_label[index]
                F = [float(entry) for entry in areas_savgol[index]]
                df.loc[j] = [A,B,D,C,E,F]
            df.to_csv(os.path.join(parent_dir,folder,prots[i] + ' Full Output.csv'),index=False)

def inclusion_lists(data_dir,files_full,RTs_graphing,mzs_graphing,ps_graphing,z_graphing,
                    png_dir,prots,es_graphing,peps_sec,parent_dir,folder,feat_RT_cf_combined):
    get_RTs = MSExperiment()
    file_path = os.path.join(data_dir,files_full[0])
    MzMLFile().load(file_path,get_RTs)

    RT_start = get_RTs[0].getRT()
    RT_end = get_RTs[get_RTs.getNrSpectra()-1].getRT()
    n_bins = int(np.ceil(RT_end-RT_start))
    del get_RTs

    for i,RT_prot in enumerate(RTs_graphing):
        mzs = mzs_graphing[i]
        ps = ps_graphing[i]
        zs = z_graphing[i]
        plt.rcParams["figure.figsize"] = 10,5
        hist,bins = np.histogram(RT_prot,bins=np.linspace(RT_start,RT_end,n_bins),density=False)
        plt.hist(RT_prot,bins=np.linspace(RT_start,RT_end,n_bins),density=False)
        plt.xlim(np.linspace(RT_start,RT_end,n_bins)[0],np.linspace(RT_start,RT_end,n_bins)[-1])
        plt.savefig(os.path.join(png_dir,prots[i]+' specific features histogram.png'))   # bins too thin to see in output lmao
        plt.show()
        plt.rcParams["figure.figsize"] = plt.rcParamsDefault["figure.figsize"]
        
        fig,axs = plt.subplots(1,2,figsize=(14,14))
        score = es_graphing[i]
        bin_members = np.digitize(RT_prot,bins=np.linspace(RT_start,RT_end,n_bins))
        unique_members = np.unique(bin_members)
        new_score_array = []  
        for member in unique_members:
            repeats = np.where(bin_members == member)[0]
            scores_sub = [score[v] for v in repeats]
            if len(repeats) > peps_sec:
                scores_sub_sorted = np.sort(scores_sub)   # sort from lowest to highest
                while len(scores_sub) > peps_sec:
                    scores_sub = np.delete(scores_sub,np.where(scores_sub==scores_sub_sorted[0])[0][0])
                    scores_sub_sorted = np.delete(scores_sub_sorted,0)
            new_score_array = np.append(new_score_array,scores_sub)
        new_RT_array = [RT_prot[np.where(score==score_red)[0][0]] for score_red in new_score_array]
        new_mz_array = [mzs[np.where(score==score_red)[0][0]] for score_red in new_score_array]
        new_z_array = [zs[np.where(score==score_red)[0][0]] for score_red in new_score_array] 
        new_p_array = [ps[np.where(score==score_red)[0][0]] for score_red in new_score_array]
        axs[0].plot(score,RT_prot,'k.',markersize=10,label='Features removed')
        axs[0].plot(new_score_array,new_RT_array,'y.',label='Features remaining')
        axs[0].legend(loc='best')
        axs[1].scatter(score,[-np.log10(p) for p in ps],label='Original features',linewidths=6)
        axs[1].scatter(new_score_array,[-np.log10(p) for p in new_p_array],label='Features remaining')
        axs[1].set_xlabel('Enrichment Score',fontsize=20)
        axs[1].set_ylabel('-Log10(P value)',fontsize=20)
        plt.savefig(os.path.join(png_dir,prots[i] + ' specific features removed plot.png'))
        plt.show()
        
        df = pd.DataFrame(columns=('Compound', 'Formula','Adduct','m/z','z'))
        for j in range(len(new_score_array)):
            A = 'EScore: ' + str(np.round(new_score_array[j],2)) + ' p value: ' + "{:.3e}".format(new_p_array[j][0]) + ' RTime: ' + str(np.round(new_RT_array[j],3))
            B = None
            C = '(no adduct)'
            D = np.round(new_mz_array[j],4)
            E = int(new_z_array[j])
            df.loc[j] = [A,B,C,D,E]
        df.to_csv(os.path.join(parent_dir,folder,prots[i] + ' Inclusion List.csv'),index=False)
        
    plt.rcParams["figure.figsize"] = 10,5
    hist = np.histogram(feat_RT_cf_combined,bins=np.linspace(RT_start,RT_end,n_bins),density=False)
    plt.hist(feat_RT_cf_combined,bins=np.linspace(RT_start,RT_end,n_bins),density=False)
    plt.xlim(np.linspace(RT_start,RT_end,n_bins)[0],np.linspace(RT_start,RT_end,n_bins)[-1])
    plt.savefig(os.path.join(png_dir,'Consensus Map Histogram.png'))
    plt.rcParams["figure.figsize"] = plt.rcParamsDefault["figure.figsize"]
    plt.show()

def compare_known_binders(known_binders,prots,parent_dir,folder,feat_mz_combined_filtered,feat_mz_combined,feat_RT_cf_combined):
    if known_binders:
        known_mzs_found = {}
        known_rts_found = {}
        known_mzs_notfound = {}
        known_rts_notfound = {}
        for i in range(len(prots)):
            known_mzs_prot_found = []
            known_rts_prot_found = []
            known_mzs_prot_notfound = []
            known_rts_prot_notfound = []
            is_feature_in = 0
            df_compare1 = pd.DataFrame(columns = ('mz','rt range low','rt range high','rt ave','rt ave sec','z','Area Intensity'))
            df_compare2 = pd.DataFrame(columns = ('mz','rt range low','rt range high','rt ave','rt ave sec','z','Area Intensity'))
            with open(os.path.join(parent_dir,'FeatureList'+prots[i]+'.csv'),mode='r',encoding='utf-8-sig') as csv_file:
                csv_reader = csv.DictReader(csv_file,delimiter=',')
                j = 0
                k = 0
                for row in csv_reader:
                    if np.round(float(row['mz']),2) in np.round(feat_mz_combined_filtered,2):
                        known_mzs_prot_found.append(float(row['mz']))
                        known_rts_prot_found.append(float(row['rt ave sec']))
                        idx = np.where(np.round(feat_mz_combined,2) == np.round(float(row['mz']),2))[0]
                        rt_list = [feat_RT_cf_combined[i] for i in idx]
                        df_compare1.loc[j] = row
                        is_feature_in -= -1
                    else:
                        known_mzs_prot_notfound.append(float(row['mz']))
                        known_rts_prot_notfound.append(float(row['rt ave sec']))
                        df_compare2.loc[k] = row
                        k -= -1
                    j -= -1
                percent_feat_ID = is_feature_in/j*100
            known_mzs_found[i] = known_mzs_prot_found
            known_rts_found[i] = known_rts_prot_found
            known_mzs_notfound[i] = known_mzs_prot_notfound
            known_rts_notfound[i] = known_rts_prot_notfound
            df_compare1.to_csv(os.path.join(parent_dir,folder,prots[i] + ', no time req ' + str(np.round(percent_feat_ID,2)) + '% new list features found, cf total = ' + str(len(feat_mz_combined)) + '.csv'))
            df_compare2.to_csv(os.path.join(parent_dir,folder,prots[i] + ', no time req ' + str(np.round(100-percent_feat_ID,2)) + '% new list features not found.csv'))

def plot_PRTC(m_z_list,feature_RT,RT_orig_list,areas,directory,savedir,weight=True,
              show_plot=False,time=5,peak_range=30,decimal=2,reps=3,baseline=0):
    RTs_full = []
    ints_full = []
    mzs_full = []
    names = []
    for file in sorted(os.listdir(directory)):
        if file.endswith('.mzML'):
            names.append(file[:-5])
            exp = MSExperiment()
            MzMLFile().load(os.path.join(directory,file),exp)
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
    for j,mz in enumerate(m_z_list):       # loop through each feature
        RT_feature = feature_RT[j]
        RT_origs = RT_orig_list[j]
        area_feature = areas[j]
        xs = []
        xs_zoom = []
        ys = []
        max_array = []
        for k in range(len(RTs_full)):             # loop through each replicate
            RTs_rep = RTs_full[k]
            RT_orig = RT_origs[k]
            if RT_orig != 0:
                RT_use = RT_orig
            else:
                RT_use = RT_feature
            ints_rep = ints_full[k]
            mzs_rep = mzs_full[k]
            
            RT_idx_low = find_nearest(RTs_rep,RT_use-time/2*60)[1]
            RT_idx_high = find_nearest(RTs_rep,RT_use+time/2*60)[1] 

            RTs_slice = RTs_rep[RT_idx_low:RT_idx_high+1]
            
            max_idx_low = find_nearest(RTs_slice,RT_use-peak_range/2)[1]
            max_idx_high = find_nearest(RTs_slice,RT_use+peak_range/2)[1]

            ints_slice = ints_rep[RT_idx_low:RT_idx_high+1]    # should take a slice of arrays for both of these
            mzs_slice = mzs_rep[RT_idx_low:RT_idx_high+1]
            mz_idx = [find_nearest_tol(entry,mz,10**(-decimal)/2)[1] for entry in mzs_slice] # find idx where the feature mz is in each scan

            ints_EIC = [array[mz_idx[o]] if mz_idx[o] > 0 else baseline for o,array in enumerate(ints_slice)]

            try:
                max_found = np.amax(ints_EIC[max_idx_low:max_idx_high])
            except ValueError:
                max_found = np.amax(ints_EIC)
                print(ints_EIC[max_idx_low:max_idx_high],ints_EIC)
            max_array.append(max_found)
            xs.append(RTs_slice)
            xs_zoom.append(RTs_slice[max_idx_low:max_idx_high])
            ys.append(ints_EIC)
        fig, axs = plt.subplots(int(len(names)/reps),reps,figsize=(14,16),sharey=True)
        for l in range(int(len(names)/reps)):
            for m in range(reps):
                x = [t/60 for t in xs[reps*l+m]]
                x_zoom = [t/60 for t in xs_zoom[reps*l+m]]
                y = ys[reps*l+m]

                axs[l,m].plot(x,y, label = f'Feature mean RT = {np.round(RT_feature/60,1)} min')
                axs[l,m].plot(x_zoom,np.ones(len(x_zoom))*max_array[reps*l+m],'k--',label=f'Area of {area_feature[reps*l+m]}')
                axs[l,m].set_title(f"EIC of m/z = {np.round(mz,decimal)} in {names[reps*l+m][:-5]}",fontsize=14)
                axs[l,m].set_xlabel('Retention time (min)',fontsize=12)
                axs[l,m].set_ylabel('Intensity',fontsize=12)
                yfmt = mticker.ScalarFormatter(useMathText=True)
                yfmt.set_powerlimits((3, 4))
                axs[l,m].yaxis.set_major_formatter(yfmt)
                axs[l,m].ticklabel_format(style='sci', axis='y', scilimits=(0,0))
                axs[l,m].get_yaxis().get_offset_text().set_visible(False)
                ax_max = max(axs[l,m].get_yticks())
                exponent_axis = np.floor(np.log10(ax_max)).astype(int)
                axs[l,m].annotate(r'$\times$10$^{%i}$'%(exponent_axis),
                             xy=(.01, .8), xycoords='axes fraction')
                axs[l,m].legend(loc='best',fontsize=10)
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        if weight:
            figname = f"mz of {np.round(mz,decimal)} at {np.round(RT_feature/60,1)}.png"
        else:
            figname = f"NORMALIZED mz of {np.round(mz,decimal)} at {np.round(RT_feature/60,1)}.png"
        plt.savefig(os.path.join(savedir,figname))
        if show_plot:
            plt.show()
        else:
            plt.close()

def plot_EIC(m_z_list,feature_RT,RT_orig_list,score_list,p_scores,areas,directory,savedir,
             weight=True,show_plot=False,time=5,peak_range=30,decimal=2,reps=3,baseline=0):
    RTs_full = []
    ints_full = []
    mzs_full = []
    names = []
    for file in sorted(os.listdir(directory)):
        if file.endswith('.mzML'):
            names.append(file[:-5])
            exp = MSExperiment()
            MzMLFile().load(os.path.join(directory,file),exp)
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
    for j,mz in enumerate(m_z_list):       # loop through each feature
        RT_feature = feature_RT[j]
        RT_origs = RT_orig_list[j]
        area_feature = areas[j]
        e_score = score_list[j]
        p_val = p_scores[j]
        xs = []
        xs_zoom = []
        ys = []
        max_array = []
        for k in range(len(RTs_full)):             # loop through each replicate
            RTs_rep = RTs_full[k]
            RT_orig = RT_origs[k]
            if RT_orig != 0:
                RT_use = RT_orig
            else:
                RT_use = RT_feature
            ints_rep = ints_full[k]
            mzs_rep = mzs_full[k]
            
            RT_idx_low = find_nearest(RTs_rep,RT_use-time/2*60)[1]
            RT_idx_high = find_nearest(RTs_rep,RT_use+time/2*60)[1] 

            RTs_slice = RTs_rep[RT_idx_low:RT_idx_high+1]
            
            max_idx_low = find_nearest(RTs_slice,RT_use-peak_range/2)[1]
            max_idx_high = find_nearest(RTs_slice,RT_use+peak_range/2)[1]

            ints_slice = ints_rep[RT_idx_low:RT_idx_high+1]    # should take a slice of arrays for both of these
            mzs_slice = mzs_rep[RT_idx_low:RT_idx_high+1]
            mz_idx = [find_nearest_tol(entry,mz,10**(-decimal)/2)[1] for entry in mzs_slice] # find idx where the feature mz is in each scan

            ints_EIC = [array[mz_idx[o]] if mz_idx[o] > 0 else baseline for o,array in enumerate(ints_slice)]

            try:
                max_found = np.amax(ints_EIC[max_idx_low:max_idx_high])
            except ValueError:
                max_found = np.amax(ints_EIC)
                print(ints_EIC[max_idx_low:max_idx_high],ints_EIC)
            max_array.append(max_found)
            xs.append(RTs_slice)
            xs_zoom.append(RTs_slice[max_idx_low:max_idx_high])
            ys.append(ints_EIC)
        fig, axs = plt.subplots(int(len(names)/reps),reps,figsize=(14,16),sharey=True)
        for l in range(int(len(names)/reps)):
            for m in range(reps):
                x = [t/60 for t in xs[reps*l+m]]
                x_zoom = [t/60 for t in xs_zoom[reps*l+m]]
                y = ys[reps*l+m]

                axs[l,m].plot(x,y, label = f'Feature mean RT = {np.round(RT_feature/60,1)} min')
                axs[l,m].plot(x_zoom,np.ones(len(x_zoom))*max_array[reps*l+m],'k--',label=f'Area of {area_feature[reps*l+m]}')
                axs[l,m].set_title(f"EIC of m/z = {np.round(mz,decimal)} in {names[reps*l+m][:-5]}",fontsize=14)
                axs[l,m].set_xlabel('Retention time (min)',fontsize=12)
                axs[l,m].set_ylabel('Intensity',fontsize=12)
                yfmt = mticker.ScalarFormatter(useMathText=True)
                yfmt.set_powerlimits((3, 4))
                axs[l,m].yaxis.set_major_formatter(yfmt)
                axs[l,m].ticklabel_format(style='sci', axis='y', scilimits=(0,0))
                axs[l,m].get_yaxis().get_offset_text().set_visible(False)
                ax_max = max(axs[l,m].get_yticks())
                exponent_axis = np.floor(np.log10(ax_max)).astype(int)
                axs[l,m].annotate(r'$\times$10$^{%i}$'%(exponent_axis),
                             xy=(.01, .8), xycoords='axes fraction')
                axs[l,m].legend(loc='best',fontsize=10)
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        if weight:
            figname = f"mz of {np.round(mz,decimal)} at {np.round(RT_feature/60,1)}, escore of {np.round(e_score,2)}, pscore of {p_val:.3E}.png"
        else:
            figname = f"NORMALIZED mz of {np.round(mz,decimal)} at {np.round(RT_feature/60,1)}, escore of {np.round(e_score,2)}, pscore of {p_val:.3E}.png"
        plt.savefig(os.path.join(savedir,figname))
        if show_plot:
            plt.show()
        else:
            plt.close()






