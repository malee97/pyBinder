import os
import platform
import csv
import math
import shutil 
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

def plot_spectra_2D_overview(exp,title,other):
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

def feature_int_extractor(m_z_feature_list,RT_feature_list,RTs_orig_list,mzs,RTs,ints,LOD=1.5E4,noise_level=150,peak_range=120,decimal=2,percent=1,subtract=False):
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

def feature_int_extractor_start(m_z_feature_list,RT_feature_list,mzs,RTs,ints,LOD=1.5E4,noise_level=150,peak_range=120,decimal=2,percent=1,subtract=False):
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

def peak_finder_savgol(rt_windows,ints_windows,names,plot=False,reps=3,default_thresh=1.5E5,kernal_size=10,width=5,prominence=2,threshold=2,rel_height=0.5):
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

def feature_area_extractor_savgol(rt_windows_filt,int_windows_filt,check=False,width_start=5,prominence=2,threshold=2,rel_height=0.5):
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

def plot_PRTC(m_z_list,feature_RT,RT_orig_list,areas,directory,savedir,weight=True,show_plot=False,time=5,peak_range=30,decimal=2,reps=3,baseline=0):
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

def plot_EIC(m_z_list,feature_RT,RT_orig_list,score_list,p_scores,areas,directory,savedir,weight=True,show_plot=False,time=5,peak_range=30,decimal=2,reps=3,baseline=0):
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






