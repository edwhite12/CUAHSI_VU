import scipy
from scipy.io import loadmat
from scipy import special
import math
import numpy as np
import matplotlib as mp
import matplotlib.pyplot as plt


def m8035(L7H_cm,a,b,c):
    # curve-fitting Model 8035 from Henszey et al., 2004
    L7H = -L7H_cm     # models were fit with groundwater elevation as positive values
    if c > 0:
        X = ((L7H-b)/c)**2
        f = max(0,a*scipy.special.erfc(X))
        return(f)
    else:
        return('error: c <= 0') 
    

def m8036(L7H_cm,a,b,c):
    # curve-fitting Model 8036 from Henszey et al., 2004
    L7H = -L7H_cm     # models were fit with groundwater elevation as positive values
    if c > 0:
       # if abs(L7H_cm) >= abs(b):
        X = (L7H-b)/c
        eX = np.e**(-X)  
        f =  max(0,4*a*eX*(1-eX))
        return(f)
      #  else:
      #S      return('error: L7H < b')
    else:
        return('error: c <= 0') 


def m8089(L7H_cm,a,b,c,d):
    # curve-fitting Model 8089 from Henszey et al., 2004
    L7H = -L7H_cm     # models were fit with groundwater elevation as positive values
    if c != 0:
        if d != 0:
            c3 = L7H - c*(math.log( (2**(1/d)) - 1 )) - b
            c4 = (np.e**(-1*c3/c))
            f =  max(0,a / ((1 + c4)**d))
            return(f)
        else:
            return('error: c = 0')
    else:
        return('error: d = 0') 

wells = ['GroundwaterDataA','GroundwaterDataB','GroundwaterDataC']
sp1 = 'Symphyotrichum lanceolatum'
sp2 = 'Sorghastrum nutans'
sp3 = 'Ambrosia psilostachya'
spp = [sp1,sp2,sp3] 

# plot plant frequency curves
L7H_list = []
pf1_list = []
pf2_list = []
pf3_list = []

for L7H_cm in range(-200,50,1):
    # coefficients for each species from Table 2 in Henszey et al., 2004
    pf1_list.append(m8036(L7H_cm,2.3,-17.2,24.1))
    pf2_list.append(m8035(L7H_cm,3.8,70.2,61.8))
    pf3_list.append(m8089(L7H_cm,4.7,113.4,0.6,0.1))
    L7H_list.append(L7H_cm)

pf_lists = [pf1_list,pf2_list,pf3_list]
for sp in spp:
    pf_list = pf_lists[spp.index(sp)]
    sp_fig = r'D:\OneDrive - Tulane University\RCSE-6900 CUAHSI VU\ecohydrology_groundwater\HW_2\%s.png' % sp
    fig = plt.figure(figsize=(7,5))
    ax = fig.add_subplot(111,facecolor='whitesmoke')
    ax.plot(L7H_list,pf_list,color='black')
    ax.set_xlim(50,-200)
    ax.set_ylim(0,5)
    ax.set_xlabel('7-day Moving Average High Water (cm)')
    ax.set_ylabel('Plant Frequency (%)')
    plt.title('%s' % sp,style='italic')
    #plt.show()
    plt.savefig(sp_fig)


for well in wells:
    # assume ground elevation to convert depth of water table to GW 
    ground_z = 0.0

    # files to save
    well_fig = r'D:\OneDrive - Tulane University\RCSE-6900 CUAHSI VU\ecohydrology_groundwater\HW_2\%s_GWelev.png' % well

    #files to read
    file = r'D:\OneDrive - Tulane University\RCSE-6900 CUAHSI VU\ecohydrology_groundwater\HW_1\%s.mat' % well
    key = well

    # import .mat file using loadmat
    gwts = loadmat(file)[key]     

    # read data into dictionaries with integer day as the dict key
    # during the first timestep of any day the 'try' statement will fail
    # the 'except' statement will initialize an empty list for the day
    # once the empty list exists for the day, the individual datapoints and timesteps are appended 
    day_dt = {}
    day_gwz = {}

    for row in gwts:
        day = int(row[0])
        dc_day = row[0]-day
        gw_d = row[1]
        try:
            checknum = float(gw_d) # will check if input values are strings
            if np.isnan(float(gw_d)) == False:
                try:
                    day_dt[day].append(dc_day)
                    day_gwz[day].append(ground_z - gw_d)    # convert depth to groundwater to watertable elevation
                except:                 
                    day_dt[day] = []
                    day_gwz[day] = []
                    day_dt[day].append(dc_day)
                    day_gwz[day].append(ground_z - gw_d)
            else:
                msg = 'NaN data found at timestep: %0.4f. Skipping.' % (day+dc_day)
        except:
            msg = 'string-formatted data found at timestep: %0.4f. Skipping.' % (day+dc_day)
         
    del(gwts)

    day0 = min(day_dt.keys())
    dayn = max(day_dt.keys())

    # calculate 7-day moving window average groundwater depth
    t_all = []
    z_all = []
    t_day = []
    z_7d = []

    for day in day_dt.keys():
        for nt in range(0,len(day_dt[day])):
            dec_day = day_dt[day][nt]
            gwz_nt = day_gwz[day][nt]    

            t_all.append(day+dec_day)
            z_all.append(gwz_nt)
        if day > day0 + 6:
            z_sum = 0
            z_cnt = 0
            for daymn in range(0,7):
                z_sum += np.mean(day_gwz[day-daymn])
                z_cnt += 1
            if z_cnt > 0:
                t_day.append(day)
                z_7d.append(z_sum/z_cnt)

    L7L_cm = min(z_7d)*100.
    L7H_cm = max(z_7d)*100.
    footnote = 'L7L = %0.2f cm; L7H = %0.2f cm' % (L7L_cm,L7H_cm)
    fig = plt.figure(figsize=(7,5))
    ax = fig.add_subplot(111,facecolor='whitesmoke')
    ax.plot(t_all,z_all,linestyle='',marker='.',markersize=2,color='darkgray')
    ax.plot(t_day,z_7d,color='black')
    ax.set_xlabel('Day of year')
    ax.set_ylabel('Water table elevation (relative to ground surface) (m)')
    plt.title('%s' % key)
    plt.figtext(0.99,0.01,footnote,horizontalalignment='right',fontsize='small')
    #plt.show()
    plt.savefig(well_fig)
    

    pf1 = m8036(L7H_cm,2.3,-17.2,24.1)
    pf2 = m8035(L7H_cm,3.8,70.2,61.8)
    pf3 = m8089(L7H_cm,4.7,113.4,0.6,0.1)

    if well == 'GroundwaterDataA':
        print('well,%s,%s,%s,L7H (cm)' % (sp1,sp2,sp3))
    print('%s,%0.4f,%0.4f,%0.4f,%0.4f' % (well,pf1,pf2,pf3,L7H_cm))






    
