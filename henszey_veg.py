from scipy.io import loadmat
from scipy import special
import numpy as np
import matplotlib as mp
import matplotlib.pyplot as plt


def m8035(L7h,a,b,c):
    if c > 0:
        X = ((L7h-b)/c)**2
        f = a*scipy.special.erfc(X)
        return(f)
    else:
        return('error: c <= 0') 
    

def m8036(L7h,a,b,c):
    if c > 0:
        if L7h >= b:
            X = (L7h-b)/c
            eX = np.e**(-X)  
            f = 4*a*eX*(1-eX)
            return(f)
        else:
            return('error: L7h < b')
    else:
        return('error: c <= 0') 


def m8036(L7h,a,b,c,d):
    if c != 0:
        if d != 0:
            numer = L7h - c*ln(2**(1/d)-1) - b
            f = a / ( (1 + np.e**(numer/c))**d )
            return(f)
        else:
            return('error: c = 0')
    else:
        return('error: d = 0') 
        

spp = ['Symphyotrichum lanceolatum','Sorghastrum nutans','Ambrosia psilostachya']
wells = ['GroundwaterDataA','GroundwaterDataB','GroundwaterDataC']

for well in wells:
    
    # assume ground elevation to convert depth of water table to GW 
    ground_z = 0.0

    # files to save
    well_fig = r'D:\OneDrive - Tulane University\RCSE-6900 CUAHSI VU\ecohydrology_groundwater\HW_1\%s_GWelev.png' % well

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
                print('NaN data found at timestep: %0.4f. Skipping.' % (day+dc_day))
        except:
            print('string-formatted data found at timestep: %0.4f. Skipping.' % (day+dc_day))
         
    del(gwts)

    day0 = min(day_dt.keys())
    dayn = max(day_dt.keys())

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
        if day < dayn - 6:
            z_sum = 0
            z_cnt = 0
            for daymn in range(0,7):
                z_sum += np.mean(day_gwz[day+daymn])
                z_cnt += 1
            if z_cnt > 0:
                t_day.append(day)
                z_7d.append(z_sum/z_cnt)



    fig = plt.figure(figsize=(7,5))
    ax = fig.add_subplot(111,facecolor='whitesmoke')
    ax.plot(t_all,z_all,linestyle='',marker='.',color='darkgray')
    ax.plot(t_day,z_7d,color='black')


#ymin = max(0.0,min(y2plt)*1.1)
#ymax = max(y2plt)*0.9
#ax.set_ylim(ymin,ymax)

    ax.set_xlabel('Day of year')
    ax.set_ylabel('Water table elevation (relative to ground surface) [L]')
    plt.title('%s' % key)
    plt.show()
    #plt.savefig(well_fig)

 


    
