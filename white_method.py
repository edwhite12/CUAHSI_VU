from scipy.io import loadmat
import numpy as np
import matplotlib as mp
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

# assume ground elevation to convert depth of water table to GW 
ground_z = 0.0

# files to save
well_fig = r'D:\OneDrive - Tulane University\RCSE-6900 CUAHSI VU\ecohydrology_groundwater\HW_1\GroundwaterDataC_GWelev.png'
ETg_fig = r'D:\OneDrive - Tulane University\RCSE-6900 CUAHSI VU\ecohydrology_groundwater\HW_1\GroundwaterDataC_ETg.png'

# files to read
file = r'D:\OneDrive - Tulane University\RCSE-6900 CUAHSI VU\ecohydrology_groundwater\HW_1\GroundwaterDataC.mat'
key = 'GroundwaterDataC'

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

# plot well hydrograph
x2plt = []
y2plt = []

for day in day_dt.keys():
    for nt in range(0,len(day_dt[day])):
        dec_day = day_dt[day][nt]
        gwz_nt = day_gwz[day][nt]
        x2plt.append(day+dec_day)
        y2plt.append(gwz_nt)


fig = plt.figure(figsize=(7,5))
ax = fig.add_subplot(111,facecolor='whitesmoke')
ax.plot(x2plt,y2plt,color='blue')


#ymin = max(0.0,min(y2plt)*1.1)
#ymax = max(y2plt)*0.9
#ax.set_ylim(ymin,ymax)

ax.set_xlabel('Day of year')
ax.set_ylabel('Water table elevation (relative to ground surface) [L]')
plt.title('%s' % key)
#plt.show()
plt.savefig(well_fig)

 
# determine change in storage and recharge rate to estimate ETg from method outlined in White (1932)
# ETg = SY(dS + t*R)
# ETg will be calculated on day0 up to dayn-1
day0 = min(day_dt.keys())
dayn = max(day_dt.keys())

# check if there are any missing days in the record
if len(day_dt.keys()) < (dayn-day0):
       print(' Length of record check failed. Looks like there may be some missing days in the record. Script may run but could potentially produce errors.')
       input(' Press <ENTER> to acknowledge and continue run.')

# check if last day has a full daily record (here defined as recording at least until 11 pm)
# if not, exclude from daily looping below
last_hr = max(day_dt[dayn])
if last_hr < 23./24.:
    dayn -= 1



# start daily loop to calculate dS and R; in Python, range is not inclusive on upper end, so last iteration will be dayn-1
dS = {}
R = {}
for di in range(day0,dayn):
    # calculate change in storage for each day
    max_di = max(day_gwz[di])
    max_di1 = max(day_gwz[di+1])
    dS[di] = max_di - max_di1       # dS>0 indicates lowering water table

    # build empty lists for the day that includes only data from 0:00 through 4:00, inclusive 
    gw4R = []
    dt4R = []
    
    for nt in range(0,len(day_dt[di])):
        dec_day = day_dt[di][nt]
        gwz_nt = day_gwz[di][nt]
        gw4R.append(gwz_nt)
        dt4R.append(dec_day)

        # check if timestep is the closest to 4:00 - if not exactly measured at 4:00, it will go to the first data point recorded after 4:00
        if dec_day >= 4./24.:
            break

    # restructure lists as 1xN and Nx1 matrices for use in regression model
    gw4R = np.asarray(gw4R)
    dt4R = np.asarray(dt4R).reshape(-1,1)
    LRmod = LinearRegression()  

    # fit the model and save slope to daily recharge rate dict
    LRmod.fit(dt4R,gw4R)
    R[di] = LRmod.coef_[0]          # R>0 indicates groundwater recharge
    #intcpt = LRmod.intercept_

# calculate ETg for a variety of Specific Yield values
# ETg = SY(dS + t*R)
t = 1.0 # day
SY_md = 0.1
SYu = 0.1       #uncertainty range (+/-) to apply to SY estimates
SY_lo = (1.0 - SYu)*SY_md
SY_hi = (1.0 + SYu)*SY_md

footnote = 'Gray bars indicate range in ETg from a +/- %s range in specific yield.' % ('{:.0%}'.format(SYu))

ETg_lo = []
ETg_md = []
ETg_hi = []
x2plt = []

for day in dS.keys():
    if dS[day] > 0:                                         # check that storage has gone down from one day to next
        if R[day] > 0:                                      # check that recharge has occurred over night
            x2plt.append(day)
            ETg_lo.append( SY_lo*(dS[day] + t*R[day]) )    
            ETg_md.append( SY_md*(dS[day] + t*R[day]) )
            ETg_hi.append( SY_hi*(dS[day] + t*R[day]) )
        else:
            print('Negative (or zero) recharge rate on morning of day %s. No ETg calculated.' % day)
    else:
        print('Groundwater storage increased from day %s to day %s. No ETg calculated.' % (day,day+1))

fig = plt.figure(figsize=(7,5))
ax = fig.add_subplot(111,facecolor='whitesmoke')
ax.plot(x2plt,ETg_lo,linestyle='',marker='_',color='darkgray')
ax.plot(x2plt,ETg_hi,linestyle='',marker='_',color='darkgray')
ax.plot(x2plt,ETg_md,linestyle='',marker='.',color='black')

ax.set_xlabel('Day of year')
ax.set_ylabel('ETg [L]/day')
plt.title('%s' % key)
plt.figtext(0.99,0.01,footnote,horizontalalignment='right',fontsize='xx-small')
#plt.show()
plt.savefig(ETg_fig)

    
