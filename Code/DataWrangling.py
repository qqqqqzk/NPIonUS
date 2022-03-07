import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#1 ILI rate
ILI= pd.read_csv('ILINet.csv', na_values='X', skiprows=1)
# exclude before 2015
ILI = ILI.loc[13886:, ['REGION', 'YEAR', 'WEEK', '%UNWEIGHTED ILI']].loc[ILI['%UNWEIGHTED ILI'].notna()].reset_index(drop=True)
ILI.rename(columns={'%UNWEIGHTED ILI':'ILI'}, inplace=True)

#2 Lab 
lab = pd.read_csv('WHO_NREVSS_Clinical_Labs.csv', na_values='X', skiprows=1)
# exclude missing value ratio > 0.5
selected = ['Idaho', 'Maine', 'North Dakota', 'Delaware', 'Mississippi', 'Vermont', 'New Mexico', 'Oregon', 'Tennessee',
            'South Carolina', 'Utah', 'Kansas', 'Maryland', 'North Carolina',
            'Hawaii', 'Iowa', 'Arkansas', 'Nebraska', 'Pennsylvania', 'Connecticut', 'Oklahoma', 
            'Texas', 'Missouri', 'California', 'South Dakota', 'Massachusetts', 
            'Montana', 'New York', 'Wisconsin', 'Virginia', 'Louisiana', 'Georgia', 
            'Colorado', 'Washington', 'Minnesota', 'Ohio', 'Michigan', 'West Virginia', 
            'Illinois', 'Kentucky', 'Arizona', 'Indiana', 'Alabama']
lab = lab.loc[lab.REGION.isin(selected), ['REGION', 'YEAR', 'WEEK', 'PERCENT POSITIVE', 'PERCENT A', 'PERCENT B']]

#3 Incidence
inc = pd.merge(lab, ILI, on=['REGION', 'YEAR', 'WEEK'], how='left')
inc['incidence'] = inc['PERCENT POSITIVE'] * inc['ILI']

#4 NPI
npi = pd.read_excel('usnpi.xlsx')
npi_l2 = npi[['State', 'Date', 'Measure_L2']].drop_duplicates() # select L2 level
npi_l2['va'] = 1 # prepare for pivot
# pivot transform to matrix
npi_l2 = npi_l2.pivot(index=['State','Date'], columns='Measure_L2', values='va').fillna(0).reset_index()
# only retain the state has both the npi and flu data
npi_l2 = npi_l2.loc[npi_l2.State.isin(set(npi_l2.State.unique())&set(inc.REGION.unique()))]
# set datetime as index
npi_l2.set_index(pd.to_datetime(npi_l2.Date), drop=False, inplace=True) 

# fill the full time series
npi_l2_ts = pd.DataFrame()
for i in npi_l2.State.unique():
    temp_df = npi_l2.loc[npi_l2.State==i]
    temp_df = temp_df.reindex(pd.date_range(temp_df.index.min(), temp_df.index.max())) # reindex补时间序列缺失的值
    temp_df.State.fillna(i, inplace=True)
    npi_l2_ts = pd.concat([npi_l2_ts, temp_df.fillna(0).drop('Date', axis=1)])
# delete the state without npi data, greater than 2020
inc = inc.loc[(inc.REGION.isin(set(npi_l2.State.unique())&set(inc.REGION.unique())))&(inc.YEAR<2021)]
# spline fill the incidence
inc_s = pd.DataFrame()
from scipy import interpolate
fig, ax = plt.subplots(21,1,figsize=(5,5*21)); n=0
for i in inc.REGION.unique():
    temp_var = inc.loc[inc.REGION==i].reset_index()
    inter = interpolate.CubicSpline(np.delete(np.arange(3, temp_var.shape[0]*7, 7), temp_var.loc[(temp_var.incidence==0)|(temp_var.incidence.isna())].index), 
                                    temp_var.loc[(temp_var.incidence!=0)&(temp_var.incidence.notna()), 'incidence'].values, 
                                    bc_type='natural') 
    xx = np.arange(0, temp_var.shape[0]*7, 1)
    interxx = inter(xx)
    interxx[interxx<0]=np.nan; interxx[interxx>500]=np.nan
    ax[n].plot(xx, interxx, label=i); ax[n].legend()
    ax[n].scatter(np.delete(np.arange(2, temp_var.shape[0]*7, 7), temp_var.loc[(temp_var.incidence==0)|(temp_var.incidence.isna())].index), 
                temp_var.loc[(temp_var.incidence!=0)&(temp_var.incidence.notna()), 'incidence'].values)
    n+=1
    inc_s = pd.concat([inc_s, pd.DataFrame(interxx[:1820], columns=['inc_splines'])])
inc_ts = pd.DataFrame()
# fill the full time series
for i in inc.REGION.unique():
    temp_df = inc.loc[inc.REGION==i].reindex(pd.date_range('2015-09-28', '2020-09-27'))
    temp_df.REGION.fillna(i, inplace=True)
    inc_ts = pd.concat([inc_ts, temp_df])
inc_ts = inc_ts.loc[~inc_ts.index.isin(pd.date_range('2015-12-28', '2016-1-3'))]
# concat dataset
case = pd.concat([inc_ts.reset_index(), inc_s.reset_index(drop=True)], axis=1)[['index', 'REGION', 'YEAR', 'WEEK', 'incidence', 'inc_splines']]
case['week'] = case['index'].dt.week
case['year'] = case['index'].dt.year
case = case.loc[:, ['index','REGION','week','year','inc_splines']]
# transform dataset
case_avg = pd.DataFrame()
for i in case.REGION.unique():
    region_df = pd.DataFrame()
    #for j, k in zip(np.array([0, 91,364,364,364,364]).cumsum(), np.array([91,364,364,364,364,131]).cumsum()):
    for j, k in zip(np.array([0,364,364,364,364]).cumsum(), np.array([364,364,364,364,364]).cumsum()):
        region_df = pd.concat([region_df, case.loc[case.REGION==i].reset_index(drop=True).iloc[j:k].reset_index(drop=True)], axis=1)
    case_avg = pd.concat([case_avg, region_df])
# average over previous five years
dif = pd.concat([case_avg.iloc[:, [-9,-8,-5,-2,-1]], 
                 pd.DataFrame(case_avg.iloc[:,np.arange(4,5*5,5)].T.mean(), columns=['avg'])], axis=1)
dif['dif'] = dif['avg'] - dif['inc_splines']
dif['ratio'] = (dif['avg'] - dif['inc_splines'])/dif['inc_splines']
# calculte Rt to R
Rt_2020 = pd.read_csv('rt_2020.csv', names=['REGION','Rt'], skiprows=1, usecols=[1,2])
Rt_p = pd.read_csv('rt_p.csv', names=['REGION','Rt'], skiprows=1, usecols=[1,2])
# concat dataset
df = pd.DataFrame()
for i in dif.REGION.unique():
    state = pd.concat([dif.loc[(dif.REGION==i)&(dif.index>6)].reset_index(drop=True), 
                       Rt_2020.loc[(Rt_2020.REGION==i), 'Rt'].reset_index(drop=True), 
                       Rt_p.loc[(Rt_p.REGION==i), 'Rt'].reset_index(drop=True)], axis=1) 
    df = pd.concat([df, state])
df.columns = ['REGION', 'week', 'date', 'year', 'inc_2020', 'inc_p', 'dif_inc', 'ratio_inc', 'Rt_2020', 'Rt_p']
df['dif_Rt'] = df.Rt_p - df.Rt_2020
df['ratio_Rt'] = (df.Rt_p - df.Rt_2020) / df.Rt_2020
# set start time
truncate = pd.DataFrame()
for i in np.arange(dif.REGION.nunique()):
    temp_v = dif.loc[dif.REGION==dif.REGION.unique()[i], 'inc_splines']
    truncate = pd.concat([truncate, 
                          pd.DataFrame([dif.REGION.unique()[i], 
                                        temp_v.loc[temp_v<10].loc[temp_v.loc[temp_v<10].index>100].index[0]], index=['REGION', 'point']).T])
truncate.reset_index(drop=True, inplace=True)
truncate['date'] = pd.Series(pd.date_range('2019-09-30', '2020-05-08'))[truncate.point].reset_index(drop=True)
truncate.date = truncate.date.astype('str')
# truncate time range
df_Rt = pd.DataFrame()
for i in df.REGION.unique():
    df_Rt = pd.concat([df_Rt, df.loc[(df.REGION==i)&(df.date.isin(pd.date_range('2020-02-27', truncate.loc[truncate.REGION==i, 'date'].values[0])))]])
dlnm = pd.merge(df_Rt.loc[:,['REGION', 'date', 'inc_2020','inc_p','Rt_2020','Rt_p','dif_Rt']], npi_l2_ts.reset_index(), left_on=['REGION','date'], right_on=['State','index']).drop(['index','REGION'],axis=1)
# delete the state time range lower than 14 days
dlnm = dlnm.loc[~dlnm.State.isin(['Michigan','Delaware'])]
# lag 2 weeks
fc = pd.DataFrame()
for j in dlnm.State.unique():
    fcc = dlnm.loc[dlnm.State==j]
    for i in np.arange(1, 15):
        fcc = pd.concat([fcc, pd.DataFrame(fcc.shift(i).iloc[:, 3:45]).add_prefix('Lag_{}_'.format(i))], axis=1) # 滞后1~7天 # 45-3=42个政策
    fc = pd.concat([fc, fcc])
# adjust columns order
lag = fc.iloc[:, :3]
for i in np.arange(3, 45):
    lag = pd.concat([lag, fc.iloc[:, i::42]], axis=1)

