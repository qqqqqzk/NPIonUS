import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

fit = pd.read_csv('matfit0528.csv')
se = pd.read_csv('matse0528.csv')
cumfit = pd.read_csv('cumfit0528.csv')
cumse = pd.read_csv('cumse0528.csv')
dlnmall = pd.read_csv('all0528.csv')

npi = pd.DataFrame(['Activate or establish emergency response',
                    'Actively communicate with healthcare professionals',
                    'Actively communicate with managers',
                    'Adapt procedures for patient management',
                    'Closure of educational institutions', 'Crisis management plans',
                    'Educate and actively communicate with the public',
                    'Enhance laboratory testing capacity', 'Increase availability of PPE',
                    'Increase healthcare workforce',
                    'Increase in medical supplies and equipment',
                    'Increase patient capacity', 'Indoor and outdoor gathering restriction',
                    'Indoor gathering restriction', 'Measures for special populations',
                    'Measures to ensure security of supply', 'National lockdown',
                    'Outdoor gathering restriction', 'Police and army interventions',
                    'Special measures for certain establishments', 'Work safety protocols'], columns=['NPI'])

fit = pd.concat([npi, fit.iloc[:,1:]], axis=1)
se = pd.concat([npi, se.iloc[:,1:]], axis=1)
cumfit = pd.concat([npi, cumfit.iloc[:,1:]], axis=1)
cumse = pd.concat([npi, cumse.iloc[:,1:]], axis=1)
dlnmall = pd.concat([npi, dlnmall.iloc[:,1:]], axis=1)

fit.columns = ['NPI'] + fit.iloc[:,1:].columns.str[3:].tolist()
se.columns = ['NPI'] + se.iloc[:,1:].columns.str[3:].tolist()
cumfit.columns = ['NPI'] + cumfit.iloc[:,1:].columns.str[3:].tolist()
cumse.columns = ['NPI'] + cumse.iloc[:,1:].columns.str[3:].tolist()

lacklag = dict(zip(['National lockdown','Outdoor gathering restriction','Work safety protocols','Increase patient capacity'], [7,9,13,15]))
for k,v in lacklag.items():
    fit.loc[fit.NPI==k, [str(i) for i in np.arange(v-1, 15)]]=np.nan
    cumfit.loc[cumfit.NPI==k, [str(i) for i in np.arange(v-1, 15)]]=np.nan

fit['overall'] = fit.iloc[:, 1:16].T.sum()
se['overall'] = fit.iloc[:, 1:16].T.sum()
fit = fit.sort_values('overall')[::-1].reset_index(drop=True)
se = se.sort_values('overall')[::-1].reset_index(drop=True)
fit = fit.drop('overall', axis=1)
se = se.drop('overall', axis=1)

# Figure 1
fig, ax = plt.subplots(figsize=(10,6))
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
plt.grid(alpha=0.3, axis='y', ls='--')
plt.barh(cumfit.loc[:, 'NPI'][-10:], cumfit.loc[:,'14'][-10:], color='#4b4b4b')
plt.rcParams['font.sans-serif'] = 'Times new roman'
for x,y in zip(np.arange(10), cumfit.loc[:,'14'][-10:]):
    plt.text(y-0.14, x, '%.2f'%y, color='w', fontsize=22, va='center', ha='center')
plt.xticks(np.arange(0, 3, 0.5), fontsize=25)
ax.set_yticklabels([str(i)+' ' for i in np.arange(10,0,-1)] + cumfit.loc[:, 'NPI'][-10:].values, fontdict={'fontsize':25})
plt.xlabel('Cumulative change of      (14 Days)', {'size':28}); plt.ylabel('NPIs', {'size':28}, labelpad=20)
plt.text(1.66,-2.45,'Rt', fontstyle='italic', fontsize=28, va='center')
plt.gca().yaxis.set_ticks_position('none');
plt.savefig('Figure 1.pdf', bbox_inches='tight');

# Figure 2
# GAM # <5
plt.rcParams['font.sans-serif'] = 'Times new roman'
fig, ax = plt.subplots(5,2, figsize=(25/2,30/2)); n, m=0, 0 
plt.subplots_adjust(left=None, bottom=None, right=None, top=None,wspace=None, hspace=0.5) 
lacklag = dict(zip(['National lockdown','Outdoor gathering restriction','Work safety protocols','Increase patient capacity'], [7,9,13,15]))
for i, sign in zip([8,2,5,9,3,6,0,1,4,7], [chr(i).upper() for i in range(97,97+10)]): # 31 18 13 
    if fit.iloc[i,0] in lacklag.keys(): 
        ax[n, m].plot(fit.iloc[i,1:lacklag[fit.iloc[i,0]]], c='black')#, label=fit.iloc[i,0]); ax[n,m].legend(loc='upper left', bbox_to_anchor=(-0.02,1.15), fontsize=14, markerfirst=False, markerscale=0.01)
        ax[n, m].fill_between(np.arange(lacklag[fit.iloc[i,0]]-1), 
                              (fit.iloc[i, 1:lacklag[fit.iloc[i,0]]]-(1.96*se.iloc[i, 1:lacklag[fit.iloc[i,0]]])).astype('float64'), 
                              (fit.iloc[i, 1:lacklag[fit.iloc[i,0]]]+(1.96*se.iloc[i, 1:lacklag[fit.iloc[i,0]]])).astype('float64'), 
                              color='gray', alpha=0.1)        
        ax[n, m].text(-lacklag[fit.iloc[i,0]]/5.5, (fit.iloc[i, 1:lacklag[fit.iloc[i,0]]]+(1.96*se.iloc[i, 1:lacklag[fit.iloc[i,0]]])).max(), sign, fontsize=25, fontweight='bold', va='center_baseline', ha='center')
    else:
        ax[n, m].plot(fit.iloc[i,1:], c='black')#, label=fit.iloc[i,0]); ax[n,m].legend(loc='upper left', bbox_to_anchor=(-0.02,1.15), fontsize=14, markerfirst=False)
        ax[n, m].fill_between(np.arange(15), 
                              (fit.iloc[i, 1:]-(1.96*se.iloc[i, 1:])).astype('float64'), 
                              (fit.iloc[i, 1:]+(1.96*se.iloc[i, 1:])).astype('float64'), 
                              color='gray', alpha=0.1)
        ax[n, m].text(-3, (fit.iloc[i, 1:]+(1.96*se.iloc[i, 1:])).max(), sign, fontsize=25, fontweight='bold', va='center_baseline', ha='center')
    ax[n, m].axhline(0, color='gray', ls='--')
    ax[n, m].tick_params(labelsize=18, rotation=0)
    ax[n, m].set_xlabel(fit.iloc[i,0], fontsize=20); ax[n, m].set_ylabel('ΔRt', fontsize=20, fontstyle='italic');
    #ax[n, m].grid(axis='x', alpha=0.1)
    ax[n, m].axvline(fit.iloc[[i],1:].T.idxmax().values[0], alpha=0.8, color='gray', ls='--') 
    take_effect = (fit.iloc[i,1:] - 1.96*se.iloc[i,1:]).dropna()
    ax[n, m].axvline(take_effect.loc[take_effect>0].index[0], alpha=0.8, color='gray', ls='--')
    m+=1
    if m==2: m=0; n+=1
plt.tight_layout()
plt.savefig('Figure 2.pdf', bbox_inches='tight'); 

f = pd.read_csv('state0520.csv')
f.date = pd.to_datetime(f.date)

f = f.loc[:, ['date', 'inc_2020', 'inc_p', 'Rt_2020', 'Rt_p', 'dif_Rt', 'State', 
          'Increase patient capacity', 'Police and army interventions',        'Outdoor gathering restriction',
          'Special measures for certain establishments',        'Increase in medical supplies and equipment',
          'Measures for special populations',        'Actively communicate with healthcare professionals',
          'Adapt procedures for patient management', 'National lockdown',        'Activate or establish emergency response']]
f.iloc[:, 7:] = f.iloc[:, 7:].T.cumsum().T * f.iloc[:, 7:]
f = f.replace(0, np.nan)

import matplotlib as mpl
mpl.rcParams['legend.handlelength']=0
plt.rcParams['font.sans-serif'] = 'Times new roman'

fig, ax = plt.subplots(10, 2, figsize=(18/2,30/2*2.5))
plt.subplots_adjust(left=None, bottom=None, right=None, top=None,wspace=None, hspace=None)
i=0; ii=0
for iii, num, sign in zip(f.groupby('State').dif_Rt.sum().sort_values().index, 
                     f.groupby('State').dif_Rt.sum().sort_values().values.round(2), 
                     [chr(i).upper() for i in range(97,97+19)]):
    temp_df = f.loc[f.State==iii]
    ax[i,ii].spines['right'].set_color('none')
    ax[i,ii].spines['top'].set_color('none')
    ax[i,ii].plot(temp_df.loc[:,'date'], temp_df.loc[:,'dif_Rt'], label=iii+' ({})'.format(num), c='black');
    ax[i,ii].legend(loc='upper right', fontsize=14)
    for j in f.columns[7:]:
        for k,l in zip(temp_df.loc[:,'date'][temp_df.loc[:,j].notna().values], temp_df.loc[temp_df[j].notna(),j]):
            ax[i,ii].text(k, 
                          min(temp_df.loc[:,'dif_Rt'].min(), 0) + (temp_df.loc[:,'dif_Rt'].max() - min(temp_df.loc[:,'dif_Rt'].min(), 0))/16*(l-1), 
                          f.columns[7:].tolist().index(j)+1, 
                          rotation=0, fontsize=15, va='center', ha='center', fontweight='bold')
    ax[i,ii].axhline(0, color='gray', ls='--')
    ax[i,ii].set_xticks(pd.date_range(temp_df.date.iloc[0], temp_df.date.iloc[-1], 5))
    month_eng = dict(zip([2,3,4], ['Feb','Mar','Apr']))
    ax[i,ii].set_xticklabels([str(i.day)+'-'+str(month_eng[i.month]) for i in pd.date_range(temp_df.date.iloc[0], temp_df.date.iloc[-1], 5)], fontsize=14)
    ax[i,ii].tick_params(labelsize=14)
    ax[i,ii].set_xlabel('Date', fontsize=14); ax[i,ii].set_ylabel('ΔRt', fontsize=14, fontstyle='italic')
    ax[i,ii].text(temp_df.iloc[0, 0] - pd.tseries.offsets.Day(int(temp_df.shape[0]/3.2)), temp_df.dif_Rt.max(), sign, fontsize=26, fontweight='bold')
    ax[i,ii].fill_between(temp_df.loc[:,'date'], temp_df.loc[:,'dif_Rt'], 0, color='gray', alpha=0.15)
    ii+=1
    if ii==2: ii=0; i+=1
ax[i, ii].axis('off')
plt.tight_layout()
plt.savefig('Figure 3.pdf', bbox_inches='tight'); 