import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
import numpy as np

#abrir todos los archivos necesarios
clin_miao = pd.read_csv('/home/vant/Miao/Clinical_Miao.csv',sep='\t') #Datos clínicos recogidos de Miao et al., 2018
clin_forde = pd.read_csv('/home/vant/Forde/Clinical_Forde.csv',sep='\t') #Datos clínicos recogidos de Forde et al., 2018
clin_hugo = pd.read_csv('/home/vant/Hugo/Clinical_Hugo.csv',sep='\t') #Datos clínicos recogidos de Hugo et al., 2016
clin_van = pd.read_csv('/home/vant/Van_allen/Clinical_Van_allen.csv',sep='\t') #Datos clínicos recogidos de Van allen et al., 2015
clin_miao2 = pd.read_csv('/home/vant/Miao2/Clinical_Miao2.csv',sep='\t') #Datos clínicos recogidos de Miao et al., 2018

df_forde = pd.read_csv('/home/vant/Forde/new_NMD_Forde.csv',sep='\t') #dataframe con variantes somáticas con NMD score de Forde et al., 2018
df_hugo = pd.read_csv('/home/vant/Hugo/new_NMD_Hugo.csv',sep='\t') #dataframe con variantes somáticas con NMD score de Hugo et al., 2016
df_van = pd.read_csv('/home/vant/Van_allen/new_NMDB_Van_Allen.csv',sep='\t') #dataframe con variantes somáticas con NMD score de Van allen et al., 2015
df_miao2 = pd.read_csv('/home/vant/Miao2/new_NMD_Miao2.csv',sep='\t') #dataframe con variantes somáticas con NMD score de Miao et al., 2018
df_miao = pd.read_csv('/home/vant/Miao/new_NMD_Miao2.0.csv',sep='\t') #dataframe con variantes somáticas con NMD score de Miao et al., 2018

df_miao = df_miao.loc[df_miao['Variant_Classification'].str.contains('Frame_Shift')]
df_miao = df_miao.loc[df_miao['NMDB'] != '.']
df_forde = df_forde.loc[df_forde['Consequence'].str.contains('Frameshift')]
df_hugo = df_hugo.loc[df_hugo['MutType'].str.contains('Frame_Shift')]
df_hugo = df_hugo.loc[df_hugo['NMDB'] != '.']
df_van = df_van.loc[df_van['Variant_Classification'].str.contains('Frame_Shift')]
df_van = df_van.loc[df_van['NMDB'] != '.']
df_miao2 = df_miao2.loc[df_miao2['pair_id'].str.contains('_T_N')]
df_miao2 = df_miao2.loc[df_miao2['Variant_Classification'].str.contains('Frame_Shift')]
df_miao2 = df_miao2.loc[df_miao2['NMDB'] != '.']

def trigger_nmd(nmd):
    if nmd == '.':
        return 0
    else:
        nmd = float(nmd)
        if nmd >= 0.65:
            return 1
        else:
            return 0

def no_trigger_nmd(nmd):
    if nmd == '.':
        return 0
    else:
        nmd = float(nmd)
        if nmd < 0.65:
            return 1
        else:
            return 0

df_forde['triggersNMD'] = df_forde.apply(lambda x: trigger_nmd(x['NMDB']), axis=1)
df_forde['notriggersNMD'] = df_forde.apply(lambda x: no_trigger_nmd(x['NMDB']), axis=1)
df_hugo['triggersNMD'] = df_hugo.apply(lambda x: trigger_nmd(x['NMDB']), axis=1)
df_hugo['notriggersNMD'] = df_hugo.apply(lambda x: no_trigger_nmd(x['NMDB']), axis=1)
df_van['triggersNMD'] = df_van.apply(lambda x: trigger_nmd(x['NMDB']), axis=1)
df_van['notriggersNMD'] = df_van.apply(lambda x: no_trigger_nmd(x['NMDB']), axis=1)
df_miao2['triggersNMD'] = df_miao2.apply(lambda x: trigger_nmd(x['NMDB']), axis=1)
df_miao2['notriggersNMD'] = df_miao2.apply(lambda x: no_trigger_nmd(x['NMDB']), axis=1)
df_miao['triggersNMD'] = df_miao.apply(lambda x: trigger_nmd(x['NMDB']), axis=1)
df_miao['notriggersNMD'] = df_miao.apply(lambda x: no_trigger_nmd(x['NMDB']), axis=1)

###Miao
df_miao['patient_id'] = df_miao['Tumor_Sample_Barcode'].str.split('-').str[0]
d = dict(zip(list(clin_miao['patient_id']),list(clin_miao['response_category'])))
df_miao['response_category'] = df_miao['patient_id'].map(d)
colors = dict(zip(['clinical benefit','intermediate benefit','no clinical benefit'],['lightgreen', 'lightgreen', 'green']))
df_miao['response_categoryColors'] = df_miao['response_category'].map(colors)
plotdf_miao = df_miao.loc[(df_miao['Variant_Classification'].str.contains('Frame_Shift'))].groupby('Tumor_Sample_Barcode')['triggersNMD'].agg('sum').sort_values(ascending=False).reset_index()
cols = ['response_category', 'response_categoryColors']
for col in cols:
    d = dict(zip(list(df_miao['Tumor_Sample_Barcode']),list(df_miao[col])))
    plotdf_miao[col] = plotdf_miao['Tumor_Sample_Barcode'].map(d)
plotdf2_miao = df_miao.loc[(df_miao['Variant_Classification'].str.contains('Frame_Shift'))].groupby('Tumor_Sample_Barcode')['notriggersNMD'].agg('sum').sort_values(ascending=False).reset_index()
cols = ['response_category', 'response_categoryColors']
for col in cols:
    d = dict(zip(list(df_miao['Tumor_Sample_Barcode']),list(df_miao[col])))
    plotdf2_miao[col] = plotdf2_miao['Tumor_Sample_Barcode'].map(d)

pat = list(set(list(df_miao['Tumor_Sample_Barcode'])))
dfpat = pd.DataFrame(pat)
dfpat['patient_id'] = dfpat[0].str.split('-').str[0]
d = dict(zip(list(clin_miao['patient_id']),list(clin_miao['response_category'])))
dfpat['response_category'] = dfpat['patient_id'].map(d)
dfpat['response_categoryColors'] = dfpat['response_category'].map(colors)
dfpat.rename(columns={0:'Tumor_Sample_Barcode'}, inplace=True)
dfpat_aux = dfpat.copy()
dfpat_aux_2 = dfpat.copy()
dfpatno = dfpat_aux
dfpatyes = dfpat_aux_2
dfpatno['notriggersNMD'] = 0
dfpatyes['triggersNMD'] = 0
mega_pat = list(set(list(plotdf2_miao['Tumor_Sample_Barcode'])))
mega_pat_2 = list(set(list(plotdf_miao['Tumor_Sample_Barcode'])))
dfpatno = dfpatno.loc[~dfpatno['Tumor_Sample_Barcode'].isin(mega_pat)]
dfpatyes = dfpatyes.loc[~dfpatyes['Tumor_Sample_Barcode'].isin(mega_pat_2)]
dfpatyes = dfpatyes.drop(['patient_id'], axis =1)
dfpatno = dfpatno.drop(['patient_id'], axis =1)
plotdf_miao_aux = pd.concat([plotdf_miao, dfpatyes], ignore_index=True)
plotdf2_miao_aux = pd.concat([plotdf2_miao, dfpatno], ignore_index=True)

###Miao2
clin_miao2 = clin_miao2.loc[clin_miao2['pair_id'].str.contains('_T_N')]
d = dict(zip(list(clin_miao2['pair_id']), list(clin_miao2['recist_response'])))
df_miao2['response_category'] = df_miao2['pair_id'].map(d)
colors = dict(zip(['clinical benefit','stable disease','no clinical benefit'],['lightgreen', 'green', 'green']))
df_miao2['response_categoryColors'] = df_miao2['response_category'].map(colors)
plotdf_miao2 = df_miao2.loc[(df_miao2['Variant_Classification'].str.contains('Frame_Shift'))].groupby('pair_id')['triggersNMD'].agg('sum').sort_values(ascending=False).reset_index()
cols = ['response_category', 'response_categoryColors']
for col in cols:
    d = dict(zip(list(df_miao2['pair_id']),list(df_miao2[col])))
    plotdf_miao2[col] = plotdf_miao2['pair_id'].map(d)
plotdf2_miao2 = df_miao2.loc[(df_miao2['Variant_Classification'].str.contains('Frame_Shift'))].groupby('pair_id')['notriggersNMD'].agg('sum').sort_values(ascending=False).reset_index()
cols = ['response_category', 'response_categoryColors']
for col in cols:
    d = dict(zip(list(df_miao2['pair_id']),list(df_miao2[col])))
    plotdf2_miao2[col] = plotdf2_miao2['pair_id'].map(d)

pat = list(set(list(df_miao2['pair_id'])))
dfpat = pd.DataFrame(pat)
dfpat.rename(columns={0:'pair_id'}, inplace=True)
# dfpat['patient_id'] = dfpat[0].str.split('-').str[0]
d = dict(zip(list(clin_miao2['pair_id']),list(clin_miao2['recist_response'])))
dfpat['response_category'] = dfpat['pair_id'].map(d)
dfpat['response_categoryColors'] = dfpat['response_category'].map(colors)
dfpat_aux = dfpat.copy()
dfpat_aux_2 = dfpat.copy()
dfpatno = dfpat_aux
dfpatyes = dfpat_aux_2
dfpatno['notriggersNMD'] = 0
dfpatyes['triggersNMD'] = 0
mega_pat = list(set(list(plotdf2_miao2['pair_id'])))
mega_pat_2 = list(set(list(plotdf_miao2['pair_id'])))
dfpatno = dfpatno.loc[~dfpatno['pair_id'].isin(mega_pat)]
dfpatyes = dfpatyes.loc[~dfpatyes['pair_id'].isin(mega_pat_2)]
# dfpatyes = dfpatyes.drop(['pair_id'], axis =1)
# dfpatno = dfpatno.drop(['pair_id'], axis =1)
plotdf_miao2_aux = pd.concat([plotdf_miao2, dfpatyes], ignore_index=True)
plotdf2_miao2_aux = pd.concat([plotdf2_miao2, dfpatno], ignore_index=True)

###Van_allen
d = dict(zip(list(clin_van['patient']),list(clin_van['group'])))
df_van['response_category'] = df_van['patient'].map(d)
colors = dict(zip(['long-survival','response','nonresponse'],['lightgreen','lightgreen', 'green']))
df_van['response_categoryColors'] = df_van['response_category'].map(colors)
plotdf_van = df_van.loc[(df_van['Variant_Classification'].str.contains('Frame_Shift'))].groupby('patient')['triggersNMD'].agg('sum').sort_values(ascending=False).reset_index()
cols = ['response_category', 'response_categoryColors']
for col in cols:
    d = dict(zip(list(df_van['patient']),list(df_van[col])))
    plotdf_van[col] = plotdf_van['patient'].map(d)
plotdf2_van = df_van.loc[(df_van['Variant_Classification'].str.contains('Frame_Shift'))].groupby('patient')['notriggersNMD'].agg('sum').sort_values(ascending=False).reset_index()
cols = ['response_category', 'response_categoryColors']
for col in cols:
    d = dict(zip(list(df_van['patient']),list(df_van[col])))
    plotdf2_van[col] = plotdf2_van['patient'].map(d)

pat = list(set(list(df_van['patient'])))
dfpat = pd.DataFrame(pat)
d = dict(zip(list(clin_van['patient']),list(clin_van['group'])))
dfpat.rename(columns={0:'patient'}, inplace=True)
dfpat['response_category'] = dfpat['patient'].map(d)
dfpat['response_categoryColors'] = dfpat['response_category'].map(colors)
dfpat_aux = dfpat.copy()
dfpat_aux_2 = dfpat.copy()
dfpatno = dfpat_aux
dfpatyes = dfpat_aux_2
dfpatno['notriggersNMD'] = 0
dfpatyes['triggersNMD'] = 0
mega_pat = list(set(list(plotdf2_van['patient'])))
mega_pat_2 = list(set(list(plotdf_van['patient'])))
dfpatno = dfpatno.loc[~dfpatno['patient'].isin(mega_pat)]
dfpatyes = dfpatyes.loc[~dfpatyes['patient'].isin(mega_pat_2)]
plotdf_van_aux = pd.concat([plotdf_van, dfpatyes], ignore_index=True)
plotdf2_van_aux = pd.concat([plotdf2_van, dfpatno], ignore_index=True)

###Hugo
d = dict(zip(list(clin_hugo['Patient ID']),list(clin_hugo['Response'])))
df_hugo['response_category'] = df_hugo['Sample'].map(d)
colors = dict(zip(['R','NR'],['lightgreen', 'green']))
df_hugo['response_categoryColors'] = df_hugo['response_category'].map(colors)
plotdf_hugo = df_hugo.loc[(df_hugo['MutType'].str.contains('Frame_Shift'))].groupby('Sample')['triggersNMD'].agg('sum').sort_values(ascending=False).reset_index()
cols = ['response_category', 'response_categoryColors']
for col in cols:
    d = dict(zip(list(df_hugo['Sample']),list(df_hugo[col])))
    plotdf_hugo[col] = plotdf_hugo['Sample'].map(d)
plotdf2_hugo = df_hugo.loc[(df_hugo['MutType'].str.contains('Frame_Shift'))].groupby('Sample')['notriggersNMD'].agg('sum').sort_values(ascending=False).reset_index()
cols = ['response_category', 'response_categoryColors']
for col in cols:
    d = dict(zip(list(df_hugo['Sample']),list(df_hugo[col])))
    plotdf2_hugo[col] = plotdf2_hugo['Sample'].map(d)

pat = list(set(list(df_hugo['Sample'])))
dfpat = pd.DataFrame(pat)
d = dict(zip(list(clin_hugo['Patient ID']),list(clin_hugo['Response'])))
dfpat.rename(columns={0:'Sample'}, inplace=True)
dfpat['response_category'] = dfpat['Sample'].map(d)
dfpat['response_categoryColors'] = dfpat['response_category'].map(colors)
dfpat.rename(columns={0:'Sample'}, inplace=True)
dfpat_aux = dfpat.copy()
dfpat_aux_2 = dfpat.copy()
dfpatno = dfpat_aux
dfpatyes = dfpat_aux_2
dfpatno['notriggersNMD'] = 0
dfpatyes['triggersNMD'] = 0
mega_pat = list(set(list(plotdf2_hugo['Sample'])))
mega_pat_2 = list(set(list(plotdf_hugo['Sample'])))
dfpatno = dfpatno.loc[~dfpatno['Sample'].isin(mega_pat)]
dfpatyes = dfpatyes.loc[~dfpatyes['Sample'].isin(mega_pat_2)]
plotdf_hugo_aux = pd.concat([plotdf_hugo, dfpatyes], ignore_index=True)
plotdf2_hugo_aux = pd.concat([plotdf2_hugo, dfpatno], ignore_index=True)

###Forde
d = dict(zip(list(clin_forde['Sample ID']),list(clin_forde['MPR'])))
df_forde['response_category'] = df_forde['Sample ID'].map(d)
colors = dict(zip(['yes','no '],['lightgreen', 'green']))
df_forde['response_categoryColors'] = df_forde['response_category'].map(colors)
plotdf_forde = df_forde.loc[(df_forde['Consequence'].str.contains('Frameshift'))].groupby('Sample ID')['triggersNMD'].agg('sum').sort_values(ascending=False).reset_index()
cols = ['response_category', 'response_categoryColors']
for col in cols:
    d = dict(zip(list(df_forde['Sample ID']),list(df_forde[col])))
    plotdf_forde[col] = plotdf_forde['Sample ID'].map(d)
plotdf2_forde = df_forde.loc[(df_forde['Consequence'].str.contains('Frameshift'))].groupby('Sample ID')['notriggersNMD'].agg('sum').sort_values(ascending=False).reset_index()
cols = ['response_category', 'response_categoryColors']
for col in cols:
    d = dict(zip(list(df_forde['Sample ID']),list(df_forde[col])))
    plotdf2_forde[col] = plotdf2_forde['Sample ID'].map(d)

pat = list(set(list(df_forde['Sample ID'])))
dfpat = pd.DataFrame(pat)
d = dict(zip(list(clin_forde['Sample ID']),list(clin_forde['MPR'])))
dfpat.rename(columns={0:'Sample'}, inplace=True)
dfpat['response_category'] = dfpat['Sample'].map(d)
dfpat['response_categoryColors'] = dfpat['response_category'].map(colors)
dfpat.rename(columns={0:'Sample'}, inplace=True)
dfpat_aux = dfpat.copy()
dfpat_aux_2 = dfpat.copy()
dfpatno = dfpat_aux
dfpatyes = dfpat_aux_2
dfpatno['notriggersNMD'] = 0
dfpatyes['triggersNMD'] = 0
mega_pat = list(set(list(plotdf2_forde['Sample ID'])))
mega_pat_2 = list(set(list(plotdf_forde['Sample ID'])))
dfpatno = dfpatno.loc[~dfpatno['Sample'].isin(mega_pat)]
dfpatyes = dfpatyes.loc[~dfpatyes['Sample'].isin(mega_pat_2)]
plotdf_forde_aux = pd.concat([plotdf_forde, dfpatyes], ignore_index=True)
plotdf2_forde_aux = pd.concat([plotdf2_forde, dfpatno], ignore_index=True)

###test estadísitco
U_miao, p_miao = mannwhitneyu(plotdf2_miao_aux['notriggersNMD'].loc[plotdf2_miao_aux['response_category'].isin(['intermediate benefit','clinical benefit'])], plotdf2_miao_aux['notriggersNMD'].loc[plotdf2_miao_aux['response_category'].isin(['no clinical benefit'])], alternative = 'greater')
U_miao2, p_miao2 = mannwhitneyu(plotdf2_miao2_aux['notriggersNMD'].loc[plotdf2_miao2_aux['response_category'].isin(['clinical benefit'])], plotdf2_miao2_aux['notriggersNMD'].loc[plotdf2_miao2_aux['response_category'].isin(['stable disease','no clinical benefit'])], alternative = 'greater')
U_van, p_van = mannwhitneyu(plotdf2_van_aux['notriggersNMD'].loc[plotdf2_van_aux['response_category'].isin(['long-survival','response'])], plotdf2_van_aux['notriggersNMD'].loc[plotdf2_van_aux['response_category'].isin(['nonresponse'])], alternative = 'greater')
U_hugo, p_hugo = mannwhitneyu(plotdf2_hugo_aux['notriggersNMD'].loc[plotdf2_hugo_aux['response_category'].isin(['R'])], plotdf2_hugo_aux['notriggersNMD'].loc[plotdf2_hugo_aux['response_category'].isin(['NR'])], alternative = 'greater')
U_forde, p_forde = mannwhitneyu(plotdf2_forde_aux['notriggersNMD'].loc[plotdf2_forde_aux['response_category'].isin(['yes'])], plotdf2_forde_aux['notriggersNMD'].loc[plotdf2_forde_aux['response_category'].isin(['no '])], alternative = 'greater')

###BOXPLOT
data_a = [plotdf2_van_aux['notriggersNMD'].loc[plotdf2_van_aux['response_category'].isin(['nonresponse'])], plotdf2_hugo_aux['notriggersNMD'].loc[plotdf2_hugo_aux['response_category'].isin(['NR'])], plotdf2_miao_aux['notriggersNMD'].loc[plotdf2_miao_aux['response_category'].isin(['no clinical benefit'])], plotdf2_miao2_aux['notriggersNMD'].loc[plotdf2_miao2_aux['response_category'].isin(['stable disease','no clinical benefit'])], plotdf2_forde_aux['notriggersNMD'].loc[plotdf2_forde_aux['response_category'].isin(['no '])]]
data_b = [plotdf2_van_aux['notriggersNMD'].loc[plotdf2_van_aux['response_category'].isin(['long-survival','response'])], plotdf2_hugo_aux['notriggersNMD'].loc[plotdf2_hugo_aux['response_category'].isin(['R'])], plotdf2_miao_aux['notriggersNMD'].loc[plotdf2_miao_aux['response_category'].isin(['intermediate benefit','clinical benefit'])], plotdf2_miao2_aux['notriggersNMD'].loc[plotdf2_miao2_aux['response_category'].isin(['clinical benefit'])], plotdf2_forde_aux['notriggersNMD'].loc[plotdf2_forde_aux['response_category'].isin(['yes'])]]

ticks = ['Van Allen et al.(2015)', 'Hugo et al.(2016)', 'Miao et al.(2018)', 'Miao et al.(2018)', 'Forde et al.(2018)']

def set_box_color(bp, color):
    plt.setp(bp['boxes'], color=color)
    plt.setp(bp['whiskers'], color='black')

    plt.setp(bp['caps'], color='black')
    plt.setp(bp['medians'], color='black')

plt.figure()

bpl = plt.boxplot(data_a, positions=np.array(range(len(data_a)))*2.0-0.4, sym='', widths=0.6, patch_artist=True)
bpr = plt.boxplot(data_b, positions=np.array(range(len(data_b)))*2.0+0.4, sym='', widths=0.6, patch_artist=True)
set_box_color(bpl, 'green')
set_box_color(bpr, 'lightgreen')
plt.ylabel('Number of fsINDELs that evade NMD', fontsize=20)
plt.yticks(fontsize=20)
plt.xticks(range(0, len(ticks) * 2, 2), ticks,fontsize=18)
plt.xlim(-2, len(ticks)*2)
plt.tight_layout()
plt.show()


#triggersNMD
#p values
tU_miao, tp_miao = mannwhitneyu(plotdf_miao_aux['triggersNMD'].loc[plotdf_miao_aux['response_category'].isin(['intermediate benefit','clinical benefit'])], plotdf_miao_aux['triggersNMD'].loc[plotdf_miao_aux['response_category'].isin(['no clinical benefit'])], alternative = 'greater')
tU_miao2, tp_miao2 = mannwhitneyu(plotdf_miao2_aux['triggersNMD'].loc[plotdf_miao2_aux['response_category'].isin(['clinical benefit'])], plotdf_miao2_aux['triggersNMD'].loc[plotdf_miao2_aux['response_category'].isin(['stable disease','no clinical benefit'])], alternative = 'greater')
tU_van, tp_van = mannwhitneyu(plotdf_van_aux['triggersNMD'].loc[plotdf_van_aux['response_category'].isin(['long-survival','response'])], plotdf_van_aux['triggersNMD'].loc[plotdf_van_aux['response_category'].isin(['nonresponse'])], alternative = 'greater')
tU_hugo, tp_hugo = mannwhitneyu(plotdf_hugo_aux['triggersNMD'].loc[plotdf_hugo_aux['response_category'].isin(['R'])], plotdf_hugo_aux['triggersNMD'].loc[plotdf_hugo_aux['response_category'].isin(['NR'])], alternative = 'greater')
tU_forde, tp_forde = mannwhitneyu(plotdf_forde_aux['triggersNMD'].loc[plotdf_forde_aux['response_category'].isin(['yes'])], plotdf_forde_aux['triggersNMD'].loc[plotdf_forde_aux['response_category'].isin(['no '])], alternative = 'greater')

#boxplot
data_a = [plotdf_van_aux['triggersNMD'].loc[plotdf_van_aux['response_category'].isin(['nonresponse'])], plotdf_hugo_aux['triggersNMD'].loc[plotdf_hugo_aux['response_category'].isin(['NR'])], plotdf_miao_aux['triggersNMD'].loc[plotdf_miao_aux['response_category'].isin(['no clinical benefit'])], plotdf_miao2_aux['triggersNMD'].loc[plotdf_miao2_aux['response_category'].isin(['stable disease','no clinical benefit'])], plotdf_forde_aux['triggersNMD'].loc[plotdf_forde_aux['response_category'].isin(['no '])]]
data_b = [plotdf_van_aux['triggersNMD'].loc[plotdf_van_aux['response_category'].isin(['long-survival','response'])], plotdf_hugo_aux['triggersNMD'].loc[plotdf_hugo_aux['response_category'].isin(['R'])], plotdf_miao_aux['triggersNMD'].loc[plotdf_miao_aux['response_category'].isin(['intermediate benefit','clinical benefit'])], plotdf_miao2_aux['triggersNMD'].loc[plotdf_miao2_aux['response_category'].isin(['clinical benefit'])], plotdf_forde_aux['triggersNMD'].loc[plotdf_forde_aux['response_category'].isin(['yes'])]]

ticks = ['Van Allen et al.(2015)', 'Hugo et al.(2016)', 'Miao et al.(2018)', 'Miao et al.(2018)', 'Forde et al.(2018)']

def set_box_color(bp, color):
    plt.setp(bp['boxes'], color=color)
    plt.setp(bp['whiskers'], color='black')

    plt.setp(bp['caps'], color='black')
    plt.setp(bp['medians'], color='black')

plt.figure()

bpl = plt.boxplot(data_a, positions=np.array(range(len(data_a)))*2.0-0.4, sym='', widths=0.6, patch_artist=True)
bpr = plt.boxplot(data_b, positions=np.array(range(len(data_b)))*2.0+0.4, sym='', widths=0.6, patch_artist=True)
set_box_color(bpl, 'green')
set_box_color(bpr, 'lightgreen')

plt.ylabel('Number of fsINDELs that trigger NMD', fontsize=20)
plt.yticks(fontsize=20)
plt.xticks(range(0, len(ticks) * 2, 2), ticks,fontsize=18)
plt.xlim(-2, len(ticks)*2)
plt.tight_layout()
plt.show()
