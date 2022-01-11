import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
import numpy as np
from scipy import stats

clin = pd.read_csv('/home/vant/Braun/Clinical_and_Immune_Data.csv',sep='\t') #Datos clínicos recogidos de Braun et al., 2020
df = pd.read_csv('/home/vant/Braun/new_NMDB_NoUCSC_Braun2.0_HD_remix_5.csv',sep='\t') #dataframe con variantes somáticas con NMD score
rna = pd.read_csv('/home/vant/Braun/RNA_Braun.csv', sep='\t') #Datos de expresión recogidos de Braun et al., 2020
cs = pd.read_csv('/home/vant/Braun/Braun_CIBERSORTX.csv', sep='\t') #Datos de deconvolution inmune recogidos de Braun et al., 2020

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

df = df.loc[df['NMDB'] != '.']
df['NMDB'] = df['NMDB'].astype(float)
df['triggersNMD'] = df.apply(lambda x: trigger_nmd(x['NMDB']), axis=1)
df['notriggersNMD'] = df.apply(lambda x: no_trigger_nmd(x['NMDB']), axis=1)

clinaux = clin.loc[~clin['MAF_Tumor_ID'].isnull()]
d = dict(zip(list(clinaux['MAF_Tumor_ID']),list(clinaux['RNA_ID'])))
df['RNA_ID'] = df['Tumor_Sample_Barcode'].map(d)
gene_list = list(rna['gene_name'])
df = df[df['Hugo_Symbol'].isin(gene_list)]

d = dict(zip(list(clinaux['MAF_Tumor_ID']),list(clinaux['ImmunoPhenotype'])))
df['InmunoPheno'] = df['Tumor_Sample_Barcode'].map(d)
d = dict(zip(list(clinaux['MAF_Tumor_ID']),list(clinaux['Arm'])))
df['Arm'] = df['Tumor_Sample_Barcode'].map(d)
d = dict(zip(list(clinaux['MAF_Tumor_ID']),list(clinaux['Cohort'])))
df['Cohort'] = df['Tumor_Sample_Barcode'].map(d)

###CIBERSORT
cells = cs['CIBERSORTX_NORM'].tolist()
cs.transpose()
cs = cs.transpose()
cs['sample'] = cs.index
colu = cells
colu.append('RNA_ID')
cs.columns = colu
cs = cs.drop('CIBERSORTX_NORM')

for i in cells:
	d = dict(zip(list(cs['RNA_ID']),list(cs[i])))
	df[i] = df['RNA_ID'].map(d)

plotdf = df.loc[(df['Variant_Classification'].str.contains('Frame_Shift'))].groupby('Tumor_Sample_Barcode')['triggersNMD'].agg('sum').sort_values(ascending=False).reset_index()
cols = ['Arm','Cohort','b_cells_naive',
       'b_cells_memory', 'plasma_cells', 't_cells_cd8', 't_cells_cd4_naive',
       't_cells_cd4_memory_resting', 't_cells_cd4_memory_activated',
       't_cells_follicular_helper', 't_cells_regulatory_tregs',
       't_cells_gamma_delta', 'nk_cells_resting', 'nk_cells_activated',
       'monocytes', 'macrophages_m0', 'macrophages_m1', 'macrophages_m2',
       'dendritic_cells_resting', 'dendritic_cells_activated',
       'mast_cells_resting', 'mast_cells_activated', 'eosinophils',
       'neutrophils', 'leukocyte_total', 'InmunoPheno']
for col in cols:
    d = dict(zip(list(df['Tumor_Sample_Barcode']),list(df[col])))
    plotdf[col] = plotdf['Tumor_Sample_Barcode'].map(d)

plotdf = plotdf.loc[plotdf['Arm'] == 'NIVOLUMAB']
plotdf['Tumor_Sample_Barcode_aux'] = plotdf['Tumor_Sample_Barcode']
plotdf['Tumor_Sample_Barcode'] =plotdf['Tumor_Sample_Barcode'].str.replace('-','_')
plotdf[['aux1', 'aux2', 'aux3','aux4','aux5','aux6','aux7','aux8']] = plotdf['Tumor_Sample_Barcode'].str.split('_', expand=True)
plotdf['ID'] = plotdf['aux3'] + '_' + plotdf['aux4']+ '_' + plotdf['aux5']
plotdf['ID'] = plotdf['ID'].fillna('0')
plotdf['ID'] = plotdf.apply(lambda x: x['aux1'] if x['ID'] == '0' else x['ID'], axis=1)
plotdf.drop(columns=['aux1', 'aux2', 'aux3','aux4','aux5','aux6','aux7','aux8'],inplace=True, axis=1)

plotdf = plotdf.set_index('ID')

plotdf2 = df.loc[(df['Variant_Classification'].str.contains('Frame_Shift'))].groupby('Tumor_Sample_Barcode')['notriggersNMD'].agg('sum').sort_values(ascending=False).reset_index()
cols = ['Arm','Cohort','b_cells_naive',
       'b_cells_memory', 'plasma_cells', 't_cells_cd8', 't_cells_cd4_naive',
       't_cells_cd4_memory_resting', 't_cells_cd4_memory_activated',
       't_cells_follicular_helper', 't_cells_regulatory_tregs',
       't_cells_gamma_delta', 'nk_cells_resting', 'nk_cells_activated',
       'monocytes', 'macrophages_m0', 'macrophages_m1', 'macrophages_m2',
       'dendritic_cells_resting', 'dendritic_cells_activated',
       'mast_cells_resting', 'mast_cells_activated', 'eosinophils',
       'neutrophils', 'leukocyte_total', 'InmunoPheno']
for col in cols:
    d = dict(zip(list(df['Tumor_Sample_Barcode']),list(df[col])))
    plotdf2[col] = plotdf2['Tumor_Sample_Barcode'].map(d)

plotdf2 = plotdf2.loc[plotdf2['Arm'] == 'NIVOLUMAB']
plotdf2['Tumor_Sample_Barcode_aux'] = plotdf2['Tumor_Sample_Barcode']
plotdf2['Tumor_Sample_Barcode'] =plotdf2['Tumor_Sample_Barcode'].str.replace('-','_')
plotdf2[['aux1', 'aux2', 'aux3','aux4','aux5','aux6','aux7','aux8']] = plotdf2['Tumor_Sample_Barcode'].str.split('_', expand=True)
plotdf2['ID'] = plotdf2['aux3'] + '_' + plotdf2['aux4']+ '_' + plotdf2['aux5']
plotdf2['ID'] = plotdf2['ID'].fillna('0')
plotdf2['ID'] = plotdf2.apply(lambda x: x['aux1'] if x['ID'] == '0' else x['ID'], axis=1)
plotdf2.drop(columns=['aux1', 'aux2', 'aux3','aux4','aux5','aux6','aux7','aux8'],inplace=True, axis=1)
plotdf2 = plotdf2.set_index('ID')

pat = list(set(list(df['Tumor_Sample_Barcode'])))
dfpat = pd.DataFrame(pat)
dfpat.rename(columns={0:'Tumor_Sample_Barcode'}, inplace=True)
d = dict(zip(list(clinaux['MAF_Tumor_ID']),list(clinaux['ImmunoPhenotype'])))
dfpat['InmunoPheno'] = dfpat['Tumor_Sample_Barcode'].map(d)
d = dict(zip(list(clinaux['MAF_Tumor_ID']),list(clinaux['Cohort'])))
dfpat['Cohort'] = dfpat['Tumor_Sample_Barcode'].map(d)
d = dict(zip(list(clinaux['MAF_Tumor_ID']),list(clinaux['Arm'])))
dfpat['Arm'] = dfpat['Tumor_Sample_Barcode'].map(d)
dfpat_aux = dfpat.copy()
dfpat_aux_2 = dfpat.copy()
dfpatno = dfpat_aux
dfpatyes = dfpat_aux_2
dfpatno['notriggersNMD'] = 0
dfpatyes['triggersNMD'] = 0
mega_pat = list(set(list(plotdf2['Tumor_Sample_Barcode_aux'])))
mega_pat_2 = list(set(list(plotdf['Tumor_Sample_Barcode_aux'])))
dfpatno = dfpatno.loc[~dfpatno['Tumor_Sample_Barcode'].isin(mega_pat)]
dfpatyes = dfpatyes.loc[~dfpatyes['Tumor_Sample_Barcode'].isin(mega_pat_2)]
plotdf_aux = pd.concat([plotdf, dfpatyes], ignore_index=True)
plotdf_aux = plotdf_aux.loc[plotdf_aux['Arm'] == 'NIVOLUMAB']
plotdf2_aux = pd.concat([plotdf2, dfpatno], ignore_index=True)
plotdf2_aux = plotdf2_aux.loc[plotdf2_aux['Arm'] == 'NIVOLUMAB']

plotdf_aux['Tumor_Sample_Barcode'] =plotdf_aux['Tumor_Sample_Barcode'].str.replace('-','_')
plotdf_aux[['aux1', 'aux2', 'aux3','aux4','aux5','aux6','aux7','aux8']] = plotdf_aux['Tumor_Sample_Barcode'].str.split('_', expand=True)
plotdf_aux['ID'] = plotdf_aux['aux3'] + '_' + plotdf_aux['aux4']+ '_' + plotdf_aux['aux5']
plotdf_aux['ID'] = plotdf_aux['ID'].fillna('0')
plotdf_aux['ID'] = plotdf_aux.apply(lambda x: x['aux1'] if x['ID'] == '0' else x['ID'], axis=1)
plotdf_aux.drop(columns=['aux1', 'aux2', 'aux3','aux4','aux5','aux6','aux7','aux8'],inplace=True, axis=1)
plotdf_aux = plotdf_aux.set_index('ID')

plotdf2_aux['Tumor_Sample_Barcode'] =plotdf2_aux['Tumor_Sample_Barcode'].str.replace('-','_')
plotdf2_aux[['aux1', 'aux2', 'aux3','aux4','aux5','aux6','aux7','aux8']] = plotdf2_aux['Tumor_Sample_Barcode'].str.split('_', expand=True)
plotdf2_aux['ID'] = plotdf2_aux['aux3'] + '_' + plotdf2_aux['aux4']+ '_' + plotdf2_aux['aux5']
plotdf2_aux['ID'] = plotdf2_aux['ID'].fillna('0')
plotdf2_aux['ID'] = plotdf2_aux.apply(lambda x: x['aux1'] if x['ID'] == '0' else x['ID'], axis=1)
plotdf2_aux.drop(columns=['aux1', 'aux2', 'aux3','aux4','aux5','aux6','aux7','aux8'],inplace=True, axis=1)
plotdf2_aux = plotdf2_aux.set_index('ID')

desert = plotdf2_aux['notriggersNMD'].loc[plotdf2_aux['InmunoPheno'] == 'Desert']
excluded =  plotdf2_aux['notriggersNMD'].loc[plotdf2_aux['InmunoPheno'] == 'Excluded']
infiltrated =  plotdf2_aux['notriggersNMD'].loc[plotdf2_aux['InmunoPheno'] == 'Infiltrated']
desert2 = plotdf_aux['triggersNMD'].loc[plotdf_aux['InmunoPheno'] == 'Desert']
excluded2 =  plotdf_aux['triggersNMD'].loc[plotdf_aux['InmunoPheno'] == 'Excluded']
infiltrated2 =  plotdf_aux['triggersNMD'].loc[plotdf_aux['InmunoPheno'] == 'Infiltrated']

#p value en distinta categoría
U1, p = mannwhitneyu(plotdf2_aux['notriggersNMD'].loc[plotdf2_aux['InmunoPheno'] == 'Desert'], plotdf2_aux['notriggersNMD'].loc[plotdf2_aux['InmunoPheno'] == 'Infiltrated'])
U3, p3 = mannwhitneyu(plotdf_aux['triggersNMD'].loc[plotdf_aux['InmunoPheno'] == 'Desert'], plotdf_aux['triggersNMD'].loc[plotdf_aux['InmunoPheno'] == 'Infiltrated'])

#p value en la misma categoría
U5, p5 = mannwhitneyu(plotdf2_aux['notriggersNMD'].loc[plotdf2_aux['InmunoPheno'] == 'Desert'], plotdf_aux['triggersNMD'].loc[plotdf_aux['InmunoPheno'] == 'Desert'])
U6, p6 = mannwhitneyu(plotdf2_aux['notriggersNMD'].loc[plotdf2_aux['InmunoPheno'] == 'Infiltrated'], plotdf_aux['triggersNMD'].loc[plotdf_aux['InmunoPheno'] == 'Infiltrated'])

#boxplot
pos = [1,2,4,5]
fig, ax = plt.subplots()
parts = ax.violinplot([desert, desert2, infiltrated, infiltrated2], pos)
ticks = ['Desert','Infiltrated']
plt.xticks([1.5,4.5], ticks,fontsize=20)
ax.set_ylabel('Nm.FS which evade/trigger NMD',fontsize=25)
plt.yticks(fontsize=20)

parts['bodies'][0].set_facecolor('blue')
parts['bodies'][0].set_edgecolor('black')
parts['bodies'][1].set_facecolor('red')
parts['bodies'][1].set_edgecolor('black')
parts['bodies'][2].set_facecolor('blue')
parts['bodies'][2].set_edgecolor('black')
parts['bodies'][3].set_facecolor('red')
parts['bodies'][3].set_edgecolor('black')


plt.show()

###CIBERSORT
cells2 = cells[0:-1]
plotdf2_aux = plotdf2_aux.dropna()
#plotear todo
for COHORT in list(set(list(plotdf2_aux['Cohort']))):
	plotaux =  plotdf2_aux.loc[plotdf2_aux['Cohort'] == COHORT]
	for i in cells2:
		plotaux.sort_values(by=i,ascending=True).plot.scatter(i,'notriggersNMD')
		plt.text(0.01, 4, stats.pearsonr(plotaux['notriggersNMD'], plotaux[i]))
		plt.show()
#para saber cuales son significativos
for COHORT in list(set(list(plotdf2_aux['Cohort']))):
	plotaux =  plotdf2_aux.loc[plotdf2_aux['Cohort'] == COHORT]
	for i in cells2:
		print(str(COHORT) + ' - '+ i + ' - ' + str(stats.pearsonr(plotaux['notriggersNMD'], plotaux[i])))

for i in cells:
	print(i + ' - ' + str(stats.pearsonr(plotdf2_aux['notriggersNMD'], plotdf2_aux[i])))

#plotear los que han salido significativos
plotdf2_aux.sort_values(by='mast_cells_resting',ascending=True).plot.scatter('mast_cells_resting','notriggersNMD')
plt.xlabel('Abundance of resting mast cells', fontsize=20)
plt.ylabel('Number of fsINDELs that evade NMD', fontsize=20)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.show()

stats.pearsonr(plotdf2_aux['mast_cells_resting'], plotdf2_aux['notriggersNMD'])

plotaux =  plotdf2_aux.loc[plotdf2_aux['Cohort'] == 'CM-010']
plotaux2 =  plotdf2_aux.loc[plotdf2_aux['Cohort'] == 'CM-025']

plotaux.sort_values(by='t_cells_cd4_memory_activated',ascending=True).plot.scatter('t_cells_cd4_memory_activated','notriggersNMD')
plt.xlabel('Abundance of CD4 activated memory T cells', fontsize=20)
plt.ylabel('Number of fsINDELs that evade NMD', fontsize=20)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.show()

stats.pearsonr(plotaux['t_cells_cd4_memory_activated'], plotaux['notriggersNMD'])

plotaux2.sort_values(by='t_cells_cd4_memory_resting',ascending=True).plot.scatter('t_cells_cd4_memory_resting','notriggersNMD')
plt.xlabel('Abundance of CD4 resting memory T cells', fontsize=20)
plt.ylabel('Number of fsINDELs that evade NMD', fontsize=20)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.show()

stats.pearsonr(plotaux2['t_cells_cd4_memory_resting'], plotaux2['notriggersNMD'])

plotaux2.sort_values(by='mast_cells_resting',ascending=True).plot.scatter('mast_cells_resting','notriggersNMD')
plt.xlabel('Abundance of resting mast cells', fontsize=20)
plt.ylabel('Number of fsINDELs that evade NMD', fontsize=20)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.show()

stats.pearsonr(plotaux2['mast_cells_resting'], plotaux2['notriggersNMD'])
