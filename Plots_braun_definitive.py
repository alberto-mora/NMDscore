import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
import numpy as np

clin = pd.read_csv('/home/vant/Braun/Clinical_and_Immune_Data.csv',sep='\t') #Datos clínicos recogidos de Braun et al., 2020
df = pd.read_csv('/home/vant/Braun/new_NMDB_NoUCSC_Braun2.0_HD_remix_5.csv',sep='\t') #dataframe con variantes somáticas con NMD score
rna = pd.read_csv('/home/vant/Braun/RNA_Braun.csv', sep='\t') #Datos de expresión recogidos de Braun et al., 2020


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

d = dict(zip(list(clinaux['MAF_Tumor_ID']),list(clinaux['ORR'])))
df['ORR'] = df['Tumor_Sample_Barcode'].map(d)
d = dict(zip(list(clinaux['MAF_Tumor_ID']),list(clinaux['Arm'])))
df['Arm'] = df['Tumor_Sample_Barcode'].map(d)
d = dict(zip(list(clinaux['MAF_Tumor_ID']),list(clinaux['Cohort'])))
df['Cohort'] = df['Tumor_Sample_Barcode'].map(d)
colors = dict(zip(['EARLY DISCONTINUATION DUE TO TOXICITY', 'PD', 'PR', 'NEVER TREATED', 'SD', 'CRPR', 'CR', 'NE'],['gray', 'green', 'lightgreen', 'gray', 'lightgreen', 'lightgreen', 'lightgreen', 'lightgreen']))
df['ORRColors'] = df['ORR'].map(colors)

# crear dataframes para plotear fsINDELs clasificados

plotdf = df.loc[(df['Variant_Classification'].str.contains('Frame_Shift'))].groupby('Tumor_Sample_Barcode')['triggersNMD'].agg('sum').sort_values(ascending=False).reset_index()
cols = ['Arm','Cohort','ORR','ORRColors']
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
cols = ['Arm','Cohort','ORR','ORRColors']
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
d = dict(zip(list(clinaux['MAF_Tumor_ID']),list(clinaux['ORR'])))
dfpat['ORR'] = dfpat['Tumor_Sample_Barcode'].map(d)
d = dict(zip(list(clinaux['MAF_Tumor_ID']),list(clinaux['Cohort'])))
dfpat['Cohort'] = dfpat['Tumor_Sample_Barcode'].map(d)
d = dict(zip(list(clinaux['MAF_Tumor_ID']),list(clinaux['Arm'])))
dfpat['Arm'] = dfpat['Tumor_Sample_Barcode'].map(d)
colors = dict(zip(['EARLY DISCONTINUATION DUE TO TOXICITY', 'PD', 'PR', 'NEVER TREATED', 'SD', 'CRPR', 'CR', 'NE'],['gray', 'green', 'lightgreen', 'gray', 'lightgreen', 'lightgreen', 'lightgreen', 'lightgreen']))
dfpat['ORRColors'] = dfpat['ORR'].map(colors)
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

# Statistics (p_value)

for COHORT in list(set(list(plotdf_aux['Cohort']))):
    plotaux = plotdf_aux.loc[plotdf_aux['Cohort'] == COHORT]
    U1, p = mannwhitneyu(plotaux['triggersNMD'].loc[plotaux['ORR'].isin(['CRPR', 'PR', 'CR', 'SD'])], plotaux['triggersNMD'].loc[plotaux['ORR'].isin(['PD'])])
    print(COHORT + '  -  ' + str(p))

for COHORT in list(set(list(plotdf2_aux['Cohort']))):
    plotaux2 = plotdf2_aux.loc[plotdf2_aux['Cohort'] == COHORT]
    U1, p = mannwhitneyu(plotaux2['notriggersNMD'].loc[plotaux2['ORR'].isin(['CRPR', 'PR', 'CR', 'SD'])], plotaux2['notriggersNMD'].loc[plotaux2['ORR'].isin(['PD'])])
    print(COHORT + '  -  ' + str(p))

U1, p = mannwhitneyu(plotdf_aux['triggersNMD'].loc[plotdf_aux['ORR'].isin(['CRPR', 'PR', 'CR', 'SD'])], plotdf_aux['triggersNMD'].loc[plotdf_aux['ORR'].isin(['PD'])])
U1, p2 = mannwhitneyu(plotdf2_aux['notriggersNMD'].loc[plotdf2_aux['ORR'].isin(['CRPR', 'PR', 'CR', 'SD'])], plotdf2_aux['notriggersNMD'].loc[plotdf2_aux['ORR'].isin(['PD'])])


### para fs en general
df = df.loc[df['Arm'] == 'NIVOLUMAB']
df = df.loc[df['Variant_Classification'].str.contains('Frame')]
df['count'] = 1
plotdf = df.loc[(df['Variant_Classification'].str.contains('Frame_Shift'))].groupby('Tumor_Sample_Barcode')['count'].agg('sum').sort_values(ascending=False).reset_index()
cols = ['Arm','Cohort','ORR','ORRColors']
for col in cols:
    d = dict(zip(list(df['Tumor_Sample_Barcode']),list(df[col])))
    plotdf[col] = plotdf['Tumor_Sample_Barcode'].map(d)

U1, p = mannwhitneyu(plotdf['count'].loc[plotdf['ORRColors'].isin(['lightgreen'])], plotdf['count'].loc[plotdf['ORRColors'].isin(['green'])], alternative = 'greater')

#boxplot
data = [plotdf['count'].loc[plotdf['ORRColors'].isin(['lightgreen'])], plotdf['count'].loc[plotdf['ORRColors'].isin(['green'])]]
fig = plt.figure(figsize =(40, 28))
ax = fig.add_subplot(111)

# Creating axes instance
bp = ax.boxplot(data, patch_artist = True, widths = (0.4,0.4))

colors = ['lightgreen', 'green']

for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)

# changing color and linewidth of
# whiskers
for whisker in bp['whiskers']:
    whisker.set(color ='black',
                linewidth = 1.5,
                linestyle =":")

# changing color and linewidth of
# caps
for cap in bp['caps']:
    cap.set(color ='#8B008B',
            linewidth = 2)

# changing color and linewidth of
# medians
for median in bp['medians']:
    median.set(color ='black',
               linewidth = 3)

# changing style of fliers
for flier in bp['fliers']:
    flier.set(marker ='D',
              color ='#e7298a',
              alpha = 0.5)

# x-axis labels
ax.set_xticklabels(['Responders', 'Nonresponders'],fontsize=25)
ax.set_ylabel('Number of fsINDELs', fontsize=25)
plt.yticks(fontsize=20)

# show plot
plt.show()

####fsINDELs evade
data = [plotdf2_aux['notriggersNMD'].loc[plotdf2_aux['ORRColors'].isin(['lightgreen'])], plotdf2_aux ['notriggersNMD'].loc[plotdf2_aux['ORRColors'].isin(['green'])]]
fig = plt.figure(figsize =(40, 28))
ax = fig.add_subplot(111)

# Creating axes instance
bp = ax.boxplot(data, patch_artist = True, widths = (0.4,0.4))

colors = ['lightgreen', 'green']

for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)

# changing color and linewidth of
# whiskers
for whisker in bp['whiskers']:
    whisker.set(color ='black',
                linewidth = 1.5,
                linestyle =":")

# changing color and linewidth of
# caps
for cap in bp['caps']:
    cap.set(color ='#8B008B',
            linewidth = 2)

# changing color and linewidth of
# medians
for median in bp['medians']:
    median.set(color ='black',
               linewidth = 3)

# changing style of fliers
for flier in bp['fliers']:
    flier.set(marker ='D',
              color ='#e7298a',
              alpha = 0.5)

# x-axis labels
ax.set_xticklabels(['Responders', 'Nonresponders'],fontsize=25)
ax.set_ylabel('Number of fsINDELs that evade NMD', fontsize=25)
plt.yticks(fontsize=20)

# show plot
plt.show()

####fsINDELs trigger
data = [plotdf_aux['triggersNMD'].loc[plotdf_aux['ORRColors'].isin(['lightgreen'])], plotdf_aux ['triggersNMD'].loc[plotdf_aux['ORRColors'].isin(['green'])]]
fig = plt.figure(figsize =(40, 28))
ax = fig.add_subplot(111)

# Creating axes instance
bp = ax.boxplot(data, patch_artist = True, widths = (0.4,0.4))

colors = ['lightgreen', 'green']

for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)

# changing color and linewidth of
# whiskers
for whisker in bp['whiskers']:
    whisker.set(color ='black',
                linewidth = 1.5,
                linestyle =":")

# changing color and linewidth of
# caps
for cap in bp['caps']:
    cap.set(color ='#8B008B',
            linewidth = 2)

# changing color and linewidth of
# medians
for median in bp['medians']:
    median.set(color ='black',
               linewidth = 4)

# changing style of fliers
for flier in bp['fliers']:
    flier.set(marker ='D',
              color ='#e7298a',
              alpha = 0.5)

# x-axis labels
ax.set_xticklabels(['Responders', 'Nonresponders'],fontsize=25)
ax.set_ylabel('Number of fsINDELs that trigger NMD', fontsize=25)
plt.yticks(fontsize=20)

# show plot
plt.show()
#####CM009
plotdf2_aux2 = plotdf2_aux.loc[plotdf2_aux['Cohort']=='CM-009']
plotdf_aux2 = plotdf_aux.loc[plotdf_aux['Cohort']=='CM-009']
####fsINDELs evade
data = [plotdf2_aux2['notriggersNMD'].loc[plotdf2_aux2['ORRColors'].isin(['lightgreen'])], plotdf2_aux2['notriggersNMD'].loc[plotdf2_aux2['ORRColors'].isin(['green'])]]
fig = plt.figure(figsize =(40, 28))
ax = fig.add_subplot(111)

# Creating axes instance
bp = ax.boxplot(data, patch_artist = True, widths = (0.4,0.4))

colors = ['lightgreen', 'green']

for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)

# changing color and linewidth of
# whiskers
for whisker in bp['whiskers']:
    whisker.set(color ='black',
                linewidth = 1.5,
                linestyle =":")

# changing color and linewidth of
# caps
for cap in bp['caps']:
    cap.set(color ='#8B008B',
            linewidth = 2)

# changing color and linewidth of
# medians
for median in bp['medians']:
    median.set(color ='black',
               linewidth = 3)

# changing style of fliers
for flier in bp['fliers']:
    flier.set(marker ='D',
              color ='#e7298a',
              alpha = 0.5)

# x-axis labels
ax.set_xticklabels(['Responders', 'Nonresponders'],fontsize=25)
ax.set_ylabel('Number of fsINDELs that evade NMD', fontsize=25)
plt.yticks(fontsize=20)

# show plot
plt.show()

####fsINDELs trigger
data = [plotdf_aux2['triggersNMD'].loc[plotdf_aux2['ORRColors'].isin(['lightgreen'])], plotdf_aux2['triggersNMD'].loc[plotdf_aux2['ORRColors'].isin(['green'])]]
fig = plt.figure(figsize =(40, 28))
ax = fig.add_subplot(111)

# Creating axes instance
bp = ax.boxplot(data, patch_artist = True, widths = (0.4,0.4))

colors = ['lightgreen', 'green']

for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)

# changing color and linewidth of
# whiskers
for whisker in bp['whiskers']:
    whisker.set(color ='black',
                linewidth = 1.5,
                linestyle =":")

# changing color and linewidth of
# caps
for cap in bp['caps']:
    cap.set(color ='#8B008B',
            linewidth = 2)

# changing color and linewidth of
# medians
for median in bp['medians']:
    median.set(color ='black',
               linewidth = 4)

# changing style of fliers
for flier in bp['fliers']:
    flier.set(marker ='D',
              color ='#e7298a',
              alpha = 0.5)

# x-axis labels
ax.set_xticklabels(['Responders', 'Nonresponders'],fontsize=25)
ax.set_ylabel('Number of fsINDELs that trigger NMD', fontsize=25)
plt.yticks(fontsize=20)

# show plot
plt.show()
#####CM010
plotdf2_aux3 = plotdf2_aux.loc[plotdf2_aux['Cohort']=='CM-010']
plotdf_aux3 = plotdf_aux.loc[plotdf_aux['Cohort']=='CM-010']
####fsINDELs evade
data = [plotdf2_aux3['notriggersNMD'].loc[plotdf2_aux3['ORRColors'].isin(['lightgreen'])], plotdf2_aux3['notriggersNMD'].loc[plotdf2_aux3['ORRColors'].isin(['green'])]]
fig = plt.figure(figsize =(40, 28))
ax = fig.add_subplot(111)

# Creating axes instance
bp = ax.boxplot(data, patch_artist = True, widths = (0.4,0.4))

colors = ['lightgreen', 'green']

for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)

# changing color and linewidth of
# whiskers
for whisker in bp['whiskers']:
    whisker.set(color ='black',
                linewidth = 1.5,
                linestyle =":")

# changing color and linewidth of
# caps
for cap in bp['caps']:
    cap.set(color ='#8B008B',
            linewidth = 2)

# changing color and linewidth of
# medians
for median in bp['medians']:
    median.set(color ='black',
               linewidth = 3)

# changing style of fliers
for flier in bp['fliers']:
    flier.set(marker ='D',
              color ='#e7298a',
              alpha = 0.5)

# x-axis labels
ax.set_xticklabels(['Responders', 'Nonresponders'],fontsize=25)
ax.set_ylabel('Number of fsINDELs that evade NMD', fontsize=25)
plt.yticks(fontsize=20)

# show plot
plt.show()

####fsINDELs trigger
data = [plotdf_aux3['triggersNMD'].loc[plotdf_aux3['ORRColors'].isin(['lightgreen'])], plotdf_aux3['triggersNMD'].loc[plotdf_aux3['ORRColors'].isin(['green'])]]
fig = plt.figure(figsize =(40, 28))
ax = fig.add_subplot(111)

# Creating axes instance
bp = ax.boxplot(data, patch_artist = True, widths = (0.4,0.4))

colors = ['lightgreen', 'green']

for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)

# changing color and linewidth of
# whiskers
for whisker in bp['whiskers']:
    whisker.set(color ='black',
                linewidth = 1.5,
                linestyle =":")

# changing color and linewidth of
# caps
for cap in bp['caps']:
    cap.set(color ='#8B008B',
            linewidth = 2)

# changing color and linewidth of
# medians
for median in bp['medians']:
    median.set(color ='black',
               linewidth = 4)

# changing style of fliers
for flier in bp['fliers']:
    flier.set(marker ='D',
              color ='#e7298a',
              alpha = 0.5)

# x-axis labels
ax.set_xticklabels(['Responders', 'Nonresponders'],fontsize=25)
ax.set_ylabel('Number of fsINDELs that trigger NMD', fontsize=25)
plt.yticks(fontsize=20)

# show plot
plt.show()
#####CM025
plotdf2_aux4 = plotdf2_aux.loc[plotdf2_aux['Cohort']=='CM-025']
plotdf_aux4 = plotdf_aux.loc[plotdf_aux['Cohort']=='CM-025']
####fsINDELs evade
data = [plotdf2_aux4['notriggersNMD'].loc[plotdf2_aux4['ORRColors'].isin(['lightgreen'])], plotdf2_aux4['notriggersNMD'].loc[plotdf2_aux4['ORRColors'].isin(['green'])]]
fig = plt.figure(figsize =(40, 28))
ax = fig.add_subplot(111)

# Creating axes instance
bp = ax.boxplot(data, patch_artist = True, widths = (0.4,0.4))

colors = ['lightgreen', 'green']

for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)

# changing color and linewidth of
# whiskers
for whisker in bp['whiskers']:
    whisker.set(color ='black',
                linewidth = 1.5,
                linestyle =":")

# changing color and linewidth of
# caps
for cap in bp['caps']:
    cap.set(color ='#8B008B',
            linewidth = 2)

# changing color and linewidth of
# medians
for median in bp['medians']:
    median.set(color ='black',
               linewidth = 3)

# changing style of fliers
for flier in bp['fliers']:
    flier.set(marker ='D',
              color ='#e7298a',
              alpha = 0.5)

# x-axis labels
ax.set_xticklabels(['Responders', 'Nonresponders'],fontsize=25)
ax.set_ylabel('Number of fsINDELs that evade NMD', fontsize=25)
plt.yticks(fontsize=20)

# show plot
plt.show()

####fsINDELs trigger
data = [plotdf_aux4['triggersNMD'].loc[plotdf_aux4['ORRColors'].isin(['lightgreen'])], plotdf_aux4['triggersNMD'].loc[plotdf_aux4['ORRColors'].isin(['green'])]]
fig = plt.figure(figsize =(40, 28))
ax = fig.add_subplot(111)

# Creating axes instance
bp = ax.boxplot(data, patch_artist = True, widths = (0.4,0.4))

colors = ['lightgreen', 'green']

for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)

# changing color and linewidth of
# whiskers
for whisker in bp['whiskers']:
    whisker.set(color ='black',
                linewidth = 1.5,
                linestyle =":")

# changing color and linewidth of
# caps
for cap in bp['caps']:
    cap.set(color ='#8B008B',
            linewidth = 2)

# changing color and linewidth of
# medians
for median in bp['medians']:
    median.set(color ='black',
               linewidth = 4)

# changing style of fliers
for flier in bp['fliers']:
    flier.set(marker ='D',
              color ='#e7298a',
              alpha = 0.5)

# x-axis labels
ax.set_xticklabels(['Responders', 'Nonresponders'],fontsize=25)
ax.set_ylabel('Number of fsINDELs that trigger NMD', fontsize=25)
plt.yticks(fontsize=20)

# show plot
plt.show()

####Bar plot (if it is a better reprsentation)
colores = {'Responders':'lightgreen','Nonresponders':'green'}
for COHORT in list(set(list(plotdf_aux['Cohort']))):
	fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(10,9))
	plotaux =  plotdf_aux.loc[plotdf_aux['Cohort'] == COHORT]
	plotaux = plotaux.loc[plotaux['ORR'].isin(['CRPR', 'PR', 'CR', 'PD', 'SD'])]
	plotaux2 =  plotdf2_aux.loc[plotdf2_aux['Cohort'] == COHORT]
	plotaux2 = plotaux2.loc[plotaux2['ORR'].isin(['CRPR', 'PR', 'CR', 'PD', 'SD'])]
	ax1 = plotaux.plot(ax=axes[1],kind='bar', y='triggersNMD',color=plotaux['ORRColors'])
	ax2 = plotaux2.plot(ax=axes[0],kind='bar', y='notriggersNMD',color=plotaux2['ORRColors'], legend=None)
	fig.suptitle('Braun et al -- ' + COHORT)
	ax1.set_ylabel('Nm.FS which trigger NMD')
	ax2.set_ylabel('Nm.FS which evade NMD')
	ax1.set_yticks([1, 6, 11])
	ax2.set_yticks([1, 6, 11])
	ax1.set_yticklabels([0, 5, 10])
	ax2.set_yticklabels([0, 5, 10])
	ax1.set_xticklabels([])
	ax2.set_xticklabels([])
	# plt.xlabel('Patients')
	# plt.ylabel('Nr.FS which evade NMD')
	labels = list(colores.keys())
	handles = [plt.Rectangle((0,0),1,1, color=colores[label]) for label in labels]
	plt.legend(handles,labels,bbox_to_anchor=(1, 1.5), loc='upper left', prop={'size': 12})
	#fig.text(0.02, 0.5, 'Nm.FS which trigger (down) and evade (up) NMD', va='center', rotation='vertical', size=12)

	#for i, j, k in zip(list(np.arange(0,len(plotaux),1)), list(plotaux['triggersNMD']), list(plotaux['ImmunoPhenotype'].replace('Excluded','').replace('No_IF','').replace('Desert','').replace('Infiltrated','*'))):
	# for i, j, k in zip(list(np.arange(0,len(plotaux),1)), list(plotaux['triggersNMD']), list(plotaux['Deletion_6p21.32'].astype(str).replace('WT','').replace('MUT','*').replace('nan','')) ):
	# 	fig.text(i, float(j) + 0.5, str(k), fontsize = 5)
	# for i, j, k in zip(list(np.arange(0,len(plotaux2),1)), list(plotaux2['notriggersNMD']), list(plotaux2['Deletion_6p21.32'].astype(str).replace('WT','').replace('MUT','*').replace('nan','')) ):
	# 	fig.text(i, float(j) + 0.5, str(k), fontsize = 5)
	plt.subplots_adjust(left=.080, bottom=.215, right=.87, top=0.940, wspace=.2, hspace=.80)
	plt.show()

fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(10,7))
plotaux = plotdf_aux.loc[plotdf_aux['ORR'].isin(['CRPR', 'PR', 'CR', 'PD', 'SD'])]
plotaux2 = plotdf2_aux.loc[plotdf2_aux['ORR'].isin(['CRPR', 'PR', 'CR', 'PD', 'SD'])]
ax1 = plotaux.plot(ax=axes[1],kind='bar', y='triggersNMD',color=plotdf['ORRColors'])
ax2 = plotaux2.plot(ax=axes[0],kind='bar', y='notriggersNMD',color=plotdf2['ORRColors'], legend=None)
fig.suptitle('Braun et al, all cohorts')
ax1.set_ylabel('Nm.FS which trigger NMD')
ax2.set_ylabel('Nm.FS which evade NMD')
ax1.set_yticks([1, 6, 11])
ax2.set_yticks([1, 6, 11])
ax1.set_yticklabels([0, 5, 10])
ax2.set_yticklabels([0, 5, 10])
ax1.set_xticklabels([])
ax2.set_xticklabels([])
labels = list(colores.keys())
handles = [plt.Rectangle((0,0),1,1, color=colores[label]) for label in labels]
plt.legend(handles, labels,bbox_to_anchor=(1, 1.5), loc='upper left', prop={'size': 12})
plt.subplots_adjust(left=.080, bottom=.215, right=.87, top=0.940, wspace=.2, hspace=.80)
plt.show()
#barplot para todos los fsINDELs
ax1 = plotdf.plot(kind='bar', y='count',color=plotdf['ORRColors'], legend=None, width=0.8)
ax1.set_xlabel('Patients', fontsize=15)
ax1.set_ylabel('Number of fsINDELs', fontsize=16)
ax1.set_xticklabels([])
plt.show()

####grupo control
###indel sin clasificar
df = df.loc[df['Arm'] == 'EVEROLIMUS']
df = df.loc[df['Variant_Classification'].str.contains('Frame')]
df['count'] = 1
plotdf = df.loc[(df['Variant_Classification'].str.contains('Frame_Shift'))].groupby('Tumor_Sample_Barcode')['count'].agg('sum').sort_values(ascending=False).reset_index()
cols = ['Arm','Cohort','ORR','ORRColors']
for col in cols:
    d = dict(zip(list(df['Tumor_Sample_Barcode']),list(df[col])))
    plotdf[col] = plotdf['Tumor_Sample_Barcode'].map(d)

U1, p = mannwhitneyu(plotdf['count'].loc[plotdf['ORR'].isin(['CRPR', 'PR', 'CR', 'SD'])], plotdf['count'].loc[plotdf['ORR'].isin(['PD'])])
#####indel clasificados
plotdf = df.loc[(df['Variant_Classification'].str.contains('Frame_Shift'))].groupby('Tumor_Sample_Barcode')['triggersNMD'].agg('sum').sort_values(ascending=False).reset_index()
cols = ['Arm','Cohort','ORR','ORRColors']
for col in cols:
    d = dict(zip(list(df['Tumor_Sample_Barcode']),list(df[col])))
    plotdf[col] = plotdf['Tumor_Sample_Barcode'].map(d)

plotdf = plotdf.loc[plotdf['Arm'] == 'EVEROLIMUS']


plotdf2 = df.loc[(df['Variant_Classification'].str.contains('Frame_Shift'))].groupby('Tumor_Sample_Barcode')['notriggersNMD'].agg('sum').sort_values(ascending=False).reset_index()
cols = ['Arm','Cohort','ORR','ORRColors']
for col in cols:
    d = dict(zip(list(df['Tumor_Sample_Barcode']),list(df[col])))
    plotdf2[col] = plotdf2['Tumor_Sample_Barcode'].map(d)

plotdf2 = plotdf2.loc[plotdf2['Arm'] == 'EVEROLIMUS']

U1, p = mannwhitneyu(plotdf['triggersNMD'].loc[plotdf['ORR'].isin(['CRPR', 'PR', 'CR', 'SD'])], plotdf['triggersNMD'].loc[plotdf['ORR'].isin(['PD'])])
U1, p2 = mannwhitneyu(plotdf2['notriggersNMD'].loc[plotdf2['ORR'].isin(['CRPR', 'PR', 'CR', 'SD'])], plotdf2['notriggersNMD'].loc[plotdf2['ORR'].isin(['PD'])])

####fsINDELs evade
data = [plotdf2['notriggersNMD'].loc[plotdf2['ORR'].isin(['CRPR', 'PR', 'CR', 'SD'])], plotdf2['notriggersNMD'].loc[plotdf2['ORR'].isin(['PD'])]]
fig = plt.figure(figsize =(40, 28))
ax = fig.add_subplot(111)

# Creating axes instance
bp = ax.boxplot(data, patch_artist = True, widths = (0.4,0.4))

colors = ['lightgreen', 'green']

for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)

# changing color and linewidth of
# whiskers
for whisker in bp['whiskers']:
    whisker.set(color ='black',
                linewidth = 1.5,
                linestyle =":")

# changing color and linewidth of
# caps
for cap in bp['caps']:
    cap.set(color ='#8B008B',
            linewidth = 2)

# changing color and linewidth of
# medians
for median in bp['medians']:
    median.set(color ='black',
               linewidth = 3)

# changing style of fliers
for flier in bp['fliers']:
    flier.set(marker ='D',
              color ='#e7298a',
              alpha = 0.5)

# x-axis labels
ax.set_xticklabels(['Responders', 'Nonresponders'],fontsize=25)
ax.set_ylabel('Number of fsINDELs that evade NMD', fontsize=25)
plt.yticks(fontsize=20)

# show plot
plt.show()

####fsINDELs trigger
data = [plotdf['triggersNMD'].loc[plotdf['ORRColors'].isin(['lightgreen'])], plotdf['triggersNMD'].loc[plotdf['ORRColors'].isin(['green'])]]
fig = plt.figure(figsize =(40, 28))
ax = fig.add_subplot(111)

# Creating axes instance
bp = ax.boxplot(data, patch_artist = True, widths = (0.4,0.4))

colors = ['lightgreen', 'green']

for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)

# changing color and linewidth of
# whiskers
for whisker in bp['whiskers']:
    whisker.set(color ='black',
                linewidth = 1.5,
                linestyle =":")

# changing color and linewidth of
# caps
for cap in bp['caps']:
    cap.set(color ='#8B008B',
            linewidth = 2)

# changing color and linewidth of
# medians
for median in bp['medians']:
    median.set(color ='black',
               linewidth = 4)

# changing style of fliers
for flier in bp['fliers']:
    flier.set(marker ='D',
              color ='#e7298a',
              alpha = 0.5)

# x-axis labels
ax.set_xticklabels(['Responders', 'Nonresponders'],fontsize=25)
ax.set_ylabel('Number of fsINDELs that trigger NMD', fontsize=25)
plt.yticks(fontsize=20)

# show plot
plt.show()
