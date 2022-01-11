import sys
import re
import argparse
import pandas as pd
import numpy as np
import numbers
from gtfparse import read_gtf
import subprocess
import gzip
from Bio import SeqIO
from Bio.Seq import Seq

##Abre el archivo csv
filename = '/home/vant/Miao/Miao_2.csv' #Dataframe con las variantes somáticas, se puede descargar de Miao et al., 2018

## Carga archivos de referencia  Archivo GTF se puede descargar de https://figshare.com/articles/dataset/NMDetective/7803398

gtf = read_gtf("/home/vant/NMDetective/CDS_hg19_NMDetectiveB_Lindeboom_et_al.v2.gtf") #GTF que contiene las zonas CDS de los trasncritos
gtf_all = read_gtf("/home/vant/NMDetective/hg19_NMDetectiveB_Lindeboom_et_al.v2.gtf") #GTF que contiene toda la información del NMDetectiveB
gtf_all = gtf_all.loc[(gtf_all['feature'] == 'exon') | (gtf_all['feature'] == 'CDS')] #Ahora contiene zonas CDS y exones enteros de los trasncritos, para sacar la zona utr
gtfA = read_gtf("/home/vant/NMDetective/CDS_hg19_NMDetectiveA_Lindeboom_et_al.v2.gtf")

gtf['Feature'] = gtf['exon_id'].str.split('.').str[0]
# Genoma de referencia
refgen = '/home/vant/Braun/ucsc.hg19.fasta'

def nmdxtract(ID, PTCposition, mega_gtf, consequence):
	##Sacamos los cromosomas y las posiciones
	CHR = ID.split('_')[0]
	if 'Nonsense' in consequence:
		POS = int(ID.split('_')[1])
		mega_gtf = mega_gtf.loc[(mega_gtf['seqname'] == CHR) & (mega_gtf['start'] <= POS) & (mega_gtf['end'] >= POS)]
		if mega_gtf.empty == False:
			start = list(mega_gtf['start'].values)[0]
			end = list(mega_gtf['end'].values)[0]
			nmd = mega_gtf['nmd_score'].head(1).values[0].split(',')
			strand = str(mega_gtf['strand'])
			r = np.arange(start, end + 1, 1)
			r = r.tolist()
			d = dict(zip(r, nmd))
			return d[POS]
		else:
			return '.'
	else:
		if (PTCposition != '.') & (PTCposition != 'utr') & (PTCposition != 'out'):
			POS = int(PTCposition)
			mega_gtf = mega_gtf.loc[(mega_gtf['seqname'] == CHR) & (mega_gtf['start'] <= POS) & (mega_gtf['end'] >= POS)]
			if mega_gtf.empty == False:
        	##Extraemos del GTF las filas que coinciden (tendría que salir 1, mirar esto)
        	#mega_gtf.loc[(mega_gtf['seqname'] == CHR) & (mega_gtf['start'] <= POS) & (mega_gtf['end'] >= POS)]
        	##Buscar el primer número del intervalo
				start = list(mega_gtf['start'].values)[0]
				##Buscar el segundo número del intervalo
				end = list(mega_gtf['end'].values)[0]
				##Extrae los valores de NMD del intervalo y lo transforma a lista separando por comas
				nmd = mega_gtf['nmd_score'].head(1).values[0].split(',')
				# Extrae la strand para ver si es forward o reverse
				strand = str(mega_gtf['strand'])
				##Crea el intervalo, se pone +1 para que coincida con len(nmd)
				r = np.arange(start, end + 1, 1)
				r = r.tolist()
				##Crea un diccionario que une posición en el intervalo con score de NMD
				d = dict(zip(r, nmd))
				return d[POS]
			else:
				return '.'
		elif PTCposition  == '.':
			return '.'
		elif PTCposition == 'utr':
			return '0.00'
		else:
			return '0.666'

def getPTCpositionFromFS(gtf, ID, refgen, consequence, ID_specific):
	if 'Nonsense' in consequence:
		return '.'
	else:
		ID = ID.split('_')
		CHR = ID[0]
		POS = int(ID[1])
		REF = ID[2]
		ALT = ID[3]
		ID_aux = '_'.join([CHR,str(POS)])

		if gtf.loc[(gtf['seqname'] == CHR) & (gtf['start'] <= POS) & (gtf['end'] >= POS)].empty == False:
			feature = gtf['Feature'].loc[(gtf['seqname'] == CHR) & (gtf['start'] <= POS) & (gtf['end'] >= POS)].values[0] #Coges el primer tránscrito aunque en esa posición genómica haya varios tránscritos
			exons = gtf.loc[gtf['exon_id'].str.contains(feature)]
			strand = exons['strand'].values[0]
			ucsc_id = exons['gene_id'].values[0]
			exons['exon_number'] = exons['exon_number'].astype(int)
			gtf_aux = gtf_all.loc[gtf_all['gene_id'] == ucsc_id]
			gtf_aux['exon_number'] = gtf_aux['exon_number'].astype(int)
			if strand == '+': #Para añadir zona UTR
				if gtf_aux['exon_number'].values[-1] == exons['exon_number'].values[-1]: #mirar porque igual no tiene sentido
					row_utr = gtf_aux.loc[(gtf_aux['feature'] == 'exon') & (gtf_aux['start'] == (exons['start'].values[-1]))]
					row_utr['start'] = (exons['end'].values[-1]) + 1
					row_utr['exon_number'] = row_utr['exon_number'].astype(int)
					row_utr['exon_number'] = row_utr['exon_number'] + 1
					exons = pd.concat([exons, row_utr], ignore_index = True, axis = 0)
				else:
					row_utr = gtf_aux.loc[(gtf_aux['feature'] == 'exon') & (gtf_aux['exon_number'] > (exons['exon_number'].values[-1]))]
					exons = pd.concat([exons, row_utr], ignore_index = True, axis = 0)
			else:
				if gtf_aux['exon_number'].values[0] == exons['exon_number'].values[0]:
					row_utr = gtf_aux.loc[(gtf_aux['feature'] == 'exon') & (gtf_aux['end'] == (exons['end'].values[0]))]
					row_utr['end'] = (exons['start'].values[0]) - 1
					row_utr['exon_number'] = row_utr['exon_number'].astype(int)
					row_utr['exon_number'] = row_utr['exon_number'] - 1
					exons = pd.concat([row_utr, exons], ignore_index = True, axis = 0)
				else:
					row_utr = gtf_aux.loc[(gtf_aux['feature'] == 'exon') & (gtf_aux['exon_number'] < (exons['exon_number'].values[0]))]
					exons = pd.concat([row_utr, exons], ignore_index = True, axis = 0)
			#if strand == '-': #probar quitandolo para ver si funciona
				#exons['exon_number'] = list(exons['exon_number'])[::-1]
			mutExon = exons.loc[(exons['seqname'] == CHR) & (exons['start'] <= POS) & (exons['end'] >= POS)]
			exons['ID'] = exons['seqname'] + ':' + exons['start'].astype(str) + '-' + exons['end'].astype(str)
			filein = ID_aux + '.exons.csv'
			exons['ID'].to_csv(filein, sep='\t', index = None, header = None)
			#Get nucleotide sequences from all exons
			fileout = ID_aux + '.outExons.fa'
			command = ' '.join(['samtools faidx', refgen, '--region-file', filein, '>', fileout])
			subprocess.run(command, shell=True)
			exonsSeqsdict = {rec.id : str(rec.seq) for rec in SeqIO.parse(fileout, "fasta")}
			exons['sequence'] = exons['ID'].map(exonsSeqsdict)
			wild_exons = exons
			wild_exons['range'] = wild_exons.apply(lambda x: list(range(x['start'], (x['end'] + 1))), axis = 1)
			wild_exons['range'] = [','.join(map(str, l)) for l in wild_exons['range']]
			wild_sequence = ''.join(list(wild_exons['sequence']))
			wild_numbers = ','.join(list(wild_exons['range']))
			strand = wild_exons['strand'].values[0]
			base = list(wild_sequence)
			coordinate = wild_numbers.split(',')
			wild_df = pd.DataFrame({'base': base, 'coordinate': coordinate})
			if 'Del' in consequence:
				base_mut = list(REF[1:])
				index = wild_df.index
				condition = wild_df['coordinate'] == str(POS)
				wild_df_indices = index[condition]
				wild_df_indices = wild_df_indices.tolist()
				dele = list(range((wild_df_indices[0] + 1), (wild_df_indices[0] + len(base_mut) + 1)))
				mut_df = wild_df.drop(dele)
				offset = len(''.join(list(mut_df['base']))) % 3
				if strand == '+':
					if offset == 0:
						sequence = Seq(''.join(list(mut_df['base'])))
					else:
						sequence = Seq(''.join(list(mut_df['base']))[:-offset])
					sequence_PTC = str(sequence.translate(to_stop=True))
					PTCposinTrans = (len(sequence_PTC)*3)
					if PTCposinTrans >= (len(sequence) - offset):
						return 'out'
					else:
						PTCposition = mut_df['coordinate'].iloc[PTCposinTrans]
						PTCposition = int(PTCposition)
						check = wild_exons.loc[(wild_exons['start'] <= PTCposition) & (wild_exons['end'] >= PTCposition)]
						if check['feature'].values[0] == 'exon':
							return 'utr'
						else:
							return PTCposition
				else:
					if offset == 0:
						sequence = Seq(''.join(list(mut_df['base']))).reverse_complement()
					else:
						sequence = Seq(''.join(list(mut_df['base']))[offset:]).reverse_complement()
					sequence_PTC = str(sequence.translate(to_stop=True))
					PTCposinTrans = (len(sequence_PTC)*3) + 1
					if PTCposinTrans >= (len(sequence) - offset):
						return 'out'
					else:
						last = mut_df.tail(1).index.item()
						PTCposinTrans2 = (last + 1) -  PTCposinTrans
						PTCposition = mut_df['coordinate'].iloc[PTCposinTrans2 - (len(REF)-1)]
						PTCposition = int(PTCposition)
						check = wild_exons.loc[(wild_exons['start'] <= PTCposition) & (wild_exons['end'] >= PTCposition)]
						if check['feature'].values[0] == 'exon':
							return 'utr'
						else:
							return PTCposition
			else:
				base_mut = list(ALT[1:])
				coordinate_mut = ['x'] * len(base_mut)
				variant_df = pd.DataFrame({'base': base_mut, 'coordinate': coordinate_mut})
				index = wild_df.index
				condition = wild_df['coordinate'] == str(POS)
				wild_df_indices = index[condition]
				wild_df_indices = wild_df_indices.tolist()
				wild_df1 = wild_df.iloc[:(wild_df_indices[0] + 1)]
				wild_df2 = wild_df.iloc[(wild_df_indices[0] + 1):]
				mut_df = pd.concat([wild_df1, variant_df, wild_df2], ignore_index = True, axis = 0)
				offset = len(''.join(list(mut_df['base']))) % 3
				if strand == '+':
					if offset == 0:
						sequence = Seq(''.join(list(mut_df['base'])))
					else:
						sequence = Seq(''.join(list(mut_df['base']))[:-offset])
					sequence_PTC = str(sequence.translate(to_stop=True))
					PTCposinTrans = (len(sequence_PTC)*3)
					if PTCposinTrans >= (len(sequence) - offset):
						return 'out'
					else:
						PTCposition = mut_df['coordinate'].iloc[PTCposinTrans]
						if PTCposition == 'x':
							PTCposition = wild_df1['coordinate'].values[-1]
							return PTCposition
						else:
							PTCposition = int(PTCposition)
							check = wild_exons.loc[(wild_exons['start'] <= PTCposition) & (wild_exons['end'] >= PTCposition)]
							if check['feature'].values[0] == 'exon':
								return 'utr'
							else:
								return PTCposition
				else:
					if offset == 0:
						sequence = Seq(''.join(list(mut_df['base']))).reverse_complement()
					else:
						sequence = Seq(''.join(list(mut_df['base']))[offset:]).reverse_complement()
					sequence_PTC = str(sequence.translate(to_stop=True))
					PTCposinTrans = (len(sequence_PTC)*3) + 1
					if PTCposinTrans >= (len(sequence) - offset):
						return 'out'
					else:
						last = mut_df.tail(1).index.item()
						PTCposinTrans2 = (last + 1) -  PTCposinTrans
						PTCposition = mut_df['coordinate'].iloc[PTCposinTrans2]
						if PTCposition == 'x':
							PTCposition = wild_df2['coordinate'].values[0]
							return PTCposition
						else:
							PTCposition = int(PTCposition)
							check = wild_exons.loc[(wild_exons['start'] <= PTCposition) & (wild_exons['end'] >= PTCposition)]
							if check['feature'].values[0] == 'exon':
								return 'utr'
							else:
								return PTCposition
		else:
			return '.'

# Read Braun's mutation file
# OUTPUT list stored into pandas df


df = pd.read_csv(filename, sep='\t')
REF = 'Tumor_Seq_Allele1'
ALT = 'Tumor_Seq_Allele2'
chr = 'Chromosome'
position = 'Start_position'

# Filtras para quedarte con las PTC en los trasncript canonicos del vcf
df = df.loc[df['Variant_Classification'].str.contains('Nonsense|Frame_Shift')]


df['VEP_allele'] = df[REF] + '/' + df[ALT]
df['ID_specific'] = df[chr].astype(str) + '_' + df[position].astype(str) + '_' + df['VEP_allele'] + '_' + df['Tumor_Sample_Barcode']
df[chr] = 'chr' + df[chr].astype(str)

########## MODIFICATION OF INDEL VARIANTS IN AGREEMENT WITH VCF_init
# Insertions, adding the reference
def addREF2Insertions(CHR, POS, refgen):
    region = CHR + ':' + str(POS) + '-' + str(POS)
    fileout = CHR + '_' + str(POS) + '_' + str(POS) + '.fa'
    command = ' '.join(['samtools faidx', refgen, region, '>', fileout])
    subprocess.run(command, shell=True)
    Seqsdict = {rec.id : str(rec.seq) for rec in SeqIO.parse(fileout, "fasta")}
    command = ' '.join(['rm', fileout])
    subprocess.run(command, shell=True)
    return Seqsdict[region]

df['aux'] = df[REF]
df[REF] = df.apply(lambda x: addREF2Insertions(x[chr], x[position], refgen).upper() if x['aux'] == '-' else x[REF], axis = 1)
df[ALT] = df.apply(lambda x:  x[REF] + x[ALT] if x['aux'] == '-' else x[ALT], axis=1)

# Deletions, adding the previous position as reference
def addREF2deletions(CHR, POS, refgen):
    region = CHR + ':' + str(POS) + '-' + str(POS)
    fileout = CHR + '_' + str(POS) + '_' + str(POS) + '.fa'
    command = ' '.join(['samtools faidx', refgen, region, '>', fileout])
    subprocess.run(command, shell=True)
    Seqsdict = {rec.id : str(rec.seq) for rec in SeqIO.parse(fileout, "fasta")}
    command = ' '.join(['rm', fileout])
    subprocess.run(command, shell=True)
    return Seqsdict[region]

df['aux'] = df[ALT]
# Position -1 for deletions
df[position] = df.apply(lambda x: x[position] - 1 if x['aux'] == '-' else x[position], axis=1)
df[ALT] = df.apply(lambda x: addREF2deletions(x[chr], x[position], refgen).upper() if x['aux'] == '-' else x[ALT], axis = 1)
df[REF] = df.apply(lambda x: x[ALT] + x[REF] if x['aux'] == '-' else x[REF], axis = 1)

df.drop(columns='aux',inplace = True)
####################################################################
####################################################################
df['ID'] = df[chr].astype(str) + '_' + df[position].astype(str) + '_' + df[REF] + '_' + df[ALT]

#VCF_init = '/mnt/64716603-5b56-4f9a-b195-c11560647a3a/data/TFM_NMD/Braun/sorted.chr.Annot_input_VEP_Braun.vcf.gz'
#VCF_init = '/home/vant/Braun/Braun.vcf.gz'
df['PTCposition'] = df.apply(lambda x: getPTCpositionFromFS(gtf, x['ID'], refgen, x['Variant_Classification'], x['ID_specific']), axis=1)

df = df.reset_index()
df.drop('index',1, inplace = True)


df['NMDB'] = df.apply(lambda x: nmdxtract(x['ID'], x['PTCposition'], gtf, x['Variant_Classification']), axis=1)
df['NMDA'] = df.apply(lambda x: nmdxtract(x['ID'], x['PTCposition'], gtfA, x['Variant_Classification']), axis=1)

# Save NMD-annotated df (NMD_Braun.csv)

df.to_csv('new_NMD_Miao2.0.csv',sep='\t',index = None)
