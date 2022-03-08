#From a results file, groups by phages and extracts different stats

import json
import os
from numpy import delete
import pandas as pd
from Bio import SeqIO, motifs
from Bio.SeqUtils import GC

def readJSON(filename):
    f=open(filename)
    print('Reading '+filename)
    json_obj=json.load(f)
    resultsFile=json_obj['Results_filename']
    resultsPath=json_obj['Results_path']
    threshold=json_obj['Threshold']
    return resultsPath, resultsFile, threshold

def getResults(path, filename):
    input_filepath=os.path.join(path, filename)
    csvfile=open(input_filepath,newline='')
    results=pd.read_csv(csvfile,header=0)
    return results

def readGB(filename):
    return SeqIO.read(filename,'genbank')

def getGB_objects(results):
    genomes=results['Genome ID'].unique()
    list_of_gb=[]
    for genome in genomes:
        gnome_filepath = os.path.join('genomes', genome+'.gb')
        gb=readGB(gnome_filepath)
        list_of_gb.append(gb)
    return list_of_gb,genomes

def write_genomes_stats(list_of_gbs, gen_IDs):
    lengths=[]
    gc_content=[]
    num_of_genes=[]
    names=[]
    genIDs=gen_IDs.tolist()
    for gb in list_of_gbs:
        lengths.append(len(gb.seq))
        genes=0
        names.append(gb.annotations['organism'])
        for f in gb.features:
            if f.type=='CDS':
                genes+=1
        num_of_genes.append(genes)
        
        
        gc_content.append(GC(gb.seq))
    
    df=pd.DataFrame(
        {'Genome ID': genIDs,
        'Organism': names,
        'Length': lengths,
        'GC Content': gc_content,
        '# of Genes': num_of_genes
        })
    with open('results/analysis/analysis.csv','w',newline='') as f:
        df.to_csv(f)
    return df

def count(results,genomes):
    start=0
    numOps=[]
    numIntras=[]
    numInters=[]
    for index, row in genomes.iterrows():
        end=start+row['Hits']
        subset=results[start:end]
        numOps.append(len(subset[subset['Region']=='OPERATOR']))
        numIntras.append(len(subset[subset['Region']=='INTRAGENETIC']))
        numInters.append(len(subset[subset['Region']=='INTERGENETIC']))
        start=end
    
    genomes['# of Operators']=numOps
    genomes['# of Intragenic']=numIntras
    genomes['# of Intergenic']=numInters
    return genomes

def find_tandems(results, genomes):
    start=0
    tandems=[]
    for index, row in genomes.iterrows():
        end=start+row['Hits']
        subset=results[start:end]
        abspos=0
        tandem=0
        for index2,row2 in subset.iterrows():
            if abspos!=0:
                if row2['Region']=='OPERATOR':
                    if abs(abspos-row2['Abs. Position'])<=100:
                        if valid==1:
                            tandem+=1
                    valid=1
                else:
                    valid=0
            abspos=row2['Abs. Position']
            
            

        tandems.append(tandem)
        start=end
    genomes['Tandem']=tandems
    return genomes
        
def compute_avgs(results, genomes, threshold):
    start=0
    
    avgSc=[]
    avgScOP=[]
    avgScINTRA=[]
    avgScINTER=[]
    
    avgPos=[]
    avgPosOP=[]
    avgPosINTRA=[]
    avgPosINTER=[]
    
    avgNormSc=[]
    avgNormScOP=[]
    avgNormScINTRA=[]
    avgNormScINTER=[]

    sumNormSc=[]
    sumNormScOP=[]
    sumNormScINTRA=[]
    sumNormScINTER=[]

    for index, row in genomes.iterrows():
        end=start+row['Hits']

        phage=results[start:end]
        phageOP=phage[phage['Region']=='OPERATOR']
        phageINTRA=phage[phage['Region']=='INTRAGENETIC']
        phageINTER=phage[phage['Region']=='INTERGENETIC']

        avgSc.append(phage['Score'].mean())
        avgNormSc.append(phage['Score'].mean()-threshold)
        sumNormSc.append((phage['Score']-threshold).sum())

        if "undefined" not in phage["Position relative to gene"].unique():
            phage["Position relative to gene"] = pd.to_numeric(phage["Position relative to gene"], downcast="float")
            phageOP["Position relative to gene"] = pd.to_numeric(phageOP["Position relative to gene"], downcast="float")
            phageINTRA["Position relative to gene"] = pd.to_numeric(phageINTRA["Position relative to gene"], downcast="float")
            phageINTER["Position relative to gene"] = pd.to_numeric(phageINTER["Position relative to gene"], downcast="float")

            
            avgScOP.append(phageOP['Score'].mean())
            avgScINTRA.append(phageINTRA['Score'].mean())
            avgScINTER.append(phageINTER['Score'].mean())

            avgPos.append(phage['Position relative to gene'].mean())
            avgPosOP.append(phageOP['Position relative to gene'].mean())
            avgPosINTRA.append(phageINTRA['Position relative to gene'].mean())
            avgPosINTER.append(phageINTER['Position relative to gene'].mean())
            
            
            avgNormScOP.append(phageOP['Score'].mean()-threshold)
            avgNormScINTRA.append(phageINTRA['Score'].mean()-threshold)
            avgNormScINTER.append(phageINTER['Score'].mean()-threshold)

            
            sumNormScOP.append((phageOP['Score']-threshold).sum())
            sumNormScINTRA.append((phageINTRA['Score']-threshold).sum())
            sumNormScINTER.append((phageINTER['Score']-threshold).sum())
        else:
            
            avgScOP.append('undefined')
            avgScINTRA.append('undefined')
            avgScINTER.append('undefined')

            avgPos.append('undefined')
            avgPosOP.append('undefined')
            avgPosINTRA.append('undefined')
            avgPosINTER.append('undefined')
            
            
            avgNormScOP.append('undefined')
            avgNormScINTRA.append('undefined')
            avgNormScINTER.append('undefined')

            
            sumNormScOP.append('undefined')
            sumNormScINTRA.append('undefined')
            sumNormScINTER.append('undefined')

        start=end

        
    
    genomes['Avg. Score']=avgSc
    genomes['Avg. Score (OPERATOR)']=avgScOP
    genomes['Avg. Score (INTRAGENETIC)']=avgScINTRA
    genomes['Avg. Score (INTERGENETIC)']=avgScINTER

    genomes['Avg. Normalized Score']=avgNormSc
    genomes['Avg. Normalized Score (OPERATOR)']=avgNormScOP
    genomes['Avg. Normalized Score (INTRAGENETIC)']=avgNormScINTRA
    genomes['Avg. Normalized Score (INTERGENETIC)']=avgNormScINTER

    genomes['Sum of Normalized Scores']=sumNormSc
    genomes['Sum of Normalized Scores (OPERATOR)']=sumNormScOP
    genomes['Sum of Normalized Scores (INTRAGENETIC)']=sumNormScINTRA
    genomes['Sum of Normalized Scores (INTERGENETIC)']=sumNormScINTER

    genomes['Avg. Relative Position']=avgPos
    genomes['Avg. Relative Position (OPERATOR)']=avgPosOP
    genomes['Avg. Relative Position (INTRAGENETIC)']=avgPosINTRA
    genomes['Avg. Relative Position (INTERGENETIC)']=avgPosINTER

    return genomes



json_filename='results.json'
results_path,results_filename, threshold=readJSON(json_filename)
results=getResults(results_path,results_filename)
if os.path.isfile('results/analysis/analysis.csv'):
    csvfile=open('results/analysis/analysis.csv',newline='')
    gbExtra=pd.read_csv(csvfile,header=0)
else:
    gbs,gen_names=getGB_objects(results)
    gbExtra=write_genomes_stats(gbs,gen_names)
dups_results=results['Genome ID'].value_counts(sort=False).to_frame().reset_index()
gbExtra['Hits']=dups_results.iloc[:,1]
gbExtra=count(results, gbExtra)

gbExtra=find_tandems(results, gbExtra)
gbExtra=compute_avgs(results, gbExtra, threshold)
with open('results/analysis/analysis.csv','w',newline='') as f:
        gbExtra.to_csv(f)