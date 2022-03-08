from contextlib import nullcontext
import os
from os import listdir
from ete3 import NCBITaxa
from math import log2
import json
from Bio import SeqIO
from Bio import motifs
from Bio import Entrez
import csv
from random import randint
ncbi=NCBITaxa()
#UNCOMMENT if desired to update taxonomy database (time consuming)
#ncbi.update_taxonomy_database()
Entrez.email = 'user@example.com'

class HitResult:
    #Object that contains all the information of interest about a hit in its atributes

    def __init__(self, gen_ID,organism, tf_ID, location, score, dist_to_gene, type, gene):
        self.tf_ID=tf_ID
        self.organism=organism
        self.genome_ID=gen_ID
        self.score=score
        self.location=location
        if gene:
            if 'locus_tag' in gene.qualifiers.keys():
                self.gene_locus_tag=gene.qualifiers['locus_tag']
            elif 'note' in gene.qualifiers.keys():
                self.gene_locus_tag=gene.qualifiers['note']
            else:
                self.gene_locus_tag='undefined'
            
            if 'product' in gene.qualifiers.keys():
                self.gene_product=gene.qualifiers['product']
            else:
                self.gene_product='undefined'

            self.gene_start=gene.location.start
            self.gene_stop=gene.location.end
            self.gene_strand=gene.strand
            self.dist_to_gene=dist_to_gene
            self.type=type
        else:
            self.gene_product='undefined'
            self.gene_locus_tag='undefined'
            self.gene_start='undefined'
            self.gene_stop='undefined'
            self.gene_strand='undefined'
            self.dist_to_gene='undefined'
            self.type='undefined'

        


def readJSON(filename):
    #Reads the configuration JSON file that specifies all the parameters necessary

    f=open(filename)
    print('Reading '+filename)
    json_obj=json.load(f)
    out_path=json_obj['Output_path']
    gnome_path=json_obj['genomes_path']
    motifs_path=json_obj['motifs_path']
    motifs_file=json_obj['TF_filename']
    gen_file=json_obj['Genome_filename']
    clade=json_obj['clade']
    up_dist=json_obj['upstream_dist']
    down_dist=json_obj['downstream_dist']
    pseudo=json_obj['pseudocounts']
    return out_path,motifs_path, gnome_path, motifs_file,gen_file, clade, up_dist, down_dist,pseudo

def readTF(filename, pseudo):
    #Reads a TF Motif file and obtains threshold (patser), pssm and reverse pssm

    motifs_parsed = SeqIO.parse(filename,"fasta")
    list_seq=[]
    for record in motifs_parsed:
        list_seq.append(record.seq.upper())
    motif=motifs.create(list_seq)
    
    pwm=motif.counts.normalize(pseudocounts=pseudo)
    pssm=pwm.log_odds()

    #UNCOMMENT following line for random permutation of pssm
    #pssm=permute(pssm,1000)
    
    reverse_pssm=pssm.reverse_complement()
    distribution = pssm.distribution(background=None,precision=10**4)
    threshold=distribution.threshold_patser()

    return motif, threshold, pssm, reverse_pssm

def readGB(filename):
    #Reads genbank file into an object and extracts some information (eg. genetic sequence)

    gb_object = SeqIO.read(filename,'genbank')
    seq=gb_object.seq
    organism = gb_object.annotations['organism']
    gen_ID=gb_object.name

    return gb_object,seq, organism, gen_ID

def permute(pssm, times):
    #Swaps columns of pssm the number of times specified

    for it in range(times):
        col1=randint(0,pssm.length-1)
        col2=randint(0,pssm.length-1)
        abc=['A','C','G','T']
        for letter in abc:
            val1=pssm[letter][col1]
            pssm[letter][col1]=pssm[letter][col2]
            pssm[letter][col2]=val1
    return pssm

def passThreshold(seq,motif,pssm,reverse_pssm,threshold):
    #Calculates scores of each position and saves those positions that pass the threshold

    scores={}
    hits=[]
    i=0
    sc=pssm.calculate(seq)
    rc_scores=reverse_pssm.calculate(seq)
    sc = [log2(2**score + 2**rc_score) for score, rc_score in zip(sc, rc_scores)]
    for score in sc:
        if score>=threshold:
            scores[i]=score
            hits.append(i)
        i+=1
    
    return hits,scores

def check_overlap(hits,motif,scores):
    #Checks if there are hits separated by less than half of the length of the motif, in which case the hit with a lower score is dismissed

    if scores:
        to_delete=[]
        for i in range(len(hits)-1):
            if (hits[i]+len(motif)/2)>=hits[i+1]:
                delete1=hits[i]
                delete2=hits[i+1]
                if scores[delete1]>=scores[delete2]:
                    to_delete.append(i+1)
                else:
                    to_delete.append(i)
        for i in sorted(to_delete,reverse=True):
            del scores[hits[i]]
            del hits[i]



    
    return hits, scores

def dist_classify(hits, seq, gb_object, up_dist, down_dist):
    dist_to_genes={}
    gene_locations={}
    genes={}
    classification={}

    for pos in hits:
        dist_final=len(seq)

        #-----Calculate distance to closest gene----------
        for f in gb_object.features:
            if f.type=='CDS':
                if f.strand is 1:
                    dist1=pos-f.location.start
                elif f.strand is -1:
                    dist1=pos-f.location.end
                dist_aux = abs(dist1)
                if dist_aux<dist_final:
                    dist_final=dist_aux
                    dist_to_genes[pos]=dist1
                    gene_locations[pos]=f.location
                    genes[pos]=f

    #------Classification of hit by their relative location to closest gene (if it exists)
        if genes:
            if genes[pos].strand is 1:
                if dist_to_genes[pos]>=-(up_dist) and dist_to_genes[pos]<=down_dist:
                    classification[pos]="OPERATOR"
                elif dist_to_genes[pos]>down_dist and dist_to_genes[pos]<gene_locations[pos].end:
                    classification[pos]="INTRAGENETIC"
                elif dist_to_genes[pos]<-(up_dist) or dist_to_genes[pos]>gene_locations[pos].end:
                    classification[pos]="INTERGENETIC"
            elif genes[pos].strand is -1:
                if dist_to_genes[pos] <= up_dist and dist_to_genes[pos] >= -(down_dist):
                    classification[pos]="OPERATOR"
                elif dist_to_genes[pos] < -(down_dist) and dist_to_genes[pos] > gene_locations[pos].start:
                    classification[pos]="INTRAGENETIC"
                elif dist_to_genes[pos] > up_dist or dist_to_genes[pos] < gene_locations[pos].start:
                    classification[pos]="INTERGENETIC"
    
    return dist_to_genes,genes,classification

def save_to_csv(results, path, threshold):
    #Writes results to a csv file in the prefered format
    #Creates json file for the analysis script that specifies the path of the csv file ant the threshold (useful for analysis script)

    header=[
        "Motif",
        "Genome",
        "Genome ID",
        "Score",
        "Abs. Position",
        "Gene locus_tag",
        "Gene Product",
        "Gene start",
        "Gene stop",
        "Gene Strand",
        "Position relative to gene",
        "Region"
        ]
    filename='results_hits.csv'
    
    filepath=os.path.join(path, filename)
    
    json_filepath=os.path.join(path, "results.json")
    data={
        "Results_filename": filename,
        "Results_path": path,
        "Threshold": threshold
    }
    with open(json_filepath, 'w') as f:
        json.dump(data, f)
    
    with open(filepath,'w',newline='') as f:
        f.truncate(0)
        writer=csv.DictWriter(f,fieldnames=header, delimiter=',')
        writer.writeheader()
        for motif_dict in results.items():
            for genome_results in motif_dict[1].items():
                for hit in genome_results[1]:
                    row={
                        'Motif':hit.tf_ID,
                        'Genome': hit.organism,
                        'Genome ID':hit.genome_ID,
                        'Score':hit.score,
                        'Abs. Position':hit.location,
                        'Gene locus_tag':hit.gene_locus_tag,
                        'Gene Product':hit.gene_product,
                        'Gene start':hit.gene_start,
                        'Gene stop':hit.gene_stop,
                        'Gene Strand':hit.gene_strand,
                        'Position relative to gene':hit.dist_to_gene,
                        'Region':hit.type
                    }
                    writer.writerow(row)

def main():
    #Main algorithm that finds matches for a list of phage genomes and TF motifs and saves them in a dictionary of dictionaries


    result_dic_GLOBAL={}
    for motif_file in tf_filenames:
        print('===== MOTIF FILE: '+ motif_file,' ================' )
        motif_filepath = os.path.join(motifs_path, motif_file)
        motif, threshold, pssm, reverse_pssm=readTF(motif_filepath,pseudo)
        
        result_dic={}
        i=0
        organisms=[]
        for gnome_file in gen_filenames:
            i=i+1
            print('===== GEN FILE: '+ gnome_file,' ================' )
            print('File ',i,'/',len(gen_filenames))
            
            gnome_filepath = os.path.join(gnome_path, gnome_file)

            gb_object,seq,organism,gen_ID=readGB(gnome_filepath)

            if organism not in organisms:
                organisms.append(organism)
                print("LENGTH IN BP: ", len(seq))

                hits, scores=passThreshold(seq,motif,pssm,reverse_pssm,threshold)
                    
                hits, scores=check_overlap(hits,motif,scores)

                print("Hits found: ",len(hits))

                dist_to_genes,genes,classification=dist_classify(hits, seq, gb_object,up_dist,down_dist)

                results=[]
                for hit in hits:
                    if genes:
                        result=HitResult(gen_ID,organism,motif_file,hit,scores[hit],dist_to_genes[hit],classification[hit],genes[hit])
                    else:
                        result=HitResult(gen_ID,organism,motif_file,hit,scores[hit],dist_to_genes,classification,genes)
                    results.append(result)
                
                result_dic[gnome_file]=results
            else:
                print("***Organism already analyzed***")
        result_dic_GLOBAL[motif_file]=result_dic
    return result_dic_GLOBAL, threshold

def get_genomes():
    #Function that finds the genomes in taxonomy under clade of interest and downloads their content using accessions from MillardLab

    descendants = ncbi.get_descendant_taxa(clade,intermediate_nodes=True,rank_limit='genus')

    descendants = ncbi.translate_to_names(descendants)

    for i in range(len(descendants)):
        if 'Candidatus' in descendants[i]:
            descendants[i]=descendants[i].replace('Candidatus ', '')

    input_filepath = os.path.join(gnome_path, gen_file)
    
    tsvfile = open(input_filepath, newline='')
    phages_hosts_df = list(csv.reader(tsvfile, delimiter='\t'))

    fields=phages_hosts_df[0]
    for i in range(len(fields)):
        if fields[i]=='Host':
            host_column=i
            break
        else:
            i+=1

    accessions=[]
    descriptions=[]
    for row in phages_hosts_df:
        if row[host_column] in descendants:
            accessions.append(row[0])
            descriptions.append(row[1])
    
    dir_contents=listdir(gnome_path)
    filenames=[]
    for i in range(len(accessions)):
        filenames.append(accessions[i]+'.gb')
        print("Processing: "+accessions[i])
        if not ((accessions[i]+'.gb') in dir_contents):
            print('--> Downloading: '+accessions[i])
            net_handle = Entrez.efetch(db='nuccore', id=accessions[i],rettype='gbwithparts',retmode='txt')
            genome_record=net_handle.read()
            out_handle=open(gnome_path+'/'+accessions[i]+".gb","w")
            out_handle.write(genome_record)
    
    return filenames

def get_motifs():
    #Reads the directory of TF Motifs files and saves the filenames

    input_filepath = os.path.join(motifs_path, motifs_file)
    csvfile = open(input_filepath, newline='')
    tf_motifs = list(csv.reader(csvfile, delimiter=','))
    
    fields=tf_motifs[0]
    for i in range(len(fields)):
        if fields[i]=='Motif_file':
            filename_column=i
            break
        else:
            i+=1
    
    filenames=[]
    for row in tf_motifs:
        if row[filename_column] != 'Motif_file':
            if row[filename_column] not in filenames:
                filenames.append(row[filename_column])
    return filenames





#---OPEN and read json file with parameters of interest
json_filename='config.json'
out_path,motifs_path, gnome_path, motifs_file,gen_file, clade, up_dist, down_dist,pseudo=readJSON(json_filename)
#---

#---Get files of genomes(online) and motifs(local)
gen_filenames=get_genomes()
tf_filenames=get_motifs()
#---

print("========== FINDING MATCHES ==========")
results, threshold=main()
save_to_csv(results,out_path, threshold)
print("========== FINISHED ==========")