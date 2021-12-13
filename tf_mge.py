from contextlib import nullcontext
import os
from os import listdir
from ete3 import NCBITaxa
import math
import json
from Bio import SeqIO
from Bio import motifs
from Bio import Entrez
import csv
ncbi=NCBITaxa()
#ncbi.update_taxonomy_database()

class HitResult:
    def __init__(self, gen_ID, tf_ID, location, site_seq, score, dist_to_gene, type, gene):
        self.genome_ID=gen_ID
        self.tf_ID=tf_ID
        self.location=location
        self.sequence=site_seq
        self.score=score
        self.dist_to_gene=dist_to_gene
        self.type=type
        self.gene=gene


def readJSON(filename):
    f=open(filename)
    print(filename)
    json_obj=json.load(f)
    gnome_path=json_obj['genomes_path']
    motifs_file=json_obj['TF_filename']
    gen_file=json_obj['Genome_filename']
    clade=json_obj['clade']
    up_dist=json_obj['upstream_dist']
    down_dist=json_obj['downstream_dist']
    pseudo=json_obj['pseudocounts']
    return motifs_file,clade, gen_file,up_dist,down_dist,pseudo,gnome_path

def readTF(filename, pseudo):
    motifs_parsed = SeqIO.parse(filename,'fasta')
    
    
    #------Convert into motif object---------
    list_seq=[]
    for record in motifs_parsed:
        list_seq.append(record.seq)
    motif=motifs.create(list_seq)
    tf_id='prueba'
    
    #------PSSM----------
    pwm=motif.counts.normalize(pseudocounts=pseudo)
    pssm=pwm.log_odds()
    reverse_pssm=pssm.reverse_complement()
    distribution = pssm.distribution(background=None,precision=10**4)
    threshold=distribution.threshold_patser()

    return motif, threshold, pssm, reverse_pssm, tf_id

def readGB(filename):
    gb_object = SeqIO.read(filename,'genbank')
    seq=gb_object.seq
    organism = gb_object.annotations['organism']
    gen_ID=gb_object.name

    return gb_object,seq, organism, gen_ID

def passThreshold(seq,motif,pssm,reverse_pssm,threshold):
    hits=[]
    scores={}
    sites={}
    #---------Obtain the sites that pass the threshold using combined scores----------
    for i in range(len(seq)-len(motif)):
        site=seq[i:i+len(motif)]
        score=pssm.calculate(site)
        reverse_score=reverse_pssm.calculate(site)
        final_score=math.log2(2**score+2**reverse_score)
        if final_score>=threshold:
            scores[i]=final_score
            hits.append(i)
            sites[i]=site
    
    return hits,scores,sites

def check_overlap(hits,motif,scores):
    for i in range(len(hits)-1):
        if (hits[i]+len(motif)/2)>=hits[i+1]:
            delete1=hits[i]
            delete2=hits[i+1]
            if scores[delete1]>=scores[delete2]:
                hits.pop(i)
                del scores[delete2]
            else:
                hits.pop(i+1)
                del scores[delete1]
    
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
                elif dist_aux>dist_final:
                    break

    #------Classification of hit by their relative location to genes
        if genes[pos].strand is 1:
            if dist_to_genes[pos]>=-(up_dist) and dist_to_genes[pos]<=down_dist:
                classification[pos]="OPERATOR"
            elif dist_to_genes[pos]>down_dist and dist_to_genes[pos]<gene_locations[pos].end:
                classification[pos]="INTRA"
            elif dist_to_genes[pos]<-(up_dist) or dist_to_genes[pos]>gene_locations[pos].end:
                classification[pos]="INTER"
        elif genes[pos].strand is -1:
            if dist_to_genes[pos] <= up_dist and dist_to_genes[pos] >= -(down_dist):
                classification[pos]="OPERATOR"
            elif dist_to_genes[pos] < -(down_dist) and dist_to_genes[pos] > gene_locations[pos].start:
                classification[pos]="INTRA"
            elif dist_to_genes[pos] > up_dist or dist_to_genes[pos] < gene_locations[pos].start:
                classification[pos]="INTER"
    
    return dist_to_genes,gene_locations,genes,classification

def main():
    #===========READING FILES===================

    #---Read specs from JSON file----------
    json_filename='jsonPrueba.json'
    motifs_file,clade, gen_file,up_dist,down_dist,pseudo,gnome_path=readJSON(json_filename)


    #-----Read TF file---------
    motif, threshold, pssm, reverse_pssm, tf_id=readTF(motifs_file,pseudo)

    #--------Read gb file-----------
    gb_object,seq,organism,gen_ID=readGB(gen_file)

    
    #==========PROCESSING================#
    hits,scores,sites=passThreshold(seq,motif,pssm,reverse_pssm,threshold)
    

    #---------Check for overlapping-------------------
    hits, scores=check_overlap(hits,motif,scores)
    
    #---------Classification of sites by relative position to genes--
    dist_to_genes,gene_locations,genes,classification=dist_classify(hits, seq, gb_object,up_dist,down_dist)


    #===============RESULTS=====================#

    #----------Saving results in list of HitResult objects--------------
    results=[]
    for hit in hits:
        result=HitResult(gen_ID,tf_id,hit,sites[hit],scores[hit],dist_to_genes[hit],classification[hit],genes[hit])
        results.append(result)

    #TODO: saving in file? in correct path. Also, find TF_ID and how to manage clade atribute and paths
    #Right now, the HitResult object contains all the info we want, 
    # including full info about the closest gene to a hit


    #--TESTING----
    print(results[1].gene.location.start)
    print(results[1].sequence)


def analysis():
    
    
    #===========Match hosts with descendants of clade===============

    descendants = ncbi.get_descendant_taxa(clade,intermediate_nodes=True,rank_limit='genus')

    print("===================================================")
    descendants = ncbi.translate_to_names(descendants)

    for i in range(len(descendants)):
        if 'Candidatus' in descendants[i]:
            descendants[i]=descendants[i].replace('Candidatus ', '')
    #print(descendants)
    #descendants.append('Unclassified')

    input_filepath = os.path.join(gnome_path, "1Nov2021_data.tsv")
    
    tsvfile = open(input_filepath, newline='')
    phages_hosts_df = list(csv.reader(tsvfile, delimiter='\t'))

    fields=phages_hosts_df[0]
    for i in range(len(fields)):
        if fields[i]=='Host':
            host_column=i
            break
        else:
            i+=1
    
    print("===================================================")
    accessions=[]
    descriptions=[]
    for row in phages_hosts_df:
        #print(row[host_column])
        if row[host_column] in descendants:
            accessions.append(row[0])
            descriptions.append(row[1])

    

    #print(accessions)
    #print(descriptions)

    #===============DOWNLOAD GENOMES============
    
    dir_contents=listdir(gnome_path)
    #print(dir_contents)
    for i in range(len(accessions)):
        print("Processing: "+accessions[i])
        if not ((accessions[i]+'.gb') in dir_contents):
            print('--> Downloading: '+accessions[i])
            net_handle = Entrez.efetch(db='nuccore', id=accessions[i],rettype='gbwithparts',retmode='txt')
            genome_record=net_handle.read()
            out_handle=open(gnome_path+'/'+accessions[i]+".gb","w")
            out_handle.write(genome_record)

    

    



#main()
json_filename='jsonPrueba.json'
motifs_file,clade, gen_file,up_dist,down_dist,pseudo,gnome_path=readJSON(json_filename)
analysis()