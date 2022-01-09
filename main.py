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
from tqdm import trange
ncbi=NCBITaxa()
#ncbi.update_taxonomy_database()
Entrez.email = 'pedro.garciasu@autonoma.cat'

class HitResult:
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
    f=open(filename)
    print(filename)
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
    motifs_parsed = SeqIO.parse(filename,"fasta")
    
    
    #------Convert into motif object---------

    list_seq=[]
    for record in motifs_parsed:
        list_seq.append(record.seq.upper())
    motif=motifs.create(list_seq)
    tf_id=filename
    
    #------PSSM----------
    pwm=motif.counts.normalize(pseudocounts=pseudo)
    pssm=pwm.log_odds()
    reverse_pssm=pssm.reverse_complement()
    distribution = pssm.distribution(background=None,precision=10**4)
    threshold=distribution.threshold_patser()

    return motif, threshold, pssm, reverse_pssm

def readGB(filename):
    gb_object = SeqIO.read(filename,'genbank')
    seq=gb_object.seq
    organism = gb_object.annotations['organism']
    gen_ID=gb_object.name

    return gb_object,seq, organism, gen_ID

def passThreshold(seq,motif,pssm,reverse_pssm,threshold):
    
    scores={}
    sites={}
    #---------Obtain the sites that pass the threshold using combined scores----------
    hits=[i for i in trange(len(seq)-len(motif),desc='Finding Hits') if log2(2**pssm.calculate(seq[i:i+len(motif)])+2**reverse_pssm.calculate(seq[i:i+len(motif)]))>=threshold]
    #hits=pssm.search(pssm.calculate(seq),threshold)
    if hits:
        for i in hits:
            scores[i]=log2(2**pssm.calculate(seq[i:i+len(motif)])+2**reverse_pssm.calculate(seq[i:i+len(motif)]))
            sites[i]=seq[i:i+len(motif)]
    #for i in trange(len(seq)-len(motif),desc='Finding Hits'):
        #site=seq[i:i+len(motif)]
        #score=pssm.calculate(seq[i:i+len(motif)])
        #reverse_score=reverse_pssm.calculate(seq[i:i+len(motif)])
        #final_score=math.log2(2**pssm.calculate(seq[i:i+len(motif)])+2**reverse_pssm.calculate(seq[i:i+len(motif)]))
        #if math.log2(2**pssm.calculate(seq[i:i+len(motif)])+2**reverse_pssm.calculate(seq[i:i+len(motif)]))>=threshold:
            #scores[i]=math.log2(2**pssm.calculate(seq[i:i+len(motif)])+2**reverse_pssm.calculate(seq[i:i+len(motif)]))
            #hits.append(i)
            #sites[i]=seq[i:i+len(motif)]
    
    return hits,scores

def check_overlap(hits,motif,scores):
    if hits:
        to_delete=[]
        for i in range(len(hits)-1):
            if (hits[i]+len(motif)/2)>=hits[i+1]:
                delete1=hits[i]
                delete2=hits[i+1]
                if scores[delete1]>=scores[delete2]:
                    to_delete.append(i+1)
                    #hits.pop(i)
                    #del scores[delete2]
                else:
                    to_delete.append(i)
                    #hits.pop(i+1)
                    #del scores[delete1]
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
                    #dist2=pos-f.location.start
                dist_aux = abs(dist1)
                if dist_aux<dist_final:
                    dist_final=dist_aux
                    dist_to_genes[pos]=dist1
                    gene_locations[pos]=f.location
                    genes[pos]=f
                #elif dist_aux>dist_final:
                    #break

    #------Classification of hit by their relative location to genes
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

def save_to_csv(results, path):
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
    filepath=os.path.join(path, "results.csv")
    with open(filepath,'a',newline='') as f:
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
    for motif_file in tf_filenames:
        print('=====MOTIF FILE: '+ motif_file,' ================' )
        motif_filepath = os.path.join(motifs_path, motif_file)
        motif, threshold, pssm, reverse_pssm=readTF(motif_filepath,pseudo)
        result_dic={}
        #gen_filenames=["MW546072.gb","KJ206559.gb","JX080301.gb"]
        i=0
        organisms=[]
        for gnome_file in gen_filenames:
            i=i+1
            print('=====GEN FILE: '+ gnome_file,' ================' )
            print('File ',i,'/',len(gen_filenames))
            gnome_filepath = os.path.join(gnome_path, gnome_file)
            #--------Read gb file-----------
            gb_object,seq,organism,gen_ID=readGB(gnome_filepath)
            

            if organism not in organisms:
                organisms.append(organism)
                print("LENGTH IN BP: ", len(seq))
                #==========PROCESSING================#
                hits,scores=passThreshold(seq,motif,pssm,reverse_pssm,threshold)
                    

                #---------Check for overlapping-------------------
                hits, scores=check_overlap(hits,motif,scores)

                print("Hits found: ",len(hits))
                    
                #---------Classification of sites by relative position to genes--
                dist_to_genes,genes,classification=dist_classify(hits, seq, gb_object,up_dist,down_dist)


                #===============RESULTS=====================#

                #----------Saving results in list of HitResult objects--------------
                results=[]
                for hit in hits:
                    if genes:
                        result=HitResult(gen_ID,organism,motif_file,hit,scores[hit],dist_to_genes[hit],classification[hit],genes[hit])
                    else:
                        result=HitResult(gen_ID,organism,motif_file,hit,scores[hit],dist_to_genes,classification,genes)
                    results.append(result)
                
                result_dic[gnome_file]=results
                #if i==10:
                    #break
        
        result_dic_GLOBAL[motif_file]=result_dic
        save_to_csv(result_dic_GLOBAL,out_path)




def get_genomes():
    
    
    #===========Match hosts with descendants of clade===============

    #get descendants of clade
    descendants = ncbi.get_descendant_taxa(clade,intermediate_nodes=True,rank_limit='genus')

    print("===================================================")
    descendants = ncbi.translate_to_names(descendants)

    for i in range(len(descendants)):
        if 'Candidatus' in descendants[i]:
            descendants[i]=descendants[i].replace('Candidatus ', '')

    input_filepath = os.path.join(gnome_path, gen_file)
    
    #open tsv with accessions
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
    #find descendants in hosts column and save accessions for each match
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
    filenames=[]
    #print(dir_contents)
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
            filenames.append(row[filename_column])
    return filenames


json_filename='jsonExample.json'
out_path,motifs_path, gnome_path, motifs_file,gen_file, clade, up_dist, down_dist,pseudo=readJSON(json_filename)
gen_filenames=get_genomes()
tf_filenames=get_motifs()

result_dic_GLOBAL={}
print("=====FINDING MATCHES========")
main()
print("==========FINISHED======")