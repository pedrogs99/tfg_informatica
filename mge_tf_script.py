from contextlib import nullcontext
import math
import Bio
import json
from Bio import SeqIO
from Bio import motifs



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


def read_from_JSON(filename):
    f=open(filename)
    json_obj=json.load(f)
    motifs_file=json_obj['TF_filename']
    gen_file=json_obj['Genome_filename']
    up_dist=json_obj['upstream_dist']
    down_dist=json_obj['downstream_dist']
    pseudo=json_obj['pseudocounts']
    return motifs_file, gen_file,up_dist,down_dist,pseudo


#===========READING FILES===================

#---Read specs from JSON file----------
json_filename='home/user/jsonfilename'
motifs_file, gen_file,up_dist,down_dist,pseudo=read_from_JSON(json_filename)


#-----Read TF file---------
motifs_parsed = SeqIO.parse(motifs_file,'fasta')

#------Convert into motif object---------
list_seq=[]
for record in motifs_parsed:
    list_seq.append(record.seq)
motif=motifs.create(list_seq)
site_length=len(motif)
tf_id='prueba'

#------PSSM----------
pwm=motif.counts.normalize(pseudocounts=pseudo)
pssm=pwm.log_odds()
reverse_pssm=pssm.reverse_complement()
distribution = pssm.distribution(background=None,precision=10**4)
threshold=distribution.threshold_patser()

#--------Read gb file-----------
gb_object = SeqIO.read(gen_file,'genbank')
seq=gb_object.seq
organism = gb_object.annotations['organism']
gen_ID=gb_object.name



#==========PROCESSING================#

hits=[]
scores={}
sites={}
#---------Obtain the sites that pass the threshold using combined scores----------
for i in range(len(seq)-site_length):
    site=seq[i:i+site_length]
    score=pssm.calculate(site)
    reverse_score=reverse_pssm.calculate(site)
    final_score=math.log2(2**score+2**reverse_score)
    if final_score>=threshold:
        scores[i]=final_score
        hits.append(i)
        sites[i]=site

#---------Check for overlapping-------------------
for i in range(len(hits)-1):
    if (hits[i]+site_length/2)>=hits[i+1]:
        delete1=hits[i]
        delete2=hits[i+1]
        if scores[delete1]>=scores[delete2]:
            hits.pop(i)
            del scores[delete2]
        else:
            hits.pop(i+1)
            del scores[delete1]


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
for result in results:
    print(result.location)
    print(result.score)
    print(result.sequence)
    print(result.type)

