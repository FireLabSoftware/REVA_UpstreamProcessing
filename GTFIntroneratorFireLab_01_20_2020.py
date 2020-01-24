#!/usr/bin/env python -i

## input is gtf annotation file including exons 
## output is gtf file with introns annotated
## input files can be gzip compressed (make sure the suffix is '.gz')
## each combination of exon-start and exon-stop is listed just once.
## Intron id are set based on the exon numbers and transcript names (e.g. "MyTranscript.1.2a.e3")
## Command Line Syntax is "python GTFIntroneratorFire.py <path_of_gtf_file> IntronsOnly=False"
## There is one user-settable parameter, "IntronsOnly".
## Setting IntronsOnly=True yields an output with only introns
## Setting IntronsOnly=False) yields output with introns embedded in the original GTF
## This is designed to work with Python 2.7 and GTF files but should work with Python 3 and may work with some GTF files
## Multiple filenames can be specific in the command line, and wildcards with * should work on Linux/MacOS

from sys import argv
import gzip
InFileList1 = []
KeepCurrentAnnotation1=True
for x in argv[1:]:
    if x.lower()=='intronsonly' or x.lower()=='intronsonly=true':
        KeepCurrentAnnotation1=False
    elif x.lower()=='intronsonly=false':
        KeepCurrentAnnotation1=True
    else:
        InFileList1.append(x)

if not(InFileList1):
    print('No Input File Specified')
    print('Sytax is "'+argv[0]+ ' <path_of_gtf_file>"')
    print('There is one user-settable parameter, "IntronsOnly"')
    print('Setting IntronsOnly=True yields an output file with only introns')
    print('Setting IntronsOnly=False yields an output file with introns embedded in the original GTF')
    exit()          

for fn1 in InFileList1:
    if fn1.endswith('.gz'):
        F = gzip.open(fn1,mode='r')
    else:
        try:
            F = open(fn1,mode='rU')  ## May need to change this to mode='r' for later Python3 versions
        except:
            F = open(fn1,mode='r')
    if '.gtf' in fn1:
        gn0 = fn1.split('.gtf',1)[0]
    else:
        gn0 = fn1.split('.gff',1)[0]
    if KeepCurrentAnnotation1:
        Nubbin1 = "_PlusImpliedIntrons"
    else:
        Nubbin1 = "_JustImpliedIntrons"
    gn1 = gn0+Nubbin1+fn1[len(gn0):]
    if gn1.endswith('.gz'):
        G =gzip.open(gn1, mode='w')        
    else:
        G =open(gn1, mode='w')        
    D0 = {}
    D1 = {}
    for L in F:
        L1 = L.split('\t')
        if len(L1)>=9 and L1[0]!='#' and L1[2]=='exon':
            pos1 = (L1[0],int(L1[3]),int(L1[4]),L1[6])
            AList1 = L1[8].split(';')
            D2 = {}
            for a1 in AList1:
                a1 = a1.strip()
                key1 = a1.split(' ')[0].replace('"','').strip()
                value1= a1.split(' ')[-1].replace('"','').strip()
                D2[key1] = value1
            if ('gene_id' in D2) and ('gene_version' in D2) and ('exon_number' in D2):
                gi1 = D2['gene_id']
                gv1 = D2['gene_version']        
                ti1 = D2['transcript_id']
                en1 = int(D2['exon_number'])
                D1[(gi1,gv1,ti1,en1)] = pos1
                if (gi1,gv1,ti1,en1+1) in D1:
                    pos0 = D1[(gi1,gv1,ti1,en1+1)]
                    posL1 = tuple([pos1[0],]+sorted([pos0[1],pos0[2],pos1[1],pos1[2]])[1:3]+[pos1[3],])
                    if not posL1 in D0:
                        L2 = L1[:]
                        L2[2] = 'intron'
                        L2[3] = str(posL1[1]+1)
                        L2[4] = str(posL1[2]-1)
                        L2[7] = '.'
                        L2[8] = L1[8].replace('exon','intron')
                        IntronID1 = ti1+'.e'+str(en1)
                        L2[8] = L2[8].replace(D2['exon_id'],IntronID1)
                        D0[posL1]='\t'.join(L2)
                        G.write(D0[posL1])        
                if (gi1,gv1,ti1,en1-1) in D1:
                    pos0 = D1[(gi1,gv1,ti1,en1-1)]
                    posL1 = tuple([pos1[0],]+sorted([pos0[1],pos0[2],pos1[1],pos1[2]])[1:3]+[pos1[3],])
                    if not posL1 in D0:
                        L2 = L1[:]
                        L2[2] = 'intron'
                        L2[3] = str(posL1[1]+1)
                        L2[4] = str(posL1[2]-1)
                        L2[7] = '.'
                        L2[8] = L1[8].replace('exon_number "'+str(en1),
                                              'exon_number "'+str(en1-1))
                        L2[8] = L2[8].replace('exon','intron')
                        IntronID1 = ti1+'.e'+str(en1-1)
                        L2[8] = L2[8].replace(D2['exon_id'],IntronID1)
                        D0[posL1]='\t'.join(L2)
                        G.write(D0[posL1])
        if KeepCurrentAnnotation1 or L1[0]=='#':
            G.write(L)
    G.close()

            
## The file assumes the following fields in the ninth column of the gtf file
##gene_id "WBGene00198386"
        ##gene_version "1"
        ##transcript_id, e.g. "cTel3X.3"
        ##exon_number, e.g.,  "1"
        ##gene_name, e.g.,  "cTel3X.3"
        ##gene_source, e.g., "WormBase"
        ##gene_biotype, e.g., "ncRNA"
        ##transcript_name, e.g., "cTel3X.3"
        ##transcript_source , e.g.,"WormBase"
        ##transcript_biotype, e.g., "ncRNA"
        ##exon_id, e.g., "cTel3X.3.e1"

##User sets the names of the input and output file and an additional value "KeepCurrentAnnotation1"
## KeepCurrentAnnotation1=True embeds the intron annotations in the original file
## KeepCurrentAnnotation1=False yields a file with just intron annotations
## Enjoy
## Copywrite 2020, FireLabInformatics, Stanford University.  No guarantees (really).  This is provisional
##    software and likely to cause heartburn in all but the most intrepid sequencer.
## Thanks to Jim and Al for coining the term "intronerator" (and to Phil for discovering introns)

    
    
        
        
    
