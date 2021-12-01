#!/usr/bin/env python -i

## Version PreREVA_zb_12_050821
PreREVA1 = True

## Modules to import
import array
from time import time, strftime, localtime, asctime
try:
    import cPickle
except:
    import pickle as cPickle
import numpy as np
from itertools import chain
try:
    from itertools import imap, izip, izip_longest 
except:
    imap=map
    izip=zip
    from itertools import zip_longest as izip_longest
try:
    temprange1 = xrange(10)
except:
    xrange = range
import gzip
import zipfile
from operator import ne
from sys import version,argv,platform,getsizeof
import os
import subprocess
import re
from glob import glob
from random import choice
def ham1(a,b):
    return(sum(imap(ne,a,b))+abs(len(a)-len(b)))  ## quick hamming distance (with penalties for any mismatch in length

## Constants relevant to the execution of the program
klen1=25  ## the Kmer length used in searches.  Recommended that this be an odd number to avoid any individual k-mer being its own antisense (which can lead to ambiguities in counting) 
BitsIndexed1=25 ## REVA uses a combination of indexing and sorting to keep track of k-mers.  BitsIndexed determines how many bits are accounted for by each.  Default is 25.
RefSeqFile1='No_File_Specified_For_Reference_Please_Specify_Your_FastA_Reference_Sequence_File_In_Command_Line' ##'miminal_hg38_Nsreplaced.fa' #### /Users/firelab08/Desktop/BigData/CircleFinder/ws220.fa' ##'OP50Mock.fa' ## Reference sequence in FastA format
IndexFile1='default' ###'/Users/firelab08/Desktop/BigData/CircleFinder/ws220_kIndex_Kmer25_IndexedBits25_FirstOccurencePlusMultiplicities.hqa2' ## A premade index file that will speed things up
UCSCAssembly1 = 'default'  ## this is, in particular the species identifier for the UCSC browser.  Current values are ce11 for worm and hg# for human
UCSCLinks1 = True ## Setting this to False skips UCSC link output
CircleMax1=100000 ##Maximum Circle Length that we'll look for (everything larger is assumed to be a structural rearrangement)
DeletionMax1=1000000 ##Maximum deletion Length that we'll look for (everything larger is assumed to be a structural rearrangement)
InsertionMin1=4 ## Minimum Circle Length: operationally this is a sequencing error filter.  Although it gets set in principal to the smallest circle.  The smallest simple deletion that could come from sequencing error is fine(4 is pretty much okay) 
Tn5DupMax1=16 ## Maximum tagmentation overlap for detecting singly-tagmented circles (16 should be fine) 
Tn5DupMin1=0  ## Minimum tagmentation overlap for detecting singly-tagmented circles (0 should be fine) 
SeparationGranularity1=10000 ## Granularity of storage for ciclular events (histogram or array)
ReadSeparationMax1=3000 ## Largest fragment we expect to be able to capture and sequence on the flow cell
ReportInterval1=20000 ## Lines between reporting intervals (text messages that indicate progress... not relevant to output of the progrma)
LinesToProcess1= -1 ## 100000 ## number of lines to process in the input files (-1 for everything)
VerboseOutput1=True  #Several extra parameters reported for potential rearranged reads and read pairs
DomainLength1=3000000 ##Length limit for all instances of an individual k-mer.  All instances of the k-mer must be within this limit and on the same chromosome for the repeat to be considered "focal".  Otherwise, the repeat is considered "dispersed".
Short1,Long1 = 100,200 ##Cutoffs in total span for a read pair to be considered short, medium, or long.  Default values are 100 for Short1, 100 for Long1.
MinRepeatEvalue1=float('1.0E-6') ## Minimum e-value for repeat elements (from a DFAM database file) to be considered as identified
CircularChromosomeLimit1=13500000 ## Above this size all chromosomes are considered linear
OtherSVCandidates1 = False ## Report translocation and inversion candidates (junctions between chromosomes or between sense and antisense strands respectively)
test1='TestSVsUnc13R1.fastq.gz'
test2='TestSVsUnc13R2.fastq.gz'
Feature1D={}
FeatureTags1=[]
Feature1 = False
FiL0=[] ##'/Users/firelab08/Desktop/BigData/ReadEvaluation/VC2010SingleWorms/VC2010-1-PCR4_S4_L001_R1_001.fastq.gz']## aardvark ##[]  ## A list of fastq Read files to match
FiD0=[] ##'R1' or 'R2' Indicates that this is read1 or read2
FiN0=[] ## Potential 'Narratives' for each read file.  Default is the file name but can include additional information after a '#' indicator
GFF1 = [] ## List of GTF/GFF data files for coverage calculation
t0=time()
MyCode='default'
Mnemonic1=[]
chrbase1='default'
DFAM1=''
DFAMCount1=False
ai1=0
FullCoverage1=True
BriefCoverage1=True
CoverageAccumulationInterval1=200000 ## Number of read pairs between updates of cumulative coverage zero for no coverage updates
FullCoverageInterval1=0 ## Zero instructs REVA to only calculate coverage once at the end of the cycle, otherwise the program will output aggretgate coverage after each CoverageInterval1 reads
FullCoverageByFile1=False ## True outputs the cumulative coverage after each examined file
KmerBinCoverageUnique1 = True
KmerBinCoverageAll1 = True

KeepDoubleExtension1 = False  ## when looking for 5' and 3' extensions, setting this to false throws away any sequences extended at both ends.  This is to avoid contaminating sequences that might match in one or a few central k-mers
SnpFilter1 = 2 ## this will filter out sequences that have an extended match followed by mismatch followed by match as likely not having extensions.  Setting this to 2 will filter out any apparent extension (5' or 3') which has 2 or more matches following the first mismatch
R1Buffer5 = 0 ## automatically trim this number of bases from the 5' end of each R1 read (after barcode elimination)
R1Buffer3 = 0 ## automatically trim this number of bases from the 3' end of each R1 read (after linker elimination)
R2Buffer5 = 0 ## automatically trim this number of bases from the 5' end of each R2 read (after barcode elimination)
R2Buffer3 = 0 ## automatically trim this number of bases from the 3' end of each R2 read (after linker elimination)

## Reporting of relative read coverage/starts/ends for regions and features 
## What are we counting
ReportReads1 = True ## Report the number of reads with start/end/k-mer in a given region
ReportPositions1 = True ## Report the number of positions in the sequence with start/end/kmer mapping to that position

ReportStarts1 = True ## Report "Starts" (literally extrapolated first match positions identified by first mapping kmer in each read) 
ReportEnds1 = True ## Report "Ends" (literally extrapolated last match positions identified by first mapping kmer in each read)
ReportCoverage1 = True ## Report coverage by every observed k-mer in each read  (not relevant for PreReva)  

##Seq Index Tools (including plurality-- which finds the most common index combination and looks at only reads with that index combination)
SeqIndexMode1 = 0  # 0 for no filtering on sequence index, 1 for plurality [most common], 2 for minority [least common], eventually holds the sequences of one or both required indices as a tuple
SeqIndex1 = ''  ## Narrative list of seq indices to either include (SeqIndexMode=1) or exclude (SeqIndexMode=2)
SeqIndexD1 = {} ## Dictionary of sequence indices to include or exclude (depending on SeqIndexMode1)
def GetSeqIndexFromReadID1(id1):  ## Input is description line in FastQ file, output is the index sequence (upper case)
    return id1.strip().split(':')[-1].upper()  ## assumes index is the bases following the terminal colon in the 

## How to handle two strands
ReportSense1 = True ## Report sense reads separately
ReportAntisense1 = True ## Report antisense reads separately
ReportBothStrands1 = False ## Report sum of sense and antisense reads

## How to handle repeated sequences
ReportUnique1 = True ## Report unique regions
ReportRepeats1 = True ## Report Local, Chromosomal, and Dispersed repeated regions (not relevant for PreReva)

FindStructuralAnomalies1=True ## Output a file of possible structural anomalies either from paired end or split read inference
CombinationList1=0 ## Set CombinationList1 to 1 to make a raw list of read pair start combinations and counts for all positions,
                   ## set CombinationList1 to 2 for just sequences in Base-ByBase output
                   ## set CombinationList1 to 3 for just paired end reads that meet criteria for a Tn5 tagmentation (8-10 base overlap with circular topology) 
RequireKMerOnly1=True ## Setting this to true outputs potential structural variants based only on k-mer matches not imposing a secondary requirement that all other bases on either side of the potential variant match.  Setting this to False is mostly of use with highly accurate reference genomes, accurate sequencing, and relatively short read lengths
MatePairSubstitutionMax1=1 ## When aligning sequences and calling variants using mate pairs how much variation (substitution tolerance) is allowed between read and genome
SplitReadSubstitutionMax1=1 ## when aligning sequences at the ends of possible variants, how many substitutions are allowed
MatePairIndelCountMax1=1 ## When aligning sequences and calling variants using mate pairs how manuy indels are allowed between read and genome
MatePairIndelLengthMax1=3 ## When aligning sequences and calling variants using mate pairs what is the maximum length of indels is allowed between read and genome
SplitReadIndelCountMax1=1 ## when aligning sequences at the ends of possible variants, how many indels are allowed (maximum value, inclusive of the number given)
SplitReadIndelLengthMax1=3 ## when aligning sequences at the ends of possible variants, how many bases of indels are allowed
CoverageDType1=np.uint16
CoverageSEMax1=255   ## maximum coverage that will be recorded in looking for starts and ends in basebybase output.  Setting to 255 conserves memory, setting to higher number will eat up a bit of memory but allow accurate numbers in high coverage regions
MaxCoverage1=2**16-1
MaxCoverHist1=100
MaxMultiplicityReported1=63 ## If multiplicity is reported, should we max out at some value (usually 63).  Setting at 63 results in best memmory usage.  For larger numbers set at 2**(2*n)-1 or 16383 to use two bytes per position, 1073741823 for four bytes and larger for 8 bytes
MultiplicityDataType1=np.uint8
IndexConstructionBins1 = 0
Debug1 = False ## For debug purposes only read and operate with the first chromosome/FastA entry from the reference file
TempFileUID1 = choice('QWERTYUIOPASDFGHJKLZXCVBNMqwertyuiopasdfghjklzxcvbnm')  ## This will be a seven character unique ID for temp files

BarcodeSeqR1 = ''
BarcodeRequireR1 = False
BarcodeSeqR2 = ''
BarcodeRequireR2 = False
BarcodeLenR1 = 0
BarcodeLenR2 = 0
AnyBarcodeR1 = False
AnyBarcodeR2 = False
LinkerR1 = ''
LinkerRequireR1 = False
LinkerR2 = ''
LinkerRequireR2 = False
LinkerLenR1 = 0
LinkerLenR2 = 0

BinByBin1 = True
ChromosomeByChromosome1 = True

BaseByBase1 = False
BaseByBaseOutFile1 = 'default'
BaseByBaseRegionsFiles1 = [] ## will take from infile if available... if a GTF/GFF is specified in the BaseByBaseCall, will use that; if there is no GTF/GFF file in the BaseByBase call but one used in the command line, will use that, otherwise will use all positions.
BaseByBaseCategories1 = {}
BaseByBaseTags1 = []
BaseByBasePositions1 = []
BaseByBaseColumns1 = 'CPFOBRMKSET'
for vL in BaseByBaseColumns1:
    vars()['bbb'+vL+'1'] = True

MaxExtension1 = 4  ## the maximun length of recorded extension for base-by-base analysis
UCSCBuffer1 = 1000 ## Buffer to add to UCSC display around any gene that is shown

FivePrimeExtensionDisallowed1 = False  ## Filtering based on homology match for R1.  This eliminates any sequence where the first base doesn't match
MinHomology1 = 0  ## Allows option to filter based on minimum homology length for R1 default is zero (no filtering)
MaxHomology1 = 999999999 ## Allows option to filter based on maximum homology length for R1 default is 999999999 (essentially no filtering)
StartHomology1 = ''  ## Allows option to filter based on start of homology sequence (e.g. 'G' insists first base of homology is a 'G') default is no filtering

SRRList1 = [] ## A List of SRR IDs for analysis from the short read archive
SRRTempList1 = [] ## A local file name list for SRR (list of files to eventually be deleted)
MetaData1 = True ## True instructs REVA to display all available metadata prominently in each output file.

DeleteNs1 = False ## Seting this to true does a prefiltering of experimental sequence (not reference genome in which all Ns are deleted.  Default is to convert all Ns to Gs.
FirstKTable1 = False ## Report a table with the first position of a unique k-mer in each sequence
FirstKByStrand1 = False ## Report separate firstK values for sense and antisense matches
InterimKTable1 = False ## Do interim K table reports
klen1d = 0 ## This will hold the max value for unmatched 5' or 3' end length (default will be 2*klen1)

R1Only = False ## Report only data from the R1 file?
R2Only = False ## Report only data from the R2 file?

FastQDumpProgram1 = 'fastq-dump' ## Location of a program that will download files from SRA if needed, can be reset by setting FastQDump=<path to fastq-dump executable>.
for i in range(8):
    TempFileUID1 += choice('QWERTYUIOPASDFGHJKLZXCVBNMqwertyuiopasdfghjklzxcvbnm')
ScratchPath1 = os.path.join(os.getcwd(), 'REVATempFiles')
if not(os.path.isdir(ScratchPath1)):
    os.mkdir(ScratchPath1)

def R2FileTest1(n):  ## Look for any occurence of 'r2' in a name with no digit afterwards
    n0 = n.lower()
    if n0.endswith('r2'):
        return True
    if 'r2' in n0:
        n1 = n0.split('r2')
        for n2 in n1[1:]:
            if not(n2[0].isdigit()):
                return True
    return False

## parse command line arguments

CircleMaxSet1 = False
DeletionMaxSet1 = False ## Keep track of whether these have been set and if so don't reset them automatically if OtherSV is set to True
KeepAllSense1 = False ## Use PreREVA to keep every read that has a sense k-mer
KeepAllAntisense1 = False ## Use PreREVA to keep every read that has an antisense k-mer
KeepAllBoth1 = False ## Keep all reads that have both a sense and an antisense k-mer 
KeepAllEither1 = False ## Keep all reads that have either a sense or an antisense k-mer 
KeepAllNeither1 = False ## Keep all reads that have neither a sense nor an antisense k-mer 
while ai1<len(argv):
    a1=argv[ai1]
    ai1+=1
    a1=a1.replace('"','').replace("'","")
    if a1=='choose' or a1.lower()=='-c' or a1.lower()=='c':
        from Tkinter import Tk
        root=Tk()
        root.attributes("-topmost", True)
        from tkFileDialog import askopenfilenames
        R1File1=askopenfilenames(title='Files with R1 reads',initialdir=os.getcwd(), filetypes=[("gz files","*.gz"),("Fastq files","*.fastq")])
        root.destroy()
        FiL0.extend(R1File1.split('#')[0])
        FiN0.extend(R1File1)
        FiD0.extend('R1'*len(R1File1))
        continue
    if a1.lower()=='exon' or a1.lower()=='-e' or a1.lower()=='e':
        Feature1D['exon'] = len(Feature1D)+1
        continue
    if a1.lower()=='gene' or a1.lower()=='-g' or a1.lower()=='g':
        Feature1D['gene'] = len(Feature1D)+1
        continue
    if not('=' in a1) and a1.lower().endswith('.py'):
        MyCode1=a1.strip()
        continue
    if not('=' in a1) and a1.lower().startswith('deleten'):
        DeleteNs1 = True
        continue
    if a1.lower().startswith('r1only') and not('false'  in a1.lower()):
        R1Only = True
        continue
    if a1.lower().startswith('r2only') and not('false'  in a1.lower()):
        R2Only = True
        continue
    if a1.lower().startswith('keepallsense') and not('false'  in a1.lower()):
        KeepAllSense1 = True
        continue
    if a1.lower().startswith('keepallantisense') and not('false'  in a1.lower()):
        KeepAllAntisense1 = True
        continue
    if a1.lower().startswith('keepalleither') and not('false'  in a1.lower()):
        KeepAllEither1 = True
        continue
    if a1.lower().startswith('keepallneither') and not('false'  in a1.lower()):
        KeepAllNeither1 = True
        continue
    if a1.lower().startswith('keepallboth') and not('false'  in a1.lower()):
        KeepAllBoth1 = True
        continue
    if not('=' in a1) and a1.lower().startswith('firstkt'):
        FirstKTable1 = True
        continue
    if not('=' in a1) and (a1.lower().startswith('firstkstr') or a1.lower().startswith('firstkbystr')):
        FirstKByStrand1 = True
        continue
    if not('=' in a1) and a1.lower().startswith('interimkt'):
        InterimKTable1 = True
        continue
    if not('=' in a1) and (a1.lower().endswith('.fa') or a1.lower().endswith('.fasta') or a1.lower().endswith('.fa.gz') or a1.lower().endswith('.fasta.gz')  or a1.lower().endswith('.fa.zip') or a1.lower().endswith('.fasta.zip')):
        RefSeqFile1=a1.strip()
        continue
    if not('=' in a1) and (a1.lower().endswith('.gtf') or a1.lower().endswith('.gtf.gz') or a1.lower().endswith('.gff') or a1.lower().endswith('.gff.gz')  or a1.lower().endswith('.gff3') or a1.lower().endswith('.gff3.gz')):
        GFF1.append(a1.strip())
        continue
    if not('=' in a1) and (a1.lower().endswith('files') or a1.lower().endswith('files.txt')):
        R1File1=a1
        Mnemonic1.append(a1)
        for fx1 in R1File1.read().replace(',',' ').split():
            if os.path.isfile(fx1):
                FiL0.append(fx1.split('#')[0])
                FiN0.append(fx1)
                iD0 = 'R1'
                if R2FileTest1(fx1):
                    if os.path.isfile(fx1.replace('R2','R1')) or os.path.isfile(fx1.replace('r2','r1')):
                        iD0 = 'R2'
                    if os.path.isfile(fx1.replace('R2','R1',1)) or os.path.isfile(fx1.replace('r2','r1',1)):
                        iD0 = 'R2'
                    if os.path.isfile(fx1[::-1].replace('2R','1R',1)[::-1]) or os.path.isfile(fx1[::-1].replace('2r','1r',1)[::-1]):
                        iD0 = 'R2'
                FiD0.append(iD0)
        continue
    if not('=' in a1) and (a1.lower().startswith('srr')):
        for a222 in a1.split(','):
            SRRList1.append(a222)
        continue
    if  not('=' in a1) and a1.lower().startswith('reportread'):
        ReportReads1 = True
        continue
    if  not('=' in a1) and a1.lower().startswith('reportposition'):
        ReportPositions1 = True
        continue
    if  not('=' in a1) and a1.lower().startswith('reportstart'):
        ReportStarts1 = True
        continue
    if  not('=' in a1) and a1.lower().startswith('reportend'):
        ReportEnds1 = True
        continue
    if  not('=' in a1) and a1.lower().startswith('reportcover'):
        ReportCoverage1 = True
        continue
    if  not('=' in a1) and a1.lower().startswith('reportsens'):
        ReportSense1 = True
        continue
    if  not('=' in a1) and a1.lower().startswith('reportantisens'):
        ReportAntisense1 = True
        continue
    if  not('=' in a1) and a1.lower().startswith('reportboth'):
        ReportBothStrands1 = True
        continue
    if  not('=' in a1) and a1.lower().startswith('reportuniqu'):
        ReportUnique1 = True
        continue
    if  not('=' in a1) and a1.lower().startswith('reportrepeat'):
        ReportRepeats1 = True
        continue
    if  not('=' in a1) and a1.lower().startswith('basebybase'):
        BaseByBase1 = True
        continue
    if  not('=' in a1) and a1.lower().startswith('binbybin'):
        BinByBin1 = True
        continue
    if  not('=' in a1) and a1.lower().startswith('chromosomebychromosome'):
        ChromosomeByChromosome1 = True
        continue
    if  not('=' in a1) and a1.lower().startswith('ucsc'):
        UCSCLinks1 = True
        continue
    if not('=' in a1) and (a1.lower().endswith('fastq.gz') or a1.lower().endswith('.fastq')):
        if 'R1' in a1 and not('R2' in a1):
            R1File1=a1.strip().replace(',',' ').split()
            for R1File1a1 in R1File1:
                FiL0.append(R1File1a1.split('#')[0])
                FiN0.append(R1File1a1)
                FiD0.append('R1')
            continue
        if 'R2' in a1 and not('R1' in a1) and not(('SRR2' in a1) and (a1.count('R2')==1)):
            R2File1=a1.strip().replace(',',' ').split()
            for R1File1a1 in R1File1:
                FiL0.append(R2File1a1.split('#')[0])
                FiN0.append(R2File1a1)
                FiD0.append('R2')
            continue
        else:
            R1File1=a1.strip().replace(',',' ').split()
            for R1File1a1 in R1File1:
                FiL0.append(R1File1a1.split('#')[0])
                FiN0.append(R1File1a1)
                FiD0.append('')
            continue
    if not('=' in a1) and (a1.lower().endswith('dfam.hits.gz') or a1.lower().endswith('dfam.hits') or a1.lower()=='dfam' or a1.lower()=='-dfam'):
        if not('.hits' in a1):
            if 'ws220' in RefSeqFile1.lower():
                DFAM1='ce10_dfam.hits.gz'
            elif 'hg38' in RefSeqFile1.lower():
                DFAM1='hg38_dfam.hits.gz'
        else:
            DFAM1=a1
        DFAMCount1=True
        continue
    if not('=' in a1) and a1.lower().endswith('.9rf'):
        IndexFile1=a1.strip()
        continue
    if a1.lower()=='-1' or a1.isdigit():
        LinesToProcess1=int(a1.strip())
        continue
    if a1.lower().startswith('nocover'):
        FullCoverage1=False
        FullCoverageInterval1=0 
        FullCoverageByFile1=False 
        continue
    if a1[0]=='-':
        a11=a1.strip()[1:].lower()
        a22=argv[ai1].strip().replace('"','').replace("'","")
        ai1+=1
    else:
        a11=a1.split('=')[0].strip().lower()
        a22=a1.split('=')[-1].strip()
    if a11.startswith('dfam'):
        DFAM1=a22
        DFAMCount1=True
    elif a11.startswith('gfffile') or a11.startswith('featurefile') or a11.startswith('gtffile') or a11.startswith('gffdatafile') or a11.startswith('gtfdatafile'):
        for a222 in a22.split(','):
            GFF1.append(a222)
    elif a11.startswith('srr') or a11.startswith('sra'):
        for a222 in a22.split(','):
            SRRList1.append(a222)
        continue
    elif a11.startswith('featuretag'):
        for a222 in a22.split(','):
            FeatureTags1.extend(a222.lower().split(','))
    elif a11.startswith('featurecategory') or a11.startswith('gfffeature') or a11.startswith('gtffeature') or a11.startswith('feature'):
        for a222 in a22.split(','):
            Feature1D[a222] = len(Feature1D)+1
    elif a11.startswith('r1buffer5'):
        R1Buffer5=int(a22)
    elif a11.startswith('r1buffer3'):
        R1Buffer3=int(a22)
    elif a11.startswith('r2buffer5'):
        R2Buffer5=int(a22)
    elif a11.startswith('r2buffer3'):
        R2Buffer3=int(a22)
    elif a11.startswith('snpfilter'):
        SnpFilter1=int(a22)
    elif a11.startswith('fastqdump') or a11.startswith('fastq-dump'):
        if os.path.isfile(a22) and not(os.path.isdir(a22)):
            FastQDumpProgram1 = a22
        elif os.path.isdir(a22):
            FastQDumpProgram1 = os.path.join(a22,'fastq-dump')
            if not(os.path.isfile(FastQDumpProgram1)):
                FastQDumpProgram1 = os.path.join(a22,'bin','fastq-dump')
                if not(os.path.isfile(FastQDumpProgram1)):
                    os.environ["PATH"]=os.getenv("PATH")+':'+a22
    elif a11.startswith('r1'):        
        for R1File1 in a22.strip().replace(',',' ').split():
            if '*' in R1File1:
                Mnemonic1.append(a22.replace('*','_all_'))
                for fx1 in glob(R1File1):
                    FiL0.append(fx1)
                    FiN0.append(fx1)
                    FiD0.append('R1')
                    if R2FileTest1(fx1):
                        if os.path.isfile(fx1.replace('R1','R2')):
                            Mnemonic1[-1] = Mnemonic1[-1].replace(('R1',''))
                        if os.path.isfile(fx1.replace('R1','R2'),1):
                            Mnemonic1[-1] = Mnemonic1[-1].replace(('R1',''),1)
                        if os.path.isfile(fx1[::-1].replace('1R','2R')[::-1],1):
                            Mnemonic1[-1] = Mnemonic1[-1][::-1].replace('1R','')[::-1]
            elif R1File1.lower().endswith('files') or R1File1.lower().endswith('files.txt'):
                Mnemonic1.append(a22)
                for fx1 in R1File1.read().replace(',',' ').split():
                    if os.path.isfile(fx1):
                        FiL0.append(fx1.split('#')[0])
                        FiN0.append(fx1)
                        FiD0.append('R1')
            else:           
                FiL0.append(R1File1.split('#')[0])
                FiN0.append(R1File1)
                FiD0.append('R1')
    elif a11.startswith('r2'):
        for R2File1 in a22.strip().replace(',',' ').split():
            if '*' in R2File1:
                if not(a22.replace('*','_all_').replace('R1','R2') in Mnemonic1):
                    Mnemonic1.append(a22.replace('*','_all_'))
                for fx1 in glob(R2File1):
                    FiL0.append(fx1)
                    FiN0.append(fx1)
                    FiD0.append('R2')
            elif R2File1.lower().endswith('files') or R2File1.lower().endswith('files.txt'):
                Mnemonic1.append(a22)
                for fx1 in R2File1.read().replace(',',' ').split():
                    if os.path.isfile(fx1):
                        FiL0.append(fx1.split('#')[0])
                        FiN0.append(fx1)
                        FiD0.append('R1')
            else:           
                FiL0.append(R2File1.split('#')[0])
                FiN0.append(R2File1)
                FiD0.append('R2')
    elif a11.startswith('ref'):
        RefSeqFile1=a22
    elif a11.startswith('indexfile'):
        IndexFile1=a22
    elif a11.startswith('seqindex'):
        if a22.lower().startswith('p'):
            SeqIndexMode1 = 1   ## plurality -- keeps only read pairs with the most prevalent combination of 5' and 3' indices
        elif a22.lower().startswith('m'):
            SeqIndexMode1 = 2   ## minority -- keeps only read pairs with something other than the most prevalent combination of 5' and 3' indices
        else:
            for SeqIndex1s in a22.split('),('):
                SeqIndex01 = SeqIndex1s.split(',')[0].strip().replace('"','').replace("'",'').replace(")",'').replace("(",'')
                if len(a22.split(','))>1:
                    SeqIndex02 = SeqIndex1s.split(',')[1].strip().replace('"','').replace("'",'').replace(")",'').replace("(",'')
                else:
                    SeqIndex02 = ''
                SeqIndexD1[(SeqIndex01.upper(),SeqIndex02.upper())] = 0        
    elif a11.startswith('bit'):
        BitsIndexed1=int(a22)
    elif a11.startswith('minre'):
        MinRepeatEvalue1=int(a22)
    elif a11.startswith('minhomol'):
        MinHomology1=int(a22)
    elif a11.startswith('maxhomol'):
        MaxHomology1=int(a22)
    elif a11.startswith('coveragesemax'):
        if int(a22)==0 or int(a22)==-1:
            CoverageSEMax1=2**32-1
        elif int(a22)<=64:
            CoverageSEMax1=2**int(a22)-1        
        else:
            CoverageSEMax1=int(a22)
    elif a11.startswith('circlemax'):
        CircleMax1=int(a22)
    elif a11.startswith('deletionmax'):
        DeletionMax1=int(a22)
    elif a11.startswith('circlemin'):
        InsertionMin1=int(a22)
    elif a11.startswith('insertionmin'):
        InsertionMin1=int(a22)
    elif a11.startswith('tn5dupmax'):
        Tn5DupMax1=int(a22)
    elif a11.startswith('tn5dupmin'):
        Tn5DupMin1=int(a22)
    elif a11.startswith('ucscbuffer'):
        UCSCBuffer1=int(a22)
    elif a11.startswith('report'):
        ReportInterval1=int(a22)
    elif a11.startswith('separation') or a11.startswith('gran'):
        SeparationGranularity1=int(a22)
    elif a11.startswith('readsep'):
        ReadSeparationMax1=int(a22)
    elif a11.startswith('domain'):
        DomainLength1=int(a22)
    elif a11.startswith('short'):
        Short1=int(a22)
    elif a11.startswith('long'):
        Long1=int(a22)
    elif a11.startswith('chrbase'):
        chrbase1=a22
    elif a11.startswith('starthomology'):
        StartHomology1 = a22.strip('"').strip("'").upper()
    elif a11.startswith('basebybasecolumns'):
        BaseByBase1 = True
        bbOptionList1 = ['C','P','F','O','B','R','M','K','S','E','T']
        for bbOs in bbOptionStrings1:
            if not('.' in bbOs):
                for vL in bbOptionList1:
                    vars()['bbb'+vL+'1'] = False
                for vS in bbOptionString1:
                    vars()['bbb'+vS+'1'] = True
    elif a11.startswith('basebybaseregions') or a11.startswith('basebybasefile') or a11.startswith('basebybasefea'):
        BaseByBase1 = True
        BaseByBaseRegionsFiles1.extend(a22.split(','))
    elif a11.startswith('basebybasecategor'):
        BaseByBase1 = True
        for a222 in a22.split(','):
            BaseByBaseCategories1[a222] = len(BaseByBaseCategories1)+1
    elif a11.startswith('basebybasetag'):
        BaseByBase1 = True
        for a222 in a22.split(','):
            BaseByBaseTags1.append(a222.strip('"').strip("'"))
    elif a11.startswith('basebybaseposition'):
        BaseByBase1 = True
        BaseByBasePositions1 =  a22.split(',')                     
    elif a11.startswith('basebybase') :
        if a22.lower().startswith('f'):
            BaseByBase1=False
        else:
            BaseByBase1=True
    elif a11.startswith('binbybin') :
        if a22.lower().startswith('f'):
            BinByBin1=False
        else:
            BinByBin1=True
    elif a11.startswith('othersv') :
        if a22.lower().startswith('f'):
            OtherSVCandidates1=False
        else:
            OtherSVCandidates1=True
            if not(CircleMaxSet1):
                CircleMax1 = 2**63
            if not(DeletionMaxSet1):
                DeletionMax1 = 2**63
    elif a11.startswith('chromosomebychromosome') :
        if a22.lower().startswith('f'):
            ChromosomeByChromosome1=False
        else:
            ChromosomeByChromosome1=True
    elif a11.startswith('reportread') :
        if a22.lower().startswith('f'):
            ReportReads1=False
        else:
            ReportReads1=True
    elif a11.startswith('reportposition') :
        if a22.lower().startswith('f'):
            ReportPositions1=False
        else:
            ReportPositions1=True
    elif a11.startswith('reportstart') :
        if a22.lower().startswith('f'):
            ReportStarts1=False
        else:
            ReportStarts1=True
    elif a11.startswith('reportend') :
        if a22.lower().startswith('f'):
            ReportEnds1=False
        else:
            ReportEnds1=True
    elif a11.startswith('reportcover') :
        if a22.lower().startswith('f'):
            ReportCoverage1=False
        else:
            ReportCoverage1=True
    elif a11.startswith('reportsens') :
        if a22.lower().startswith('f'):
            ReportSense1=False
        else:
            ReportSense1=True
    elif a11.startswith('reportanti') :
        if a22.lower().startswith('f'):
            ReportAntisense1=False
        else:
            ReportAntisense1=True
    elif a11.startswith('reportboth') :
        if a22.lower().startswith('f'):
            ReportBothStrands1=False
        else:
            ReportBothStrands1=True
    elif a11.startswith('reportuniqu') :
        if a22.lower().startswith('f'):
            ReportUnique1=False
        else:
            ReportUnique1=True
    elif a11.startswith('reportrepeat') :
        if a22.lower().startswith('f'):
            ReportRepeats1=False
        else:
            ReportRepeats1=True
    elif a11.startswith('fullcover') :
        if a22.lower().startswith('f'):
            FullCoverage1=False
        else:
            FullCoverage1=True
            FullCoverageByFile1=False 
    elif a11.startswith('gffstart') or a11.startswith('featurestart') or a11.startswith('gtfstart') :
        ReportStarts1=True
        if a22.lower().startswith('f'):
            ReportStarts1=False
    elif a11.startswith('featureend') or a11.startswith('gffend') or a11.startswith('gtfend') :
        ReportEnds1=True
        if a22.lower().startswith('f'):
            ReportEnds1=False
    elif a11.startswith('featurerepeat') or a11.startswith('gffrepeat') or a11.startswith('gtfrepeat') :
        ReportRepeats1=True
        if a22.lower().startswith('f'):
            ReportRepeats1=False
    elif a11.startswith('briefcover'):
        BriefCoverage1=True
        if a22.lower().startswith('f'):
            BriefCoverage1=False
    elif a11.startswith('metadata'):
        MetaData1=True
        if a22.lower().startswith('f'):
            MetaData1=False
    elif a11.startswith('deleten'):
        DeleteNs1=True
        if a22.lower().startswith('f'):
            DeleteNs1=False
    elif a11.startswith('firstkt'):
        FirstKTable1=True
        if a22.lower().startswith('f'):
            FirstKTable1=False
    elif a11.startswith('interimkt'):
        InterimKTable1=True
        if a22.lower().startswith('f'):
            InterimKTable1=False
    elif a11.startswith('firstkstr') or a11.startswith('firstkbystr') :
        FirstKByStrand1=True
        if a22.lower().startswith('f'):
            FirstKByStrand1=False
    elif a11.startswith('fiveprimeex'):
        FivePrimeExtensionDisallowed1=True
        if a22.lower().startswith('f'):
            FivePrimeExtensionDisallowed1=False
    elif a11.startswith('keepdouble'):
        KeepDoubleExtension1=True
        if a22.lower().startswith('f'):
            KeepDoubleExtension1=False
    elif a11.startswith('findstructure') :
        FindStructuralAnomalies1=True
        if a22.lower().startswith('f'):
            FindStructuralAnomalies1=False
    elif a11.startswith('ucsclinks') :
        UCSCLinks1=True
        if a22.lower().startswith('f'):
            UCSCLinks1=False
    elif a11.startswith('combination'):
        if a22.lower().startswith('a') or a22.lower().startswith('c') or a22.lower().startswith('1'):  ## 'all' or no entry -- record all paired read locations
            CombinationList1=1
        elif a22.lower().startswith('b') or a22.lower().startswith('i') or a22.lower().startswith('2'):  ## use base-by-base positions to record
            CombinationList1=2
            BaseByBase1 = True
        elif a22.lower().startswith('t') or a22.lower().startswith('3'):  ## only record for potential Tn5 insertions with 8-10b duplication
            CombinationList1=3
    elif a11.startswith('requirek') or a11.startswith('kmerrequire') :
        RequireKMerOnly1=True
        if a22.lower().startswith('f'):
            RequireKMerOnly1=False
    elif a11.startswith("verbose"):
        VerboseOutput1=True
        if a22.lower().startswith('f'):
            VerboseOutput1=False
    elif a1.lower().startswith('coveragebyfile'):        
        FullCoverageByFile1=True 
        if a22.lower().startswith('f'):
            FullCoverageByFile1=False
    elif a1.lower().startswith('kmerbincoverageunique'):
        KmerBinCoverageUnique1= True
        if a22.lower().startswith('f'):
            KmerBinCoverageUnique1=False
    elif a1.lower().startswith('kmerbincoverageall'):
        KmerBinCoverageAll1 = True
        if a22.lower().startswith('f'):
            KmerBinCoverageAll1=False
    elif a11.startswith('barcoderequire') or a11.startswith('requirebarcode') :
        if '2' in a11:
            if 'f' in a22.lower():
                BarcodeRequireR2=False
            else:
                BarcodeRequireR2=True
        else:
            if 'f' in a22.lower():
                BarcodeRequireR1=False
            else:
                BarcodeRequireR1=True
    elif a11.startswith('barcode'):
        if '2' in a11:
            BarcodeR2 = a22.upper()
            if BarcodeR2.count('N')==len(BarcodeR2):
                AnyBarcodeR2 = True
            BarcodeLenR2 = len(BarcodeR2)
        else:
            BarcodeR1 = a22.upper()
            if BarcodeR1.count('N')==len(BarcodeR1):
                AnyBarcodeR1 = True
            BarcodeLenR1 = len(BarcodeR1)
    elif a11.startswith('linkerrequire') or a11.startswith('requirelinker')  :
        if '2' in a11:
            if 'true' in a22.lower():
                LinkerRequireR2=True
            elif 'f' in a22.lower():
                LinkerRequireR2=False
        else:
            if 'true' in a22.lower():
                LinkerRequireR1=True
            if 'f' in a22.lower():
                LinkerRequireR1=False
    elif a11.startswith('linker')  :
        if '2' in a11:
            LinkerR2 = a22.upper()
            LinkerLenR2 = len(LinkerR2)
        else:
            LinkerR1 = a22.upper()
            LinkerLenR1 = len(LinkerR1)
    elif a11.startswith('matepairsub') :
        MatePairSubstitutionMax1=int(a22)
    elif a11.startswith('matepairindelcount') :
        MatePairIndelCountMax1=int(a22)
    elif a11.startswith('matepairindellength') :
        MatePairIndelLengthMax1=int(a22)
    elif a11.startswith('splitreadsub') :
        SplitReadSubstitutionMax1=int(a22)
    elif a11.startswith('splitreadindelcount') :
        SplitReadIndelCountMax1=int(a22)
    elif a11.startswith('splitreadindellength') :
        SplitReadIndelLengthMax1=int(a22)
    elif a11.startswith('maxexten') :
        MaxExtension1=int(a22)
    elif a11.startswith('lines'):
        LinesToProcess1=int(a22)
    elif a11.startswith('firstkmax') or a11.startswith('maxfirstk'):
        klen1d=int(a22)
    elif a11.startswith('ucscassembl'):
        UCSCAssembly1=a22
    elif a11.startswith('k'):
        klen1=int(a22)
    elif a11.startswith('multiplicity') :
        MaxMultiplicityReported1=int(a22)
    elif a11.startswith('coverageaccumulationinterval'):
        CoverageAccumulationInterval1=int(a22)
    elif a11.startswith('fullcoverageinterval'):
        FullCoverageInterval1=int(a22)
        FullCoverageByFile1=True
        FullCoverage1=True
    elif a11.startswith('maxcoverhist'):
        MaxCoverHist1=int(a22)
    elif a11.startswith('fullcoveragebyfile'):
        CoverageAccumulationInterval1=int(a22)
        FullCoverageByFile1=True
        FullCoverage1=True
    elif a11.startswith('circularchromosomelimit'):
        CircularChromosomeLimit1=int(a22)
    elif a11.startswith('coveragedtype'):
        dt11=np.dtype(a22.split('.')[-1])
        CoverageDType1=dt11
    elif a11.startswith('indexcon'):
        IndexConstructionBins1=int(a22)
    elif a11.startswith('debug'):
        Debug1=True

if PreREVA1:
    ReportRepeats1 = False  ## No repeat data to report for PreREVA, so repeat columns are not useful
    ReportCoverage1 = False

if len(Feature1D)>0:
    Feature1All = False
else:
    Feature1All = True
RefAbbrev1=os.path.basename(RefSeqFile1).split('.')[0]    
if IndexFile1=='default':
    IndexFile0='PreREVAIndex_'+RefAbbrev1+'_k'+str(klen1)+'_b'+str(BitsIndexed1)+'_wv2.9rf'
    if Debug1:
        IndexFile0 = 'Debug'+IndexFile0
    if os.path.isfile(IndexFile0):
        IndexFile1=IndexFile0
    else:
        IndexFile1=os.path.join( os.path.dirname(RefSeqFile1) , IndexFile0 ) ## try to find/make the index in the folder with the Reference Sequence File
        if not ( os.path.isfile(IndexFile1) ):
            try:
                TryIndexOpen1=open(IndexFile1,mode='wb')  ##make sure we have the ability to open a file at that location, otherwise default to the current folder
                TryIndexOpen1.close()
            except:
                IndexFile1=IndexFile0
                    
if chrbase1=='default':
    chrbase1='chr'

pypy1=False
if 'pypy' in version.lower():
    pypy1=True

if SRRList1:
    for srr1 in SRRList1:
        FiL0.append(srr1.split('#')[0])
        FiN0.append(srr1)
        FiD0.append('R1')
        
if len(FiL0)>0 and not(Mnemonic1):
    for fil0,fin0 in zip(FiL0,FiN0):
        mn0 = os.path.basename(fil0).split('.')[0]
        if '#' in fin0:
            mn0 += fin0.split('#',1)[1]
        mn1 = mn0.replace('R1','R2')
        mn2 = mn0.replace('R2','R1')
        if mn1 in Mnemonic1:
            mi1 = Mnemonic1.index(mn1)
            Mnemonic1[mi1] = mn0.replace('R1','R1R2')
        elif mn2 in Mnemonic1:
            mi2 = Mnemonic1.index(mn2)
            Mnemonic1[mi2] = mn0.replace('R1','R1R2')
        elif not(mn0 in Mnemonic1):
            Mnemonic1.append(mn0)
if Mnemonic1:
    Mnemonic1='_'.join(Mnemonic1)
else:
    Mnemonic1='Diagnostic'
OutFileNameBase='PreREVA_'+Mnemonic1+'_'+RefAbbrev1+'_'+strftime("D_%m_%d_%y_T_%H_%M_%S",localtime())

LogFile1="LogSummary_"+OutFileNameBase+'.tdv'
def LogNote1(note,LogFileName):
    LogFile=open(LogFileName,mode='a')
    LogFile.write(note+'\t'+'; Time='+"{0:.2f}".format(time()-t0)+' sec'+'\t'+strftime("D_%m_%d_%y_T_%H_%M_%S",localtime())+' \n')
    LogFile.close()
    print(note.split('#')[0].replace('\t',' ').strip(',') + '; Time='+"{0:.2f}".format(time()-t0)+' sec')

def HeaderTranspose(hT1):
    hT2 = hT1.split('\t')
    hT0 = '<!--\tOutput_Key\t\t-->\n'
    hT0 += '<!--\tColumnNumber\tColumnHeader\t-->\n'
    for iT2,nT2 in enumerate(hT2):
        nT2 = nT2.strip()
        if not(nT2.startswith('<!')):
            hT0 += '<!--\t'+str(iT2)+'\t'+nT2+'\t-->\n'
    hT3 = argv
    if MetaData1:
        hT4 = '<!--\tRunning REVA With Parameters\t-->\n'
        for hT5 in hT3:
            if hT5:
                hT4 += '<!--\t'+hT5+'\t\t-->\n'
        hT4 += '<!--\t\t\t-->\n'
    else:
        hT4 =''
    return hT4+hT0+'<!--\t\t\t-->\n'
   
def FileInfo1(FileID):
    if type(FileID)==str:
        Fn1=FileID
    else:
        Fn1=FileID.name
    s11=os.stat(Fn1)
    return ','.join([Fn1,
                     '#',
                      'Path='+os.path.abspath(Fn1),
                        'Size='+str(s11[6]),
                        'Accessed='+strftime("%m_%d_%y_at_%H_%M_%S",localtime(s11[7])),
                        'Modified='+strftime("%m_%d_%y_at_%H_%M_%S",localtime(s11[8])),
                        'Created='+strftime("%m_%d_%y_at_%H_%M_%S",localtime(s11[9])),
                        'FullInfo='+str(s11)])
LogNote1("Running PreREVA with parameters:"+' '.join(argv)+' #Python Version'+version,LogFile1)

if not 'pypy' in version.lower():
    LogNote1('********',LogFile1)
    LogNote1('WARNING PreREVA is designed to run with the PyPy interpreter.',LogFile1)
    LogNote1('Now running with version "'+version+'", which may be much slower.',LogFile1)
    LogNote1('Please see REVA documentation for instructions on installing PyPy and NumPy.',LogFile1)
    LogNote1('********',LogFile1)

def LogOpeningFile1(FileID):
    LogNote1("Opening "+FileInfo1(FileID),LogFile1)
def LogClosingFile1(FileID):
    LogNote1("Closed "+FileInfo1(FileID),LogFile1)
def LogRunningFile1(FileID):
    LogNote1("Running "+FileInfo1(FileID),LogFile1)
LogRunningFile1(MyCode1)

VirtualFPath1 = {}
VirtualDtype1 = {}
def npush(Variable0):
    if type(Variable0) ==str:
        VList = [Variable0,]
    else:
        VList = Variable0
    for Variable1 in VList:
        if not(Variable1 in VirtualFPath1):
            fp1 = os.path.join(ScratchPath1,Variable1+'_'+TempFileUID1+'_reva.tmp')
            VirtualFPath1[Variable1] = fp1
            if type(globals()[Variable1]) == np.ndarray:
                VirtualDtype1[Variable1] = globals()[Variable1].dtype
                open(fp1,mode='wb').write(globals()[Variable1])                
            else:
                pf0 = open(fp1, mode='wb')
                pp0 = cPickle.Pickler(pf0)
                pp0.dump(globals()[Variable1])
                pf0.close()
                VirtualDtype1[Variable1] ='pickle'
        globals()[Variable1] = None

def npull(Variable0):
    if type(Variable0) ==str:
        VList = [Variable0,]
    else:
        VList = Variable0
    for Variable1 in VList:
        fp1 = VirtualFPath1[Variable0]
        dt1 = VirtualDtype1[Variable0]
        if  type(dt1)==str:
            pf0 = open(fp1, mode='rb')
            pp0 = cPickle.Unpickler(pf0)
            globals()[Variable1] = pp0.load()
            pf0.close()
        else:
            globals()[Variable1] =  np.frombuffer(open(fp1,mode='rb').read(),dtype=VirtualDtype1[Variable1])

if 'darwin' in platform:
    LogNote1('Trying to Run \'caffeinate\' on MacOSX to prevent the system from sleeping',LogFile1)
    try:
        Coffee_process1=subprocess.Popen('caffeinate')
    except:
        LogNote1("Couldn't start 'caffeinate', you may need to manually set your mac to avoid dozing (System Preferences, Energy)",LogFile1)
else:
    LogNote1('Trying to Run \'caffeine\' to prevent the system from sleeping',LogFile1)
    try:
        Coffee_process1=subprocess.Popen('caffeine')
    except:
        LogNote1("Couldn't start 'caffeine', you may need to manually set your mac to avoid dozing (System Preferences, Energy)",LogFile1)
        
BitMask0=4**klen1-1   
BitsSorted1=klen1*2-BitsIndexed1
Bins1=2**BitsIndexed1
BitMask1=2**BitsSorted1-1

AllBase1=['G','A','T','C']
Numbase1={AllBase1[0]:0,AllBase1[1]:1,AllBase1[2]:2,AllBase1[3]:3,'N':0}
## e4 is a set of 4**x exponents
e4=[]
ne4=[]
for i in range(32):
    e4.append(4**i)
    ne4.append(4**i-1)
ex2=np.array([2**i for i in range(64)],dtype=np.uint64)
## Purge Files from a directory
def PurgeREVATempFiles(dir1):
    LogNote1('Starting Purge of Unused Files from Temporary Directory '+ScratchPath1,LogFile1)
    try:
        for f1 in os.listdir(dir1):
            if f1.endswith('_'+TempFileUID1+'_reva.tmp') or f1.endswith('_'+TempFileUID1+'.fasta.gz'):
                os.remove(os.path.join(ScratchPath1,f1))
        LogNote1('Finished Purge of Unused Files from Temporary Directory '+ScratchPath1,LogFile1)
    except Exception as e:
        LogNote1('Purge of temporary files may have failed.  Reason='+str(e)+' --- you may want to do this manually when program is finished.  Directory='+ScratchPath1,LogFile1)

## Open and unpack the index of unique sequences (three sorted lists of integers with 2*klen1 bit list of values-converted-to binary, sorted, then two lists that are the source sequence (chromosome for C. elegna) and position
## fastfind is the routine that returns the index for any given binary sequence k-mer representation
AntisenseB1={'A':'T','T':'A','G':'C','C':'G'}
Tr2=''  ## Tr2 allows a quick translation of sequence to a numerical array.
for i in range(256):
    if chr(i) in 'AGCTagct':
        Tr2+=chr(Numbase1[chr(i).upper()])
    else:  ##in case there are unusual characters in the sequence, they will be converted to Gs
        Tr2+=chr(Numbase1['N'])
TrN1=''  ## Tr2 allows a quick translation of sequence to a numerical array.
for i in range(256):
    if chr(i) in 'AGCTagct':
        TrN1+=chr(1)
    else:  ##in case there are unusual characters in the sequence, they will be converted to Gs
        TrN1+=chr(0)
def seqnum(s):
    ''' converts an input sequence s (upper case string) into a numpy array of values from 0 to 3.  Translation between sequence and numbers from AllBase1'''
    return np.frombuffer(s.translate(Tr2).encode('ascii'),dtype=np.uint8)
def seqNotN(s):
    ''' converts an input sequence s (upper case string) into a numpy array of values from 0 to 3.  Translation between sequence and numbers from AllBase1'''
    return np.frombuffer(s.translate(TrN1).encode('ascii'),dtype=np.bool)

## antisense-- returns the reverse complement of argument string 
filterminus=''
ASB11="AaCcNn*nNgGtT"
for i in range(256):
    if chr(i) in ASB11:
        filterminus+=ASB11[12-ASB11.find(chr(i))]
    else:
        filterminus+='N'
def antisense(s):
    '''return an antisense and filtered version of any sequence)'''
    return s.translate(filterminus)[::-1]

BitDepthD1={np.dtype('bool'):1,
    np.dtype('int8'):1,
    np.dtype('int16'):2,
    np.dtype('int32'):4,
    np.dtype('int64'):8,
    np.dtype('uint8'):1,
    np.dtype('uint16'):2,
    np.dtype('uint32'):4,
    np.dtype('uint64'):8,
    np.dtype('float16'):2,
    np.dtype('float32'):4,
    np.dtype('float64'):8,
    np.dtype('complex64'):8,
    np.dtype('complex128'):16}

def AltFromFile(F_file,dtype,count):  ## PyPy Numpy has no .fromfile method for arrays, so I made this very simple equivalent for 1-D arrays.  Could fail if the index was written on another platform 
    ## F_file needs to be open with 'wb' mode, 'dtype is a numpy data type, and count is the number of entries to read)
    h=F_file.read( BitDepthD1[dtype] * count )
    return np.frombuffer(h,dtype=dtype)
def AltToFile (Array1,F_file):  ## PyPy Numpy has no .fromfile method for arrays, so I made this very simple equivalent for 1-D arrays.  Could fail if the index was written on another platform 
    ## F_file needs to be open with 'wb' mode, Array1 is a numpy array.  May fail if file is written on one platform and read on another)
    F_file.write(Array1)

dtype1=np.uint16
if klen1>8:
    dtype1=np.uint32
if klen1>16:
    dtype1=np.uint64

b_a1 = 12  ## bit length for individual counts of events.  2**b1 is the maximum sequence length for analysis
         ## for anything like 'current' technology (2018) the implicit maximum of 4096 bases per read seems quite adequate
b_a2 = b_a1*2
b_a3 = b_a1*3
b_a4 = b_a1*4
b_a5 = b_a1*5
   
LogNote1('starting fastA file read',LogFile1)
LogOpeningFile1(RefSeqFile1)
if RefSeqFile1.endswith('.gz'):
    if version.startswith('2.'):
        F0=gzip.open(RefSeqFile1,mode='r')
    else:
        F0=gzip.open(RefSeqFile1,mode='rt')
elif RefSeqFile1.endswith('.zip'):
    ZipFile1 = zipfile.ZipFile(RefSeqFile1)  ##.readline()
    F0 = []
    for zfn1 in ZipFile1.namelist():
        if zfn1[0]!='_':
            F0.extend(ZipFile1.open(zfn1,mode='r').read().splitlines())
else:
    if version.startswith('2.'):
        F0=open(RefSeqFile1,mode='rU')
    else:
        F0=open(RefSeqFile1,mode='r')
SD1=[] ## Sense Sequences
AD1=[] ## AntiSense Sequences
LD1=[] ## Length of each sequence
NameA1=[] ## Names of each reference DNA entity by number
NameD1={} ## Numbers for each reference DNA entity by name
CircD1=[] ## Cicular status of each input DNA entity by number if anything less than 13.5MB is cir
for L0 in F0:
    L1=L0.strip()
    if L1.startswith('>'):
        if len(SD1)>0 and Debug1:  ## DEbug mode-- only look at first chromosomal unit in FastA reference file
            break
        NameA1.append(L1[1:].strip().split()[0])
        if "circular" in L1.lower():
            CircD1.append(True)
        elif "linear" in L1.lower():
            CircD1.append(False)
        else:
            CircD1.append(-1)            
        sn1=NameA1[-1]
        NameD1[sn1]=len(NameA1)-1
        SD1.append([])
    else:
        SD1[-1].append(L1.upper().replace('U','T'))
for sn1 in range(len(SD1)):
    SD1[sn1]=''.join(SD1[sn1])   ##SD1[chr][1] is the sense sequence of chr1, SD1.  N's will be converted to "G" to preserve absolute position.
    AD1.append(antisense(SD1[sn1]))
    LD1.append(len(SD1[sn1]))
    if LD1[sn1]<CircularChromosomeLimit1:
        if CircD1[sn1] == -1: 
            CircD1[sn1] = True
    else:
        if CircD1[sn1] == -1: 
            CircD1[sn1] = False
    if CircD1[sn1]:
        SD1[sn1]=SD1[sn1]+SD1[sn1][:klen1-1]
        AD1[sn1]=AD1[sn1]+AD1[sn1][:klen1-1]
        
TotalRefBases1=sum(LD1)
LongestContig1=max(LD1)
dtype2=np.int32  ## data type (signed) for any reference to positions (sense + antisense-) in the individual chromosomes 
if LongestContig1>=2**31-2*klen1-2*Tn5DupMax1:  ## 2*klen1+2*Tn5DupMax1 is a "safety factor to make sure there is no overflow
    dtype2=np.int64
dtype4=np.uint16
if BitsSorted1>16:
    dtype4=np.uint32
if BitsSorted1>32:
    dtype4=np.uint64
dtype5=np.uint16
if BitsIndexed1>16:
    dtype5=np.uint32
if BitsIndexed1>32:
    dtype5=np.uint64

NumSeq1=len(NameA1)
NameARange0=range(NumSeq1)
NameARange1=range(NumSeq1+1)
if not(RefSeqFile1.endswith('zip')):
    F0.close()
    LogClosingFile1(F0)

if IndexConstructionBins1 == 0:
    if TotalRefBases1<120000000:
        IndexConstructionBins1 = 1
    elif TotalRefBases1<400000000:
        IndexConstructionBins1 = 4
    elif TotalRefBases1<1600000000:
        IndexConstructionBins1 = 16
    else:
        IndexConstructionBins1 = 64        
qssD = 1+ (4**klen1//IndexConstructionBins1)

LogNote1('finished fastA file read',LogFile1)
LogNote1('Getting ready to try unpickle/unpack of Index File: '+IndexFile1,LogFile1)
try:
    pf1=open(IndexFile1,mode='rb')
    LogOpeningFile1(pf1)
    pp1=cPickle.Unpickler(pf1)
    RefSeqFileR1=pp1.load()
    klenR1=pp1.load()
    ulenR1=pp1.load()
    V1dtype=pp1.load()
    P1dtype=pp1.load()
    I1dtype=pp1.load()
    NameA1=pp1.load()
    FastIndexDataType1=pp1.load()
    if pypy1:
        Sort_V1dedup=AltFromFile(pf1,dtype=V1dtype,count=ulenR1)
        Sort_P1dedup=AltFromFile(pf1,dtype=P1dtype,count=ulenR1+1)
        Sort_I1dedup=AltFromFile(pf1,dtype=I1dtype,count=ulenR1+1)
        FastIndex1=AltFromFile(pf1,dtype=FastIndexDataType1,count=2**BitsIndexed1+1)
    else:
        Sort_V1dedup=np.fromfile(pf1,dtype=V1dtype,count=ulenR1)
        Sort_P1dedup=np.fromfile(pf1,dtype=P1dtype,count=ulenR1+1)
        Sort_I1dedup=np.fromfile(pf1,dtype=I1dtype,count=ulenR1+1)
        FastIndex1=np.fromfile(pf1,dtype=FastIndexDataType1,count=2**BitsIndexed1+1)
    pf1.close()
    LogClosingFile1(pf1)
    LogNote1('Finished unpicle/unpack',LogFile1)
except Exception as e:
    LogNote1('Exception encountered in finding/reading the index file from disk: '+str(e),LogFile1)
    LogNote1('Creating the index from scratch (takes extra time but not a problem otherwise)',LogFile1)
    ## IndexConstructionBins1 = How many bins to divide data in for sorting
    ## The code below is a very rough attempt to avoid suboptimal setting.  It could be tweaked in various ways,
    ## or IndexConstructionBins1 can be set from the command line.
    qssT = 0
    if IndexConstructionBins1>1:
        IdF1 = [open(os.path.join(ScratchPath1,'Id_'+str(i)+'_'+TempFileUID1+'_reva.tmp'),mode='wb') for i in xrange(IndexConstructionBins1)]
        ValF1 = [open(os.path.join(ScratchPath1,'Val_'+str(i)+'_'+TempFileUID1+'_reva.tmp'),mode='wb') for i in xrange(IndexConstructionBins1)]
        PosF1 = [open(os.path.join(ScratchPath1,'Pos_'+str(i)+'_'+TempFileUID1+'_reva.tmp'),mode='wb') for i in xrange(IndexConstructionBins1)]
    else:
        IdTx1 = []
        ValTx1 = []
        PosTx1 = []
    if NumSeq1<256:
        IdAdtype=np.uint8
    elif NumSeq1<65536:
        IdAdtype=np.uint16
    elif NumSeq1<4294967296:
        IdAdtype=np.uint32
    else:
        IdAdtype=np.uint64 
    for Ni1 in NameARange0:
        sL1=len(SD1[Ni1]) ## number of bases for analysis: featureLen+klen1-1 for a circular feature, featureLen for a linear feature
        sL2=sL1-klen1+1 ##number of k-mers that can be extracted
        sL3=sL1-(CircD1[Ni1]*(klen1-1))  ## the actual number of bases in the feature
        s=np.zeros(sL1,dtype=dtype1)+seqnum(SD1[Ni1])
        a=np.zeros(sL1,dtype=dtype1)+seqnum(AD1[Ni1])
        n=seqNotN(SD1[Ni1])
        if sL2>=0:
            sA1=np.zeros(sL2,dtype=dtype1)
            sN1=np.ones(sL2,dtype=np.bool)
            aA1=np.zeros(sL2,dtype=dtype1)
            sP1=np.arange(1,sL2+1,dtype=dtype2)        
            aP1=np.arange(-sL3,sL2-sL3,dtype=dtype2)
            i1=0
            for kin1 in (1,2,4,8,16,32):
                if kin1>klen1:
                    break
                if klen1 & kin1:                   
                    if i1==0:
                        sA1 += s[klen1-kin1:]*e4[i1]
                        aA1 += a[klen1-kin1:]*e4[i1]
                        sN1 &= n[klen1-kin1:]
                    else:
                        sA1 += s[klen1-kin1-i1:-i1]*e4[i1]
                        aA1 += a[klen1-kin1-i1:-i1]*e4[i1]
                        sN1 &= n[klen1-kin1-i1:-i1]
                    i1+=kin1                
                s = s[:-kin1]*e4[kin1]+s[kin1:]
                a = a[:-kin1]*e4[kin1]+a[kin1:]
                n = n[:-kin1] & n[kin1:]
            NFilteredLength1 = np.sum(sN1)
            ValA1x = np.empty(NFilteredLength1*2, dtype=dtype1)
            PosA1x = np.empty(NFilteredLength1*2, dtype=dtype2)                
            ValA1x[0::2] = sA1[sN1]
            PosA1x[0::2] = sP1[sN1]
            if CircD1[Ni1]:  ## added to properly permute the junctional k-mers in the antisense orientation for proper merger with sense elements
                aA2 = np.concatenate( [aA1[-(klen1-1):] , aA1[:-(klen1-1 )]] )
                aP2 = np.concatenate( [aP1[-(klen1-1):] , aP1[:-(klen1-1 )]] )
                ValA1x[1::2] = aA2[::-1][sN1]
                PosA1x[1::2] = aP2[::-1][sN1]
            else:
                ValA1x[1::2] = aA1[::-1][sN1]
                PosA1x[1::2] = aP1[::-1][sN1]
            IdA1x=np.full(NFilteredLength1*2,Ni1,dtype=IdAdtype)
            if IndexConstructionBins1>1:
                for qssi in xrange(IndexConstructionBins1):
                    if qssi==0:
                        qssMask = np.nonzero( ValA1x<(qssi+1)*qssD )[0]
                    elif qssi<IndexConstructionBins1-1:
                        qssMask = np.nonzero( (ValA1x>=qssi*qssD) & (ValA1x<(qssi+1)*qssD) )[0]
                    else:
                        qssMask = np.nonzero( (ValA1x>=qssi*qssD) )[0]
                    IdF1[qssi].write(IdA1x[qssMask])
                    ValF1[qssi].write(ValA1x[qssMask])
                    PosF1[qssi].write(PosA1x[qssMask])
            else:
                IdTx1.append(IdA1x)
                ValTx1.append(ValA1x)
                PosTx1.append(PosA1x)         
            qssT += 2*len(sP1)
    s=None;  sA1=None; a=None; aA1=None; sP1=None; aP1=None ##Manual Garbage Collection here and below may help avoid memory overuse
    IdA1x = None; ValA1x = None; PosA1x =  None
    LogNote1('Initial Conversion completed; '+str(klen1)+'-mer positions='+str(qssT),LogFile1)
    if IndexConstructionBins1>1:
        ValA1x=None; PosA1x=None; IdA1x=None
        for qssf1 in ValF1+PosF1+IdF1:
            qssf1.close()
        IdF2 = open(os.path.join(ScratchPath1,'IdAll'+'_'+TempFileUID1+'_reva.tmp'),mode='wb')
        ValUF2 = open(os.path.join(ScratchPath1,'ValU'+'_'+TempFileUID1+'_reva.tmp'),mode='wb')
        ValLF2 = open(os.path.join(ScratchPath1,'ValL'+'_'+TempFileUID1+'_reva.tmp'),mode='wb')
        PosF2 = open(os.path.join(ScratchPath1,'PosAll'+'_'+TempFileUID1+'_reva.tmp'),mode='wb')
    else:
        IdTx1 = np.concatenate(IdTx1)
        ValTx1 = np.concatenate(ValTx1)
        PosTx1 =  np.concatenate(PosTx1)                
    TotalKmers1 = 0
    for qssi in xrange(IndexConstructionBins1):
        if IndexConstructionBins1>1:
            IdTx1 = np.frombuffer(open(os.path.join(ScratchPath1,'Id_'+str(qssi)+'_'+TempFileUID1+'_reva.tmp'),mode='rb').read(),dtype=IdAdtype)
            ValTx1 = np.frombuffer(open(os.path.join(ScratchPath1,'Val_'+str(qssi)+'_'+TempFileUID1+'_reva.tmp'),mode='rb').read(),dtype=dtype1)
            PosTx1 = np.frombuffer(open(os.path.join(ScratchPath1,'Pos_'+str(qssi)+'_'+TempFileUID1+'_reva.tmp'),mode='rb').read(),dtype=dtype2)
        Sort_i1 = np.argsort(ValTx1,kind='mergesort')
        LogNote1('Sort completed: '+str(qssi),LogFile1)
        ValTx1 = ValTx1[Sort_i1]
        LogNote1('Array ValTx1 rearranged: '+str(qssi),LogFile1)
        PosTx1 = PosTx1[Sort_i1]
        LogNote1('Array PosTx1 rearranged: '+str(qssi),LogFile1)
        IdTx1 = IdTx1[Sort_i1]
        LogNote1('Array IdTx1 rearranged: '+str(qssi),LogFile1)
        u1=np.concatenate([[True],np.not_equal(ValTx1[1:],ValTx1[:-1])])  ##true at first occurence of each unique k-mer    
        u2=np.concatenate([np.not_equal(ValTx1[1:],ValTx1[:-1]),[True]])  ##true at last occurence of each unique k-mer    
        u1 = u1  & u2 ## true only for unique k-mers
        Sort_V1dedup=ValTx1[u1]
        ValTx1=None
        Sort_P1dedup=PosTx1[u1]
        Sort_I1dedup=IdTx1[u1]
        LogNote1('Deduplication completed: '+str(qssi),LogFile1)
        Sort_I1=None
        u1=None; Sort_P1=None
        LogNote1('Span test completed: '+str(qssi),LogFile1)
        TotalKmers1 += len(Sort_V1dedup)
        if IndexConstructionBins1>1:
            IdF2.write(Sort_I1dedup)
            PosF2.write(Sort_P1dedup)
            ValUF2.write((Sort_V1dedup>>BitsSorted1).astype(dtype5))
            ValLF2.write((Sort_V1dedup & BitMask1).astype(dtype4))
            Sort_I1dedup = None
            Sort_P1dedup = None
            Sort_V1dedup = None
            LogNote1('TempFile writes completed: '+str(qssi),LogFile1)
    LogNote1('All deduplication completed; '+str(klen1)+'-mer positions='+str(TotalKmers1),LogFile1)
    LogNote1('TerminalCatenateStart',LogFile1)
    if IndexConstructionBins1>1:
        IdF2.write(np.array([NumSeq1],dtype=IdAdtype))
        PosF2.write(np.zeros(1,dtype=dtype2))
        IdF2.close(); ValUF2.close(); ValLF2.close(); PosF2.close()
    else:
        Sort_I1dedup = np.concatenate((Sort_I1dedup, np.array([NumSeq1],dtype=IdAdtype)))
        Sort_P1dedup = np.concatenate((Sort_P1dedup, np.zeros(1,dtype=dtype2)))
    LogNote1('TerminalCatenateEnd',LogFile1)
    
    LogNote1('FastIndex assembly starting',LogFile1)
    FastIndexDataType1=np.uint32  ## data type (unsigned) for any reference to positions in the sorted arrays (FastIndex) 
    if TotalKmers1>=2**32-2:
        FastIndexDataType1=np.uint64
    if IndexConstructionBins1>1:
        u3 = np.frombuffer(open(os.path.join(ScratchPath1,'ValU_'+TempFileUID1+'_reva.tmp'),mode='rb').read(),dtype=dtype5)
    else:
        u3 = (Sort_V1dedup>>BitsSorted1).astype(dtype5)
    FastIndex1=np.zeros(2**BitsIndexed1,dtype=FastIndexDataType1) ## FastIndex1[ix1] is 1+the first position in the Sort_V1dedup array where Sort_V1dedup[iy1] >> BitsSorted1>=ix1  (zero if never ==ix1)
    u4=np.concatenate(([True],np.not_equal(u3[:-1],u3[1:]),[True]))
    u5=np.nonzero(u4)[0]
    u6=(u5[1:]-u5[:-1]).astype(FastIndexDataType1)
    u5=None
    FastIndex1[u3[u4[1:]]]=u6
    u3=None;u4=None;u6=None
    FastIndex1=np.cumsum(np.concatenate([np.zeros(1,dtype=FastIndex1.dtype),FastIndex1]))
    LogNote1('FastIndex completed',LogFile1)
    if IndexConstructionBins1>1:
        Sort_V1dedup=np.frombuffer(open(os.path.join(ScratchPath1,'ValL_'+TempFileUID1+'_reva.tmp'),mode='rb').read(),dtype=dtype4)
        Sort_P1dedup=np.frombuffer(open(os.path.join(ScratchPath1,'PosAll_'+TempFileUID1+'_reva.tmp'),mode='rb').read(),dtype=dtype2)
        Sort_I1dedup=np.frombuffer(open(os.path.join(ScratchPath1,'IdAll_'+TempFileUID1+'_reva.tmp'),mode='rb').read(),dtype=IdAdtype)
        PurgeREVATempFiles(ScratchPath1)
    else:
        Sort_V1dedup = (Sort_V1dedup & BitMask1).astype(dtype4)
    V1dtype=Sort_V1dedup.dtype
    P1dtype=Sort_P1dedup.dtype
    I1dtype=Sort_I1dedup.dtype
    
    try:
        ## the following lines store the index as a disk file (making startup somewhat faster)
        LogNote1('Will try to write the index to disk',LogFile1)
        pf1=open(IndexFile1,mode='wb')
        LogOpeningFile1(pf1)
        pp1=cPickle.Pickler(pf1)
        pp1.dump(RefSeqFile1)
        pp1.dump(klen1)
        pp1.dump(len(Sort_V1dedup))
        pp1.dump(V1dtype)
        pp1.dump(P1dtype)
        pp1.dump(I1dtype)
        pp1.dump(NameA1)
        pp1.dump(FastIndex1.dtype)
        if pypy1:
            AltToFile(Sort_V1dedup,pf1)
            AltToFile(Sort_P1dedup,pf1)
            AltToFile(Sort_I1dedup,pf1)
            AltToFile(FastIndex1,pf1) 
        else:
            Sort_V1dedup.tofile(pf1)
            Sort_P1dedup.tofile(pf1)
            Sort_I1dedup.tofile(pf1)
            FastIndex1.tofile(pf1)
        pf1.close()
        LogClosingFile1(pf1)
        LogNote1('Success in writing index to disk',LogFile1)
    except Exception as e:
        LogNote1('Exception Encountered writing index to disk: '+str(e),LogFile1)
        LogNote1('****IMPORTANT: DUE TO INCOMPLETE WRITE, YOU WILL NEED TO DELETE FILE '+IndexFile1+' BEFORE RERUNNING PROGRAM',LogFile1)
        LogNote1('Unable to write indexes to disk (program may run but will need to assemble indexes next time',LogFile1)
ulen1=len(Sort_V1dedup)
NameA1.append('NotFound')
LD1.append(0)
SD1.append('')
AD1.append('')
CircD1.append(False)

CoverageSEDType1 = np.uint8
if CoverageSEMax1>255:
    CoverageSEDType1 = np.uint16
if CoverageSEMax1>65535:
    CoverageSEDType1 = np.uint32
if CoverageSEMax1>2**32-1:
    CoverageSEDType1 = np.uint64
if FullCoverage1:
    CoverageK1=np.zeros(ulen1+1,dtype=CoverageDType1)
    if ReportStarts1:
        CoverageS1=np.zeros(ulen1+1,dtype=CoverageSEDType1)
    if ReportEnds1:
        CoverageE1=np.zeros(ulen1+1,dtype=CoverageSEDType1)

def FastFindS1(aV,N1):
    lastL=int(FastIndex1[N1>>BitsSorted1])
    lastH=int(FastIndex1[(N1>>BitsSorted1)+1])
    N2=N1 & BitMask1
    if lastL==lastH:
        return ulen1
    if lastL+1==lastH:
        if aV[lastL]==N2:
            return lastL
        else:
            return ulen1
    while lastL<lastH:
        i=(lastL+lastH) >> 1
        if aV[i]<N2:
            if lastL==i:
                return ulen1
            lastL=i
        elif aV[i]>N2:
            lastH=i
        else:
            return i
    return ulen1        


filelinenum=0
Dashes1000='-'*1000

eList16=(1,2,4,8,16)
def strplus1(num1,lim1):
    if num1==lim1:
        return str(num1)+'+'
    else:
        return str(num1)

FiL1=[]
FiD1=[]
FiN1=[]
FiLR1=[]
FiLR2=[]
FiLN1=[]
def versioner1(v88):
    s88 = ['0']
    for c1 in v88:
        if c1.isdigit():
            s88[-1]+=c1
        else:
            if s88[-1]!='0':
                s88.append('0')
    return(map(int,s88))
fqdVersion1 = ['0']

def Moment1(a):
    s1 = sum(a)
    s2 = sum(a[i]*i for i in range(len(a)))
    if s1==0:
        s1=1
    return float(s2)/s1

for f0,d0,n0 in zip(FiL0,FiD0,FiN0):
    if os.path.isfile(f0) or f0.endswith('.fastq') or f0.endswith('.fastq.gz') or f0.endswith('.fasta') or f0.endswith('.fasta.gz'):
        FiL1.append(f0)
        FiD1.append(d0)
        FiN1.append(n0)
    else:
        LogNote1(f0+" looks like a non-fasta, non-fastq filename; will assume it's an NCBI SRA link and try to download",LogFile1)
        LogNote1("Preparing to download sequence read set "+f0+" from NCBI",LogFile1)
        try:
            TryFastQDump1 = subprocess.check_output([FastQDumpProgram1,'-version'])
        except:
            LogNote1("Searching for a version of fastq-dump that will run; if this fails, you may need to redownload the program and unzip the archive, also add FastQDumpProgram=<path to program> to command line",LogFile1)
            os.environ["PATH"]=os.getenv("PATH")+':./:/opt/local/bin:/opt/local/sbin:/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin:/usr/X11/bin:/Applications:~/Downloads'
            for fqd1 in os.getenv("PATH").split(':'):
                if os.path.isdir(fqd1):
                    for fqd2 in os.listdir(fqd1):
                        if fqd2.startswith('sratoolkit'):
                            fqd3 = os.path.join(fqd1,fqd2)
                            if os.path.isdir(fqd3):
                                fqd4 = os.path.join(fqd3,'bin','fastq-dump')
                                if os.path.isfile(fqd4):
                                    if versioner1(fqd2) > versioner1(fqdVersion1):
                                        fqdVersion1 = fqd2
                                        FastQDumpProgram1 = fqd4
                                        TryFastQDump1 = subprocess.check_output([FastQDumpProgram1,'-version'])
        LogNote1("Trying presumed fastq-dump program file located at "+FastQDumpProgram1,LogFile1)
        if LinesToProcess1<1:
            TryFastQDump1 = subprocess.check_output([FastQDumpProgram1,
                                                 '--fasta',
                                                 '0',
                                                 '--origfmt',
                                                 '--split-files',
                                                 '--gzip',
                                                 '--outdir',
                                                 ScratchPath1,
                                                 f0])
        else:
            TryFastQDump1 = subprocess.check_output(['fastq-dump',
                                                 '--fasta',
                                                 '0',
                                                 '--origfmt',
                                                 '--split-files',
                                                 '--gzip',
                                                 '-X',
                                                 str(LinesToProcess1),
                                                 '--outdir',
                                                 ScratchPath1,
                                                 f0])
        LogNote1("Result of "+f0+" NCBI Download " + TryFastQDump1,LogFile1)
        PresumptiveRead1FilePath1 = os.path.join(ScratchPath1,f0+'_1.fasta.gz')
        PresumptiveRead1FilePath2 = os.path.join(ScratchPath1,f0+'_2.fasta.gz')
        if os.path.isfile(PresumptiveRead1FilePath1):
            SRRReadPath1 = os.path.join(ScratchPath1,f0.lower()+'R1_'+TempFileUID1+'.fasta.gz')
            SRRReadPath2 = os.path.join(ScratchPath1,f0.lower()+'R2_'+TempFileUID1+'.fasta.gz')
            os.rename(PresumptiveRead1FilePath1, SRRReadPath1)
            SRRTempList1.append(SRRReadPath1)
            FiL1.append(SRRReadPath1)
            FiD1.append('R1')
            FiN1.append(n0)
            if os.path.isfile(PresumptiveRead1FilePath2):
                os.rename(PresumptiveRead1FilePath2, SRRReadPath2)
                SRRTempList1.append(SRRReadPath2)

## The following conversion to Python-style structured arrays may help with some versions of Python but for the moment no advantage in terms of speed on MACOS with pypy7.3 (1/14/20)
##TypeConvertD1 = {'int8':'b', 'uint8':'B', 'unicode':'u', 'int16':'h', 'uint16':'H','int32':'l','uint32':'L','int64':'q','uint64':'Q','float32':'f','float64':'d'}
##def NumpyAsPythonArray(a):
##    return array.array(TypeConvertD1[str(a.dtype)],a.tolist()) ## as a structured python list
##    ##return a.tolist() ##[alternative as a regular python list]
##Sort_V1dedup = NumpyAsPythonArray(Sort_V1dedup)

for f0,d0,n0 in zip(FiL1,FiD1,FiN1):
    if d0=='R1':
        rDom1 = 'R1'
    elif  d0=='R2':
        rDom1 = 'R2'
    else:
        rDom1='R1'
        if 'R1' in f0 and not('R2' in f0):
            rDom1='R1'
        elif 'R2' in f0 and not('R1' in f0) and not(('SRR2' in f0) and f0.count('R2')==1):
            rDom1='R2'
        elif 'R1' in f0 and 'R2' in f0:
            if f0.rfind('R1')>f0.rfind('R2'):
                rDom1='R1'
            else:
                rDom1='R2'
    if rDom1=='R1':
        f0r1=f0[:]
        f0r2=f0[::-1].replace('1R','2R',1)[::-1]
    else:
        f0r2=f0[:]
        f0r1=f0[::-1].replace('2R','1R',1)[::-1]
    if not(f0r1 in FiLR1):
        if os.path.isfile(f0r1):
            FiLR1.append(f0r1)
        else:
            FiLR1.append('')
        if (f0r1 != f0r2) and os.path.isfile(f0r2):
            FiLR2.append(f0r2)
        else:
            FiLR2.append('')
        FiLN1.append(n0)
        if FiLR1[-1]=='' and FiLR2[-1]=='':
            LogNote1("***Warning*** Input Files "+f0r1+' and '+f0r2+' not found.  Trying to proceed with any other files that are found',LogFile1)
if FiLR1==[]:
    LogNote1("NO READ INPUT FILES (e.g., fastq's) SPECIFIED- RUNNING IN DIAGNOSTIC MODE WITH A FEW TEST SEQUENCES",LogFile1)
    FiLR1.append(test1)
    FiLR2.append(test2)
    
if any(FiLR2):
    Read2Input = True
else:
    Read2Input = False
if any(FiLR1):
    Read1Input = True
else:
    Read1Input = False

TaskHeader1 = '<!--PreREVA_Task_Header: '+OutFileNameBase+'-->\n'+\
    '<!--Reference_File: '+RefSeqFile1+'-->\n'
for fin1,(fi1,fi2) in enumerate(zip(FiLR1,FiLR2)):
    if fi1:
        TaskHeader1+='<!--Read_1_File_'+str(fin1)+': '+FileInfo1(fi1)+'-->\n'
    if fi2:
        TaskHeader1+='<!--Read_2_File_'+str(fin1)+': '+FileInfo1(fi2)+'-->\n'
for fin1,GFF11 in enumerate(GFF1):
    TaskHeader1+='<!--Feature_File_'+str(fin1)+': '+FileInfo1(GFF11)+'-->\n'
for fin1,GFF11 in enumerate(BaseByBaseRegionsFiles1):
    TaskHeader1+='<!--BaseByBaseRegions_File_'+str(fin1)+': '+FileInfo1(GFF11)+'-->\n'
TaskHeader1 +=  '<!--DFAM_File: '+DFAM1+'-->\n'+\
    '<!--K-mer_Length: '+str(klen1)+'-->\n'+\
    '<!--Command_Line: '+' '.join(argv)+'-->\n'+\
    '<!--PythonVersion: '+','.join(version.splitlines())+'-->\n'+\
    '<!--PreREVA_Version: '+FileInfo1(argv[0])+'-->\n'
TaskHeader1+='\n'
AbbrevHeader1 = ''.join(TaskHeader1.splitlines()[:-1])+'<!--PreREVATableHeader-->'  ##ending with ':REVATableHeader' identifies a line as a row of table headers

NameA2 = [na1+'__'+RefAbbrev1 for na1 in NameA1]

krange1 = range(klen1)
xe4=np.array([4**i for i in range(klen1)],dtype=dtype1)[::-1]

exptlinenum1=0
klen2 = klen1//2
klen1dm2 = 2*klen1-2
if not(klen1d):
    klen1d = 2*klen1
SingleReadMode1 = 0  ## 0 means both R1 and R2, 1 for just R1, 2 for just R2. SingleReadMode1 is set for each input file
SingleReadMode0 = 0  ## 0 means both R1 and R2, 1 for just R1, 2 for just R2. SingleReadMode0 is set for the entire REVA run
FastAProject0 = False
for fnn1 in FiLR1+FiLR2:
    if fnn1.lower().endswith('fasta') or fnn1.lower().endswith('fasta.gz'):
        FastAProject0 = True

if ''.join(FiLR1):
    SingleReadMode0+=1
    if FastAProject0:
        FastF1 = OutFileNameBase+'_R1_PreFilter.fasta'
    else:
        FastF1 = OutFileNameBase+'_R1_PreFilter.fastq'
    F11 = open(FastF1, mode='w')    
if ''.join(FiLR2):
    SingleReadMode0+=2
    if FastAProject0:
        FastF2 = OutFileNameBase+'_R2_PreFilter.fasta'
    else:
        FastF2 = OutFileNameBase+'_R2_PreFilter.fastq'
    F12 = open(FastF2, mode='w')    

if FirstKTable1 or InterimKTable1:
    J0A  = [0]*(klen1d+1); J1A  = [0]*(klen1d+1); J2A  = [0]*(klen1d+1); J3A  = [0]*(klen1d+1)
    JS0A = [0]*(klen1d+1); JS1A = [0]*(klen1d+1); JS2A = [0]*(klen1d+1); JS3A = [0]*(klen1d+1)
    JA0A = [0]*(klen1d+1); JA1A = [0]*(klen1d+1); JA2A = [0]*(klen1d+1); JA3A = [0]*(klen1d+1)        
if FirstKTable1 and InterimKTable1:
    LogNote1('FirstKTable and InterimKTable options are incompatible, running InterimKTable Only.  Run without InterimKTable to obtain global KTable',LogFile1)        
    FirstKTable1 = False

SpecialKeep1 = KeepAllSense1 or KeepAllAntisense1 or KeepAllBoth1 or KeepAllEither1 or KeepAllNeither1  ## Any special instructions on keeping reads
for (f1,f2,n00) in zip(FiLR1,FiLR2,FiLN1):
    if R1Only:
        f2 = ''
    if R2Only:
        f1 = ''
    if f1:
        if f1.endswith('.gz'):
            if version.startswith('2.'):
                F1=gzip.open(f1,mode='r')
            else:
                F1=gzip.open(f1,mode='rt')
        else:
            if version.startswith('2.'):
                F1=open(f1,mode='rU')
            else:
                F1=open(f1,mode='r')                
        LogOpeningFile1(F1)
    else:
        F1 = iter(())
        SingleReadMode1 = 2
        if not(R2Only):
            LogNote1('no Read1 file for '+f2+' was found, running PreREVA in single read mode with just read file '+f1,LogFile1)
    if f2:
        if f2.endswith('.gz'):
            if version.startswith('2.'):
                F2=gzip.open(f2,mode='r')
            else:
                F2=gzip.open(f2,mode='rt')
        else:
            if version.startswith('2.'):
                F2=open(f2,mode='rU')
            else:
                F2=open(f2,mode='r')
        LogOpeningFile1(F2)
    else:
        F2 = iter(())
        SingleReadMode1 = 1
        if not(R1Only):
            LogNote1('no Read2 file for '+f1+' was found, running PreREVA in single read mode with just read file '+f1,LogFile1)
    filelinenum=0
    if f1:
        Mnemonic2=os.path.basename(f1)[::-1].split('1R',1)[-1][::-1].strip('_')  # file base for a given run
    else:
        Mnemonic2=os.path.basename(f2)[::-1].split('2R',1)[-1][::-1].strip('_')
    if '#' in n00:
        Mnemonic2 += 'ie'+n00.split('#',1)[1]
    FastAFile1 = False
    LineDensity1 = 4
    HotLine0 = 1
    if f1.lower().endswith('fasta') or  f1.lower().endswith('fasta.gz') or f2.lower().endswith('fasta') or  f2.lower().endswith('fasta.gz'):
        FastAFile1 = True
        LineDensity1 = 2
        HotLine0 = 1
    LinesToProcessA=LinesToProcess1*LineDensity1
    ReportIntervalA=LineDensity1*ReportInterval1
    FullCoverageIntervalA=LineDensity1*FullCoverageInterval1
    CoverageAccumulationIntervalA=LineDensity1*CoverageAccumulationInterval1
    PreFilters1 = 0
    PassFilters1 = 0
    BuffR1 = ['','','','']
    BuffR2 = ['','','','']
    Roton0 = 0
    if SeqIndexMode1>0:
        SeqIndexD1 = {}
        for rn1,(L0,M0) in enumerate(izip_longest(F1,F2,fillvalue='')):
            if rn1%LineDensity1 == 0:
                bcL1 = GetSeqIndexFromReadID1(L0)
                bcM1 = GetSeqIndexFromReadID1(M0)
                if not((bcL1,bcM1)) in SeqIndexD1:
                    SeqIndexD1[(bcL1,bcM1)] = 0
                SeqIndexD1[(bcL1,bcM1)] +=  1
        SeqIndexMax1 = max(SeqIndexD1.values())
        for key1 in list(SeqIndexD1):
            if not SeqIndexD1[key1]==SeqIndexMax1:
                del(SeqIndexD1[key1])
        if f1:
            F1.seek(0)
        if f2:
            F2.seek(0)
    for L0,M0 in izip_longest(F1,F2,fillvalue=''):
        if Roton0 == 0:
            GoodSeqIndex1 = True
            if SeqIndexD1:
                bcL1 = L0.strip().split(':')[-1]
                bcM1 = M0.strip().split(':')[-1]
                if not((bcL1,bcM1) in SeqIndexD1):
                    GoodSeqIndex1 = False
        BuffR1[Roton0] = L0
        BuffR2[Roton0] = M0            
        Roton0 += 1
        if Roton0==LineDensity1:
            Roton0 = 0
        filelinenum += 1
        if filelinenum % ReportIntervalA==0:
            LogNote1('finished '+Mnemonic2+
                     ' read: '+str(filelinenum//LineDensity1)+
                     ' PreFilters: '+str(PreFilters1)+
                     ' PassFilters: '+str(PassFilters1),LogFile1)
            PreFilters1 = 0
            PassFilters1 = 0
            if InterimKTable1:
                if SingleReadMode1==0:
                    LogNote1('    Average_FirstKValues_Cumulative: StartR1=' + '{0:.2f}'.format(Moment1(J0A))+
                                                 ', EndR1='   + '{0:.2f}'.format(Moment1(J1A))+
                                                 ', StartR2=' + '{0:.2f}'.format(Moment1(J2A))+
                                                 ', EndR2='   + '{0:.2f}'.format(Moment1(J3A)),LogFile1)
                elif SingleReadMode1==1:
                    LogNote1('    Average_FirstKValues_Cumulative: StartR1=' + '{0:.2f}'.format(Moment1(J0A))+
                                                 ', EndR1='   + '{0:.2f}'.format(Moment1(J1A)),LogFile1)
                elif SingleReadMode1==2:
                    LogNote1('    Average_FirstKValues_Cumulative: StartR1=' + '{0:.2f}'.format(Moment1(J2A))+
                                            ', EndR1='   + '{0:.2f}'.format(Moment1(J3A)),LogFile1)
                J0A  = [0]*(klen1d+1); J1A  = [0]*(klen1d+1); J2A  = [0]*(klen1d+1); J3A  = [0]*(klen1d+1)
                JS0A = [0]*(klen1d+1); JS1A = [0]*(klen1d+1); JS2A = [0]*(klen1d+1); JS3A = [0]*(klen1d+1)
                JA0A = [0]*(klen1d+1); JA1A = [0]*(klen1d+1); JA2A = [0]*(klen1d+1); JA3A = [0]*(klen1d+1)        
        if filelinenum == LinesToProcessA:
            LogNote1('Hit pre-specified line maximum of '+str(LinesToProcess1)+' reads.  Proceeding after '+Mnemonic2,LogFile1)
            break
        if Roton0==0 and ((GoodSeqIndex1 and SeqIndexMode1<2) or (not(GoodSeqIndex1) and SeqIndexMode1==2)):    ##only pay attention to the second line of each 4 (fastq file structure)
            PreFilters1 += 1
            exptlinenum1+=1
            L0=BuffR1[HotLine0].strip()
            M0=BuffR2[HotLine0].strip()
            if R1Buffer5>0:
                L0 = L0[R1Buffer5:]
            if R1Buffer3>0:
                L0 = L0[:-R1Buffer3]
            if R2Buffer5>0:
                M0 = M0[R2Buffer5:]
            if R2Buffer3>0:
                M0 = M0[:-R2Buffer3]
            if DeleteNs1:
                L0=L0.replace('N','')
                M0=M0.replace('N','')
            lL0=len(L0)
            lM0=len(M0)
            MinCirDiff1 = max(lL0-Tn5DupMax1,InsertionMin1)
            MinCirDiff2 = max(lM0-Tn5DupMax1,InsertionMin1)
            if lL0<klen1 and lM0<klen1:
                continue
            x0=ulen1
            x2=ulen1
            IHaveSense1 = False
            IHaveAntisense1 = False
            if lL0>=klen1:
                b0 = array.array('B',L0.translate(Tr2).encode())
                k0 = 0
                for j0 in xrange(klen1):
                    k0 = (k0 << 2) + b0[j0]
                x0 = FastFindS1(Sort_V1dedup,k0)
                if x0==ulen1:
                    for j0 in xrange(klen1,lL0):
                        k0 = ((k0 << 2) + b0[j0]) & BitMask0
                        x0 = FastFindS1(Sort_V1dedup,k0)
                        if x0<ulen1:
                            break
                if x0<ulen1:
                    k1 = 0
                    for j1 in xrange(lL0-klen1,lL0):
                        k1 = (k1 << 2) + b0[j1]
                    x1 = FastFindS1(Sort_V1dedup,k1)
                    if x1==ulen1:
                        for j1 in xrange(lL0-1,0,-1):
                            k1 = (k1>>2) + (b0[j1-klen1]<<klen1dm2)
                            x1 = FastFindS1(Sort_V1dedup,k1)
                            if x1<ulen1:
                                break
            if lM0>=klen1:
                b2 = array.array('B',M0.translate(Tr2).encode())
                k2 = 0
                j2 = 0
                for j2 in xrange(klen1):
                    k2 = (k2 << 2) + b2[j2]
                x2 = FastFindS1(Sort_V1dedup,k2)
                if x2==ulen1:
                    for j2 in xrange(klen1,lM0):
                        k2 = ((k2 << 2) + b2[j2]) & BitMask0
                        x2 = FastFindS1(Sort_V1dedup,k2)
                        if x2<ulen1:
                            break
                if x2<ulen1:
                    k3 = 0
                    for j3 in xrange(lM0-klen1,lM0):
                        k3 = (k3 << 2) + b2[j3]                   
                    x3 = FastFindS1(Sort_V1dedup,k3)
                    if x3==ulen1:
                        for j3 in xrange(lM0-1,0,-1):
                            k3 = (k3>>2) + (b2[j3-klen1]<<klen1dm2)
                            x3 = FastFindS1(Sort_V1dedup,k3)
                            if x3<ulen1:
                                break
            RClass0 = 0 ## Classification of possible Read1 junction. 0 is junk, 1 is normal, 2 is potential categorizable rearrangment, 3 is other mapped rearrangement
            RClass2 = 0 ## Classification of possible Read2 junction.
            RClass02 = 0 ## Classification of possible Read12 junction.  0 is "Not applicable (one of more reads doesn;t mal), 1 is normal, 2 is potential categorizable rearrangment, 3 is other mapped rearrangement

            if x0<ulen1:
                c0 = Sort_I1dedup[x0]; c1 = Sort_I1dedup[x1]
                p0 = Sort_P1dedup[x0]; p1 = Sort_P1dedup[x1]
                if FirstKTable1 or InterimKTable1:
                    J0A[min(klen1d,j0-klen1+1)] += 1
                    J1A[min(klen1d,lL0-j1-1)] += 1
                    if FirstKByStrand1:
                        if p0>0:
                            JS0A[min(klen1d,j0-klen1+1)] += 1
                        elif p0<0:
                            JA0A[min(klen1d,j0-klen1+1)] += 1
                        if p1>0:
                            JS1A[min(klen1d,lL0-j1-1)] += 1
                        elif p1<0:
                            JA1A[min(klen1d,lL0-j1-1)] += 1                           
                if p0<0 or p1<0: IHaveAntisense1 = True
                if p0>0 or p1>0: IHaveSense1 = True
                Es01 = p0-p1-j0+j1
                if c0==c1 and abs(Es01)<MatePairIndelLengthMax1:
                    RClass0 = 1
                elif c0==c1 and (CircleMax1>=Es01>=MinCirDiff1) or (DeletionMax1>=-Es01>=MinCirDiff1) or (MinCirDiff1>Es01>=InsertionMin1):
                    RClass0 = 2
                elif p0!=0 and p1!=0:
                    RClass0 = 3
                else:
                    RClass0 = 4
                if ReportStarts1:
                    CoverageS1[x0]= min(CoverageS1[x0]+1,CoverageSEMax1)  ## if there is a unique k-mer assign the start to the start of that k-mer; this results in a potential anomaly where the interface between repeated and unique sequence will serve as the effective start point for any read spanning the junction.  So caveat emptor!
                if ReportEnds1 and (SingleReadMode1==1):
                    CoverageE1[x1]= min(CoverageE1[x1]+1,CoverageSEMax1)
            if x2<ulen1:
                c2 = Sort_I1dedup[x2]; c3 = Sort_I1dedup[x3]
                p2 = Sort_P1dedup[x2]; p3 = Sort_P1dedup[x3]
                if FirstKTable1 or InterimKTable1:
                    J2A[min(klen1d,j2-klen1+1)] += 1
                    J3A[min(klen1d,lM0-j3-1)] += 1
                    if FirstKByStrand1:
                        if p2>0:
                            JA2A[min(klen1d,j2-klen1+1)] += 1
                        elif p2<0:
                            JS2A[min(klen1d,j2-klen1+1)] += 1
                        if p3>0:
                            JA3A[min(klen1d,lM0-j3-1)] += 1
                        elif p3<0:
                            JS3A[min(klen1d,lM0-j3-1)] += 1                           
                if p2<0 or p3<0: IHaveSense1 = True
                if p2>0 or p3>0: IHaveAntisense1 = True
                Es23 = p2-p3-j2+j3
                if c2==c3 and abs(Es23)<MatePairIndelLengthMax1:
                    RClass2 = 1
                elif c2==c3 and (CircleMax1>=Es23>=MinCirDiff1) or (DeletionMax1>=-Es23>=MinCirDiff1) or (MinCirDiff1>Es23>=InsertionMin1):
                    RClass2 = 2
                elif p0!=0 and p1!=0:
                    RClass2 = 3
                else:
                    RClass2 = 4
                if ReportEnds1:
                    CoverageE1[x2]= min(CoverageE1[x2]+1,CoverageSEMax1)
            if RClass0>0 and RClass2>0:
                EsPair1 = j1+j3-p1-p3+2*klen1
                if c1==c3 and klen1<EsPair1<=ReadSeparationMax1:
                    RClass02 = 1
                elif c1==c3 and Tn5DupMin1<=EsPair1<=Tn5DupMax1 or DeletionMax1>=EsPair1>ReadSeparationMax1:
                    RClass02 = 2
                elif p0!=0 and p2!=0:
                    RClass02 = 3
                else:
                    RClass02 = 4
            StrandKeep1 = False
            if SpecialKeep1:
                RClass0 = 0
                RClass2 = 0
                RClass02 = 0
                if (
                    (IHaveSense1 and KeepAllSense1) or 
                    (IHaveAntisense1 and KeepAllAntisense1) or 
                    (IHaveSense1 and IHaveAntisense1 and KeepAllBoth1) or 
                    ((IHaveSense1 or IHaveAntisense1) and KeepAllEither1) or 
                    (not(IHaveSense1 or IHaveAntisense1) and KeepAllNeither1)):
                        StrandKeep1 = True  ## Tells PreREVA to keep this particular read pair
            if RClass0==2 or RClass2==2 or RClass02==2 or (OtherSVCandidates1 and ( RClass0==3 or RClass2==3 or RClass02==3 )) or StrandKeep1:
                PassFilters1 += 1
                if SingleReadMode0 & 1:
                    if FastAProject0:
                        F11.write('>PreREVApick:'+Mnemonic2+'#'+str(filelinenum//LineDensity1)+"_"+BuffR1[0][1:])
                        F11.write(BuffR1[1])
                    else:
                        F11.write('>PreREVApick:'+Mnemonic2+'#'+str(filelinenum//LineDensity1)+"_"+BuffR1[0][1:])
                        for i in range(1,4):
                            F11.write(BuffR1[i])
                if SingleReadMode0 & 2:
                    if FastAProject0:
                        F12.write('>PreREVApick:'+Mnemonic2+'#'+str(filelinenum//LineDensity1)+"_"+BuffR2[0][1:])
                        F12.write(BuffR2[1])
                    else:
                        F12.write('>PreREVApick:'+Mnemonic2+'#'+str(filelinenum//LineDensity1)+"_"+BuffR2[0][1:])
                        for i in range(1,4):
                            F12.write(BuffR2[i])

    if 'file' in str(type(F1)).lower():
        F1.close()
        LogClosingFile1(F1)
    if 'file' in str(type(F2)).lower():
        F2.close()
        LogClosingFile1(F2)         

if SingleReadMode0 & 1:
    F11.close()
    LogClosingFile1(F11)
if SingleReadMode0 & 2:
    F12.close()
    LogClosingFile1(F12)


OutputHeaders1 = []
## Make generalized header list for bin-by-bin or feature-by feature output summaries
if ReportUnique1:
    OutputHeaders1.append('Reference_SingleCopyKMerCount___'+RefAbbrev1)
if ReportRepeats1:
    if ReportSense1:
        OutputHeaders1.append('Reference_LocalRepeatKMerCount___'+RefAbbrev1)
    if ReportSense1:
        OutputHeaders1.append('Reference_ChromosomalRepeatKMerCount___'+RefAbbrev1)
    if ReportSense1:
        OutputHeaders1.append('Reference_DispersedRepeatKMerCount___'+RefAbbrev1)
MultiplicityTypes1 = []
if ReportUnique1:
    MultiplicityTypes1.append('SingleCopyKmers')
if ReportRepeats1:
    MultiplicityTypes1.append('LocalRepeatKmers')
if ReportRepeats1:
    MultiplicityTypes1.append('ChromosomalRepeatKmers')
if ReportRepeats1:
    MultiplicityTypes1.append('DispersedRepeatKmers')
EventTypes1 = []
if ReportCoverage1:
    EventTypes1.append('Covering')
if ReportStarts1:
    EventTypes1.append('Start')
if ReportEnds1:
    EventTypes1.append('End')
Strands1 = []
if ReportSense1:
    Strands1.append('Sense')
if ReportAntisense1:
    Strands1.append('Antisense')
if ReportBothStrands1:
    Strands1.append('')
DuplicityTypes1 = []
if ReportPositions1:
    DuplicityTypes1.append('Positions')
if ReportReads1:
    DuplicityTypes1.append('Counts')
for ri3 in MultiplicityTypes1:
    for ri1 in EventTypes1:
        for ri4 in Strands1:
            for ri2 in DuplicityTypes1:
                CurHead1 = ''.join((ri4,ri1,ri2,'_',ri3))
                CurHead1 = CurHead1.replace("CoveringPositions","CoveredPositions")
                CurHead1 = CurHead1.replace("CoveringCounts","SummedKMerCoverage")
                if PreREVA1:
                    CurHead1 = CurHead1+"_PreREVA"
                OutputHeaders1.append(CurHead1+'__'+Mnemonic1)
if UCSCLinks1:
    OutputHeaders1.append('UCSC_Link__'+RefAbbrev1)
OutputHeaders1 = '\t'.join(OutputHeaders1)+'\t'+AbbrevHeader1
if GFF1:
    F7Buff1 = []
    Dft1 = [] # list of feature types
    Dfl1 = [] # list of stripped lines in GFF file(s)
    Dfc1 = [] # list of ordinal chromosome numbers for gffs
    Dfs1 = [] # list of GFF feature starts
    Dfe1 = [] # list of GFF feature ends
    Dfo1 = [] # orientation of GFF feature (+/-/.)  ## not used for now but could be used
    HeadWritten1 = False
    for GFF11 in GFF1:
        LogOpeningFile1(GFF11)
        if GFF11.endswith('gz'):
            if version.startswith('2.'):
                GFFile1=gzip.open(GFF11,mode='r')
            else:
                GFFile1=gzip.open(GFF11,mode='rt')
        else:
            if version.startswith('2.'):
                GFFile1=open(GFF11,mode='rU')
            else:
                GFFile1=open(GFF11,mode='r')
        for L0 in GFFile1:
            L1=L0.strip().split('\t')
            if len(L1)<5 and not(HeadWritten1):
                F7Buff1.append(L0.strip())
                continue
            if not(HeadWritten1) and not(L1[3].isdigit() and L1[4].isdigit()):
                F7Buff1.append(L0.strip()+'\t'+OutputHeaders1+'\t'+AbbrevHeader1+'\n')
                HeadWritten1 = True
                continue
            if not(HeadWritten1):
                GFFHead0 = ['Seqname','Source','Feature','Start','End','Score','Strand','Frame','Attribute']
                if len(L1)>len(GFFHead0):
                    GFFHead0.extend(['GFF_ColumnUnlabeled']*(len(L1)-len(GFFHead0)))
                F7Buff1.append('\t'.join([GFFt1+'__'+GFF11 for GFFt1 in GFFHead0[:len(L1)]])+'\t'+OutputHeaders1+'\t'+AbbrevHeader1+'\n')
                HeadWritten1 = True
            ch = L1[0]
            if ch.lower().startswith('m') and not(ch in NameD1) and not('chr'+ch in NameD1):  ## very cumbersome temporary code trying to deal with different ways of naming MtDNA
                if ('M' in NameD1) or ('chrM' in NameD1):
                    ch = 'M'
            if not ch in NameD1:
                ch = 'chr'+ch
            if (ch in NameD1) and ((L1[2] in Feature1D) or Feature1All):
                if FeatureTags1:
                    Keep1 = False
                    for ft1 in FeatureTags1:
                        if ft1 in L0.lower():
                            Keep1 = True
                            break
                    if not(Keep1):
                        continue                       
                Dfc1.append(NameD1[ch])
                Dft1.append(L1[2])
                Dfs1.append(int(L1[3]))
                Dfe1.append(int(L1[4]))
                Dfl1.append(L0.strip())
                Dfo1.append(L1[6])
        GFFile1.close()
    F7=open("FeatureSummary_"+OutFileNameBase+'.tdv',mode='w')
    LogOpeningFile1(F7)
    F7.write('<!--REVA_FeatureByFeatureCounts-->\n')
    F7.write(TaskHeader1+'\n')
    F7.write(HeaderTranspose(F7Buff1[-1])+'\n')
    F7.write('\n'.join(F7Buff1))


if BinByBin1:
    BinHead1  = '\t'.join(['Chromosome__'+RefAbbrev1,
                          'BinStart__'+RefAbbrev1,
                          'BinLength__'+RefAbbrev1])+'\t'+OutputHeaders1

    F10=open("BinByBinReadCountSummary_"+OutFileNameBase+'.tdv',mode='w')
    LogOpeningFile1(F10)
    F10.write('<!--REVA_BinByBinReadCountSummary-->\n')
    F10.write(TaskHeader1+'\n')
    F10.write(HeaderTranspose(BinHead1)+'\n')
    F10.write(BinHead1+'\n')

if ChromosomeByChromosome1:
    ChromosomeHead1  = '\t'.join(['Chromosome__'+RefAbbrev1,
                          'BinStart__'+RefAbbrev1,
                          'BinLength__'+RefAbbrev1])+'\t'+OutputHeaders1

    F11=open("ChromosomeByChromosomeReadCountSummary_"+OutFileNameBase+'.tdv',mode='w')
    LogOpeningFile1(F11)
    F11.write('<!--REVA_ChromosomeByChromosomeReadCountSummary-->\n')
    F11.write(TaskHeader1+'\n')
    F11.write(HeaderTranspose(ChromosomeHead1)+'\n')
    F11.write(ChromosomeHead1+'\n')


if GFF1 or BaseByBase1 or BinByBin1 or ChromosomeByChromosome1:
    for z in range(NumSeq1):        
        PA7s = []
        PA7a = []
        PA7k = []
        allchr = np.where(Sort_I1dedup==z)[0]
        allchrS = allchr[Sort_P1dedup[allchr]>0]
        allchrA = allchr[Sort_P1dedup[allchr]<0]
        allqS = (Sort_P1dedup[allchrS]-1) % LD1[z]
        allqA = (-Sort_P1dedup[allchrA]-klen1) % LD1[z]
        KzSU = np.zeros(LD1[z],np.uint8)
        KzSU[(allqS+klen2)%LD1[z]] = 1
        PA7k.append(KzSU)
        if ReportStarts1:
            StartzSU = np.zeros(LD1[z],CoverageS1.dtype)
            StartzAU = np.zeros(LD1[z],CoverageS1.dtype)
            StartzSU[allqS] = CoverageS1[allchrS]
            StartzAU[(allqA+klen1-1)%LD1[z]] = CoverageS1[allchrA]
            PA7s.append(StartzSU)
            PA7a.append(StartzAU)
        if ReportEnds1:
            EndzSU = np.zeros(LD1[z],CoverageS1.dtype)
            EndzAU = np.zeros(LD1[z],CoverageS1.dtype)
            EndzSU[(allqS+klen1-1)%LD1[z]] = CoverageE1[allchrS]
            EndzAU[allqA] = CoverageE1[allchrA]
            PA7s.append(EndzSU)
            PA7a.append(EndzAU)
        if GFF1:
            for t7,l7,c7,s7,e7,o7 in izip(Dft1,Dfl1,Dfc1,Dfs1,Dfe1,Dfo1):
                if c7==z:
                    if s7<=0:
                        s7=1
                    if e7>LD1[z]:
                        e7=LD1[z]
                    F7Buffer1= [l7,]
                    for pak7 in PA7k:
                        F7Buffer1.append(np.sum(pak7[s7-1:e7]))
                    for pa7s,pa7a in izip(PA7s,PA7a):
                        if o7=='-':
                            pa7s,pa7a = pa7a,pa7s
                        if ReportSense1:
                            if ReportPositions1:
                                F7Buffer1.append(np.count_nonzero(pa7s[s7-1:e7]))
                            if ReportReads1:
                                F7Buffer1.append(np.sum(pa7s[s7-1:e7]))
                        if ReportAntisense1:
                            if ReportPositions1:
                                F7Buffer1.append(np.count_nonzero(pa7a[s7-1:e7]))
                            if ReportReads1:
                                F7Buffer1.append(np.sum(pa7a[s7-1:e7]))
                        if ReportBothStrands1:
                            if ReportPositions1:
                                F7Buffer1.append(np.count_nonzero(pa7a[s7-1:e7])+np.count_nonzero(pa7s[s7-1:e7]))
                            if ReportReads1:
                                F7Buffer1.append(np.sum(pa7a[s7-1:e7])+np.sum(pa7s[s7-1:e7]))
                    if UCSCLinks1:
                        try:
                            F7Buffer1.append(UCSClink1(c7,s7-UCSCBuffer1,e7+UCSCBuffer1))
                        except:
                            pass
                    ##Wormbaselink1(n7),
                    F7.write('\t'.join(map(str,F7Buffer1))+'\n')
        if BinByBin1:
            for bi1 in xrange(0,LD1[z],SeparationGranularity1):
                end1 = min(bi1+SeparationGranularity1,LD1[z])
                F10Buffer1= [NameA1[z],]
                F10Buffer1.append(bi1+1)
                F10Buffer1.append(end1-bi1)
                for pak10 in PA7k:
                    F10Buffer1.append(np.sum(pak10[bi1:end1]))
                for pa7s,pa7a in izip(PA7s,PA7a):
                    if ReportSense1:
                        if ReportPositions1:
                            F10Buffer1.append(np.count_nonzero(pa7s[bi1:end1]))
                        if ReportReads1:
                            F10Buffer1.append(np.sum(pa7s[bi1:end1]))
                    if ReportAntisense1:
                        if ReportPositions1:
                            F10Buffer1.append(np.count_nonzero(pa7a[bi1:end1]))
                        if ReportReads1:
                            F10Buffer1.append(np.sum(pa7a[bi1:end1]))
                    if ReportBothStrands1:
                        if ReportPositions1:
                            F10Buffer1.append(np.count_nonzero(pa7a[bi1:end1])+np.count_nonzero(pa7s[bi1:end1]))
                        if ReportReads1:
                            F10Buffer1.append(np.sum(pa7a[bi1:end1])+np.sum(pa7s[bi1:end1]))
                F10.write('\t'.join(map(str,F10Buffer1))+'\n')
        if ChromosomeByChromosome1:
            F11Buffer1= [NameA1[z],]
            F11Buffer1.append(1)
            F11Buffer1.append(LD1[z])
            for pak11 in PA7k:
                F11Buffer1.append(np.sum(pak11))
            for pa7s,pa7a in izip(PA7s,PA7a):
                if ReportSense1:
                    if ReportPositions1:
                        F11Buffer1.append(np.count_nonzero(pa7s))
                    if ReportReads1:
                        F11Buffer1.append(np.sum(pa7s))
                if ReportAntisense1:
                    if ReportPositions1:
                        F11Buffer1.append(np.count_nonzero(pa7a))
                    if ReportReads1:
                        F11Buffer1.append(np.sum(pa7a))
                if ReportBothStrands1:
                    if ReportPositions1:
                        F11Buffer1.append(np.count_nonzero(pa7a)+np.count_nonzero(pa7s))
                    if ReportReads1:
                        F11Buffer1.append(np.sum(pa7a)+np.sum(pa7s))
            F11.write('\t'.join(map(str,F11Buffer1))+'\n')
    if GFF1:
        F7.close()
        LogClosingFile1(F7)
    if BinByBin1:
        F10.close()
        LogClosingFile1(F10)
    if ChromosomeByChromosome1:
        F11.close()
        LogClosingFile1(F11)

if FirstKTable1:
    J0As = float(max(1,sum(J0A)))
    J1As = float(max(1,sum(J1A)))
    J2As = float(max(1,sum(J2A)))
    J3As = float(max(1,sum(J3A)))
    if FirstKByStrand1:
        JS0As = float(max(1,sum(JS0A)))
        JS1As = float(max(1,sum(JS1A)))
        JS2As = float(max(1,sum(JS2A)))
        JS3As = float(max(1,sum(JS3A)))
        JA0As = float(max(1,sum(JA0A)))
        JA1As = float(max(1,sum(JA1A)))
        JA2As = float(max(1,sum(JA2A)))
        JA3As = float(max(1,sum(JA3A)))
    FirstKHead1 = 'FirstKOffset__'+RefAbbrev1
    if SingleReadMode1==0 or (J0As>0 and J2As>0):
        AveNote1 = '    Average_FirstKValues_Cumulative: StartR1=' + '{0:.2f}'.format(Moment1(J0A))+\
                                     ', EndR1='   + '{0:.2f}'.format(Moment1(J1A))+\
                                     ', StartR2=' + '{0:.2f}'.format(Moment1(J2A))+\
                                     ', EndR2='   + '{0:.2f}'.format(Moment1(J3A))
        FirstKHead1 += '\tR1Start_FirstKOffsetCounts__'+RefAbbrev1
        FirstKHead1 += '\tR1End_FirstKOffsetCounts__'+RefAbbrev1
        FirstKHead1 += '\tR2Start_FirstKOffsetCounts__'+RefAbbrev1
        FirstKHead1 += '\tR2End_FirstKOffsetCounts__'+RefAbbrev1
        FirstKHead1 += '\tR1Start_FirstKOffsetFraction__'+RefAbbrev1
        FirstKHead1 += '\tR1End_FirstKOffsetFraction__'+RefAbbrev1
        FirstKHead1 += '\tR2Start_FirstKOffsetFraction__'+RefAbbrev1
        FirstKHead1 += '\tR2End_FirstKOffsetFraction__'+RefAbbrev1
        if FirstKByStrand1:
            FirstKHead1 += '\tR1StartSense_FirstKOffsetCounts__'+RefAbbrev1
            FirstKHead1 += '\tR1EndSense_FirstKOffsetCounts__'+RefAbbrev1
            FirstKHead1 += '\tR2StartSense_FirstKOffsetCounts__'+RefAbbrev1
            FirstKHead1 += '\tR2EndSense_FirstKOffsetCounts__'+RefAbbrev1
            FirstKHead1 += '\tR1StartSense_FirstKOffsetFraction__'+RefAbbrev1
            FirstKHead1 += '\tR1EndSense_FirstKOffsetFraction__'+RefAbbrev1
            FirstKHead1 += '\tR2StartSense_FirstKOffsetFraction__'+RefAbbrev1
            FirstKHead1 += '\tR2EndSense_FirstKOffsetFraction__'+RefAbbrev1
            FirstKHead1 += '\tR1StartAntiSense_FirstKOffsetCounts__'+RefAbbrev1
            FirstKHead1 += '\tR1EndAntiSense_FirstKOffsetCounts__'+RefAbbrev1
            FirstKHead1 += '\tR2StartAntiSense_FirstKOffsetCounts__'+RefAbbrev1
            FirstKHead1 += '\tR2EndAntiSense_FirstKOffsetCounts__'+RefAbbrev1
            FirstKHead1 += '\tR1StartAntiSense_FirstKOffsetFraction__'+RefAbbrev1
            FirstKHead1 += '\tR1EndAntiSense_FirstKOffsetFraction__'+RefAbbrev1
            FirstKHead1 += '\tR2StartAntiSense_FirstKOffsetFraction__'+RefAbbrev1
            FirstKHead1 += '\tR2EndAntiSense_FirstKOffsetFraction__'+RefAbbrev1
    elif SingleReadMode1==1 or J0As>0:
        AveNote1 = '    Average_FirstKValues_Cumulative: StartR1=' + '{0:.2f}'.format(Moment1(J0A))+\
                                     ', EndR1='   + '{0:.2f}'.format(Moment1(J1A)) 
        FirstKHead1 += '\tR1Start_FirstKOffsetCounts__'+RefAbbrev1
        FirstKHead1 += '\tR1End_FirstKOffsetCounts__'+RefAbbrev1
        FirstKHead1 += '\tR1Start_FirstKOffsetFraction__'+RefAbbrev1
        FirstKHead1 += '\tR1End_FirstKOffsetFraction__'+RefAbbrev1
        if FirstKByStrand1:
            FirstKHead1 += '\tR1StartSense_FirstKOffsetCounts__'+RefAbbrev1
            FirstKHead1 += '\tR1EndSense_FirstKOffsetCounts__'+RefAbbrev1
            FirstKHead1 += '\tR1StartSense_FirstKOffsetFraction__'+RefAbbrev1
            FirstKHead1 += '\tR1EndSense_FirstKOffsetFraction__'+RefAbbrev1
            FirstKHead1 += '\tR1StartAntiSense_FirstKOffsetCounts__'+RefAbbrev1
            FirstKHead1 += '\tR1EndAntiSense_FirstKOffsetCounts__'+RefAbbrev1
            FirstKHead1 += '\tR1StartAntiSense_FirstKOffsetFraction__'+RefAbbrev1
            FirstKHead1 += '\tR1EndAntiSense_FirstKOffsetFraction__'+RefAbbrev1
    elif SingleReadMode1==2 or J2As>0:
        AveNote1 = '    Average_FirstKValues_Cumulative: StartR2=' + '{0:.2f}'.format(Moment1(J2A))+\
                                     ', EndR2='   + '{0:.2f}'.format(Moment1(J3A))
        FirstKHead1+= '\tR2Start_FirstKOffsetCounts__'+RefAbbrev1
        FirstKHead1+= '\tR2End_FirstKOffsetCounts__'+RefAbbrev1
        FirstKHead1 += '\tR2Start_FirstKOffsetFraction__'+RefAbbrev1
        FirstKHead1 += '\tR2End_FirstKOffsetFraction__'+RefAbbrev1
        if FirstKByStrand1:
            FirstKHead1 += '\tR2StartSense_FirstKOffsetCounts__'+RefAbbrev1
            FirstKHead1 += '\tR2EndSense_FirstKOffsetCounts__'+RefAbbrev1
            FirstKHead1 += '\tR2StartSense_FirstKOffsetFraction__'+RefAbbrev1
            FirstKHead1 += '\tR2EndSense_FirstKOffsetFraction__'+RefAbbrev1
            FirstKHead1 += '\tR2StartAntiSense_FirstKOffsetCounts__'+RefAbbrev1
            FirstKHead1 += '\tR2EndAntiSense_FirstKOffsetCounts__'+RefAbbrev1
            FirstKHead1 += '\tR2StartAntiSense_FirstKOffsetFraction__'+RefAbbrev1
            FirstKHead1 += '\tR2EndAntiSense_FirstKOffsetFraction__'+RefAbbrev1
    CapNote1 = "Note that all offset values above k-len cutoff (in this case cutoff=" +str(len(J0A)-1)+ ") have been assigned offset value at cutoff"    
    LogNote1(AveNote1,LogFile1)
    LogNote1(CapNote1,LogFile1)
    F13=open("FirstKTable_"+OutFileNameBase+'.tdv',mode='w')
    LogOpeningFile1(F13)
    F13.write('<!--REVA_FirstMatchedKMerSummary-->\n')
    F13.write(TaskHeader1+'\n')
    FirstKHead1 += '\t'+AbbrevHeader1
    F13.write(HeaderTranspose(FirstKHead1)+'\n')
    F13.write('<!--'+AveNote1+'\n')
    F13.write('<!--'+CapNote1+'\n')
    F13.write('\n')
    F13.write(FirstKHead1+'\n')
    for i in range(len(J0A)):
        F13.write(str(i)+'\t')
        if SingleReadMode1==0 or (J0As>0 and J2As>0):
            F13.write(str(J0A[i])+'\t')
            F13.write(str(J1A[i])+'\t')
            F13.write(str(J2A[i])+'\t')
            F13.write(str(J3A[i])+'\t')
            F13.write('{0:.6f}'.format(J0A[i]/J0As)+'\t')
            F13.write('{0:.6f}'.format(J1A[i]/J1As)+'\t')
            F13.write('{0:.6f}'.format(J2A[i]/J2As)+'\t')
            F13.write('{0:.6f}'.format(J3A[i]/J3As)+'\t')
            if FirstKByStrand1:
                F13.write(str(JS0A[i])+'\t')
                F13.write(str(JS1A[i])+'\t')
                F13.write(str(JS2A[i])+'\t')
                F13.write(str(JS3A[i])+'\t')
                F13.write('{0:.6f}'.format(JS0A[i]/JS0As)+'\t')
                F13.write('{0:.6f}'.format(JS1A[i]/JS1As)+'\t')
                F13.write('{0:.6f}'.format(JS2A[i]/JS2As)+'\t')
                F13.write('{0:.6f}'.format(JS3A[i]/JS3As)+'\t')
                F13.write(str(JA0A[i])+'\t')
                F13.write(str(JA1A[i])+'\t')
                F13.write(str(JA2A[i])+'\t')
                F13.write(str(JA3A[i])+'\t')
                F13.write('{0:.6f}'.format(JA0A[i]/JA0As)+'\t')
                F13.write('{0:.6f}'.format(JA1A[i]/JA1As)+'\t')
                F13.write('{0:.6f}'.format(JA2A[i]/JA2As)+'\t')
                F13.write('{0:.6f}'.format(JA3A[i]/JA3As)+'\t')

        elif SingleReadMode1==1 or J0As>0:
            F13.write(str(J0A[i])+'\t')
            F13.write(str(J1A[i])+'\t')
            F13.write('{0:.6f}'.format(J0A[i]/J0As)+'\t')
            F13.write('{0:.6f}'.format(J1A[i]/J1As)+'\t')
            if FirstKByStrand1:
                F13.write(str(JS0A[i])+'\t')
                F13.write(str(JS1A[i])+'\t')
                F13.write('{0:.6f}'.format(JS0A[i]/JS0As)+'\t')
                F13.write('{0:.6f}'.format(JS1A[i]/JS1As)+'\t')
                F13.write(str(JA0A[i])+'\t')
                F13.write(str(JA1A[i])+'\t')
                F13.write('{0:.6f}'.format(JA0A[i]/JA0As)+'\t')
                F13.write('{0:.6f}'.format(JA1A[i]/JA1As)+'\t')
        elif SingleReadMode1==2 or J2As>0:
            F13.write(str(J2A[i])+'\t')
            F13.write(str(J3A[i])+'\t')
            F13.write('{0:.6f}'.format(J2A[i]/J2As)+'\t')
            F13.write('{0:.6f}'.format(J3A[i]/J3As)+'\t')
            if FirstKByStrand1:
                F13.write(str(JS2A[i])+'\t')
                F13.write(str(JS3A[i])+'\t')
                F13.write('{0:.6f}'.format(JS2A[i]/JS2As)+'\t')
                F13.write('{0:.6f}'.format(JS3A[i]/JS3As)+'\t')
                F13.write(str(JA2A[i])+'\t')
                F13.write(str(JA3A[i])+'\t')
                F13.write('{0:.6f}'.format(JA2A[i]/JA2As)+'\t')
                F13.write('{0:.6f}'.format(JA3A[i]/JA3As)+'\t')
        F13.write('\n')
    F13.close()
    LogClosingFile1(F13)    

if 'darwin' in platform:
    try:
        Coffee_process1.kill()
    except:
        LogNote1("Couldn't kill 'caffeinate' process, you may need to manually reset your energy saving preferences (System Preferences, Energy)",LogFile1)

PurgeREVATempFiles(ScratchPath1)
LogNote1('Finished running '+' '.join(argv),LogFile1)
try:
    LogNote1(open(MyCode1,mode='r').read(),LogFile1)
except:
    pass


