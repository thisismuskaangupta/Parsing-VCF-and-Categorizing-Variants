#TO-DO
'''
problems to report in report
1) the pseudogenes thing
2) the one protein that is invalid??? perhaps the annotation is incorrect.
'''

#setting up the logger first and foremost
#importing 'logging', prompting the user to install it if not available.
try:
    import logging
except ModuleNotFoundError:
    print(f"Hello, the logging module was not found, please make sure it is installed, and run the script again.")
    raise SystemExit(1)

#setting up two loggers - one to print serious errors in the terminal, the other to generate the log file.
file_logger = logging.getLogger('file_logger') #getting the logger. the name is passed because multiple loggers are being used.
file_logger.setLevel(logging.INFO) #setting the level to info as we want to write information into the log file.
fhandler = logging.FileHandler('log.log') #setting up filehandler to produce an output log file called 'log.log'.
fhandler.setFormatter(logging.Formatter('%(levelname)s - %(asctime)s - {%(pathname)s:%(lineno)d} - %(message)s')) #formatting the logger messages
file_logger.addHandler(fhandler) #adding the handler to the logger

#repeating the same for a terminal logger
stream_logger = logging.getLogger('stream_logger')
stream_logger.setLevel(logging.ERROR)
shandler = logging.StreamHandler()
shandler.setFormatter(logging.Formatter('%(levelname)s - %(asctime)s - {%(pathname)s:%(lineno)d} - %(message)s'))
stream_logger.addHandler(shandler)

#importing all the necessary modules, showing an error and exiting the program if they're not installed.
try:
    from Bio.Seq import Seq
    import vcf
    import os
    import gffutils
    import math
    import Levenshtein
    import argparse
    import pandas as pd
    import seaborn as sns 
    import matplotlib.pyplot as plt
except ModuleNotFoundError as e:
    stream_logger.error(f"Hello, the following module was not found, please make sure it is installed and run the script again.\n{e}")
    raise SystemExit(1)

#setting up the command line
#parser is created
parser = argparse.ArgumentParser(prog='Biopython Assignment',description='this program takes a vcf file, a gff file, and a genome fasta file. for every variant listed in the vcf file, it uses the genome fasta files (sequence(s)), and the genome annotation file (gff), to calculate if the variant is non-coding, synonymous, or non-synonymous, and reports this information and other relevant data in an output file. it also creates a barplot of each variant type, and a log file showing progress of the program.',epilog='thanks and regards',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
#adding the arguments, all filenames are taken as positional arguments
parser.add_argument('vcfFile',help='please write the path of the vcfFile in this argument. if the vcf file is gunzipped, please make sure that the corresponding tabix index exists in the same directory with the same filename.')
parser.add_argument('gffFile',help='please write the path of the gffFile here.')
parser.add_argument('genomeFastaFile',help='please write the path of the genome fasta file here.')
#a object containing the arguments passed in the command line is created.
args = parser.parse_args()

#these are assigned to variables
gffFile = args.gffFile
fastaFile = args.genomeFastaFile
vcfFile = args.vcfFile
#the output table is also named and assigned to a variable
outputfile = 'outputTable.tsv'

#the input file names are written in the log file
file_logger.info(f"The following files were passed in the command line.\n{args.vcfFile}\n{args.gffFile}\n{args.genomeFastaFile}\n")

#checking that these files exist (except for the output file)
for file in [gffFile,fastaFile,vcfFile]:
    if not os.path.isfile(file): #raising an error and a SystemExit if they don't, as the program cannot run without them.
        stream_logger.error(f"Hello, the provided {file} does not exist in the current working directory. Please check and run the script again.")
        raise SystemExit(1)

#setting up counts
low_qual_count = 0 #this will count the variants that have a QUAL of 20 or lesser.
non_coding_count = 0 #this will count the number of non-coding variants, including those in pseudogenes.
syn_count = 0 #this will count the number of synonymous variants.
non_syn_count = 0 #this will count the number of non-synonymous variants.

#defining the functions used

#this one takes in a sequence and returns its biological complement sequence
def complement_generator(sequence:str):
    #it converts the sequence to upper case, and that manually replaces every base with its respective complement base
    sequence = sequence.upper()
    complementseq = sequence.replace('A','t')
    complementseq = complementseq.replace('C','g')
    complementseq = complementseq.replace('T','a')
    complementseq = complementseq.replace('G','c')
    complementseq = complementseq.upper()
    return complementseq

#test suite 
assert complement_generator('AATTGGCC') == 'TTAACCGG'
assert complement_generator('aattcc') == 'TTAAGG'

#this one mutates a sequence
'''
the arguments are:
1) the sequence - python str type
2) the mutation base - python str type 
3) the position where the mutation occurs in the sequence - python int type
'''
def mutation_creator(sequence:str,mutation:str,position:int): #position 1-indexed, so can directly use the POS attribute from the vcf file
    #the sequence and mutation are converted to upper case
    sequence = sequence.upper()
    mutation = mutation.upper()
    #two strings are extracted from the sequence - the one to the left of the mutation, and the one to the right. indexing is used with the assumption that the position supplied is 1-indexed.
    #these strings are then joined, first the left, then the mutation, and then the right.
    split_sequence_left = sequence[:(position-1)]
    split_sequence_right = sequence[(position):]
    return split_sequence_left+mutation+split_sequence_right

#test suite 
#also clarifies the working of the function further
assert mutation_creator('ABCDE','F',3) == 'ABFDE'
assert mutation_creator('ABCDE','F',1) == 'FBCDE'
assert mutation_creator('ABCDE','F',5) == 'ABCDF'

#this one builds and returns the sequence of interest containing the CDS of interest, and also returns the position in the sequence of interest where the variant occurs.
#the sequence of interest is a subset of the sequence of genomic DNA which encodes the open reading frame (i.e. it is coding). 
'''
the arguments are (in chronological order):
1) the transcript from the gff file which contains the variant
2) name of the database object contructed from the gff file
3) the CDS from the gff file which contains the variant
4) the vcf file record containing the variant in question
5) the fasta file containing the reference genome sequence(s)
'''
#the function returns the position of the variant (python type int) and the sequence of interest (python type str), in a list format, in that order
def seq_and_pos_calculator(gff_transcript,database_name,gff_CDS_of_interest,vcf_record,genome_fasta_filename):
    #this function iterates over the CDSes of a transcript and builds the sequence, and calculates the position of the variant.
    length = 0
    sequence_of_interest = ''
    #creating conditionals for which strand it is as the calculation is different
    if gff_transcript.strand == '+':
        CDSes = database_name.children(gff_transcript,featuretype='CDS',order_by='start') #all the CDS children of the transcript are found, and are ordered by their start position.
        for CDS in CDSes:
            if CDS != gff_CDS_of_interest: #if CDS is not the one that contains the variant, the length of the CDS is added to the length.
                length = int(CDS.end) - int(CDS.start) + 1 + length
                sequence_of_interest = sequence_of_interest + str(CDS.sequence(genome_fasta_filename,use_strand=True)) #use_strand = True ensures that the coding sequence is read with respect to the strand directionality.
            elif CDS == gff_CDS_of_interest: #if the CDS is the one that contains the variant, then the variant position is calculated.
                variant_pos = int(vcf_record.POS) - int(CDS.start) + 1 + length
                sequence_of_interest = sequence_of_interest + str(CDS.sequence(genome_fasta_filename,use_strand=True))
                #break
    
    elif gff_transcript.strand == '-': #similar logic is followed if the strand is negative.
        CDSes = database_name.children(gff_transcript,featuretype='CDS',order_by='start',reverse=True) #all the CDS children of the transcript are found, and are ordered by their start position, and then reversed because the strand is negative.
        for CDS in CDSes:
            if CDS != gff_CDS_of_interest:
                length = int(CDS.end) - int(CDS.start) + 1 + length
                sequence_of_interest = sequence_of_interest + str(CDS.sequence(genome_fasta_filename,use_strand=True))
            elif CDS == gff_CDS_of_interest:
                variant_pos = int(CDS.end) - int(vcf_record.POS) + 1 + length
                sequence_of_interest = sequence_of_interest + str(CDS.sequence(genome_fasta_filename,use_strand=True))
                #break
    '''
    initially, I chose to break the loop when the CDS of interest was found, however this was leading to errors down the line when Biopython was trying to translate the sequence, in particular, the following error:
    
    BiopythonWarning: Partial codon, len(sequence) not a multiple of three. Explicitly trim the sequence or add trailing N before translation. This may become an error in future.

    due to this, I assumed the CDSes are not multiples of three and that some of the sequences were not getting translated reliably, e.g. if the sequence is 'AAATC' and the variant is the T or a C, then we do not have enough information to calculate if the corresponding amino acid changes. for this reason, I chose to use the full sequence for every transcript here, at the cost of extra computational resources, so that the final data may be correct and reliable.
    '''
    return [variant_pos,sequence_of_interest]
#a test suite cannot be built to test this but it was tested on toy data.

#this function takes a protein sequence and tests if it is valid or not. it runs on three assumptions - that a valid protein starts with an M, ends with a *, and that there are no internal *.
#it returns the validity as a boolean value
def prot_validity_checker(sequence: str):
    sequence = sequence.upper()
    if sequence.startswith('M') and sequence.endswith('*') and sequence.count('*') == 1:
        return True
    else:
        return False
    
#test suite
assert prot_validity_checker('MLLL*') == True
assert prot_validity_checker('MMM**') == False
assert prot_validity_checker('LLL*') == False
assert prot_validity_checker('MMM*M') == False

#creating the vcfReader object
vcfReader = vcf.Reader(filename=vcfFile)

#creating the genome annotation database from the gff file
db = gffFile.replace('.gff','.db')
#if the database already exists, this is loaded. else, a new database is created.
if not os.path.isfile(db):
    db = gffutils.create_db(gffFile,dbfn=db,force=True,keep_order=True)
else:
    db = gffutils.FeatureDB(db, keep_order=True)

#opening the output file
with open(outputfile,'w') as out_file:
    #writing the header
    out_file.write(f"Chrom\tPos\tRef\tAlt\tType\tTranscript\tProtein Location\tRef AA\tAlt AA\n")
    #iterating through the vcfReader object
    for record in vcfReader:
        #filtering out low quality records
        if record.QUAL <= 20:
            #counting this
            low_qual_count = low_qual_count + 1
        else: #the records that pass the quality filter are worked with.
            #checking if the given record exists in a coding region
            if len(list(db.region(seqid=record.CHROM,start=record.POS,end=record.POS,featuretype='CDS'))) == 0:
                #if this conditional is true, that means that the record is non-coding. this is counted, and the data is added to the output file.
                non_coding_count = non_coding_count + 1
                out_file.write(f"{record.CHROM}\t{record.POS}\t{record.REF}\t{str((record.ALT))[1]}\tNon-Coding\tNA\tNA\tNA\tNA\n") #str((record.ALT))[1] is used because record.ALT returns base in '[base]' format, e.g. '[A]'. this is converted to a string and only the character at index 1, i.e. A (in this example) is written.
            else: 
                #we now know that the record is found in coding region(s)
                #these could either be functional genes, or pseudogenes
                #we iterate through every CDS that the variant occurs in
                for CDS_of_interest in db.region(seqid=record.CHROM,start=record.POS,end=record.POS,featuretype='CDS'):
                    #we find the parent transcripts for each CDS of interest
                    transcripts = db.parents(CDS_of_interest)
                    #we iterate through the transcripts
                    for transcript in transcripts:
                        #we count those that occur in pseudogenes as non-coding
                        if transcript.featuretype == 'pseudogenic_transcript':
                            non_coding_count = non_coding_count + 1
                            out_file.write(f"{record.CHROM}\t{record.POS}\t{record.REF}\t{str((record.ALT))[1]}\tNon-Coding\tNA\tNA\tNA\tNA\n") 
                        elif transcript.featuretype == 'mRNA':
                            #we filter for transcripts that encode proteins
                            #for every transcript, we find where the variant occurs and create a sequence containing all the CDSes of the respective transcript 
                            variant_pos, sequence_of_interest = seq_and_pos_calculator(transcript,db,CDS_of_interest,record,fastaFile)
                            #the function returns the position and the sequence in int and str format but to make sure they're as expected, they're converted.
                            variant_pos = int(variant_pos)
                            sequence_of_interest = str(sequence_of_interest)
                            #calculating the protein position
                            protein_pos = math.ceil(variant_pos/3) 
                            '''
                            if the sequence is AAATTT, the indexes are 1,2,3,4,5,6 (1 indexed), and the protein is Lys-Phe. 
                            the above calculation will divide this by 3 and round it up, so the indexes are 1,1,1,2,2,2 respectively, which do correspond correctly to the protein product.
                            '''
                            #now we mutate the sequence of interest
                            #if transcript is on the positive strand, then the record.ALT can directly be used.
                            if transcript.strand == '+':
                                mutated_sequence_of_interest = mutation_creator(sequence_of_interest,(str(record.ALT)[1]),variant_pos)
                            #if transcript is on the negative strand, then we need to use the complement of record.ALT as this is provided with respect to the positive strand
                            elif transcript.strand == '-':
                                ALT_base = complement_generator((str(record.ALT))[1])
                                mutated_sequence_of_interest = mutation_creator(sequence_of_interest,ALT_base,variant_pos)
                            
                            #now we have the reference DNA sequence, and the mutated DNA sequence. the next step is to translate these and compare the protein products to check if the variant is synonymous or non-synonymous.
                            #creating Seq objects.
                            try:
                                original_DNA_sequence = Seq(sequence_of_interest)
                                mutated_DNA_sequence = Seq(mutated_sequence_of_interest)
                            except:
                                file_logger.warning(f"The script could not create Seq objects for the sequences of the following transcript.\n{transcript.id}\nSkipping these.\n")
                                continue
                            #translating these
                            try:
                                original_prot_sequence = original_DNA_sequence.translate()
                                mutated_prot_sequence = mutated_DNA_sequence.translate()
                            except:
                                file_logger.warning(f"The script could not translate the sequences of the following transcript.\n{transcript.id}\nSkipping this.\n")
                                continue

                            #if everything ran correctly, then three things must be true - that the reference (original) protein sequence is valid (starts with M, ends with *, with no internal stops (*)) (we don't want to test this for the mutated sequence as the mutation may have rendered it invalid), the lengths of both the protein sequences is the same and that there is a maximum of one difference between them (the mutation).
                            #checking these

                            #checking that the protein is valid
                            try:
                                assert prot_validity_checker(original_prot_sequence) == True
                            except AssertionError:
                                file_logger.warning(f'The reference protein sequence for the transcript {transcript.id} is not valid. Following is the sequence:\n{original_prot_sequence}\nSkipping this.\n')
                                #print(original_prot_sequence)
                                continue

                            #checking that the lengths are equal
                            try:
                                assert len(original_prot_sequence) == len(mutated_prot_sequence)
                            except AssertionError:
                                file_logger.warning(f"There was an error in translating the sequences for the following transcript.\n{transcript.id}\nSkipping this.\n")
                                continue

                            #for the last assertion, we use the concept of Levenshtein distance. this is the number of edits that are required to convert one string to another. for example the Levenshtein distance between 'ABCDE' and 'ABCDF' is 1, and between 'ABCDE' and 'ABCFG', it is 2. for identical strings, it is zero.
                            try:
                                assert Levenshtein.distance(original_prot_sequence,mutated_prot_sequence) <= 1
                            except AssertionError:
                                file_logger.warning(f"There was an error in mutating the sequences for the following transcript.\n{transcript.id}\nSkipping this.\n")
                                continue

                            #if the Levenshtein distance is 0, then the variant must be non-synonymous as the mutation in the DNA sequence did not cause a change in the protein sequence, and so the two sequences are identical. this is the same as testing sequence1 == sequence2.
                            if Levenshtein.distance(original_prot_sequence,mutated_prot_sequence) == 0:
                                #this is counted and written into the output file.
                                syn_count = syn_count + 1
                                out_file.write(f"{record.CHROM}\t{record.POS}\t{record.REF}\t{(str(record.ALT))[1]}\tSynonymous\t{transcript.id}\t{protein_pos}\t{original_prot_sequence[protein_pos-1]}\tNA\n") #the index protein_pos - 1 is used because protein_pos is 1-indexed but python is zero-indexed.
                            
                            #if the Levenshtein distance is 1, then the variant must be non-synonymous.
                            elif Levenshtein.distance(original_prot_sequence,mutated_prot_sequence) == 1:
                                #this is counted and written into the output file.
                                non_syn_count = non_syn_count + 1
                                out_file.write(f"{record.CHROM}\t{record.POS}\t{record.REF}\t{(str(record.ALT))[1]}\tNon-Synonymous\t{transcript.id}\t{protein_pos}\t{original_prot_sequence[protein_pos-1]}\t{mutated_prot_sequence[protein_pos-1]}\n") 

#the main output file has been written.
                                
#now creating the barplot of the counts.
#first a pandas dataframe is created
dataframe = pd.DataFrame({'Variant Type': ['non-coding', 'synonymous', 'non-synonymous'], 'Counts': [non_coding_count, syn_count, non_syn_count]})
#print(dataframe)
#then the barplot is created
sns.barplot(data=dataframe,x='Variant Type',y='Counts')
#the plot is titled
plt.title('Types of Variants and Their Counts')
#the plot is saved
plt.savefig('TypesOfVariants.png')
#the plot object is erased.
plt.clf()

#the count of low quality records in the vcf fike is reported in the log file
file_logger.info(f"The number of records in the vcf file with QUAL <= 20 is {low_qual_count}.\n")

#if the computer reads this line, that means that the program has run successfully with no major errors. now, the list of output files can be reported in the log file
file_logger.info(f"Following are the output files, all located in the current working directory:\nThe tab-separated table is 'outputTable.tsv',\nthe bar plot is named 'TypesOfVariants.png',\nand the log file is 'log.log'.")