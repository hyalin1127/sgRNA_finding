import argparse,textwrap
import os,sys,re,collections,math,twobitreader,glob,time,operator
from collections import defaultdict
import subprocess as sp
from multiprocessing import Pool
import itertools
from argparse import RawTextHelpFormatter
from itertools import chain
from Bio.Seq import Seq

def Info(infoStr):
    print("[%s] %s" %(time.strftime('%H:%M:%S'), infoStr))

def nucleotide_code(list_of_PAM_motif):
    #http://www.bioinformatics.org/sms/iupac.html
    return_list=[]
    coding_bank=dict()
    coding_bank["N"]=["A","T","C","G"]
    coding_bank["R"]=["A","G"]
    coding_bank["Y"]=["C","T"]
    coding_bank["A"]=["A"]
    coding_bank["T"]=["T"]
    coding_bank["C"]=["C"]
    coding_bank["G"]=["G"]
    for PAM_motif in list_of_PAM_motif:
        temp=[]
        for i in PAM_motif:
            temp.append(coding_bank[i])
        temp= list(itertools.product(*temp))
        for j in temp:
            return_list.append("".join(j))
    return(return_list)

def bed_format_trial(XX):
    '''
    Separate the input file into respective bed file with corresponding chromosome and GG/CC!
    '''
    (input_file,output_file_1,output_file_2,selection_range) = XX
    data=open(input_file,'rt')

    file_dict=open(output_file_1,'wt')
    file_dict_2bit=open(output_file_2,'wt')

    for line in data:
        elements=line.strip().split("\t")
        chr,start,end,gene=[elements[i] for i in range(4)]

        combined_start=int(start)-30
        combined_end=int(end)+30
        try:
            file_dict.write("%s\n" %"\t".join([chr,str(combined_start),str(combined_end),gene]))
            file_dict_2bit.write("%s\n" %"\t".join([chr,str(combined_start),str(combined_end)]))
        except:
            pass

    file_dict.close()
    file_dict_2bit.close()

    os.remove(input_file)

def twobitreader_trial(XX):
    (bit_path,bit_version,target_path,input_file) = XX
    import twobitreader
    os.chdir(bit_path)
    genome=twobitreader.TwoBitFile("%s.2bit" %bit_version)
    os.chdir(target_path)

    inputfile = open(input_file, 'rt')
    output_name=input_file[2:input_file.index(".bed")]

    orig_stdout = sys.stdout
    o_f = open("5_%s_sequence.txt" %output_name,'wt');sys.stdout = o_f

    twobitreader.twobit_reader(genome, inputfile)

    sys.stdout = orig_stdout
    inputfile.close()
    o_f.close()
    os.remove(input_file)

def sequence_reformat_trial(XX):
    (ref,sequence_data_from_2bit) = XX
    original_data=open(ref,'rt')
    if "temp_2" in locals():
        del temp_2

    output_name=sequence_data_from_2bit[2:sequence_data_from_2bit.index(".txt")]
    output=open("6_%s.txt" %output_name,'wt')

    data=open(sequence_data_from_2bit,'rt')
    i=0
    for line in data:
        elements=line.strip().split("\t")
        if elements[0].startswith(">"):
            i=i+1
            if "temp_2" in locals():
                gene_symbol= original_line[3]
                insert=[str(chr),str(start),str(end),str(temp_2),str(gene_symbol)];s = "\t".join(insert);output.write("%s\n" %s)

            temp_1=elements[0];position_1=temp_1.index(":");position_2=temp_1.index("-")
            chr=temp_1[1:position_1];start=temp_1[position_1+1:position_2];end=temp_1[position_2+1:len(temp_1)]
            temp_2=""
            original_line=original_data.readline().strip().split("\t")
        else:
            temp_2=temp_2+str(elements[0])

    if "original_line" in locals():
        gene_symbol=original_line[3]
        insert=[str(chr),str(start),str(end),str(temp_2),str(gene_symbol)];s = "\t".join(insert);output.write("%s\n" %s)
    os.remove(ref)
    os.remove(sequence_data_from_2bit)
    output.close()

def bp_finding(XX):
    (sequence_input,cutting_range,PAM_motifs,stream,spacer_length) = XX
    for PAM_motif in PAM_motifs:
        reverse_PAM_motif=Seq(PAM_motif).reverse_complement()

        data=open(sequence_input,'rt')
        output_name=sequence_input[2:sequence_input.index("_2bit")]
        output=dict()

        output["positive"]=open("7_%s_positive_%s_sgRNA.txt" %(output_name,PAM_motif),'w')
        output["negative"]=open("7_%s_negative_%s_sgRNA.txt" %(output_name,PAM_motif),'w')

        for line in data:
            elements=line.strip().split("\t")
            chr,old_start,old_end,original_sequence,gene=[elements[i] for i in range(5)]

            if isinstance(original_sequence, float)==True:
                print ("float")
                original_sequence="###"

            original_sequence=original_sequence.upper()
            for target_bp in [PAM_motif,reverse_PAM_motif]:
                if target_bp==PAM_motif:
                    strand="+"
                    if stream=="-":
                        searching_start=30-cutting_range[0]+1
                        searching_end=len(original_sequence)-30-cutting_range[1]+len(PAM_motif)-1
                    else:
                        searching_start=30-cutting_range[0]-len(PAM_motif)+1
                        searching_end=len(original_sequence)-30-cutting_range[1]-1

                    searching_sequence=original_sequence[searching_start:searching_end]
                    #searching_start=(30+selection_range[0]+1)
                elif target_bp==reverse_PAM_motif:
                    strand="-"
                    if stream=="-":
                        searching_start=30+cutting_range[1]-len(PAM_motif)+1
                        searching_end=len(original_sequence)-30+cutting_range[0]-1
                    else:
                        searching_start=30+cutting_range[1]+1
                        searching_end=len(original_sequence)-30+cutting_range[0]+len(PAM_motif)-1
                    searching_sequence=original_sequence[searching_start:searching_end]
                    #searching_start=+(30-selection_range[1]-3)

                for m in re.finditer("(?=%s)" %target_bp,searching_sequence):
                    if target_bp==PAM_motif:
                        if stream=="-":
                            relative_long_start=searching_start+m.start()-30
                            relative_short_start=searching_start+m.start()-spacer_length
                            absolute_long_start=int(old_start)+relative_long_start
                            absolute_short_start=int(old_start)+relative_short_start
                            long_segment_sequence=original_sequence[relative_long_start:(relative_long_start+40)]
                            short_segment_sequence=long_segment_sequence[40-10-spacer_length:30]
                        else:
                            relative_long_start=searching_start+m.start()+len(PAM_motif)-10
                            relative_short_start=searching_start+m.start()+len(PAM_motif)
                            absolute_long_start=int(old_start)+relative_long_start
                            absolute_short_start=int(old_start)+relative_short_start
                            long_segment_sequence=original_sequence[relative_long_start:(relative_long_start+40)]
                            short_segment_sequence=long_segment_sequence[10:10+spacer_length]
                    if target_bp==reverse_PAM_motif:
                        if stream=="-":
                            relative_long_start=searching_start+m.start()+len(PAM_motif)-10
                            relative_short_start=searching_start+m.start()+len(PAM_motif)
                            absolute_long_start=int(old_start)+relative_long_start
                            absolute_short_start=int(old_start)+relative_short_start
                            long_segment_sequence=original_sequence[relative_long_start:(relative_long_start+40)]
                            long_segment_sequence=Seq(long_segment_sequence).reverse_complement()
                            short_segment_sequence=long_segment_sequence[40-10-spacer_length:30]
                        else:
                            relative_long_start=searching_start+m.start()-30
                            relative_short_start=searching_start+m.start()-spacer_length
                            absolute_long_start=int(old_start)+relative_long_start
                            absolute_short_start=int(old_start)+relative_short_start
                            long_segment_sequence=original_sequence[relative_long_start:(relative_long_start+40)]
                            long_segment_sequence=Seq(long_segment_sequence).reverse_complement()
                            short_segment_sequence=long_segment_sequence[10:10+spacer_length]


                    insert=[str(i) for i in [chr,absolute_long_start,absolute_long_start+40,absolute_short_start,absolute_short_start+spacer_length,long_segment_sequence.upper(),short_segment_sequence.upper(),gene,strand]]
                    s = "\t".join(insert)
                    if target_bp==PAM_motif:
                        output["positive"].write("%s\n" %s)
                    if target_bp==reverse_PAM_motif:
                        output["negative"].write("%s\n" %s)
        data.close()
        for i in output.values():
            i.close()
    os.remove(sequence_input)
