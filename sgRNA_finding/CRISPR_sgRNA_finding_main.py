'''
Created on Nov 27, 2016
@author: Chen-Hao Chen
'''

import argparse,textwrap
import os,sys,re,collections,math,twobitreader,glob,time,operator
from collections import defaultdict
import subprocess as sp
from multiprocessing import Pool
import itertools
from argparse import RawTextHelpFormatter
from itertools import chain
from sgRNA_finding.CRISPR_sgRNA_finding_supp import *
from Bio.Seq import Seq
'''
def prepare_optparser():
    parser=argparse.ArgumentParser(description='Input a bed file containing search region, output all sgRNAs within.', formatter_class=RawTextHelpFormatter)

    reqgroup=parser.add_argument_group(title='Required arguments',description='')
    reqgroup.add_argument("-f","--input_file",required=True,help="Input the input file in which the sgRNA is designed, the correct format is chromosome, start, end, gene_symbol")
    reqgroup.add_argument("-b","--bit_path",required=True,help="Directory where 2bit file is stored")
    reqgroup.add_argument("-v","--bit_version",required=True,help="Bit version used to extract the sequence. Please choose from:"
                                                                       "\n 1) hg38"
                                                                       "\n 2) hg19"
                                                                       "\n 3) mm10"
                                                                       "\n 4) mm9")
    reqgroup.add_argument("-l","--spacer_length",required=True,type=int,help="Spacer length. Ex: Common SpCas9 is 20.")
    reqgroup.add_argument("-n","--PAM_name",required=True,choices=['SpCas9','SaCas9','Cpf1'],help="The name for PAM motif.")

    iogroup=parser.add_argument_group(title='Optional arguments',description='')
    iogroup.add_argument("-t","--target_path",help="Directory in which all the output stored, default is current path/input_file_sgRNA_folder",default="%s" %os.getcwd())
    iogroup.add_argument("-p","--processing",help="Number of multiprocessing! Default=100",default=100)
    iogroup.add_argument("-m","--PAM_motif",help="PAM motif. Default motifs include:"
                                            "\n 1) SpCas9: NGG"
                                            "\n 2) SaCas9: NNGRR"
                                            "\n 3) Cpf1: TTN"
                                            "\n You could modify it if you have updated motif.")

    advancegroup=parser.add_argument_group(title='Example',
    description='python2 CRISPR_sgRNA_finding.py -f demo_hg38_refseq_CDS.txt -v hg38 -l 19 -b /path/of/2bit/file/ -n SpCas9')

    args=parser.parse_args()
    return(args)

def postargs(args):
    if args.PAM_name == "SpCas9":
        args.cutting_range = "-4,-3"
        if args.PAM_motif is None:
            args.PAM_motif = "NGG"
        args.stream = "-"
    if args.PAM_name == "SaCas9":
        args.cutting_range = "-4,-3"
        if args.PAM_motif is None:
            args.PAM_motif = "NNGRRT"
        args.stream = "-"
    if args.PAM_name == "Cpf1":
        args.cutting_range = "18,23"
        if args.PAM_motif is None:
            args.PAM_motif = "TTN"
        args.stream = "+"

    if ".txt" in args.input_file:
        args.simplified_file_name=args.bit_version+"_"+args.input_file[:args.input_file.index(".txt")]
        if not os.path.isdir("%s/%s_%s_%s_sgRNA_folder/" %(args.target_path,args.input_file[:args.input_file.index(".txt")],args.spacer_length,args.PAM_name)):
            os.mkdir("%s/%s_%s_%s_sgRNA_folder/" %(args.target_path,args.input_file[:args.input_file.index(".txt")],args.spacer_length,args.PAM_name))
    else:
        args.simplified_file_name=args.bit_version+"_"+args.input_file
        if not os.path.isdir("%s/%s_%s_%s_sgRNA_folder/" %(args.target_path,args.input_file[:options.input_file],args.spacer_length,args.PAM_name)):
            os.mkdir("%s/%s_%s_%s_sgRNA_folder/" %(args.target_path,args.input_file[:args.input_file],args.spacer_length,args.PAM_name))

    if os.path.isfile("/%s/%s.2bit" %(args.bit_path,args.bit_version))==False:
        Info("No %s.2bit file in %s" %(args.bit_version,args.bit_path))
        sys.exit(1)

    args.PAM_motif=args.PAM_motif.split(",")
    args.PAM_motif=nucleotide_code(args.PAM_motif)
    args.cutting_range=[int(i) for i in args.cutting_range.split(",")]

    Info("Argument List: ")
    Info("input file = " + args.input_file)
    Info("Bit version = " + args.bit_version)
    Info("Bit path = " + args.bit_path)
    Info("PAM motif = " + ",".join(args.PAM_motif))
    Info("Selection range = " + ",".join([str(args.cutting_range[i]) for i in [0,1]]))
    Info("Number of process = " + str(args.processing))
    return(args)
'''
class sgRNA():
    def __init__(self, options):
        self.input_file = options.input_file
        self.simplified_file_name=options.simplified_file_name
        self.bit_version = options.bit_version
        self.bit_path=options.bit_path

        self.cutting_range = options.cutting_range
        self.PAM_motif=options.PAM_motif
        self.stream=options.stream
        self.PAM_name=options.PAM_name

        self.target_path="%s/%s_%s_%s_sgRNA_folder/" %(options.target_path,options.input_file[:options.input_file.index(".txt")],options.spacer_length,options.PAM_name)
        self.opts_string = "# Argument List:\n" +\
                           "# Name = %s\n" %options.input_file +\
                           "# bit version = %s\n" %options.bit_version
        self.cut_limit=1000
        self.processing=options.processing
        self.spacer_length=options.spacer_length

    def data_reformat(self):
        '''
        check the correct data format as standard input file: chromosome, start, end, gene_symbol
        '''
        data=open(self.input_file,"rt")
        os.chdir(self.target_path)
        output=open("%s/3_%s_original.txt" %(self.target_path,self.simplified_file_name),"wt")
        for line in data:
            elements=line.strip().split("\t")
            gene=elements[3]
            if gene=="n/a":
                continue
            try:
                com=gene.index(",")
                if (com+1)<len(gene):
                    pass
                else:
                    insert=elements;s = "\t".join(map(str,insert));s=s.replace(",","");output.write("%s\n" %s)
            except:
                insert=elements;s = "\t".join(map(str,insert));s=s.replace(",","");output.write("%s\n" %s)
        output.close()
        data.close()
        Info("Finished file cleaning! ")

        os.chdir(self.target_path)
        data=open("3_%s_original.txt" %self.simplified_file_name,'rt')
        output=open("3_%s_secondary.txt" %(self.simplified_file_name),'wt')
        for line in data:
            chr,start,end,symbol=line.strip().split("\t")
            start=int(start);end=int(end)
            if (end-start)<self.cut_limit:
                insert=[chr,str(start),str(end),symbol]
                output.write("%s\n" %"\t".join(insert))
            else:
                for j in range(int((end-start)/self.cut_limit)):
                    if j==0:
                        insert=[chr,str(start+(j*self.cut_limit)),str(start+(j+1)*self.cut_limit+20),symbol]
                    else:
                        insert=[chr,str(start+(j*self.cut_limit)-20),str(start+(j+1)*self.cut_limit+20),symbol]
                    output.write("%s\n" %"\t".join(insert))
                insert=[chr,str(start+((j+1)*self.cut_limit)-20),str(end),symbol]
                output.write("%s\n" %"\t".join(insert))
        data.close()
        output.close()

        data=open("3_%s_secondary.txt" %self.simplified_file_name,'rt')
        count=0;line_number=0
        for line in data:
            if line_number % 100 ==0:
                try:
                    output.close()
                except:
                    pass
                output=open("3_%s_%s.txt" %(self.simplified_file_name,count),'wt')
                count=count+1
            output.write("%s" %line)
            line_number=line_number+1

        self.counts=count
        os.remove("3_%s_original.txt" %self.simplified_file_name)
        os.remove("3_%s_secondary.txt" %self.simplified_file_name)

    def multitask(self,task):
        if task=="bed_format_trial":
            print("Totally %s iterations are going to be processed!" %self.counts)
            file_list=[];fun=eval(task)
            for count in range(self.counts):
                input_file="3_%s_%s.txt" %(self.simplified_file_name,count)
                output_file_1="4_%s_%s.txt" %(self.simplified_file_name,count)
                output_file_2="4_%s_%s_2bit.bed" %(self.simplified_file_name,count)
                file_list.append((input_file,output_file_1,output_file_2,self.cutting_range))
            if len(file_list)!=self.counts:
                Info("Wrong implement in data separation!")
                sys.exit(1)
            process_ratio=1

        elif task=="twobitreader_trial":
            file_list=[];fun=eval(task)
            for input_file in glob.glob("4_*2bit.bed"):
                file_list.append((self.bit_path,self.bit_version,self.target_path,input_file))
            process_ratio=1

        elif task=="sequence_reformat_trial":
            file_list=[];fun=eval(task)
            for count in range(self.counts):
                file_list.append(("4_%s_%s.txt" %(self.simplified_file_name,count),"5_%s_%s_2bit_sequence.txt" %(self.simplified_file_name,count)))
            process_ratio=1

        elif task=="bp_finding":
            file_list=[];fun=eval(task)
            for input_file in glob.glob("6_*"):
                file_list.append((input_file,self.cutting_range,self.PAM_motif,self.stream,self.spacer_length))
            process_ratio=1

        else:
            Info("Wrong multi-task assignment!")
            sys.exit(1)

        if self.counts>self.processing:
            for j in range(int(self.counts/self.processing)):
                Info("%s of %s has been processed!" %((j+1),(int(self.counts/self.processing))+1))
                pool=Pool(processes=self.processing*process_ratio)
                pool.map(fun,file_list[j*self.processing*process_ratio:(j+1)*self.processing*process_ratio])
                pool.close()
                pool.terminate()
                pool.join()
            pool=Pool(processes=(self.counts % self.processing)*process_ratio)
            Info("%s of %s has been processed!" %((j+2),(int(self.counts/self.processing))+1))
            pool.map(fun,file_list[(j+1)*self.processing*process_ratio:])
            pool.close()
            pool.terminate()
            pool.join()
        else:
            pool=Pool(processes=self.counts*process_ratio)
            pool.map(fun,file_list)
            pool.close()
            pool.terminate()
            pool.join()

    def file_merge_and_remove_duplicates(self):
        lines_seen=set()
        self.sequence_list=[]
        with open("8_%s_merged_sgRNA.txt" %self.simplified_file_name, 'wt') as outfile:
            for filename in glob.glob("7_*"):
                with open(filename) as readfile:
                    for line in readfile:
                        if line not in lines_seen:
                            outfile.write(line)
                            lines_seen.add(line)
                            elements=line.strip().split("\t")
                            self.sequence_list.append(elements[6].upper())
                readfile.close()
                os.remove(filename)
        outfile.close()

    def finalize(self):
        D=defaultdict(list)
        for i, item in enumerate(self.sequence_list):
            D[item].append(i)
        single_D = {k:v for k,v in D.items() if len(v)==1}

        single_index=list(chain(*(single_D.values())))
        single_index.sort(reverse=True)
        non_single_D = {k:v for k,v in D.items() if len(v)!=1}
        non_single_index=list(chain(*(non_single_D.values())))
        non_single_index.sort(reverse=True)
        output_1=open("Non_repeat_40bp_sgRNA_%s_%s" %(self.simplified_file_name,self.PAM_name),'wt')
        output_3=open("Non_repeat_%sbp_sgRNA_%s_%s" %(self.spacer_length,self.simplified_file_name,self.PAM_name),'wt')

        output_2=open("Repeat_40bp_sgRNA_%s_%s" %(self.simplified_file_name,self.PAM_name),'wt')
        output_4=open("Repeat_%sbp_sgRNA_%s_%s" %(self.spacer_length,self.simplified_file_name,self.PAM_name),'wt')

        count=0
        data=open("8_%s_merged_sgRNA.txt" %self.simplified_file_name,'rt')

        single_index_set=set(single_index)
        non_single_index_set=set(non_single_index)

        for line in data:
            elements=line.strip().split("\t")
            short_insert = [elements[i] for i in [0,3,4,6,7,8]]
            long_insert = [elements[i] for i in [0,1,2,5,7,8]]
            seq = elements[6]
            long_insert[3] = (long_insert[3]).lower()
            long_insert[3] = (long_insert[3]).replace(seq.lower(),seq)
            if count in single_index_set:
                output_3.write("%s\n" %"\t".join(short_insert))
                output_1.write("%s\n" %"\t".join(long_insert))
            elif count in non_single_index_set:
                output_4.write("%s\n" %"\t".join(short_insert))
                output_2.write("%s\n" %"\t".join(long_insert))

            count=count+1

        output_1.close()
        output_2.close()
        output_3.close()
        output_4.close()
        data.close()
        os.remove("8_%s_merged_sgRNA.txt" %self.simplified_file_name)
'''
def main():
    initial_args=prepare_optparser()
    opts=postargs(initial_args)
    g = sgRNA(opts)
    g.data_reformat()
    g.multitask("bed_format_trial");Info("Finished bed file construction!")
    g.multitask("twobitreader_trial");Info("Finished twobitreader!")
    g.multitask("sequence_reformat_trial");Info("Finished sequence reformat!")
    g.multitask("bp_finding");Info("Finished sgRNAs finding!")
    g.file_merge_and_remove_duplicates();Info("Files merged!")
    g.finalize();Info("Finished! Congratulations! :))))")

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me! ;-) Bye!\n")
        sys.exit(0)
'''
