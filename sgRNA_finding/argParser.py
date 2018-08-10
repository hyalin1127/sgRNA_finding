
import argparse,textwrap
import os,sys,re,collections,math,twobitreader,glob,time,operator
from collections import defaultdict
import subprocess as sp
from multiprocessing import Pool
import itertools
from argparse import RawTextHelpFormatter
from itertools import chain
from sgRNA_finding.CRISPR_sgRNA_finding_supp import *

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
    description='sgRNA_finding.py -f demo_hg38_refseq_CDS.txt -v hg38 -l 19 -b /path/of/2bit/file/ -n SpCas9')

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
