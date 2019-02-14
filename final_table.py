# A script allowing to modify a gap file according the modification
# apply to a genome. Modification information saved in an other file

##############################IMPORTS##############################

import os
import sys
import getopt
import os.path
import numpy as np
import pandas as pd
from Bio import BiopythonWarning
from Bio import SeqIO

def usage():
    print('Usage:\n')
    print('\tpython ' +
          sys.argv[0] + ' -r <directory path> [-o <output file>]')
    print('\t\t-h or --help : display this help')
    print('\t\t-r or --directory_path : The directory of the peak_locus files')
    print('\t\t-o or --output_file : the output file')

def main(argv):

    directory_path = ""
    output_file = "final_table.output.csv"

    try:
        opts, args = getopt.getopt(sys.argv[1:], 'r:o:', [
                                   'directory_path=', 'output_file=', 'help'])
    except getopt.GetoptError:
        usage()
        sys.exit(2)

##############################OPTIONS##############################

    for opt, arg in opts:
        if opt in ('-h', '--help'):
            usage()
            sys.exit(2)
        elif opt in ('-r', '--directory_path'):
            directory_path = arg
        elif opt in ('-o', '--output_file'):
            output_file = arg
        else:
            print("Error : Bad option -> " + opt)
            usage()
            sys.exit(2)

##############################CHECK UP/SET UP##############################

    # CHECK DIRECTORY PATH
    if not os.path.exists(directory_path):
        print("Error : You have to set a directory path!\n")
        usage()
        sys.exit(2)
    else:
        if directory_path[-1] != "/":
            directory_path += "/"

##############################PRINTS##############################

    print('\n-----------------------------------------')
    print('Directory path : ' + directory_path)
    print('Output file : ' + output_file)
    print('-----------------------------------------\n')

##############################PROGRAM##############################

    # WRITE DATAFRAME ON ONE LINE
    pd.set_option('expand_frame_repr', False)

    pd.options.mode.chained_assignment = None  # default='warn'

    dict_final = {}
    #READ FILE BY FILE
    for filename in os.listdir(directory_path):
        if filename.endswith(".bdg"):
            condition = filename.split("_")[2]
            histone = filename.split("_")[4]
            locus = filename.split("_")[5].split(".")[0]
            if histone not in dict_final:
                dict_final[histone] = {}
            if locus not in dict_final[histone] :
                dict_final[histone][locus] = {}
            if condition not in dict_final[histone][locus] :
                dict_final[histone][locus][condition] = {}
            with open(directory_path+"/"+filename) as f:
                for line in f:
                    chrom = line.split("\t")[0]
                    start = line.split("\t")[1]
                    end = line.split("\t")[2]
                    likelihood = line.split("\t")[3][:-1]
                    dict_final[histone][locus][condition][chrom+"-"+start+"-"+end] = likelihood
    #print(dict_final)

    output = open(output_file, 'w')
    output.write("\t\taB\t\t\t\trB\t\t\t\tcommon\n")
    output.write("\t\tChr\tStart\tEnd\tLikelihood_ratio\tChr\tStart\tEnd\tLikelihood_ratio\tChr\tStart\tEnd\tLikelihood_ratio\n")
    for key_histone, value_histone in dict_final.items():
        output.write(key_histone+"\t")
        for key_locus, value_locus in dict_final[key_histone].items():
            output.write(key_locus+"\t")
            if len(dict_final[key_histone][key_locus]['aB']) != 0 and len(dict_final[key_histone][key_locus]['rB']) != 0 and len(dict_final[key_histone][key_locus]['common']) != 0:
                count = 0
                max_dict = max(len(dict_final[key_histone][key_locus]['aB']), len(dict_final[key_histone][key_locus]['rB']), len(dict_final[key_histone][key_locus]['common']))
                while count < max_dict:
                    if count < len(dict_final[key_histone][key_locus]['aB']):
                        output.write(list(dict_final[key_histone][key_locus]['aB'].keys())[count].split("-")[0]+"\t"+list(dict_final[key_histone][key_locus]['aB'].keys())[count].split("-")[1]+"\t"+list(dict_final[key_histone][key_locus]['aB'].keys())[count].split("-")[2]+"\t"+list(dict_final[key_histone][key_locus]['aB'].values())[count]+"\t")
                    else:
                        output.write("\t\t\t\t")
                    if count < len(dict_final[key_histone][key_locus]['rB']):
                        output.write(list(dict_final[key_histone][key_locus]['rB'].keys())[count].split("-")[0]+"\t"+list(dict_final[key_histone][key_locus]['rB'].keys())[count].split("-")[1]+"\t"+list(dict_final[key_histone][key_locus]['rB'].keys())[count].split("-")[2]+"\t"+list(dict_final[key_histone][key_locus]['rB'].values())[count]+"\t")
                    else:
                        output.write("\t\t\t\t")
                    if count < len(dict_final[key_histone][key_locus]['common']):
                        output.write(list(dict_final[key_histone][key_locus]['common'].keys())[count].split("-")[0]+"\t"+list(dict_final[key_histone][key_locus]['common'].keys())[count].split("-")[1]+"\t"+list(dict_final[key_histone][key_locus]['common'].keys())[count].split("-")[2]+"\t"+list(dict_final[key_histone][key_locus]['common'].values())[count])
                    output.write("\n\t")
                    count+=1
                    if count < max_dict:
                        output.write("\t")

            elif len(dict_final[key_histone][key_locus]['aB']) != 0 and len(dict_final[key_histone][key_locus]['rB']) != 0:
                count = 0
                max_dict = max(len(dict_final[key_histone][key_locus]['aB']), len(dict_final[key_histone][key_locus]['rB']))
                while count < max_dict:
                    if count < len(dict_final[key_histone][key_locus]['aB']):
                        output.write(list(dict_final[key_histone][key_locus]['aB'].keys())[count].split("-")[0]+"\t"+list(dict_final[key_histone][key_locus]['aB'].keys())[count].split("-")[1]+"\t"+list(dict_final[key_histone][key_locus]['aB'].keys())[count].split("-")[2]+"\t"+list(dict_final[key_histone][key_locus]['aB'].values())[count]+"\t")
                    else:
                        output.write("\t\t\t\t")
                    if count < len(dict_final[key_histone][key_locus]['rB']):
                        output.write(list(dict_final[key_histone][key_locus]['rB'].keys())[count].split("-")[0]+"\t"+list(dict_final[key_histone][key_locus]['rB'].keys())[count].split("-")[1]+"\t"+list(dict_final[key_histone][key_locus]['rB'].keys())[count].split("-")[2]+"\t"+list(dict_final[key_histone][key_locus]['rB'].values())[count])
                    output.write("\n\t")
                    count+=1
                    if count < max_dict:
                        output.write("\t")

            elif len(dict_final[key_histone][key_locus]['aB']) != 0 and len(dict_final[key_histone][key_locus]['common']) != 0:
                count = 0
                max_dict = max(len(dict_final[key_histone][key_locus]['aB']), len(dict_final[key_histone][key_locus]['common']))
                while count < max_dict:
                    if count < len(dict_final[key_histone][key_locus]['aB']):
                        output.write(list(dict_final[key_histone][key_locus]['aB'].keys())[count].split("-")[0]+"\t"+list(dict_final[key_histone][key_locus]['aB'].keys())[count].split("-")[1]+"\t"+list(dict_final[key_histone][key_locus]['aB'].keys())[count].split("-")[2]+"\t"+list(dict_final[key_histone][key_locus]['aB'].values())[count]+"\t\t\t\t\t")
                    else:
                        output.write("\t\t\t\t\t\t\t\t")
                    if count < len(dict_final[key_histone][key_locus]['common']):
                        output.write(list(dict_final[key_histone][key_locus]['common'].keys())[count].split("-")[0]+"\t"+list(dict_final[key_histone][key_locus]['common'].keys())[count].split("-")[1]+"\t"+list(dict_final[key_histone][key_locus]['common'].keys())[count].split("-")[2]+"\t"+list(dict_final[key_histone][key_locus]['common'].values())[count])
                    output.write("\n\t")
                    count+=1
                    if count < max_dict:
                        output.write("\t")

            elif len(dict_final[key_histone][key_locus]['rB']) != 0 and len(dict_final[key_histone][key_locus]['common']) != 0:
                count = 0
                max_dict = max(len(dict_final[key_histone][key_locus]['rB']), len(dict_final[key_histone][key_locus]['common']))
                while count < max_dict:
                    output.write("\t\t\t\t")
                    if count < len(dict_final[key_histone][key_locus]['rB']):
                        output.write(list(dict_final[key_histone][key_locus]['rB'].keys())[count].split("-")[0]+"\t"+list(dict_final[key_histone][key_locus]['rB'].keys())[count].split("-")[1]+"\t"+list(dict_final[key_histone][key_locus]['rB'].keys())[count].split("-")[2]+"\t"+list(dict_final[key_histone][key_locus]['rB'].values())[count]+"\t")
                    else:
                        output.write("\t\t\t\t")
                    if count < len(dict_final[key_histone][key_locus]['common']):
                        output.write(list(dict_final[key_histone][key_locus]['common'].keys())[count].split("-")[0]+"\t"+list(dict_final[key_histone][key_locus]['common'].keys())[count].split("-")[1]+"\t"+list(dict_final[key_histone][key_locus]['common'].keys())[count].split("-")[2]+"\t"+list(dict_final[key_histone][key_locus]['common'].values())[count])
                    output.write("\n\t")
                    count+=1
                    if count < max_dict:
                        output.write("\t")

            elif len(dict_final[key_histone][key_locus]['aB']) != 0:
                for key_condition_aB,value_condition_aB in dict_final[key_histone][key_locus]['aB'].items():
                    output.write(key_condition_aB.split("-")[0]+"\t"+key_condition_aB.split("-")[1]+"\t"+key_condition_aB.split("-")[2]+"\t"+value_condition_aB+"\n\t")
                    if key_condition_aB != list(dict_final[key_histone][key_locus]['aB'].keys())[-1]:
                        output.write("\t")
            elif len(dict_final[key_histone][key_locus]['rB']) != 0:

                for key_condition_rB,value_condition_rB in dict_final[key_histone][key_locus]['rB'].items():
                    output.write("\t\t\t\t"+key_condition_rB.split("-")[0]+"\t"+key_condition_rB.split("-")[1]+"\t"+key_condition_rB.split("-")[2]+"\t"+value_condition_rB+"\n\t")
                    if key_condition_rB != list(dict_final[key_histone][key_locus]['rB'].keys())[-1]:
                        output.write("\t")
            elif len(dict_final[key_histone][key_locus]['common']) != 0:

                for key_condition_common,value_condition_common in dict_final[key_histone][key_locus]['common'].items():
                    output.write("\t\t\t\t\t\t\t\t"+key_condition_common.split("-")[0]+"\t"+key_condition_common.split("-")[1]+"\t"+key_condition_common.split("-")[2]+"\t"+value_condition_common+"\n\t")
                    if key_condition_common != list(dict_final[key_histone][key_locus]['common'].keys())[-1]:
                        output.write("\t")
            else:
                output.write("\n\t")
        output.write("\n")
        #sys.exit()

if __name__ == '__main__':
    main(sys.argv[1:])
