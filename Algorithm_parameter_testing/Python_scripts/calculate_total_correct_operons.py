###############################################################################################################################

# script name: calculate_total_correct_operons.py

# calculates the percentage of correctly called full-length operons from multiple files after finding the number of exact matches to operons from literature
###############################################################################################################################

# Requirements:

#1. path to csv files containing operons e.g. /home/tracey/Desktop/PhD_Project/data/software/Operon_algorithms/Python_new/tb_operon_detection/output/SRR_correct_files/SRR*

# 2. path to operon list, e.g. /home/tracey/Desktop/PhD_Project/data/Python_scripts/Combined_operon_list.txt

###############################################################################################################################
import string
import glob
import os.path

message =  '\n Please note that this script works on the basis that your input file generated from COSMO, has the file name format:\n "filename_parameter1_parameter2_parameter3_parameter4.csv" \n E.g. StrainS1_D8_d3_F5_f5.csv \n'
print(message)

dir_pred_operons = input ('\nEnter the path to the file with file name, eg. /home/my_name/Desktop/my_files*: \n')
#/home/tracey/Desktop/PhD_Project/data/software/Operon_algorithms/Python_new/tb_operon_detection/output/SRR_correct_files/SRR*

operon_list_dir = input("\n Enter path to list of operons: \n")

list_files_pred_operons = glob.glob(dir_pred_operons)
outfile_strain = dir_pred_operons.split("/")[-1] # strain with its cut-off paramaters
#sample_name = outfile_strain.split("_")[0] # the actual name of sample or strain e.g. M2, M3
outfile_name = outfile_strain + "_COSMO_total_correct_operons.txt"
out_file = open(outfile_name, "w")

################################################
# add the lines from each file to a list
################################################

operon_list = []
operons = open(operon_list_dir)
for operon_line in operons:
    if operon_line.startswith("Rv"):
        operon_line = operon_line.strip().split("\t")
        operon_list.append(operon_line)

for pred_op_file in list_files_pred_operons:
    pred_operon_list = []
    op_file_name = os.path.basename(pred_op_file)
    strain = op_file_name.strip(".csv")
    print(strain)
    CDS_cut_off = strain.split("_")[-4].strip("D")
    IGR_cut_off = strain.split("_")[-3].strip("d")
    CDS_cov_diff = strain.split("_")[-2].strip("F")
    IGR_cov_diff = strain.split("_")[-1].strip("f")
    pred_operons = open(pred_op_file)
    for row in pred_operons:
        # if not row.startswith("gene_id"): # skip the header lines
        #     continue
        if not "NOT EXPRESSED" in row and row.startswith("Rv"): # skip genes that are not expressed
            row = row.strip().split(",")
            pred_operon_list.append(row)
    #print(len(pred_operon_list), strain)

# # ###########################################################
# # # Loop through the list containing the real operons and
# # # check if they were predicted by algorithm
# # ###########################################################
# #
    all_operons = []
    TP_list = []

    for data in operon_list:
        operon = data[0]
        all_operons.append(operon)

        for elem in pred_operon_list:
            if ("-") in elem[0]:
                pred_operon_start = elem[0].split(" - ")[0]#.strip('"')
                pred_operon_end = elem[0].split(" - ")[1]#.strip('"')
                if pred_operon_start.startswith("EBG"):
                    continue
                else:
                    pred_operon = pred_operon_start + "-" + pred_operon_end


 ##################################################################
# # # # if the real operon matches the predicted operon, print it
# # # ##################################################################
                if operon == pred_operon:
                    #print (operon)
                    TP_list.append(operon)
                    #len_TP_list = len(TP_list)
    len_TP_list = len(TP_list)
    Total_operons = len_TP_list/len(operon_list)
    out_line = f"{strain} {CDS_cut_off} {IGR_cut_off} {CDS_cov_diff} {IGR_cov_diff} {len_TP_list} {Total_operons}\n"
    #out_line2 = str(TP_list)
    out_file.write(out_line)

print("\n Success!! \n Your files have been printed to the directory where this script is stored.\n")
