import os
import sys
#import re
#directory = input("Enter directory: ")
#os.chdir(directory)

#infile="regions_variants_ref.fasta"
#outfile="regions_variants_ref_oneline.fasta"

infile = sys.argv[1]
outfile = sys.argv[1]+'oneline.fasta'

def multi_2_one_fa(input_file, output_file):
    """
        This will convert multiline fasta file to one line fasta file 
    """
    print(os.getcwd())
    #output_file = re.match(r'(.*\.\d+)\.fa$',input_file).group(1) + ".fasta"
    #print(output_file)
    with open (input_file, 'r') as infile, open(output_file, 'w') as outfile:
        block = []
        #[print(line) for line in infile]
        for line in infile:
            if line.startswith('>'):
                if block:
                    outfile.write(''.join(block) + '\n')
                    block = []
                #outfile.write(line)
                foo=line.strip() + "\t"
                block.append(foo)
            else:
                block.append(line.strip())
            
        if block:
            outfile.write(''.join(block) + '\n')
            

multi_2_one_fa(infile, outfile)
# for filename in os.listdir(directory):
#     if re.match(r'.*\.\d+\.fa$',filename) is not None:
#         infile = directory + "/" + filename
#         outfile = re.match(r'(.*\.\d+)\.fa$',filename).group(1) + ".fasta"
#         multi_2_one_fa(infile, outfile)
