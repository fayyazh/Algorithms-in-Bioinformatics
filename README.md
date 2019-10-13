# Algorithms-in-Bioinformatics
The procedure to run the algorithms:

1. In command prompt, be in same directory where files are placed

2. simple python run command. i.e python filename.py fastafile (without extension)

3. It will run with default settings

4. You can enter additional values like gap penalty, scoring matrix (PAM250/BLOSUM62 without .txt), tree (upgma/wpgma)

5. These options can be entered with any order.

6. Example commmand: python NW.py fasta -g -8 -s PAM250 -t upgma
