elliotl9.seqviewer.py ReadMe

Welcome to elliotl9.seqviewer.py! This readme.txt file provides information on setting up and using elliotl9.seqviewer.py.
Please follow the steps below to get started.

Overview
elliotl9.seqviewer.py is a SeqViewer. It is written in Python3 and designed to build a sequence viewer that helps
visualize sequence data, such as: sequence name, length, ATGC Counts, GC content, CpG islands, detected homopolymers,
Motif Search, and codon profiling. Users can upload a FASTA file and select features to change what data is displayed.
A database is embedded in the program and displays data immediately upon running the program. The program accepts two
file types: .FASTA and .Fa, and these files can contain multiple sequences.

Prerequisites
Ensure that you have the following prerequisites installed on your system:

Python (version 2.8.2): Download Python
tkinter
os

Files:
sequences.txt
sequences.sql

Additional Features
elliotl9.seqviwer.py contains two classes: ModifiedFASTAViewer and UpdatedFASTAViewer. ModifiedFASTAViewer sets up the
GUI to display sequence data and allow the user to make selections. UpdatedFASTAViewer checks the GUI variables in order
to display the correct data.

Users are able to upload multiple .FASTA .FA files without the data clearing, but users can also click the clear button
to remove data and start with a clean slate.

Additionally, elliotl9.seqviewer.py contains many functions that aid in calculating sequence data, like Motif Search,
Codon Profiling, ATGC Counts, etc. These functions are called when needed, and are called from the UpdatedFASTAViewer
class.

*When getting the ATGC counts, the function takes a sequence and gets sequence length and nucleotide count. A dictionary
of these values were returned.

*When getting the GC counts, the function takes a sequence and returns the GC Content (G counts +C counts)/ Total
Sequence Length in float value.

*A spacer function was utilized to display the sequence in non-fasta format. The printout has ruler bars to aid in
sequence viewing.

*When getting the homopolymer counts, the function takes a sequence and calculates the homopolymers by looping through
the sequence. In this case, a homopolymer was considered a minimum of 10 repeated nucleotides, with no mismatch. The
homopolymer data included the starting and ending indexes with the nucleotide.

*When getting the CpG islands, the function takes a sequence and loops through the sequence. The data was added to a
dictionary, containing the starting and ending indexes with the length. In this case, a CpG island is CG repeated a
minimum of 6 times, having a length of at least 12 nucelotides.

*The codon profiling function calculates the codon profiling by taking a sequence and utilizing a sequence of values.
If a codon match the value, the count is upped. An output was formatted using the dictionary.

*The motif search takes the sequence and user intput from the text box and loops through the function, checking
if the target is in the sequence. If it is, the starting and ending index is displayed, along with the length.

Troubleshooting
Users cannot upload files that are not .FASTA or .FA to prevent the program from crashing. Users can download .fa or
.FASTA files from a public database to use the program. Additionally, the program runs with data pre-uploading to
let users test without needing a compatible file.

Acknowledgments
The python GUI prototype was provided by Prof. Dr. Chun Liang at Miami University.

Thank you for using elliotl9.seqviewer.py!




