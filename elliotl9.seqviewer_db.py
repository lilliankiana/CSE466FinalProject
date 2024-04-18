import os
import tkinter as tk
from collections import OrderedDict
from tkinter import filedialog, Text, Scrollbar, Toplevel, Label, Button, ttk


root = tk.Tk()
root.configure(background='lightblue3')


def is_empty_file(filepath):
    return os.path.getsize(filepath) == 0


class ModifiedFASTAViewer:
    """Creation of GUI and reads displays initial database information"""
    def __init__(self, master):
        """Creates database elements"""
        # Main frame
        self.master = master
        self.master.title("FASTA Viewer")
        self.seq_header = ["Sequence", "Length", "ATGC Counts", "GC Content"]

        # Create a button to upload the FASTA file
        self.upload_button = Button(self.master, text="Upload FASTA", command=self.upload_file)
        self.upload_button.pack(pady=10)
        self.upload_button.configure(background='lightblue1')

        # Label to display the total number of sequences
        self.sequence_count_label = Label(self.master, text="")
        self.sequence_count_label.pack(pady=10)

        # Clear button to reset selections and clear the sequence display
        self.clear_button = Button(self.master, text="Clear", command=self.clear_widget)
        self.clear_button.pack(pady=10, padx=10)
        self.clear_button.configure(background='lightblue1')

        # Listbox  frame to display sequence headers with a Scrollbar
        self.listbox_frame = ttk.Frame(self.master)
        self.listbox_frame.pack(pady=10, padx=20, fill=tk.BOTH, expand=True)

        # Column by Treeview to separate data, with scrollbar
        self.tree_scroll = Scrollbar(self.listbox_frame, orient=tk.VERTICAL)
        self.tree = ttk.Treeview(self.listbox_frame, columns=self.seq_header, show="headings", yscrollcommand=self.tree_scroll.set)
        self.tree_scroll.config(command=self.tree.yview)
        self.tree_scroll.pack(side=tk.RIGHT, fill=tk.Y)
        self.tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self.tree.bind("<<TreeviewSelect>>", self.show_sequence)

        # Check boxes frame
        self.tool_bar = tk.Frame(self.master)
        self.tool_bar.pack(pady=10, padx=20)

        # Spacer Checkbox
        self.spacer = tk.IntVar()
        self.spacerCheck = tk.Checkbutton(self.tool_bar, text="Spacer", variable=self.spacer, onvalue=1, offvalue=0)
        self.spacerCheck.pack(side=tk.LEFT)

        # Homopolymer Checkbox
        self.homopolymer = tk.IntVar()
        self.homopolymerCheck = tk.Checkbutton(self.tool_bar, text="Homopolymer", variable=self.homopolymer, onvalue=1,
                                               offvalue=0)
        self.homopolymerCheck.pack(side=tk.LEFT)

        # CpG Island Checkbox
        self.cpgisland = tk.IntVar()
        self.cgpislandCheck = tk.Checkbutton(self.tool_bar, text="CpG Island", variable=self.cpgisland, onvalue=1,
                                             offvalue=0)
        self.cgpislandCheck.pack(side=tk.LEFT)

        # Codon Profile Checkbox
        self.codonprofile = tk.IntVar()
        self.codonprofilecheck = tk.Checkbutton(self.tool_bar, text="Codon Profile", variable=self.codonprofile,
                                                onvalue=1,
                                                offvalue=0)
        self.codonprofilecheck.pack(side=tk.LEFT)

        # Motif Seach Checkbox and text entry
        self.motifsearch = tk.IntVar()
        self.motifsearchcheck = tk.Checkbutton(self.tool_bar, text="Motif Search", variable=self.motifsearch, onvalue=1, offvalue=0)
        self.motifsearchcheck.pack(side=tk.LEFT)
        self.motifentry = tk.Entry(self.tool_bar, width=20)
        self.motifentry.focus_set()
        self.motifentry.pack(side=tk.LEFT)

        # Text box to display the sequence
        self.display_frame = tk.Frame(self.master)
        self.display_frame.pack(pady=10, padx=20, fill=tk.BOTH, expand=True)

        # Scrollbar to scroll through sequence output
        self.display_scrollbar = Scrollbar(self.display_frame, orient=tk.VERTICAL)
        self.sequence_display = Text(self.display_frame, yscrollcommand=self.display_scrollbar.set, height=20, width=200)
        self.display_scrollbar.config(command=self.sequence_display.yview)
        self.display_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        self.sequence_display.pack(pady=10, padx=20)

        # Dictionary to store sequences from the FASTA file
        self.sequences = OrderedDict()

        # Uploads tab-deliminated text file containing sequence database
        with open('sequences.txt', 'r') as sequence_file:
            # bypass header line
            next(sequence_file)
            sequence_data = sequence_file.readlines()

        for col in self.seq_header:
            self.tree.heading(col, text=col.title())
            self.tree.column(column=col, width=100, minwidth=0)
        sequence_list = []
        for line in sequence_data:
            line = line.strip()
            sequence = ""
            header = ""
            index = 0
            seqindex = 0
            for character in line:
                if character.isspace() and len(header) > 0 and seqindex == 0:
                    seqindex = index + 1
                elif character.isspace() and seqindex == 0:
                    header = line[0:index]
                elif character.isspace():
                    sequence = line[seqindex:index].strip()
                index += 1

            self.sequences[header] = sequence
            atgc_count = getATGC(sequence)
            gc_content = getGC(sequence)
            sequence_list.append((header, len(sequence), atgc_count, gc_content))

        for item in sequence_list:
            self.tree.insert('', 'end', values=item)
            # Update the sequence count label


    def upload_file(self):
        """Method to upload the fasta file and display subsequent data"""
        # Open a file dialog to select the FASTA file
        fasta_file = filedialog.askopenfilename(title="Select a FASTA file",
                                                filetypes=[("FASTA files", "*.fa*"), ("All files", "*.*")])
        print("================", fasta_file)

        # Display a loading message for large files
        file_size = os.path.getsize(fasta_file) / (1024 * 1024)  # File size in MB

        if file_size > 10:  # Threshold set to 10MB, can be adjusted
            loading_window = Toplevel(self.master)
            loading_label = Label(loading_window, text="Processing large file. Please wait...")
            loading_label.pack(pady=20, padx=20)
            self.master.update()
            loading_window.destroy()
        # Check if the selected file is empty
        if is_empty_file(fasta_file):
            tk.messagebox.showerror("Error", "The selected file is empty!")
            return
        if not fasta_file:
            return

        with open(fasta_file, 'r') as file:
            content = file.readlines()

        for col in self.seq_header:
            self.tree.heading(col, text=col.title())
            self.tree.column(column=col, width=100, minwidth=0)
        header = None
        sequence = ""
        sequence_list = []
        for line in content:
            line = line.strip()
            if line.startswith('>'):
                if header:
                    self.sequences[header] = sequence
                    sequence = sequence.upper()
                    atgc_count = getATGC(sequence)
                    gc_content = getGC(sequence)
                    sequence_list.append((header, len(sequence), atgc_count, gc_content))
                header = line
                sequence = ""
            else:
                sequence += line.upper()
        if header:
            self.sequences[header] = sequence
            sequence = sequence.upper()

        for item in sequence_list:
            self.tree.insert('', 'end', values=item)
            # Update the sequence count label

        # Check if the FASTA file contains any sequences
        if not self.sequences:
            tk.messagebox.showerror("Error", "The selected FASTA file contains no sequences!")
            return
        self.sequence_count_label.config(text=f"Total Sequences: {len(self.sequences)}")


def getATGC(sequence):
    """Takes a sequence and gets sequence length and nucleotide count. Returns a dictionary"""
    a_count = sequence.count("A")
    t_count = sequence.count("T")
    g_count = sequence.count("G")
    c_count = sequence.count("C")
    n_count = sequence.count("N")

    nucleotide_dict = {"A": a_count, "T": t_count, "G": g_count, "C": c_count, "N": n_count}
    return nucleotide_dict


def getGC(sequence):
    """Takes a sequence and returns the GC Content
    (G counts +C counts)/ Total Sequence Length in float value."""
    # keys of the detected CG nums
    c_count = sequence.count("C")
    g_count = sequence.count("G")

    gc_total = g_count + c_count

    gc_content = gc_total / float(len(sequence))

    return gc_content


# Modifying the ModifiedFASTAViewer class to display long sequences in a new window
def addspacerfun(sequence_display, sequence):
    """Displays the sequence in non-fasta format(ruler format)."""
    spacer_message = "\nNon-Fasta Format\n"
    # Sets up ruler
    spacer_message += "              1          2          3          4          5\n"
    spacer_message += "Line 123456890 1234567890 1234567890 1234567890 1234567890\n"
    count = 1
    # prints the sequence fragments
    for i in range(0, len(sequence), 50):
        line = sequence[i: i + 50]
        spacer_message += "  "
        spacer_message += (str(count) + " ")
        for fragment in range(0, len(line), 10):
            spacer_message += (line[fragment: fragment + 10] + ' ')

        count = count + 1
        spacer_message += "\n"

    sequence_display.insert(tk.END, spacer_message)


def addhomopolymerfun(sequence_display, sequence):
    """Calculates the homopolymers and displays in sequence text box"""
    homoList = []
    count = 1
    homopolymer = ""
    polymerCount = 0
    for i in range(1, len(sequence) - 1):
        if sequence[i] == sequence[i - 1]:
            polymerCount = polymerCount + 1
            homopolymer = sequence[i]
        elif polymerCount >= 10:
            endBase = i
            startBase = i - polymerCount
            format = str(startBase) + "-" + str(endBase) + "_" + str(endBase - startBase + 1) + "_" + homopolymer
            homoList.append(format)
            count = count + 1
            homopolymer = ""
            polymerCount = 0
        else:
            homopolymer = ""
            polymerCount = 0

    sequence_display.insert(tk.END, "\nHomopolymer Search\n")
    sequence_display.insert(tk.END, homoList)


def addcpgfun(sequence_display, sequence):
    """Calculates the CpG islands and displays in the sequence_display"""
    cpg_dict = {}
    cpgIsland = ""
    count = 1
    for i in range(0, len(sequence) - 1, 2):
        if sequence[i: i + 2] == "CG":
            cpgIsland = cpgIsland + sequence[i: i + 2]
        else:
            if len(cpgIsland) >= 12:
                format_output = str(i - len(cpgIsland)) + "-" + str(i) + "_" + str(len(cpgIsland))
                cpg_dict[count] = format_output
                count = count + 1
            cpgIsland = ""

    sequence_display.insert(tk.END, "\nCpG Islands\n")
    sequence_display.insert(tk.END, cpg_dict)


def addcodonprofilefun(sequence_display, sequence):
    """Calculates the Codon Profiling and displays in the sequence_display"""
    codon_dict = {}
    index = range(64)

    # List of dictionary values
    values = ["TTT", "TCT", "TAT", "TGT", "TTC", "TCC", "TAC", "TGC",
              "TTA", "TCA", "TAA", "TGA", "TTG", "TCG", "TAG", "TGG",
              "CTT", "CCT", "CAT", "CGT", "CTC", "CCC", "CAC", "CGC",
              "CTA", "CCA", "CAA", "CGA", "CTG", "CCG", "CAG", "CGG",
              "ATT", "ACT", "AAT", "AGT", "ATC", "ACC", "AAC", "AGC",
              "ATA", "ACA", "AAA", "AGA", "ATG", "ACG", "AAG", "AGG",
              "GTT", "GCT", "GAT", "GGT", "GTC", "GCC", "GAC", "GGC",
              "GTA", "GCA", "GAA", "GGA", "GTG", "GCG", "GAG", "GGG"]

    # Creates all dictionary values and sets them to 0
    for i in index:
        for x in values:
            codon_dict[i] = x
            codon_dict[x] = 0

    # Looping through sequence to count codons
    for i in range(len(sequence)):
        codon = sequence[i:i + 3]
        if len(codon) == 3:
            previousvalue = codon_dict[codon]
            codon_dict[codon] = previousvalue + 1

    codon_message = "\nCodon Profile\n"

    codon_message += "       -------------2nd-------------       \n"
    codon_message += "1st    T       C       A       G       3rd \n"

    # Section T
    codon_message += ("T      TTT=" + str(codon_dict["TTT"]) + "   TCT=" + str(codon_dict["TCT"]) + "   TAT=" +
                      str(codon_dict["TAT"]) + "   TGT=" + str(codon_dict["TGT"]) + "   T   \n")
    codon_message += ("       TTC=" + str(codon_dict["TTC"]) + "   TCC=" + str(codon_dict["TCC"]) + "   TAC=" +
                      str(codon_dict["TAC"]) + "   TGC=" + str(codon_dict["TGC"]) + "   C   \n")
    codon_message += ("       TTA=" + str(codon_dict["TTA"]) + "   TCA=" + str(codon_dict["TCA"]) + "   TAA=" +
                      str(codon_dict["TAA"]) + "   TGA=" + str(codon_dict["TGA"]) + "   A   \n")
    codon_message += ("       TTG=" + str(codon_dict["TTG"]) + "   TCG=" + str(codon_dict["TCG"]) + "   TAG=" +
                      str(codon_dict["TAG"]) + "   TGG=" + str(codon_dict["TGG"]) + "   G   \n")

    # Section C
    codon_message += "\n"
    codon_message += ("C      CTT=" + str(codon_dict["CTT"]) + "   CCT=" + str(codon_dict["CCT"]) + "   CAT=" +
                      str(codon_dict["CAT"]) + "   CGT=" + str(codon_dict["CGT"]) + "   T   \n")
    codon_message += ("       CTC=" + str(codon_dict["CTC"]) + "   CCC=" + str(codon_dict["CCC"]) + "   CAC=" +
                      str(codon_dict["CAC"]) + "   CGC=" + str(codon_dict["CGC"]) + "   C   \n")
    codon_message += ("       CTA=" + str(codon_dict["CTA"]) + "   CCA=" + str(codon_dict["CCA"]) + "   CAA=" +
                      str(codon_dict["CAA"]) + "   CGA=" + str(codon_dict["CGA"]) + "   A   \n")
    codon_message += ("       CTG=" + str(codon_dict["CTG"]) + "   CCG=" + str(codon_dict["CCG"]) + "   CAG=" +
                      str(codon_dict["CAG"]) + "   CGG=" + str(codon_dict["CGG"]) + "   G   \n")

    # Section A
    codon_message += "\n"
    codon_message += ("A      ATT=" + str(codon_dict["ATT"]) + "   ACT=" + str(codon_dict["ACT"]) + "   AAT=" +
                      str(codon_dict["AAT"]) + "   AGT=" + str(codon_dict["AGT"]) + "   T   \n")
    codon_message += ("       ATC=" + str(codon_dict["ATC"]) + "   ACC=" + str(codon_dict["ACC"]) + "   AAC=" +
                      str(codon_dict["AAC"]) + "   AGC=" + str(codon_dict["AGC"]) + "   C   \n")
    codon_message += ("       ATA=" + str(codon_dict["ATA"]) + "   ACA=" + str(codon_dict["ACA"]) + "   AAA=" +
                      str(codon_dict["AAA"]) + "   AGA=" + str(codon_dict["AGA"]) + "   A   \n")
    codon_message += ("       ATG=" + str(codon_dict["ATG"]) + "   ACG=" + str(codon_dict["ACG"]) + "   AAG=" +
                      str(codon_dict["AAG"]) + "   AGG=" + str(codon_dict["AGG"]) + "   G   \n")

    # Section G
    codon_message += "\n"
    codon_message += ("G      GTT=" + str(codon_dict["GTT"]) + "   GCT=" + str(codon_dict["GCT"]) + "   GAT=" +
                      str(codon_dict["GAT"]) + "   GGT=" + str(codon_dict["GGT"]) + "   T   \n")
    codon_message += ("       GTC=" + str(codon_dict["GTC"]) + "   GCC=" + str(codon_dict["GCC"]) + "   GAC=" +
                      str(codon_dict["GAC"]) + "   GGC=" + str(codon_dict["GGC"]) + "   C   \n")
    codon_message += ("       GTA=" + str(codon_dict["GTA"]) + "   GCA=" + str(codon_dict["GCA"]) + "   GAA=" +
                      str(codon_dict["GAA"]) + "   GGA=" + str(codon_dict["GGA"]) + "   A   \n")
    codon_message += ("       GTG=" + str(codon_dict["GTG"]) + "   GCG=" + str(codon_dict["GCG"]) + "   GAG=" +
                      str(codon_dict["GAG"]) + "   GGG=" + str(codon_dict["GGG"]) + "   G    \n")
    codon_message += "\n"

    sequence_display.insert(tk.END, codon_message)


def addMotifSearch(sequence_display, sequence, target):
    """Retrieves the user;s entry from the text box Calculates the motif search and displays in the sequence_display"""
    if len(target) == 0:
        sequence_display.insert(tk.END, "\nNo target entered for Motif Search!!!")
        return
    motifList = []
    count = 1
    length = len(target)
    for i in range(0, len(sequence) - 1, length):
        if (sequence[i: i + length] == target):
            format = str(i) + "-" + str(i + length - 1) + "_" + str(length)
            motifList.append(format)
            count = count + 1

    sequence_display.insert(tk.END, "\nMotif Search:\n")
    sequence_display.insert(tk.END, motifList)



class UpdatedFASTAViewer(ModifiedFASTAViewer):
    """Updates the GUI components based on user selections"""
    def show_sequence(self, sequence):
        """Displays the selected sequence data in the sequence_display"""
        selected = self.tree.focus()
        selectedvalue = self.tree.item(selected)
        if not selected:
            return

        # Removes previous values from the display
        self.sequence_display.delete(1.0, tk.END)
        header = selectedvalue["values"][0]
        sequence = self.sequences[header]

        # Checks the checkbox values
        addspacer = self.spacer.get()
        addhomopolymer = self.homopolymer.get()
        addcpg = self.cpgisland.get()
        addcodonprofile = self.codonprofile.get()
        addmotifsearch = self.motifsearch.get()

        # If sequence length exceeds 5000, open it in a new window
        if len(sequence) > 5000:
            self.show_in_new_window(sequence)
        else:
            if addspacer == 1:
                self.sequence_display.delete(1.0, tk.END)
                addspacerfun(self.sequence_display, sequence)
            else:
                self.sequence_display.insert(tk.END, "FASTA Format\n")
                self.sequence_display.insert(tk.END, sequence + '\n')

            if addhomopolymer == 1:
                addhomopolymerfun(self.sequence_display, sequence)
            if addcpg == 1:
                addcpgfun(self.sequence_display, sequence)
            if addmotifsearch == 1:
                motif = self.motifentry.get()
                addMotifSearch(self.sequence_display, sequence, motif)
            if addcodonprofile == 1:
                addcodonprofilefun(self.sequence_display, sequence)

    def show_in_new_window(self, sequence):
        """Opens the sequence data in a new window if necessary"""
        new_window = Toplevel(self.master)
        new_window.title("Detailed Sequence View")

        sequence_display = Text(new_window, wrap=tk.NONE, height=40, width=100)  # Increased dimensions
        sequence_scrollbar_x = Scrollbar(new_window, orient=tk.HORIZONTAL, command=sequence_display.xview)
        sequence_scrollbar_y = Scrollbar(new_window, orient=tk.VERTICAL, command=sequence_display.yview)

        sequence_display.config(xscrollcommand=sequence_scrollbar_x.set, yscrollcommand=sequence_scrollbar_y.set)

        sequence_scrollbar_x.pack(side=tk.BOTTOM, fill=tk.X)
        sequence_scrollbar_y.pack(side=tk.RIGHT, fill=tk.Y)
        sequence_display.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        addspacer = self.spacer.get()
        addhomopolymer = self.homopolymer.get()
        addcpg = self.cpgisland.get()
        addcodonprofile = self.codonprofile.get()

        if addspacer == 1:
            self.sequence_display.delete(1.0, tk.END)
            addspacerfun(self.sequence_display, sequence)
        else:
            self.sequence_display.insert(tk.END, "FASTA Format\n")
            sequence_display.insert(tk.END, sequence)

        if addhomopolymer == 1:
            addhomopolymerfun(self.sequence_display, sequence)
        if addcpg == 1:
            addcpgfun(self.sequence_display, sequence)
        if addcodonprofile == 1:
            addcodonprofilefun(self.sequence_display, sequence)

    def clear_widget(self):
        """Clears file data and GUI display when clear button is clicked"""
        self.sequence_display.delete(1.0, tk.END)
        for i in self.tree.get_children():
            self.tree.delete(i)
        self.sequence_count_label.config(text=f"")
        self.sequences.clear()


# Display GUI
app = UpdatedFASTAViewer(root)
root.mainloop()
