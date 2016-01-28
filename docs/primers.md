# Primer Design 101

*Note that this tutorial is in development.  If you are not associated with the GERMS lab or the DARTE-QM research project, this tutorial is not meant for you.

This documentation guides how to design primers from a list of genes.  It should be proceeded by a good understanding of our in-house primer design pipeline until someone gets around to putting together the documentation of the workflow or better yet publishes it.

This tutorial is also written specifically for performing on an ubuntu-based server.  But...the caveat is...I have NO IDEA where the software is installed but you will need to in order to run this code.  So if you don't understand paths...I suggest you go find your favorite expert, buy them a coffee, and get a 5 minute explanation.

### Choose your favorite genes.

The first thing you need to do is go get a FASTA formatted sequencing file of genes of interest.  How do to this is entirely up to you -- but I recommend that you pick fairly well-curated genes.  Junk in - Junk out...  You will need at minimum the nucleotide sequences and ideally (and for this tutorial I am assuming) both nucleoitide and amino acids.

If you want to try this out with tutorial files, they live [here](https://github.com/germs-lab/Primer_Design/tree/master/tutorial_seqs).  These can be obtained with:  

    git clone https://github.com/germs-lab/Primer_Design.git

And are in the *tutorial_seqs* directory.  This is private github repository and I will need to grant you access.  However, access to this repository needs to be granted in order for you.  

### Generate a HMM for your gene of interest.

If you are using RDP Fungene, you will likely have model that you can use.  If not, you will need to read the HMMER manual and generate a model based on a sequence alignment. Here's the [manual](http://hmmer.org/).

If you want to align using ClustalO:

    clustlo -i <gene.aa.fa> > <gene.aa_aln.fa>
### Align your genes to a model

Prior to this step, you need to make sure that the names of your genes in your FASTA file make sense to you and that they don't contain funny characters like @!*#! because they will @!*#! your results.

I am also assuming at this point that we are working with a protein sequence, our preferred route for now.

To align genes, we will use the HMMER aligner.

    hmmalign -o <gene.stk> <gene.hmm> <gene.faa>

Then, you need to do some file parsing and fixing which involves multiple steps and you'll realize why you want to learn automation in the shell.

This command will convert your somewhat unreadable file to a readable file that can be read by the RDP tools:

    esl-reformat --informat stockholm afa <gene.stk> > <gene.aln_aa.fa>

*Note that the input is a HMMER alignment and the output is an aligned FASTA file.

Next, we will retrieve the alignment of the nucleotides that match these encoded and now aligned amino acids using RDP Tools. The inputs are the aligned protein sequences and the un-aligned nucleotide sequences.  The output is some statistical summary and and aligned nucleotide sequence file.

    java -jar AlignmentTools.jar align-nucl-to-prot <gene.aln_aa.fa> <gene.nuc.fa> <gene.aln_nuc.fa> <gene.aln.stat>

### Making a tree (optional)

This is an optional step for the primer design.  A bit of a phylogentic guided of an approach.  To build a tree of your aligned sequence that can be inputted into the Primer Design pipeline (for a nucleotide sequence file):

    fasttree -nt -gtr <gene.aln_nuc.fa> > <gene.aln_nuc.nwk>

### Getting your primers!

The input into this script is going to be

With tree:

    java -jar PrimerDesign.jar  select-primer-pipeline -i gene.aln_na.fa -t gene.aln_na.nwk --productLengthMin 150 --isHenikoffWeightNeeded false

Without tree:

    java -jar PrimerDesign.jar  select-primer-pipeline -i gene.aln_na.fa --productLengthMin 150 --isHenikoffWeightNeeded false --isTreeWeightNeeded false

To see the alignment on a graph:

    java -jar PrimerDesign.jar screen-oligo-pipeline -i gene.aln_na.fa -o gene_graphs

Other options can be explored using:

    java -jar PrimerDesign_v0.0.1/PrimerDesign.jar --help
