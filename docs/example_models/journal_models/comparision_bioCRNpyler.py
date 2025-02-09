from biocrnpyler import *

# Define a set of DNA parts
# this is a promoter repressed by tetR and has a leak reaction
ptet = RegulatedPromoter("ptet", ["tetr"], leak=True)
# constitutive promoter
pconst = Promoter("pconst")
# the Combinations A and B or just A or just B be transcribed
pcomb = CombinatorialPromoter("pcomb", ["arac","laci"], leak=False, tx_capable_list=[["arac"], ["laci"]])
# regular RBS
utr1 = RBS("UTR1")
# regular RBS
utr2 = RBS("UTR1")
# a CDS has a name and a protein name. so this one is called GFP and the protein is also called GFP
gfp = CDS("GFP", "GFP")
# you can say that a protein has no stop codon.

# This is a little different from a fusion protein, because in this case you are saying that the ribosome reads
# through two proteins but still produces two distinct proteins, rather than one fused protein.
# This can happen in the case of the ta peptide which causes a peptide bond not to be formed while making a protein.
fusrfp = CDS("fusRFP", "RFP", no_stop_codons=["forward"])

# regular RFP
rfp = CDS("RFP", "RFP")
# cfp
cfp = CDS("CFP", "CFP")
# a terminator stops transcription
t16 = Terminator("t16")

# Combine the parts together in a DNA_construct with their directions
construct = DNA_construct([[ptet, "forward"], [utr1, "forward"], [gfp, "forward"], [t16, "forward"],
                           [t16, "reverse"], [rfp, "reverse"], [utr1, "reverse"], [pconst, "reverse"]])

# some very basic parameters are defined - these are sufficient for the whole model to compile!
parameters = {"cooperativity": 2, "kb": 100, "ku": 10, "ktx": .05, "ktl": .2, "kdeg": 2, "kint": .05}

# Place the construct in a context (TxTlExtract models a bacterial lysate
# with machinery like Ribosomes and Polymerases modelled explicitly)
myMixture = TxTlExtract(name="txtl", parameters=parameters, components=[construct])

# Compile the CRN
myCRN = myMixture.compile_crn()
print(myCRN.pretty_print(show_rates=True, show_keys=True))

# plotting not shown -
# but BioCRNpyler automatically produces interactive reaction network graphs to help visualize and debug complex CRNs!