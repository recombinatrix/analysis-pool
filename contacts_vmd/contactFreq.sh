# use this script to define the contact frequency calculations

# first define atom selection macros

atomselect macro GlyT2 {protein and not resname SUB}
atomselect macro ligand {resname OP48 and resid 1234}   

# you can select a specific functional group
atomselect macro head {resname OP48 and resid 1234 and name CA CB CG HA HB HG }  
atomselect macro tail {resname OP48 and resid 1234 and name C1 C2 C3 H1 H2 H3  }  

# then run the calculation 
#
# syntax: 
#   contactFreq $sel1 $sel2 $dist $threshold $output_filename.dat
#
#   residues are in contact in a given frame if the minimum distance between any atom of sel1 and any atom of sel2 is less than or equal to $dist
#   contact frequency is only reported for contacts present in at least $threshold percent of the simulation
#   output is a whitespace delimited text file

contactFreq {GlyT2} {ligand} 4 0 contact-GlyT2_ligand.dat
contactFreq {GlyT2} {head} 4 0 contact-GlyT2_ligand-head.dat
contactFreq {GlyT2} {tail} 4 0 contact-GlyT2_ligand-tail.dat

