# QMS
Quick Marker Selection for Species Trees

# This is not ready to use yet

You need to install:
ete3
pandas


Steps:

1. Create gene families of **homologous** sequences (you can use MCL or Silix, for instance)
2. Select gene families that:
    ** Are represented in most species (a rule of thumb, use gene families that have more than 80% of your species)
    ** The number of sequences is similar to the number of species
3. Compute gene trees for those gene families (using for instance IQ tree)
4. Sequences should be named: nameOfSpecies_sequenceIdentifier
5. Put all the resulting trees in the same folder
6. Launch QMS:

python QMS.py folder_with_candidate_trees output_folder_name




