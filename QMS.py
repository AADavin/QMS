import seaborn as sns
import itertools
import os
import random
import pandas as pd
import ete3
import sys
import numpy as np

_, folder_with_trees, output_folder = sys.argv


def QMS():

    good_families = list()
    bad_families = list()


    # First we read the trees

    tree_files = [x for x in os.listdir(folder_with_trees) if x.endswith("treefile")]

    # Second we go tree by tree

    for tree_file in tree_files:

        # Now, we remove obvious duplications taking place within one species. We select the
        # copy that is closer to the root
        with open(os.path.join(folder_with_trees, tree_file)) as f:
            tree = ete3.Tree(f.readline().strip(), format=1)

        for l in tree:
            l.add_feature("txid", l.name.split("_")[0])

        # Need to root cleverly here

        tree = filter_repeated_nodes(tree)
        problematic_species = return_problematic_species(tree)
        if len(problematic_species) == 0:
            good_families.append((tree, tree_file, len(tree)))
        else:
            bad_families.append((tree, tree_file, len(tree)))


    # Now remove families under x (in this case 15, you NEED to modify this

    bad_families = [x for x in bad_families if x[2] > 15]
    good_families = [x for x in good_families if x[2] > 15]

    # For each family I am going to recover the ML distance matrix. I use the genes I have already selected above
    # Then, I compute what I call the Global matrix. In here, I compute the mean ML distance between one species and the others. I take into account the fact that some species might have several copies
    # Once I have the Global Matrix, I go through the bad families and I select the genes whose ML distance has the highest Pearson correlation with the Global distance matrix
    # Since the ML distance from one family to another can vary, I think I should normalize this

    # I create also a dict to track the number of times that a species have appeared

    k = 0
    fam_to_df_mldists = dict()

    for mtree, iqtree, ngenes in itertools.chain(good_families, bad_families):
        # This reads the ml dists and returns a df # Need to modify this using the patristic distances
        if ngenes < 40:
            continue
        k += 1
        mldist_file = iqtree.replace("treefile", "mldist")
        with open(os.path.join(mainpath, mldist_file)) as f:
            n = int(f.readline().strip())
            mldist = np.zeros((n, n))
            mindex = list()
            for i, l in enumerate(f):
                gene = l.strip().split(" ")[0]
                mindex.append(gene)
                dists = [x.strip() for x in l.strip().split(" ")[1:] if x != ""]
                for j, vl in enumerate(dists):
                    mldist[i][j] = vl
            df_mldist = pd.DataFrame(mldist, index=mindex, columns=mindex)
        # Now I want to drop rows and columns we already discarded
        removed_leaves = set(df_mldist.columns) - set([x.code for x in mtree.get_leaves()])
        df_mldist = df_mldist.drop(removed_leaves, axis=0)
        df_mldist = df_mldist.drop(removed_leaves, axis=1)

        fam_to_df_mldists[iqtree.split(".")[0]] = df_mldist

    print("Done!")

    # I obtained the MML dists. I obtain the mean values grouping per species

    fam_to_mml = dict()

    k = 0

    for fam, df_mldist in fam_to_df_mldists.items():

        species = list()
        renames = dict()

        # We add a new column with the species

        df_n_mldist = pd.DataFrame.copy(df_mldist)

        new_col = list()

        for name in df_mldist.index:
            sp = txid_to_sp[name.split("_")[0].replace("CYANO", "").replace("EGT", "").replace("PLASTID", "")][1]
            new_col.append(sp)

        df_n_mldist.insert(0, "species", new_col)
        firstmean = df_n_mldist.groupby('species').mean().T

        new_col = list()

        for name in firstmean.index:
            sp = txid_to_sp[name.split("_")[0].replace("CYANO", "").replace("EGT", "").replace("PLASTID", "")][1]
            new_col.append(sp)

        firstmean.insert(0, "species2", new_col)
        df_mmldist = firstmean.groupby("species2").mean()

        # We make zeros the main diagonal

        for sp1 in df_mmldist:
            df_mmldist[sp1][sp1] = 0

        # secondmean = firstmean.groupby('species2').mean()
        # df_mmldist = secondmean/np.max(secondmean.values)   Actually we don't have to normalize

        fam_to_mml[fam] = df_mmldist

        k += 1

    print("Done!")

    # I create the global dataframe

    print(len(species_retained))

    global_dist = np.zeros((len(species_retained), len(species_retained)))
    mindex = [txid_to_sp[x][1] for x in list(species_retained)]
    df_global_dist = pd.DataFrame(global_dist, index=mindex, columns=mindex)

    print("Done!")

    # I obtain the global distances

    species_to_count = dict()

    sps = [txid_to_sp[x][1] for x in species_retained]

    for sp1 in sps:
        if sp1 not in species_to_count:
            species_to_count[sp1] = dict()
        for sp2 in sps:
            if sp2 not in species_to_count[sp1]:
                species_to_count[sp1][sp2] = 0

    for fam, ml in fam_to_mml.items():
        for sp1 in ml.columns:
            for sp2 in ml.columns:
                if sp1 == sp2:
                    continue
                species_to_count[sp1][sp2] += 1
                species_to_count[sp2][sp1] += 1

        df_global_dist = df_global_dist.add(ml, fill_value=0)
    print("Done!")

    # We correct for multiple values

    for sp1, msps in species_to_count.items():

        for sp2, vl in msps.items():

            if sp1 not in df_global_dist.columns or sp2 not in df_global_dist.columns:
                continue

            if vl == 0:
                continue

            df_global_dist.loc[sp1, sp2] = df_global_dist.loc[sp1, sp2] / vl

    print("Done!")

    # We normalize

    # df_global_dist = df_global_dist/np.max(df_global_dist.values)

    # print("Done!")

    # Now I iterate over the normalized PER COLUMN. And I get, for each species with several copies of the gene,
    # the copy with the better Pearson correlation with the global matrix

    best_genes = dict()
    k = 0

    for fam, df_mldist in fam_to_df_mldists.items():

        best_genes[fam] = list()

        # First I do the mean per species per column

        species_to_cols = dict()
        for cindex, column in df_mldist.iteritems():
            cname = txid_to_sp[cindex.split("_")[0].replace("CYANO", "").replace("EGT", "").replace("PLASTID", "")][1]
            if cname not in species_to_cols:
                species_to_cols[cname] = list()
            species_to_cols[cname].append(cindex)

        df_meancols_mldist = pd.DataFrame(columns=species_to_cols.keys())

        for sp, cols in species_to_cols.items():
            df_meancols_mldist[sp] = df_mldist[cols].mean(axis=1)

        df_meancols_mldist = df_meancols_mldist.T

        species_to_cols = dict()

        for cindex, column in df_meancols_mldist.iteritems():
            cname = txid_to_sp[cindex.split("_")[0].replace("CYANO", "").replace("EGT", "").replace("PLASTID", "")][1]
            if cname not in species_to_cols:
                species_to_cols[cname] = list()
            species_to_cols[cname].append(cindex)

        for sp, cols in species_to_cols.items():

            if len(cols) != 1:
                mcorrvalues = list()
                for col in cols:
                    df_meancols_mldist[col][sp] = 0
                    md = pd.DataFrame(columns=[col, sp])
                    md[col] = df_global_dist[sp]
                    md[sp] = df_meancols_mldist[col]
                    mcorrvalues.append((col, (md[col].corr(md[sp]))))

                mychoice = max(mcorrvalues, key=lambda x: x[1])[0]
                best_genes[fam].append(mychoice)

            else:
                best_genes[fam].append(cols[0])
    print("Done!")


def return_problematic_species(tree):
    problematic_species = list()
    species_to_count = dict()

    for n in tree:
        if n.name not in species_to_count:
            species_to_count[n.name] = 0
        species_to_count[n.name] += 1

    for k, v in species_to_count.items():
        if v > 1:
            problematic_species.append((k, v))
    return problematic_species

def filter_repeated_nodes(tree):

    species2number = dict() # Stores the number of copies that a species has
    poly_nodes = list() # Stores genes containing one species with several copies
    mroot = tree.get_tree_root()

    for n in tree.traverse():
        if not n.is_leaf():
            descendant_txids = [x.txid for x in n.get_leaves()]
            if len(descendant_txids) > 1 and len(set(descendant_txids)) == 1:
                if species2number[n.get_leaves()[0].name] == len(descendant_txids):
                    poly_nodes.append(n)

    # We need to select the gene closest to the root

    def checkEqual1(iterator):
        iterator = iter(iterator)
        try:
            first = next(iterator)
        except StopIteration:
            return True
        return all(first == rest for rest in iterator)

    remove = set()

    for n in poly_nodes:
        dists = [(l.get_distance(mroot), l.number) for l in n.get_leaves()]

        # If all dists are equal, we pick at random

        if checkEqual1([x[0] for x in dists]):
            nmin = random.sample(dists, 1)[0][1]
        else:
            nmin = min(dists, key=lambda x: x[0])[1]
        for _, number in dists:
            if number != nmin:
                remove.add(number)

    for l in tree.get_leaves():
        if l.number in remove:
            l.delete(prevent_nondicotomic=True, preserve_branch_length=True)

    return tree











