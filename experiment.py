import os
from tabulate import tabulate
import matplotlib.pyplot as plt
from ete3 import Tree
import numpy as np
import math

def run_inference(msa_path, model, prefix, args = ""):
    if not os.path.isfile(msa_path):
        print("MSA " + msa_path + " does not exist")
        return
    prefix_dir = "/".join(prefix.split("/")[:-1])
    if not os.path.isdir(prefix_dir):
        os.makedirs(prefix_dir)
    if not os.path.isfile(prefix + ".raxml.bestTree"):
        args = args + " --redo"
    command = "./bin/raxml-ng"
    command += " --msa " + msa_path
    command += " --model " + model
    command += " --prefix " + prefix
    command += " --threads auto --seed 2 "
    command += " " + args
    os.system(command)


def gq_distance(tree_name1, tree_name2):
    if tree_name1 is None or tree_name2 is None:
        return float('nan')
    if tree_name1 != tree_name1 or tree_name2 != tree_name2:
        return float("nan")
    os.system("./bin/qdist " + tree_name1 + " " + tree_name2 + " >out.txt")
    lines = open("out.txt").readlines()
    if len(lines) < 2: #error occurred
        return float('nan')
    res_q = float(lines[1].split("\t")[-3])
    qdist = 1 - res_q
    os.remove("out.txt")
    return qdist


def rf_distance(tree_name1, tree_name2):
    try:
        t1 = Tree(tree_name1)
        t2 = Tree(tree_name2)
    except:
        return float('nan')
    if t1 is None or t2 is None:
        return float('nan')
    if t1 != t1 or t2 != t2:
        return float("nan")
    rf, max_rf, common_leaves, parts_t1, parts_t2,discard_t1, discart_t2 = t1.robinson_foulds(t2, unrooted_trees = True)
    if max_rf == 0:
        return float('nan')
    return rf/max_rf


def plot_distribution(data, label):
    plt.hist(data, 20)
    plt.xlabel(label)
    plt.ylabel('Number of datasets')
    plt.savefig(os.path.join(plots_dir, "hist_" + label +  ".png"))
    plt.clf()
    plt.close()

num_samples = 100
plots_dir = "data/plots"
for dataset in os.listdir("data/msa"):
    if dataset in ["abvdoceanic", "bowernpny"]:
        continue
    glottolog_tree_path = os.path.join("data/glottolog_trees", dataset, "glottolog.tree")
    if not os.path.isfile(glottolog_tree_path):
        continue
    bin_msa_path = os.path.join("data/msa/", dataset, "bin.phy")
    #run_inference(bin_msa_path, "BIN+G", os.path.join("data/raxml/", dataset, "bin"))
    dolgo_msa_path = os.path.join("data/sound_msa/", dataset, "dolgo.phy")
    #run_inference(dolgo_msa_path, "MULTI14_MK+M{VKPHJMNSRTW+1_}{-}", os.path.join("data/raxml/", dataset, "dolgo"))
    dolgo_catg_msa_path = os.path.join("data/sound_msa/", dataset, "dolgo.catg")
    #run_inference(dolgo_catg_msa_path, "MULTI15_MK+M{VKPHJMNSRTW+1_~}{-}", os.path.join("data/raxml/", dataset, "dolgo_catg"), args = "--prob-msa on")
    bin_samples_dir = os.path.join("data/msa/", dataset, "samples")
    bin_prefix = os.path.join("data/raxml/", dataset, "bin_samples")
    for i in range(num_samples):
        bin_msa_path = os.path.join(bin_samples_dir, "sample" + str(i) + "_bin.phy")
        #run_inference(bin_msa_path, "BIN+G", os.path.join(bin_prefix, "sample" + str(i) + "_bin"))
    dolgo_samples_dir = os.path.join("data/sound_msa/", dataset, "samples")
    dolgo_prefix = os.path.join("data/raxml/", dataset, "dolgo_samples")
    for i in range(num_samples):
        dolgo_msa_path = os.path.join(dolgo_samples_dir,  "sample" + str(i) + "_dolgo.phy")
       # run_inference(dolgo_msa_path, "MULTI14_MK+M{VKPHJMNSRTW+1_}{-}", os.path.join(dolgo_prefix, "sample" + str(i) + "_dolgo"))
results = []
bin_dists = []
dolgo_dists = []
dolgo_catg_dists = []
bin_sample_dists = []
dolgo_sample_dists = []
rf_bin = []
rf_dolgo = []
for dataset in os.listdir("data/raxml"):
    glottolog_tree_path = os.path.join("data/glottolog_trees", dataset, "glottolog.tree")
    best_bin_tree_path = os.path.join("data/raxml/", dataset, "bin.raxml.bestTree")
    best_dolgo_tree_path = os.path.join("data/raxml/", dataset, "dolgo.raxml.bestTree")
    best_dolgo_catg_tree_path = os.path.join("data/raxml/", dataset, "dolgo_catg.raxml.bestTree")
    bin_prefix = os.path.join("data/raxml/", dataset, "bin_samples")
    gq_bin = gq_distance(glottolog_tree_path, best_bin_tree_path)
    bin_sample_dists.append([gq_distance(glottolog_tree_path, os.path.join(bin_prefix, "sample" + str(i) + "_bin.raxml.bestTree")) for i in range(num_samples)])
    rf_bin.append([rf_distance(best_bin_tree_path, os.path.join(bin_prefix, "sample" + str(i) + "_bin.raxml.bestTree")) for i in range(num_samples)])
    gq_dolgo = gq_distance(glottolog_tree_path, best_dolgo_tree_path)
    gq_dolgo_catg = gq_distance(glottolog_tree_path, best_dolgo_catg_tree_path)
    dolgo_prefix = os.path.join("data/raxml/", dataset, "dolgo_samples")
    dolgo_sample_dists.append([gq_distance(glottolog_tree_path, os.path.join(dolgo_prefix, "sample" + str(i) + "_dolgo.raxml.bestTree")) for i in range(num_samples)])
    rf_dolgo.append([rf_distance(best_dolgo_catg_tree_path, os.path.join(dolgo_prefix, "sample" + str(i) + "_dolgo.raxml.bestTree")) for i in range(num_samples)])
    results.append([dataset, gq_bin, gq_dolgo, gq_dolgo_catg])
    bin_dists.append(gq_bin)
    dolgo_dists.append(gq_dolgo)
    dolgo_catg_dists.append(gq_dolgo_catg)
print(tabulate(results, tablefmt="pipe", floatfmt=".3f", headers = ["dataset", "gq bin", "gq dolgo", "gq dolgo catg"]))


plot_distribution([sum(dists) / len(dists) for dists in rf_bin], "mean_rf_bin")
plot_distribution([sum(dists) / len(dists) for dists in rf_dolgo], "mean_rf_dolgo")
plot_distribution([max(dists) for dists in rf_bin] , "max_rf_bin")
plot_distribution([max(dists) for dists in rf_dolgo] , "max_rf_dolgo")
plot_distribution([np.std(dists) for dists in bin_sample_dists] , "std_gqd_bin")
plot_distribution([np.std(dists) for dists in dolgo_sample_dists] , "std_gqd_dolgo")

plt.axline([0, 0], slope=1, color = 'lightgray', linewidth = 1, linestyle = "--")
plt.scatter(bin_dists, [sum(dists) / len(dists) for dists in bin_sample_dists], s=10)
plt.xlabel("gq_full_bin")
plt.ylabel("avg_gq_samples_bin")
plt.savefig(os.path.join(plots_dir, "scatter_samples_bin.png"))
plt.clf()
plt.close()

plt.axline([0, 0], slope=1, color = 'lightgray', linewidth = 1, linestyle = "--")
plt.scatter(dolgo_catg_dists, [sum(dists) / len(dists) for dists in dolgo_sample_dists], s=10)
plt.xlabel("gq_full_dolgo_catg")
plt.ylabel("avg_gq_samples_dolgo")
plt.savefig(os.path.join(plots_dir, "scatter_samples_dolgo.png"))
plt.clf()
plt.close()


plt.axline([0, 0], slope=1, color = 'lightgray', linewidth = 1, linestyle = "--")
plt.scatter(bin_dists, dolgo_catg_dists, s = 10)
plt.xlabel("bin")
plt.ylabel("dolgo")
plt.savefig(os.path.join(plots_dir, "scatter_bin_dolgo_catg.png"))
plt.clf()
plt.close()

plt.axline([0, 0], slope=1, color = 'lightgray', linewidth = 1, linestyle = "--")
plt.scatter(bin_dists, [sum(dists) / len(dists) for dists in dolgo_sample_dists], s=10)
plt.xlabel("gq_full_bin")
plt.ylabel("avg_gq_samples_dolgo")
plt.savefig(os.path.join(plots_dir, "scatter_samplesi_dolgo_bin.png"))
plt.clf()
plt.close()
