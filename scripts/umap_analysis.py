import numpy as np
import umap
import pandas as pd
import sys
import time

from itertools import product
from pathlib import Path

from sklearn.preprocessing import StandardScaler

# load params
data_input = snakemake.params.data_input
output_folder_tab = f"{snakemake.params.results_folder}/{data_input['rna_class']}_{data_input['detection_rate']}/results_{data_input['data_sub_set']}/matrices/umap/coords"

# stuff for umap to loop over
top_list = snakemake.params.top_list
n_neighbors_s = snakemake.params.n_neighbors
min_dist_s = snakemake.params.min_dist
metric_s = snakemake.params.metric
init_s = snakemake.params.init
seed_s = snakemake.params.seed
norm_s = snakemake.params.norm

# load expr matrix
expr_filename = "{}_{}/{}_{}_quantification_{}.detection_rate_{}_per_{}.csv".format(data_input["raw_data_folder"],
                                                                                                        data_input["data_sub_set"],
                                                                                                        data_input["rna_class"],
                                                                                                        data_input["feature_filtering"],
                                                                                                        data_input["norm"],
                                                                                                        data_input["detection_rate"],
                                                                                                        data_input["detection_group"])

for top_feature in top_list:
    # loading happens in the inner loop to load a "fresh" expression matrix for every umap
    for n_neighbors, min_dist, metric, init, seed, norm in product(n_neighbors_s, min_dist_s, metric_s, init_s, seed_s, norm_s):
        #print(expr_filename)
        expr = pd.read_csv(expr_filename, sep='\t', index_col=0)

        # filter for top features
        if top_feature != "all":
            features_list = pd.read_csv(f"{snakemake.params.results_folder}/{data_input['rna_class']}_{data_input['detection_rate']}/results_{data_input['data_sub_set']}/matrices/most_expressed/top_{top_feature}.csv", sep='\t', index_col=None, header=0)
            expr = expr.loc[features_list[data_input['rna_class']]]
        #else:
        #    feature_list = pd.read_csv(f"{snakemake.params.results_folder}/matrices/most_expressed/{top_feature}.csv", sep='\t', index_col=None, header=0)

        # transpose matrix (Tobias did this in the snakefile)
        tbl = expr.transpose()
        ids = list(tbl.index)
        
        # according to sncrna-analysis pipeline
        output_folder_with_norm = Path(output_folder_tab) 

        if top_feature != "all":
            output_file = output_folder_with_norm / f"umap_top_{top_feature}_expressed_{data_input['rna_class']}s_neigh{n_neighbors}_mdist{min_dist}_{metric}_init{init}_seed{seed}_norm_{norm}.csv"
        else:
            output_file = output_folder_with_norm / f"umap_neigh{n_neighbors}_mdist{min_dist}_{metric}_init{init}_seed{seed}_norm_{norm}.csv"
       
        # print(output_file)
        if norm == "std":
            tbl = StandardScaler().fit_transform(tbl)
        elif norm != "raw":
            print("Unknown normalization chosen! ({})".format(norm))
            sys.exit(1)
        if False:
            print(tbl.shape)
            print(n_neighbors)
            print(min_dist)
            print(metric)
            print(init)
            print(seed)
        emb = umap.UMAP(
            n_neighbors=n_neighbors,
            min_dist=min_dist,
            metric=metric,
            init=init,
            random_state=seed
        ).fit_transform(tbl)

        result = pd.DataFrame(emb, columns=["D1", "D2"], index=ids)
        #print(result.shape)
        # create path
        output_folder_with_norm.mkdir(parents=True, exist_ok=True)
        # save
        result.to_csv(path_or_buf=output_file, sep='\t', header=True, index=True, index_label='ID')

    # create .done file
    with open(snakemake.output[0], "w") as file:
        file.write("")
