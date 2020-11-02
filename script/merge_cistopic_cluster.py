import pandas as pd
import sys


def add_new_cluster_information(input, cluster, output):
    pdf = pd.read_csv(input, header=0, index_col=0)
    cdf = pd.read_csv(cluster, header=0)
    df = pdf.merge(cdf, how='left', left_index=True, right_index=True)




if __name__ == "__main__":
    assert len(sys.argv) >= 5
    input   = sys.argv[1]
    cluster = sys.argv[2]
    add_new_cluster_information(input, cluster, output)