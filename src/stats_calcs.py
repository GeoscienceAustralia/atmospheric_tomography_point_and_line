import numpy as np

def post_stats(x, name):
    values = (name, np.mean(x), np.std(x), np.var(x), np.percentile(x, 2.5), np.percentile(x, 97.5))
    return values


