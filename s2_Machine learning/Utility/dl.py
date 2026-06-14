import numpy as np
import scipy
from scipy.stats import entropy

def delong_test(pred_proba, true_labels, pred_labels, num_classes):

    pred_entropy = entropy(pred_proba)

    true_entropy = entropy(true_labels)

    pred_labels_entropy = entropy(pred_labels, base=num_classes)

    delong_stat = np.mean((pred_entropy - true_entropy) ** 2 - (pred_labels_entropy - true_entropy) ** 2)

    p_value = 1 - scipy.stats.chi2.cdf(delong_stat, df=2 * num_classes)
    return p_value
