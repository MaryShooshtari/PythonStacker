import numpy as np

class bsm_reweighter:
    def __init__(self, order):
        # Set up matrix and invert based on order we go to. Constant component should always be 0!
        self.order = order
        if self.order == 2:
            self.weightvariations = np.array((0.5, 2.))
        else:
            self.weightvariations = np.array((0.1, 0.2, 0.5, 2., 3.))[:order]
        # now build self.matrix. Each row is one order higher of weight var
        # make sure this is not hardcoded:
        self.matrix = np.empty((order, order), dtype=np.float64)
        for i in range(order):
            for j in range(order):
                self.matrix[j][i] = self.weightvariations[j] ** (i * 2 + 2)
        self.matrix_inv = np.linalg.inv(self.matrix)
    def transform_weights(self, weights):
        # Need to swap order of column and weights: weights{2, 1] must become weights[1,2]:
        return np.dot(self.matrix_inv, np.transpose(weights[:, :self.order]))


# Constant component of polynomial is always 0.
def getBSMVariations(order=2):
    return ["Quad", "Quartic"][:order]


def getBSMVariationsGroomed(order=2):
    variations = getBSMVariations(order)
    ret = []
    for bsm_var in variations:
        if len(bsm_var) == 0:
            continue
        bsm_var_name = "BSM_" + bsm_var
        ret.append(bsm_var_name)
    return ret
