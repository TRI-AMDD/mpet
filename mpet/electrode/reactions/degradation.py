def degradation(eta,k0):
    Rate = -k0 * np.exp(eta)
    # lithium plating/stripping cannot happen if there is no Li
    return Rate
