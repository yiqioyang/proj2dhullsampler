def compute_mi_pair(col_x, col_y):
    """Function to compute MI between two Series."""
    x = col_x.values.reshape(-1, 1)
    y = col_y.values
    mi = mutual_info_regression(x, y, discrete_features=False)
    return mi[0]
