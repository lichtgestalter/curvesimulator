import scipy.stats as stats


def chi_square_pvalue(chi_square, n_measurements, n_parameters):
    """
    Calculate the p-value for a chi-square test.
    This is the probability of observing a chi-square value >= your observed value.
    chi_square :     The chi-square test statistic
    n_measurements : Number of measurements/observations
    n_parameters :   Number of free parameters in the model
    """
    df = n_measurements - n_parameters  # degrees of freedom
    p_value = 1 - stats.chi2.cdf(chi_square, df)  # cumulative distribution function
    print(f"Chi-square value:       {chi_square}")
    print(f"Number of measurements: {n_measurements}")
    print(f"Number of parameters:   {n_parameters}")
    print(f"Degrees of freedom:     {df}")
    print(f"P-value:                {p_value:.4f}\n")
    if p_value > 0.05:
        print(f"The fit is acceptable (p > 0.05). There's no significant evidence that your model is inconsistent with the data.")
    else:
        print(f"The fit is poor (p < 0.05). Your model may not adequately describe the data.")
    return p_value


p_value_uli     = chi_square_pvalue(chi_square=31.57, n_measurements=39, n_parameters=13)  # Uli T062
p_value_vitkova = chi_square_pvalue(chi_square=38.74, n_measurements=39, n_parameters=13)  # Vitkova
