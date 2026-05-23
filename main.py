import pandas as pd
import numpy as np
from scipy.optimize import linprog
import matplotlib.pyplot as plt


def filter_df(f_p, n_p, fdf):
    relevant_nutrient_ids = [n_id for n_ids in n_p.values() for n_id in n_ids]
    return fdf[
        (fdf["fdc_id"].isin(f_p)) & (fdf["nutrient_id"].isin(relevant_nutrient_ids))]


def build_matrix(f_p, n_p, fdf):
    matrix = {}
    for n, n_ids in n_p.items():
        temp = []
        for f_id in f_p:
            i = 0
            found = False
            while i < len(n_ids) and not found:
                n_id = n_ids[i]
                amount = fdf[(fdf["fdc_id"] == f_id) & (
                        fdf["nutrient_id"] == n_id)]["amount"]
                if amount.empty:
                    i += 1
                else:
                    temp.append(amount.iloc[0])
                    found = True
            if not found:
                temp.append(0)
        matrix[n] = temp
    return matrix


def optimise(f_o, n_o, cost_v, ctrl_v, lower_v, upper_v, matrix):
    c = np.array([cost_v[f_id] for f_id in f_o])
    A_ub_rows = []
    b_ub_values = []
    for n in n_o:
        if n in lower_v:
            A_ub_rows.append([-amount for amount in matrix[n]])
            b_ub_values.append(-lower_v[n])
        if n in upper_v:
            A_ub_rows.append(matrix[n])
            b_ub_values.append(upper_v[n])
    A_ub = np.array(A_ub_rows)
    b_ub = np.array(b_ub_values)
    bounds = np.array([(0, ctrl_v[f_id]) for f_id in f_o])
    return linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method="highs")


def plot_diet(ax, f_names, r):
    names = list(f_names.values())
    ax.barh(names, r.x * 100)
    ax.set_xlabel("grams/day")
    return


def create_labels(lower_v, upper_v, n_o):
    labels_array = []
    for n in n_o:
        if n in lower_v:
            labels_array.append(f"{n} >= {lower_v[n]}")
        if n in upper_v:
            labels_array.append(f"{n} <= {upper_v[n]}")
    return labels_array


def main():
    food_nutrients = pd.read_csv('food_nutrient.csv', low_memory=False)
    prices = pd.read_csv('price.csv')

    food_pool = prices["fdc_id"].values
    nutrient_pool = {"Calories": [2047, 1008, 2048], "Protein": [1053, 1003], "Carbs": [2039, 1005, 1050],
                     "Fat": [1004, 1085]}

    food_nutrients_filtered = filter_df(food_pool, nutrient_pool, food_nutrients)
    cost_vector = {}
    control_vector = {}
    for food_id in food_pool:
        info = prices[prices["fdc_id"] == food_id]
        cost_vector[food_id] = (info["pack_price_gbp"].iloc[0] / info["pack_size_g"].iloc[0]) * 100
        control_vector[food_id] = info["max_serving_100g"].iloc[0]

    nutrient_matrix = build_matrix(food_pool, nutrient_pool, food_nutrients_filtered)

    # Locking the order
    food_order = list(food_pool)
    nutrient_order = list(nutrient_pool.keys())
    A = np.array([nutrient_matrix[nutrient] for nutrient in nutrient_order])  # nutrient_matrix but in order
    food_names = {f: prices[prices["fdc_id"] == f]["description"].iloc[0].split(",")[0] for f in food_order}

    lower_vector = {"Protein": 150, "Calories": 3300, "Fat": 75}  # lower limits
    upper_vector = {"Calories": 3400, "Fat": 100}  # upper limits
    inequality_labels = create_labels(lower_vector, upper_vector, nutrient_order)

    result = optimise(food_order, nutrient_order, cost_vector, control_vector, lower_vector, upper_vector,
                      nutrient_matrix)
    totals = A @ result.x  # individual nutrient totals

    # result.x = [a1, a2, a3, ... ax]
    individual_totals = {nutrient: (A * result.x)[i] for i, nutrient in enumerate(nutrient_order)}

    print(f"Daily cost: £{result.fun:.2f}")
    for label, price in zip(inequality_labels, result.ineqlin.marginals):
        if abs(price) > 1e-6:
            print(f"{label}: shadow price = {price:.4f}£/unit")
    print(f"Total food mass: {sum(result.x) * 100:.0f}g/day")
    print()

    for i, n in enumerate(nutrient_order):
        L = lower_vector.get(n, "-")
        U = upper_vector.get(n, "-")
        print(f"{n}: {L} <= {totals[i]:.1f} <= {U}")

    # to display  x: calories, ... ->
    for i, f in enumerate(food_order):
        if result.x[i] < 1e-6:
            continue
        name = food_names[f]
        print(f"{name}: {result.x[i] * 100:.0f}g")
        for n in nutrient_order:
            print(f"{n}: {individual_totals[n][i]:.1f}")

    fig, axes = plt.subplots()
    plot_diet(axes, food_names, result)
    fig.tight_layout()
    plt.show()
    return


if __name__ == "__main__":
    main()
