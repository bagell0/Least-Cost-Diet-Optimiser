import pandas as pd
import numpy as np
from scipy.optimize import linprog
import plotly.graph_objects as go
from plotly.subplots import make_subplots

def load_constraints(path):
    # ids column is a comma-separated string of USDA nutrient ids, most-reliable first (fallback chain)
    df = pd.read_csv(path)
    nutrient_pool, lower_vector, upper_vector = {}, {}, {}
    for _, row in df.iterrows():
        if pd.isna(row["nutrient"]):  # skip blank/trailing rows
            continue
        name = row["nutrient"]
        nutrient_pool[name] = [int(x) for x in str(row["ids"]).split(",")]
        if pd.notna(row["lower"]):
            lower_vector[name] = row["lower"]
        if pd.notna(row["upper"]):
            upper_vector[name] = row["upper"]
    return nutrient_pool, lower_vector, upper_vector


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
                    # USDA's "by difference" nutrients (e.g. carbs = 100 - moisture - protein - fat - ash)
                    # can come out slightly negative on near-zero foods due to rounding error elsewhere;
                    # negative nutrient content is physically impossible, so floor it at 0.
                    temp.append(max(amount.iloc[0], 0))
                    found = True
            if not found:
                temp.append(0)
        matrix[n] = temp
    return matrix


def display_name(f_id, p):
    row = p[p["fdc_id"] == f_id].iloc[0]
    if pd.notna(row["product_name"]) and str(row["product_name"]).strip():
        return row["product_name"]
    # no retail product_name set yet - fall back to a short USDA description
    # (first word alone is too lossy, e.g. "Nuts, brazilnuts, raw" -> "Nuts"; keep two)
    return ", ".join(str(row["usda_dscr"]).split(",")[:2]).strip()


# Sub-nutrients whose absence is expected, not a gap, once their parent is ~0 for a food (e.g. no
# amino acid panel on a near-zero-protein oil). Maps each id to its parent's name so
# audit_missing_data can suppress those automatically instead of needing a per-food allowlist.
AMINO_ACID_IDS = {1210, 1211, 1212, 1213, 1214, 1215, 1216, 1217, 1218, 1219, 1220, 1221, 1222,
                   1223, 1224, 1225, 1226, 1227, 1228}
OMEGA3_IDS = {1404, 1278, 1272}
SUBNUTRIENT_PARENT = {**{i: "Protein" for i in AMINO_ACID_IDS}, **{i: "Fat" for i in OMEGA3_IDS}}


def audit_missing_data(n_p, f_o, fdf, p):
    # build_matrix defaults to 0 when no fallback id has data for a food - indistinguishable from
    # genuinely having none, unless flagged. Grouped one line per food (not per food-nutrient pair)
    # so a food with a long list stands out as having a thin USDA panel worth reconsidering.
    gaps_by_food = {}
    for n, n_ids in n_p.items():
        parent_name = next((SUBNUTRIENT_PARENT[i] for i in n_ids if i in SUBNUTRIENT_PARENT), None)
        parent_ids = n_p.get(parent_name, []) if parent_name else []
        for f in f_o:
            has_data = not fdf[(fdf["fdc_id"] == f) & (fdf["nutrient_id"].isin(n_ids))].empty
            if has_data:
                continue
            if parent_ids:
                parent_amt = fdf[(fdf["fdc_id"] == f) & (fdf["nutrient_id"].isin(parent_ids))]["amount"]
                if not parent_amt.empty and parent_amt.iloc[0] < 0.5:
                    continue
            gaps_by_food.setdefault(f, []).append(n)
    for f, missing in gaps_by_food.items():
        name = display_name(f, p)
        print(f"WARNING: {name} has no recorded data for: {', '.join(missing)} "
              f"(defaulting to 0, which may be inaccurate).")


def audit_calories(n_matrix, f_o, p, threshold=0.3):
    # A missing energy value defaults to 0 (see build_matrix), which the LP treats as free calories.
    # Flags foods whose reported kcal doesn't roughly match the Atwater estimate (4*protein +
    # 4*carbs + 9*fat) - usually a sign the record has no energy value at all. Threshold is wide
    # (30%) since high-fibre foods miss the estimate by ~20% on their own; a genuinely missing
    # value is off by ~100%.
    if not {"Calories", "Protein", "Carbs", "Fat"} <= n_matrix.keys():
        return
    for i, f in enumerate(f_o):
        protein, carbs, fat, kcal = (n_matrix["Protein"][i], n_matrix["Carbs"][i],
                                      n_matrix["Fat"][i], n_matrix["Calories"][i])
        estimate = 4 * protein + 4 * carbs + 9 * fat
        if estimate > 0 and abs(kcal - estimate) / estimate > threshold:
            name = display_name(f, p)
            print(f"WARNING: {name} reports {kcal:.0f} kcal/100g but its macros imply "
                  f"~{estimate:.0f} kcal/100g - USDA data may be incomplete for this food.")


def optimise(f_o, n_o, cost_v, min_v, max_v, lower_v, upper_v, matrix):
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
    bounds = np.array([(min_v[f_id], max_v[f_id]) for f_id in f_o])
    return linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method="highs")


def create_labels(lower_v, upper_v, n_o):
    labels_array = []
    for n in n_o:
        if n in lower_v:
            labels_array.append(f"{n} >= {lower_v[n]}")
        if n in upper_v:
            labels_array.append(f"{n} <= {upper_v[n]}")
    return labels_array


def display_results(r, p, n_matrix, l_v, u_v, n_o, f_o):
    A = np.array([n_matrix[nutrient] for nutrient in n_o])  # nutrient_matrix but in order
    totals = A @ r.x  # individual nutrient totals
    individual_totals = {nutrient: (A * r.x)[i] for i, nutrient in enumerate(n_o)}  # per-food, indexed like f_o
    inequality_labels = create_labels(l_v, u_v, n_o)

    def fmt(v):
        return v if isinstance(v, str) else f"{v:g}"

    def decimal_places(v):
        s = f"{v:g}"
        return len(s.split(".")[1]) if "." in s else 0

    used = [i for i, f in enumerate(f_o) if r.x[i] >= 1e-6]
    bar_names = [display_name(f_o[i], p) for i in used]
    bar_grams = [r.x[i] * 100 for i in used]
    bar_stores = [p[p["fdc_id"] == f_o[i]]["store"].iloc[0] for i in used]
    bar_hover = [
        "<b>%s</b><br>%.0fg<br>%s" % (
            display_name(f_o[i], p), r.x[i] * 100,
            "<br>".join(f"{n}: {individual_totals[n][i]:.1f}" for n in n_o),
        )
        for i in used
    ]

    diet_rows = list(zip(bar_names, [f"{g:.0f}g" for g in bar_grams]))

    # Actual is shown to whichever nutrient's own Lower/Upper decimal precision is (e.g. 3 d.p. for
    # a 0.125g floor, 0 d.p. for a whole-number one) - a flat precision either hides the digit that
    # decides if a tight gram-scale floor is met, or over-precises a mg-scale total for no reason.
    nutrient_rows = []
    for i, n in enumerate(n_o):
        bound_decimals = [decimal_places(v) for v in (l_v.get(n), u_v.get(n)) if v is not None]
        precision = max(bound_decimals) if bound_decimals else 1
        nutrient_rows.append((n, fmt(l_v.get(n, "-")), f"{totals[i]:.{precision}f}", fmt(u_v.get(n, "-"))))

    shadow_rows = [(label, f"{price:.4f}")
                   for label, price in zip(inequality_labels, r.ineqlin.marginals) if abs(price) > 1e-6]
    if not shadow_rows:
        shadow_rows = [("(none binding)", "")]

    # per-nutrient pie data, switched via the buttons below (see button block for why buttons).
    pie_values_by_nutrient = {n: [individual_totals[n][i] for i in used] for n in n_o}
    first_n = n_o[0]
    pie_steps = [
        dict(
            label=n, method="restyle",
            args=[{"values": [pie_values_by_nutrient[n]], "title.text": f"{n}<br>total: {totals[i]:.1f}"}, [4]],
        )
        for i, n in enumerate(n_o)
    ]

    # row heights are sized off actual row counts, not fixed fractions - the nutrient/shadow-price
    # tables grow as more nutrients are added to constraints.csv, and a fixed layout clips rows
    # silently (the data is still there, it's just not visible) once they outgrow it.
    bar_h = max(300, 45 * len(used) + 130)
    diet_h = max(150, 34 * len(used) + 90)
    nutrient_h = max(150, 34 * len(n_o) + 90)
    shadow_h = max(120, 34 * len(shadow_rows) + 90)
    buttons_h = 44 * len(n_o) + 30  # always-visible button column beside the pie, one per nutrient
    pie_h = max(600, buttons_h)  # row needs to fit whichever of the two is taller
    # make_subplots renormalizes row_heights by their own sum, then shrinks the result by
    # (rows-1) vertical_spacing gaps - so bar_h/total_h doesn't put bar_h pixels on screen, it puts
    # fewer. Solving total_h to exactly cancel that shrink keeps every row at its intended size.
    vertical_spacing = 0.06
    row_h_sum = bar_h + diet_h + nutrient_h + shadow_h + pie_h
    margin_px = 100  # top (90) + bottom (10) margins set in update_layout below
    total_h = row_h_sum / (1 - 4 * vertical_spacing) + margin_px

    fig = make_subplots(
        rows=5, cols=1,
        row_heights=[bar_h / total_h, diet_h / total_h, nutrient_h / total_h, shadow_h / total_h, pie_h / total_h],
        specs=[[{"type": "bar"}], [{"type": "table"}], [{"type": "table"}], [{"type": "table"}], [{"type": "domain"}]],
        subplot_titles=("Grams per day (hover for full breakdown)", "Diet summary", "Nutrient totals",
                         "Shadow prices", "Nutrient breakdown by food (buttons on the left switch nutrient)"),
        vertical_spacing=vertical_spacing,
    )
    fig.add_trace(go.Bar(
        y=bar_names, x=bar_grams, orientation="h",
        marker_color="#5B8DEF", hovertext=bar_hover, hoverinfo="text",
    ), row=1, col=1)
    for name, grams, store in zip(bar_names, bar_grams, bar_stores):
        fig.add_annotation(
            x=grams, y=name, text=f"Store: {store}", showarrow=False,
            xanchor="left", xshift=8, font=dict(size=10, color="#9AA0A6"),
            row=1, col=1,
        )
    fig.add_trace(go.Table(
        header=dict(values=["Product", "Grams"], fill_color="#2A2A2A", font=dict(color="white")),
        cells=dict(values=list(zip(*diet_rows)), fill_color="#161616", font=dict(color="white"), height=26),
    ), row=2, col=1)
    fig.add_trace(go.Table(
        header=dict(values=["Nutrient", "Lower", "Actual", "Upper"], fill_color="#2A2A2A", font=dict(color="white")),
        cells=dict(values=list(zip(*nutrient_rows)), fill_color="#161616", font=dict(color="white"), height=26),
    ), row=3, col=1)
    fig.add_trace(go.Table(
        header=dict(values=["Constraint", "Shadow price (£/unit)"], fill_color="#2A2A2A", font=dict(color="white")),
        cells=dict(values=list(zip(*shadow_rows)), fill_color="#161616", font=dict(color="white"), height=26),
    ), row=4, col=1)
    fig.add_trace(go.Pie(
        labels=bar_names, values=pie_values_by_nutrient[first_n], hole=0.4,
        title=dict(text=f"{first_n}<br>total: {totals[0]:.1f}"),
        hovertemplate="%{label}<br>%{value:.1f} (%{percent})<extra></extra>",
    ), row=5, col=1)

    # Always-visible buttons, not a dropdown (scroll-hijacks once its open list overflows the
    # viewport) or a slider (hides labels once there are too many ticks). Buttons sit in a left
    # column, pie fills the rest - both spanning the row's full height. (An earlier version shrunk
    # the pie's *vertical* span to make room above it, but the buttons are only as wide as their
    # labels, so that just left a tall blank strip beside them instead of stacking cleanly.)
    pie_domain = fig.data[4].domain
    buttons_frac = 0.22  # width reserved for the button column, wide enough for the longest labels
    col_width = pie_domain.x[1] - pie_domain.x[0]
    fig.data[4].update(domain=dict(x=[pie_domain.x[0] + buttons_frac * col_width, pie_domain.x[1]], y=list(pie_domain.y)))

    fig.update_layout(
        template="plotly_dark",
        paper_bgcolor="#0E0E0E", plot_bgcolor="#0E0E0E",
        title=f"Daily cost: £{r.fun:.2f}   ·   "
              f"Total food mass (as-purchased/raw weight): {sum(r.x) * 100:.0f}g/day",
        showlegend=False,
        height=total_h,
        margin=dict(l=10, r=10, t=90, b=10),
        updatemenus=[dict(
            # Plotly hardcodes the *active* button's fill to a near-white color no matter what
            # bgcolor is set here - rather than fight that, every button is styled light with dark
            # text, so the active one (forced light) and inactive ones (set light) both stay readable.
            type="buttons", direction="down", buttons=pie_steps, active=0,
            bgcolor="#E8E8E8", bordercolor="#888", font=dict(color="#111111", size=14),
            pad=dict(t=6, b=6, l=10, r=10),
            x=pie_domain.x[0], xanchor="left",
            y=pie_domain.y[1], yanchor="top",
        )],
    )
    fig.update_xaxes(title_text="grams/day", row=1, col=1)
    fig.show()


def main():
    food_nutrients = pd.read_csv('food_nutrient.csv', low_memory=False)
    prices = pd.read_csv('price.csv')

    food_pool = prices["fdc_id"].values

    nutrient_pool, lower_vector, upper_vector = load_constraints('constraints.csv')

    food_nutrients_filtered = filter_df(food_pool, nutrient_pool, food_nutrients)
    cost_vector = {}
    min_vector = {}
    max_vector = {}
    for food_id in food_pool:
        info = prices[prices["fdc_id"] == food_id]
        cost_vector[food_id] = (info["pack_price_gbp"].iloc[0] / info["pack_size_g"].iloc[0]) * 100
        # a blank min_serving_100g must default to 0, not NaN - scipy treats a NaN bound as
        # unbounded, which would let the LP select a negative (physically meaningless) serving
        min_vector[food_id] = info["min_serving_100g"].iloc[0] if pd.notna(info["min_serving_100g"].iloc[0]) else 0
        max_vector[food_id] = info["max_serving_100g"].iloc[0]

    nutrient_matrix = build_matrix(food_pool, nutrient_pool, food_nutrients_filtered)
    audit_missing_data(nutrient_pool, food_pool, food_nutrients_filtered, prices)
    audit_calories(nutrient_matrix, food_pool, prices)

    food_order = list(food_pool)
    nutrient_order = list(nutrient_pool.keys())

    result = optimise(food_order, nutrient_order, cost_vector, min_vector, max_vector, lower_vector, upper_vector,
                       nutrient_matrix)

    if result.success:
        display_results(result, prices, nutrient_matrix, lower_vector, upper_vector, nutrient_order, food_order)
    else:
        print(result.message)


if __name__ == "__main__":
    main()
