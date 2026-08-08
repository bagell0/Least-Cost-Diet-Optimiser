"""
Look up USDA ids: a food's fdc_id for price.csv, or a nutrient's id for constraints.csv.

Just run this file (e.g. `python find.py`, or hit Run in your IDE). Pick food or nutrient
search when prompted, then type search terms - it prints every match. Type "food" or
"nutrient" at any point to switch modes. Press Enter with nothing typed to quit.

Food search only looks at SR Legacy foods - Foundation Foods records are frequently missing
vitamin data, so they're excluded here even though they still exist in the full USDA download.
"""

import pandas as pd


def find_food(query, foods):
    foods = foods[foods["data_type"] == "sr_legacy_food"]
    # match if every word in the query appears somewhere in the description, in any order -
    # USDA descriptions follow a fixed "Category, qualifier..." convention that rarely matches
    # how a person phrases the same food (e.g. "Salt, table", not "table salt").
    mask = pd.Series(True, index=foods.index)
    for word in query.split():
        mask &= foods["description"].str.contains(word, case=False, na=False, regex=False)
    return foods[mask]


def find_nutrient(query, nutrients):
    return nutrients[nutrients["name"].str.contains(query, case=False, na=False)]


if __name__ == "__main__":
    print("Loading USDA food list (food.csv) and nutrient list (nutrient.csv)...")
    foods = pd.read_csv("food.csv")
    nutrients = pd.read_csv("nutrient.csv")

    mode = "food"
    while True:
        query = input(f"\nSearch for a {mode} (blank to quit, or type 'food'/'nutrient' to switch): ").strip()
        if not query:
            break
        if query.lower() in ("food", "nutrient"):
            mode = query.lower()
            continue
        if mode == "food":
            matches = find_food(query, foods)
            columns = ("fdc_id", "description")
        else:
            matches = find_nutrient(query, nutrients)
            columns = ("id", "name", "unit_name")
        if matches.empty:
            print("No matches found.")
        else:
            for _, row in matches.iterrows():
                print("\t".join(str(row[c]) for c in columns))
