import csv
import os

def plantpad_annotate(d):
    try:
        pred_path = os.path.join(d, "plantpad_annotate.csv")
        gold_path = "benchmark/gold_results/plantpad_annotate.csv"

        if not os.path.exists(pred_path):
            print(f"Prediction file not found: {pred_path}")
            return False

        if not os.path.exists(gold_path):
            print(f"Gold file not found: {gold_path}")
            return False

        def load_csv_as_set(path):
            rows_set = set()
            with open(path, newline='', encoding="utf-8") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    # normalize: lowercase keys and values, strip spaces
                    normalized = {
                        k.strip().lower(): v.strip().lower()
                        for k, v in row.items()
                    }
                    # convert dict to a sorted tuple so order doesn't matter
                    rows_set.add(tuple(sorted(normalized.items())))
            return rows_set

        pred_rows = load_csv_as_set(pred_path)
        gold_rows = load_csv_as_set(gold_path)

        return pred_rows == gold_rows

    except Exception as e:
        print(f"Error processing CSV files: {e}")
        return False
