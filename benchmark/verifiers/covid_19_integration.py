import csv

def covid_19_integration(d):
    try:
        pred_path = f"{d}/covid_19_integration.csv"
        gold_path = "benchmark/gold_results/covid_19_integration.csv"

        def load_csv_as_set(path):
            """
            Load CSV and convert it into a canonical set representation
            that is insensitive to row/column order and case.
            """
            with open(path, newline='', encoding='utf-8') as f:
                reader = csv.reader(f)
                rows = list(reader)

            if not rows:
                return set()

            header = rows[0]
            col_names = [c.strip().lower() for c in header[1:]]

            result = set()
            for row in rows[1:]:
                row_key = row[0].strip().lower()
                values = row[1:]

                for col, val in zip(col_names, values):
                    v = val.strip().lower()
                    # normalize numeric formatting
                    try:
                        v = str(float(v))
                    except:
                        pass
                    result.add((row_key, col, v))

            return result

        pred_set = load_csv_as_set(pred_path)
        gold_set = load_csv_as_set(gold_path)

        return pred_set == gold_set

    except Exception as e:
        print(f"Error processing CSV files: {e}")
        return False
