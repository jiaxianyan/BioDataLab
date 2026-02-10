import csv

def cellstar_integration(d):
    try:
        pred_path = f"{d}/cellstar_integration.csv"
        gold_path = "benchmark/gold_results/cellstar_integration.csv"

        def load_csv_as_set(path):
            rows = set()
            with open(path, "r", newline="", encoding="utf-8") as f:
                reader = csv.reader(f)
                header = next(reader)

                # 严格校验表头
                expected_header = ["species", "tissue", "cell type"]
                header = [h.strip().lower() for h in header]
                if header != expected_header:
                    raise ValueError(f"Header mismatch in {path}: {header}")

                for row in reader:
                    if len(row) != 3:
                        raise ValueError(f"Invalid row length in {path}: {row}")

                    species, tissue, cell_type = row
                    rows.add((
                        species.strip().lower(),
                        tissue.strip().lower(),
                        cell_type.strip().lower()
                    ))
            return rows

        pred_rows = load_csv_as_set(pred_path)
        gold_rows = load_csv_as_set(gold_path)

        return pred_rows == gold_rows

    except Exception as e:
        print(f"Error evaluating cellstar_integration: {e}")
        return False
