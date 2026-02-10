import csv

def amdb_extract(d):
    try:
        pred_path = f'{d}/amdb_extract.csv'
        gold_path = 'benchmark/gold_results/amdb_extract.csv'

        def load_csv_as_set(path):
            with open(path, newline='', encoding='utf-8') as f:
                reader = csv.DictReader(f)
                rows = []
                for row in reader:
                    # normalize: lowercase + strip spaces
                    normalized_row = {
                        k.strip().lower(): (v.strip().lower() if v is not None else "")
                        for k, v in row.items()
                    }
                    rows.append(normalized_row)

                # convert each row dict into a frozenset of (key, value)
                return set(
                    frozenset(item.items())
                    for item in rows
                )

        pred_set = load_csv_as_set(pred_path)
        gold_set = load_csv_as_set(gold_path)

        return pred_set == gold_set

    except Exception as e:
        print(f"Error processing CSV files: {e}")
        return False
