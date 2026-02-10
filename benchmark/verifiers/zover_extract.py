import csv

def zover_extract(d):
    try:
        pred_path = f'{d}/zover_extract.csv'
        gold_path = 'benchmark/gold_results/zover_extract.csv'

        def load_and_normalize_csv(path):
            """
            Load CSV and normalize it:
            - lower case
            - strip whitespace
            - treat each row as a frozenset of (column, value)
            - final result is a set of rows, order-insensitive
            """
            with open(path, newline='', encoding='utf-8') as f:
                reader = csv.DictReader(f)
                rows = set()
                for row in reader:
                    normalized_row = frozenset(
                        (k.strip().lower(), (v or '').strip().lower())
                        for k, v in row.items()
                    )
                    rows.add(normalized_row)
                return rows

        pred_rows = load_and_normalize_csv(pred_path)
        gold_rows = load_and_normalize_csv(gold_path)

        return pred_rows == gold_rows

    except Exception as e:
        print(f"Error processing CSV files: {e}")
        return False
