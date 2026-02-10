import os


def plabdab_annotate_1(d) -> bool:
    output_csv_path = f"{d}/plabdab_annotate_1.csv"
    return os.path.exists(output_csv_path)