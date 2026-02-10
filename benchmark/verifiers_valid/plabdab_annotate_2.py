import os


def plabdab_annotate_2(d) -> bool:
    output_csv_path = f"{d}/plabdab_annotate_2.csv"
    return os.path.exists(output_csv_path)