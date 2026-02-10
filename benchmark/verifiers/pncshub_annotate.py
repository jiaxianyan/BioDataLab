import os

def pncshub_annotate(d):
    try:
        pred_path = os.path.join(d, "pncshub_annotate.txt")
        gold_path = "benchmark/gold_results/pncshub_annotate.txt"

        def load_mechanisms(txt_path):
            with open(txt_path, 'r') as f:
                mechanisms = f.read().strip()
            return mechanisms

        pred_mechanisms = load_mechanisms(pred_path)
        gold_mechanisms = load_mechanisms(gold_path)

        return sorted(pred_mechanisms) == sorted(gold_mechanisms)

    except Exception as e:
        print(f"Error processing pncshub_annotate text files: {e}")
        return False
