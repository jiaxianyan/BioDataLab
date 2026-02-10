def cancermirnome_annotate(d):
    try:
        # Load predicted result
        with open(f'{d}/cancermirnome_annotate.txt', 'r') as f:
            pred_list = [line.strip() for line in f if line.strip()]

        # Load ground truth
        with open('benchmark/gold_results/cancermirnome_annotate.txt', 'r') as f:
            gold_list = [line.strip() for line in f if line.strip()]

        # Compare lists: case-sensitive, order-insensitive
        return sorted(pred_list) == sorted(gold_list)

    except Exception as e:
        print(f"Error processing cancermirnome_annotate results: {e}")
        return False
