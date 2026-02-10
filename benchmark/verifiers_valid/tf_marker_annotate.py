import os

def tf_marker_annotate(d):
    try:
        # check if generated file exists
        return os.path.exists(f'{d}/tf_marker_annotate.json')
    except Exception as e:
        print(f"Error processing tf_marker_annotate.json: {e}")
        return False