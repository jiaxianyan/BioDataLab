import json

def eval_list(file1_path: str, file2_path: str, order_matters: bool = False) -> bool:
    """
    Reads two JSON files, extracts the lists, and compares them.

    Args:
        file1_path (str): The path to the first JSON file.
        file2_path (str): The path to the second JSON file.
        order_matters (bool): If True, performs a direct comparison where element
                              order is important. If False, compares the lists
                              irrespective of element order. Defaults to False.

    Returns:
        bool: True if the lists are considered the same, False otherwise.
    """
    try:
        # Open and load the first JSON file
        with open(file1_path, 'r') as f1:
            list1 = json.load(f1)

        # Open and load the second JSON file
        with open(file2_path, 'r') as f2:
            list2 = json.load(f2)

    except FileNotFoundError as e:
        print(f"Error: File not found. Details: {e}")
        return False
    except json.JSONDecodeError as e:
        print(f"Error: Could not decode JSON. Make sure files are valid. Details: {e}")
        return False
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return False

    # Ensure that the loaded data from both files are actually lists
    if not isinstance(list1, list) or not isinstance(list2, list):
        print("Error: One or both of the JSON files do not contain a list at the top level.")
        return False

    # --- Comparison Logic ---
    if order_matters:
        # Direct comparison: checks length, elements, and order.
        # This is the fastest and simplest comparison.
        # Example: [1, 2, 3] == [1, 2, 3] -> True
        # Example: [1, 2, 3] == [3, 2, 1] -> False
        return list1 == list2
    else:
        # Order-insensitive comparison: checks if they have the same elements,
        # regardless of their order.
        # First, a quick check on length. If lengths differ, they can't be the same.
        if len(list1) != len(list2):
            return False
            
        # We sort both lists and then compare.
        # NOTE: This will fail with a TypeError if the lists contain complex,
        # non-sortable objects like dictionaries.
        # Example: [1, 2, 3] vs [3, 1, 2] -> sorted versions are both [1, 2, 3] -> True
        try:
            return sorted(list1) == sorted(list2)
        except TypeError:
            print("Error: The lists contain non-sortable items (like dictionaries).")
            print("Cannot perform an order-insensitive comparison using the default sort method.")
            # For complex cases like lists of dictionaries, a more advanced comparison
            # would be needed. For this function, we'll consider them not the same.
            return False
