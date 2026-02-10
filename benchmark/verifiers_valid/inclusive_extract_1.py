import os


def inclusive_extract_1(d) -> bool:
    print('Checking if generated files exist...')
    
    result_1 = os.path.exists(f"{d}/inclusive_extract_1/paper_data_2126_2127.csv")
    print('Retrieve ncAAs metadata from Journal of Translational Medicine...')
    print('True/False:', result_1)
    print('-' * 50)

    result_2 = os.path.exists(f"{d}/inclusive_extract_1/paper_data_1748.csv")
    print('Retrieve ncAAs metadata from ACS Chemical Biology....')
    print('True/False:', result_2)
    print('-' * 50)

    result_3 = os.path.exists(f"{d}/inclusive_extract_1/paper_data_947.csv")
    print('Retrieve ncAAs metadata from PNAS....')
    print('True/False:', result_3)
    print('-' * 50)

    result_4 = os.path.exists(f"{d}/inclusive_extract_1/paper_data_966_1011.csv")
    print('Retrieve ncAAs metadata from Nature Chemical Biology and its SI....')
    print('True/False:', result_4)
    print('-' * 50)

    result_5 = os.path.exists(f"{d}/inclusive_extract_1/paper_data_17_21.csv")
    print('Retrieve ncAAs metadata from Science....')
    print('True/False:', result_5)
    print('-' * 50)

    print('Evaluation is done!')
    return all([result_1, result_2, result_3, result_4, result_5])