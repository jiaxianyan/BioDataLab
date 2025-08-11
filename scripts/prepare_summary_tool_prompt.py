import os 
names = os.listdir('MinerU_output/NAR_cetergory')
with open('prompts/help_summary_one_paper_database_and_tools.txt', 'r') as f:
    prompt = f.read()
for name in names:
    for d in os.listdir(f'MinerU_output/NAR_cetergory/{name}'):
        with open(f'MinerU_output/NAR_cetergory/{name}/{d}/auto/{d}_cn_summary_v2.md', 'w') as f:
            f.write('\n')

        with open(f'MinerU_output/NAR_cetergory/{name}/{d}/auto/{d}_cn.md', 'r') as f:
            paper_cn = f.read()

        paper_prompt = prompt + f"{paper_cn}\n"

        with open(f'MinerU_output/NAR_cetergory/{name}/{d}/auto/{d}_cn_summary_v2_input_prompt.md', 'w') as f:
            f.write(paper_prompt)