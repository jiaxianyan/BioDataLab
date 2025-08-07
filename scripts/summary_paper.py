import os 
from deep_data_research.commons.utils import pmap_multi
from openai import OpenAI
from langgraph.prebuilt import create_react_agent
from langchain.chat_models import init_chat_model

def translate_english_into_chinese(text):
    client = OpenAI(base_url=BASE_URL, api_key=API_SECRET_KEY)
    response = client.chat.completions.create(messages=[{"role": "user", "content": 'Please translate the following english text into chinese:\n' + text}], 
                                                        model='gemini-2.0-flash', temperature=0)
    output = response.choices[0].message.content
    return output

def summary_chinese_content(text):
    client = OpenAI(base_url=BASE_URL, api_key=API_SECRET_KEY)
    response = client.chat.completions.create(messages=[{"role": "user", "content": '下面是一个关于数据集的论文，请你帮我总结论文，以下面的形式\n任务目标：<构建这个数据集的目的是什么，为什么要构建这个数据集>\n详细解释：<详细解释任务目标，帮助我们更加的理解，包括对一些会出现的科学概念的解释>\n构建过程：<具体的构建步骤>，目的是什么，输入是什么，输出是什么，如何实现的(用的什么工具，什么操作)>\n下面是论文具体的内容：' + text}], 
                                                        model='gemini-2.5-flash-lite-preview-06-17', temperature=0)
    output = response.choices[0].message.content
    return output

API_SECRET_KEY = "sk-di6EQiNDxyGzzoFa93389fA755B14eAfAf31B863F2E6C369"
BASE_URL = "https://aihubmix.com/v1"
os.environ["OPENAI_API_KEY"] = API_SECRET_KEY
os.environ["OPENAI_BASE_URL"] = BASE_URL
os.environ["GOOGLE_API_KEY"] = API_SECRET_KEY
os.environ["GOOGLE_BASE_URL"] = BASE_URL

root_dir = 'MinerU_output/NAR_cetergory'
dirs = os.listdir(root_dir)

for d in dirs:
    data_dir = os.path.join(root_dir, d)
    pdf_dirs = os.listdir(data_dir)

    texts = []
    for pdf_dir in pdf_dirs:
        text_md = os.path.join(data_dir, pdf_dir, 'auto', f'{pdf_dir}.md')
        with open(text_md, 'r') as f:
            text = f.read()
        texts.append(text)

    outputs = pmap_multi(translate_english_into_chinese, zip(texts), n_jobs=256, desc='translating chinese to english')

    for pdf_dir, output in zip(pdf_dirs, outputs):
        text_cn_md = os.path.join(data_dir, pdf_dir, 'auto', f'{pdf_dir}_cn.md')
        with open(text_cn_md, 'w') as f:
            f.write(output)

    texts = []
    for pdf_dir in pdf_dirs:
        text_md = os.path.join(data_dir, pdf_dir, 'auto', f'{pdf_dir}_cn.md')
        with open(text_md, 'r') as f:
            text = f.read()
        texts.append(text)

    outputs = pmap_multi(summary_chinese_content, zip(texts), n_jobs=256, desc='summary chinese content')

    for pdf_dir, output in zip(pdf_dirs, outputs):
        text_cn_md = os.path.join(data_dir, pdf_dir, 'auto', f'{pdf_dir}_cn_summary.md')
        with open(text_cn_md, 'w') as f:
            f.write(output)