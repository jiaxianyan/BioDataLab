### 更新./assistant/tool_biodatalab:
1. database.py 复用Biomni，数据库查询工具
2. database_tools.py 数据库下载工具和pubmed query，输入为database.py的query结果
3. structure_tools.py structure中用到的各种需要安装的工具，qiime封装了四个工具，工具默认路径可能需要安装完成后更改
   
### 更新./operation_env/tool_lake
提供了12个工具的安装脚本，其中annovar只能去官网下载zip包，网址已在bash脚本中给出，这些工具的python封装已在structure_tools.py中实现

### 更新./data/biodatalab_data/benchmark/tasks
十个structure benchmark任务，任务已经定义好，任务对应的工具也已经提供完成。输入数据，refdata和eval函数是占位符，数据体量太大

### 更新./tool_agent
自动化生成python代码和description的脚本，输入配置在tool_builder.yml中，运行时cd到此子文件夹直接运行./tool_generate即可

### 8.22 更新./tool_agent/tool_generate_wrap.py
将tool生成自动化脚本封装为每个外部工具整体封装为一个tool函数，运行方式仍然是更改tool_builder.yml中的配置，直接在此子文件夹下运行python3 ./tool_generate_wrap，增加bash脚本生成：在dependencies.md中要求llm给出安装对应工具的bash脚本


### 8.22 更新./operation_env/tool_lake
1. cactus, qiime, ncbi_datasets, exonerate 原下载链接失效，更改为conda安装，bash脚本中无内容，仅有conda安装的命令，直接通过命令在conda中安装即可
2. InterProScan 依赖java环境，可参考 https://www.jianshu.com/p/1e6da8a5ead3
3. annovar wget链接不稳定，可能需要从官网手动下载安装包并解压
4. 其余脚本直接运行，已确认无问题
5. 解压完成后annovar, arts, bowtie2, fastqc, InterProScan, mash, sra_toolkit, Trimmomatic 需要在代码中更新默认地址参数（也可以在说明文档中标明让agent自己传入此参数）