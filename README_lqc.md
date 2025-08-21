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