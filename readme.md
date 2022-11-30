# 矩阵分析与应用

### 文件说明

- /config：存放预设的配置文件，提供不通过输入参数直接加载配置文件的选项
- /data：存放测试数据，包括矩阵和用于测试方程组求解的向量
- /output：存放已经执行的程序结果
- /src：存放功能实现源代码
  - utils.py：包含从文件中读取矩阵和向量的函数实现
  - Matrix.py：包含矩阵分解以及行列式、Ax=b求解函数实现
  - main.py：主函数
  - main_config.py：通过配置文件执行的主函数

- requiremens.txt：程序的依赖包列表
- Readme.md：本文件

### 依赖安装

在项目主目录下执行`pip install -r requirements.txt`即可，依赖如下：

> numpy
>
> argparse
>
> configparser

### 使用方法

#### 1.通过命令行参数

​				在src目录下，执行`python main.py --m=[分解种类] --input=[矩阵文件目录] --solve=[是否求解Ax=b] --vector=[向量文件目录]`即可。示例如下：

```shell
python main.py --m=LU --input=../data/LU.txt --solve=1 --vector=../data/V_LU.txt
python main.py --m=GS --input=../data/GS.txt --solve=1 --vector=../data/V_GS.txt
python main.py --m=H --input=../data/H.txt --solve=1 --vector=../data/V_H.txt
python main.py --m=G --input=../data/G.txt --solve=1 --vector=../data/V_G.txt
python main.py --m=URV --input=../data/URV.txt --solve=1 --vector=../data/V_URV.txt
```

​				上述五个命令可以完成五种分解的测试

#### 2.通过加载配置文件

​				为使用方便，附带了通过加载配置文件进行测试的方法。在`src`目录下执行`python main_config.py`即可完成五种分解的测试，具体内容与按序执行方法一中五个命令相同。