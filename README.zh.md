# hetfamily — 异质性处理效应（jamovi 模块）

在随机实验中找出谁获益更多：通用机器学习推断（Chernozhukov、Demirer、Duflo、Fernandez-Val）、优先级曲线、效应轮廓和处理规则。 英文完整说明见 [README.md](README.md)。

## 分析（jamovi 菜单：Heterogeneity）

| 分析 | 内容 |
|---|---|
| **异质性效应（BLP、GATES、CLAN）** | 最佳线性预测、按预测效应排序的分组平均效应、最受影响与最不受影响组的协变量均值 |
| **优先级曲线（RATE）** | TOC 曲线、AUTOC 与 Qini，bootstrap 标准误 |
| **效应轮廓（CATE 回归）** | 双重稳健得分对调节变量回归：平均效应及每个标准差的变化 |
| **处理规则** | 预测效应超过成本时才处理：相对不处理和全部处理的收益 |

每个分析的输出部分都有"方法说明"和"数据与设置指南"。

学习器：lasso、岭、弹性网、随机森林、线性回归；多次随机分样，用中位数合并。处理必须随机分配（或倾向得分已知）。第二阶段与 `GenericML`、`grf` 对照一致（`tests/reference.R`、`tests/profile.R`、`tests/rate.R`）。示例数据：`sim_het_trial.csv`、`sim_homogeneous_trial.csv`。

## 安装
- **jamovi 模块库**：在 jamovi 右上角打开模块按钮，进入 *jamovi library*，搜索并安装（模块上架后）。
- **侧载**：下载 `.jmo` 文件，在同一菜单选 *Install from file*。`.jmo` 只适用于构建它的操作系统、处理器类型和 jamovi 版本。
- **从源码**：需要 R 和 `jmvtools`：`jmvtools::install(home = "C:/Program Files/jamovi 28.17.0.0")`（路径按本机 jamovi 修改）。

## 引用与许可
引用：Sun, M. *hetfamily* (jamovi module). https://github.com/mingchao-sun/hetfamily 。请同时引用所用方法的原始文献（输出中列出）。

许可：见 `DESCRIPTION`（GPL）；第三方代码与数据见 `NOTICE`。问题与建议：https://github.com/mingchao-sun/hetfamily/issues
