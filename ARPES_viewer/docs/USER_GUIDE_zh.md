# ARPES viewer 用户指南

*[English version / 英文版](USER_GUIDE_en.md) · [Version française / 法文版](USER_GUIDE_fr.md)*

*最后更新：2026 年 10 月 25 日——对应 ARPES viewer v1.0.0*

程序的版本号和更新日期显示在主窗口标题栏和主窗口最下方。报告 bug、提问或获取最新版本：
<https://github.com/DingPei1995/ARPES-data-browser>（报告 bug 时请注明版本号）。

本指南按照日常使用的顺序，配合截图介绍程序的主要功能：载入数据、浏览数据、
转换到动量空间、分析，以及制作论文图。每一步背后的技术细节（公式、测试结果、
设计取舍）见 [README](../README.md) 和 [CHANGELOG](../CHANGELOG.md)。

所有截图都来自实际运行的程序，所用数据是 `devtools/demo_data.py` 生成的模拟
数据，因此没有实验数据也可以跟着操作一遍（见[用演示数据试用](#2-用演示数据试用)）。

---

## 目录

1. [安装与启动](#1-安装与启动)
2. [用演示数据试用](#2-用演示数据试用)
3. [主窗口](#3-主窗口)
4. [载入数据](#4-载入数据)
5. [各查看窗口的共同操作](#5-各查看窗口的共同操作)
6. [Map：等能面及其两个切面](#6-map等能面及其两个切面)
7. [Cut：从原始谱到能带参数](#7-cut从原始谱到能带参数)
8. [Map 的角度→动量转换与布里渊区](#8-map-的角度动量转换与布里渊区)
9. [任意路径切面与费米面校正](#9-任意路径切面与费米面校正)
10. [数据处理：平滑、求导、曲率](#10-数据处理平滑求导曲率)
11. [堆叠图与三维视图](#11-堆叠图与三维视图)
12. [Cut 运算：二色性与参考谱归一](#12-cut-运算二色性与参考谱归一)
13. [光子能量扫描（kz）](#13-光子能量扫描kz)
14. [去除探测器网格](#14-去除探测器网格)
15. [曲线：EDC、MDC 与自旋分辨数据](#15-曲线edcmdc-与自旋分辨数据)
16. [空间分辨扫描（nano-ARPES）](#16-空间分辨扫描nano-arpes)
17. [作图与导出](#17-作图与导出)
18. [自动保存、内存与恢复](#18-自动保存内存与恢复)
19. [技巧与常见问题](#19-技巧与常见问题)

---

## 1. 安装与启动

需要 Python 3.9 或更高版本，以及五个依赖包：

```bash
pip install PyQt5 pyqtgraph numpy h5py scipy
```

然后在 `ARPES_viewer` 目录下运行：

```bash
python ARPES_viewer.py
```

支持读取的数据：

| 线站 / 来源 | 文件 |
|---|---|
| SOLEIL **ANTARES** | NeXus/HDF5 `*.nxs`（cut、偏转器 map、空间扫描） |
| SOLEIL **CASSIOPEE** | Scienta SES 文本文件——单张谱，以及按编号组成的文件夹（自动拼成 map 或光子能量扫描） |
| CASSIOPEE **自旋端站**（MBS A-1） | `*.krx` 与 `*.txt`——cut、偏转器 map 和自旋分辨 EDC |
| 本程序 | 本程序自己保存的 `*.nxs` 格式 |

## 2. 用演示数据试用

```bash
python devtools/demo_data.py            # 生成 ./demo_data/*.nxs
```

会生成六个模拟数据文件：一张费米面 map、一对 LH / LV 偏振下的 cut、一个光子能量
扫描、一条自旋分辨 EDC、一个小的实空间扫描，以及一张叠加了探测器网格的 map。
它们来自一个简单的能带模型（Γ 点的六角电子口袋与空穴带、K 点的空穴口袋），
包含费米截断、分析器分辨率与泊松噪声，打开方式与实测数据完全相同。

修改窗口之后，如需重新生成本指南的全部截图：

```bash
python devtools/guide_screenshots.py    # 写入 docs/images/*.png
```

---

## 3. 主窗口

![主窗口](images/main_window_info.png)

主窗口是一个**数据浏览器**：每行对应一个*数据集*，而不是一个文件（一个 `.nxs`
文件常常包含多个测量，每个测量单独一行）。名字后面标明数据类型：`[Cut]`、
`[Map]`、`[kz map]`、`[k-map]`、`[SPEM]`、`[EDC]`、`[Spin EDC]` 等。

| 操作 | 方法 |
|---|---|
| 在独立窗口中打开数据集 | **双击**该行（或右键 → **Open**） |
| 把多张 cut 叠成一个 map（温度、门电压、光子能量 → kz map 等），或把分段采集的 map / kz map 拼接起来 | 选中它们 → **Combine datasets...**（或右键 → **Combine N datasets...**） |
| 其余所有针对数据集的操作 | **右键**该行（或选中的多行） |
| 查看各轴的首末值、点数、步长 | 右键 → **Show information**（填入 *Data Information* 表） |
| 改变数据的网格或标度——裁剪范围、合并像素、每张谱自归一——可多个数据集同时进行 | **Crop / bin / normalise...** |
| 突出能带的滤波——平滑、二阶导数、曲率、扣背景、对称化 | **Smooth / derivative...** |
| 查看并释放内存 | 底部的 **Memory** 区 |
| 查看当前版本、报告 bug、获取最新版本 | 标题栏，以及窗口最后一行（GitHub 链接） |

列表下方的三个按钮位于 **With the selected datasets**（对选中的数据集操作）框中：
它们都作用于上方选中的行，结果都作为新行加入列表——原数据不会被修改。

### 各名称对应的功能

v1.0.0 中若干按钮改了名，使名称直接说明它做什么。新旧名称对照：

| 现在 | 以前 | 作用 |
|---|---|---|
| **Combine datasets...** | Combine data... | 把 cut 叠成 map，或把分段采集的数据拼接起来 |
| **Crop / bin / normalise...** | Data operations... | 改变网格或标度：裁剪（truncate）、合并像素（compress）、自归一 |
| **Smooth / derivative...** | Process... | 突出能带的滤波：平滑、求导、曲率、扣背景、对称化 |
| Functions 菜单中的 **Transform** 分组 | Data operations | 查看窗口中的转换与校正 |
| **Fermi edge fit...** | Fermi level... | 在 cut 上拟合 E_F（及分辨率） |
| **Convert to momentum (k)...** | Cut / Map k conversion... | 角度 → k∥ |
| **Align Fermi levels (kz map)...** | kz map processing... | 拟合并对齐每个光子能量的 E_F |
| **Convert to momentum (kz, k∥)...** | kz -> momentum... | 光子能量扫描 → (k_z, k∥) |
| **Straighten bands (FS correction)...** | FS correction... | 把弯曲的费米边或特征拉直 |
| **Remove detector grid...** | De-grid... / De-grid map... | 去除 MCP / 网栅图案 |
| **Multi-slice figure...** | Slice figure... | 由 map 的一系列切片组成的一页图 |
| **Freeze slice in a new window** | Open slice in a new panel | 当前切片的冻结副本 |
| **Plot tools**（查看窗口顶部一行） | —（新增） | 在实时图上画辅助线和文字 |
| **Figure composer...**（冻结切片窗口）、**Open in figure composer...**（列表右键） | Plot tools、Plot as a figure... | 制作论文图 |
| **Figure tools...**（作图器内） | Plot tools | 该论文图的尺寸、字体、刻度、色标、标注等 |

右键菜单只列出适用于当前选择的操作。灰色的条目在鼠标悬停时会说明原因——
例如 *Cut arithmetic* 需要恰好两个 cut，*3D view* 需要 map。

**Show information** 还会在独立窗口中显示该数据集的全部元数据，带过滤框，
并可 **Export metadata (.csv)**：

![数据集信息](images/information_window.png)

**Clear list** 会先确认；如果本次会话中有计算得到但尚未保存的数据，会提示保存。
磁盘上的文件永远不会被删除。

### 合并数据集：叠加（Stack）与拼接（Join）

**Combine datasets...**（列表下方，或右键 → **Combine N datasets...**）可以做两件
不同的事，并会让你选择：

- **Stack（叠加）**——*仅限 cut*。每张 cut 成为新 map 的一个切片，沿一个由你命名的
  轴：温度序列、门电压序列、每个光子能量一张的 cut。结果可以像普通 map 一样查看——
  在固定能量处切片，或沿新轴切片。
- **Join（拼接）**——*cut、map、kz map 或 k-map，但须同一类型*。把一次分段采集的
  测量按其自身的某个轴首尾拼接起来：偏转角范围分两次扫描的 map、光子能量范围分几次
  采集的 kz map、能量窗口分两半采集的 cut。

**叠加 cut**

1. 在列表中选中这些 cut（Ctrl/Shift 多选），点 **Combine datasets...**。叠加要求
   各 cut 的坐标轴逐点一致（点数、数值、单位都相同）；若不一致，*Stack* 选项变灰并
   指出第一处差异，此时可改用 *Join*。
2. 程序会询问是否**按列表中的顺序**合并。选 *Yes* 保持该顺序；选 *No* 则可自行
   调整顺序。
3. 选择 **Stack**，然后设置新轴：**New axis is**（温度、门电压、光子能量、延迟、
   位置、其他——会自动填入 *Name* 和 *Unit*，均可修改）、**Start**（起始值）和
   **Step**（间隔）。每张 cut 旁边会显示它对应的数值。如果选择了调整顺序，可
   **拖动**各行（或用 **Move up** / **Move down** / **Reverse**）。点 **Combine**。

![合并 cut](images/combine_cuts_dialog.png)

由角度 cut 组成的光子能量序列会成为 `kz map`，可直接使用 *Align Fermi levels (kz
map)* 和 *Convert to momentum (kz, k∥)*（见[第 13 节](#13-光子能量扫描kz)）；
其他轴得到的是 `map`，其第一个轴被标记为非角度，因此不会提供 k 转换。

**拼接分段采集的数据**

1. 选中各段（全是 map，或全是 kz map，或全是 cut……），点 **Combine datasets...**。
   map 和 kz map 直接进入 *Join*。
2. **Join along**：测量被分段的那个轴——map 默认是扫描（偏转角）轴，kz map 默认是
   光子能量；对 cut 可选择角度轴或能量轴。
3. **Where pieces overlap**（重叠部分）：*Average them*（取平均），或按下方列表的
   顺序保留较早 / 较晚的一段（拖动各行可改变顺序）。
4. 对话框列出每一段的范围，并给出 **Result**（拼接后的范围和点数），以及处理过程中
   所做的说明。点 **Combine**。

![拼接分两段采集的 map](images/combine_dialog.png)

坐标轴的处理方式：

- 各段在同一网格上（步长相同，范围相接或按整点重叠——通常情况）时，数据**按原样
  复制**，不做插值。
- 某段步长不同或有半步偏移时，会被**插值**到公共轴上（采用第一段的步长），对话框
  会说明。
- 若**其他轴**不同——几次采集的能量窗口相差几 meV、slit 通道数不同——所有段都会
  重采样到它们**共同覆盖**的范围，对话框会说明。
- 没有任何一段覆盖的区间是**空缺**，填为零并给出提示。
- 拼接后的 kz map 仍是 `kz map`：请在拼接结果上运行 *Align Fermi levels (kz
  map)*，而不是在各段上（各段按光子能量的拟合结果不会被带入）。

---

## 4. 载入数据

**Load data...** 打开载入窗口；**Quick add...** 跳过该窗口，直接以自动识别的
设置加入所选文件。

![载入窗口](images/loader_dialog.png)

1. **Files**——选择一个或多个文件（*Select file(s)...*）。文件类型过滤器由已安装的
   读取器自动生成。
2. **Reader**——自动识别，并显示识别结果（如 *Detected: 3 × This program's own
   format*）。只有在文件被误读时才需要手动指定。
3. **Datasets found**——每个文件里有哪些数据集。这里只读取文件*结构*，即使是几 GB
   的扫描也是瞬间完成；在真正打开某一行之前不会读取数据本身。
4. **Axes**——*Axis order* 在载入时一次性调整数组维度顺序；*First axis of a map
   is* 说明扫描的是什么：偏转器角度（默认）、**光子能量**、温度、栅压、延时……
   这样后续每一步都使用正确的单位，也不会对非角度轴提供 k 转换。

**CASSIOPEE 文件夹**：选中编号序列（`name_1_ROI1_.txt`、`name_2_ROI1_.txt`……）
中的任意一个文件，整个文件夹会作为一个拼好的数据立方体列出一次。它是角度 map
还是 kz map，由谱旁边的 `_i` 参数文件自动判断。列表中的名称带有从这些文件读出的
光子能量和偏振——如 `FineMap 117eV LH`，光子能量扫描则为 `KzScan 40-120eV LH`——
因此同一扫描在不同能量或偏振下的重复测量可以区分开。

**ANTARES 光子能量（PhD）文件**：PhD `.nxs` 文件（每个光子能量一个 `Phd_NNNN`
条目）作为一个 `kz map` 列出：各帧按光子能量排序，每帧的动能轴用其自身的 hν
换算到 E_F（功函数取 4.5 eV，k_z 转换也会用到）。每一帧仍可单独作为 cut 打开。

### 数据类型

| 类型 | 轴 | 含义 |
|---|---|---|
| `cut` | 角度（或 k）× E | 一张谱 |
| `map` | 角度 × 角度 × E | 偏转器或极角扫描 |
| `k_map` | kx × ky × E | 已转换到动量空间的 map |
| `kz_map` | hν × 角度 × E | 光子能量扫描 |
| `kz_map_k` | k_z × k∥ × E | 转换到动量空间的光子能量扫描 |
| `spem_4d`、`spem_1d` | x（× y）× 角度 × E | 实空间扫描 |
| `edc`、`mdc`、`spin_edc` | E（或 k）× 通道 | 曲线，可带 σ 误差通道 |

---

## 5. 各查看窗口的共同操作

每个数据集在自己的窗口中打开，因此可以把 map、它的切面以及另一个文件并排摆放。

![带读数光标的 cut 窗口](images/cut_viewer.png)

- **顶部一行**：**Functions** 菜单（本窗口能做的所有事）、**Plot tools**（在本图上
  画辅助线和文字，见下文）、色表（每个窗口独立）及 **Flip**（反转），以及视图控制——坐标范围、**Auto**、**Inv**（反向）、
  **Reset**、**Grid**。
- **图像上方**：**Levels** 色阶条（拖动两个滑块，或输入 **Min%/Max%** 百分位，
  **Clip** 应用）以及 **Gamma**。您设置的最小/最大值在移动切片位置时（contour 的
  能量、cut 的角度、空间 map 的光标）会保持不变，使各位置处于同一色标；把滑块拖回
  两端，或设为 0%–100%，即恢复为跟随每一帧自身的范围。
- **图像下方**：**Smooth**（仅影响显示，不改变数据）和 **Lock x/y**（锁定纵横比）。
- **读数光标**——在图像上右键 → *Readout cursor*。它显示光标处的数值；在 cut 上
  还会画出经过该点的 **EDC**（右侧）和 **MDC**（下方）。**EDC ±** / **MDC ±**
  设置积分窗口，并在图上以半透明带显示。在 **Cursor** 行输入坐标并回车，可把
  光标精确放到该位置。
- **EDC → list** / **MDC → list** 把当前显示的曲线作为独立数据集加入主列表
  （见[第 15 节](#15-曲线edcmdc-与自旋分辨数据)）。
- **图像右键 → Export** 导出该面板当前显示的内容——数据（`.nxs`）、图片，或
  矢量坐标轴（`.eps`），见[第 17 节](#17-作图与导出)。
- **Functions → Slice → Save slice to the main list / Freeze slice in a new
  window** 保存当前切片的副本。冻结窗口（以及空间 map 上的 *Open this position in
  a new window*）会沿用原面板的最小/最大值和 gamma，便于直接对比。它有自己的
  **Plot tools**，以及用来制作论文图的 **Figure composer...** 按钮
  （见[第 17 节](#17-作图与导出)）。

**Functions** 菜单按类别分组：

| 窗口 | Analysis（分析） | Transform（转换与校正） | Visualization（可视化） |
|---|---|---|---|
| map / kz map / k-map | | Arbitrary cut、Convert to momentum (k)、Align Fermi levels (kz map)、Convert to momentum (kz, k∥)、Remove detector grid | Brillouin zone（k-map）、Multi-slice figure |
| cut | Fermi edge fit、MDC / EDC fit | Straighten bands (FS correction)、Convert to momentum (k)、Cut arithmetic、Remove detector grid | |
| 曲线 | Curve fit、Spin analysis | Crop、Bin、Normalise、Background、Shift axis、Add √N errors | Figure composer |

不适用于当前数据的条目不会显示；适用但暂时无法运行的条目显示为灰色，
悬停提示会说明原因。

### Plot tools：在图上画辅助线和文字

每个查看窗口（cut、map、空间扫描、曲线、冻结切片）顶部一行都有 **Plot tools**
按钮，打开一个小窗口，**直接在眼前这幅图上**作画：

![cut 上的辅助线](images/plot_tools_viewer.png)

| ![Plot tools](images/plot_tools.png) | |
|---|---|

- **Type**（类型）：*Vertical line*（竖线）、*Horizontal line*（横线）、*Line at an
  angle*（过某点的斜线；角度按图本身的单位计，即斜率为 tan(角度) 个 y 单位 / x 单位）、
  *Text label*（文字）。
- **Position**（x, y）、**Angle**（角度）、可选的 **Label**（显示在线上，如 `E_F`、
  `Γ`）、**Colour**（九种一键颜色，或 **More...** 自选）、**Width**（线宽）、
  **Style**（实线、虚线、点线、点划线）和 **Text size**（字号）。
- **Add** 在输入的位置画出（**Centre of view** 可填入当前视图中心）；**Place by
  clicking** 则在图上点击处画出。
- 在图上**拖动**线条即可移动，窗口中的数值同步更新。在列表中选中某条可**编辑**
  （修改立即生效）或 **Remove**（删除）；**New guide** 结束编辑；**Remove all**
  清空。
- 有两个面板时（空间扫描），用 **Draw on** 选择画在哪一个上。
- 移动切片时辅助线保持不变（浏览 map 时 E_F 线始终在 E_F 处），关闭 Plot tools
  窗口后也保留。面板的 **Export → Image** 也会画出它们（*Include the Plot tools
  guides*）。

*Plot tools* 用于看数据；要制作论文图，请使用**作图器（Figure composer）**及其
**Figure tools**（见[第 17 节](#17-作图与导出)）。

---

## 6. Map：等能面及其两个切面

map 打开时显示**等能面（constant-energy contour）**。顶部滑块选择能量（**Pos**
以 eV 为单位，**Ind** 为序号；**±** 设置积分范围，可选 **mean** 或 **sum**）。
**Deflector cut** 和 **Slit cut** 分别在新窗口中打开两个正交切面。

![map 的等能面](images/map_contour.png)

| Slit cut（狭缝角 vs 能量） | Deflector cut（偏转器角 vs 能量） |
|---|---|
| ![Slit cut](images/map_slit_cut.png) | ![Deflector cut](images/map_deflector_cut.png) |

**三个窗口共用一个光标。** 在任一窗口中打开读数光标，其余窗口同时打开；在任一
窗口中移动光标，其余窗口会沿新位置重新切片。每个轴在所有窗口中都用同一种颜色：
**偏转器为红色、狭缝为绿色、能量为蓝色**。光标线的颜色是它保持不变的那个轴的
颜色，它所对应的 EDC 或 MDC 也用同一颜色绘制。积分宽度（± 框）属于各自窗口，
不会在窗口之间共享。

也可以直接从列表中打开 map 的切面（右键 → **Open slit cut** / **Open deflector
cut**），此时等能面窗口保持隐藏，需要时用 *Functions → Show the map* 调出。

---

## 7. Cut：从原始谱到能带参数

单张 cut 的常规流程：**费米边拟合 → 平移能量轴 → 转换到动量 → MDC/EDC 拟合**。

### 7.1 费米边拟合

*Functions → Fermi edge fit...* 对选框内求和得到的 EDC 拟合费米边（费米-狄拉克
分布 × 线性态密度，再与高斯分辨率卷积）。

![费米能级对话框](images/fermi_dialog.png)

1. 点击 **Show the selection box**，在 cut 上把选框拖到费米边处（或直接输入角度
   与能量范围）。角度范围决定哪些通道被求和，能量范围就是拟合窗口。

   ![cut 上的拟合选框](images/fermi_cut_box.png)

2. 先 **Estimate**，再 **Fit**。勾选 **Hold** 可固定某个参数（温度默认固定为
   文件记录值）。每个参数都给出不确定度。
3. **Offset energy axis** 生成一份能量轴为 **E − E_F** 的副本并加入列表。

*Advanced* 区域可以除去费米截断（用于能隙与参考样品），或逐通道拟合 E_F。

### 7.2 把 cut 转换到动量空间

*Functions → Convert to momentum (k)...* 把角度轴转换为 k∥（Å⁻¹）。一张 cut 是 (kx, ky)
平面上的一条线，通常不经过 Γ 点，所以对话框需要知道 Γ 的位置：手动输入角度，
或用 **From a converted map...** 从已转换的 map 继承。

![Cut k 转换](images/cut_k_dialog.png)

如果该 cut 已经平移到 E − E_F，对话框会自动把扣除的 E_F 填入 **Energy offset**
（下方的提示会说明）：转换时使用真实的动能，而结果的能量轴仍保持 E − E_F。

![转换到动量空间的 cut](images/cut_in_k.png)

### 7.3 MDC / EDC 拟合

*Functions → MDC / EDC fit...* 逐行拟合峰，再从拟合出的峰位读出能带。它只接受
横轴已经是 Å⁻¹ 的 cut——以 eV·deg 为单位的费米速度没有意义。

![MDC 拟合面板](images/mdc_fit_panel.png)

1. **Bands** 页：每条能带分支点一次 **Add band**；选中一条能带，按 **Place a
   seed**，再在右侧曲线上点击峰的位置。沿能带放四五个种子点就够了；**Line at**
   切换当前查看的行。
2. **Try this line** 显示当前行的拟合结果（虚线为初始猜测）。
3. **Fit settings** 页：选择 MDC 或 EDC、峰形（Lorentzian、Gaussian、Voigt）、
   背景、步长，然后 **Fit the series**。拟合出的峰位画在图上，点的大小代表振幅。
4. **Dispersion...** 拟合能带：**v_F**（直线）、**m\***（抛物线），并检查拟合
   在多大范围内稳定。**Self-energy...** 给出 Re Σ 与 Im Σ。

![色散窗口](images/dispersion_window.png)

**Save table...** 保存所有行的拟合结果；**Band to list** 把一条能带加入主列表。

---

### 7.4 芯能级检查（Core-level check）

宽范围的芯能级扫描也以 cut 的形式列出。在主列表中右键点击它 →
**Core-level check...**，即可判断各个峰的来源。

1. 输入样品所含元素——`NbSCl` 或 `Nb, S, Cl`——然后按 **Next**。
   *Match within*（0.5 eV）是峰与手册中结合能的最大偏差，在此范围内即认为是样品信号。

   ![元素](images/core_level_elements.png)

2. 检查能量参考。若窗口内包含费米边，则拟合得到 E_F；否则取 hν − W。
   如光子能量、功函数或 E_F 有误可直接修改（绿色竖线可拖动）；若 hν − E_F
   明显偏离 4–5 eV，说明光子能量不准。按 **Analyse**。

   ![能量参考](images/core_level_reference.png)

3. 结果：以结合能为横轴的积分 EDC，每个峰都有标注——绿色 *sample*
   （标出轨道），橙色 *uncertain*（2 eV 内有样品能级：化学位移或荷电，由您判断），
   红色 *not the sample*（列出可能的常见元素：污染、高次光、俄歇峰）。
   点击表格某行可查看全部候选；在 *Assignment* 列可填写您自己的判读。

   ![芯能级检查](images/core_level_check.png)

**弱峰与手动添加峰**：勾选 **Log counts** 可在强峰旁看清弱峰；**Min height**
（默认为最强峰的 0.2%）设定可被找到的峰的最小高度。要添加搜索遗漏的峰，按下
**Click to add peak** 后在图上点击（双击始终有效），或在 *or at* 后输入结合能并点
**Add**：其位置由局部拟合精修，标记为 ◆，并列出 ±2 eV 内所有表列能级。
**Remove selected** 删除选中的峰；**Reset peaks** 恢复自动找到的峰并清除手动添加的
峰。您手动填写的判读在 **Re-analyse** 后仍会保留。

**Peak table...** 列出所有表列能级，可按元素和结合能范围筛选：查找某个峰、在图上
显示某能级，或以该能级为判读将其添加为峰。**Select range** 在图上显示一个可拖动的
区间；**Zoom window...** 把该区间（或当前视图）在独立窗口中打开，窗口带光标、光标处
的高斯拟合，以及除以强度的负二阶导数——强峰旁的肩峰或弱峰也能和强峰一样清楚。

若所有峰整体偏移相同的量，备注中会给出提示，**Apply suggested shift**
可一键修正（样品荷电或 E_F 有误）。**Add EDC to list** 将结合能标度下的 EDC
加入列表；**Copy report** 复制全部结果文本。仍有疑问时，换一个光子能量再测：
芯能级的结合能不变，俄歇峰和高次光会移动。

## 8. Map 的角度→动量转换与布里渊区

### 8.1 Map 的 k 转换

在等能面窗口中选择 *Functions → Convert to momentum (k)...*。把读数光标拖到应成为
**k = (0, 0)** 的点（对话框会跟随光标），然后点 **Convert**，结果作为新的
`k-map` 行加入列表。

![Map k 转换](images/map_k_dialog.png)

- **Sample rotation** 绕表面法线旋转 map；**Set rotation from contour...** 通过在
  等能面上点击一个方向（例如沿镜面线的两点）读出旋转角。
- **Output grid**：按点数或按分辨率（Å⁻¹、meV）设置输出网格。
- **Energy offset** 在转换前加到能量轴上，用于能量轴不是真实动能的情况
  （对 E − E_F 的 map 会自动填好）。

![转换后的 k-map](images/k_map_contour.png)

### 8.2 布里渊区

在 k-map 上：*Functions → Visualization → Brillouin zone...*。

![布里渊区对话框](images/brillouin_dialog.png)

1. 先选**空间群**（1–230）；被空间群约束的晶格参数显示为灰色并自动跟随。也可
   切换到 **2D layer** 使用二维表面晶格。
2. **Cut plane**：表面法线 (hkl) 及沿法线的偏移——画出的是三维布里渊区被该平面
   截出的截面，适用于任意晶面，而不只是 (001)。**Through point** 列出平面经过高对称点
   （菱方晶体为 Γ、Z、L、F……）时的偏移量，直接选择即可，无需自己计算。
   **Show cut plane in 3D...** 可预览，图中标出高对称点，并圈出该平面经过的点；
   **Label high-symmetry points** 会把这些点标在等能面图上，方便设置方位角。
   对话框会自动沿用上一次使用的晶体参数。
3. 选择 **Conventional**（常规）或 **Irreducible**（不可约）区；设置 **Azimuth**
   （或 **Pick direction from contour...**）；**Repeat across the field of view**
   在视野内重复绘制。
4. **Moiré** 加入第二层（扭转或晶格失配）及其莫尔布里渊区。
5. **Plot**。

![叠加布里渊区的 k-map](images/k_map_with_bz.png)

---

## 9. 任意路径切面与费米面校正

### 9.1 任意路径切面

在任意等能面上 *Functions → Arbitrary cut...*：点击（或输入）最多六个点，切面
按顺序经过这些点——这里是 k-map 上的 Γ → M → K → Γ。勾选 **One dataset per
segment** 可让每一段分别成为一个数据集。

| 等能面上的路径 | 路径点 |
|---|---|
| ![任意切面路径](images/arbitrary_cut_path.png) | ![任意切面对话框](images/arbitrary_cut_dialog.png) |

结果是以沿路径距离为横轴的 cut；路径的拐点用绿色虚线标出，避免把拐点处的折线
误认为物理特征。

![沿 Γ–M–K–Γ 的切面](images/arbitrary_cut_result.png)

### 9.2 费米面（FS）校正

在 cut（或 map 的 slit cut）上 *Functions → Straighten bands (FS correction)...*，把本应平直的特征
拉直——通常是随分析器角度弯曲的费米边。沿该特征点击若干点，程序拟合多项式
（默认 2 阶），并把每一列在能量方向上平移。**From a reference (Au)...** 改用列表
中的金参考谱确定曲线。在 map 的 slit cut 上操作时，同样的平移会应用到整个数据
立方体。

![费米面校正](images/fs_correction_dialog.png)

---

## 10. 数据处理：平滑、求导、曲率

在列表中选中一个或多个 cut（或 map），点击 **Smooth / derivative...**。每个结果都是新的
数据集，并记录生成它的处理步骤。

![二维处理](images/process_dialog.png)

- 各页：**Smooth**（高斯、方框等，可先去尖峰）、**Derivative**（求导）、
  **Curvature**（一维与二维曲率）、**Background**（背景）、**Symmetrise**（对称化）。
- **Source**（原图）与 **Result**（结果）并排显示、色阶一致，便于看出参数是否把
  能带"吃掉"了；两者同步平移和缩放。
- **Shape**：*Square pixels*（每个数据点为正方形——对角度-能量 cut 不失真）、
  *Equal axes*（两轴单位长度相同——适合角度-角度或 k-k 图，两轴单位相同时自动
  选用）、*Fill the box*（填满窗口）。
- **Live preview** 在修改参数时实时重算；**Apply** 把结果加入列表。

对 map，处理按平面进行，并首先询问**作用于哪组平面**——等能面还是色散切面，
因为等能面的曲率与 cut 的曲率是不同的物理量。**Preview at** 选择预览的平面；
**Apply to the whole cube** 对所有平面执行。

![三维处理](images/volume_process_dialog.png)

主窗口的 **Crop / bin / normalise...** 可同时改变一个或多个数据集的网格或标度——
前提是它们类型相同，且各轴点数一致。各页：**Crop (truncate)** 裁剪到某一范围、
**Normalise** 每张谱自归一、**Bin (compress)** 合并像素。

![裁剪 / 合并 / 归一](images/data_operations.png)

---

## 11. 堆叠图与三维视图

**Stack plot**（右键 cut）：EDC 或 MDC 的瀑布图，可设置曲线条数、合并、归一化、
偏移、错切与配色；**Send to figure composer** 把它作为一个图面板送入作图器。

![堆叠图](images/stack_plot.png)

**3D view**（右键 map）首先询问区域与采样密度，因为每一帧都要对数据立方体重新
采样——*Fast*、*Balanced* 或 *Fine* 决定能否跟上鼠标操作。

![三维视图：区域与采样](images/volume_range_dialog.png)

之后拖动图像即可旋转。**Slices** 在指定位置显示三个正交平面，**Notched cube**
切掉立方体的一个角，**Projection** 沿某轴求和。右侧可设置色表、色阶、gamma、
方位角、仰角和深度比例；**Save the view...** 导出图像。

![三维视图](images/volume_view.png)

---

## 12. Cut 运算：二色性与参考谱归一

把两张 cut 组合成第三张：线二色性（LH − LV）、圆二色性、除以金参考谱，或普通的
差、比、和。在列表中选中两张 cut → 右键 → **Cut arithmetic on the two...**
（先选的为 A）；或在 cut 窗口中 *Functions → Cut arithmetic...*，然后在列表中点击
另一张 cut。

![Cut 运算](images/cut_arithmetic.png)

- **顶部表格**对比两次测量记录的条件（偏振、光子能量、角度、样品位置、温度、
  通能、透镜模式）。橙色行表示差异超过容差——在样品上两个不同位置之间做的
  二色性图，反映的其实是*位置*差异。
- **Purpose** 预设：*Linear dichroism*（A − B，B 按总强度缩放）、*Circular
  dichroism*（(A − B)/(A + B)）、*Divide by a reference*（除以参考谱）、*Custom*。
  若偏振看起来颠倒了会给出警告（可用 **Swap A ↔ B** 交换）。
- **Scale B to A** 按总强度或某一区域的强度缩放；**Hide where A+B below** 屏蔽
  背景区域。
- A、B 使用同一色阶；结果用以零为中心对称的红-白-蓝色阶。若两者都是原始计数，
  会传播泊松误差，并在图下方给出统计。
- **Result to list**、**Uncertainty to list**、**All three as a figure**。

---

## 13. 光子能量扫描（kz）

光子能量扫描在列表中显示为 `kz map`；其"等能面"是某一能量下光子能量对角度的图。

![kz map](images/kz_map_contour.png)

### 13.1 对齐 kz map 的费米能级

文件中并没有记录每张谱真正的费米能级位置。*Functions → Align Fermi levels (kz map)...*
直接从谱本身测量它。

1. slit cut 随之打开，上面有一个选框。把选框拖到费米边处——角度范围决定求和的
   通道，能量范围是拟合窗口。

   ![slit cut 上的选框](images/kz_processing_cut.png)

2. **Fit and assemble**：对每张谱使用同一个选框（按探测器像素序号），把每张谱
   平移使其自身的费米边位于零点，把数据立方体裁剪到所有谱共有的能量范围，并
   （默认）**按每张谱的总强度归一**，以消除光通量的变化。

![对齐费米能级 (kz map)](images/kz_processing_dialog.png)

图中给出 E_F 随光子能量的变化。曲线平滑说明拟合可靠；拟合不可信的点会被标记，
并用相邻点插值补上。

### 13.2 把 kz map 转换到动量空间

*Functions → Convert to momentum (kz, k∥)...* 把扫描转换为 (k_z, k∥)。

![kz → 动量](images/kz_conversion_dialog.png)

- **内势 V0** 是未知量：修改它时预览会在毫秒内重新转换一个能量切片。
- 给定**晶格**（空间群、a、b、c）、**表面法线**和面内方向后，预览上会画出布里渊
  区边界。
- **Scan V0** 在一系列 V0 下转换，报告 k_z 周期与晶格匹配时的 V0，以及 V0 能被
  确定到什么精度。
- **Measure a k_z period...**：点击相隔一个布里渊区的同一特征上的两点，程序列出
  与之匹配的倒格子晶面。
- 黄色横幅提示费米边仍随分析器角度弯曲——应先在 slit cut 上做 Straighten bands (FS correction)。
- **Convert the whole cube** 生成 `kz map (k)`。

![转换后的 kz map](images/kz_map_in_k.png)

---

## 14. 去除探测器网格

MCP（六角）和网栅（方形）图案会印在每一幅图像上。**请最先做这一步**：只有仍在
探测器原始像素上的数据才能去网格，因此对已转换、校正或平滑过的数据，该按钮会
拒绝并说明原因。

在 map（或其 slit cut）或 cut 上 *Functions → Remove detector grid...*。
对 map 不需要参考测量——网格在每个切片中都位于同样的像素上，而光电子信号会移动。
点击 **Find and remove the grid**：

![去网格](images/degrid_dialog.png)

- 所选切片的 **Before / After**，使用同一色阶（勾选 **Show −∂²I/∂E²** 可让残余
  网格一目了然）。
- 找到的**网格图案**、其在 **k 空间的区域**（红色），以及每个切片拟合的**对比度
  与平移**。
- 网格区域中剩余的功率（1 表示网格已完全去除）。
- **Also list the grid pattern** 会保留一个 `[grid]` 数据集；之后用相同透镜模式和
  通能测得的 cut 会自动找到并使用它。

---

## 15. 曲线：EDC、MDC 与自旋分辨数据

任意图像下方的 **EDC → list** / **MDC → list**，或自旋分辨文件，都会得到曲线
数据集，在**曲线窗口**中打开：

![曲线窗口](images/curve_viewer.png)

可以隐藏通道、用 **Offset** 做成瀑布图、显示 **Error bars**（来自 σ 通道）；
**Range** 标记一个 x 范围供各操作使用。*Functions → Transform* 提供
**Crop**（裁剪）、**Bin**（合并）、**Normalise**（归一）、**Subtract a
background**（扣背景：常数、线性、Shirley）、**Shift the axis**（平移横轴）和
**Add counting errors (√N)**（添加计数误差）；每个操作都生成新数据集，并正确
传播 σ。

### 15.1 曲线拟合

*Functions → Curve fit...*：**Peaks on a background**（Lorentzian、Gaussian 或
Voigt 峰，点击放置或 **Find peaks** 自动寻峰；无 / 常数 / 线性 / Shirley 背景；
仪器分辨率；可选 **× Fermi edge** 费米截断），或者直接拟合 **Fermi edge**（费米边）。
权重可取通道自身的 σ、√N 或均匀权重。每个参数都给出误差，残差以 σ 为单位显示。

![曲线拟合](images/curve_fit_dialog.png)

**Fit to list** 保存数据、拟合、残差及各分量；**Copy results** 复制制表符分隔的
结果表。费米边拟合可用 **List a copy with E_F = 0** 生成平移后的副本。

### 15.2 自旋分析

自旋 EDC 打开时每个自旋通道一条曲线：

![自旋 EDC](images/spin_edc_viewer.png)

*Functions → Spin analysis...*：

![自旋分析](images/spin_analysis.png)

- 表格读取每个通道的标签（样品台设置、靶磁化方向、计为"上"的轴和符号）；
  **Axis** 与 **Counts as** 可手动修改。
- **Cross ratio**（交叉比，通道齐全时为默认）可精确消除透过率和靶反射率；也可用
  **单对通道**，但它带有仪器不对称性，程序会测量并报告该不对称性（ε）。
- **S_eff**（Sherman 函数，默认 0.2）、**Bin**、可选的 **Zero reference** 窗口。
- 图：原始通道、带误差棒的极化率 P，以及带误差的 I↑ / I↓。
- **P to list**、**I↑ / I↓ to list**、**As a figure**（打开作图器）。

---

## 16. 空间分辨扫描（nano-ARPES）

实空间扫描打开时左侧是**空间分布图**（每个位置的总强度），右侧是**光标处的
E-k 谱**。拖动红色靶标即可查看另一点的谱。在面板上右键可进行两种积分：

- 在空间图上，**Integrate selection → E vs k panel** 把样品某一区域内的谱求和；
- 在谱上，**Integrate selection → spatial panel** 画出所选 (k, E) 窗口内的强度
  在样品上的分布——用于寻找畴区或样品边缘。

![空间扫描](images/spatial_scan.png)

---

## 17. 作图与导出

### 17.1 作图器（Figure composer）

在列表中选中若干 cut（或曲线）→ 右键 → **Open in figure composer...**（或在冻结切片
窗口中点 **Figure composer...**）。作图器按论文格式
排版：同一列面板尺寸一致、只在外侧边缘标注坐标、自动添加 (a)/(b)/(c) 标号，
并可选择共用色阶。

![作图器](images/figure_composer.png)

**Figure tools...** 可设置所有细节（毫米为单位的尺寸、字体、刻度、色标、线条、箭头、
文字、插图、比例尺）——注意不要与查看窗口的 **Plot tools** 混淆，后者是在实时图上
画辅助线；**Add panel...** 添加数据集；**PDF**、**PNG**、**SVG**
导出，**Copy** 复制到剪贴板。

### 17.2 切片组图（Multi-slice figure）

在 map 上 *Functions → Multi-slice figure...*：把数据立方体沿某一轴的一系列切片——
等能面或 cut——排成一幅图，可等间距或指定位置，可积分，可共用色阶并自动加标题。

| 设置 | 结果 |
|---|---|
| ![切片组图设置](images/slice_figure_dialog.png) | ![切片组图](images/slice_figure.png) |

### 17.3 导出单个面板

在任意图像上右键 → **Export** → **Data (.nxs)...**、**Image...** 或 **Axes only
(.eps)...**。导出的总是该面板当前显示的内容（包括显示范围）：数据导出为本程序
可重新读取的数据集；图片应用当前色表（PNG、TIFF 或 JPEG，可设定尺寸和 dpi，
可带或不带坐标轴，如有 Plot tools 辅助线也会一并画出）；或者把边框、刻度和标签导出为矢量 PostScript，便于在
Illustrator 中放入图片。

![导出](images/export_dialog.png)

在列表中右键 → **Save** 把数据集写入你自己的 `.nxs` 文件；选中多行时会写入同一个文件。

---

## 18. 自动保存、内存与恢复

- **自动保存**：每个计算得到的数据集（k-map、处理后的 cut、拟合结果……）一经生成，
  就以**不压缩**的方式写入本次会话的文件夹 `~/.arpes_viewer/sessions/`。实空间
  扫描除外（它们就是测量数据本身，而且很大）。
- **正常关闭即删除自动保存文件**：关闭主窗口并确认后（所需数据已保存，或选择
  *Close anyway*），程序会关闭所有查看窗口，并立即删除本次会话的文件夹。某个数据集
  一旦被保存到你自己的文件，它的自动保存副本也会立即删除。
- **崩溃后的恢复**：只有程序没有正常关闭（崩溃、被强制结束）时，会话文件夹才会
  保留：下次启动时会提示恢复其中的数据集，三天后自动删除。
- **操作日志**：`~/.arpes_viewer/sessions/operations.log`，每次存储 / 保存 /
  丢弃记录一行——用来查找数据的去向。
- **Memory 区**：显示本程序与系统的内存占用、打开的窗口数与 HDF5 文件数，以及
  内存中的计算数据量。**Free memory** 丢弃计算数据的缓存副本（需要时从自动保存的
  文件重新读取）；**Close all viewers** 关闭所有窗口并释放内存。
- **Open auto-save folder**（Memory 区下方）在系统文件管理器中打开自动保存文件夹，
  可以手动删除旧的临时副本。当前会话的文件夹（状态栏中显示其名称）正在使用，请勿删除。
- **耗时操作**（读取大 map、转换、去网格、拼接文件夹）在后台线程运行，带进度条，
  且可以真正取消。
- **关闭程序**时如有未保存的计算数据，会先询问，并列出这些数据集的名字：
  **Save to a file...**（保存）、**Close anyway**（不保存直接关闭，自动保存副本
  随之删除）或 **Cancel**（取消关闭）。

## 19. 技巧与常见问题

- **Functions 菜单里显示 "Not available here"**——把鼠标悬停在该条目上，提示会说明
  原因（例如 MDC/EDC 拟合要求 cut 已经是 Å⁻¹ 单位）。
- **没有提供 k 转换**——该 map 的第一个轴不是角度（而是光子能量、温度等）。请在
  载入窗口检查 *First axis of a map is*。
- **对 E − E_F 的 cut 做 k 转换**——能量偏移会根据费米拟合自动填好；如果是用其他
  方式平移的能量轴，请在 *Energy offset* 中填入 E_F 对应的动能（光子能量 − 功函数）。
- **文件被误读**——在载入窗口中手动指定读取器，并参阅 README 中的 *Assumptions
  worth knowing about*（SOLEIL 文件的轴顺序）。
- **在大 map 上拖动光标卡顿**——从文件惰性读取的 map 每个窗口重新切片约需 0.1 s，
  因此其他窗口会等鼠标停下再更新；map 读入内存后即可实时跟随。
- **输入数字**——无论系统语言如何，所有数字框都使用小数点（`0.05`），不接受逗号。
- **系统开始使用交换内存**——留意 Memory 进度条（超过内存三分之一变为琥珀色，
  超过一半变为红色），使用 **Free memory** 或 **Close all viewers**。
