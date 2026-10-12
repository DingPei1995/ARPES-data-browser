# ARPES viewer 用户指南

*[English version / 英文版](USER_GUIDE_en.md) · [Version française / 法文版](USER_GUIDE_fr.md)*

**程序版本 2026-10-12** · 源代码：
<https://github.com/DingPei1995/ARPES-data-browser>

本指南分为四个部分：

- **第一部分　安装与设置**：版本号、源码地址、支持的数据类型，以及 Python
  环境的配置与程序启动。
- **第二部分　实验中的数据可视化**：主窗口、右键菜单，以及 SPEM、map、cut
  等数据的查看窗口。这些是在线站上快速判断"下一步测哪里"的工具。
- **第三部分　高级数据处理（用于论文）**：去网格、费米能级、k 转换、布里渊区、
  kz、曲率、二色性、峰拟合与作图，按论文数据通常的处理顺序排列。
- **第四部分　常见问题与排错**：自动保存与崩溃恢复、内存，以及常见问题的解答。

所有截图都来自实际运行的程序，所用数据是 `devtools/demo_data.py` 生成的模拟数据，
因此没有实验数据也可以跟着操作一遍（见[第 3 节](#3-用演示数据试用)）。每一步背后的
公式、测试和设计取舍见 [README](../README.md) 和 [CHANGELOG](../CHANGELOG.md)。

---

## 目录

**第一部分　安装与设置**

1. [关于本程序](#1-关于本程序)
2. [安装 Python 并启动程序](#2-安装-python-并启动程序)
3. [用演示数据试用](#3-用演示数据试用)

**第二部分　实验中的数据可视化**

4. [主窗口与列表右键菜单](#4-主窗口与列表右键菜单)
5. [载入数据](#5-载入数据)
6. [图像窗口：基本操作与右键菜单](#6-图像窗口基本操作与右键菜单)
7. [SPEM 实空间扫描：寻找测量位置](#7-spem-实空间扫描寻找测量位置)
8. [Map：等能面及其两个切面](#8-map等能面及其两个切面)
9. [Cut：快速查看](#9-cut快速查看)
10. [光子能量扫描、曲线与自旋 EDC：快速查看](#10-光子能量扫描曲线与自旋-edc快速查看)
11. [Map 的三维视图](#11-map-的三维视图)

**第三部分　高级数据处理（用于论文）**

12. [推荐的处理顺序](#12-推荐的处理顺序)
13. [去除探测器网格](#13-去除探测器网格)
14. [费米能级与费米面校正](#14-费米能级与费米面校正)
15. [从角度到动量](#15-从角度到动量)
16. [布里渊区与任意路径切面](#16-布里渊区与任意路径切面)
17. [光子能量扫描：kz 处理与 kz → 动量](#17-光子能量扫描kz-处理与-kz--动量)
18. [数据处理：平滑、求导、曲率](#18-数据处理平滑求导曲率)
19. [Cut 运算：二色性与参考谱归一](#19-cut-运算二色性与参考谱归一)
20. [MDC / EDC 拟合：色散、v_F、m\*、自能](#20-mdc--edc-拟合色散v_fm自能)
21. [曲线拟合与自旋分析](#21-曲线拟合与自旋分析)
22. [作图、导出与保存](#22-作图导出与保存)

**第四部分　常见问题与排错**

23. [自动保存、崩溃恢复与内存](#23-自动保存崩溃恢复与内存)
24. [常见问题](#24-常见问题)

---

# 第一部分　安装与设置

## 1. 关于本程序

| | |
|---|---|
| **版本** | 2026-10-12（最后修改于 2026 年 10 月 12 日）。版本号也写在 [README](../README.md) 顶部；[CHANGELOG.md](../CHANGELOG.md) 记录了每一处改动。 |
| **源代码** | <https://github.com/DingPei1995/ARPES-data-browser>（程序位于 `ARPES_viewer/` 目录） |
| **编写语言** | Python 3、PyQt5 与 pyqtgraph；数值计算与文件读写使用 numpy、scipy、h5py |
| **运行平台** | Windows、macOS、Linux——凡是能运行 PyQt5 的系统 |

### 支持读取的数据

| 线站 / 来源 | 文件 | 打开为 |
|---|---|---|
| SOLEIL **ANTARES** | NeXus/HDF5 `*.nxs` | cut、偏转器 map、实空间（SPEM）扫描 |
| SOLEIL **CASSIOPEE** | Scienta SES 文本文件：单张谱，以及按编号组成的文件夹 | cut；一个文件夹拼成一张 map 或一个光子能量（kz）扫描 |
| CASSIOPEE **自旋端站**（MBS A-1） | `*.krx` 与 `*.txt` | cut、偏转器 map、自旋分辨 EDC |
| 本程序 | 本程序自己保存的 `*.nxs` | 本程序保存的一切数据，并带有生成它的处理步骤 |

增加一个线站，只需在 `loader/` 中加一个读取器文件（见 README 中的 *Adding a beamline*）。

### 数据类型

程序中的每个数据集都属于下列类型之一。类型决定用哪个窗口打开它、可以对它使用
哪些功能。

| 类型 | 列表中显示为 | 轴 | 含义 |
|---|---|---|---|
| `cut` | `[Cut]` | 角度（或 k）× E | 一张谱 |
| `map` | `[Map]` | 角度 × 角度 × E | 偏转器或极角扫描 |
| `k_map` | `[k-map]` | kx × ky × E | 已转换到动量空间的 map |
| `kz_map` | `[kz map]` | hν × 角度 × E | 光子能量扫描 |
| `kz_map_k` | `[kz map (k)]` | k_z × k∥ × E | 转换到动量空间的光子能量扫描 |
| `spem_4d`、`spem_1d` | `[SPEM]` | x（× y）× 角度 × E | 实空间（nano-ARPES）扫描 |
| `edc`、`mdc`、`spin_edc` | `[EDC]`、`[MDC]`、`[Spin EDC]` | E（或 k）× 通道 | 曲线，可带 σ 误差通道 |

## 2. 安装 Python 并启动程序

**Python 版本**：需要 Python **3.9 或更高**（已测试到 3.13）。任何发行版都可以：
python.org 官方版、Anaconda/Miniconda，或 Linux 自带的 Python。

**依赖包**：共五个：

```bash
pip install PyQt5 pyqtgraph numpy h5py scipy
```

建议装在独立的环境中，这样程序既不依赖、也不会影响本机上其他 Python 工作：

```bash
# 使用 venv
python -m venv arpes-env
source arpes-env/bin/activate          # Windows：arpes-env\Scripts\activate
pip install PyQt5 pyqtgraph numpy h5py scipy

# 或使用 conda
conda create -n arpes python=3.12 pyqt pyqtgraph numpy h5py scipy -c conda-forge
conda activate arpes
```

最近一次测试所用版本：PyQt5 5.15、pyqtgraph 0.14（0.13 也可）、numpy 2.3、
h5py 3.16、scipy 1.18。

**获取代码**：在 GitHub 页面上 *Code → Download ZIP* 下载，或者：

```bash
git clone https://github.com/DingPei1995/ARPES-data-browser.git
cd ARPES-data-browser/ARPES_viewer
```

**启动**：在 `ARPES_viewer` 目录下运行：

```bash
python ARPES_viewer.py
```

**设置与工作文件**保存在用户主目录下的 `.arpes_viewer` 文件夹中：

| 文件 | 内容 |
|---|---|
| `~/.arpes_viewer/config.json` | 上次使用的数据文件夹，以及新窗口默认的色表 |
| `~/.arpes_viewer/sessions/` | 计算结果的自动保存副本（见[第 23 节](#23-自动保存崩溃恢复与内存)） |
| `~/.arpes_viewer/sessions/operations.log` | 每次自动保存、保存、删除各记一行 |

删除 `config.json` 即可恢复默认设置；程序不保存其他任何东西。

**检查安装是否正常**：测试套件不需要显示器即可运行：

```bash
python -m pytest            # 在 ARPES_viewer 目录下运行；需要 pip install pytest
```

## 3. 用演示数据试用

```bash
python devtools/demo_data.py            # 生成 ./demo_data/*.nxs
```

会生成六个模拟数据文件：一张费米面 map、一对 LH / LV 偏振下的 cut、一个光子能量
扫描、一条自旋分辨 EDC、一个含两个畴区的薄片的小实空间扫描，以及一张叠加了探测器
网格的 map。它们来自一个简单的能带模型（Γ 点的六角电子口袋与空穴带、K 点的空穴
口袋），包含费米截断、分析器分辨率与泊松噪声，打开方式与实测数据完全相同。把主
窗口的文件夹框指向 `demo_data/` 并加入它们即可。

修改窗口之后，如需重新生成本指南的全部截图：

```bash
python devtools/guide_screenshots.py    # 写入 docs/images/*.png
python devtools/build_guide_pdfs.py     # 重新生成 PDF（需要 pandoc 和 playwright）
```

---

# 第二部分　实验中的数据可视化

## 4. 主窗口与列表右键菜单

![主窗口](images/main_window_info.png)

主窗口是一个**数据浏览器**：每行对应一个*数据集*，而不是一个文件。一个 `.nxs`
文件常常包含多个测量，每个测量单独一行。名字后面标明数据类型：`[Cut]`、`[Map]`、
`[kz map]`、`[k-map]`、`[SPEM]`、`[EDC]`、`[Spin EDC]` 等。

| 要做的事 | 方法 |
|---|---|
| 在独立窗口中打开数据集 | **双击**该行，或点 **Open selected** |
| 对一个或多个数据集操作 | 选中后**右键**（菜单见下） |
| 查看各轴的首末值、点数、步长 | 右键 → **Show information**（填入 *Data Information* 表） |
| 同时对多个数据集截取、归一或合并像素 | **Data operations...**（[第 18 节](#18-数据处理平滑求导曲率)） |
| 平滑、求导、曲率、扣背景、对称化 | **Process...**（[第 18 节](#18-数据处理平滑求导曲率)） |
| 查看并释放内存 | 底部的 **Memory** 区（[第 23 节](#23-自动保存崩溃恢复与内存)） |

### 列表的右键菜单

选中一行或多行后右键。菜单作用于**整个选择**（按住 Shift 或 Ctrl 点击可多选）：

![列表右键菜单（选中两张 cut 时）](images/list_menu.png)

| 菜单项 | 适用于 | 作用 |
|---|---|---|
| **Open** | 一个数据集 | 在查看窗口中打开——与双击相同 |
| **Open slit cut** / **Open deflector cut** | 一个或多个 map | 只打开每张 map 的这一切面，不打开等能面 |
| **Show information** | 一个数据集 | 在信息窗口中显示类型、各轴与全部元数据 |
| **Rename...** | 一个数据集 | 修改它在列表中显示和保存时使用的名字（磁盘上的文件名不变） |
| **Plot as a figure...** | cut 或曲线 | 在作图器中每个数据集一个面板（[第 22 节](#22-作图导出与保存)） |
| **Stack plot (EDC / MDC)...** | 一张 cut | EDC 或 MDC 瀑布图（[第 9 节](#9-cut快速查看)） |
| **MDC / EDC fit...** | 一张已转换为 k 的 cut | 峰拟合面板（[第 20 节](#20-mdc--edc-拟合色散v_fm自能)） |
| **3D view...** | 一张 map | 正交切片、缺角立方体与投影（[第 11 节](#11-map-的三维视图)） |
| **Cut arithmetic on the two...** | 恰好两张 cut | 二色性、比值、差值（[第 19 节](#19-cut-运算二色性与参考谱归一)）；先选的为 A |
| **Save dataset(s)...** | 任意 | 把所选数据写入你自己的一个 `.nxs` 文件 |
| **Remove from list** | 任意 | 确认后从列表中移除；磁盘上的文件永远不会被删除 |
| **Session log...** | — | 本次会话自动保存、保存、删除了什么，以及文件在哪里（[第 23 节](#23-自动保存崩溃恢复与内存)） |

菜单只根据所选行的类型决定，不读取任何数据，所以无论文件多大都能立即弹出。不适用
于当前选择的条目显示为**灰色**；鼠标悬停时提示会说明原因（例如 *Cut arithmetic*
需要恰好两张 cut，*3D view* 需要 map）。

**Show information** 还会在独立窗口中显示该数据集的全部元数据，带过滤框，
并可 **Export metadata (.csv)**：

![数据集信息](images/information_window.png)

**Clear list** 会先确认；如果本次会话中有计算得到但尚未保存的数据，会提示保存。
磁盘上的文件永远不会被删除。

## 5. 载入数据

**Load data...** 打开载入窗口；**Quick add...** 跳过该窗口，直接以自动识别的设置
加入所选文件——实验期间最快的方式。

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
还是 kz map，由谱旁边的 `_i` 参数文件自动判断。

## 6. 图像窗口：基本操作与右键菜单

每个数据集在自己的窗口中打开，因此可以把 map、它的切面以及另一个文件并排摆放。

![带读数光标的 cut 窗口](images/cut_viewer.png)

- **顶部一行**：**Functions** 菜单（本窗口能做的所有事）、色表（每个窗口独立）
  及 **Flip**（反转），以及视图控制——坐标范围、**Auto**、**Inv**（反向）、
  **Reset**、**Grid**。
- **图像上方**：**Levels** 色阶条（拖动两个滑块，或输入 **Min%/Max%** 百分位后按
  **Clip**；在 **Clip** 按钮上*右键*可选预设：0.5–99.5 %、1–99 %、2–98 %、
  5–95 %、0–100 %）以及 **Gamma**。
- **图像下方**：**Smooth**（仅影响显示，不改变数据）和 **Lock x/y**（锁定纵横比）。
- **鼠标**：左键拖动平移，滚轮缩放，**右键拖动**分别缩放两个轴。

### 图像右键菜单

在任意图像上右键：

![图像右键菜单](images/image_menu.png)

| 菜单项 | 作用 |
|---|---|
| **View All** | 把视图恢复到显示整幅图像 |
| **Readout cursor**（读数光标） | 一个可拖动的十字光标。图像下方的 **Cursor** 行显示其位置和该处的数值；输入坐标并回车可把光标精确放到该位置。在 cut 上还会画出经过该点的 **EDC**（右侧）和 **MDC**（下方）；**EDC ±** / **MDC ±** 设置积分窗口，并在图上以半透明带显示。 |
| **Selection box**（选框） | 一个带两个角点手柄的矩形。角点坐标显示在图像下方的 **Selection** 行，可直接输入（再按 **Refresh**）。在 SPEM 窗口中选框就是积分区域（[第 7 节](#7-spem-实空间扫描寻找测量位置)）；在所有窗口中它还决定堆叠图的范围。 |
| **Line profile...**（线剖面） | 沿任意方向直线的强度，可设积分宽度（[见下](#线剖面)） |
| **Stack plot of the selection...** | 选框范围内的 EDC 或 MDC 瀑布图——没有选框时为整幅图。菜单项上会显示将要使用的范围。 |
| **Export → Data (.nxs) / Image / Axes only (.eps)** | 导出该面板当前显示的内容（[第 22 节](#22-作图导出与保存)） |

SPEM 窗口还有三个专用菜单项，见[第 7 节](#7-spem-实空间扫描寻找测量位置)。
pyqtgraph 自带的菜单项（鼠标模式、坐标轴子菜单、Plot Options、它自己的导出对话框）
已被有意去掉：其中有用的功能都放在了可见的视图控制栏上。

曲线旁的 **EDC → list** / **MDC → list** 把当前显示的曲线作为独立数据集加入主列表
（[第 10 节](#10-光子能量扫描曲线与自旋-edc快速查看)）。**Functions → Slice → Save
slice to the main list / Open slice in a new panel** 保存当前切片的副本。

### 线剖面

**右键 → Line profile...** 在图像上放一条橙色直线，并打开一个窗口显示沿线的强度。
拖动直线或它的两个端点即可移动；**Width** 设置直线两侧的积分宽度（以纵轴单位计），
**Samples** 设置采样点数，**To list** 把剖面作为数据集保存。关闭窗口时直线随之
移除。当能带斜着穿过图像、EDC 和 MDC 都无法沿着它取时，就用线剖面。

| cut 上的直线 | 沿线的强度 |
|---|---|
| ![cut 上的直线](images/line_profile_cut.png) | ![线剖面](images/line_profile.png) |

### Functions 菜单

| 窗口 | Analysis（分析） | Data operations（数据操作） | Visualization（可视化） |
|---|---|---|---|
| map / kz map / k-map | | Arbitrary cut、Map k conversion、kz map processing、kz → momentum、De-grid map | Brillouin zone（k-map）、Slice figure |
| cut | Fermi level、MDC / EDC fit | FS correction、Cut k conversion、Cut arithmetic、De-grid | |
| 曲线 | Curve fit、Spin analysis | Crop、Bin、Normalise、Background、Shift axis、Add √N errors | As a figure |

不适用于当前数据的条目不会显示；适用但暂时无法运行的条目显示为灰色，悬停提示会
说明原因。

## 7. SPEM 实空间扫描：寻找测量位置

实空间扫描（nano-ARPES，`[SPEM]`）打开时并排显示两个面板：左侧是**空间分布图**
（每个样品位置的总强度），右侧是 **E-k 谱**。空间图上的红色靶标选择显示哪个像素的
谱，右侧面板的标题会写明是哪个像素（如 *E vs k at x=20, y=24*）。

![SPEM 扫描：空间图与红色靶标处的谱](images/spatial_scan.png)

**拖动红色靶标**即可查看样品上任意位置的谱。这是判断一个薄片有没有能带、能带有多
锐利的最快方法。

### 在空间图上右键

![空间图的右键菜单](images/spem_spatial_menu.png)

| 菜单项 | 作用 |
|---|---|
| **Selection box** | 在空间图上画一个选框 |
| **Integrate selection -> E vs k panel** | 把**选框内所有像素的谱求和**：计数更多，得到样品某一区域的平均谱 |
| **Spatial map over the whole spectrum** | 在把谱上的选框积分到空间图之后出现（见下）：恢复为对所有 k、E 求和的空间图 |
| **Open this position in a new window** | 把当前显示的谱复制到独立窗口——可以同时打开多个，比较不同位置 |
| **Swap X and Y** | 转置显示空间图，适用于习惯另一个样品台方向为横轴的实验室。光标和选框仍停留在同样的样品像素上，只改变显示方式。 |
| Readout cursor、Line profile、Stack plot、Export | 与所有图像相同（[第 6 节](#6-图像窗口基本操作与右键菜单)） |

**对样品某一区域积分**：在空间图上右键 → **Selection box**，把选框拖到该区域
（或在 **Selection** 行输入角点坐标后按 **Refresh**），再右键 → **Integrate
selection -> E vs k panel**。此时右侧的谱就是选框内的总和，标题也会写明（如
*E vs k, summed over x 14 to 26, y 14 to 30 (63 pixels)*）：

![薄片上选框内求和的谱](images/spem_box_to_spectrum.png)

### 在谱上右键：某个特征来自样品的哪里？

![谱的右键菜单](images/spem_spectrum_menu.png)

这是反方向的操作。在谱上右键 → **Selection box**，在 (k, E) 中框住一个特征——
一条能带、一个肩峰、某一相的费米边——再右键 → **Integrate selection -> spatial
panel**。空间图会重新绘制为**仅该窗口内**的强度在每个位置的分布，标题写明所用的
窗口。下图中选框框住的是只有薄片的一个畴区才有的能带，空间图上恰好亮起了这个畴区：

![选框内能带的空间分布：一个畴区亮起](images/spem_box_to_map.png)

畴区、薄片边缘、台阶或第二相就是这样找到的。要恢复完整的空间图，在空间图上右键 →
**Spatial map over the whole spectrum**。

**线站上的实用流程**：

1. 打开总览扫描，拖动红色靶标扫过薄片，观察各处的谱。
2. 在谱上框住一个特征 → 积分到空间图，看这个特征在哪里最强。
3. 在空间图上框住最好的区域 → 积分到谱，以较好的统计确认能带。
4. 对候选位置使用 **Open this position in a new window**，并排比较，再从光标行读出
   坐标，用于下一次测量。

线扫描（`spem_1d`）也以同样方式打开，左侧为位置 × 角度；选框积分需要二维空间图，
因此只对二维扫描提供。

## 8. Map：等能面及其两个切面

map 打开时显示**等能面（constant-energy contour）**。顶部滑块选择能量（**Pos**
以 eV 为单位，**Ind** 为序号；**±** 设置积分范围，可选 **mean** 或 **sum**）。
**Deflector cut** 和 **Slit cut** 分别在新窗口中打开两个正交切面。

![map 的等能面](images/map_contour.png)

| Slit cut（狭缝角 vs 能量） | Deflector cut（偏转器角 vs 能量） |
|---|---|
| ![Slit cut](images/map_slit_cut.png) | ![Deflector cut](images/map_deflector_cut.png) |

**三个窗口共用一个光标。** 在任一窗口中打开读数光标（右键 → *Readout cursor*），
其余窗口同时打开；在任一窗口中移动光标，其余窗口会沿新位置重新切片。每个轴在所有
窗口中都用同一种颜色：**偏转器为红色、狭缝为绿色、能量为蓝色**。光标线的颜色是它
保持不变的那个轴的颜色，它所对应的 EDC 或 MDC 也用同一颜色绘制。积分宽度（± 框）
属于各自窗口，不会在窗口之间共享。

也可以直接从列表中打开 map 的切面（右键 → **Open slit cut** / **Open deflector
cut**），此时等能面窗口保持隐藏，需要时用 *Functions → Show the map* 调出。

在实验中，这一步回答的是：*样品是否对准、Γ 点在哪里、切面将沿哪个方向*——
在任何转换之前就能看出来。

## 9. Cut：快速查看

cut 在单个图像面板中打开（[第 6 节](#6-图像窗口基本操作与右键菜单)）。右键工具可以
快速给出答案：

- **读数光标**配合 **EDC ±** / **MDC ±**：任意一点处的曲线，可在窗口内积分。
- **Line profile...**：沿任意直线的强度。
- **堆叠图**：在 cut 上（或在列表中它的那一行）右键 → **Stack plot of the
  selection...**。若有选框，曲线就取在选框范围内：选框高大于宽时给出其角度范围内
  的 EDC，宽大于高时给出其能量范围内的 MDC。窗口中可设置曲线条数、合并、归一化、
  偏移、错切与配色；**Send to figure composer** 把它作为一个图面板送入作图器。

![堆叠图](images/stack_plot.png)

- **费米能级**：*Functions → Fermi level...* 几秒钟就能给出 E_F 和分辨率
  （[第 14 节](#14-费米能级与费米面校正)）——在改变光子能量之后检查能量标定很有用。

## 10. 光子能量扫描、曲线与自旋 EDC：快速查看

**光子能量扫描**在列表中显示为 `kz map`；其"等能面"是某一能量下光子能量对角度的
图，切面的打开方式与 map 相同。处理之前（第 17 节）各张谱的费米能级尚未对齐，
因此只作快速查看之用。

![kz map](images/kz_map_contour.png)

**曲线**：任意图像下方的 **EDC → list** / **MDC → list**，或自旋分辨文件，都会得到
曲线数据集，在**曲线窗口**中打开：

![曲线窗口](images/curve_viewer.png)

可以隐藏通道、用 **Offset** 做成瀑布图、显示 **Error bars**（来自 σ 通道）；
**Range** 标记一个 x 范围供各操作使用。*Functions → Data operations* 提供
**Crop**（裁剪）、**Bin**（合并）、**Normalise**（归一）、**Subtract a
background**（扣背景：常数、线性、Shirley）、**Shift the axis**（平移横轴）和
**Add counting errors (√N)**（添加计数误差）；每个操作都生成新数据集，并传播 σ。

**自旋 EDC**：打开时每个自旋通道一条曲线，可以在测量过程中检查各通道的计数率。
极化率本身在 *Functions → Spin analysis...* 中计算（[第 21 节](#21-曲线拟合与自旋分析)）。

![自旋 EDC](images/spin_edc_viewer.png)

## 11. Map 的三维视图

**3D view**（在列表中右键 map）首先询问区域与采样密度，因为每一帧都要对数据立方体
重新采样——*Fast*、*Balanced* 或 *Fine* 决定能否跟上鼠标操作。

![三维视图：区域与采样](images/volume_range_dialog.png)

之后拖动图像即可旋转。**Slices** 在指定位置显示三个正交平面，**Notched cube**
切掉立方体的一个角，**Projection** 沿某轴求和。右侧可设置色表、色阶、gamma、
方位角、仰角和深度比例；**Save the view...** 导出图像。

![三维视图](images/volume_view.png)

---

# 第三部分　高级数据处理（用于论文）

## 12. 推荐的处理顺序

每一步都会在列表中生成一个**新的数据集**，并记录生成它的步骤，因此原始数据永远不会
被修改，每个结果都可以追溯。顺序很重要，因为有些步骤只能用于探测器原始数据：

1. 若有网栅或 MCP 图案，先对原始 map 或 cut **去网格**（[第 13 节](#13-去除探测器网格)）。
2. **费米能级**；若费米边随分析器角度弯曲，再做 **FS 校正**（[第 14 节](#14-费米能级与费米面校正)）。
3. 对 cut 或 map 做 **k 转换**（[第 15 节](#15-从角度到动量)）；光子能量扫描则先做
   **kz 处理**，再做 **kz → 动量**（[第 17 节](#17-光子能量扫描kz-处理与-kz--动量)）。
4. 叠加**布里渊区**，并沿高对称路径取**任意切面**（[第 16 节](#16-布里渊区与任意路径切面)）。
5. 用**曲率 / 二阶导数**锐化能带；**二色性**（[第 18 节](#18-数据处理平滑求导曲率)、[第 19 节](#19-cut-运算二色性与参考谱归一)）。
6. **峰拟合**：色散、v_F、m\*、自能（[第 20 节](#20-mdc--edc-拟合色散v_fm自能)）；曲线拟合与自旋分析（[第 21 节](#21-曲线拟合与自旋分析)）。
7. **作图与导出**（[第 22 节](#22-作图导出与保存)）。

## 13. 去除探测器网格

MCP（六角）和网栅（方形）图案会印在每一幅图像上。**请最先做这一步**：只有仍在
探测器原始像素上的数据才能去网格，因此对已转换、校正或平滑过的数据，该功能会
拒绝并说明原因。

在 map（或其 slit cut）上 *Functions → De-grid map...*，在 cut 上 *De-grid...*。
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

## 14. 费米能级与费米面校正

### 14.1 费米能级

*Functions → Fermi level...* 对选框内求和得到的 EDC 拟合费米边（费米-狄拉克
分布 × 线性态密度，再与高斯分辨率卷积）。

![费米能级对话框](images/fermi_dialog.png)

1. 点击 **Show the selection box**，在 cut 上把选框拖到费米边处（或直接输入角度
   与能量范围）。角度范围决定哪些通道被求和，能量范围就是拟合窗口。

   ![cut 上的拟合选框](images/fermi_cut_box.png)

2. 先 **Estimate**，再 **Fit**。勾选 **Hold** 可固定某个参数（温度默认固定为
   文件记录值）。每个参数都给出不确定度。
3. **Offset energy axis** 生成一份能量轴为 **E − E_F** 的副本并加入列表。

*Advanced* 区域可以除去费米截断（用于能隙与参考样品），或逐通道拟合 E_F。

### 14.2 费米面（FS）校正

在 cut（或 map 的 slit cut）上 *Functions → FS correction...*，把本应平直的特征
拉直——通常是随分析器角度弯曲的费米边。沿该特征点击若干点，程序拟合多项式
（默认 2 阶），并把每一列在能量方向上平移。**From a reference (Au)...** 改用列表
中的金参考谱确定曲线。在 map 的 slit cut 上操作时，同样的平移会应用到整个数据
立方体。

![费米面校正](images/fs_correction_dialog.png)

## 15. 从角度到动量

### 15.1 Cut 的 k 转换

*Functions → Cut k conversion...* 把角度轴转换为 k∥（Å⁻¹）。一张 cut 是 (kx, ky)
平面上的一条线，通常不经过 Γ 点，所以对话框需要知道 Γ 的位置：手动输入角度，
或用 **From a converted map...** 从已转换的 map 继承。

![Cut k 转换](images/cut_k_dialog.png)

如果该 cut 已经平移到 E − E_F，对话框会自动把扣除的 E_F 填入 **Energy offset**
（下方的提示会说明）：转换时使用真实的动能，而结果的能量轴仍保持 E − E_F。

![转换到动量空间的 cut](images/cut_in_k.png)

### 15.2 Map 的 k 转换

在等能面窗口中选择 *Functions → Map k conversion...*。把读数光标拖到应成为
**k = (0, 0)** 的点（对话框会跟随光标），然后点 **Convert**，结果作为新的
`k-map` 行加入列表。

![Map k 转换](images/map_k_dialog.png)

- **Sample rotation** 绕表面法线旋转 map；**Set rotation from contour...** 通过在
  等能面上点击一个方向（例如沿镜面线的两点）读出旋转角。
- **Output grid**：按点数或按分辨率（Å⁻¹、meV）设置输出网格。
- **Energy offset** 在转换前加到能量轴上，用于能量轴不是真实动能的情况
  （对 E − E_F 的 map 会自动填好）。

![转换后的 k-map](images/k_map_contour.png)

## 16. 布里渊区与任意路径切面

### 16.1 布里渊区

在 k-map 上：*Functions → Visualization → Brillouin zone...*。

![布里渊区对话框](images/brillouin_dialog.png)

1. 先选**空间群**（1–230）；被空间群约束的晶格参数显示为灰色并自动跟随。也可
   切换到 **2D layer** 使用二维表面晶格。
2. **Cut plane**：表面法线 (hkl) 及沿法线的偏移——画出的是三维布里渊区被该平面
   截出的截面，适用于任意晶面，而不只是 (001)。**Show cut plane in 3D...** 可预览
   （拖动旋转，右键拖动缩放）。
3. 选择 **Conventional**（常规）或 **Irreducible**（不可约）区；设置 **Azimuth**
   （或 **Pick direction from contour...**）；**Repeat across the field of view**
   在视野内重复绘制。
4. **Moiré** 加入第二层（扭转或晶格失配）及其莫尔布里渊区。
5. **Plot**。

![叠加布里渊区的 k-map](images/k_map_with_bz.png)

### 16.2 任意路径切面

在任意等能面上 *Functions → Arbitrary cut...*：点击（或输入）最多六个点，切面
按顺序经过这些点——这里是 k-map 上的 Γ → M → K → Γ。勾选 **One dataset per
segment** 可让每一段分别成为一个数据集。

| 等能面上的路径 | 路径点 |
|---|---|
| ![任意切面路径](images/arbitrary_cut_path.png) | ![任意切面对话框](images/arbitrary_cut_dialog.png) |

结果是以沿路径距离为横轴的 cut；路径的拐点用绿色虚线标出，避免把拐点处的折线
误认为物理特征。

![沿 Γ–M–K–Γ 的切面](images/arbitrary_cut_result.png)

## 17. 光子能量扫描：kz 处理与 kz → 动量

### 17.1 kz map 处理

文件中并没有记录每张谱真正的费米能级位置。*Functions → kz map processing...*
直接从谱本身测量它。

1. slit cut 随之打开，上面有一个选框。把选框拖到费米边处——角度范围决定求和的
   通道，能量范围是拟合窗口。

   ![slit cut 上的选框](images/kz_processing_cut.png)

2. **Fit and assemble**：对每张谱使用同一个选框（按探测器像素序号），把每张谱
   平移使其自身的费米边位于零点，把数据立方体裁剪到所有谱共有的能量范围，并
   （默认）**按每张谱的总强度归一**，以消除光通量的变化。

![kz map 处理](images/kz_processing_dialog.png)

图中给出 E_F 随光子能量的变化。曲线平滑说明拟合可靠；拟合不可信的点会被标记，
并用相邻点插值补上。

### 17.2 kz → 动量

*Functions → kz -> momentum...* 把扫描转换为 (k_z, k∥)。

![kz → 动量](images/kz_conversion_dialog.png)

- **内势 V0** 是未知量：修改它时预览会在毫秒内重新转换一个能量切片。
- 给定**晶格**（空间群、a、b、c）、**表面法线**和面内方向后，预览上会画出布里渊
  区边界。
- **Scan V0** 在一系列 V0 下转换，报告 k_z 周期与晶格匹配时的 V0，以及 V0 能被
  确定到什么精度。
- **Measure a k_z period...**：点击相隔一个布里渊区的同一特征上的两点，程序列出
  与之匹配的倒格子晶面。
- 黄色横幅提示费米边仍随分析器角度弯曲——应先在 slit cut 上做 FS 校正。
- **Convert the whole cube** 生成 `kz map (k)`。

![转换后的 kz map](images/kz_map_in_k.png)

## 18. 数据处理：平滑、求导、曲率

在列表中选中一个或多个 cut（或 map），点击 **Process...**。每个结果都是新的
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

主窗口的 **Data operations...** 可同时对一个或多个数据集进行截取（truncate）、
自归一（self-normalise）或压缩（合并像素）——前提是它们类型相同，且各轴点数一致。

![数据操作](images/data_operations.png)

## 19. Cut 运算：二色性与参考谱归一

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

## 20. MDC / EDC 拟合：色散、v_F、m\*、自能

*Functions → MDC / EDC fit...*（或在列表中右键 → **MDC / EDC fit...**）逐行拟合峰，
再从拟合出的峰位读出能带。它只接受横轴已经是 Å⁻¹ 的 cut——以 eV·deg 为单位的
费米速度没有意义——所以请先转换 cut（[第 15 节](#15-从角度到动量)）。

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

## 21. 曲线拟合与自旋分析

### 21.1 曲线拟合

在曲线上 *Functions → Curve fit...*：**Peaks on a background**（Lorentzian、
Gaussian 或 Voigt 峰，点击放置或 **Find peaks** 自动寻峰；无 / 常数 / 线性 /
Shirley 背景；仪器分辨率；可选 **× Fermi edge** 费米截断），或者直接拟合
**Fermi edge**（费米边）。权重可取通道自身的 σ、√N 或均匀权重。每个参数都给出
误差，残差以 σ 为单位显示。

![曲线拟合](images/curve_fit_dialog.png)

**Fit to list** 保存数据、拟合、残差及各分量；**Copy results** 复制制表符分隔的
结果表。费米边拟合可用 **List a copy with E_F = 0** 生成平移后的副本。

### 21.2 自旋分析

在自旋 EDC 上 *Functions → Spin analysis...*：

![自旋分析](images/spin_analysis.png)

- 表格读取每个通道的标签（样品台设置、靶磁化方向、计为"上"的轴和符号）；
  **Axis** 与 **Counts as** 可手动修改。
- **Cross ratio**（交叉比，通道齐全时为默认）可精确消除透过率和靶反射率；也可用
  **单对通道**，但它带有仪器不对称性，程序会测量并报告该不对称性（ε）。
- **S_eff**（Sherman 函数，默认 0.2）、**Bin**、可选的 **Zero reference** 窗口。
- 图：原始通道、带误差棒的极化率 P，以及带误差的 I↑ / I↓。
- **P to list**、**I↑ / I↓ to list**、**As a figure**。

## 22. 作图、导出与保存

### 22.1 作图器（Figure composer）

在列表中选中若干 cut（或曲线）→ 右键 → **Plot as a figure**。作图器按论文格式
排版：同一列面板尺寸一致、只在外侧边缘标注坐标、自动添加 (a)/(b)/(c) 标号，
并可选择共用色阶。

![作图器](images/figure_composer.png)

**Plot tools** 可设置所有细节（毫米为单位的尺寸、字体、刻度、色标、线条、箭头、
文字、插图、比例尺）；**Add panel...** 添加数据集；**PDF**、**PNG**、**SVG**
导出，**Copy** 复制到剪贴板。堆叠图（**Send to figure composer**）和二色性结果
（**All three as a figure**）也可以送到这里。

### 22.2 切片组图（Slice figure）

在 map 上 *Functions → Slice figure...*：把数据立方体沿某一轴的一系列切片——
等能面或 cut——排成一幅图，可等间距或指定位置，可积分，可共用色阶并自动加标题。

| 设置 | 结果 |
|---|---|
| ![切片组图设置](images/slice_figure_dialog.png) | ![切片组图](images/slice_figure.png) |

### 22.3 导出单个面板

在任意图像上右键 → **Export** → **Data (.nxs)...**、**Image...** 或 **Axes only
(.eps)...**。导出的总是该面板当前显示的内容（包括显示范围）：数据导出为本程序
可重新读取的数据集；图片应用当前色表（PNG、TIFF 或 JPEG，可设定尺寸和 dpi，
可带或不带坐标轴）；或者把边框、刻度和标签导出为矢量 PostScript，便于在
Illustrator 中放入图片。在 SPEM 窗口中，在要导出的那个面板上右键——空间图和谱
分别导出。

![导出](images/export_dialog.png)

### 22.4 保存数据集

在列表中右键 → **Save dataset(s)...** 把数据集写入你自己的 `.nxs` 文件；选中多行时
写入同一个文件，每个数据集使用列表中的名字（可先用 **Rename...** 改名）。保存的
文件再次打开时，每个数据集都带有生成它的全部步骤；保存操作也会记入会话日志。

---

# 第四部分　常见问题与排错

## 23. 自动保存、崩溃恢复与内存

- **自动保存**：每个计算得到的数据集（k-map、处理后的 cut、拟合结果……）一经生成，
  就以**不压缩**的方式写入本次会话的文件夹 `~/.arpes_viewer/sessions/`。实空间
  扫描除外（它们就是测量数据本身，而且很大）。
- **正常关闭即删除自动保存文件**：关闭主窗口并确认后（所需数据已保存，或选择
  *Close anyway*），程序会关闭所有查看窗口，并立即删除本次会话的文件夹。某个数据集
  一旦被保存到你自己的文件，它的自动保存副本也会立即删除。
- **崩溃后的恢复**：只有程序没有正常关闭（崩溃、被强制结束、断电）时，会话文件夹
  才会保留：下次启动时会提示恢复其中的数据集，三天后自动删除。
- **操作日志**：`~/.arpes_viewer/sessions/operations.log`（或在列表中右键 →
  **Session log...**），每次自动保存、保存、删除各记一行并带路径——用来查找数据的
  去向。
- **关闭程序**时如有未保存的计算数据，会先询问，并列出这些数据集的名字：
  **Save to a file...**（保存）、**Close anyway**（不保存直接关闭，自动保存副本
  随之删除）或 **Cancel**（取消关闭）。
- **Memory 区**：显示本程序与系统的内存占用、打开的窗口数与 HDF5 文件数，以及
  内存中的计算数据量。**Free memory** 丢弃计算数据的缓存副本（需要时从自动保存的
  文件重新读取）；**Close all viewers** 关闭所有窗口并释放内存。
- **耗时操作**（读取大 map、转换、去网格、拼接文件夹）在后台线程运行，带进度条，
  且可以真正取消。

## 24. 常见问题

### 安装与启动

- **`ModuleNotFoundError: No module named 'PyQt5'`**（或 pyqtgraph、h5py 等）——
  依赖包装在了另一个 Python 中，而不是运行程序的这个。激活对应的环境
  （[第 2 节](#2-安装-python-并启动程序)），并用
  `python -c "import PyQt5, pyqtgraph, h5py, scipy"` 检查。
- **`qt.qpa.plugin: Could not load the Qt platform plugin "xcb"`**（Linux）——缺少
  系统的 Qt 库；用包管理器安装（Ubuntu 上为
  `sudo apt install libxcb-xinerama0 libxkbcommon-x11-0`），或改用 conda 安装。
- **字体或窗口过大 / 过小**（高分辨率屏幕）——程序跟随系统缩放；启动前在环境变量
  中设置 `QT_SCALE_FACTOR=0.9`（或其他数值）。
- **设置看起来不对**——删除 `~/.arpes_viewer/config.json`。

### 读取数据

- **文件被误读**——在载入窗口中手动指定读取器，并参阅 README 中的 *Assumptions
  worth knowing about*（SOLEIL 文件的轴顺序）。
- **没有提供 k 转换**——该 map 的第一个轴不是角度（而是光子能量、温度等）。请在
  载入窗口检查 *First axis of a map is*。
- **CASSIOPEE 文件夹被列成一张张单独的谱**——谱旁边的 `_i` 参数文件缺失或被改名；
  正是这些文件说明了扫描的是什么。

### 右键菜单

- **列表菜单中的条目是灰色的**——鼠标悬停在上面：提示会说明当前选择缺了什么
  （一行、两张 cut、一张 map、一张 k 空间的 cut……）。
- **"Integrate selection" 是灰色的**——先画选框（右键 → **Selection box**）。
- **没有 "Integrate selection" 这一项**——它只存在于带二维空间图的 SPEM 窗口；在线
  扫描、cut 或 map 上，选框用于堆叠图和坐标读取。
- **SPEM 空间图变暗或只显示样品的一部分**——它正在显示谱上某个选框的强度（标题会
  写明）。在空间图上右键 → **Spatial map over the whole spectrum**。
- **SPEM 空间图看起来被转置或镜像了**——**Swap X and Y**（在空间图上右键）和视图栏
  中的 **Inv** 只改变显示方式；ANTARES 扫描默认按线站的习惯方向打开。
- **线剖面的直线不见了**——关闭线剖面窗口时直线会随之移除；从右键菜单重新打开即可。
- **MDC / EDC 拟合拒绝某张 cut**——它的横轴仍是角度；请先转换（[第 15 节](#15-从角度到动量)）。

### 数据处理

- **Functions 菜单里显示 "Not available here"**——把鼠标悬停在该条目上，提示会说明
  原因。
- **去网格拒绝某个数据集**——它已被转换、校正或平滑过；请先对原始数据去网格
  （[第 13 节](#13-去除探测器网格)）。
- **对 E − E_F 的 cut 做 k 转换**——能量偏移会根据费米拟合自动填好；如果是用其他
  方式平移的能量轴，请在 *Energy offset* 中填入 E_F 对应的动能（光子能量 − 功函数）。

### 速度与内存

- **在大 map 上拖动光标卡顿**——从文件惰性读取的 map 每个窗口重新切片约需 0.1 s，
  因此其他窗口会等鼠标停下再更新；map 读入内存后即可实时跟随。
- **系统开始使用交换内存**——留意 Memory 进度条（超过内存三分之一变为琥珀色，
  超过一半变为红色），使用 **Free memory** 或 **Close all viewers**。
- **程序崩溃之后**——重新启动程序：它会提示恢复自动保存的数据集。日志
  （**Session log...**）写明了每个数据集在哪里。
