# 量子散射计算程序集 | Quantum Scattering Programs Collection

这个仓库包含多个用于量子散射计算的独立程序，每个目录对应一个完整的项目。以下是对各个项目的简要介绍。

This repository contains several independent programs for quantum scattering calculations. Each directory corresponds to a complete project. Below is a brief introduction to each project.

## 项目列表 | Project List

### 1. AAXV - 原子-双原子系统振动非弹性碰撞截面的评估程序
### 1. AAXV - A Program for Evaluating Vibrationally Inelastic Collisional Cross Sections of Atom-Diatom Systems

位置: `./A program to evaluate vibrationally inelastic collisional cross sections of atom-diatom systems/`

Location: `./A program to evaluate vibrationally inelastic collisional cross sections of atom-diatom systems/`

本程序用于计算原子与双原子分子之间的振动非弹性碰撞截面。程序基于扭曲波（Distorted Wave）或指数型扭曲波（Exponential Distorted Wave）近似方法，用于量子力学中的散射计算。

This program calculates vibrationally inelastic collisional cross sections between atoms and diatomic molecules. It is based on the Distorted Wave (DW) or Exponential Distorted Wave (EDW) approximation methods for scattering calculations in quantum mechanics.

原始参考文献：M.M. Novak, Computer Physics Communications 46 (1987) 417

Original reference: M.M. Novak, Computer Physics Communications 46 (1987) 417

### 2. FERM3D - 有限元R矩阵电子-分子散射代码
### 2. FERM3D - A Finite Element R-matrix Electron-Molecule Scattering Code

位置: `./FERM3D A finite element R-matrix electron molecule scattering code/`

Location: `./FERM3D A finite element R-matrix electron molecule scattering code/`

FERM3D是一个用于计算电子与分子之间散射过程的程序。它使用有限元方法实现R矩阵理论来求解散射问题。该程序可以处理弹性散射和光致电离过程。

FERM3D is a program for calculating scattering processes between electrons and molecules. It implements the R-matrix theory using finite element methods to solve scattering problems. The program can handle elastic scattering and photoionization processes.

---

## AAXV 详细说明 | AAXV Detailed Description

### 功能特点 | Features

- 计算原子-双原子分子系统的振动非弹性碰撞截面
- 支持扭曲波和指数型扭曲波近似方法
- 提供WKBJ近似方法选项
- 可调整总角动量和能量范围
- 计算S矩阵和T矩阵
- 灵活的输出控制选项

---

- Calculates vibrationally inelastic collisional cross sections for atom-diatom systems
- Supports both Distorted Wave and Exponential Distorted Wave approximation methods
- Offers WKBJ approximation method option
- Adjustable total angular momentum and energy ranges
- Calculates S-matrices and T-matrices
- Flexible output control options

### 安装与使用 | Installation and Usage

#### 系统要求 | System Requirements

- Fortran编译器（如gfortran、ifort等）
- BLAS和LAPACK库（数值计算）

---

- Fortran compiler (e.g., gfortran, ifort)
- BLAS and LAPACK libraries (for numerical calculations)

#### 编译 | Compilation

使用Fortran编译器编译主程序文件`aaxv_v1_0.f`：

Compile the main program file `aaxv_v1_0.f` using a Fortran compiler:

```bash
gfortran -o aaxv aaxv_v1_0.f -lblas -llapack
```

#### 输入文件 | Input Files

程序需要一个名为`DIA.DAT`的输入文件，包含所有计算参数。主要参数包括：

The program requires an input file named `DIA.DAT` containing all calculation parameters. Main parameters include:

- `IPRNT`：输出详细程度（1-4）
- `JTMIN`/`JTMAX`：总角动量的最小/最大值
- `EMIN`/`ETOTAL`：能量范围
- `HBAR`/`RSCALE`/`VSCALE`：能量和距离的缩放因子
- `RMASS`：约化质量
- `RBEG`/`REND`：积分范围
- `IEDW`：选择扭曲波(-2)或指数型扭曲波(-1)方法
- `LWKBJ`：是否使用WKBJ近似（0或1）
- `NVMIN`/`NVMAX`：振动量子数范围

---

- `IPRNT`: Output detail level (1-4)
- `JTMIN`/`JTMAX`: Minimum/maximum values of total angular momentum
- `EMIN`/`ETOTAL`: Energy range
- `HBAR`/`RSCALE`/`VSCALE`: Energy and distance scaling factors
- `RMASS`: Reduced mass
- `RBEG`/`REND`: Integration range
- `IEDW`: Select distorted (-2) or exponential distorted (-1) wave method
- `LWKBJ`: Whether to use WKBJ approximation (0 or 1)
- `NVMIN`/`NVMAX`: Vibrational quantum number range

#### 运行程序 | Running the Program

编译后，执行程序：

After compilation, execute the program:

```bash
./aaxv
```

确保`DIA.DAT`文件位于当前目录。

Make sure the `DIA.DAT` file is in the current directory.

### 输入文件示例 | Sample Input File

```
示例计算：He-H2系统
&EDW
  IPRNT=2, IDWINTG=1, IDISC=0,
  JTMIN=0, JTMAX=20, JTSTEP=2,
  ETOTAL=0.1, NESTEP=0,
  HBAR=1.0, RSCALE=1.0, VSCALE=1.0,
  RMASS=1.819, DFACTOR=30.0,
  RBEG=0.5, REND=12.0, EPSIL=0.002,
  IEDW=-1, LWKBJ=1, STWKBJ=0.01,
  VIBQUAN=0.5159, NVMIN=0, NVMAX=1,
  NGAM=8, NROINT=148,
  RVDW=3.0, DVDW=0.003, ALVDW=1.5,
  RDIA=1.4, DDIA=4.7, ALDIA=1.9, DIMASS=0.5
&END
```

## FERM3D 详细说明 | FERM3D Detailed Description

### 安装说明 | Installation Instructions

FERM3D提供了三种不同的安装脚本：
1. 适用于Intel编译器（文件"install.adv.pl.Intel"），适用于Intel或AMD架构
2. 适用于Tru64（文件"install.adv.pl.Tru64"），在HP Alpha机器上测试，使用HP Fortran编译器
3. 适用于Tru64（文件"install.adv.pl.Tru64.SuperLU"），使用SuperLU代替Pardiso

FERM3D provides three different installation scripts:
1. One for Intel compilers (file "install.adv.pl.Intel"), suitable for Intel or AMD architectures
2. One for Tru64 (file "install.adv.pl.Tru64"), tested on HP Alpha machines with HP Fortran compiler
3. Another for Tru64 (file "install.adv.pl.Tru64.SuperLU"), using SuperLU instead of Pardiso

### 依赖库 | Required Libraries

- LAPACK/BLAS（必需）
- PARDISO库（或可选的SuperLU库）
- CERNLIB和ARPACK（可选）

- LAPACK/BLAS (required)
- PARDISO library (or optional SuperLU library)
- CERNLIB and ARPACK (optional)

### 安装步骤 | Installation Steps

1. 安装所有必需的库
2. 切换到`<FERM3D-dir>/make`目录
3. 在`install.adv.pl<system>`中修改编译器和库的位置（在脚本开头定义）
4. 运行Perl安装脚本：`perl install.adv.pl<system>`
5. 所有二进制文件将安装在`<FERM3D-dir>/bin/`目录中

1. Install all required libraries
2. Change to the `<FERM3D-dir>/make` directory
3. Modify compiler and library locations in "install.adv.pl<system>" (defined at beginning of script)
4. Run the Perl installation script: `perl install.adv.pl<system>`
5. All binaries will be installed in the `<FERM3D-dir>/bin/` directory

### 运行示例 | Running Examples

示例位于"examples"目录中。提供了三个示例：CO2中性（弹性散射）、CO2+（光致电离）和苯（弹性散射）。

Examples are in the "examples" directory. Three examples are provided: neutral CO2 (elastic scattering), CO2+ (photoionization), and benzene (elastic scattering).

## 通用参考公式 | Common Reference Formulas

对于扭曲波近似，程序求解以下方程：

For the Distorted Wave approximation, the program solves the following equation:

$$\left(\frac{d^2}{dr^2} + k_j^2 - \frac{l_j(l_j+1)}{r^2} - V_{jj}(r)\right)\psi_j(r) = 0$$

其中$k_j$是通道波数，$l_j$是角动量，$V_{jj}$是势能。

where $k_j$ is the channel wave number, $l_j$ is the angular momentum, and $V_{jj}$ is the potential energy.

散射截面由T矩阵计算：

The scattering cross sections are calculated from the T-matrix:

$$\sigma_{j \rightarrow j'} = \frac{\pi}{k_j^2} \sum_l (2l+1) |T_{j'j}^l|^2$$

## 许可证 | License

这个仓库中的代码基于原始开发者的工作，请参考相应的原始文献了解使用条款。

The code in this repository is based on the work of the original developers. Please refer to the corresponding original literature for terms of use.#   s c a t t e r i n g 1 2  
 