# BluePebble 最终分析流程：完整顺序、目的与运行方式

本文档描述如何按最终项目逻辑重建 Step01–Step17。`scripts/` 保留 Step17 最终脚本的分析逻辑、参数和依赖关系，但发布前已把个人项目路径、home 路径、集群用户名和 Slurm 账户替换为通用占位符。首次执行前必须先阅读 `PORTABILITY_AND_SECURITY.md`，并替换所有配置占位符。

## 一、运行前准备

### 1. 设置项目结构

GitHub 脱敏版本统一使用以下项目根目录占位符：

```text
/path/to/mollusc_project
```

运行前须把脚本和文档中的 `/path/to/mollusc_project`、`/path/to/miniconda3`、`/path/to/hpc_home`、`/path/to/eggnog.taxa.db`、`HPC_USER` 和 `YOUR_SLURM_ACCOUNT` 替换为目标集群上的实际授权值；也可以通过环境变量 `EGGNOG_TAXDB` 指定 eggNOG taxonomy 数据库。不要把含个人路径或凭据的定制版本再次提交到公共仓库。

至少应准备以下目录：

```text
mollusc_project/
├── scripts/
├── logs/
├── orthofinder_proteomes/
├── genome_fastas/
├── busco_downloads/lineages/metazoa_odb10/
├── conda_envs/
├── databases/gene_ontology/
├── results_or_reports/
└── result/
```

将仓库脚本复制到项目脚本目录后赋予执行权限：

```bash
cd /path/to/mollusc_project
cp /path/to/repository/scripts/* scripts/
chmod u+x scripts/*
mkdir -p logs results_or_reports result
```

### 2. 输入数据

活跃蛋白组目录 `orthofinder_proteomes/` 应只包含 14 个非空 `.faa` 文件，文件名必须与脚本中的物种名一致。最终类群设计为：

- Bivalvia（6）：`Acanthocardia_echinata`、`Americardia_media`、`Cerastoderma_edule`、`Fragum_sueziense`、`Tridacna_gigas`、`Tridacna_maxima`。
- Gastropoda（8）：`Aplysia_californica`、`Berghia_stephanieae`、`Bullacta_exarata`、`Elysia_chlorotica`、`Elysia_crispata`、`Elysia_marginata`、`Onchidella_celtica`、`Plakobranchus_ocellatus`。

不要同时把 `Fragum_fragum.faa` 和 `Fragum_sueziense.faa` 留在活跃输入目录。前者可移到单独的历史输入目录，但不应进入 14 物种最终分析。

### 3. 软件和数据库

最终项目涉及：BUSCO、seqkit、OrthoFinder、DIAMOND、FAMSA、FastTree、IQ-TREE、CAFE5、eggNOG-mapper、HMMER、GO OBO、MAFFT、trimAl、HyPhy、Salmon、Python、R 和 ggplot2。已明确记录的版本包括 IQ-TREE 3.1.3、HyPhy 2.5.101、eggNOG-mapper 2.1.12、Pfam 38.2、Salmon 1.10.2、Python 3.13.11、R 4.5.2 和 ggplot2 4.0.3；详见 `metadata/SOFTWARE_AND_DATABASE_VERSIONS.tsv`。

### 4. 资源原则

类群分析及 Step04 之后的作业均不超过 12 CPUs/16 GB。最早的原始蛋白组 BUSCO 和 14 物种 combined OrthoFinder 脚本形成于该限制之前，归档中仍保留其实际请求（BUSCO 8 CPUs/32 GB；combined OrthoFinder 20 CPUs/80 GB）。若必须把所有重跑作业都压到 12 CPUs/16 GB，需要先在目标集群小规模验证内存，而不是直接假设 OrthoFinder 可以无损降配。

## 二、Step01–Step17 的执行顺序

### Step 01 — Proteome Quality Assessment

**目的：**检查每个蛋白组的 FASTA 合法性、蛋白数量、重复 ID、异常字符、BUSCO 完整度、单拷贝/重复/缺失比例；选择代表转录本，并以 BUSCO complete ≥80% 定义核心物种集。

**主要软件：**seqkit、BUSCO（`metazoa_odb10`）、Python。

**1A. 重建最终 *Fragum* 替代数据（仅在尚未准备 *F. sueziense* 时运行）：**

```bash
cd /path/to/mollusc_project
bash scripts/12_submit_fragum_replacement.sh
```

该链依次运行：

1. `10_download_fragum_sueziense.sh`：从 Ensembl Rapid Release 下载并核验 `GCA_963680895.1` 的 genome、peptide 和 GFF3；解压到对应目录；筛选代表蛋白。
2. `11_busco_fragum_sueziense.sh`：对代表蛋白运行 BUSCO 并刷新汇总。

**1B. 原始蛋白组 QC：**

```bash
bash scripts/03_submit_busco_qc.sh
```

该 wrapper 自动提交 BUSCO 数组 `01_proteome_busco_qc.sh`，成功后运行 `02_summarise_busco.sh`。

**1C. 代表蛋白与代表蛋白 BUSCO：**

Step01 与 Step02 的官方 wrapper 是连续链：

```bash
bash scripts/09_submit_primary_to_orthofinder.sh
```

其前 3 个依赖任务依次为 `04_prepare_primary_proteomes.sh` → `05_primary_busco_qc.sh` 数组 → `06_summarise_primary_busco.sh`。最后一个任务进入 Step02 的 combined OrthoFinder。

**主要输出：**

- `proteome_busco_summary.tsv`
- `primary_proteome_selection.tsv`
- `primary_proteome_busco_summary.tsv`
- `orthofinder_species_decisions.tsv`
- `orthofinder_primary_proteomes/*.faa`
- `orthofinder_core_complete80/*.faa`

---

### Step 02 — Orthogroup Inference

**目的：**在 14 物种 combined 基线和两个纲内数据集中推断 orthogroups、同源/旁系同源关系、基因拷贝数矩阵、单拷贝直系同源基因和 OrthoFinder 基础树。

**主要软件：**OrthoFinder、DIAMOND、FAMSA、FastTree。

**2A. Combined 14 物种基线：**

由上一步的 `09_submit_primary_to_orthofinder.sh` 自动以 `08_run_orthofinder.sh` 提交。若 Step01 已单独完成，也可直接提交：

```bash
sbatch scripts/08_run_orthofinder.sh
```

**2B. 分纲运行 Bivalvia 与 Gastropoda：**

```bash
bash scripts/15_submit_class_analyses.sh
```

`14_run_class_orthofinder_iqtree.sh` 会为每个纲建立软链接输入目录，运行 OrthoFinder，并继续生成初始全 ModelFinder IQ-TREE。两个纲可并行运行。

**主要输出：**

- `Orthogroups/Orthogroups.tsv`
- `Orthogroups/Orthogroups.GeneCount.tsv`
- `Orthogroups/Orthogroups_SingleCopyOrthologues.txt`
- `Species_Tree/SpeciesTree_rooted.txt`
- OrthoFinder 的 gene trees、duplication tables 和 hierarchical orthogroups。

---

### Step 03 — Species Tree Reconstruction

**目的：**使用单拷贝 orthologue 拼接蛋白比对构建高支持度物种树，为 CAFE5 和 gene-tree/species-tree reconciliation 提供系统发育框架。

**主要软件：**IQ-TREE 3.1.3；1,000 次 UFBoot 和 1,000 次 SH-aLRT。

最终选择如下：

- Bivalvia：保留 `14_run_class_orthofinder_iqtree.sh` 完成的纲内 MFP 树。
- Gastropoda：使用 restricted-MFP 树。
- Combined 14 物种：使用 restricted-MFP 基线树。
- `13_run_iqtree_species_tree.sh` 的 combined 全 MFP 入口已被替代，不属于最终运行顺序。

提交最终 restricted-MFP 两项：

```bash
combined_iq=$(sbatch --parsable \
  --job-name=iqtree_restricted_combined \
  --export=ALL,ANALYSIS=combined \
  scripts/21_run_restricted_mfp_iqtree.slurm)

gastropoda_iq=$(sbatch --parsable \
  --job-name=iqtree_restricted_gastropoda \
  --export=ALL,ANALYSIS=gastropoda \
  scripts/21_run_restricted_mfp_iqtree.slurm)

printf 'combined=%s\ngastropoda=%s\n' "$combined_iq" "$gastropoda_iq"
```

restricted ModelFinder 的候选集合为 `Q.INSECT,Q.YEAST,JTT,LG,WAG`，频率模式 `FU,F`，速率模式 `E,I,G,I+G,R`，最大 6 个速率类别，随机种子 12345。

**主要输出：**`SpeciesTree_IQTREE.treefile`、`.contree`、`.iqtree`、`.log` 和重命名后的拼接比对。

---

### Step 04 — Orthogroup Distribution Analysis

**目的：**量化 core、soft-core、shell、lineage-specific、species-specific orthogroups，统计蛋白分配率和未分配蛋白，描述两个纲内的 orthogroup 组成。

**主要工具：**Python，自定义统计脚本。

```bash
step04=$(sbatch --parsable scripts/16_step04_orthogroup_distribution.slurm)
echo "$step04"
```

此步骤只依赖 Step02，可与 Step03 后半段、Step05 的准备工作并行，但最终论文图表要使用权威运行 `run_18269662`。

**主要输出目录：**`results_or_reports/orthogroup_distribution/run_<JOB_ID>`。

---

### Step 05 — Gene-Family Expansion and Contraction Analysis

**目的：**把纲内 gene-count matrix 与物种树转换为 CAFE5 输入，估计家族演化率 λ，并识别分支上的显著扩张和收缩。

**主要软件：**CAFE5；相对根年龄固定为 100，仅用于探索性相对时间筛查。

首次使用时安装环境：

```bash
cafe_setup=$(sbatch --parsable scripts/18_step05_install_cafe5.slurm)
```

可选地验证 Bivalvia 输入准备：

```bash
cafe_check=$(sbatch --parsable --dependency=afterok:"$cafe_setup" \
  scripts/20_step05_validate_prepare.slurm)
```

两个纲可在各自物种树完成后并行提交：

```bash
biv_cafe=$(sbatch --parsable \
  --dependency=afterok:"$cafe_setup" \
  --job-name=step05_cafe_bivalvia \
  --export=ALL,GROUP=bivalvia \
  scripts/19_step05_run_cafe_group.slurm)

gas_cafe=$(sbatch --parsable \
  --dependency=afterok:"$cafe_setup":"$gastropoda_iq" \
  --job-name=step05_cafe_gastropoda \
  --export=ALL,GROUP=gastropoda \
  scripts/19_step05_run_cafe_group.slurm)
```

最终解释以 Base/error-model 输出为主。Gastropoda Gamma 模型虽然运行结束，但诊断显示全部输入家族的失败率超过 20%，因此不能把其显著性结果作为主要证据。

**主要输出：**`cafe_input_primary.tsv`、超度量相对时间树、`cafe_base_error/`、`cafe_gamma_k2/`、`cafe_input_summary.tsv` 和 manifest。

---

### Step 06 — Candidate Gene-Family Identification

**目的：**从 CAFE5 Base 模型提取显著扩张、收缩及终端分支变化家族，结合变化幅度和家族层级显著性建立候选集合与蛋白 FASTA。

**主要工具：**Python。

```bash
step06=$(sbatch --parsable \
  --dependency=afterok:"$biv_cafe":"$gas_cafe" \
  scripts/22_step06_candidate_families.slurm)
```

关键阈值：family FDR 0.05、high-confidence branch p 0.01、candidate branch p 0.05、high-confidence minimum absolute copy-number change 2。

**主要输出目录：**`results_or_reports/candidate_gene_families/base_model_v1/run_<JOB_ID>`。

---

### Step 07 — Functional Annotation

**目的：**为候选蛋白分配 eggNOG seed orthologues、功能描述、preferred name、GO、KEGG、EC、COG、CAZy 和 Pfam 转移注释，并生成 orthogroup 层级汇总。

**主要软件：**eggNOG-mapper 2.1.12 / DIAMOND。

两个纲并行，完成后合并：

```bash
biv_ann=$(sbatch --parsable --dependency=afterok:"$step06" \
  --job-name=step07_eggnog_bivalvia \
  scripts/23_step07_functional_annotation.slurm bivalvia)

gas_ann=$(sbatch --parsable --dependency=afterok:"$step06" \
  --job-name=step07_eggnog_gastropoda \
  scripts/23_step07_functional_annotation.slurm gastropoda)

step07_merge=$(sbatch --parsable \
  --dependency=afterok:"$biv_ann":"$gas_ann" \
  scripts/23_step07_combine.slurm)
```

**主要输出：**`protein_functional_annotations.tsv`、`candidate_family_functional_annotations.tsv`、GO/KEGG/Pfam 长表、coverage summary 和 combined summary。

---

### Step 08 — Functional Enrichment Analysis

**目的：**在候选家族内部背景下，对不同候选集合进行 Fisher exact test，并用 Benjamini–Hochberg 控制 FDR；分析 GO BP/MF/CC、KEGG、Pfam、CAZy 和 EC。

**主要工具：**Python，Fisher exact test，BH FDR。

```bash
step08=$(sbatch --parsable \
  --dependency=afterok:"$step07_merge" \
  scripts/24_step08_functional_enrichment.slurm)
```

注意：背景是“具有相应注释的 CAFE 显著候选 orthogroups”，不是全蛋白组，因此结果只能解释为候选集合内部的相对富集。

**主要输出目录：**`results_or_reports/functional_enrichment/orthogroup_level_v1/run_<JOB_ID>`。

---

### Step 09 — Candidate Sequence Validation

**目的：**用 Pfam 38.2 的直接 HMM 命中、蛋白长度、序列完整性、精确重复、分类学来源和高风险功能标记验证候选序列；区分 `pass`、`review`、`review_high` 和 `exclude_recommended`。

**主要软件：**HMMER 3.4、Pfam 38.2、Python。

先下载、校验并 `hmmpress` Pfam：

```bash
pfam_job=$(sbatch --parsable scripts/25_step09_prepare_pfam.slurm)
```

两个纲并行验证，随后合并：

```bash
biv_val=$(sbatch --parsable --dependency=afterok:"$pfam_job":"$biv_ann" \
  --job-name=step09_validate_bivalvia \
  scripts/25_step09_validate_group.slurm bivalvia)

gas_val=$(sbatch --parsable --dependency=afterok:"$pfam_job":"$gas_ann" \
  --job-name=step09_validate_gastropoda \
  scripts/25_step09_validate_group.slurm gastropoda)

step09_merge=$(sbatch --parsable \
  --dependency=afterok:"$biv_val":"$gas_val" \
  scripts/25_step09_combine.slurm)
```

**主要输出：**`sequence_validation.tsv`、`family_validation_summary.tsv`、`direct_pfam_domains.tsv`、`exact_duplicate_groups.tsv`、`high_risk_sequence_flags.tsv`、`validated_pass.faa`、`manual_review.faa` 和 `exclude_recommended.faa`。

---

### Step 10 — Candidate Gene-Family Phylogenetics

**目的：**从验证后的候选中为每个纲选择最多 30 个家族，进行 MAFFT 比对、trimAl 修剪和 IQ-TREE 基因树重建，从而检查候选家族内部复制和谱系特异分化。

**主要软件：**MAFFT、trimAl、IQ-TREE。

```bash
step10_prep=$(sbatch --parsable \
  --dependency=afterok:"$step09_merge" \
  scripts/26_step10_prepare.slurm)

biv_trees=$(sbatch --parsable \
  --dependency=afterok:"$step10_prep" \
  --array=0-29%3 \
  --job-name=step10_tree_bivalvia \
  scripts/26_step10_gene_tree_array.slurm bivalvia)

gas_trees=$(sbatch --parsable \
  --dependency=afterok:"$step10_prep" \
  --array=0-29%3 \
  --job-name=step10_tree_gastropoda \
  scripts/26_step10_gene_tree_array.slurm gastropoda)

step10_merge=$(sbatch --parsable \
  --dependency=afterany:"$biv_trees":"$gas_trees" \
  scripts/26_step10_summarise.slurm)
```

数组每个任务 4 CPUs/8 GB，两个数组最大并发均为 3。支持树使用 1,000 UFBoot + 1,000 SH-aLRT；若支持计算失败，脚本允许 ML-only 软降级并显式标记，汇总时必须统计。

**主要输出：**每家族 `trimmed.faa`、`.treefile`、`.iqtree`、`.contree`（ML-only 可能无 `.contree`），以及 `gene_tree_summary.tsv` 和 `failed_or_missing_trees.tsv`。

---

### Step 11 — Gene-Tree/Species-Tree Reconciliation

**目的：**把候选基因树与对应纲的物种树进行 duplication-loss LCA reconciliation，定位复制与丢失事件；避免直接凭基因树拓扑推断复制时间。

**主要工具：**自定义 Python DL-LCA reconciliation；duplication cost = 1，loss cost = 1。

```bash
biv_rec=$(sbatch --parsable --dependency=afterok:"$step10_merge" \
  --job-name=step11_reconcile_bivalvia \
  scripts/27_step11_reconcile_group.slurm bivalvia)

gas_rec=$(sbatch --parsable --dependency=afterok:"$step10_merge" \
  --job-name=step11_reconcile_gastropoda \
  scripts/27_step11_reconcile_group.slurm gastropoda)

step11_merge=$(sbatch --parsable \
  --dependency=afterok:"$biv_rec":"$gas_rec" \
  scripts/27_step11_combine.slurm)
```

**主要输出：**`family_reconciliation_summary.tsv`、event tables、per-family reconciled trees 和 combined summary。

完成 Step11 后，可把 Step01–Step11 的权威结果复制到中央 `result/` 目录；这是复制而非移动，不删除原始结果：

```bash
archive_01_11=$(sbatch --parsable \
  --dependency=afterok:"$step11_merge" \
  scripts/28_archive_results.slurm)
```

---

### Step 12 — Molecular Evolution Analysis

**目的：**将候选蛋白映射回 CDS，生成 codon-aware alignments，并用 HyPhy BUSTED、aBSREL 和 RELAX 检验基因/分支层面的正选择、分支特异选择及选择强度变化。

**主要软件：**HyPhy 2.5.101，加上脚本内的 CDS/蛋白映射和 codon alignment 处理。

Step12 的 `29_submit_step12.sh` 只提交环境检查和输入准备。推荐显式建立完整依赖链：

```bash
hyphy_setup=$(sbatch --parsable scripts/29_step12_setup_hyphy.slurm)
step12_prep=$(sbatch --parsable \
  --dependency=afterok:"$step11_merge" \
  scripts/29_step12_prepare.slurm)

step12_launch=$(sbatch --parsable \
  --dependency=afterok:"$hyphy_setup":"$step12_prep" \
  scripts/29_step12_launch.slurm "$step12_prep")
```

`29_step12_launch.slurm` 会验证 HyPhy 参数帮助、两个纲各 30 行 manifest，然后自动提交两个数组、汇总任务和归档任务。

**主要输出目录：**`results_or_reports/molecular_evolution/hyphy_busted_absrel_relax_v1/run_<PREP_JOB_ID>`；主要表为 family-level molecular-evolution results 和 method-specific summaries。

---

### Step 13 — Optional Transcriptomic Integration

**目的：**用单个公开 paired-end library 为部分物种的候选转录本提供“是否检测到/库内 TPM”支持；不进行 differential expression，不跨物种直接比较绝对 TPM。

**主要软件：**Salmon 1.10.2、Python。

该步骤依赖 Step12 生成的 candidate CDS，同时要求 `30_step13_transcriptomics.py` 中列出的 6 个 cleaned read pairs 已存在：

```bash
bash scripts/30_submit_step13.sh
```

wrapper 自动提交 prepare → Salmon 数组（0–5，最大并发 2）→ summarise → archive。

最终 6 个 run accessions 为 `ERR15985856`、`ERR15696446`、`SRR34876314`、`ERR12342469`、`ERR10378018` 和 `SRR28022337`。`Fragum_fragum` 的 `ERR12245519` 不可作为核心物种 `Fragum_sueziense` 的表达代理。

**主要输出：**`candidate_expression.tsv`、`species_expression_summary.tsv`、`orthogroup_species_expression.tsv`、`orthogroup_transcriptomic_integration.tsv` 和明确的 unavailable-species table。

---

### Step 14 — Integrated Candidate Prioritisation

**目的：**整合 gene-family change、功能注释、序列质量、基因树、reconciliation、分子进化和可选表达证据，生成 family 和 gene 两级的优先级排名。

**主要工具：**Python，自定义 evidence-weighted score。

确保 Step06–Step13 的 `latest/COMPLETE` 均指向目标结果后运行：

```bash
bash scripts/31_submit_step14.sh
```

family 分数的最大构成为：family change 25、functional annotation 10、sequence quality 20、gene tree 10、reconciliation 10、molecular evolution 20、expression support 5。缺失 HyPhy/RNA evidence 时按可用最大分数归一化，不能把缺失当作 0 生物学效应。Step08 enrichment 只作上下文，不重复加入分数。

**主要输出：**family ranking、gene ranking、evidence matrix、priority tiers、methods/limitations 和 source manifest。

---

### Step 15 — Sensitivity and Bias Analysis

**目的：**在 13 个替代情景下测试候选排序对蛋白组质量、重复度、注释来源、缺失证据和权重选择的敏感性，量化排名稳定性和潜在偏差。

**主要工具：**Python。

不要使用历史运行 `18428504`。最终权威逻辑对应修正后的脚本和运行 `18431248`。新重跑时显式传入新的 Step14 目录：

```bash
STEP14_RUN=/path/to/mollusc_project/results_or_reports/integrated_candidate_prioritisation/evidence_weighted_v1/latest
bash scripts/32_submit_step15.sh "$STEP14_RUN"
```

**主要输出：**scenario rankings、rank correlations、top-candidate stability、quality/annotation-source sensitivity tables 和 summary。

---

### Step 16 — Statistical Analysis and Visualisation

**目的：**汇总 proteome quality、orthogroup composition、family changes、系统发育、功能富集、验证、分子进化、表达支持、综合排名和敏感性结果，生成论文用表和图。

**主要软件：**R 4.5.2、ggplot2 4.0.3。

```bash
STEP14_RUN=/path/to/mollusc_project/results_or_reports/integrated_candidate_prioritisation/evidence_weighted_v1/latest
STEP15_RUN=/path/to/mollusc_project/results_or_reports/sensitivity_bias_analysis/ranking_sensitivity_v1/latest
bash scripts/33_submit_step16.sh "$STEP14_RUN" "$STEP15_RUN"
```

wrapper 自动提交统计/绘图任务和中央归档任务。

**主要输出：**统计 TSV、9 张 PNG、PDF 图册、figure captions、综合报告和 manifest。最终权威历史运行是 `18451973`。

---

### Step 17 — Reproducibility and Result Archiving

**目的：**冻结 Step01–Step16 的最终运行选择，生成软件/数据库版本、Slurm ledger、数据 lineage、分析决策、结果 manifest、SHA-256、敏感信息扫描和可解压验证；制作紧凑交付包。

**主要工具：**Python、tar/gzip、SHA-256 和中央结果同步脚本。

先确认 Step12–Step16 的 archive jobs 已完成，中央目录 `result/` 包含 Step01–Step16，再运行：

```bash
bash scripts/34_submit_step17.sh
```

wrapper 自动提交 `34_step17_reproducibility_archive.slurm`，成功后再运行 `34_step17_archive.slurm` 把 Step17 结果复制到中央目录。

**主要输出：**

- `FINAL_RUN_SELECTION.tsv`
- `DATA_LINEAGE.tsv`
- `ANALYSIS_DECISIONS.tsv`
- `SOFTWARE_AND_DATABASE_VERSIONS.tsv`
- `SLURM_JOB_LEDGER.psv`
- `MASTER_RESULT_MANIFEST.tsv`
- `CENTRAL_RESULT_SHA256SUMS`
- `SECURITY_SCAN.tsv`
- `ARCHIVE_VALIDATION_REPORT.tsv`
- reproducibility scripts、全部 TSV 和 visualisation compact archives。

## 三、依赖关系和可并行部分

严格依赖主链为：

```text
Step01 → Step02 → Step03 → Step05 → Step06 → Step07 → Step09 → Step10 → Step11
                                                └────→ Step08
Step10 + Step11 → Step12 → Step13
Step06–Step13 → Step14 → Step15 → Step16 → Step17
```

可安全并行的部分：

- Step02 的 Bivalvia 与 Gastropoda 类群作业并行；
- Step03 combined 与 Gastropoda restricted-MFP 并行；
- Step04 在 Step02 完成后可与后续树分析并行；
- Step05 两个纲并行；
- Step07 两个纲并行；
- Step09 两个纲并行；
- Step10 两个 30-task 数组并行，各自限制最大并发 3；
- Step11 两个纲并行；
- Step12 两个 30-task HyPhy 数组并行，各自限制最大并发 3。

不能提前的关键点：Step05 必须使用最终类群树；Step10 必须等待 Step09；Step11 必须等待 Step10；Step13 需要 Step12 的 candidate CDS；Step14 必须等待所需的 Step06–Step13 evidence；Step17 必须等待所有归档完成。

## 四、作业监控与完成判据

```bash
squeue -u HPC_USER
sacct -j JOB_ID --format=JobIDRaw,JobName,State,ExitCode,Elapsed,AllocCPUS,ReqMem,MaxRSS
```

不要只以 Slurm `COMPLETED` 判断生物信息学步骤成功。每个步骤还应核对：

1. 对应 `COMPLETE` 或 `PREPARED` 标记存在；
2. 主要 TSV/树/figure 文件非空；
3. `.err` 中无 traceback、out-of-memory、quota/no-space、missing input 或工具级 fatal error；
4. `latest` 符号链接指向本次权威运行；
5. 归档完成后原始结果目录仍存在，确认执行的是复制而不是移动。

历史权威结果路径应以 `metadata/FINAL_RUN_SELECTION.tsv` 为准，而不是根据目录日期自行选择。
