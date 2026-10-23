# REGSTA-script

**Single-cell-resolution spatial transcriptomics for plant tissue: cell-wall segmentation → cell-level count matrix → Seurat analysis.**

Scripts and a trained segmentation model used in

> Yao J†, Chu Q†, Guo X†, Shao W†, et al., Liu H\*, Fan L\*. **Spatiotemporal transcriptomic landscape of rice embryonic cells during seed germination.** *Developmental Cell* (2024). [doi:10.1016/j.devcel.2024.05.016](https://doi.org/10.1016/j.devcel.2024.05.016)

Plant cells have walls and vary widely in size, so bin-based spatial transcriptomics (e.g. Stereo-seq bin50) does not map onto real cells. This pipeline segments individual cells from a stained tissue image with a **custom Cellpose model trained on plant cell walls**, assigns every spatial barcode (DNB) to the cell it falls in, and produces a standard 10x-style count matrix that can be analysed as single-cell data.

> 面向植物组织的单细胞分辨率空间转录组流程：用针对植物细胞壁训练的 Cellpose 模型分割单个细胞，将 Stereo-seq 的空间条形码归属到细胞，生成可直接用 Seurat 分析的细胞级表达矩阵。配套 *Developmental Cell* 2024 水稻胚萌发时空转录组论文。

## Pipeline

| Step | File | Input → Output |
|---|---|---|
| 1. Cell segmentation | `01cellposePrediction.txt` | tissue image (`data4test/img/`) → cell outlines (`txt_outlines/*_cp_outlines.txt`) using `model/cellpose_*` |
| 2. Barcodes → cells | `02outline2scResGemCsv.py` | outlines + image + Stereo-seq `.gem` → `*_scResGem.csv` (each DNB labelled with its cell ID) |
| 3. Build matrix | `03scResGemCsv2mtx.sh` | `*_scResGem.csv` → MatrixMarket folder (`matrix4test/`: `matrix.mtx`, `genes.tsv`, `barcodes.tsv`) + cell-centroid table |
| 4. Analysis | `04scST_SCT4test.R` | matrix + centroids → QC (mt/pt genes), SCTransform, clustering, spatial plots (`scSTplot/`) |

`colorTools.py` is a small RGB⇄int helper used by step 2. `data4test/` holds a down-sized rice-embryo example so every step can be run end-to-end.

## Quick start

```bash
# 1. Segment cells (Cellpose ≥ 2.0; see 01cellposePrediction.txt for the conda env)
python -m cellpose --dir ./data4test/img --pretrained_model ./model/cellpose_residual_on_style_on_concatenation_off_CP_1_2022_11_26_02_48_20.469512 \
  --diameter 0 --flow_threshold 5 --exclude_on_edges --save_outlines --save_txt --in_folders --verbose

# 2. Assign DNBs to cells  (edit the paths at the top of the script)
python 02outline2scResGemCsv.py

# 3. Write a 10x-style matrix
bash 03scResGemCsv2mtx.sh

# 4. Analyse in R (Seurat ≥ 4)
Rscript 04scST_SCT4test.R
```

**Dependencies.** Python 3.8+: `cellpose`, `opencv-python`, `pillow`, `numpy`, `pandas`, `matplotlib`, `tqdm`. R ≥ 4.1: `Seurat`, `SeuratObject`, `dplyr`, `ggplot2`, `ggpubr`, `ggthemes`, `Matrix`.

## Model

`model/cellpose_residual_on_style_on_concatenation_off_CP_1_2022_11_26_02_48_20.469512` is a Cellpose model fine-tuned on annotated rice embryo cell-wall images. It may transfer to other plant tissues with visible cell walls; if your staining or tissue differs, fine-tune it further on a few annotated slices.

## Citation

If you use this code or model, please cite the *Developmental Cell* paper above.

## License

MIT © 2022–2026 Jie Yao
