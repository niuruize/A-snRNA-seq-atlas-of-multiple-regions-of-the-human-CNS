#scFEA
#下载安装
#git clone https://github.com/changwn/scFEA #下载scFEA
cd scFEA
conda install --file requirements
conda install pytorch torchvision -c pytorch
pip install --user magic-impute


%%bash
cd scFEA
python src/scFEA.py --data_dir data --input_dir input \
                    --test_file Melissa_full.csv \
                    --moduleGene_file module_gene_m168.csv \
                    --stoichiometry_matrix cmMat_c70_m168.csv \

python src/scFEA.py --data_dir data --input_dir input \
                    --test_file Melissa_full.csv \
                    --moduleGene_file module_gene_m168.csv \
                    --stoichiometry_matrix cmMat_c70_m168.csv
                    --output_flux_file output/Melissa_flux.csv \
                    --output_balance_file output/Melissa_balance.csv

python src/scFEA.py --data_dir data --input_dir input \
                    --test_file EC_seurat_scFEA_all_data.csv \
                    --moduleGene_file module_gene_m168.csv \
                    --stoichiometry_matrix cmMat_c70_m168.csv
                    
#最多10000个细胞
python src/scFEA.py --data_dir data --input_dir input \
                    --test_file EC_10000_data.csv \
                    --moduleGene_file module_gene_m168.csv \
                    --stoichiometry_matrix cmMat_c70_m168.csv \
                    --output_flux_file output/EC_10000_flux.csv \
                    --output_balance_file output/EC_10000_balance.csv

%%bash
python src/scFEA.py --data_dir data --input_dir input \
                    --test_file Seurat_geneExpr.csv \
                    --moduleGene_file module_gene_glutaminolysis1_m23.csv \
                    --stoichiometry_matrix cmMat_glutaminolysis1_c17_m23.csv \
                    --cName_file cName_glutaminolysis1_c17_m23.csv \
                    --output_flux_file output/Seurat_gluta_flux.csv \
                    --output_balance_file output/Seurat_gluta_balance.csv





#https://github.com/Dulab2020/singlecell-analysis/blob/master/Analysis/scFEA_analysis.py
import pandas as pd

## get data
data = pd.read_csv('./output/mouse_flux.csv')
meta = pd.read_csv("./cell_id.csv")

data['Unnamed: 0'] = meta['x']
cluster_mean = data.groupby('Unnamed: 0').mean()


data['cluster'] = meta['cell_type']
data['pid'] = meta['pid']

## mean expression for each celltype
cluster_mean = data.groupby('cluster').mean()
cluster_mean.to_csv('./cluster_mean.csv')

## mean expression for (celltype, patient)
bulk_mean = data.groupby(['pid', 'cluster']).mean()
bulk_mean.to_csv("./bulk_mean.csv")



