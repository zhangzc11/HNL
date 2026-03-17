directory = "myTupleOS1J"
mc_file = "/pool1/lhcb/zhicaiz/mc/2024/HNL_Wmumuqq/ntuple/ntuple_20260228/42912010_5GeVCtau0ps_options_20260228_111712_skim.root"
data_file = "/pool1/lhcb/zhicaiz/data/Collision24/ntuple_HNL_Wmumuqq/data_MagUp_Sprucing24r1_all_disk_options_20260228_011229_skim.root"

sort: ./csv/branch_ks.csv
	python3 weight_sort.py ${directory} ${mc_file} ${data_file}

./csv/branch_ks.csv: ./ks_andeson.py 
	python3 ks_andeson.py ${directory} ${mc_file} ${data_file}

.PHONY: clean
clean:
	rm -rf ./png/* ./csv/*