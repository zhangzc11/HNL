DECFILESROOT_custom=/home/zhicaiz/GaussDev_v56r11/Gen/DecFiles

lb-run Gauss/v56r11 gaudirun.py \
    '$APPCONFIGOPTS/Gauss/Beam6800GeV-mu100-2024.W31-nu6.3.py' \
    '$APPCONFIGOPTS/Gauss/DataType-2024.py' \
    '$GAUSSOPTS/GenStandAlone.py' \
    $DECFILESROOT_custom/options/11372001.py \
    '$LBPYTHIA8ROOT/options/Pythia8.py' \
    Gauss-Job.py
