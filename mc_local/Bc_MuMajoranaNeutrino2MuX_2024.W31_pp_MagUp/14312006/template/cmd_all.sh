echo "=================DEBUG beginning step 1====="
date


lb-run --siteroot=/cvmfs/lhcb.cern.ch/lib/ -c x86_64_v2-centos7-gcc11-opt --use=AppConfig.v3r460 --use=Gen/DecFiles.v32r41 --use=ProdConf Gauss/v56r11 gaudirun.py -T '$APPCONFIGOPTS/Gauss/Beam6800GeV-mu100-2024.W31-nu6.3.py' '$APPCONFIGOPTS/Gauss/EnableSpillover-25ns.py' '$APPCONFIGOPTS/Gauss/Run3-detector.py' '$APPCONFIGOPTS/Gauss/DataType-2024.py' '$DECFILESROOT/options/14312006.py' '$LBBCVEGPYROOT/options/BcVegPyPythia8.py' '$APPCONFIGOPTS/Gauss/G4PL_FTFP_BERT_EmOpt2.py' '$APPCONFIGOPTS/Persistency/BasketSize-10.py' '$APPCONFIGOPTS/Persistency/Compression-ZSTD-1.py' prodConf_00012345_00006789_1.py


echo "=================DEBUG beginning step 2====="
date

lb-run --siteroot=/cvmfs/lhcb.cern.ch/lib/ -c x86_64_v2-el9-gcc13+detdesc-opt --use=AppConfig.v3r460 --use=ProdConf Boole/v47r0p1 gaudirun.py -T '$APPCONFIGOPTS/Boole/Default.py' '$APPCONFIGOPTS/Boole/EnableSpillover.py' '$APPCONFIGOPTS/Boole/Boole-Upgrade-Baseline-20200616.py' '$APPCONFIGOPTS/Boole/Upgrade-RichMaPMT-NoSpilloverDigi.py' '$APPCONFIGOPTS/Boole/Boole-Upgrade-IntegratedLumi-0fb.py' '$APPCONFIGOPTS/Boole/Run3-VP-NoSpillOver.py' '$APPCONFIGOPTS/Persistency/BasketSize-10.py' '$APPCONFIGOPTS/Boole/MuonLowE-Bkg-G4.py' '$APPCONFIGOPTS/Persistency/Compression-ZSTD-1.py' prodConf_00012345_00006789_2.py


echo "=================DEBUG beginning step 3====="
date

lb-run --siteroot=/cvmfs/lhcb.cern.ch/lib/ -c x86_64_v2-el9-gcc13+detdesc-opt Moore/v55r11p10 lbexec Moore.production:hlt1 lbexec_options_00012345_00006789_3.yaml -- --sequence=hlt1_pp_matching_no_ut_1000KHz --flagging


echo "=================DEBUG beginning step 4====="
date

lb-run --siteroot=/cvmfs/lhcb.cern.ch/lib/ -c x86_64_v2-el9-gcc13+detdesc-opt Moore/v55r11p10 lbexec Moore.production:hlt2 lbexec_options_00012345_00006789_4.yaml -- --velo-source=VPRetinaCluster --settings=hlt2_pp_2024 --flagging


echo "=================DEBUG beginning step 5====="
date

lb-run --siteroot=/cvmfs/lhcb.cern.ch/lib/ -c x86_64_v2-el9-gcc13+detdesc-opt Moore/v55r11p10 lbexec Moore.production:spruce lbexec_options_00012345_00006789_5.yaml -- --velo-source=VPRetinaCluster --settings=Sprucing_production_physics_pp_Collision24c2 --flagging


echo "=================DEBUG beginning step 6====="
date

lb-run --siteroot=/cvmfs/lhcb.cern.ch/lib/ -c x86_64_v3-el9-gcc13+detdesc-opt+g --use=AppConfig.v3r460 LHCb/v58r6 lbexec GaudiConf.mergeDST:dst lbexec_options_00012345_00006789_6.yaml -- --merge-genfsr

echo "=================DEBUG end step 6====="
date
