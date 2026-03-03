from ProdConf import ProdConf

ProdConf(
    Application='Boole',
    AppVersion='v47r0p1',
    InputFiles=['LFN:00012345_00006789_1.sim'],
    OutputFilePrefix='00012345_00006789_2',
    OutputFileTypes=['digi'],
    XMLSummaryFile='summaryBoole_00012345_00006789_2.xml',
    XMLFileCatalog='pool_xml_catalog.xml',
    DDDBTag='2024-v00.04',
    CondDBTag='sim10-2024.W31-v00.00-mu100',
    NOfEvents=-1,
    TCK='',
    FirstEventNumber=0,
    NThreads=1,
)
