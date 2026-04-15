from options_Zmumu import make_algs_Zmumu
from options_Wmunu import make_algs_Wmunu
from DaVinci import make_config, Options

def main(options: Options):

    algs = make_algs_Zmumu(options, "myTupleZmumu")
    algs.update(make_algs_Wmunu(options, "myTupleWmunu"))

    return make_config(options, algs)
