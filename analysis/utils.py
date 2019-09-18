import os, ast, json
import glob
import numpy as np
import uproot
import uproot_methods
import pandas as pd
from scipy.stats import poisson

def LoadData(indir, years=[2016,2017,2018], do_phi=True, apply_mcand_mass_cut=True, apply_leppho_mass_cut=True, do_blind=True):
    # returns as dict data[YEAR][SAMPLE] = dataframe

    dfs = {}

    for y in years:
        dfs[y] = {}
        for fname in glob.glob(os.path.join(indir, str(y), "*.root")):
            print "Loading from file", fname
            samp = fname.split("/")[-1].split(".")[0]
            isSignal = samp.lower().startswith("wh")

            fid = uproot.open(fname)
            tree = fid["Events"]
            branches_to_read = tree.keys()
            to_remove = ["mcand_dR_genMeson"]
            for b in branches_to_read:
                # if b.startswith("evt"):
                    # to_remove.append(b)
                if not isSignal and b.startswith("gen"):
                    to_remove.append(b)
            for b in to_remove:
                branches_to_read.remove(b)

            df = tree.pandas.df(branches_to_read)

            # in case it's not skimmed
            df.drop(df[(df.nmcand < 1) | (df.nlep != 1) | (df.nphoton < 1)].index, inplace=True)

            # drop events failing MET filters
            df.drop(df[~df.passFilters].index, inplace=True)

            # for data, check golden JSON and triggers
            if "data" in samp.lower():
                df.drop(df[~df.isGolden].index, inplace=True)
                if "SingleMuon" in samp:
                    df.drop(df[(df.HLT_SingleMu == 0) | ((df.lepton_pdgId != 13) & (df.lepton_pdgId != -13))].index, inplace=True)
                elif "EGamma" in samp or "SingleElectron" in samp:
                    df.drop(df[(df.HLT_SingleEl == 0) | ((df.lepton_pdgId != 11) & (df.lepton_pdgId != -11))].index, inplace=True)
                else:
                    raise Exception("Unknown data sample name: "+samp)

            # lepton pT cuts from trigger
            df.drop(df[((df.lepton_pdgId == 11) | (df.lepton_pdgId == -11)) & (df.lepton_pt < 35)].index, inplace=True)
            df.drop(df[((df.lepton_pdgId == 13) | (df.lepton_pdgId == -13)) & (df.lepton_pt < 30)].index, inplace=True)

            # for signal, skim based on decay mode
            if isSignal:
                if do_phi:
                    df.drop(df[df.genHMeson_pdgId != 333].index, inplace=True)
                else:
                    df.drop(df[df.genHMeson_pdgId != 113].index, inplace=True)

            # remove phi/rho-specific branches
            if do_phi:
                df.drop("mcand_mass_pion", axis=1, inplace=True)
                df.rename(columns={"mcand_mass_kaon":"mcand_mass"}, inplace=True)
            else:
                df.drop("mcand_mass_kaon", axis=1, inplace=True)
                df.rename(columns={"mcand_mass_pion":"mcand_mass"}, inplace=True)


            # if multiple mesons candidates, select the closest to phi/rho mass
            # main_index gets the index of each "event" (level 0) that we haven't dropped
            main_index = sorted(set(df.index.levels[0][df.index.labels[0]]))
            # best_cand gives the index of the best meson candidate within each event (level 1)
            # So if event 94 has 3 meson candidates, and #1 is the best, then the index we want is (94,1)
            best_cand = df.best_phiCand_idx if do_phi else df.best_rhoCand_idx
            df = df.loc[zip(main_index, best_cand.loc[main_index,0])]

            # cut out events not in phi/rho mass peak
            if apply_mcand_mass_cut:
                if do_phi:
                    df.drop(df[(df.mcand_mass < 1.00) | (df.mcand_mass > 1.04)].index, inplace=True)
                else:
                    df.drop(df[(df.mcand_mass < 0.50) | (df.mcand_mass > 1.00)].index, inplace=True)

            # blind out higgs mass region if do_blind
            if "data" in samp.lower() and do_blind:
                df.drop(df[(df.mcand_photon_mass > 120) & (df.mcand_photon_mass < 130)].index, inplace=True)

            # re-index
            df.index = np.arange(df.shape[0])

            # compute lepton_photon mass so we can cut on it
            photon_p4 = uproot_methods.TLorentzVectorArray.from_ptetaphim(df.photon_pt, df.photon_eta, df.photon_phi, 0)
            lepton_p4 = uproot_methods.TLorentzVectorArray.from_ptetaphim(df.lepton_pt, df.lepton_eta, df.lepton_phi, 0)
            df["lepton_photon_mass"] = (photon_p4 + lepton_p4).mass

            if apply_leppho_mass_cut:
                df.drop(df[(df.lepton_photon_mass > 80) & (df.lepton_photon_mass < 95)].index, inplace=True)
            # don't actually want this branch (some have bad values, causes problems later)
            df.drop(columns=["lepton_photon_mass"], inplace=True)

            dfs[y][samp] = df

    return dfs

def fit_pow(x, par):
    N = par[0]
    a = par[1]
    return N * x**a / ((180**(a+1) - 100**(a+1)) / (a+1))

def fit_exp(x, par):
    N = par[0]
    a = par[1]
    return N * np.exp(-(x-100)/a) / (a*(1-np.exp(-80/a)))

def fit_quad(x, par):
    N = par[0]
    a = par[1]
    b = par[2]
    return N * (a*(x-100)**2 + b*(x-100) + 1) / (80.0/3*(6400*a + 120*b + 3))

def fit_crystalBall(x, par):
    norm = par[0]
    mean = par[1]
    sigma = par[2]
    alpha = par[3]
    n = par[4]

    t = (x-mean)/sigma
    if alpha < 0:
        t *= -1

    t = np.array(t)
    isfloat = len(t.shape)==0
    if isfloat:
        t = t.reshape(1,)

    res = norm * np.exp(-0.5*t*t)
    idx = t < -abs(alpha)
    a = (n/abs(alpha))**n * np.exp(-0.5*abs(alpha)**2)
    b = n/abs(alpha) - abs(alpha)
    res[idx] = norm * a / (b-t[idx])**n

    if isfloat:
        return res[0]
    else:
        return res
    
# some global parameters for controlling fits
FIT_BLINDED = True
FIT_BINS = np.arange(100,180.1,1)
FIT_BLINDBINS = (20,30)
FIT_BINCONTENTS = None  # set this before fitting
FIT_FUNC = fit_pow

def binned_likelihood(par):
    nll = 0.0
    for i in range(FIT_BINS.size-1):
        if FIT_BLINDED and i >= FIT_BLINDBINS[0] and i < FIT_BLINDBINS[1]:
            continue
        x = 0.5*(FIT_BINS[i]+FIT_BINS[i+1])
        mu = FIT_FUNC(x, par)
        k = FIT_BINCONTENTS[i]
        nll -= poisson.logpmf(k, mu)
    return nll

# chi-squared function, using sqrt(contents) as the error in each bin
def chisquared(par):
    cs = 0.0
    for i in range(FIT_BINS.size-1):
        if FIT_BLINDED and i >= FIT_BLINDBINS[0] and i < FIT_BLINDBINS[1]:
            continue
        x = 0.5*(FIT_BINS[i]+FIT_BINS[i+1])
        mu = FIT_FUNC(x, par)
        k = FIT_BINCONTENTS[i]        
        if k>0:
            cs += abs(mu - k)**2 / k
    return cs

def BDTtoJSON(fname, bst, features, labels=[]):
    """Return JSON of BDT model (Written by Nick Amin)
    Parameters
    ----------
    fname : str
        Desired name of output file
    bst : xgboost.core.Booster
        Trained BDT
    features : list[str]
        List of features used to train BDT
    labels : list[str], optional
        Class labels (default [])
    """
    buff = "[\n"
    trees = bst.get_dump()
    ntrees = len(trees)
    for itree,tree in enumerate(trees):
        prev_depth = 0
        depth = 0
        for line in tree.splitlines():
            depth = line.count("\t")
            nascending = prev_depth - depth
            (depth == prev_depth-1)
            # print ascending, depth, prev_depth
            prev_depth = depth
            parts = line.strip().split()
            padding = "    "*depth
            for iasc in range(nascending):
                buff += "{padding}]}},\n".format(padding="    "*(depth-iasc+1))
            if len(parts) == 1:  # leaf
                nodeid = int(parts[0].split(":")[0])
                leaf = float(parts[0].split("=")[-1])
                # print "leaf: ",depth,nodeid,val
                buff += """{padding}{{ "nodeid": {nodeid}, "leaf": {leaf} }},\n""".format(
                        padding=padding,
                        nodeid=nodeid,
                        leaf=leaf,
                        )
            else:
                nodeid = int(parts[0].split(":")[0])
                split, split_condition = parts[0].split(":")[1].replace("[","").replace("]","").split("<")
                split_condition = float(split_condition)
                yes, no, missing = map(lambda x:int(x.split("=")[-1]), parts[1].split(","))
                buff += """{padding}{{ "nodeid": {nodeid}, "depth": {depth}, "split": "{split}", "split_condition": {split_condition}, "yes": {yes}, "no": {no}, "missing": {missing}, "children": [\n""".format(
                        padding=padding,
                        nodeid=nodeid,
                        depth=depth,
                        split=split,
                        split_condition=split_condition,
                        yes=yes,
                        no=no,
                        missing=missing,
                        )
        for i in range(depth):
            padding = "    "*(max(depth-1,0))
            if i == 0:
                buff += "{padding}]}}".format(padding=padding)
            else:
                buff += "\n{padding}]}}".format(padding=padding)
            depth -= 1
        if itree != len(trees)-1:
            buff += ",\n"
    buff += "\n]"
    # print buff
    to_dump = {
            "trees": list(ast.literal_eval(buff)),
            "features": features,
            "labels": map(int,np.array(labels).tolist()), # numpy array not json serializable
            }
    with open(fname, "w") as fhout:
        json.dump(to_dump,fhout,indent=2)

def JSONtoC(fname_in, fname_out=None):
    """Return C function that mimics BDT output (Written by Nick Amin)
    Parameters
    ----------
    fname_in : str
        Name of input JSON file
    fname_out : str, optional
        Desired name of output C file
    """
    with open(fname_in, "r") as fhin:
        data = json.loads(fhin.read())
        trees = data["trees"]
        features = data["features"]
        labels = data["labels"]
    def get_leaf(entry, depth=0):
        if "leaf" in entry: 
            return entry["leaf"]
        splitvar = entry["split"]
        splitval = entry["split_condition"]
        yesnode = [c for c in entry["children"] if c["nodeid"] == entry["yes"]][0]
        nonode = [c for c in entry["children"] if c["nodeid"] == entry["no"]][0]
        # if "mcand" in splitvar or "vhAngles" in splitvar:
            # splitvar += "[best_phiCand_idx]"
        return "({} < {} ? {} : {})".format(splitvar, splitval, get_leaf(yesnode, depth=depth+1), get_leaf(nonode, depth=depth+1))
    buff = ""
    multi = False
    if len(labels) > 0:
        multi = True
        ntrees_per_class = len(trees) // len(labels)
        nclasses = len(labels)
    if multi:
        buff += "std::vector<float> get_prediction({}) {{\n".format(",".join(map(lambda x: "float {}".format(x), features)))
        for ic in labels:
            buff += "  float w_{} = 0.;\n".format(ic)
        for itree,j in enumerate(trees):
            iclass = int(labels[itree % nclasses])
            buff += "  w_{} += {};\n".format(iclass,get_leaf(j))
        buff += "  float w_sum = {};\n".format("+".join("exp(w_{})".format(ic) for ic in labels))
        for ic in labels:
            buff += "  w_{0} = exp(w_{0})/w_sum;\n".format(ic)
        buff += "  return {{ {} }};\n".format(",".join("w_{}".format(ic) for ic in labels))
    else:
        buff += "float get_prediction({}) {{\n".format(",".join(map(lambda x: "float {}".format(x), features)))
        buff += "  float w = 0.;\n"
        for itree,j in enumerate(trees):
            buff += "  w += {};\n".format(get_leaf(j))
            # buff += "  {}+".format(get_leaf(j))
        buff += "  return 1.0/(1.0+exp(-1.0*w));\n"
    buff += "}"
    buff = buff.strip().strip("+")
    buff += "\n"
    if fname_out:
        with open(fname_out, "w") as fhout:
            fhout.write(buff)
    else:
        return buff

