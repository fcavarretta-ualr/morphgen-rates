import neurom as nm
import matplotlib.pyplot as plt
from neurom.core.types import NeuriteType
import numpy as np
from pathlib import Path
from collections import defaultdict
import matplotlib.pyplot as plt
import pandas as pd
import os

def _load_blacklist(directory: str, filename: str = "blacklist.txt") -> list[str]:
    path = Path(directory) / filename
    if path.exists() and path.is_file():
        with path.open("r", encoding="utf-8") as f:
            return [line.strip() for line in f if line.strip()]
    return []

def list_morphologies(directory):
    bl = _load_blacklist(directory)
    return [os.path.join(directory, f) for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f)) and f.endswith(".swc") and f.replace('.CNG.swc', '') not in bl]

def process_morphology(m, dx=50.0, section_type=4):
    collections = defaultdict(int)
    segments = []

    tot = 0.0
    bw = 0.0
    bp = 0
    for s in nm.iter_sections(m):
        if s.type == section_type:
            if s.parent is None:
                print('first point', s.points[0])
                collections[0] += 1
                
            if len(s.children) == 2:
                bp += 1
            elif len(s.children) > 2 or len(s.children) == 1:
                print(len(s.children))
                raise Exception()
                  
            for i in range(1, len(s.points)):
                slen = abs(np.linalg.norm(s.points[i - 1][:3]) - np.linalg.norm(s.points[i][:3]))
                tot += slen
                
                i0 = int(np.linalg.norm(s.points[i - 1][:3]) / dx)
                i1 = int(np.linalg.norm(s.points[i][:3]) / dx)
                segments.append(np.linalg.norm(s.points[i][:3] - s.points[i - 1][:3]))
                if i1 != i0:
                    if i0 > i1:
                        print('Warning there is a segment that move backward: ', slen)
                        bw += slen
                        collections[i0 * dx] += 1
                    elif i0 < i1:
                        if i1 - i0 > 1:
                            print('Warning there is a segment that skip spatial intervals')
                        else:
                            # count intersection
                            collections[i1 * dx] += 1
    l = bw / tot
    print('segments (min/max):', min(segments), '/', max(segments), ' l=', l)
    print()
    return collections, l, bp


def process_morphologies(dirname, xmax, filename, filename2, dx=50.0, section_type=4):
    xp = np.arange(0, xmax, dx)
    yp_all = []
    l_all = []
    bp_all = []
    for fname in list_morphologies(dirname):
        try:
            p = nm.load_morphology(fname)
            collections, l, bp = process_morphology(p, section_type=section_type)
            tmp = np.array(sorted(collections.items()))
            yp = np.interp(xp, tmp[:, 0], tmp[:, 1])
            yp[xp > tmp[-1, 0]] = 0
            yp_all.append(yp)
            plt.plot(xp, yp, color='gray')
            l_all.append(l)
            bp_all.append(bp)
        except:
            print('\t\t\t\t\t\t\tWarning! Morphology ', fname, ' cannot be processed')
    plt.errorbar(xp, np.mean(yp_all, axis=0), yerr=np.std(yp_all, axis=0), color='black')
    plt.show()
    df = pd.DataFrame()
    df['Distance'] = xp
    for i, c in enumerate(yp_all):
        df['Count%d' % i] = c
    df.set_index('Distance').to_csv(filename)
    print('l +/- SE:', np.mean(l_all), np.std(l_all) / np.sqrt(len(l_all)))
    print('branch points:', bp_all)
    print('branch points +/- SD:', np.mean(bp_all), np.std(bp_all))
    pd.DataFrame(bp_all).to_csv(filename2)

if __name__ == '__main__':
    print('\nPYR')
    process_morphologies('morphologies/bathellier/PYR', 800, 'pyr_apical_sholl_plot.txt', 'pyr_bifurcations.txt')
    print('\nSL')
    process_morphologies('morphologies/bathellier/SL', 500, 'sl_apical_sholl_plot.txt', 'sl_bifurcations.txt')
    print('\nNeocortex')
    process_morphologies('morphologies/markram_dataset', 1400, 'neocortex_apical_sholl_plot.txt', 'neocortex_bifurcations.txt')
    print('\nMitral')
    process_morphologies('morphologies/mitral', 1400, 'mitral_sholl_plot.txt', 'mitral_bifurcations.txt', section_type=3)
    print('\nTufted')
    process_morphologies('morphologies/tufted', 1400, 'tufted_sholl_plot.txt', 'tufted_bifurcations.txt', section_type=3)
