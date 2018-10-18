import json
import math
try:
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
except:
    pass

#
# Statistics
#

class stats(object):
    tags = {} # dict from tag to ids
    values = {} # dict from stat to dict from id to value
    references = {} # dict from str to values

    @staticmethod
    def add(id, tags, **args):
        for k,v in args.items():
            if id in stats.values.setdefault(k, {}):
                raise Exception("Repeated id!")
            stats.values[k][id] = v
        for t in tags:
            stats.tags.setdefault(t, []).append(id)
        stats.tags.setdefault("*", []).append(id) 

    @staticmethod
    def show():
        def mean(x):
            return sum(x)/float(len(x)) if len(x) > 0 else None, len(x)

        result = {}
        for tag in sorted(stats.tags.keys()):
            ids = frozenset(stats.tags[tag])
            for k in stats.values.keys():
                v = [x for id,x in stats.values[k].items() if id in ids]
                #yield (tag, k, len(v), "sum", sum(v))
                result.setdefault(tag, {}).setdefault("itself", {})["sum " + k] = sum(v)
                result.setdefault(tag, {}).setdefault("itself", {})["n"] = len(v)
                for rname, rvalue in stats.references.items():
                    if k not in rvalue: continue
                    mu, n = mean([ (rvalue[k][id] - x)/float(rvalue[k][id]) for id,x in stats.values[k].items() if id in ids and id in rvalue[k] ])
                    var, _ = mean([ ((rvalue[k][id] - x)/float(rvalue[k][id]) - mu)**2 for id,x in stats.values[k].items() if id in ids and id in rvalue[k] ])
                    #yield (tag, k, n, "rel " + rname, mu, math.sqrt(var) if var is not None else None)
                    result.setdefault(tag, {}).setdefault("rel " + rname, {})["n"] = n
                    result.setdefault(tag, {}).setdefault("rel " + rname, {})["avg " + k] = mu
                    result.setdefault(tag, {}).setdefault("rel " + rname, {})["std error " + k] = math.sqrt(var) if var is not None else None
        return result

    @staticmethod
    def compare(r0, r1, key="mv"):
        ids = [ id for id in stats.references[r0][key].keys() if id in stats.references[r1][key] ]
        def onpick3(event):
            ind = event.ind
            print('onpick3 scatter:', ind, [ids[i] for i in ind])

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel(r0 + ' ' + key)
        ax.set_ylabel(r1 + ' ' + key)
        ax.set_title(key)

        #ax = fig.add_subplot(111, projection='3d')
        x = [ stats.references[r0][key][id] for id in ids ]
        y = [ stats.references[r1][key][id] for id in ids ]
        col = ax.scatter(x, y, picker=True)
        fig.canvas.mpl_connect('pick_event', onpick3)
        plt.show()

                   
    @staticmethod
    def dump(filename):
        f = open(filename, "w")
        json.dump(stats.values, f)
        f.close()

    @staticmethod
    def load_reference(filename):
        f = open(filename, "r")
        stats.references[filename] = json.load(f)
        f.close()


if __name__ == "__main__":
    stats.load_reference("tests-master.json")
    stats.load_reference("tests.json")
    print(stats.show())
    #stats.compare("tests-master.json", "tests.json")
   
