# Recipe: take a cut along a straight line of a k-converted map with the viewer's "Arbitrary cut".
# Example: a slice at fixed kx running along ky (the same position as a measured cut).
import numpy as np
MAP_ROW = "FineMap_s_nor [k] [E-Ef]"        # name of the k-converted map row
KX = 0.644                                    # 1/A, position of the line
NAME = "Cut_from_map [k] [E-Ef]"
key = {r["name"]: k for k, r in A.loaded_items.items()}[MAP_ROW]
m = A.load_dataset(key).scan
ky = np.asarray(m.k, float)
cont = A.map_contour(key)                     # the map's contour window (not shown if offscreen)
out = cont.run_arbitrary_cut([(KX, float(ky.min())), (KX, float(ky.max()))], name=NAME, open_window=False)
print("created", [o.source_label for o in out], flush=True)
# NOTE: the new row's x axis is the distance along the path starting at 0; add ky.min() to get ky.
