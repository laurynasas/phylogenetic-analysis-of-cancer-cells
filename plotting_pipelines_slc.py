import matplotlib
import matplotlib.pyplot as plt
import matplotlib.text as mpl_text
import numpy as np

print matplotlib.__version__


class AnyObject(object):
    def __init__(self, text, color):
        self.my_text = text
        self.my_color = color


class AnyObjectHandler(object):
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        print orig_handle
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        patch = mpl_text.Text(x=0, y=0, text=orig_handle.my_text, color=orig_handle.my_color,
                              verticalalignment=u'baseline',
                              horizontalalignment=u'left', multialignment=None,
                              fontproperties=None, rotation=0, linespacing=None,
                              rotation_mode=None)
        handlebox.add_artist(patch)
        return patch


N = 3

small_pars = (2.449489743,
              4.795831523,
              4)

medium_pars = (16.0623784,
               9.591663047)

medium_pars_isolated = (4)
medium_pars_isolated_std = (0)

large_pars = (21.88606863,
              23.83275058,
              15.06651917)

small_pars_std = (0,
                  0.,
                  0.)

medium_pars_std = (0.,
                   0.)

large_pars_std = (0.,
                  0.,
                  0.)

small_nj = (2.449489743,
            3.464101615,
            2.645751311)

medium_nj = (15.49193338,
             16.82260384,
             9.797958971)
large_nj = (15.19868415,
            24.73863375,
            22.53885534)
small_nj_std = (0,
                0.,
                0.)
medium_nj_std = (0.,
                 0.,
                 0.)
large_nj_std = (0.,
                0.,
                0.)

small_upgma = (2.449489743,
               4.795831523,
               4)
medium_upgma = (16.24807681,
                9.539392014,
                9.591663047)
large_upgma = (19.36491673,
               23.10844002,
               16.70329309)
small_upgma_std = (0,
                   0.,
                   0.)
medium_upgma_std = (0.,
                    0.,
                    0.)
large_upgma_std = (0.,
                   0.,
                   0.)


def autolabel(rects, text, extra_height=0):
    """
    Attach a text label above each bar displaying its height
    """
    for index, rect in enumerate(rects):

        height = rect.get_height()
        if extra_height != 0 and index == 2:
            extra_height = 0.5
        if extra_height != 0 and index == 0:
            extra_height = 2.5

        plt.text(rect.get_x() + rect.get_width() / 2., height + 4 + extra_height,
                 text,
                 ha='center', va='bottom')


ind = np.arange(N)  # the x locations for the groups
width = 0.35  # the width of the bars: can also be len(x) sequence
p3 = plt.bar(ind, large_pars, 0.5 * width, color='#3D6E05', yerr=large_pars_std, ecolor='b')
ind_isolated = np.array([0,2])
p2 = plt.bar(ind_isolated, medium_pars, 0.5 * width, color='#55234D', yerr=medium_pars_std, ecolor='r')
ind = np.arange(N)

p1 = plt.bar(ind, small_pars, 0.5 * width, color='#d62728', yerr=small_pars_std, ecolor='y')
p2_isolated = plt.bar(np.array([1]), medium_pars_isolated, 0.5 * width, color='#55234D', yerr=medium_pars_isolated_std, ecolor='r')

p4 = plt.bar(ind + 0.5 * width, large_nj, 0.5 * width, color='#3D6E05', yerr=large_nj_std, ecolor='b')
p5 = plt.bar(ind + 0.5 * width, medium_nj, 0.5 * width, color='#55234D', yerr=medium_nj_std, ecolor='r')
p6 = plt.bar(ind + 0.5 * width, small_nj, 0.5 * width, color='#d62728', yerr=small_nj_std, ecolor='y')

p7 = plt.bar(ind + width, large_upgma, 0.5 * width, color='#3D6E05', yerr=large_upgma_std, ecolor='b')
p8 = plt.bar(ind + width, medium_upgma, 0.5 * width, color='#55234D', yerr=medium_upgma_std, ecolor='r')
p9 = plt.bar(ind + width, small_upgma, 0.5 * width, color='#d62728', yerr=small_upgma_std, ecolor='y')

autolabel(p3, "T1")
autolabel(p4, "T2", extra_height=0)
autolabel(p7, "T3")
# autolabel(p2)
# p4 = plt.bar(ind+0.5*width, womenMeans1, 0.5*width,
#              bottom=menMeans1, yerr=womenStd)
# p5 = plt.bar(ind+0.5*width, menMeans, 0.5*width,
#              bottom=menMeans1, yerr=womenStd)
# p6 = plt.bar(ind+0.5*width, menMeans2, 0.5*width, color='#ff0000',
#              bottom=menMeans1, yerr=womenStd)


plt.ylabel('Distance between true and predicted trees')
plt.title('Performance of SLC clustering with\n different tree reconstruction methods on different-size data')
plt.xticks(ind + 0.75 * width, ('Data sets with 1% error', 'Data sets with 5% error', 'Data sets with 10% error'))
plt.yticks(np.arange(0, 30, 5))

plt.ylim(ymax=40,ymin=0)
plt.xlim(xmax=3 - 1.2 * width)

obj_0 = AnyObject("T1", "black")
obj_1 = AnyObject("T2", "black")
obj_2 = AnyObject("T3", "black")

plt.legend((p1[0], p2[0], p3[0], obj_0, obj_1, obj_2),
           ('Small data set', 'Medium data set', 'Large data set', "Parsimony", "Neighbour Joining", "UPGMA"),
           handler_map={obj_0: AnyObjectHandler(), obj_1: AnyObjectHandler(), obj_2: AnyObjectHandler()})

plt.show()
