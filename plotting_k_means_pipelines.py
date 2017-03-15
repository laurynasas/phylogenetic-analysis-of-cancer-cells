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

small_pars = (0,
              3.820418676,
              0.7656854249)

medium_pars = (4.636443966,
                   8.322263063,
               7.483760659)

large_pars = (21.44938823,
              16.72148201,
              8.190360323)

small_pars_std = (0,
                  0.7102120374,
                  0.5134065547)

medium_pars_std = (0.3215390309,
                   0.8601962076,
                   1.730123233)

large_pars_std = (3.581025624,
                  4.12214014,
                  2.901378599)

small_nj = (0,
            2.950890227,
            1.053516943)

medium_nj = (4.637803989,
             9.604848754,
             6.161548599)
large_nj = (16.13111997,
            13.13300496,
            7.623047086)
small_nj_std = (0,
                0.6262961494,
                0.7000728901)
medium_nj_std = (0.3012875069,
                 3.232163427,
                 1.932697301)
large_nj_std = (5.47603584,
                3.744887288,
                2.255915142)

small_upgma = (0,
               3.547727023,
               0.9071067812)
medium_upgma = (4.910954446,
                7.629331301,
                7.020508834)
large_upgma = (23.52050676,
               15.55883612,
               7.872690889)
small_upgma_std = (0,
                   0.4622044659,
                   0.759708686)
medium_upgma_std = (0.2872741359,
                    0.6271394595,
                    1.764215324)
large_upgma_std = (2.736012011,
                   2.867510892,
                   3.036566839)


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
p2 = plt.bar(ind, medium_pars, 0.5 * width, color='#55234D', yerr=medium_pars_std, ecolor='r')
p1 = plt.bar(ind, small_pars, 0.5 * width, color='#d62728', yerr=small_pars_std, ecolor='y')

p4 = plt.bar(ind + 0.5 * width, large_nj, 0.5 * width, color='#3D6E05', yerr=large_nj_std, ecolor='b')
p5 = plt.bar(ind + 0.5 * width, medium_nj, 0.5 * width, color='#55234D', yerr=medium_nj_std, ecolor='r')
p6 = plt.bar(ind + 0.5 * width, small_nj, 0.5 * width, color='#d62728', yerr=small_nj_std, ecolor='y')

p7 = plt.bar(ind + width, large_upgma, 0.5 * width, color='#3D6E05', yerr=large_upgma_std, ecolor='b')
p8 = plt.bar(ind + width, medium_upgma, 0.5 * width, color='#55234D', yerr=medium_upgma_std, ecolor='r')
p9 = plt.bar(ind + width, small_upgma, 0.5 * width, color='#d62728', yerr=small_upgma_std, ecolor='y')

autolabel(p3, "T1")
autolabel(p4, "T2", extra_height=2)
autolabel(p7, "T3")
# autolabel(p2)
# p4 = plt.bar(ind+0.5*width, womenMeans1, 0.5*width,
#              bottom=menMeans1, yerr=womenStd)
# p5 = plt.bar(ind+0.5*width, menMeans, 0.5*width,
#              bottom=menMeans1, yerr=womenStd)
# p6 = plt.bar(ind+0.5*width, menMeans2, 0.5*width, color='#ff0000',
#              bottom=menMeans1, yerr=womenStd)


plt.ylabel('Distance between true and predicted trees')
plt.title('Performance of k-means clustering with\n different tree reconstruction methods on different-size data')
plt.xticks(ind + 0.75 * width, ('Data sets with 1% error', 'Data sets with 5% error', 'Data sets with 10% error'))
plt.yticks(np.arange(0, 40, 5))

plt.ylim(ymax=40,ymin=0)

plt.xlim(xmax=3 - 1.2 * width)

obj_0 = AnyObject("T1", "black")
obj_1 = AnyObject("T2", "black")
obj_2 = AnyObject("T3", "black")

plt.legend((p1[0], p2[0], p3[0], obj_0, obj_1, obj_2),
           ('Small data set', 'Medium data set', 'Large data set', "Parsimony", "Neighbour Joining", "UPGMA"),
           handler_map={obj_0: AnyObjectHandler(), obj_1: AnyObjectHandler(), obj_2: AnyObjectHandler()})

plt.show()
