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

small_k_means = (1,
                 0.7738431781,
                 0.7277936963)

medium_k_means = (0.9699099121,
                  0.5306479887,
                  0.4829900944)

large_k_means = (0.908006128,
                 0.810535397,
                 0.5355520363)


small_slc = (1,
             0.7656571194,
             0.5995040298)

medium_slc = (0.99,
              0.7489166563,
              0.4778077295)
large_slc = (0.8282598476,
             0.7108345146,
             0.7502722354)


small_bmm = (0.9818337989,
             0.8779685748,
             0.6982171485)
medium_bmm = (0.840118066,
              0.6400818059,
              0.4579926847)
large_bmm = (0.8554337059,
             0.7510225512,
             0.6766839868)


def autolabel(rects, text, extra_height=0):
    """
    Attach a text label above each bar displaying its height
    """
    for index, rect in enumerate(rects):

        height = rect.get_height()
        if extra_height == 0 and index == 1:
            extra_height = 0.1
        if height >0.7 and height<0.8 and index ==2:
            extra_height = -0.1
        elif height >0.6 and height<0.7 and index ==2:
            extra_height = -0.05
        else:
            extra_height  = 0.1
        if extra_height == 0 and index == 0:
            extra_height = 0.1

        plt.text(rect.get_x() + rect.get_width() / 2., height + 0.1 + extra_height,
                 text,
                 ha='center', va='bottom')


ind = np.arange(N)  # the x locations for the groups
width = 0.35  # the width of the bars: can also be len(x) sequence
p1 = plt.bar(ind, small_k_means, 0.5 * width, color='#b0c0ad', hatch='//')
p2 = plt.bar(ind, medium_k_means, 0.5 * width, color='#20c78a', hatch='*')
p3 = plt.bar(ind, large_k_means, 0.5 * width, color='#20c7ff', hatch='O')
p1 = plt.bar([1], small_k_means[1], 0.5 * width, color='#b0c0ad', hatch='//')

p2 = plt.bar([1,2], [medium_k_means[1],medium_k_means[2]], 0.5 * width, color='#20c78a', hatch='*')


p6 = plt.bar(ind + 0.5 * width, small_slc, 0.5 * width, color='#b0c0ad', hatch='//')
p5 = plt.bar(ind + 0.5 * width, medium_slc, 0.5 * width, color='#20c78a', hatch='*')
p4 = plt.bar(ind + 0.5 * width, large_slc, 0.5 * width, color='#20c7ff', hatch='O')
p6 = plt.bar([2+0.5*width], small_slc[2], 0.5 * width, color='#b0c0ad', hatch='//')
p5 = plt.bar([2+0.5*width] , medium_slc[2], 0.5 * width, color='#20c78a', hatch='*')

p9 = plt.bar(ind + width, small_bmm, 0.5 * width, color='#b0c0ad', hatch='//')
p7 = plt.bar(ind + width, large_bmm, 0.5 * width, color='#20c7ff', hatch='O')
p8 = plt.bar(ind + width, medium_bmm, 0.5 * width, color='#20c78a', hatch='*')
threshold_l1 =1.0
thr = plt.plot([0, 2+2*width], [threshold_l1, threshold_l1], "k--")
# p8 = plt.bar([0 + width, 1+width], [medium_bmm[0],medium_bmm[1]], 0.5 * width, color='#20c78a', hatch='*')

autolabel(p3, "C1")
autolabel(p4, "C2", extra_height=2)
autolabel(p7, "C3")
# autolabel(p2)
# p4 = plt.bar(ind+0.5*width, womenMeans1, 0.5*width,
#              bottom=menMeans1, yerr=womenStd)
# p5 = plt.bar(ind+0.5*width, menMeans, 0.5*width,
#              bottom=menMeans1, yerr=womenStd)
# p6 = plt.bar(ind+0.5*width, menMeans2, 0.5*width, color='#ff0000',
#              bottom=menMeans1, yerr=womenStd)


plt.ylabel('Adjusted Rand Index')
plt.title('Similarity between predicted and true clusters')
plt.xticks(ind + 0.75 * width, ('Data sets with 1% error', 'Data sets with 5% error', 'Data sets with 10% error'))
plt.yticks(np.arange(0, 1.5, 0.2))

plt.ylim(ymax=1.5,ymin=0)

plt.xlim(xmax=3 - 1.2 * width)

obj_0 = AnyObject("C1", "black")
obj_1 = AnyObject("C2", "black")
obj_2 = AnyObject("C3", "black")
thr = AnyObject("- - -", "black")

plt.legend((p1[0], p2[0], p3[0], thr, obj_0, obj_1, obj_2),
           ('Small data set', 'Medium data set', 'Large data set',"Max threshold" ,"k-means", "SLC", "BMM"),
           handler_map={obj_0: AnyObjectHandler(), obj_1: AnyObjectHandler(), obj_2: AnyObjectHandler(), thr: AnyObjectHandler()})

plt.show()
