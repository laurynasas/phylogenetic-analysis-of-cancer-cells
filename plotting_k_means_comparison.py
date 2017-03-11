import matplotlib.pyplot as plt
import numpy as np

def autolabel(rects):
    """
    Attach a text label above each bar displaying its height
    """
    for rect in rects:
        height = rect.get_height()

        ax.text(rect.get_x() + rect.get_width()/2., 1.01 * height,
                '%s' % str(round(height,3)),
                ha='center', va='bottom')

x = (
    0.5979963153,
    0.2307576041,
    -0.06374255207,
)
y = (
1,
0.9489112699,
0.8023501378
)
N = 3
ind = np.arange(N)
nd = np.arange(N)  # the x locations for the groups
width = 0.3  # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(ind, x, width, color='r')

rects2 = ax.bar(ind + width, y, width, color='y')


# add some text for labels, title and axes ticks
ax.set_ylabel('Predicted Silhouette Score')
ax.set_title('Comparison of my k-means clustering implementation vs sklearn library')
ax.set_xticks(ind + width)
ax.set_xticklabels(('Small dataset (1% error)', 'Medium dataset (1% error)', 'Large dataset (1% error)'))
ax.set_ylim([0,1.3])
threshold_s1 = 1.0
# plt.text(1.1,2-2*width,'blah',rotation=90)

thr, = ax.plot([0., 2*width], [threshold_s1, threshold_s1], "k--")
threshold_m1 =0.925534404
ax.plot([1, 1+2*width], [threshold_m1, threshold_m1], "k--")
threshold_l1 =0.8643318013
ax.plot([2, 2+2*width], [threshold_l1, threshold_l1], "k--")

ax.legend((rects1[0], rects2[0], thr), ('k-medoids by me', 'k-means by sklearn', "Max score from true data"))

autolabel(rects1)
autolabel(rects2)

plt.show()