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


N = 3
hier_pre = (
    1,
    0.925534404,
    0.7717890624
)

k_pred = (
    1,
    0.9489112699,
    0.8023501378
)

bern_pre = (
    0.9763164525,
    0.8385854643,
    0.7744317069
)

# men_std = (2, 3, 4, 1, 2)

ind = np.arange(N)  # the x locations for the groups
width = 0.2  # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(ind, bern_pre, width, color='r')

rects2 = ax.bar(ind + width, hier_pre, width, color='y')

rects3 = ax.bar(ind + 2*width, k_pred, width, color='g')

# add some text for labels, title and axes ticks
ax.set_ylabel('Predicted Silhouette Score')
ax.set_title('Score by clustering algorithm')
ax.set_xticks(ind + 1.5*width)
ax.set_xticklabels(('Small_1', 'Medium_1', 'Large_1'))
ax.set_ylim([0,1.3])
threshold_s1 = 1.0
# plt.text(1.1,2-2*width,'blah',rotation=90)
ax.plot([0., 1-2*width], [threshold_s1, threshold_s1], "k--")
threshold_m1 =0.925534404
ax.plot([1, 2-2*width], [threshold_m1, threshold_m1], "k--")
threshold_l1 =0.8643318013
ax.plot([2, 3-2*width], [threshold_l1, threshold_l1], "k--")

ax.legend((rects1[0], rects2[0], rects3[0]), ('Hierarchical', 'k-means', 'Bernoulli mixture'))

autolabel(rects1)
autolabel(rects2)
autolabel(rects3)


plt.show()
