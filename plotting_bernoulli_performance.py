import matplotlib
import matplotlib as mpl

matplotlib.use('QT4Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

mpl.rcParams['legend.fontsize'] = 10
fig = plt.figure()
ax = fig.gca(projection='3d')

# Make data.
X = (
    1,
    2,
    3,
    4,
    5,
    6,
    7,
    8,
    9,
    10,
    11,
    12,
    13,
    14,
    15,
    16,
    17,
    18,
    19,
    20,
    21,
    22,
    23,
    24,
    25
)
Z = (
    0.2451670915,
    0.5970772489,
    0.6221022494,
    0.6839926589,
    0.743230738,
    0.7844683574,
    0.7848660102,
    0.8143552828,
    0.8141547922,
    0.8233092602,
    0.824309257,
    0.8300712744,
    0.8326346053,
    0.8334916117,
    0.826418805,
    0.8316589039,
    0.8423091746,
    0.8395929651,
    0.8446476358,
    0.8389935835,
    0.8451420274,
    0.8518168249,
    0.8463992285,
    0.8397485283,
    0.8485602733
)

# R = np.sqrt(X ** 2 + Y ** 2)
Y = (
    0.2117615795,
    0.4258240795,
    0.6691293192,
    0.8769654989,
    1.049256108,
    1.24985158,
    1.45722224,
    1.682014341,
    1.872970479,
    2.176216419,
    2.308609722,
    2.837633572,
    3.11453239,
    3.229220581,
    3.580776422,
    3.37629714,
    3.512906389,
    3.713023469,
    3.90750319,
    4.149899018,
    4.33176506,
    4.531050241,
    4.711484763,
    5.05448488,
    5.24570801
)
Z = Z[::-1]
Y = Y[::-1]
X = X[::-1]
# X, Y = np.meshgrid(X, Y)

ax.plot(X, Y, Z, label='parametric curve')
for x, y, z in zip(X, Y, Z):
    ax.plot([x, x], [y, y], [0, z], "k--", label='parametric curve')

# ax.plot([len(X)*0], Y, Z, label='parametric curve')

ax.plot(X, Y)
# ax.plot(X, np.array(Y), Z, color='g')
ax.set_xlabel('Iterations')
ax.set_ylabel('Average time to cluster (sec)')
ax.set_zlabel('Silhouette score')


# ax.legend()
# for y,z in zip(Y,Z):
# plt.plot([Z[0],Z[1]], [X[0],X[1]], color='k', linestyle='-', linewidth=2)


class MyAxes3D(axes3d.Axes3D):
    def __init__(self, baseObject, sides_to_draw):
        self.__class__ = type(baseObject.__class__.__name__,
                              (self.__class__, baseObject.__class__),
                              {})
        self.__dict__ = baseObject.__dict__
        self.sides_to_draw = list(sides_to_draw)
        self.mouse_init()

    def set_some_features_visibility(self, visible):
        for t in self.w_zaxis.get_ticklines() + self.w_zaxis.get_ticklabels():
            t.set_visible(visible)
        self.w_zaxis.line.set_visible(visible)
        self.w_zaxis.pane.set_visible(visible)
        self.w_zaxis.label.set_visible(visible)

    def draw(self, renderer):
        # set visibility of some features False
        self.set_some_features_visibility(False)
        # draw the axes
        super(MyAxes3D, self).draw(renderer)
        # set visibility of some features True.
        # This could be adapted to set your features to desired visibility,
        # e.g. storing the previous values and restoring the values
        self.set_some_features_visibility(True)

        zaxis = self.zaxis
        draw_grid_old = zaxis.axes._draw_grid
        # disable draw grid
        zaxis.axes._draw_grid = False

        tmp_planes = zaxis._PLANES

        if 'l' in self.sides_to_draw:
            # draw zaxis on the left side
            zaxis._PLANES = (tmp_planes[2], tmp_planes[3],
                             tmp_planes[0], tmp_planes[1],
                             tmp_planes[4], tmp_planes[5])
            zaxis.draw(renderer)
        if 'r' in self.sides_to_draw:
            # draw zaxis on the right side
            zaxis._PLANES = (tmp_planes[3], tmp_planes[2],
                             tmp_planes[1], tmp_planes[0],
                             tmp_planes[4], tmp_planes[5])
            zaxis.draw(renderer)

        zaxis._PLANES = tmp_planes

        # disable draw grid
        zaxis.axes._draw_grid = draw_grid_old


plt.title('Bernoulli mixture model performance on medium dataset (1 % error)')

ax = fig.add_axes(MyAxes3D(ax, 'lr'))

plt.show()
