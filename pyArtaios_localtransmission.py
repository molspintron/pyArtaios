import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
import numpy as np
from matplotlib.widgets import Slider, Button, RadioButtons


class Arrow3D(FancyArrowPatch):
    '''class for making nice arrows'''
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)


def load_t(file_localtransmission):
    '''load matrix with local transmission coefficients'''

    with open(file_localtransmission) as file:
        lines = file.readlines()

    # get number of atoms from length of second line
    n_atoms = len(lines[1].split())
    # get number of energies by looking at the length of the
    # file
    n_energies = int(len(lines)/(n_atoms + 1))

    # make iterator from list of lines
    lines = iter(lines)

    local_transmissions_alpha = {}
    local_transmissions_beta  = {}

    for i in range(n_energies):
        temp = next(lines).strip('\n').split()

        energy = round(float(temp[0]), 2)
        spin = int(temp[1])
        temp_mat = []
        for j in range(n_atoms):
            temp_mat.append(next(lines).strip('\n').split())

        if spin == 1:
            local_transmissions_alpha[energy] = np.array(temp_mat, dtype = 'float')
        else:
            local_transmissions_beta[energy] = np.array(temp_mat, dtype = 'float')

    return local_transmissions_alpha, local_transmissions_beta


def plot_local_interactive(atoms, coords, localt, tlabels = True, threshold = 0.1, threshold_method = 'relative'):
    '''interactive plot for the local transmissions
    realised with matplotlib

    Parameters:
    -----------
    atoms : list of atom symbols
    coords : numpy array with cartesian coordinates
    localt : dictionary with local transmissions between
        all atoms for all energies
    tlabels : switch for displaying the local transmission
    threshold : initial threshold value
    threshold_method : method to determine the threshold,
        either in absolute values or relative to the maximal
        local transmission of the current energy
    '''

    # some settings which control the
    # displayed molecule (size and color)
    # using molecules with atoms not listed here
    # will lead to a crash, then just add it
    size_scale = 50
    size = {'au': 10*size_scale, 'c': 5*size_scale, 's': 5*size_scale, 'si': 7*size_scale, 'n': 5*size_scale, 'h': 1*size_scale, 'fe':7*size_scale, 'o': 5*size_scale}
    color = {'au': '#ff9900', 'c': '#444444', 's': '#ffff00', 'si': '#cccc00','n': '#0000ff', 'h':  '#aaaaaa', 'fe': '#962700', 'o': '#aa0000', 'f': '#000099'}

    # initialize figure
    fig = plt.figure(1, figsize = (8, 4))
    fig.clf()
    ax = Axes3D(fig, proj_type='persp')
    ax.grid(False)
    plt.axis('off')

    # select alpha local transmissions for the start
    # localt = local_trans[0]

    # create lists without gold for local transmissions
    central_atoms = []
    central_coords = []

    for i in range(len(atoms)):
        if atoms[i].lower() != 'au':
            central_atoms.append(atoms[i])
            central_coords.append(coords[i, :])

    central_atoms = np.array(central_atoms)
    central_coords = np.array(central_coords)

    # plot atoms
    for i in range(len(atoms)):
        ax.scatter(coords[i, 2], coords[i, 1], coords[i, 0],
                   s = size[atoms[i].lower()],
                   c = color[atoms[i].lower()])

    # get initial energy
    energy = list(localt[0].keys())[0]

    # set initial spin
    spin = 0

    # build alpha list of arrows
    arrows_list = []
    for i in range(len(central_atoms)):
        for j in range(len(central_atoms)):
            # get local transmission
            loct = localt[spin][energy][i, j]

            # threshold for plotting
            if threshold_method == 'relative':
                thres = localt[spin][energy].max()*threshold
            else:
                thres = threshold

            # we want the arrows to start and stop
            # a bit after/before the atom
            distance = 0.1

            vec = [central_coords[j, 2] - central_coords[i, 2],
                   central_coords[j, 1] - central_coords[i, 1],
                   central_coords[j, 0] - central_coords[i, 0]]

            l = 0.1/(np.sqrt(vec[2]**2 + vec[1]**2 + vec[0]**2))
            p1 = [central_coords[i, 2] + l*vec[0],
                  central_coords[i, 1] + l*vec[1],
                  central_coords[i, 0] + l*vec[2],
                 ]

            p2 = [central_coords[j, 2] - l*vec[0],
                  central_coords[j, 1] - l*vec[1],
                  central_coords[j, 0] - l*vec[2],
                 ]

            # now we plot the arrow with the new points
            arrows_list.append(Arrow3D([p1[0], p2[0]],
                               [p1[1], p2[1]],
                               [p1[2], p2[2]],
                               mutation_scale=10, lw=(4*loct/localt[spin][energy].max()),
                               arrowstyle="-|>"))


            ax.add_artist(arrows_list[-1])
            ax.artists[-1].set_visible(False)

            # set visible if above initial threshold
            if loct > thres:
                ax.artists[-1].set_visible(True)

    # update function for slider
    def update(val):
        # get energy and threshold from slider
        energy = round(s_energy.val,2)
        threshold = s_thresh.val+0.00001

        counter = 0
        for i in range(len(central_atoms)):
            for j in range(len(central_atoms)):
                # get local transmission
                loct = localt[spin][energy][i, j]
                ax.artists[counter].update({'lw' : 4*loct/localt[spin][energy].max()})

                # threshold for plotting
                if threshold_method == 'relative':
                    thres = localt[spin][energy].max()*threshold
                else:
                    thres = threshold

                # make arrows visible
                if loct > thres:
                    ax.artists[counter].set_visible(True)
                else:
                    ax.artists[counter].set_visible(False)

                counter += 1

        fig.canvas.draw_idle()

    def update_spin(label):
        nonlocal spin

        if label == 'beta':
            spin = 1
        else:
            spin = 0

        update(0)

    # add slider for energy
    ax_slider_energy = plt.axes([0.25, 0.17, 0.65, 0.03])
    ax_slider_thresh = plt.axes([0.25, 0.10, 0.65, 0.03])
    s_energy = Slider(ax_slider_energy, 'Energy',
                  min(list(localt[0].keys())),
                  max(list(localt[0].keys())),
                  valinit = 0,
                  valstep = abs(list(localt[0].keys())[1] - list(localt[0].keys())[0]))

    s_thresh = Slider(ax_slider_thresh, 'Threshold',
                  0, 1,valinit = 0.1)

    # add radio button for spin
    rax = plt.axes([0.025, 0.10, 0.10, 0.1])
    if len(localt[1]) != 0:
        spin_button = RadioButtons(rax, ('alpha', 'beta'), active=0)

    # "callback"
    s_energy.on_changed(update)
    s_thresh.on_changed(update)
    spin_button.on_clicked(update_spin)

    # Create cubic bounding box to simulate equal aspect ratio
    max_range = np.array([coords[:, 2].max()-coords[:, 2].min(), coords[:, 1].max()-coords[:, 1].min(), coords[:, 0].max()-coords[:, 0].min()]).max()
    Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(coords[:, 2].max()+coords[:, 2].min())
    Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(coords[:, 1].max()+coords[:, 1].min())
    Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(coords[:, 0].max()+coords[:, 0].min())

    # Comment or uncomment following both lines to test the fake bounding box:
    for xb, yb, zb in zip(Xb, Yb, Zb):
        ax.plot([xb], [yb], [zb], 'w')

    try:
        ax.set_aspect('equal')
    except:
        pass
    fig.tight_layout()
    plt.show()


def plot_local(atoms, coords, localt, tlabels = True, threshold = 0.1, threshold_method = 'relative'):
    '''plots molecule with local transmissions'''

    # some settings which control the
    # displayed molecule (size and color)
    size_scale = 50
    size = {'au': 10*size_scale, 'c': 5*size_scale, 's': 5*size_scale, 'si': 7*size_scale, 'n': 5*size_scale, 'h': 1*size_scale, 'fe':7*size_scale, 'o': 5*size_scale}
    color = {'au': '#ff9900', 'c': '#444444', 's': '#ffff00', 'si': '#cccc00','n': '#0000ff', 'h':  '#aaaaaa', 'fe': '#962700', 'o': '#aa0000', 'f': '#000099'}

    # initialize figure
    fig = plt.figure(1, figsize = (8, 4))
    fig.clf()
    ax = Axes3D(fig, proj_type='persp')
    ax.grid(False)
    plt.axis('off')

    # create lists without gold for local transmissions
    central_atoms = []
    central_coords = []

    for i in range(len(atoms)):
        if atoms[i].lower() != 'au':
            central_atoms.append(atoms[i])
            central_coords.append(coords[i, :])

    central_atoms = np.array(central_atoms)
    central_coords = np.array(central_coords)

    # plot arrows for local transmissions
    for i in range(len(central_atoms)):
        for j in range(len(central_atoms)):
            # get local transmission
            loct = localt[i, j]

            # threshold for plotting
            if threshold_method == 'relative':
                thres = localt.max()*threshold
            else:
                thres = threshold

            if loct > thres:

                # we want the arrows to start and stop
                # a bit after/before the atom
                distance = 0.1

                vec = [central_coords[j, 2] - central_coords[i, 2],
                       central_coords[j, 1] - central_coords[i, 1],
                       central_coords[j, 0] - central_coords[i, 0]]

                l = 0.1/(np.sqrt(vec[2]**2 + vec[1]**2 + vec[0]**2))
                p1 = [central_coords[i, 2] + l*vec[0],
                      central_coords[i, 1] + l*vec[1],
                      central_coords[i, 0] + l*vec[2],
                     ]

                p2 = [central_coords[j, 2] - l*vec[0],
                      central_coords[j, 1] - l*vec[1],
                      central_coords[j, 0] - l*vec[2],
                     ]

                # now we plot the arrow with the new points
                a = Arrow3D([p1[0], p2[0]],
                            [p1[1], p2[1]],
                            [p1[2], p2[2]],
                             mutation_scale=10, lw=(4*loct/localt.max()),
                             arrowstyle="-|>", color="r", alpha = loct/localt.max())

                ax.add_artist(a)

                if tlabels:
                    # find middle of the local transmission
                    l = 0.5*np.sqrt(vec[2]**2 + vec[1]**2 + vec[0]**2)/\
                        (np.sqrt(vec[2]**2 + vec[1]**2 + vec[0]**2))
                    p3 = [central_coords[i, 2] + l*vec[0],
                          central_coords[i, 1] + l*vec[1],
                          central_coords[i, 0] + l*vec[2],
                          ]
                    #add text for the rounded local transmission
                    ax.text(p3[0], p3[1], p3[2], str(round(loct,4)))

    # plot atoms
    for i in range(len(atoms)):
        ax.scatter(coords[i, 2], coords[i, 1], coords[i, 0],
                   s = size[atoms[i].lower()],
                   c = color[atoms[i].lower()])

    # Create cubic bounding box to simulate equal aspect ratio
    max_range = np.array([coords[:, 2].max()-coords[:, 2].min(), coords[:, 1].max()-coords[:, 1].min(), coords[:, 0].max()-coords[:, 0].min()]).max()
    Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(coords[:, 2].max()+coords[:, 2].min())
    Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(coords[:, 1].max()+coords[:, 1].min())
    Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(coords[:, 0].max()+coords[:, 0].min())

    # Comment or uncomment following both lines to test the fake bounding box:
    for xb, yb, zb in zip(Xb, Yb, Zb):
        ax.plot([xb], [yb], [zb], 'w')

    try:
        ax.set_aspect('equal')
    except:
        pass
    fig.tight_layout()
    plt.show()
