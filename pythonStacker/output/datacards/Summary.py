import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use a non-interactive backend like Agg
import matplotlib.pyplot as plt
import ROOT as root
#from root_numpy import array2hist, hist2array
#import mplhep as hep
from scipy.interpolate import CubicSpline
import CombineHarvester.CombineTools.plotting as plot

def FindCrossings(x,y,value):
    crossings = []
    for i_p, p in enumerate(x):
        try:
            if (y[i_p] < value and y[i_p + 1] > value) or (y[i_p] > value and y[i_p + 1] < value):

                #Dummy mean because many points in spline
                crossings.append(0.5*(x[i_p] + x[i_p+1]))
        except:
            continue
    return np.unique(crossings)

crossings_list = {
    'cQQ1': {
        'Asimov': {
            'frozen': {
                'min': 0.0,
                'one': np.array([-0.563,  0.700]),
                'two': np.array([-0.860,  1.000]),
            },
            'profiled': {
                'min': 0.0,
                'one': np.array([-0.611,  0.709]),
                'two': np.array([-0.917,  1.023]),
            },
        },
        'data': {
            'frozen': {
                'min': 0.711,
                'one': np.array([-0.785, -0.102,  0.166,  1.006]),
                'two': np.array([-1.082,  1.260]),
            },
            'profiled': {
                'min': 0.650,
                'one': np.array([-0.789,  0.977]),
                'two': np.array([-1.131,  1.261]),
            },
        },
    },

    'cQt1': {
        'Asimov': {
            'frozen': {
                'min': 0.0,
                'one': np.array([-1.134,  0.976]),
                'two': np.array([-1.641,  1.477]),
            },
            'profiled': {
                'min': 0.0,
                'one': np.array([-1.159,  1.075]),
                'two': np.array([-1.684,  1.592]),
            },
        },
        'data': {
            'frozen': {
                'min': 0.850,
                'one': np.array([-1.539,  1.389]),
                'two': np.array([-1.991,  1.834]),
            },
            'profiled': {
                'min': 0.150,
                'one': np.array([-1.197,  1.231]),
                'two': np.array([-1.912,  1.835]),
            },
        },
    },

    'cQt8': {
        'Asimov': {
            'frozen': {
                'min': 0.0,
                'one': np.array([-1.949,  2.411]),
                'two': np.array([-2.981,  3.448]),
            },
            'profiled': {
                'min': 0.0,
                'one': np.array([-2.122,  2.417]),
                'two': np.array([-3.204,  3.466]),
            },
        },
        'data': {
            'frozen': {
                'min': -1.679,
                'one': np.array([-2.790,  3.219]),
                'two': np.array([-3.706,  4.156]),
            },
            'profiled': {
                'min': -0.300,
                'one': np.array([-2.461,  2.503]),
                'two': np.array([-3.770,  3.962]),
            },
        },
    },

    'ctHIm': {
        'Asimov': {
            'frozen': {
                'min': 0.0,
                'one': np.array([-11.847, 12.214]),
                'two': np.array([-17.196, 17.543]),
            },
            'profiled': {
                'min': -12.500,
                'one': np.array([-20.929, 14.547]),
                'two': np.array([-25.952, 21.649]),
            },
        },
        'data': {
            'frozen': {
                'min': -6.800,
                'one': np.array([-14.091, 13.646]),
                'two': np.array([-19.128, 19.098]),
            },
            'profiled': {
                'min': -12.500,
                'one': np.array([-22.000, 15.500]),
                'two': np.array([-26.000, 21.500]),
            },
        },
    },

    'ctHRe': {
        'Asimov': {
            'frozen': {
                'min': 0.0,
                'one': np.array([-2.607,  3.229]),
                'two': np.array([-5.037,  7.867]),
            },
            'profiled': {
                'min': 2.500,
                'one': np.array([-2.134, 11.129]),
                'two': np.array([-5.171, 30.277]),
            },
        },
        'data': {
            'frozen': {
                'min': -0.188,
                'one': np.array([-2.962,  3.136]),
                'two': np.array([-5.590,  8.060]),
            },
            'profiled': {
                'min': 2.500,
                'one': np.array([-3.000, 14.500]),
                'two': np.array([-6.000, 31.000]),
            },
        },
    },

    'ctt': {
        'Asimov': {
            'frozen': {
                'min': 0.0,
                'one': np.array([-0.624,  0.739]),
                'two': np.array([-0.948,  1.067]),
            },
            'profiled': {
                'min': 0.0,
                'one': np.array([-0.679,  0.753]),
                'two': np.array([-1.012,  1.090]),
            },
        },
        'data': {
            'frozen': {
                'min': 0.660,
                'one': np.array([-0.872,  1.010]),
                'two': np.array([-1.172,  1.301]),
            },
            'profiled': {
                'min': 0.250,
                'one': np.array([-0.753,  0.888]),
                'two': np.array([-1.188,  1.262]),
            },
        },
    },
}





def lookup_crossings(crossings_dict, coeff, dataset, mode):
    """
    Retrieve crossings for a given coefficient (e.g. 'ctt'),
    dataset ('data' or 'Asimov'), mode ('profiled' or 'frozen'),
    and level ('one', 'two', 'min'). Returns a numpy array.
    """
    try:
        vals = crossings_dict[coeff][dataset][mode]
    except KeyError:
        return {}
    return vals


def tree_to_asymmgraph(tree, x_branch):#, y_branch, ex_low_branch, ex_high_branch, ey_low_branch, ey_high_branch):
    # Find branches by name, assuming a common naming convention
    branch_names = [b.GetName() for b in tree.GetListOfBranches()]
    #print (branch_names)
    
    graph = root.TGraphAsymmErrors()
    
    for i in range(tree.GetEntries()):
        tree.GetEntry(i)
        quantile_value = getattr(tree, "quantileExpected")
        # Only include points where quantileExpected > -1.5
        if quantile_value > -1.5:
            x_value = getattr(tree, x_branch)
            y_value = getattr(tree, "deltaNLL")

        graph.SetPoint(i, x_value, 2*y_value)
        graph.Sort()

    return graph



from ctypes import c_double

def getGraph(fname, w):
    # open file + build graph
    with root.TFile.Open(fname) as f:
        graph = tree_to_asymmgraph(f.Get("limit"), w)

    # pull out N points into NumPy arrays
    n = graph.GetN()
    x = np.empty(n)
    y = np.empty(n)
    for i in range(n):
        xd, yd = c_double(), c_double()
        graph.GetPoint(i, xd, yd)
        x[i], y[i] = xd.value, yd.value

    # compute crossings and minimum
    crossings = {
        "one": FindCrossings(x, y, 1.0),
        "two": FindCrossings(x, y, 4.0),
        "min": x[np.argmin(y)]
    }

    print(crossings)
    return crossings



def plotConstraints(path, WC, WCtex, Factors):

    fig, ax = plt.subplots(figsize=(8, 7))
    fig.set_size_inches(25,15)
    plt.subplots_adjust(left=0.17, right=0.98, top=0.95, bottom=0.1)
    ax.vlines(0,-2,2*7, color = "grey", linestyle = "dashed", linewidth = 2)
    ax.vlines(2,-2,2*7, color = "grey", linewidth = 2)

    for i_w, w in enumerate(WC):

        filepath = "v36_IM/Allops" + "/higgsCombine_"+w+"_combinedFit_profiled_data_fixed.root"#+w.split("_")[1]+"_combinedFit_profiled.root" 
        crossings = getGraph(filepath,w)
        filepath_frozen = path + "/higgsCombine_"+w+"_combinedFit_data.MultiDimFit.mH120.root"#+w.split("_")[1]+"_combinedFit.root" 
        crossings_frozen = getGraph(filepath_frozen,w)
        crossings = lookup_crossings(crossings_list, w, 'data', 'profiled') 
        crossings_frozen = lookup_crossings(crossings_list, w, 'data', 'frozen') 
        #crossings = lookup_crossings(crossings_list, w, 'Asimov', 'profiled') 
        #crossings_frozen = lookup_crossings(crossings_list, w, 'Asimov', 'frozen') 

        offset = 2.5
        #Floating 
        ax.text(2.2, offset*(len(WC) - i_w -1 - 0.1), "["  + str(round(crossings["one"][0],3)) + "," + str(round(crossings["one"][1],3))+ "]"  , fontsize=22, color='black')
        #2 sigma
        if len(crossings["two"]) >2:
            ax.text(2.8, offset*(len(WC) - i_w -1 - 0.1), "["  + str(round(crossings["two"][0],3)) + "," + str(round(crossings["two"][1],3))+ "]U[" + str(round(crossings["two"][2],3)) + "," + str(round(crossings["two"][3],3))+ "]" , fontsize=22, color='black')
        else:
            ax.text(3.2, offset*(len(WC) - i_w -1 - 0.1), "["  + str(round(crossings["two"][0],3)) + "," + str(round(crossings["two"][1],3))+ "]"  , fontsize=22, color='black')

        ax.text(2.2,14, r'1 $\sigma$ profiled', fontsize=22, color='black')
        ax.text(3.2,14, r'2 $\sigma$ profiled', fontsize=22, color='black')

        for k,v in crossings.items():
            crossings[k] = Factors[i_w]*v
            crossings_frozen[k] = Factors[i_w]*crossings_frozen[k]

        #Floating
        if i_w == 0:
            the_label = "Best profiled fit"
        else:
            the_label = None

        ax.plot([crossings["min"]],[offset*(len(WC) - i_w - 1)],marker = "o" , markersize = 10, markerfacecolor='red', markeredgecolor = "black", color = "black", linestyle = '', label = the_label)
        
        if i_w == 0:
            the_label = "Best frozen fit"
        else:
            the_label = None
        
        ax.plot([crossings_frozen["min"]],[offset*(len(WC) - i_w - 1) -0.5],marker = "o" , markersize = 10, markerfacecolor='maroon', markeredgecolor = "black", color = "black", linestyle = '', label = the_label)

        # Define a helper function for plotting
        def plot_crossings(crossings, frozen_crossings, label_prefix, offset, i_w, WC, color_profiled, color_frozen):
            """
            Plots the crossings for profiled and frozen data.
            All lines from 'one' are solid, and all lines from 'two' are dashed.
            Arguments:
            - crossings: Dictionary containing crossing data for profiled.
            - frozen_crossings: Dictionary containing crossing data for frozen.
            - label_prefix: Prefix for the labels, e.g., "q < 1".
            - offset: Offset multiplier for line positions.
            - i_w: Current index in WC loop.
            - WC: List of Wilson Coefficients.
            - color_profiled: Color for profiled lines.
            - color_frozen: Color for frozen lines.
            """
            for key, linestyle,linewidth in [("one", "solid",5), ("two", "dashed",3)]:
                # Profiled crossings
                for j in range(len(crossings[key]) // 2):  # Use xrange for Python 2
                    if i_w == 0 and j == 0:
                        the_label = "{} (profiled)".format(label_prefix)
                    else:
                        the_label = None
        
                    ax.hlines(
                        offset * (len(WC) - i_w - 1),
                        crossings[key][2 * j],
                        crossings[key][2 * j + 1],
                        color=color_profiled,
                        label=the_label,
                        linestyle=linestyle,
                        linewidth = linewidth
                    )
        
                # Frozen crossings
                for j in range(len(frozen_crossings[key]) // 2):  # Use xrange for Python 2
                    if i_w == 0 and j == 0:
                        the_label = "{} (frozen)".format(label_prefix)
                    else:
                        the_label = None
        
                    ax.hlines(
                        offset * (len(WC) - i_w - 1) - 0.5,
                        frozen_crossings[key][2 * j],
                        frozen_crossings[key][2 * j + 1],
                        color=color_frozen,
                        label=the_label,
                        linestyle=linestyle,
                        linewidth = linewidth 
                    )
        
        # Use the helper function for plotting
        plot_crossings(
            crossings,
            crossings_frozen,
            label_prefix="q < 1",
            offset=offset,
            i_w=i_w,
            WC=WC,
            color_profiled="forestgreen",
            color_frozen="black",
        )
        #WC name and factor
        ax.text(-2.90, offset*(len(WC) - i_w -1 - 0.1)-0.25, WCtex[i_w], fontsize=35, color='black')
        ax.text(-2.55, offset*(len(WC) - i_w -1 - 0.1)-0.25,"[x" + str(Factors[i_w]) + "]" , fontsize=30, color='black')

    handles, labels = ax.get_legend_handles_labels()
    lgd = ax.legend(handles, labels, loc='upper left', bbox_to_anchor=(-0.0,0.9), ncol = 3, fontsize = 30)
    lgd.get_frame().set_facecolor('none')
    lgd.get_frame().set_linewidth(0)

    ax.tick_params(left = False, labelleft = False) 
    ax.set_xlim(-2,4)
    ax.set_xlabel("Wilson coefficient value", fontsize = 35)
    
    ax.set_ylim(-2,20)
    
    the_ticks = np.arange(-2, 2.2, 0.2)
    ticklabels = ["" for i in the_ticks]

    ticklabels[0] = round(the_ticks[0],2)
    ticklabels[5] = round(the_ticks[5],2)
    ticklabels[10] = 0
    ticklabels[15] = round(the_ticks[15],2)
    ticklabels[20] = round(the_ticks[20],2)
    ax.set_xticks(the_ticks)

    ax.set_xticklabels(ticklabels)
    ax.tick_params(axis = 'both', direction = 'in', which = 'both', length = 10, width = 2)

    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(2)

    ax.tick_params('x', labelsize = 30)  

    fig.savefig("/user/mshoosht/public_html/Interpretations/Plots/SS2L_3L_fit_v36/SummaryScan_data.pdf")
    plt.clf()
    plt.close()

    #fig.tight_layout(pad = 5)




if __name__ == "__main__":

    print("Plotting Constraints")
    f = "v36_IM/Allops"
    WC = ["ctt","cQQ1","cQt1","cQt8","ctHRe","ctHIm"] 
    Factors = [1,1,1,0.4,0.05,0.05]
    WCtex = [r'$c_{tt}$',r'$c_{QQ}^{1}$',r'$c_{Qt}^{1}$',r'$c_{Qt}^{8}$',r'$c_{tH}^{\Re}$',r'$c_{tH}^{\Im}$'] 
    plotConstraints(f, WC, WCtex = WCtex, Factors = Factors)

