import numpy as np
import matplotlib as mpl
import sys


#matplotlib.use['qt4agg']
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster.hierarchy import cophenet
from scipy.spatial.distance import pdist
from mpl_toolkits import axes_grid1
import matplotlib.ticker as ticker

from scipy import stats
from itertools import *

import seaborn as sns

## default style setting, some of them were replaced in matplotlib 2.0
plt.rcParams['patch.force_edgecolor'] = False


#legend
plt.rcParams['legend.fancybox']       = True
plt.rcParams['legend.numpoints']      = 1
plt.rcParams['legend.framealpha']     = 0.7
plt.rcParams['legend.scatterpoints']  = 1
#plt.rcParams['borderaxespad']         = 0.5 # teh pad between the axes and the legend border

plt.rcParams['axes.xmargin']          = 0

plt.rcParams['xtick.top']             = False
plt.rcParams['ytick.right']           = False
plt.style.use(u'seaborn-bright')
#plt.style.use(u'grayscale')

myfavcolors=['grey','tomato','peru','dodgerblue','brown','darkslategray','lightsalmon','limegreen','deeppink','mediumpurple','darkkhaki','aqua','coral','gold','cornflowerblue','crimson']
mat2color=['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b','#e377c2','#7f7f7f','#bcbd22','#17becf']

'''
============================== ===============================================
marker                         description
============================== ===============================================
`"."`                          point
`","`                          pixel
`"o"`                          circle
`"v"`                          triangle_down
`"^"`                          triangle_up
`"<"`                          triangle_left
`">"`                          triangle_right
`"1"`                          tri_down
`"2"`                          tri_up
`"3"`                          tri_left
`"4"`                          tri_right
`"8"`                          octagon
`"s"`                          square
`"p"`                          pentagon
`"P"`                          plus (filled)
`"*"`                          star
`"h"`                          hexagon1
`"H"`                          hexagon2
`"+"`                          plus
`"x"`                          x
'''

myfavmarkers = ['o','+','2','x','5','.','<','8','s']
myfavcmap = sns.diverging_palette(220, 10, as_cmap=True)

'''
==============================================================================
==============================================================================
# text style
text_style = dict(horizontalalignment='right', verticalalignment='center',fontsize=12, fontdict={'family': 'monospace'})

# text position
from matplotlib.transforms import blended_transform_factory
# the reference point (0 in Axes coords, y tick value in Data coords).
reference_transform = blended_transform_factory(ax.transAxes, ax.transData)
for i, (name, linestyle) in enumerate(linestyles.items()):
    ax.annotate(str(linestyle), xy=(0.0, i), xycoords=reference_transform,
                xytext=(-6, -12), textcoords='offset points', color="blue",
                fontsize=8, ha="right", family="monospace")

'''

def styleNum(number):
    colors = cycle(myfavcolors)
    markers = cycle(myfavmarkers)
    retcol = []
    retmak = []
    for i in xrange(number):
        c=colors.next()
        m=markers.next()
        retcol.append(c)
        retmak.append(m)
    return retcol,retmak


def autolabel(rects,horizontal=0):
    # attach some text labels
    ax=plt.gca()
    for rect in rects:
        if not horizontal:
            height = rect.get_height()
            ax.text(rect.get_x() + rect.get_width()/2., 1.02*height,'%d' % int(height),ha='center', va='bottom')
        else:
            height = rect.get_height()
            ax.text(rect.get_width()*1.02,rect.get_y()+height/2., '%d' % int(height),ha='center', va='center')
    return 0

def setupTick(ax):
    ax.spines['right'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.yaxis.set_major_locator(ticker.NullLocator())   # remove all tickers,it seems that it shows only major ticks by default,so this line is enough normaly
    ax.yaxis.set_minor_locator(ticker.NullLocator())
    
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5))  # show a major ticker every 0.5 space
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))

    majors = [0, 1, 5]
    ax.xaxis.set_major_locator(ticker.FixedLocator(majors))  # also you can define any ticker position, even irregularly

    ax.xaxis.set_major_locator(ticker.MaxNLocator(4))  # MaxN Locator

    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.tick_params(which='major', width=1.00)
    ax.tick_params(which='major', length=5)
    ax.tick_params(which='minor', width=0.75)
    ax.tick_params(which='minor', length=2.5)
    ax.set_xlim(0, 5)
    ax.set_ylim(0, 1)
    ax.patch.set_alpha(0.0)


def barplot(plotarr,arrlabel,xlabel="x",ylabel="y",title="barplot",picname="test",color='gray',log=0,horizontal=0,rotation=0,):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    y_pos = np.arange(len(plotarr))
    if not horizontal:
        if log:ax.set_yscale('log')
        rects=ax.bar(y_pos,plotarr,align='center',color=color,log=log)
        ax.set_xticks(y_pos)
        ax.set_xticklabels(arrlabel,rotation=rotation)
        autolabel(rects,horizontal)
        ax.set_ylim(top=ax.axis()[-1]*1.1)
    else:
        if log:ax.set_xscale('log')
        rects=ax.barh(y_pos, plotarr, align='center',color=color,log=log,)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(arrlabel,rotation=rotation)
        autolabel(rects,horizontal)
        ax.set_xlim(right=ax.axis()[1]*1.1)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    plt.savefig("%s.png"%picname,format='png',dpi=300)
    plt.savefig("%s.svg"%picname,format='svg',dpi=300)
    plt.clf()
    plt.close()
    return 0

def bargroup(matrix,samplenames,ticklabels,xlabel="x",ylabel="y",title="bargroupplot",picname="test",colors=['gray','red','yellow'],width = 0.35,log=0,horizontal=0):
    mat = np.asarray(matrix)
    r,c = mat.shape
    print r,c
    try:
        assert c == len(samplenames)
    except AssertionError,e:
        sys.stderr.write("wrong data dimension!!!\n")
    if len(colors) != c:
        colors,markers = styleNum(c)
    fig, ax = plt.subplots()
    ind = np.arange(r)
    for i in range(c):
        rects=ax.bar(ind+(i+1)*width, mat[:,i], width, color=colors[i],label=samplenames[i])
        autolabel(rects,horizontal)
    ax.set_xticks(ind + width*(r+1)/2)
    ax.set_xticklabels(ticklabels,rotation=rotation)
    ax.set_ylim(top=ax.axis()[-1]*1.1)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    plt.savefig("%s.png"%picname,format='png',dpi=300)
    plt.savefig("%s.svg"%picname,format='svg',dpi=300)
    plt.clf()
    plt.close()
    return 0

def clearaxis(ax):
    ''' remove xaxis and yaxis of the ax '''
    ax.margins(0.2)
    ax.set_axis_off()

def clearHalfaxis(ax):
    '''Only show ticks on the left and bottom spines. '''
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')


def boxplot(matrix,labels,xlabel="x",ylabel='y',title='boxplot',picname="test",rotation=0,patch_artist=True,notch=True):
    mat = np.asarray(matrix)
    r,c = mat.shape
    print r,c
    try:
        assert c == len(labels)
    except AssertionError,e:
        sys.stderr.write("wrong data dimension!!!\n")
    fig, ax = plt.subplots()
    #mat=np.asarray(matrix).T
    boxes=ax.boxplot(mat,notch=notch,patch_artist=patch_artist)
    colors,markers = styleNum(len(labels))
    for patch, color in zip(boxes['boxes'], colors):
        patch.set_facecolor(color)
    ax.set_xticks([i+1 for i in range(len(labels))],)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xticklabels(labels,rotation=rotation)
    ax.set_title(title)
    plt.tight_layout()
    plt.savefig("%s.png"%picname,format='png',dpi=300)
    plt.savefig("%s.svg"%picname,format='svg',dpi=300)
    plt.clf()
    plt.close()
    return 0

def densityplot(plotarrs,labels,xlabel="x",ylabel='y',title='density',picname="test"):
    try:
        assert len(plotarrs) == len(labels)
    except AssertionError,e:
        sys.stderr.write("inconsitent data and labels!!!\n")
    fig, ax = plt.subplots()
    for i,label in enumerate(labels):
        data=plotarrs[i]
        data=np.asarray(data)
        density = stats.kde.gaussian_kde(data)
        x=np.linspace(data.min(),data.max(), num=100)
        ax.plot(x, density(x),color=myfavcolors[i],label=label)
    if len(labels) <= 6:
        plt.legend(loc='best')
    else:
        ncol=min(int(len(labels)/2.0),5)
        plt.legend(loc='upper center',bbox_to_anchor=(.5,1.1),ncol=ncol,fancybox=True,frameon=True)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plt.savefig("%s.png"%picname,format='png',dpi=300)
    plt.savefig("%s.svg"%picname,format='svg',dpi=300)
    plt.clf()
    plt.close()
    return 0

def clusterHC(data,method='average', metric='euclidean',picname='test'):
    fig, ax = plt.subplots()
    #data=data.T
    Z=linkage(data,method=method,metric=metric)
    c, coph_dists = cophenet(Z, pdist(data))
    a=dendrogram(Z,leaf_rotation=0, leaf_font_size=8.,orientation='left') 
    #print Z[1]
    #plt.savefig("%s.png"%picname,format='png',dpi=300)
    #plt.savefig("%s.svg"%picname,format='svg',dpi=300)
    plt.clf()
    plt.close()
    return Z

def clusterheatmap(mat,samplenames,featurenames,labelsize=14,plotxlabel=1,plotylabel=1,cbarlabel="Expression",cut_tree = 3,cmap='coolwarm',picname='cluster_heatmap'):
    # mat is r * c matrix
    r,c=mat.shape
    figw=12 if c<30 else c/3
    figh=10 if r<40 else r/4.0 
    
    fig, ax = plt.subplots(figsize=(figw,figh))
    divider = axes_grid1.make_axes_locatable(ax)
    sys.stderr.write('plot size is %d*%d\n'%(figw,figh))
    try:
        assert len(samplenames) == c and len(featurenames) == r
    except AssertionError,e:
        sys.stderr.write("inconsitent data and labels!!!\n")
    if len(featurenames)>50:
        plotylabel=0
    cmap='RdYlBu_r'
    vmin = np.floor(np.min(mat))
    vmax = np.ceil(np.max(mat))
    vmax = max([vmax,abs(vmin)]) # choose larger of vmin and vmax
    my_norm = mpl.colors.Normalize(vmin, vmax)
    
    leftx=0.06
    rightx=0.85
    topy=0.9
    bottomy=0.1
    rotationx=90
    if len(featurenames[0]) > 12:
        #labelsize=8
        rightx=0.8
        leftx=0.05
    if len(samplenames[0])>10:
        bottomy=0.15
        rotationx=45

    fontdict={'fontsize':labelsize}

    leftax=divider.append_axes(position="left",size=1.3,pad=0,frameon=False)
    clearaxis(leftax)
    distmat_left=clusterHC(mat)
    featuredict=dendrogram(distmat_left,ax=leftax,leaf_rotation=0, leaf_font_size=8.,orientation='left',no_labels=True)
    orderByfeature = mat[featuredict['leaves'],:]
    plotfeaturename = []
    for i in featuredict['leaves']:
        plotfeaturename.append(featurenames[i])
    
    topax = divider.append_axes(position="top",size=0.5,pad=0,frameon=False)
    #topax.set_position([0.125, 0.50,1.2, 0.3])
    clearaxis(topax)
    distmat_top=clusterHC(orderByfeature.T)
    featuredict=dendrogram(distmat_top,ax=topax,leaf_rotation=0, leaf_font_size=8.,orientation='top',no_labels=True)
    plotdata = orderByfeature[:,featuredict['leaves']]
    plotsamplename = []
    for i in featuredict['leaves']:
        plotsamplename.append(samplenames[i])

    im=ax.imshow(plotdata, aspect='auto', cmap=cmap,interpolation='nearest',norm=my_norm,)
    ax.set_xticks([i for i in range(c)],)
    ax.set_xticklabels(plotsamplename,rotation=rotationx,fontdict=fontdict)
    ax.set_yticks([i for i in range(r)],)
    ax.tick_params(axis=u'both', which=u'both',length=0)
    ax.yaxis.set_ticks_position('right')
    if plotylabel:
        ax.set_yticklabels(plotfeaturename,fontdict=fontdict)
    
    cax=fig.add_axes([0.04,0.70,0.01,0.25],frameon=False,)
    cax.set_title(cbarlabel,)
    plt.colorbar(im,cax=cax,)
    cax.yaxis.set_ticks_position('left')
    cax.yaxis.set_label_position('left')
    plt.subplots_adjust(left  = leftx,right = rightx,bottom=bottomy,top=topy)
    #if len(featurenames[0])>15:
    #    plt.tight_layout()
    plt.savefig("%s.png"%picname,format='png',dpi=300)
    plt.savefig("%s.svg"%picname,format='svg',dpi=300)
    plt.clf()
    plt.close()
    return 0


def corrheatplot(mat,samplenames,labelsize=10,picname='corrheat'):
    r,c=mat.shape
    try:
        assert r==len(samplenames)
    except AssertionError,e:
        mat=mat.T
        sys.stderr.write("inconsitent data dimension with length of labels!!! Will try to transpose the matrix\n")
    corrMat = np.corrcoef(mat)
    r,c = corrMat.shape
    fig, ax = plt.subplots(figsize=(10,10))
    my_norm = mpl.colors.Normalize(-1, 1)
    cmap=sns.diverging_palette(220, 10, as_cmap=True)
    im=ax.imshow(corrMat,aspect='auto', cmap=cmap,interpolation='nearest',norm=my_norm,)
    ax.grid(False)
    fontdict={'fontsize':labelsize}
    rotationx=90
    ax.set_xticks([i for i in range(r)],)
    ax.set_xticklabels(samplenames,rotation=rotationx,fontdict=fontdict)
    ax.set_yticks([i for i in range(r)],)
    rotationx=0
    ax.set_yticklabels(samplenames,rotation=rotationx,fontdict=fontdict)
    ax.yaxis.set_ticks_position('right')
    ax.tick_params(axis=u'both', which=u'both',length=0)

    cax=fig.add_axes([0.04,0.70,0.01,0.25],frameon=False,)
    plt.colorbar(im,cax=cax,)
    cax.yaxis.set_ticks_position('right')
    cax.yaxis.set_label_position('right')
    plt.savefig("%s.png"%picname,format='png',dpi=300)
    plt.savefig("%s.svg"%picname,format='svg',dpi=300)
    plt.clf()
    plt.close()
    return 0


def scatterplot(xarr,yarr,data,samplenames,colors=['C0','C1','C2'],marker=0,xlabel="x",ylabel='y',title='scatter plot',cm=0,picname="scatter",alpha=0.8,figsize=(8,6)):
    c = len(xarr)
    from sklearn import datasets
    circles = datasets.make_circles()
    try:
        assert c == len(samplenames)
    except AssertionError,e:
        sys.stderr.write("inconsitent data dimension with length of samplenames!!!\n")
    if len(colors) != c:
        colors,markers = styleNum(c)
    if not marker:
        markers = ['o',]*c
    fig, ax = plt.subplots(figsize=figsize)
    xlimmin = min([min(i) for i in xarr])
    xlimmax = max([max(i) for i in xarr])
    ylimmin = min([min(i) for i in yarr])
    ylimmax = max([max(i) for i in yarr])
    cmap='RdYlBu_r'
    #cmap=myfavcmap
    vmax = np.ceil(max([max(i) for i in data]))
    vmin = np.floor(min([min(i) for i in data]))
                                                 #data = np.asarray(data,dtype=np.float64) in case of a non-matrix data
    vmax = max([vmax,abs(vmin)]) # choose larger of vmin and vmax
    my_norm = mpl.colors.Normalize(vmin, vmax)
    for i in range(c):
        if cm:
            #ax.scatter(xarr[i], yarr[i],c=(0.2,0.4,0.6,0.8), s=data[i],marker=markers[i], alpha=0.8,label=samplenames[i], norm=my_norm)
            im=ax.scatter(xarr[i], yarr[i],c=data[i], s=data[i],marker=markers[i], alpha=0.8,label=samplenames[i], cmap=cmap)
            cax=fig.add_axes([0.04,0.70,0.01,0.25],frameon=False,)
            plt.colorbar(im,cax=cax,)
            '''
            cax_yticks = np.linspace(vmin,vmax,6)
            cax.set_yticks(cax_yticks)
            cax.set_yticklabels(cax_yticks)
            '''
            tickL = cax.yaxis.get_ticklabels()
            for t in tickL:
                t.set_fontsize(t.get_fontsize() - 0.8)
            cax.yaxis.set_ticks_position('right')
            cax.yaxis.set_label_position('right')
        else:
            ax.scatter(xarr[i], yarr[i], s=data[i], c=colors[i], marker=markers[i], alpha=alpha, label=samplenames[i])
            ax.legend(loc='best')
    ax.set_title(title)
    ax.set_xlim(xlimmin*0.85,xlimmax*1.1)
    ax.set_ylim(ylimmin*0.85,ylimmax*1.1)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    #ax.legend(loc='best')
    plt.savefig("%s.png"%picname,format='png',dpi=300)
    plt.savefig("%s.svg"%picname,format='svg',dpi=300)
    plt.clf()
    plt.close()
    return 0 


def plot_chromosome():
    pass
