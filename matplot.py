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
import matplotlib.gridspec as gridspec

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

fig, axes = plt.subplots(nrows=2, ncols=2)
ax0, ax1, ax2, ax3 = axes.flatten()

font = {'family': 'serif', 'color':  'darkred', 'weight': 'normal', 'size': 16,}

arrowprop = dict(arrowprops=u'wedge,connectionstyle="arc3,rad=-0.05",pacthB=None,shrinkA=5,shrinkB=5,bbox=dict(boxstyle=bboxstylename, fc="w", ec="k")))
#patchA (or patchB) is given, the returned path is clipped so that it start (or end) from the boundary of the patch. The path is further shrunk by shrinkA (or shrinkB) which is given in points.

ax.axis([xmin, xmax, ymin, ymax])
fig.add_axes([*left*, *bottom*, *width*, *height*])


**********************************************************
Artist.set_clip_path(path,)
set_clip_path is a method of artist, you should supply a path instance or patch instance as parameter, but I don't know if we can control the direction of clipping, ie inner or outer
**********************************************************

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
            ax.text(rect.get_x() + rect.get_width()/2., rect.get_y() + height()/2.0*0.8, str(int(height)),ha='center', va='bottom')
        else:
            height = rect.get_width()
            ax.text(rect.get_x() + rect.get_width()/2., rect.get_y() + rect.get_height()/2.0*0.9, str(int(height)),ha='center', va='center')
    return 0

def setupTick(ax):
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter
    ax.spines['right'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_visible(False)
    ax.yaxis.set_major_locator(ticker.NullLocator())   # remove all tickers,it seems that it shows only major ticks by default,so this line is enough normaly
    ax.yaxis.set_minor_locator(ticker.NullLocator())
    
    ax.xaxis.set_major_locator(ticker.LinearLocator())  #evenly spaced ticks from min to max

    majorFormatter = FormatStrFormatter('%d')
    ax.xaxis.set_major_formatter(majorFormatter) # with tickers on but no tickerlabels
    ax.xaxis.set_minor_formatter(ticker.NullFormatter())

    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5))  # show a major ticker every 0.5 space
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
    ax.xaxis.set_major_formatter(ticker.FixedFormatter(majors))


    majors = [0, 1, 5]
    ax.xaxis.set_major_locator(ticker.FixedLocator(majors))  # also you can define any ticker position, even irregularly
    
    ax.xaxis.set_major_locator(ticker.MaxNLocator(4))  # MaxN Locator

    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.tick_params(which='major', width=1.00)
    ax.tick_params(which='major', length=5)
    ax.tick_params(which='minor', width=0.75)
    ax.tick_params(which='minor', length=2.5)
    ax.tick_params(axis='both', which='both', bottom='off', top='off',labelbottom='on', left='off', right='off', labelleft='on')
    
    ax.set_xlim(0, 5)
    ax.set_ylim(0, 1)
    ax.patch.set_alpha(0.0)


def barplot(plotarr,arrlabel,xlabel="x",ylabel="y",width=0.5,title="barplot",picname="test",color='gray',log=0,horizontal=0,rotation=0,):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    y_pos = np.arange(len(plotarr))
    if not horizontal:
        if log:ax.set_yscale('log')
        rects = ax.bar(y_pos,plotarr,align='center',color=color,log=log,width=width)
        ax.set_xticks(y_pos)
        ax.set_xticklabels(arrlabel,rotation=rotation)
        autolabel(rects,horizontal)
        ax.set_ylim(top=ax.axis()[-1]*1.1)
    else:
        if log:
            ax.set_xscale('log')
        rects = ax.barh(y_pos, plotarr, align='center',color=color,log=log,)
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

def bargroup(matrix,samplenames,ticklabels,xlabel="x",ylabel="y",title="bargroupplot",picname="test",colors=['gray','red','yellow'],width = 0.35,log=0,horizontal=0,rotation=0,stack=0,norm=0):
    # if you want to get same height for every bar, enable norm equals 1, will override parameters stack
    mat = np.asarray(matrix,dtype=np.float64)
    if log: 
        mat = np.log10(mat)
    r,c = mat.shape
    if norm:
        mat = np.divide(mat.T,np.sum(mat,axis=1)).T
        #mat = mat/np.sum(mat,axis=1)
        stack = 1
    try:
        assert c == len(samplenames)
    except AssertionError,e:
        sys.stderr.write("wrong data dimension!!!\n")
    if len(colors) != c:
        colors,markers = styleNum(c)
    fig, ax = plt.subplots()
    ind = np.arange(r)
    bottoms = np.zeros(r)
    lefts = np.zeros(r)
    for i in range(c):
        if stack:
            if horizontal:
                rects=ax.barh(ind, mat[:,i], height = width*r, align='edge',color=colors[i],label=samplenames[i],left=lefts)
                lefts += mat[:,i]
            else:
                rects=ax.bar(ind, mat[:,i], width*r, align='edge',color=colors[i],label=samplenames[i],bottom=bottoms)
                bottoms += mat[:,i]
            if not norm:
                autolabel(rects,horizontal)
        else:
            if horizontal:
                rects=ax.barh(ind+(i+1)*width, mat[:,i], height=width, align='edge', color=colors[i],label=samplenames[i])
                autolabel(rects,horizontal)
            else:
                rects=ax.bar(ind+(i+1)*width, mat[:,i], width, align='edge', color=colors[i],label=samplenames[i])
                autolabel(rects,horizontal)
    
    if not horizontal:
        ax.set_xticks(ind + width*r/2)
        ax.set_xticklabels(ticklabels,rotation=rotation)
        ax.set_ylim(top=ax.axis()[-1]*1.1)
    else:
        ax.set_yticks(ind + width*r/2) 
        ax.set_yticklabels(ticklabels,rotation=rotation)
        ax.set_xlim(right=ax.axis()[1]*1.1)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ncol = c/2 if c>6 else c
    plt.legend(loc='best',ncol=ncol,bbox_to_anchor=[0.6,0.95,0.4,0.1])   #The extend of the bounding box is zero - being equivalent to bbox_to_anchor=[x0, y0, 0, 0]
    plt.savefig("%s.png"%picname,format='png',dpi=300)
    plt.savefig("%s.svg"%picname,format='svg',dpi=300)
    plt.clf()
    plt.close()
    return 0



def clearaxis(ax):
    ''' remove xaxis and yaxis of the ax '''
    ax.margins(0.1)
    ax.axis('off')

def clearHalfaxis(ax):
    '''Only show ticks on the left and bottom spines. '''
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')


def boxplot(matrix,labels,xlabel="x",ylabel='y',title='boxplot',picname="test",rotation=0,patch_artist=True,notch=False):
    mat = np.asarray(matrix)
    r,c = mat.shape
    print r,c
    try:
        assert c == len(labels)
    except AssertionError,e:
        sys.stderr.write("wrong data dimension!!!\n")
    fig, ax = plt.subplots()
    #mat=np.asarray(matrix).T
    boxes=ax.boxplot(mat,sym='k+',bootstrap=5000,notch=notch,patch_artist=patch_artist)
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

def densityplot(plotarrs,labels,xlabel="x",ylabel='y',title='density plot',picname="test",fill=0):
    plotarrs = np.asarray(plotarrs)
    r,c = plotarrs.shape
    try:
        assert c == len(labels)
    except AssertionError,e:
        sys.stderr.write("inconsitent data and labels!!!\n")
    fig, ax = plt.subplots()
    for i,label in enumerate(labels):
        data=plotarrs[i]
        data=np.asarray(data)
        density = stats.kde.gaussian_kde(data)
        x=np.linspace(data.min(),data.max(), num=100)
        ax.plot(x, density(x),color=myfavcolors[i],label=label)
        if fill:
            ax.fill_between(x,0,density(x),color=myfavcolors[i],alpha=0.2)
    if len(labels) <= 6:
        plt.legend(loc='best')
    else:
        ncol=min(int(len(labels)/2.0),5)
        plt.legend(loc='upper center',bbox_to_anchor=(.5,1.1),ncol=ncol,fancybox=True,frameon=True)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
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
    fig, ax = plt.subplots(figsize=figsize)
    c = len(xarr)
    try:
        assert c == len(samplenames)
    except AssertionError,e:
        sys.stderr.write("inconsitent data dimension with length of samplenames!!!\n")
    if len(colors) != c:
        colors,markers = styleNum(c)
    if not marker:
        markers = ['o',]*c
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


def scatter3Dplot(xarr,yarr,zarr,samplenames,colors=['C0','C1','C2'],marker='o',xlabel="x",ylabel='y',zlabel='z',title='scatter plot',picname="scatter3d",alpha=0.8):
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax  = fig.add_subplot(111, projection='3d')
    c = len(xarr)
    try:
        assert c == len(samplenames)
    except AssertionError,e:
        sys.stderr.write("inconsitent data dimension with length of samplenames!!!\n")
    if len(colors) != c:
        sys.stderr.write("[WARNING] you should give correct colors parameters, will assign colors automatically!\n")
        colors,markers = styleNum(c)
    xlim = np.max(np.fabs(np.asarray(xarr)))
    ylim = np.max(np.fabs(np.asarray(yarr)))
    zlim = np.max(np.fabs(np.asarray(zarr)))
    for i in range(c):
        ax.scatter(xarr[i], yarr[i], zarr[i], c=colors[i], marker=marker, alpha=0.8,label=samplenames[i])
    ax.legend(loc='best')
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)
    ax.xaxis.set_major_locator(ticker.LinearLocator(6))
    ax.yaxis.set_major_locator(ticker.LinearLocator(6))
    ax.zaxis.set_major_locator(ticker.LinearLocator(6))
    majorFormatter = ticker.FormatStrFormatter('%d')
    ax.xaxis.set_major_formatter(majorFormatter)
    ax.yaxis.set_major_formatter(majorFormatter)
    ax.zaxis.set_major_formatter(majorFormatter)
    ax.tick_params(axis=u'both', which=u'both',length=20,color='red',pad=1)
    ax.spines['right'].set_color('none')
    plt.savefig("%s.png"%picname,format='png',dpi=300)
    plt.savefig("%s.svg"%picname,format='svg',dpi=300)
    plt.clf()
    plt.close()
    return 0


def pieplot(data,labels,cutoff=0.01,title='pie plot',pctdistance=1.1, startangle=0, picname='pie_plot', expidx=None):
    c = len(data)
    try:
        assert c == len(labels)
    except AssertionError,e:
        sys.stderr.write("inconsitent data dimension with length of laebls!!!\n")
    fig, ax = plt.subplots(figsize=(8,6))
    explode = np.zeros(c)
    if expidx:
        explode[[expidx]] = 0.1
    data   = np.asarray(data,dtype=np.float64)
    labels = np.asarray(labels)
    data   = data/np.sum(data)
    plotdata = data[np.where(data > cutoff)]
    explode  = explode[np.where(data > cutoff)]
    labels   = labels[np.where(data > cutoff)]
    colors,markers = styleNum(len(explode))
    patches,text1,text2 = ax.pie(plotdata, explode=explode, colors=colors, autopct='%1.2f%%', pctdistance=pctdistance, shadow=False, labeldistance=1.1, startangle=startangle, radius=None, counterclock=True,)
    ax.set_title(title)
    plt.legend(patches, labels, bbox_to_anchor=(0.68,0.9),bbox_transform=plt.gcf().transFigure)
    plt.gca().axis('equal')
    draw_pie(ax,X=2,Y=2)
    #plt.subplots_adjust(left=0.05, bottom=0.1, right=0.65)
    #plt.tight_layout()
    plt.savefig("%s.png"%picname,bbox_inches="tight",format='png',dpi=300)
    plt.savefig("%s.svg"%picname,bbox_inches="tight",format='svg',dpi=300)
    plt.clf()
    plt.close()
    return 0


def draw_pie(ax,ratios=[0.4,0.3,0.3], X=0, Y=0, size = 10):
    ''' referred Thomas Lecocq's blog about how to add a float pie to a particular axes'''
    N = len(ratios)
    colors,markers = styleNum(N)
    xy = []
    start = 0.
    for ratio in ratios:
        x = [0] + np.cos(np.linspace(2*np.pi*start,2*np.pi*(start+ratio), 30)).tolist()
        y = [0] + np.sin(np.linspace(2*np.pi*start,2*np.pi*(start+ratio), 30)).tolist()
        xy1 = zip(x,y)
        xy.append(xy1)
        start += ratio
    for i, xyi in enumerate(xy):
        ax.scatter([X],[Y] , marker=(xyi,0), s=size, facecolor=colors[i] )
    #a list of (x, y) pairs used for Path vertices. The center of the marker is located at (0,0) and the size is normalized.
    return 0


def scatterHistplot(xarr,yarr,colors='C1',marker='o',xlabel="x",ylabel='y',title='scatter hist plot',bins=100,picname="scatter",alpha=0.8,figsize=(8,6)):
    # support only 1-dimension plot data 
    fig = plt.figure(figsize=figsize)
    nullfmt = ticker.NullFormatter()         # no labels
    # definitions for the axes
    xarr = np.asarray(xarr)
    yarr = np.asarray(yarr)
    xlim = np.max(np.fabs(xarr))*1.1
    ylim = np.max(np.fabs(yarr))*1.1
    axscatter = plt.subplot2grid((4,4),(1,0),colspan=3,rowspan=3)
    axscatter.scatter(xarr, yarr, c=colors, marker=marker, alpha=0.8,s=14)
    histx = plt.subplot2grid((4,4),(0,0),colspan=3)
    histy = plt.subplot2grid((4,4),(1,3),rowspan=3)
    histx.xaxis.set_major_formatter(nullfmt)
    histy.yaxis.set_major_formatter(nullfmt)

    histx.hist(xarr, bins=bins, color=colors)
    histy.hist(yarr, bins=bins, color=colors, orientation='horizontal')
    axscatter.set_xlim((-xlim,xlim))
    axscatter.set_ylim((-ylim,ylim))
    histx.set_xlim(axscatter.get_xlim())
    histy.set_ylim(axscatter.get_ylim())
    histx.spines['left'].set_smart_bounds(True)
    histy.spines['bottom'].set_smart_bounds(True)
    #axscatter.set_title(title)
    axscatter.set_xlabel(xlabel)
    axscatter.set_ylabel(ylabel)
    plt.tight_layout(h_pad=0.5,w_pad=0.5)
    plt.savefig("%s.png"%picname,format='png',dpi=300)
    plt.savefig("%s.svg"%picname,format='svg',dpi=300)
    plt.clf()   
    plt.close()
    return 0


def venn2plot(listA,listB,labels=('A','B'),labelfz=20,numfz=12,picname='venn2'):
    # this function is just a template, the circle size is unchanged by number points, maybe I'll modify it to make it more verstile later 
    plt.style.use(u'classic')
    fig = plt.figure()
    ax = fig.add_axes([0 ,0 , 1, 1])
    setA = set(listA)
    setB = set(listB)
    AB = setA&setB
    Ab, AB, aB = map(list,[setA-AB,AB,setB-AB])
    lenAb, lenAB, lenaB = map(len,[Ab, AB, aB])
    import matplotlib.patches as mpatches
    #from matplotlib.collections import PatchCollection
    circle = mpatches.Circle((0.34,0.5), radius=0.3, ec="none",fc='C1',alpha=0.3)
    ax.add_patch(circle)
    circle = mpatches.Circle((0.66,0.5), radius=0.3, ec="none",fc='C2',alpha=0.3)
    ax.add_patch(circle)
    ax.text(0.2,0.5,str(lenAb),fontsize=numfz)
    ax.text(0.77,0.5,str(lenaB),fontsize=numfz)
    ax.text(0.5,0.5,str(lenAB),fontsize=numfz)
    ax.text(0.28,0.7,labels[0],fontsize=labelfz)
    ax.text(0.71,0.7,labels[1],fontsize=labelfz)
    ax.set_aspect(1)
    ax.axis('off')
    plt.savefig("%s.png"%picname,format='png',dpi=300)
    plt.savefig("%s.svg"%picname,format='svg',dpi=300)
    plt.clf()
    plt.close()
    return 0

def venn3plot(listA,listB,listC,labels=('A','B','C'),labelfz=18,numfz=12,picname='venn3'):
    plt.style.use(u'classic')
    fig = plt.figure()
    ax = fig.add_axes([0 ,0 , 1, 1])
    setA = set(listA)
    setB = set(listB)
    setC = set(listC)
    AB = setA&setB
    AC = setA&setC
    BC = setB&setC
    ABC = setA&setB&setC
    aBc = (setB - AB - BC)
    Abc = (setA - AB - AC)
    abC = (setC - BC - AC)
    import matplotlib.patches as mpatches
    circle = mpatches.Circle((0.5,0.7), radius=0.25, ec="none",fc='C1',alpha=0.3)
    ax.add_patch(circle)
    circle = mpatches.Circle((0.35,0.45), radius=0.25, ec="none",fc='C2',alpha=0.3)
    ax.add_patch(circle)
    circle = mpatches.Circle((0.65,0.45), radius=0.25, ec="none",fc='C3',alpha=0.3)
    ax.add_patch(circle)
    labelpos = [[0.8,0.85],[0.05,0.2],[0.85,0.2]]
    addvenntext(labelpos,labels,labelfz)
    textpos = [[0.2,0.4],[0.5,0.8],[0.75,0.4],[0.5,0.53]]
    texts = map(str,map(len,[aBc,Abc,abC,abC]))
    addvenntext(textpos,texts,numfz)
    # add center point
    ax.text(0.5,0.53,str(len(ABC)),fontsize=numfz)
    textpos = [[0.35,0.6],[0.475,0.4],[0.625,0.6]]
    texts = map(str,map(len,[AB-ABC,BC-ABC,AC-ABC]))
    addvenntext(textpos,texts,numfz)
    ax.set_aspect(1)
    ax.axis('off')
    plt.savefig("%s.png"%picname,format='png',dpi=300)
    plt.savefig("%s.svg"%picname,format='svg',dpi=300)
    plt.clf()
    plt.close()
    return 0

def addvenntext(labelpos,labels,fontsize):
    ax=plt.gca()
    for (x,y),label in zip(labelpos,labels):
        ax.text(x,y,label,fontsize=fontsize)
    return ax 

def venn4plot(listA,listB,listC,listD,labels=('A','B','C','D'),labelfz=18,numfz=12,r=0.28,picname='venn4'):
    plt.style.use(u'classic')
    fig = plt.figure()
    ax = fig.add_axes([0 ,0 , 1, 1])
    A = set(listA)
    B = set(listB)
    C = set(listC)
    D = set(listD)
    AB = A&B; AC = A&C; AD = A&D; BC = B&C; BD = B&D; CD = C&D
    ABC = A&B&C; ABD = A&B&D; ACD = A&C&D; BCD = B&C&D
    ABCD = A&B&C&D
    Abcd = A - AB - AC - AD
    aBcd = B - AB - BC - BD
    abCd = C - AC - BC - CD
    abcD = D - AD - BD - CD
    ABcd = AB - C - D; AbCd = AC - B - D;aBcD = BD - A - C; abCD = CD - A - B
    ABCd = ABC - D; ABcD = ABD - C; AbCD = ACD - B; aBCD = BCD - A
    import matplotlib.patches as mpatches
    circleA = mpatches.Circle((0.50,0.7), radius=r, ec="none",fc='C1',alpha=0.3)
    ax.add_patch(circleA)
    circleB = mpatches.Circle((0.30,0.5), radius=r, ec="none",fc='C2',alpha=0.3)
    ax.add_patch(circleB)
    circleC = mpatches.Circle((0.70,0.5), radius=r, ec="none",fc='C3',alpha=0.3)
    ax.add_patch(circleC)
    circleD = mpatches.Circle((0.50,0.3), radius=r, ec="none",fc='C4',alpha=0.3)
    ax.add_patch(circleD)
    # add label 
    labelpos = [[0.25,0.93],[0.1,0.2],[0.93,0.73],[0.75,0.1]]
    addvenntext(labelpos,labels,labelfz)
    # add element points
    textpos = [[0.5,0.85],[0.15,0.5],[0.85,0.5],[0.5,0.15]]
    texts = map(str,map(len,([Abcd,aBcd,abCd,abcD])))
    addvenntext(textpos,texts,numfz)
    textpos = [[0.65/2,1.35/2],[1.35/2,1.35/2],[1.35/2,0.65/2],[0.65/2,0.65/2]]
    texts = map(str,map(len,[ABcd,AbCd,abCD,aBcD]))
    addvenntext(textpos,texts,numfz)
    textpos = [[0.5,0.6],[0.5,0.38],[0.38,0.5],[0.6,0.5]]
    texts = map(str,map(len,[ABCd,aBCD,AbCD,AbCD]))
    addvenntext(textpos,texts,numfz)
    # add center
    ax.text(0.5,0.5,str(len(ABCD)))
    ax.set_aspect(1)
    ax.axis('off')
    plt.savefig("%s.png"%picname,format='png',dpi=300)
    plt.savefig("%s.svg"%picname,format='svg',dpi=300)
    plt.clf()
    plt.close()
    return 0


def plot_chromosome():
    pass
