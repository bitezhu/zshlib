import numpy as np
import matplotlib
import sys

#matplotlib.use['qt4agg']
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import stats
from itertools import *


#myfavcolors = ['grey','tomato','blue','green','gold','black','deepskyblue','olive','steelblue','brown','plum','chocolate','cyan','crimson','aqua']
myfavcolors=['grey','tomato','peru','dodgerblue','brown','darkslategray','lightsalmon','limegreen','deeppink','mediumpurple','darkkhaki','aqua','coral','gold','cornflowerblue','crimson']
Markers = Line2D.markers
Fillmarkers = Line2D.filled_markers
myfavmarkers = [Markers['+'],Markers[2],Markers['x'],Markers[5],Markers['.'],Fillmarkers[0],Fillmarkers[2],Fillmarkers[-1],Fillmarkers[-6],Fillmarkers[-5],Markers['D']]

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


def bargroup(matrix,samplenames,ticklabels,xlabel="x",ylabel="y",title="bargroupplot",picname="test",colors=['gray','red','yellow'],width = 0.35,log=0,horizontal=0):
    mat = np.asarray(matrix)
    r,c = mat.shape
    print r,c
    try:
        assert c == len(samplenames)
    except AssertionError,e:
        sys.stderr.write("wrong data dimension!!!\n")
    if len(colors) != c:
        colors,marker = styleNum(c)
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

def clearaxis(ax):
    ''' remove xaxis and yaxis of the ax '''
    ax.margins(0.2)
    ax.set_axis_off()


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
    colors,marker = styleNum(len(labels))
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

