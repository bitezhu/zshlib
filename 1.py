import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


axes=[]
ax1 = plt.subplot2grid((1,10),(0,0),colspan=4,rowspan=1)
ax1.plot([0,1],[2,3])
ax2 = plt.subplot2grid((1,10),(0,4),colspan=1,rowspan=1)
ax2.plot([1,2],[3,4])
ax3 = plt.subplot2grid((1,10),(0,5),colspan=3,rowspan=1)
ax3.plot([2,3],[4,5])
ax4 = plt.subplot2grid((1,10),(0,8),colspan=2,rowspan=1)
ax4.plot([3,4],[5,6])
axes=[ax1,ax2,ax3,ax4]

ax1.spines['right'].set_visible(False)
ax1.set_xticks([0,1])
ax1.set_xticklabels(['0','1'])
ax2.spines['right'].set_visible(False)
ax2.spines['left'].set_visible(False)
ax2.yaxis.set_major_locator(ticker.NullLocator())
ax2.set_xticks([2])
ax2.set_xticklabels(['2'])
ax3.spines['right'].set_visible(False)
ax3.spines['left'].set_visible(False)
ax3.yaxis.set_major_locator(ticker.NullLocator())
ax3.set_xticks([3])
ax3.set_xticklabels(['3'])
ax4.spines['left'].set_visible(False)
ax4.yaxis.set_major_locator(ticker.NullLocator())
ax4.set_xticks([4])
ax4.set_xticklabels(['4'])

[plt.setp(axes[i],xlim=[i+0,i+1]) for i in range(4)]
[plt.setp(axes[i],ylim=[2,6]) for i in range(4)]
plt.subplots_adjust(wspace=0,)

plt.savefig('xx.png',format='png',dpi=300)
