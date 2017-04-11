# -- load the raw grids
dgrid = np.load('../output/lwir_distgrid_raw.npy')
xgrid = np.load('../output/lwir_xgrid_raw.npy')
ygrid = np.load('../output/lwir_ygrid_raw.npy')

# -- loop through rows and columns and replace further distances
nrow, ncol = dgrid.shape
print "joining top to bottom"
for ii in range(1,nrow-1):
    dgrid[0,:] = dgrid.max() # avoid anomalous hits in top row
    for jj in range(1,ncol-1):
        if (dgrid[ii-1,jj-1] == dgrid[ii-1,jj] == dgrid[ii-1,jj+1] ==
            dgrid[ii,jj+1] == dgrid[ii+1,jj+1] == dgrid[ii+1,jj] == 
            dgrid[ii+1,jj-1] == dgrid[ii,jj-1]):
            dgrid[ii,jj] = dgrid[ii-1,jj-1]
            continue
        if dgrid[ii,jj]>dgrid[ii-1,jj]:
            dgrid[ii,jj]=dgrid[ii-1,jj]
            xgrid[ii,jj]=xgrid[ii-1,jj]
            ygrid[ii,jj]=ygrid[ii-1,jj]

#print "smoothing"
#for ii in range(nrow):
#    for jj in range(1,ncol-1):
#        if (dgrid[ii,jj-1]==dgrid[ii,jj+1]) and \
#                (dgrid[ii,jj-1]!=dgrid[ii,jj]):
#            dgrid[ii,jj]=dgrid[ii,jj-1]
#            xgrid[ii,jj]=xgrid[ii,jj-1]
#            ygrid[ii,jj]=ygrid[ii,jj-1]

np.save('../output/lwir_distgrid.npy',dgrid)
np.save('../output/lwir_xgrid.npy',xgrid)
np.save('../output/lwir_ygrid.npy',ygrid)
