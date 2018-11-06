'''This function returns the index of the grid point just below xi and the contributions to p(Xj) and p(Xj+1) for the i particule'''

def barycentre(deltax,xi):
 j = int(xi/deltax)
 alpha = (xi-j*deltax)/deltax
 return j,1-alpha,alpha
 
print(barycentre(0.5,0.75))